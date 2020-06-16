#!/usr/bin/env python
#-------------------------------------------------------------------------------------------------------------------
#  pvclient.py  -  client for the PseudoViewer web service (pseudoviewer.inha.ac.kr/)
#  Copyright Steven Busan 2014 (www.chem.unc.edu/rna/)
#  Included as part of the SHAPE-MaP analysis pipeline (ShapeMapper)
#
#   - Requires a secondary structure/list of structures in .ct format (max sequence length ~4000 nucleotides)
#   - Requests a rendered structure in postscript format from the PseudoViewer server,
#       cleans up the image, and writes to local file.  Also creates xrna files.
#   - Optional arguments for coloring by chemical probing reactivities or arbitrary nucleotide ranges
#
#  Quick usage examples -
#
#   render 1st 10 structures in .ct file and color by SHAPE reactivity file:
#        pvclient.py --ct folds.ct --shape 1M7.shape --structures 10
#
#
#   render 1st structure and color by differential reactivity file, hide title
#        pvclient.py --ct folds.ct --diff NMIA-1M6.diff --no_title
#
#
#   render all structures, highlight nucleotides 50-60 in magenta, 80-90 in yellow
#        pvclient.py --ct folds.ct --range 50 60 magenta 80 90 yellow --structures 9999
#
#
#  Full argument descriptions -
#
#  --ct <ctFilePath.ct>
#       Required file containing 1 or more secondary structures in connect-table format
#
#  --structures <number>
#       Number of structures in .ct file to render (default: 1)
#  --shape <shapeFilePath.shape>
#       Color nucleotides by reactivity. Red >= 0.85, orange >= 0.4, black < 0.4, gray <-998.5    
#  --diff <differenceFilePath.diff/.dif> [<upperColor> <lowerColor> <upperThreshold> <lowerThreshold>]
#       Color nucleotides by differential reactivity. Optionally specify the positive and
#       negative thresholds and colors (default: green blue 0.5 -0.5)
#  --range <startNumber> <endNumber> <color>
#       Color nucleotides in range (inclusive) given color.
#       Multiple ranges and colors may be specified - see example above
#  --out <destinationPath>
#       Write files to the location specified.  The structure number will be appended to
#       the end of the filename. New folders will not be created. (default: "structure_")
#  --title <title>
#       Write the given string to each postscript image in the upper-left corner instead of
#       the name pulled from the .ct file
#  --no_title
#       Do not add any title to the rendered images
#
#  Available colors -
#       red blue green cyan magenta yellow orange purple brown darkred lightblue lightgreen pink teal
#
#-------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# GPL statement:
#
# This file is part of Shapemapper.
#
# ShapeMapper is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ShapeMapper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ShapeMapper.  If not, see <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------------------------------

DEBUG = False

import re, sys, os, traceback, argparse, httplib2, RNAtools

# display the header of this file if executed with no arguments
showHelp = False
if len(sys.argv) == 1:
    showHelp = True
elif sys.argv[1] == "-h" or sys.argv[1] == "--help":
    showHelp = True
if showHelp == True:
    thisFile = open(os.path.realpath(__file__),"r")
    inBlock = False
    endBlock = False
    while endBlock == False:
        line = thisFile.readline().replace("#","").rstrip()
        if inBlock == True:
            print line
            if line.find("-------------------------------") != -1:
                endBlock = True
        elif line.find("--------------------------------") != -1:
            print line
            inBlock = True
    sys.exit(1)

class ConnectTableFile(argparse.FileType):
    def __call__(self, string):
        allowedExtensions = ['.ct','.CT']
        base, ext = os.path.splitext(string)
        if ext not in allowedExtensions:
            raise ValueError('\'%s\' is not a recognized connect-table file extension.'%ext)
        return super(ConnectTableFile, self).__call__(string)

class ShapeFile(argparse.FileType):
    def __call__(self, string):
        allowedExtensions = ['.shape','.SHAPE']
        base, ext = os.path.splitext(string)
        if ext not in allowedExtensions:
            raise ValueError('\'%s\' is not a recognized difference file extension.'%ext)
        return super(ShapeFile, self).__call__(string)

def parseDiffFilePath(arg):
    allowedExtensions = ['.diff','.DIFF','.dif','.DIF']
    base, ext = os.path.splitext(arg)
    if ext in allowedExtensions:
        return open(arg, 'r')
    else:
        return None

#class DotBracketFile(argparse.FileType):
#    def __call__(self, string):
#        allowedExtensions = ['.db','.DB','.dbn','.DBN','.dp','.DP']
#        base, ext = os.path.splitext(string)
#        if ext not in allowedExtensions:
#            raise ValueError('\'%s\' is not a recognized dot-bracket file extension.'%ext)
#        return super(DotBracketFile, self).__call__(string)

# define color names and associated postscript strings
colorDict = {
    'gray':"0.6 0.6 0.6 setrgbcolor",
    'red':"1 0 0 setrgbcolor",
    'blue':"0.3 0.3 0.9 setrgbcolor",
    'green':"0 0.75 0 setrgbcolor",
    'magenta':"0.9 0 0.9 setrgbcolor",
    
    'yellow':"0.8 0.8 0 setrgbcolor",
    'orange':"0.95 0.5 0.1 setrgbcolor",
    'lightblue':"0.6 0.6 1 setrgbcolor",
    'black':"0 0 0 setrgbcolor",
    'purple':"0.65 0.0 0.8 setrgbcolor",
    
    'teal':"0 0.8 0.6 setrgbcolor",
    'cyan':"0 0.9 0.9 setrgbcolor",
    'lightgreen':"0.4 1 .4 setrgbcolor",
    'brown':"0.7 0.55 0.2 setrgbcolor",
    'pink':"1 0.6 0.7 setrgbcolor",
    'darkred':"0.75 0 0 setrgbcolor"}

# default colors for differential reactivity coloring
defaultUpperColors = ['yellow','green']
defaultLowerColors = ['magenta','blue']

xrnaMainTemplate = """<ComplexDocument Name='#name'>
<SceneNodeGeom CenterX='#centerX' CenterY='#centerY' />
<LabelList>
s #titleX #titleY 0.0 16 2 0 "#name"
</LabelList>
<Complex Name='#name'> 
<RNAMolecule Name='#name'> 
<NucListData StartNucID='#startNuc' DataType='NucChar.XPos.YPos'>
#nucData
</NucListData>
#nucColors
<Nuc RefIDs='#startNuc-#endNuc' IsNucPath='true' NucPathColor='dddddd' NucPathLineWidth='1.5' />
#nucNumbers
</RNAMolecule>
</Complex>
</ComplexDocument>"""

xrnaHelixTemplate = """<BasePairs nucID='%i' length='%i' bpNucID='%i' />"""
xrnaNucTemplate = """%s %f %f"""
xrnaColorTemplate = """<Nuc RefIDs='%i-%i' Color='%02x%02x%02x' FontID='2' FontSize='11' />"""

def parseColor(arg):
    if arg in colorDict.keys():
        return arg
    else:
        return None

def parseDiffCutoff(arg):
    try:
        cutoff = float(arg)
        return cutoff
    except ValueError:
        return None

def parseRangeNum(arg):
    try:
        num = int(arg)
        return num
    except ValueError:
        return None
    
class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("\n\nERROR: Command-line arguments are incorrectly formatted. Run pvclient with no args for help.\n\n")
        sys.exit(2)

def parseStructArgs():
    parser = CustomParser(description="Render RNA secondary structures.")
    structPathGroup = parser.add_mutually_exclusive_group(required=True)
    #structPathGroup.add_argument("--db", type=DotBracketFile('r'), help="Dot-bracket secondary structure.")
    structPathGroup.add_argument("--ct", type=ConnectTableFile('rU'), help="Connect-table secondary structure.")
    parser.add_argument("--structures", type=int, default=1)
    parser.add_argument("--out", default="structure_")
    parser.add_argument("--title", default="")
    parser.add_argument("--no_title", default=False, action="store_true")
    # accept either a SHAPE file or 1 or more difference files
    #colorPathGroup = parser.add_mutually_exclusive_group()
    colorPathGroup = parser.add_argument_group()
    colorPathGroup.add_argument("--shape", type=ShapeFile('r'), help="Chemical-probing data.")
    colorPathGroup.add_argument("--diff", nargs='+', help="<differenceFile.diff> [<upperColor> <lowerColor> <upperThreshold> <lowerThreshold>] (repeat)")
    # color a range of nucleotides (inclusive) a given color
    colorPathGroup.add_argument("--range", nargs='+')
    #args = parser.parse_args(argString.split())
    args = parser.parse_args()

    # no easy way I can find to handle repeated groups of positional args within argparse, so do it "manually" here
    if args.range:
        args.ranges = []
        args.rangeColors = []
        i = 0
        while i < len(args.range):
            args.ranges.append([parseRangeNum(args.range[i]), parseRangeNum(args.range[i+1])])
            args.rangeColors.append(parseColor(args.range[i+2]))
            i += 3
        del args.range
    
    if args.diff:
        args.diffFiles = []
        args.diffUpperColors = []
        args.diffLowerColors = []
        args.diffUpperCutoffs = []
        args.diffLowerCutoffs = []

        for arg in args.diff:
            diffFile = parseDiffFilePath(arg)
            color = parseColor(arg)
            cutoff = parseDiffCutoff(arg) 
                
            if diffFile:
                args.diffFiles.append(diffFile)
                findUpperColor = True
                findUpperCutoff = True
            elif color:
                if findUpperColor == True:
                    args.diffUpperColors[len(args.diffFiles)-1] = color
                    findUpperColor = False
                else:
                    args.diffLowerColors[len(args.diffFiles)-1] = color
            elif cutoff:
                if findUpperCutoff == True:
                    args.diffUpperCutoffs[len(args.diffFiles)-1] = cutoff
                    findUpperCutoff = False
                else:
                    args.diffLowerCutoffs[len(args.diffFiles)-1] = cutoff
            else:
                raise ValueError('Difference file arguments are incorrectly formatted.')

            # set default colors
            if len(args.diffFiles) > len(args.diffUpperColors):
                args.diffUpperColors.append(defaultUpperColors.pop())
            if len(args.diffFiles) > len(args.diffLowerColors):
                args.diffLowerColors.append(defaultLowerColors.pop())

            # set default thresholds
            if len(args.diffFiles) > len(args.diffUpperCutoffs):
                args.diffUpperCutoffs.append(0.5)
            if len(args.diffFiles) > len(args.diffLowerCutoffs):
                args.diffLowerCutoffs.append(-0.5)
                
        del args.diff

    #print str(args)
    return args

def ct_to_dbn(lines, maxStructures):

    linesToWrite = []

    leftBracket = "([{<"
    rightBracket = ")]}>"

    titleFinder = re.compile(r'[0-9]+')

    i = 0
    structCount = 0
    titleFound = False
    while i < len(lines) and structCount < maxStructures:
        line = lines[i]
        index = 0
        if titleFound == False:
            try:
                firstWord = line.strip().split()[0]
                numberNucs = int(firstWord)
                titleFound = True
            except:
                pass
        if titleFound == True:
            structCount += 1
            newLine = ">   "+line[index:].strip()+"\n"
            linesToWrite.append(newLine)
            splitLine = line.strip().split()
            seq = list("_"*numberNucs)
            struct = list("."*numberNucs)
            pairedNuc = [0]*numberNucs
            for n in range(0,numberNucs):
                lineIndex = n+i+1
                splitLine = lines[lineIndex].strip().split()
                seq[n] = splitLine[1]
                pairedNuc[n] = int(splitLine[4])
            for n in range(0,numberNucs):
                if pairedNuc[n] == 0:
                    struct[n] = "."
                elif pairedNuc[n]-1 > n:
                    if struct[n] == ".":
                        struct[n] = "("
                        struct[pairedNuc[n]-1] = ")"
                    bracketIndex = leftBracket.find(struct[n])
                    partner = pairedNuc[n]-1
                    # search for non-nested basepairs
                    for k in range(n+1,partner):
                        if pairedNuc[k] != 0:
                            subPartner = pairedNuc[k]-1
                            if subPartner > partner:
                                # pseudoknot detected
                                struct[k] = leftBracket[bracketIndex+1]
                                struct[subPartner] = rightBracket[bracketIndex+1]
            linesToWrite.append("".join(seq)+"\n")
            linesToWrite.append("".join(struct)+"\n")
            i += n+1
            titleFound = False
        i += 1
    return linesToWrite


def parseSHAPE(file):
    readLines = file.read()
    file.seek(0) # return to beginning of file (this is an ugly kludge - I read the same file multiple times)
    lines = readLines.splitlines()

    # parse SHAPE reactivities
    reactivityDict = {}
    for line in lines:
        splitLine = line.strip().split()
        nucNum = int(splitLine[0])
        reactivity = float(splitLine[1])
        reactivityDict[nucNum] = reactivity
    print "loaded SHAPE reactivities:"
    if DEBUG==True:
        for n in reactivityDict:
            print "%i: %0.2f"%(n,reactivityDict[n])
    return reactivityDict

def parseDIFF(file, positiveThreshold, negativeThreshold):
    readLines = file.read()
    file.seek(0)
    lines = readLines.splitlines()

    # parse differential reactivities
    reactivityDict = {}
    for line in lines:
        splitLine = [field.strip() for field in line.split("\t")]
        nucNum = int(splitLine[0])
        if splitLine[1] == "":
            reactivity = -999.0
        else:
            reactivity = float(splitLine[1])
        if(reactivity >= positiveThreshold):
            reactivityDict[nucNum] = "above"
        elif reactivity <= negativeThreshold and reactivity >= -998.5:
            reactivityDict[nucNum] = "below"
        elif reactivity < -998.5:
            reactivityDict[nucNum] = "no_data"
        else:
            reactivityDict[nucNum] = "within"
    return reactivityDict

def parseDotBracket(lineGroup):
    title = lineGroup[0].strip()[1:]
    seq = lineGroup[1].strip()
    struct = lineGroup[2].strip()

    leftBracket = "([{<"
    rightBracket = ")]}>"

    pairedNuc = [0]*len(struct)
    openList = [[] for i in range(len(leftBracket))]
    for i in range(len(struct)):
        pairedNuc[i] = -1
        if struct[i] in leftBracket:
            openList[leftBracket.find(struct[i])].append(i)
        elif struct[i] in rightBracket:
            pairedNuc[i] = openList[rightBracket.find(struct[i])].pop()
            pairedNuc[pairedNuc[i]] = i

    return title, seq, struct, pairedNuc

def parseDotBracketFile(lines):
    lineGroups = []
    lineGroup = []
    for line in lines:
        if len(line.strip()) > 0:
            if len(lineGroup) < 3:
                lineGroup.append(line)
            if len(lineGroup) == 3:
                lineGroups.append(lineGroup)
                lineGroup = []
    return lineGroups


def getEPS(title, startNuc, seq, struct):
    xmlTemplate = """<?xml version="1.0" encoding="utf-8"?>
    <soap:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                    xmlns:xsd="http://www.w3.org/2001/XMLSchema"
                    xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
    <soap:Body> 
    <WSPVRun xmlns="http://wilab.inha.ac.kr/WSPseudoViewer/">
    <WSPVRequest xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns="http://wilab.inha.ac.kr/WSPseudoViewer/">
    <Option Scale="1" Output_name="test" Output_type="raw eps"/>
    <WSPVIn_file_data>
    <PV_file>
    # %s
    %i
    %s
    %s
    </PV_file>
    </WSPVIn_file_data>
    </WSPVRequest>
    </WSPVRun>
    </soap:Body>
    </soap:Envelope>"""

    xmlRequest = xmlTemplate%(title,startNuc,seq,struct)

    URL = "http://165.246.44.42/WSPseudoViewer/WSPseudoViewer.asmx?WSDL"

    headers = { "Content-type": "text/xml; charset=utf-8",
                "Content-length": "%d" % len(xmlRequest),
                "Connection": "keep-alive",
                "SOAPAction": "\"http://wilab.inha.ac.kr/WSPseudoViewer/WSPVRun\""
                }

    # httplib2 method, handles cacheing/packet splitting properly
    http = httplib2.Http(timeout=30)
    response, content = http.request(URL, 'POST', headers=headers, body=xmlRequest)
    stringResponse = str(content)
    #print stringResponse
    # pull eps file from response
    epsRegex = re.compile(r'<WSPVOut_EPS>(.*?)</WSPVOut_EPS>', re.DOTALL)
    faultRegex = re.compile(r'<faultstring>(.*?)</faultstring>', re.DOTALL)
    match = epsRegex.search(stringResponse)
    faultMatch = faultRegex.search(stringResponse)
    epsResponse = ""
    if match != None:
        epsResponse = match.group(1)
    elif faultMatch != None:
        raise Exception("Pseudoviewer error: "+faultMatch.group(1))
    else:
        raise Exception("Pseudoviewer server returned: "+stringResponse)
        
    return epsResponse

def makeOrderedNucIndices(pairedNuc):
    # PseudoViewer's postscript draws helical nucs before unpaired nucs
    orderedNucIndices = []
    pairedIndices = []
    unpairedIndices = []
    for i in range(len(pairedNuc)):
        if pairedNuc[i] != -1:
            pairedIndices.append(i+1)
        else:
            unpairedIndices.append(i+1)
    orderedNucIndices = pairedIndices+unpairedIndices
    return orderedNucIndices

def makeShapeColorStrings(seq, startNuc, reactivityDict):
    # prepare list of colors for each nucleotide
    gray = "0.6 0.6 0.6 setrgbcolor"
    red = "1 0 0 setrgbcolor"
    rcomp = 245/255.0
    gcomp = 127/255.0
    bcomp = 32/255.0
    orange = "%.3f %.3f %.3f setrgbcolor"%(rcomp, gcomp, bcomp)
    black = "0 0 0 setrgbcolor"
    colorStrings = []
    for i in range(startNuc, len(seq)+startNuc):
        if i in reactivityDict.keys():
            if reactivityDict[i] < -998.5:
                colorStrings.append(gray)
            elif reactivityDict[i] <= 0.4:
                colorStrings.append(black)
            elif reactivityDict[i] <= 0.85:
                colorStrings.append(orange)
            else:
                colorStrings.append(red)
        else:
            colorStrings.append(gray)
    return colorStrings


def makeDiffColorStrings(seq, startNuc, differentialDicts, diffUpperColors, diffLowerColors):
    colorStrings = []
    for i in range(startNuc, len(seq)+startNuc):
        colorString = colorDict["black"]
        for j in range(len(differentialDicts)):
            differentialDict = differentialDicts[j]
            diffUpperColor = diffUpperColors[j]
            diffLowerColor = diffLowerColors[j]
            if i in differentialDict.keys():
                if differentialDict[i] == "above":
                    colorString = colorDict[diffUpperColor]
                elif differentialDict[i] == "within":
                    colorString = colorDict['black']
                elif differentialDict[i] == "below":
                    colorString = colorDict[diffLowerColor]
                elif differentialDict[i] == "no_data":
                    colorString = colorDict['gray']
        colorStrings.append(colorString)
    return colorStrings

def makeRangeColorStrings(seq, startNuc, args):
    colorStrings = []
    for i in range(startNuc, len(seq)+startNuc):
        colorString = colorDict['black']
        # this loop should really go outside the outer loop for speed, but oh well
        for j in range(len(args.ranges)):
            start = args.ranges[j][0]
            end = args.ranges[j][1]
            color = args.rangeColors[j]
            if i in range(start, end+1):
                colorString = colorDict[color]
        colorStrings.append(colorString)
    return colorStrings

def addRangeColorStrings(seq, startNuc, args, colorStrings):
    for i in range(startNuc, len(seq)+startNuc):
        for j in range(len(args.ranges)):
            start = args.ranges[j][0]
            end = args.ranges[j][1]
            color = args.rangeColors[j]
            if i in range(start, end+1):
                colorString = colorDict[color]
                colorStrings[i-1] = colorString
    return colorStrings

def makeSolidColorStrings(seq, startNuc):
    colorStrings = []
    for i in range(startNuc, len(seq)+startNuc):
        colorStrings.append(colorDict["black"])
    return colorStrings

def modifyEPS(epsResponse, seq, startNuc, colorStrings, title, pairedNuc, startNumFrom=1):
    # pretty up the drawing with title, colored reactivities, bigger fonts, etc.
    orderedNucIndices = makeOrderedNucIndices(pairedNuc)
    if DEBUG == True:
        print "orderedNucIndices: "+str(orderedNucIndices)
    #nucRegex = re.compile(r'\S+\s\S+\smoveto\s\([AUGCTNaugctn]\)\sshow')
    # nucleotide line-matching regex excludes numbers but allows other chars now
    nucRegex = re.compile(r'\S+\s\S+\smoveto\s\([^0-9]\)\sshow')
    numRegex = re.compile(r'\S+\s\S+\smoveto\s\([0-9]*\)\sshow')
    nucsFound = 0
    leftX = 0.0
    topY = 0.0
    bottomY = 0.0
    rightX = 0.0
    centerX = 0.0
    centerY = 0.0
    outLines = ""
    xrnaNucStrings = [""]*len(seq)
    for line in epsResponse.splitlines():
        #if "BoundingBox" in line:
        #    leftX = -float(line.split()[3])/2
        #    topY = float(line.split()[4])/2
        # add structure title to page
        if nucsFound == len(orderedNucIndices):
            outLines += "0 0 0 setrgbcolor\n"
            outLines += "/Helvetica-Bold findfont 16 scalefont setfont\n"
            outLines += "%.3f %.3f moveto (%s) show\n"%(leftX, topY, title)
            nucsFound += 1
        match = nucRegex.match(line)
        numMatch = numRegex.match(line)
        
        if match:
            outLines += "%% nucleotide %i\n"%(orderedNucIndices[nucsFound]+startNuc)
            splitLine = line.strip().split()
            nucIndex = orderedNucIndices[nucsFound]-1
            if nucIndex in range(len(xrnaNucStrings)):
                xrnaNucStrings[nucIndex] = "%s %s %s"%(splitLine[3][1],splitLine[0],splitLine[1])
            outLines += colorStrings[orderedNucIndices[nucsFound]-startNuc]+"\n"
            #print orderedNucIndices, nucsFound
            outLines += line.strip()+"\n"
            nucsFound += 1
            # find upper-left corner of image, since bounding box can't be trusted (fix bounding box in the future)
            currentX = float(line.split()[0])
            currentY = float(line.split()[1])
            if currentX < leftX:
                leftX = currentX
            if currentY > topY:
                topY = currentY
            # find bottom-right corner of image
            if currentX > rightX:
                rightX = currentX
            if currentY < bottomY:
                bottomY = currentY

        elif numMatch:
            splitLine = line.rstrip().split()
            oldNum = int(splitLine[3].replace("(","").replace(")",""))-1 + startNumFrom
            newLine = "{0[0]} {0[1]} {0[2]} ({1}) show\n".format(splitLine, oldNum)
            outLines += newLine
            
        # make background strokes white
        elif "setrgbcolor" in line \
        and "0 0 0 setrgbcolor" not in line \
        and "0.5 0.5 0.5 setrgbcolor" not in line:
            outLines += "1 1 1 setrgbcolor\n"
        # replace default font
        elif "findfont" in line:
            outLines += "/Helvetica-Bold findfont 10.5 scalefont setfont\n"
        else:
            outLines += line.strip()+"\n"
    centerX = -(leftX+(rightX-leftX)/2)
    centerY = -(topY+(bottomY-topY)/2)

    xrnaColorStrings = [xrnaColorTemplate]*len(seq)
    for i in range(len(seq)):
        splitColor = colorStrings[i].split()
        colorInts = [int(float(splitColor[n])*255) for n in range(3)]
        xrnaColorStrings[i] = xrnaColorStrings[i]%(i+startNuc,i+startNuc,colorInts[0],colorInts[1],colorInts[2])

    xrnaHelixStrings = [""]
    inHelix = False
    helixLength = 0
    helixStart = 0
    helixPair = 0
    # pairedNuc is 0-based, with -1 indicating no basepair
    for i in range(len(pairedNuc)):
        currentNuc = i+1
        currentPaired = pairedNuc[i]+1
        if currentPaired > currentNuc:
            if inHelix == False:
                inHelix = True
                helixStart = currentNuc
                helixPair = currentPaired
                helixLength = 1
                xrnaHelixStrings.append(xrnaHelixTemplate%(helixStart,helixLength,helixPair))
            else:
                if abs(currentNuc-lastNuc) == 1 and abs(currentPaired-lastPaired) == 1:
                    helixLength = currentNuc - helixStart + 1
                    xrnaHelixStrings[-1] = xrnaHelixTemplate%(helixStart,helixLength,helixPair)
                else:
                    xrnaHelixStrings[-1] = xrnaHelixTemplate%(helixStart,helixLength,helixPair)
                    xrnaHelixStrings.append("")
                    helixStart = currentNuc
                    helixPair = currentPaired
                    helixLength = 1
        else:
            inHelix = False
        lastNuc = currentNuc
        lastPaired = currentPaired
    xrnaOutLines = xrnaMainTemplate.replace("#name",title)
    xrnaOutLines = xrnaOutLines.replace("#startNuc",str(startNuc)).replace("#endNuc",str(startNuc+len(seq)-1))
    xrnaOutLines = xrnaOutLines.replace("#titleX",str(leftX)).replace("#titleY",str(topY))
    xrnaOutLines = xrnaOutLines.replace("#nucData","\n".join(xrnaNucStrings))
    xrnaOutLines = xrnaOutLines.replace("#nucColors", "\n".join(xrnaColorStrings))
    xrnaOutLines = xrnaOutLines.replace("#nucNumbers","\n".join(xrnaHelixStrings))
    xrnaOutLines = xrnaOutLines.replace("#centerX",str(centerX)).replace("#centerY",str(centerY))

    return outLines, xrnaOutLines


def makeColorStrings(seq, startNuc, args):
    colorStrings = None
    diffFilesFound = False
    try:
        if args.diffFiles:
            diffFilesFound = True
    except:
        pass
    shapeFileFound = False
    try:
        if args.shape:
            shapeFileFound = True
    except:
        pass
    rangeFound = False
    try:
        if args.ranges:
            rangeFound = True
    except:
        pass

    # color by differential reactivities
    if diffFilesFound:
        differentialDicts = []
        # convert differential reactivity file into a dictionary of strings indicating value with respect to thresholds:
        # "above", "below", "within", or "no_data"
        for i in range(len(args.diffFiles)):
            differentialDicts.append({})
            differentialDicts[-1] = parseDIFF(args.diffFiles[i], args.diffUpperCutoffs[i], args.diffLowerCutoffs[i])
        # generate colors based on nucs passing threshold
        colorStrings = makeDiffColorStrings(seq, startNuc, differentialDicts, args.diffUpperColors, args.diffLowerColors)

    # color by SHAPE reactivity
    elif shapeFileFound:
        reactivityDict = parseSHAPE(args.shape)
        colorStrings = makeShapeColorStrings(seq, startNuc, reactivityDict)

    # color by defined ranges
    #elif rangeFound:
    #    colorStrings = makeRangeColorStrings(seq, startNuc, args)

    # no color
    else:
        colorStrings = makeSolidColorStrings(seq, startNuc)

    # overwrite existing colors with range colors
    if rangeFound:
        print "updating colors with range colors"
        colorStrings = addRangeColorStrings(seq, startNuc, args, colorStrings)

    return colorStrings
    

def python_client(CT, shape, startNumFrom=1, writeFolder=""):
    startNuc = 1
    CT.writeCT("tmpCTfile.ct")
    with open("tmpCTfile.ct") as tmp:
        tmpCT = tmp.readlines()
        dbnLines = ct_to_dbn(tmpCT,2)
    os.remove("tmpCTfile.ct")
    
    
    shapeDict = {}
    for i in range(len(shape)):
        shapeDict[i+startNuc] = shape[i]
        
    
    #print dbnLines
    lineGroups = parseDotBracketFile(dbnLines)
    structureIndex = 1
    for lineGroup in lineGroups:
        fullTitle, seq, struct, pairedNuc = parseDotBracket(lineGroup)

        #print "Asking for structure %d . . ."%structureIndex
        #print "args.out: "+args.out
        epsResponse = getEPS(fullTitle, startNuc, seq, struct)
        if DEBUG==True:
            debugEPS = open("raw_response.eps","w")
            debugEPS.write(epsResponse)
            debugEPS.close()
        #print "Coloring image %d . . ."%structureIndex
        #colorStrings = makeColorStrings(seq, startNuc, args)
        colorStrings = makeShapeColorStrings(CT.seq, startNuc, shapeDict)
        if DEBUG==True:
            print "colorStrings: "+str(colorStrings)
        outLines, xrnaOutLines = modifyEPS(epsResponse, seq, startNuc, colorStrings, fullTitle, pairedNuc, startNumFrom)
        outputPath = os.path.normpath(writeFolder)
        #print "outputPath: "+outputPath
        outFile = open(outputPath+".eps", "w")
        outFile.write(outLines)
        outFile.close()
        #print "Wrote image %d."%structureIndex
        outFile = open(outputPath+".xrna", "w")
        outFile.write(xrnaOutLines)
        outFile.close()
        #print "Wrote xrna file %d."%structureIndex
        structureIndex += 1
    #print "Done."

#---------------------------------------------------
# Main execution
def main(args):
    try:
        #ctRead = args.ct.read().replace("\r","")
        #ctLines = ctRead.split("\n")
        ctLines = args.ct.readlines()
        dbnLines = ct_to_dbn(ctLines, args.structures)
        print dbnLines
        lineGroups = parseDotBracketFile(dbnLines)
        structureIndex = 1
        for lineGroup in lineGroups:
            title, seq, struct, pairedNuc = parseDotBracket(lineGroup)
            startNuc = 1
            if not args.no_title:
                fullTitle = title
                if args.title != "":
                    fullTitle = args.title
            else:
                fullTitle = " "
            print "Asking for structure %d . . ."%structureIndex
            #print "args.out: "+args.out
            epsResponse = getEPS(fullTitle, startNuc, seq, struct)
            if DEBUG==True:
                debugEPS = open("raw_response.eps","w")
                debugEPS.write(epsResponse)
                debugEPS.close()
            print "Coloring image %d . . ."%structureIndex
            colorStrings = makeColorStrings(seq, startNuc, args)
            print len(colorStrings)
            if DEBUG==True:
                print "colorStrings: "+str(colorStrings)
            outLines, xrnaOutLines = modifyEPS(epsResponse, seq, startNuc, colorStrings, fullTitle, pairedNuc)
            outputPath = os.path.normpath(args.out + str(structureIndex))
            print "outputPath: "+outputPath
            outFile = open(outputPath+".eps", "w")
            outFile.write(outLines)
            outFile.close()
            print "Wrote image %d."%structureIndex
            outFile = open(outputPath+".xrna", "w")
            outFile.write(xrnaOutLines)
            outFile.close()
            print "Wrote xrna file %d."%structureIndex
            structureIndex += 1
        print "Done."
    except Exception as e:
        sys.stderr.write(str(e.args))
        if "Pseudoviewer" in str(e.args[0]):
            sys.stderr.write(str(e.args[0])+"\n")
        else:
            sys.stderr.write("Error:"+traceback.format_exc())
        exit(1)

if __name__ == '__main__':
    args = parseStructArgs()
    main(args)
#---------------------------------------------------
