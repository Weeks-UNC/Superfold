##################################################################################
# GPL statement:
# This file is part of SuperFold.
#
# SuperFold is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SuperFold is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SuperFold.  If not, see <http://www.gnu.org/licenses/>.
#
# 17 Nov 2014
# 1.0 build
#
# Copywrite 2014 
# Greggory M Rice
# all rights reserved
##################################################################################
# Usage: predicted CTFile, correct CTFile, Shannon file, SHAPE File
#  this script generates postscript code from RNAstructure's Circlecompare
#  but with extra plotting features, such as annotations around nucleotide numbering
#
#  The CircleCompare program is available at: http://rna.urmc.rochester.edu/rnastructure.html
#
#  This derivative work was written by Gregg Rice - gmr@unc.edu
#  as a part of the Weeks lab at UNC - Chapel Hill
#
#  the codes from this file are available under a creative commons 4.0 share-alike
#  available here:   http://creativecommons.org/licenses/by-sa/4.0/
#
#  Date: 3 August 2014
#
#  Version: 1.0
#
#  Currently prints to command line have to pipe it into a .ps file
#  SHAPE/Shannon file needs to be length of RNA


import sys
#from RNAtools import CT
import math

# CT class from RNAtools
class CT:
    def __init__(self,fIN=None):
        #if givin an input file construct the ct object automatically
        if fIN:
            self.name = fIN
            self.num,self.seq,self.ct = self.readCT(fIN)
    
    def __str__(self):
        a = '{ Name= %s, len(CT)= %s }' % (self.name, str(len(self.ct)))
        return a
    
    def readCT(self,z):
        #reads a ct file, !reqires header!
        num,seq,bp = [],[],[]
        linesToRead = int(open(z).readlines()[0].rstrip().split()[0])
        #print linesToRead
        for i in open(z).readlines()[1:linesToRead+1]:
            a = i.rstrip().split()
            num.append(int(a[0])),seq.append(str(a[1])),bp.append(int(a[4]))
        return num,seq,bp
    
    def writeCT(self,fOUT):
        #writes a ct file
        w = open(fOUT,'w')
        line = '{0:6d} {1}\n'.format(len(self.num),self.name)
        for i in range(len(self.num)-1):
        #line+='{0:5d} {1} {2:5d} {3:5d} {4:5d} {5:5d} {0:5d}\n'\
        #.format(self.num[i],self.seq[i],self.num[i]-1,self.num[i]+1,self.ct[i])
            line += '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d}\n'.format(self.num[i],self.seq[i],self.num[i]-1,self.num[i]+1,self.ct[i])
        
        #last line is different
        i = len(self.num)-1
        line += '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d}\n'.format(self.num[i],self.seq[i],self.num[i]-1,0,self.ct[i])
        w.write(line)
        w.close()
    
    def copy(self):
        """
        returns a copy of the ct object
        """
        out = CT()
        out.name = self.name[:]
        out.num = self.num[:]
        out.seq = self.seq[:]
        out.ct = self.ct[:]
        return out
    
    def pairList(self):
        #returns a list of base pairs i<j
        out = []
        for nt in range(len(self.ct)):
            if self.ct[nt] != 0 and self.ct[nt] > self.num[nt]:
                out.append((self.num[nt],self.ct[nt]))
        return out
    
    def pair2CT(self, pairs, seq, name=None, skipConflicting=False):
        """
        pairs are an array of bp tuples ( i < j )
           i.e. [(4,26),(5,25),...]
        length is implied from the given sequence
        """
        length = len(seq)
        
        self.num = range(1,length+1)
        self.seq=seq
        
        #give it a name if it has one
        if name:
            self.name=name
        else:
            self.name='RNA_'+str(length)
        
        self.ct = []
        for i in range(length):self.ct.append(0)
        
        for i,j in pairs:
            if self.ct[i-1]!=0:
                print 'Warning: conflicting pairs, (%s - %s) : (%s - %s)' % (str(i),str(j),str(self.ct[i-1]),str(i))
                if skipConflicting: continue
            if self.ct[j-1]!=0:
                print 'Warning: conflicting pairs, (%s - %s) : (%s - %s)' % (str(i),str(j),str(j),str(self.ct[j-1]))
                if skipConflicting: continue
            self.ct[i-1]=j
            self.ct[j-1]=i
    
    def cutCT(self,start,end):
        """
        inclusively cuts a ct file and removes pairs outside the cutsite
        """
        start = int(start)
        end = int(end)
        
        out = CT()
        out.seq = self.seq[start-1:end]
        out.num = range(1,end-start+2)
        out.name = self.name + '_cut_'+str(start)+'_'+str(end)
        
        out.ct = []
        temp  = self.ct[start-1:end]
        # renumber from 1
        for nt in temp:
            nt_out = nt-(start-1)
            # cut out pairings that lie outside the window
            if nt_out <= 0 or nt_out > end-(start-1):
                nt_out = 0
            out.ct.append(nt_out)
        
        return out
    
    def contactDistance(self,i,j):
        """
        calculates the contact distance between pairs i,j in
        in the RNA using the RNAtools CT object. Method for
        calculating the contact is described in Hajdin et al
        (2013).
        """    
        #correct for indexing offset
        i, j = i-1, j-1
        
        #error out if nucleotide out of range
        if max(i,j) > len(self.ct):
            print 'Error!, nucleotide {0} out of range!'.format(max(i,j)+1)
            return
        
        #i must always be less than j, correct for this
        if j < i:
            i, j = j, i
        
        
        count = 0.0
        tempBack = 10000.0
        k = int(i)
        
        def backTrace(rna,j,k):
            bcount = 2
            k -= 1
            if k-1 == j:
                return bcount
            bcount += 1
            while k > j:
                if rna.ct[k] == 0:
                    k -= 1
                    bcount += 1
                # if not single stranded, exit the backtrace
                else:
                    return None
            return bcount
        # search forward through sequence
        while k < j:
            #debuging stuff, prevent infinite loop
            #if count > 200:
            #    print i,j
            #    break
            
            #nonpaired nucleotides are single beads
            if self.ct[k] == 0:
                k += 1
                #count += 6.5
                count += 1
            
            #branches through helices can't be skipped
            #treated as single beads
            elif self.ct[k] > j + 1:
                # try backtracing a few if it is close (within 10)
                if self.ct[k] - 10 < j:
                    back = backTrace(self, j,self.ct[k])
                    # if the backtracing is able to reach
                    # your nt, add its length to the count
                    # and break
                    if back:
                        #print 'backitup'
                        # store the backtracing count for later
                        # if it ends up being lower than the final
                        # we will use it instead
                        if count + back < tempBack:
                            tempBack = count + back
                k += 1
                #count += 6.5
                count += 1
            
            # simple stepping
            elif self.ct[k] < i + 1:
                k += 1
                #count += 6.5
                count += 1
            
            #handle branching, jumping the base of a helix is only 1
            else:
                # one count for the step to the next
                count += 1
                k = self.ct[k]
                
                # if by jumping you land on your nt
                # stop counting
                if k - 1 == j:
                    break
                
                # one count for jumping the helix
                #count += 15.0
                count += 1
            #print i,k,j
        finalCount = min(count, tempBack)
        return finalCount
    
    def extractHelicies(self,fillPairs=True):
        """
        returns a list of helicies in a given CT file and the nucleotides
        making up the list as a dict object. e.g. {0:[(1,50),(2,49),(3,48)}
        
        defaults to filling 1,1 and 2,2 mismatches
        """
        # first step is to find all the helicies
        rna = self.copy()
        if fillPairs:
            rna = rna.fillPairs()
        helicies = {}
        nt = 0
        heNum = 0
        
        while nt < len(rna.ct):
            if rna.ct[nt] == 0:
                nt += 1
                continue
            else:
                # skip dups
                if nt > rna.ct[nt]:
                    nt += 1
                    continue
                
                tempPairs = []
                stillHelix = True
                
                previous = rna.ct[nt]
                
                while stillHelix:
                    # see if the helix juts up against another one
                    if abs(rna.ct[nt] - previous ) > 1:
                        break
                    #add the pairs
                    elif rna.ct[nt] != 0:
                        tempPairs.append((nt+1,rna.ct[nt]))
                        previous = rna.ct[nt]
                        nt+=1
                    #handle slip case
                    else:
                        if rna.ct[nt+1] == rna.ct[nt-1]-1:
                            nt+=1
                            #print 'slip'
                        else:
                            break
                
                # remove single bp helicies
                if len(tempPairs) <= 1:
                    continue
                helicies[heNum] = tempPairs
                heNum += 1
        return helicies
    
    def fillPairs(self):
        """
        fills 1,1 and 2,2 mismatches in an RNA structure
        """
        rna = self.copy()
        # fill in 1,1 mismatch, 2,2 mismatch
        for i in xrange(len(rna.ct)-3):
            if rna.ct[i+1] == 0:
                if rna.ct[i] - rna.ct[i+2] == 2:
                    rna.ct[i+1] = rna.ct[i] - 1
            if rna.ct[i+1] + rna.ct[i+2] == 0:
                if rna.ct[i] - rna.ct[i+3] == 3:
                    rna.ct[i+1] = rna.ct[i] - 1
                    rna.ct[i+2] = rna.ct[i] - 2
        return rna
    
    
    def extractPK(self,fillPairs=True):
        """
        returns the pk1 and pk2 pairs from a CT file. Ignores single base
        pairs. PK1 is defined as the helix crossing the most other helices.
        if there is a tie, the most 5' helix called pk1
        
        returns pk1 and pk2 as a list of base pairs e.g [(1,10),(2,9)...
        """
        
        def checkOverlap(h1,h2):
            # only need to check one set of pairs from each
            # of the helicies. Test is to see if they form
            # a cross hatching pattern
            if max(h1[0]) > min(h2[0]) > min(h1[0]):
                if max(h2[0]) > max(h1[0]):return True
            if max(h2[0]) > min(h1[0]) > min(h2[0]):
                if max(h1[0]) > max(h2[0]):return True
            else: return False
        
        
        # make a copy so we don't destroy the original object
        rna = self.copy()
        
        
        #self.writeCT('foo.ct')
        #rna.writeCT('bar.ct')
        
        # get the helicies by calling the extract helix function
        helicies = rna.extractHelicies(fillPairs=fillPairs)
        heNum = len(helicies)

        
        # do the helicies overlap? Check for a crosshatching pattern
        # append them to a list if they have it.
        overlaps = [] # stores the helix number
        
        for i in xrange(heNum):
            for j in xrange(i+1,heNum):
                if checkOverlap(helicies[i],helicies[j]):
                    overlaps.append((i,j))
        
        #print overlaps
        #print '#'*30
        # if there are no overlapping bps, return none
        if len(overlaps) == 0:
            return None, None
        
        
        # determine which is pk1
        allHelix = []
        for i,j in overlaps:
            allHelix.append(i), allHelix.append(j)
        pk1Helix = max(set(allHelix), key=allHelix.count)
        pk2Helix = filter(lambda x: x != pk1Helix, allHelix)
        
        # construct list of base pairs
        pk1 = helicies[pk1Helix]
        pk2 = []
        for i in pk2Helix:
            for j in helicies[i]:
                pk2.append(j)
        
        return pk1, pk2

#below is the stuff for PyCircleCompare

def compareRNA(predRNA, correctRNA, allowSlip=True):
    """
    Compares a predicted structure to a correct structure. Returns
    correct, missing, and extra base pairs as arrays of tuples [(1,10),(2,8)...]
    the flag allowSlip will count pairs that are within one nucleotide of the correct
    pairing as correct
    
    returns correct, missing, extra, int(numCorrect)
    """
    
    # get the base pairs from objects
    correctBPs = correctRNA.pairList()
    predBPs = predRNA.pairList()
    
    correct, missing, extra = [], [], []
    numCorrect = 0
    
    # calculate the num of correct bps, and extra base pairs
    for pair in predBPs:
        if pair in correctBPs:
            correct.append(pair)
            numCorrect += 1
        #slip case
        elif allowSlip:
            if (pair[0],pair[1]+1) in correctBPs:
                numCorrect += 1
            elif (pair[0],pair[1]-1) in correctBPs:
                numCorrect += 1
            elif (pair[0]+1,pair[1]) in correctBPs:
                numCorrect += 1
            elif (pair[0]-1,pair[1]) in correctBPs:
                numCorrect += 1
            extra.append(pair)
        else:
            extra.append(pair)
    
    # determine the missing base pairs
    for pair in correctBPs:
        if pair not in predBPs:
            missing.append(pair)
    
    return correct, missing, extra, numCorrect

def genSHAPEcolorString(shapeData, diffColor=False):
    """
    generates the colors for the SHAPE annotations in the Circle plot
    """
    
    line = ''
    for i in shapeData:
        # if the shape reactivities are differential convert them
        # back to a shapelike scale assuming standard sm and si
        if diffColor:
            if i > -500:
                i = math.e**( (math.log(i+1.0)-0.4)/1.8 ) - 1.0
        
        if i >= 0.85:
            line += "[1.00 0.00 0.00]\n"
        elif i >= 0.4:
            line += "[1.00 0.50 0.00]\n"
        elif i >= -500:
            line += "[0.00 0.00 0.00]\n"
        else:
            line += "[0.67 0.67 0.67]\n"
    return line

def genShannoncolorString(shannonData):
    """
    generates the colors for the shannon annotations in the Circle plot
    """
    # alpha and beta describing the SE distribution: alpha0.2,beta=2.63
    # the cutoffs are the ppf of the beta distribution at 0.85, 0.91, and 0.94
    # new cutoffs - 6 aug 2013
    s1, s2, s3 = 0.186,  0.277 ,  0.33207707
    line = ''
    for i in shannonData:
        if   i > s3:
            line += "[0.00 0.00 0.00]\n"
        elif i >= s2:
            line += "[0.00 0.00 1.00]\n"
        elif i >= s1:
            line += "[0.00 0.67 1.00]\n"
        else:
            line += "[1.00 1.00 1.00]\n"
    return line
    
def genPairingString(correct, missing, extra, slipped, length):
    """
    generates the pairing strings and determines whether the pair is one of the
    slipped base pairs. Also colors based on their correctness
    """
    line = ''
    count = 0
    for nt in range(length):
        
        for pair in correct:
            if pair[0] == nt+1:
                line += '[{0} {1} 0.40 0.40 0.40 {2}]\n'.format(pair[0],pair[1],int(pair in slipped))
                count += 1
        for pair in missing:
            if pair[0] == nt+1:
                line += '[{0} {1} 1.00 0.00 0.00 {2}]\n'.format(pair[0],pair[1],int(pair in slipped))
                count += 1
        for pair in extra:
            if pair[0] == nt+1:
                line += '[{0} {1} 0.40 0.00 0.60 {2}]\n'.format(pair[0],pair[1],int(pair in slipped))
                count += 1
    return line, count

def getSHAPE(fIN):
    """
    reads a .SHAPE formatted file and returns an array of the values
    """
    out = []
    for line in open(fIN).readlines():
        i = line.rstrip().split()
        out.append(float(i[1]))
    return out

def getCorrelData(fIN):
    out = {'i':[],'j':[],'correl':[]}
    for line in open(fIN).readlines()[1:]:
        k = line.rstrip().split()
        i = int(k[0])
        j = int(k[1])
        correl = float(k[2])
        out['i'].append(i), out['j'].append(j), out['correl'].append(correl)
    return out

def getScaleFactor(length):
    """
    calculates teh scaling factor. Was determined heuristically from a range or RNA sizes in CircleCompair
    """
    scaleFactor = 1.0
    if length > 76:
        scaleFactor = 74.0875 / float(length) + 0.020367
    else:
        scaleFactor = -1.06*math.log10(length) + 5.44
    
    if scaleFactor > 1.0:
        scaleFactor = 1.0
    
    return scaleFactor

def zeros_like(arr):
    """
    returns an array of zeros the length of the input
    """
    out = []
    for i in range(len(arr)):
        out.append(0.00)
    return out


def genCorrelString(correlDat):
    """
    correlDat is a dict. i, j, correl are names.
    """
    line = ""
    count = 0
    length = len(correlDat['i'])
    for num in xrange(length):
        if correlDat['correl'][num] > 0.035:
            line += '[{0} {1} 0.00 0.50 0.00 {2}]\n'.format(correlDat['i'][num],correlDat['j'][num],0)
            count += 1
        elif correlDat['correl'][num] > 0.025:
            line += '[{0} {1} 0.80 0.80 0.10 {2}]\n'.format(correlDat['i'][num],correlDat['j'][num],0)
            count += 1
    return line, count

def makeCircle(ct_pred, ct_correct, shannon, shape, correlDat, slipped, diffColor=False, offset=1):
    """
    function to make a circleplot from two structures, shannon and shape reactivity and slip status.
    diffColor uses the differential coloring scale based on a slope and intercept of 1 and -1. offset
    delineates nt start position.
    
    returns postscript file lines as a string
    """

    correct, missing, extra, bpCorrect = compareRNA(ct_pred, ct_correct)
    
    
    shannonColoring = genShannoncolorString(shannon)
    shapeColoring = genSHAPEcolorString(shape, diffColor)
    
    pairingTable, numPairs = genPairingString(correct,missing,extra,slipped,len(ct_pred.ct))
    
    pairingTableCorrel, numPairsCorrel = genCorrelString(correlDat)
    
    # add the new correlations to the pairing array
    pairingTable += pairingTableCorrel
    numPairs += numPairsCorrel
    #print pairingTable
    scaleFactor = getScaleFactor(len(ct_correct.ct))
    
    # heuristically based on circlecompare
    cirRadius = 3*len(ct_correct.ct)+14
    
    sens = bpCorrect/float(len(correct)+len(missing))
    ppv = bpCorrect/float(len(correct)+len(extra))


    Circle = """%!

% Set font size and type.
/fontSize 24 def
/quarterFont fontSize 4 div def
/halfFont fontSize 2 div def
/Courier findfont fontSize scalefont setfont
% Set variables handling number, placement of nucleotides.
/numBases {0} def
/currentBase 0 def
/basePoints numBases array def

% Set variables handling scaling, translation of circular backbone.
/scaleFactor {1} def
/scaleFactorX scaleFactor def
/scaleFactorY scaleFactor def
/translateFactorX 0 def
/translateFactorY 0 def

% Set variables handling properties of circular backbone.
/radius {2} def
/center 306 scaleFactor div def
/angle 360 numBases 2 add div -1 mul def
/labelSpace 40 def

% Create the array of nucleotides.
/bases [
""".format(len(ct_correct.ct),scaleFactor,cirRadius)

    Circle += ' '.join(map('({0})'.format, ct_correct.seq)) # lists sequence as '(A) (G)...'
    Circle +="""
] def

% Write the array of pairings.
/pairings [ """
    Circle += pairingTable
    Circle += """
] def

% Set variables handling number, placement of pairings.
/numPairings {0} def
/numPseudoknotted 0 def
/currentPairing 0 def

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write out the combined circular structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the font size and type for the descriptor.
/Courier-Bold findfont fontSize 2 div scalefont setfont

/numDescriptors 1 def
/descriptors [
    [(Predicted Structure file name: {1}) 0.00 0.00 0.00 40 740]
] def

% Write the contents of the descriptor array.
descriptors {{ aload pop moveto setrgbcolor show }} forall

% Write the array of annotation colors.
/annotations [ """.format(numPairs, ct_pred.name, ct_correct.name)
    Circle += shapeColoring
    Circle += """] def

% Create the legend array.
/legend [
    [(0.035 >= Chi-sq) 0.0 0.5 0.0]
    [(0.035 >= Chi-sq >0.025) 0.8 0.8 0.1]
] def

% Show the legend.
%/Courier-Bold findfont 8 scalefont setfont
%40 595 moveto

%legend {{ gsave aload pop setrgbcolor show grestore 0 -10 rmoveto }} forall

% Reset the font to its original size, type, and color.
/Courier findfont fontSize scalefont setfont
0 setgray

% Write the translation and scaling of the main image.
gsave
scaleFactorX scaleFactorY scale
translateFactorX translateFactorY translate

% Use repeat loop to write nucleotides on the circular path.
0 1 numBases {{
    clear

    % Move to appropriate point, angle to write next nucleotide.
    center center moveto
    currentBase angle mul rotate
    quarterFont radius rmoveto

    % Determine x and y coordinates of the nucleotide location.
    currentBase angle mul -1 mul rotate
    currentpoint
    /y exch cvi def
    /x exch cvi def
    currentBase angle mul rotate
    quarterFont -1 mul 0 rmoveto

    % Save the nucleotide location and show the nucleotide.
    /point [ x y ] def
    basePoints currentBase point put
    
    %define offset for nums
    /offset {2} def

    % Set variables for conditions where number labels are found.
    /numCond1 currentBase offset add 10 mod 0 eq def
    /numCond2 currentBase offset add 1 eq def
    /numCond3 currentBase offset add numBases eq def

    /annotation annotations currentBase get def
    annotation 0 get annotation 1 get annotation 2 get setrgbcolor

    % If a condition is met for a number label, write a label.
    numCond1 numCond2 or numCond3 or {{
        /Courier-Bold findfont fontSize scalefont setfont
        bases currentBase get show
        halfFont -1 mul labelSpace rmoveto
        /numString 7 string def
        currentBase offset add numString cvs show
        0 labelSpace -1 mul rmoveto
        /Courier findfont fontSize scalefont setfont
    }}
    {{ bases currentBase get show }} ifelse
    0 setgray

    % Return to circle center, rotate to 0, increment base.
    0 radius -1 mul rmoveto
    currentBase angle mul -1 mul rotate
    /currentBase currentBase 1 add def
}} repeat

% Write the pairing lines inside the backbone.
0 setgray
0 1 numPairings {{
    /pair pairings currentPairing get def

    % Determine coordinates of first nucleotide in pair.
    /base1 pair 0 get 1 sub def
    /point1 basePoints base1 get def
    /x1 point1 0 get def
    /y1 point1 1 get def

    % Determine coordinates of second nucleotide in pair.
    /base2 pair 1 get 1 sub def
    /point2 basePoints base2 get def
    /x2 point2 0 get def
    /y2 point2 1 get def

    % If the pair should be colored, set its appropriate color.
    pair 2 get pair 3 get pair 4 get setrgbcolor
    
    % if the pair is slipped set variable so
    /isSlipped pair 5 get def

    % Draw current pair, then increment current pairing.
    /between base2 base1 sub def
    between numBases 2 div gt {{ /between numBases between sub def }} if

    /midX x1 x2 add 2 div def
    /midY y1 y2 add 2 div def

    /gamma 0.9 def
    /centerThresh radius 2 mul 8 div 5 mul def

    /ends x2 x1 sub x2 x1 sub mul y2 y1 sub y2 y1 sub mul add sqrt def
    /distance between 2 mul numBases div radius mul gamma mul def
    /lineAngle center midY sub center midX sub atan def
    /distX lineAngle cos distance mul 2 mul def
    /distY lineAngle sin distance mul 2 mul def

    /controlX midX distX add def
    /controlY midY distY add def

    ends centerThresh ge {{
        /controlX center def
        /controlY center def
    }} if

    isSlipped 0.5 gt {{
    x1 y1 moveto x1 y1 controlX controlY x2 y2 curveto [15 5] 0 setdash stroke
    }}{{
    x1 y1 moveto x1 y1 controlX controlY x2 y2 curveto [100000 1] 0 setdash stroke
    }} ifelse
    0 setgray
    /currentPairing currentPairing 1 add def
}} repeat


 % Set variables handling number, placement of nucleotides.
 /numBasesf {1}""".format(numPairs,len(ct_correct.ct),offset)
    Circle += """ def
 /currentBasef 0 def
 /basePointsf numBasesf array def
 
 
 % Set variables handling properties of circular backbone.
 /radiusf {0} def
 /centerf 306 scaleFactor div def
 /anglef 359.75 numBases 2 add div -1 mul def
 /labelSpace 20 def
 
 /annotationsf [ """.format(int(cirRadius+10)) + shannonColoring + """
 ] def
 
 
 /basesf [ """
 
    Circle += ' '.join(map('(*)'.format, ct_correct.seq))
    Circle += """

 ] def
 
 
 
 
 
 /fontSizef 40 def
 /quarterFontf fontSizef 4 div def
 
 
 % Use repeat loop to write nucleotides on the circular path.
 0 1 numBasesf {{
     clear
 
     % Move to appropriate point, angle to write next nucleotide.
     center centerf moveto
     currentBasef anglef mul rotate
     quarterFontf radiusf rmoveto
 
     % Determine x and y coordinates of the nucleotide location.
     currentBasef anglef mul -1 mul rotate
     currentpoint
     /xf exch cvi def
     /yf exch cvi def
     currentBasef angle mul rotate
     quarterFontf -1 mul 0 rmoveto
 
     % Save the nucleotide location and show the nucleotide.
     /point [ xf yf ] def
     basePointsf currentBasef point put
     /offset {5} def
 
     % Set variables for conditions where number labels are found.
     /numCond1 currentBasef offset add 10 mod 0 eq def
     /numCond2 currentBasef offset add 1 eq def
     /numCond3 currentBasef offset add numBases eq def
 
     /annotationf annotationsf currentBasef get def
     annotationf 0 get annotationf 1 get annotationf 2 get setrgbcolor
 
     % If a condition is met for a number label, write a label.
     numCond1 numCond2 or numCond3 or {{
         /Courier-Bold findfont fontSizef scalefont setfont
         basesf currentBasef get show
         /Courier-Bold findfont fontSize scalefont setfont
         halfFont -1 mul labelSpace rmoveto
         /numString 7 string def
        % currentBasef offset add numString cvs show
         0 labelSpace -1 mul rmoveto
         /Courier findfont fontSizef scalefont setfont
     }}
     {{ basesf currentBasef get show }} ifelse
     0 setgray
 
     % Return to circle center, rotate to 0, increment base.
     0 radiusf -1 mul rmoveto
     currentBasef angle mul -1 mul rotate
     /currentBasef currentBasef 1 add def
 }} repeat
 
  % Write the array of annotation colors.
 

grestore

/Courier-Bold findfont fontSize 2 div scalefont setfont

/x1 570 (Accepted:) stringwidth pop sub def
/x2 570 (Pairs: {2}) stringwidth pop sub def
/x3 570 (Pseudoknotted Pairs: 0) stringwidth pop sub def
/x4 570 (Sensitivity: {0} / {2} = {3:.2%}) stringwidth pop sub def
/x5 570 (PPV: {0} / {1} = {4:.2%}) stringwidth pop sub def

/statistics [
    [(Predicted:) 40 70]
    [(Pairs: {1}) 40 55]
    [(Pseudoknotted Pairs: 0) 40 40]
    [(Accepted:) x1 70]
    [(Pairs: {2}) x2 55]
    [(Pseudoknotted Pairs: 0) x3 40]
    [(Sensitivity: {0} / {2} = {3:.2%}) x4 595]
    [(PPV: {0} / {1} = {4:.2%}) x5 580]
] def

%statistics {{ aload pop moveto show }} forall
showpage
""".format(len(correct), len(correct)+len(extra), len(correct)+len(missing), sens, ppv, offset)

    return Circle


if __name__ == '__main__':
    if len(sys.argv) <=3:
        print 'Usage: CircleDMS.py structure.ct correlData.txt output.ps'
        sys.exit()
    ct_pred = CT(sys.argv[1])
    
    correlDat = getCorrelData(sys.argv[2])
    # slipped base pairs are an array of tuples of the slipped pairs
    # slipped = [(1,78),(2,77)]
    slipped = []
    shannon = zeros_like(ct_pred.ct)
    shape = zeros_like(ct_pred.ct)
    #shannon = getSHAPE(sys.argv[3])
    #shape = getSHAPE(sys.argv[4])

    lines = makeCircle(ct_pred, ct_pred, shannon, shape, correlDat, slipped, offset=1)
    w = open(sys.argv[3],'w')
    w.write(lines)
    w.close()
