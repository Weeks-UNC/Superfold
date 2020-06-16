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

# 17 Nov 2014
# Copywrite 2014
# Greggory M Rice
# all rights reserved
# 1.0 build
##################################################################################
# render one or two CT secondary structures with helices shown as semicircular arcs
# This simple version of the script does not plot any additional data above the structure(s).
# Written by Steve Busan, modified in 2014 with permission for python hooks Gregg Rice
# ------------------------------------------------------------------------------------------

import sys, os, traceback, copy, re, math

def parseCTFile(lines):
    numberNucs = int(lines[0].strip().split()[0])
    seq = list("_"*numberNucs)
    pairedNuc = [0]*numberNucs
    for n in range(0,numberNucs):
        lineIndex = n+1
        splitLine = lines[lineIndex].strip().split()
        seq[n] = splitLine[1]
        pairedNuc[n] = int(splitLine[4])
        
    return seq, pairedNuc

def plotArcRibbons(pairedNuc, ax, color, alpha=1.0):
    from matplotlib.path import Path
    import matplotlib.patches as patches

    handleLengthFactor = 4*(math.sqrt(2)-1)/3
    #vert = 100.0

    i = 0
    while i < len(pairedNuc):
        if pairedNuc[i] > i+1:
            outerPair = [i+0.5,pairedNuc[i]+0.5]
            # find the right side of helix
            lastPairedNuc = pairedNuc[i]
            offset = 1
            while i+offset<len(pairedNuc) and abs(pairedNuc[i+offset]-lastPairedNuc) == 1:
                lastPairedNuc = pairedNuc[i+offset]
                offset += 1
            innerPair = [i+offset+0.5, pairedNuc[i+offset-1]-0.5] 
            i += offset-1
            outerRadius = (outerPair[1]-outerPair[0])/2.0
            innerRadius = (innerPair[1]-innerPair[0])/2.0
            #print "innerPair %s, outerPair %s"%(str(innerPair),str(outerPair))
            verts = [
            (outerPair[0], 0), # outer left

            (outerPair[0], -handleLengthFactor*outerRadius), # outer left control 1
            (outerPair[0]+outerRadius-handleLengthFactor*outerRadius, -outerRadius), # outer left control 2
            (outerPair[0]+outerRadius, -outerRadius), # outer center
            (outerPair[0]+outerRadius+handleLengthFactor*outerRadius, -outerRadius), # outer right control 1
            (outerPair[1], -handleLengthFactor*outerRadius), # outer right control 2
            
            (outerPair[1], 0), # outer right
            (innerPair[1], 0), # inner right

            (innerPair[1], -handleLengthFactor*innerRadius), # inner right control 1
            (innerPair[0]+innerRadius+handleLengthFactor*innerRadius, -innerRadius), # inner right control 2
            (innerPair[0]+innerRadius, -innerRadius), # inner center
            (innerPair[0]+innerRadius-handleLengthFactor*innerRadius, -innerRadius), # inner right control 1
            (innerPair[0], -handleLengthFactor*innerRadius), # inner right control 2

            (innerPair[0], 0), # inner left
            (outerPair[0], 0) # outer left duplicate point
            ]

            for n in xrange(len(verts)):
                verts[n] = [verts[n][0], verts[n][1]-1.2]#/vert-0.013]
            codes = [
            Path.MOVETO,

            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,

            Path.LINETO,
            Path.LINETO,
            Path.LINETO,

            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,

            Path.LINETO,
            Path.LINETO,
            Path.CLOSEPOLY,
            ]
            path = Path(verts, codes)
            patch = patches.PathPatch(path, facecolor=color, linewidth=0, edgecolor='none', alpha=alpha)
            ax.add_patch(patch)
            patch.set_clip_on(False)
        i += 1

def writePlot(outPath="arcs.pdf",title="",seq=["A"],pairedNucA=[],pairedNucB=[],alpha=0.5):
    import matplotlib as mp
    mp.use('Agg')
    mp.rcParams['xtick.major.size'] = 8
    mp.rcParams['xtick.major.width'] = 2.5
    mp.rcParams['xtick.direction'] = 'out'
    mp.rcParams['xtick.minor.size'] = 4 
    mp.rcParams['xtick.minor.width'] = 1


    import matplotlib.pyplot as plot
    import matplotlib.patches as patches
    import matplotlib.gridspec as gridspec

    num = range(1,len(seq)+1)
    scaleFactor = 0.05
    # find longest base-pair and scale height of plot to fit this arc
    maxDistance = 0
    for i in range(len(pairedNucA)):
        fromNuc = i+1
        toNuc = pairedNucA[i]    
        if toNuc == 0:
            toNuc = fromNuc
        dist = toNuc-fromNuc
        if dist > maxDistance:
            maxDistance = dist
    if pairedNucB != []:
        for i in range(len(pairedNucA)):
            fromNuc = i+1
            toNuc = pairedNucA[i]    
            if toNuc == 0:
                toNuc = fromNuc
            dist = toNuc-fromNuc
            if dist > maxDistance:
                maxDistance = dist
    figWidth = len(seq)*scaleFactor
    figHeight = maxDistance/2.0*scaleFactor
    #print "seq len: %f, maxDistance: %f, figwidth: %f, figheight: %f"%(len(seq), maxDistance, figWidth, figHeight)
    fig = plot.figure(figsize=(figWidth, figHeight)) # 500*scaleFactor
    
    #ax1 = plot.subplot(211)
    #ax1.set_frame_on(False)
    #fig.text(0.0, 0.25, title, horizontalalignment='left',size="24",weight="bold")
    plot.xlim(0,len(seq))
    #plot.yticks([threshold,-threshold],[threshold,-threshold],fontsize=9)
    #plot.yticks(fontsize=20)

    ax2 = plot.gca()

    #ax1.axes.get_xaxis().
    ax2.get_xaxis().tick_top()  
    #ax2.xaxis.labelpad = 20
    #ax2.tick_params(axis='x', direction='out')
    #ax2.tick_params(which='major', fontproperties = {'fontsize':14,'weight':'bold','rotation':30})

        
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    if len(seq) <= 500:
        majorLocator = MultipleLocator(100)
        minorLocator = MultipleLocator(10)
        interval = 10
    else:
        majorLocator = MultipleLocator(500)
        minorLocator = MultipleLocator(100)
        interval = 5
    majorFormatter = FormatStrFormatter('%i')
    minorFormatter = FormatStrFormatter('%i')

    ax2.xaxis.set_major_locator(majorLocator)
    ax2.xaxis.set_major_formatter(majorFormatter)
    ax2.xaxis.set_minor_locator(minorLocator)
    ax2.xaxis.set_minor_formatter(minorFormatter)
    
    #ax1.xaxis.set_visible(False)

    plot.subplots_adjust(hspace=0)
    plot.ylim((-float(maxDistance)/2.0,0))
    #plot.ylim((-1.0,0))

    ax2.set_frame_on(False)
    ax2.axes.get_yaxis().set_visible(False)
    #ax2.axes.get_xaxis().tick_bottom()
    xlabels = ax2.axes.get_xaxis().get_majorticklabels()
    xlabels[0].set_visible(False)
    for label in xlabels:
        label.set_weight('bold')
        label.set_size(14)
        #label.set_rotation(30)
    xlabels = ax2.axes.get_xaxis().get_minorticklabels()
    xlabels[0].set_visible(False)
    labelCount = 0
    for label in xlabels:
        label.set_size(7)
        #label.set_rotation(30)
        # also need to hide minor tick labels that overlap major tick labels
        if labelCount%interval==0:
            label.set_visible(False)
        labelCount += 1

    xticks = ax2.axes.get_xaxis().get_major_ticks()
    xticks[0].set_visible(False)
    xticks = ax2.axes.get_xaxis().get_minor_ticks()
    xticks[0].set_visible(False)

    # write title (structure names)
    #plot.text(0.0,-maxDistance/2.0, title, horizontalalignment='left',size="30",weight="bold", color="blue")

    bothColor = "green"
    aColor = "red"
    bColor = "purple"

    if pairedNucB == []:
        bothColor = "black"
        pairedNucB = pairedNucA        

    # split paired nucs into both, AnotB, BnotA
    bothPaired = [0]*len(pairedNucA)
    aOnlyPaired = [0]*len(pairedNucA)
    bOnlyPaired = [0]*len(pairedNucA)
    for i in xrange(len(pairedNucA)):
        a = pairedNucA[i]
        b = pairedNucB[i]
        if a == b:
            if a > i+1:
                bothPaired[i] = a
        else:
            if a > i+1:
                aOnlyPaired[i] = a
            if b > i+1:
                bOnlyPaired[i] = b
        
    print bothPaired
    plotArcRibbons(aOnlyPaired, ax2, aColor, alpha=alpha)
    plotArcRibbons(bOnlyPaired, ax2, bColor, alpha=alpha)
    plotArcRibbons(bothPaired, ax2, bothColor, alpha=alpha)
    
    xmax, xmin, ymin, ymax = plot.axis()

    # put nuc sequence on axis 
    if len(seq) <= 500:                                                                            
        fontProp = mp.font_manager.FontProperties(family = "monospace", 
                                              style="normal",
                                              weight="extra bold",
                                              size="4")
        for i in range(len(seq)):
            nuc = seq[i]
            if nuc == "T":
                nuc = "U"
            col = "black"
            plot.annotate(nuc, xy=(i+0.5,ymax),fontproperties=fontProp,color=col,annotation_clip=False,verticalalignment="top")

    plot.savefig(outPath,dpi=100,bbox_inches="tight")
    
def writeManyPlot(outPath="arcs.pdf",title="",seq=["A"],pairedNucArr=[], arcColors = [],alpha=[1.0], maxDistance=None):
    """
    draws arcs for many sets of nucleotides, pairedNucArr is a 2D array containing all sets of plotting
    elements. arcColor is the same length as pairedNucArr but contains color strings
    """
    
    def findMaxDistance(pairingArrays):
        """
        finds the maximum pairing distance from a 2D array of all arc elements
        """
        maxDistance = 0
        
        for pairedNucA in pairingArrays:
            for i in range(len(pairedNucA)):
                fromNuc = i+1
                toNuc = pairedNucA[i]    
                if toNuc == 0:
                    toNuc = fromNuc
                dist = toNuc-fromNuc
                if dist > maxDistance:
                    maxDistance = dist
        
        return maxDistance
    
    import matplotlib as mp
    mp.use('Agg')
    mp.rcParams['xtick.major.size'] = 8
    mp.rcParams['xtick.major.width'] = 2.5
    mp.rcParams['xtick.direction'] = 'out'
    mp.rcParams['xtick.minor.size'] = 4 
    mp.rcParams['xtick.minor.width'] = 1


    import matplotlib.pyplot as plot
    import matplotlib.patches as patches
    import matplotlib.gridspec as gridspec

    num = range(1,len(seq)+1)
    scaleFactor = 0.05
    
    # find longest base-pair and scale height of plot to fit this arc
    if not maxDistance:
        maxDistance = findMaxDistance(pairedNucArr)
    
    # set up figure dimensions
    figWidth = len(seq)*scaleFactor
    figHeight = maxDistance/2.0*scaleFactor
    
    fig = plot.figure(figsize=(figWidth, figHeight)) # 500*scaleFactor
    
    plot.xlim(0,len(seq))

    ax2 = plot.gca()

    # ticks on top
    ax2.get_xaxis().tick_top()  

    # determine tick locations based on sequece length
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    if len(seq) <= 500:
        majorLocator = MultipleLocator(100)
        minorLocator = MultipleLocator(10)
        interval = 10
    else:
        majorLocator = MultipleLocator(500)
        minorLocator = MultipleLocator(100)
        interval = 5
    majorFormatter = FormatStrFormatter('%i')
    minorFormatter = FormatStrFormatter('%i')

    ax2.xaxis.set_major_locator(majorLocator)
    ax2.xaxis.set_major_formatter(majorFormatter)
    ax2.xaxis.set_minor_locator(minorLocator)
    ax2.xaxis.set_minor_formatter(minorFormatter)
    

    plot.subplots_adjust(hspace=0)
    plot.ylim((-float(maxDistance)/2.0,0))
    #plot.ylim((-1.0,0))

    ax2.set_frame_on(False)
    ax2.axes.get_yaxis().set_visible(False)
    #ax2.axes.get_xaxis().tick_bottom()
    xlabels = ax2.axes.get_xaxis().get_majorticklabels()
    xlabels[0].set_visible(False)
    for label in xlabels:
        label.set_weight('bold')
        label.set_size(14)
        #label.set_rotation(30)
    xlabels = ax2.axes.get_xaxis().get_minorticklabels()
    xlabels[0].set_visible(False)
    labelCount = 0
    for label in xlabels:
        label.set_size(7)
        #label.set_rotation(30)
        # also need to hide minor tick labels that overlap major tick labels
        if labelCount%interval==1:
            label.set_visible(False)
        labelCount += 1

    xticks = ax2.axes.get_xaxis().get_major_ticks()
    xticks[0].set_visible(False)
    xticks = ax2.axes.get_xaxis().get_minor_ticks()
    xticks[0].set_visible(False)


    bothColor = "green"
    aColor = "red"
    bColor = "purple"

    # plot the arcs
    for arc in range(len(pairedNucArr)):
        #if arc % 500 == 0:
        #    print arc, len(pairedNucArr)
        plotArcRibbons(pairedNucArr[arc], ax2, arcColors[arc], alpha=alpha[arc])
    
    #plotArcRibbons(aOnlyPaired, ax2, aColor, alpha=alpha)
    #plotArcRibbons(bOnlyPaired, ax2, bColor, alpha=alpha)
    #plotArcRibbons(bothPaired, ax2, bothColor, alpha=alpha)
    
    xmax, xmin, ymin, ymax = plot.axis()

    # put nuc sequence on axis 
    if len(seq) <= 500:                                                                            
        fontProp = mp.font_manager.FontProperties(family = "monospace", 
                                              style="normal",
                                              weight="extra bold",
                                              size="4")
        for i in range(len(seq)):
            nuc = seq[i]
            if nuc == "T":
                nuc = "U"
            col = "black"
            plot.annotate(nuc, xy=(i+0.5,ymax),fontproperties=fontProp,color=col,annotation_clip=False,verticalalignment="top")

    plot.savefig(outPath,dpi=100,bbox_inches="tight")

if __name__=="__main__":
    if len(sys.argv) < 3:
        print "\nUsage: python drawArcRibbons_simple.py <out_name.pdf> <accepted_structure.ct> [<different_structure.ct>] [<transparency 0.0-1.0>]"
        print "\n    pdf, eps, and png are all acceptable, but png will be large\n"
        exit()
    outPath = sys.argv[1]
    ctFileA = open(sys.argv[2],"rU")
    alpha=1.0
    alphaArg3 = False
    try:
        alpha = 1.0-float(sys.argv[3])
        alphaArg3 = True
    except (ValueError, IndexError) as e:
        try:
            alpha = 1.0-float(sys.argv[4])
        except IndexError:
            pass

 
    if alphaArg3 == True:
        haveSecondCT = False
    else:
        haveSecondCT = True
        try:
            ctFileB = open(sys.argv[3],"rU")
        except IndexError:
            haveSecondCT = False
 
    #title = os.path.basename(sys.argv[2])
    #try:
    #    title += "  vs.  "+os.path.basename(sys.argv[3])
    #except IndexError:
    #    pass

    ctLines = ctFileA.readlines()
    seq, pairedNucA = parseCTFile(ctLines)
    
    if haveSecondCT == True:
        ctLines = ctFileB.readlines()
        trash, pairedNucB = parseCTFile(ctLines)
    else:
        pairedNucB = []
    
    writePlot(outPath=outPath,title="",seq=seq,pairedNucA=pairedNucA,pairedNucB=pairedNucB,alpha=alpha)

