"""
Protein Sequencing Project
Name: David Bautista
"""

project = "ProteinSeq"

#### CHECK-IN 1 ####
def readFile(filename):
    f=open(filename, "r")
    GeneInfoV1=f.read()
    GeneInfoV2=GeneInfoV1.replace("\n","")
    return GeneInfoV2
    
def dnaToRna(dna, startIndex):
    RNA=[]
    for i in range(startIndex,len(dna),3):
        Codon=dna[i:i+3]
        Codon=Codon.replace("T","U")
        RNA.append(Codon)
        if Codon=="UAA" or Codon=="UAG" or Codon=="UGA":
            return RNA
    return RNA

def makeCodonDictionary():
    import json
    f=open("codon_table.json","r")
    content=f.read()
    aminoD=json.loads(content)
    CodonD={}
    for amino in aminoD:
        for dnaTriplet in aminoD[amino]:
            Codon=dnaTriplet.replace("T","U")
            CodonD[Codon]=amino
    return CodonD

def generateProtein(codons, codonD):
    Protein=[]
    for i in range(len(codons)):
        if codons[i]=="AUG" and i==0:
            Protein.append("Start")
        else:
            Protein.append(codonD[codons[i]])
    return Protein

def synthesizeProteins(filename):
    DNA=readFile(filename)
    CodonD=makeCodonDictionary()
    Proteins=[]
    print(len(DNA))
    i=0
    while i<(len(DNA)-2):
        if (DNA[i]+DNA[i+1]+DNA[i+2])=="ATG":
            RNA=dnaToRna(DNA,i)
            Protein=generateProtein(RNA,CodonD)
            Proteins.append(Protein)
            i=i+(3*len(RNA))
        else:
            i+=1
    print("Unused-base count: "+str(i)+"\n"+"Length of Synthesized Proteins: "+str(len(Proteins))+"\n"+"Synthesized Proteins: "+str(Proteins))
    return Proteins
    
def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("human_p53.txt")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("elephant_p53.txt")

### WEEK 1 TEST FUNCTIONS ###

def testReadFile():
    print("Testing readFile()...", end="")
    text = readFile("human_p53.txt")
    # If the length is not correct, check that you're
    # removing newlines, and that you copied the whole sequence
    assert(len(text) == 19149)
    assert(text[:10] == "GATGGGATTG")
    print("... done!")

def testDnaToRna():
    print("Testing dnaToRna()...", end="")
    # Test a basic sequence
    dna = "ATGGATGGACTCTAA"
    assert(dnaToRna(dna, 0) == ["AUG", "GAU", "GGA", "CUC", "UAA"])
    # Test two mRNA strands in a row, with a random codon in between
    dna = "ATGGATGGACTCTAACTCATGCCCTTTTAG"
    assert(dnaToRna(dna, 0) == ["AUG", "GAU", "GGA", "CUC", "UAA"])
    assert(dnaToRna(dna, 18) == ["AUG", "CCC", "UUU", "UAG"])
    # Test a DNA strand that doesn't end properly
    dna = "CCTATGGACCAT"
    assert(dnaToRna(dna, 3) == ["AUG", "GAC", "CAU"])
    # Test a DNA strand with random bases in between
    dna = "ATGGATGGACTCTAACGCAATGCCCTTTTAG"
    assert(dnaToRna(dna, 0) == ["AUG", "GAU", "GGA", "CUC", "UAA"])
    assert(dnaToRna(dna, 19) == ["AUG", "CCC", "UUU", "UAG"])
    print("... done!")

def testMakeCodonDictionary():
    print("Testing makeCodonDictionary()...", end="")
    d = makeCodonDictionary()
    assert(d["AAA"] == "Lys")
    assert(d["GGA"] == "Gly")
    assert(d["AUG"] == "Met")
    assert(d["UAA"] == "Stop")
    print("... done!")

def testGenerateProtein():
    print("Testing generateProtein()...", end="")
    codonD = makeCodonDictionary()
    rna = ["AUG", "GAU", "GGA", "CUC", "UAA"]
    assert(generateProtein(rna, codonD) == ["Start", "Asp", "Gly", "Leu", "Stop"])
    rna = ["AUG", "CCC", "UUU", "UAG"]
    assert(generateProtein(rna, codonD) == ["Start", "Pro", "Phe", "Stop"])
    rna = ["AUG", "GAC", "CAU"]
    assert(generateProtein(rna, codonD) == [ "Start", "Asp", "His"])
    print("... done!")
    
def week1Tests():
    testReadFile()
    testDnaToRna()
    testMakeCodonDictionary()
    testGenerateProtein()

week1Tests()
runWeek1()

#### CHECK-IN 2 ####

def commonProteins(proteinList1, proteinList2):
    UniqueProteins=[]
    for protein1 in proteinList1:
        for protein2 in proteinList2:
            if protein1==protein2:
                if protein1 not in UniqueProteins:
                    UniqueProtein=protein1
                    UniqueProteins.append(UniqueProtein)
    return UniqueProteins

def combineProteins(proteinList):
    aminoList=[]
    for protein in proteinList:
        for amino in protein:
            aminoList.append(amino)
    return aminoList

def aminoAcidDictionary(aaList):
    aaDictionary={}
    for amino in aaList:
        if amino not in aaDictionary:
            aaDictionary[amino]=1
        else:
            aaDictionary[amino]=aaDictionary[amino]+1
    return aaDictionary

def sortAminoAcidsByFreq(aaList):
    total=len(aaList)
    aaDictionary=aminoAcidDictionary(aaList)
    FrequencyList=[]
    for amino in aaDictionary:
        if amino!="Start" and amino!="Stop":
            frequency=aaDictionary[amino]/total
            FrequencyList.append([frequency,amino])
            FrequencyList.sort()
    return FrequencyList

def findAminoAcidDifferences(proteinList1, proteinList2):
    Aminos1=combineProteins(proteinList1)
    Aminos2=combineProteins(proteinList2)
    FreqList1=sortAminoAcidsByFreq(Aminos1)
    FreqList2=sortAminoAcidsByFreq(Aminos2)
    DiffList=[]
    for Freq1 in FreqList1:
        for Freq2 in FreqList2:
            if Freq1[1]==Freq2[1]:
                diff=abs(Freq1[0]-Freq2[0])
                if FreqList1.index(Freq1)!=FreqList2.index(Freq2) \
                and diff>=0.005:
                    DiffList.append([Freq1[1],Freq1[0],Freq2[0]])
    return DiffList

def displayTextResults(commonalities, differences):
    NoRepeats=[]
    for protein in commonalities:
        if ["Start","Stop"]==protein:
            commonalities.remove(["Start","Stop"])
        else:
            for amino in protein:
                if amino not in NoRepeats and amino!="Start" and \
                amino!="Stop":
                    NoRepeats.append(amino)
    for AMINO in NoRepeats:
        print(AMINO+"\n")
    for rates in differences:
        AminoAcid=rates[0]
        Rate1=str(rates[1]*100)+"%"
        Rate2=str(rates[2]*100)+"%"
        print(AminoAcid+":"+Rate1+" in Human Gene,",Rate2+" in Elephant Gene")
    return

def runWeek2():
    humanProteins = synthesizeProteins("human_p53.txt")
    elephantProteins = synthesizeProteins("elephant_p53.txt")
    
    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins)
    displayTextResults(commonalities, differences)

#### WEEK 2 TESTS ####

def testCommonProteins():
    print("Testing commonProteins()...", end="")
    plist1 = [ [ "Start", "Pro", "Val", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "His", "Stop" ] ]
    plist2 = [ [ "Start", "Cys", "Cys", "Tyr", "Stop" ], ["Start", "Glu", "Asp", "Stop" ],
               [ "Start", "His", "Stop" ], [ "Start", "Stop" ], [ "Start", "Met", "Leu", "Stop" ] ]
    plist3 = [ [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Lys", "Stop" ], 
               [ "Start", "Asn", "Asn", "Asn", "Asn", "Stop" ] ]
    assert(commonProteins(plist1, plist2) == [ [ "Start", "His", "Stop" ] ])
    assert(sorted(commonProteins(plist1, plist3)) == [ [ "Start", "Asp", "Glu", "Stop" ], 
                                                       [ "Start", "Phe", "Stop" ] ])
    assert(commonProteins(plist2, plist3) == [ ])
    print("... done!")

def testCombineProteins():
    print("Testing combineProteins()...", end="")
    plist1 = [ [ "Start", "Pro", "Val", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "His", "Stop" ] ]
    plist2 = [ [ "Start", "Cys", "Cys", "Tyr", "Stop" ], ["Start", "Glu", "Asp", "Stop" ],
               [ "Start", "His", "Stop" ], [ "Start", "Stop" ], [ "Start", "Met", "Leu", "Stop" ] ]
    plist3 = [ [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Lys", "Stop" ], 
               [ "Start", "Asn", "Asn", "Asn", "Asn", "Stop" ] ]
    assert(combineProteins(plist1) == [ "Start", "Pro", "Val", "Stop", "Start", 
                                        "Phe", "Stop", "Start", "Asp", "Glu", 
                                        "Stop", "Start", "His", "Stop" ])
    assert(combineProteins(plist2) == [ "Start", "Cys", "Cys", "Tyr", "Stop", 
                                        "Start", "Glu", "Asp", "Stop", "Start", 
                                        "His", "Stop", "Start", "Stop", "Start", 
                                        "Met", "Leu", "Stop" ])
    assert(combineProteins(plist3) == [ "Start", "Asp", "Glu", "Stop",  "Start", 
                                        "Phe", "Stop", "Start", "Asp", "Glu", 
                                        "Stop", "Start", "Lys", "Stop", "Start", 
                                        "Asn", "Asn", "Asn", "Asn", "Stop" ])
    print("... done!")

def testAminoAcidDictionary():
    print("Testing aminoAcidDictionary()...", end="")
    aaList1 = [ "Start", "Pro", "Val", "Stop", "Start", "Phe", "Stop", "Start", 
                "Asp", "Glu", "Stop", "Start", "His", "Stop" ]
    aaList2 = [ "Start", "Cys", "Cys", "Tyr", "Stop", "Start", "Glu", "Asp", 
                "Stop", "Start", "His", "Stop", "Start", "Stop", "Start", "Met", 
                "Leu", "Stop" ]
    aaList3 = [ "Start", "Asp", "Glu", "Stop",  "Start", "Phe", "Stop", "Start", 
                "Asp", "Glu", "Stop", "Start", "Lys", "Stop", "Start", "Asn", 
                "Asn", "Asn", "Asn", "Stop" ]
    assert(aminoAcidDictionary(aaList1) == { "Start" : 4, "Pro" : 1, "Val" : 1,
                "Stop" : 4, "Phe" : 1, "Asp" : 1, "Glu" : 1, "His" : 1 })
    assert(aminoAcidDictionary(aaList2) == { "Start" : 5, "Cys" : 2, "Tyr" : 1,
                "Stop" : 5, "Glu" : 1, "Asp" : 1, "His" : 1, "Met" : 1, "Leu" : 1 })
    assert(aminoAcidDictionary(aaList3) == { "Start" : 5, "Asp" : 2, "Glu" : 2, 
                "Stop" : 5, "Phe" : 1, "Lys" : 1, "Asn" : 4 })
    print("... done!")

def testSortAminoAcidsByFreq():
    print("Testing sortAminoAcidsByFreq()...", end="")
    aaList1 = [ "Start", "Pro", "Val", "Stop", "Start", "Phe", "Stop", "Start", 
                "Asp", "Glu", "Stop", "Start", "His", "Stop" ]
    aaList2 = [ "Start", "Cys", "Cys", "Tyr", "Stop", "Start", "Glu", "Asp", 
                "Stop", "Start", "His", "Stop", "Start", "Stop", "Start", "Met", 
                "Leu", "Stop" ]
    aaList3 = [ "Start", "Asp", "Glu", "Stop",  "Start", "Phe", "Stop", "Start", 
                "Asp", "Glu", "Stop", "Start", "Lys", "Stop", "Start", "Asn", 
                "Asn", "Asn", "Asn", "Stop" ]
    result1 = sortAminoAcidsByFreq(aaList1)
    assert(len(result1) == 6)
    # All the amino acids occur once, so they should all have the same frequency
    for i in range(len(result1)):
        assert(result1[0][0] == result1[i][0])
    result2 = sortAminoAcidsByFreq(aaList2)
    assert(len(result2) == 7)
    # All but one of the amino acids occur once
    for i in range(len(result2)-1):
        assert(result2[0][0] == result2[i][0])
    # And Cys occurs twice
    assert(result2[len(result2)-1][0] == 2*result2[0][0])
    assert(result2[len(result2)-1][1] == "Cys")
    # The last one has simple frequencies, we can check it directly
    result3 = sortAminoAcidsByFreq(aaList3)
    assert((result3[0:2] == [[0.05, 'Lys'], [0.05, 'Phe']]) or \
           (result3[0:2] == [[0.05, 'Phe'], [0.05, 'Lys']]))
    assert((result3[2:4] == [[0.1, 'Asp'], [0.1, 'Glu']]) or \
           (result3[2:4] == [[0.1, 'Glu'], [0.1, 'Asp']]))
    assert(result3[4] == [0.2, 'Asn'])    
    print("... done!")

def testFindAminoAcidDifferences():
    print("Testing findAminoAcidDifferences()...", end="")
    plist1 = [ [ "Start", "Pro", "Val", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "His", "Stop" ] ]
    plist2 = [ [ "Start", "Cys", "Cys", "Tyr", "Stop" ], ["Start", "Glu", "Asp", "Stop" ],
               [ "Start", "His", "Stop" ], [ "Start", "Stop" ], [ "Start", "Met", "Leu", "Stop" ] ]
    plist3 = [ [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Lys", "Stop" ], 
               [ "Start", "Asn", "Asn", "Asn", "Asn", "Stop" ] ]
    result1 = findAminoAcidDifferences(plist1, plist2)
    assert(result1 == [])
    result2 = findAminoAcidDifferences(plist1, plist3)
    assert(len(result2) == 3)
    assert(sorted([result2[0][0], result2[1][0]]) == ["Asp", "Glu"])
    assert(result2[2][0] == "Phe")
    result3 = findAminoAcidDifferences(plist2, plist3)
    assert(len(result3) == 2)
    assert(sorted([result3[0][0], result3[1][0]]) == ["Asp", "Glu"])
    print("... done!")

def week2Tests():
    testCommonProteins()
    testCombineProteins()
    testAminoAcidDictionary()
    testSortAminoAcidsByFreq()
    testFindAminoAcidDifferences()

week2Tests()
runWeek2()

#### FULL ASSIGNMENT ####

def makeAminoAcidLabels(geneList):
    aaList2D=[]
    aaList1D=[]
    for gene in geneList:
        aaList2D.append(combineProteins(gene))
    for aminoList in aaList2D:
        for amino in aminoList:
            if amino not in aaList1D:
                aaList1D.append(amino)
                aaList1D.sort()
    return aaList1D

def setupChartData(labels, geneList):
    FreqList=[]
    for gene in geneList:
        CombinedProteins=combineProteins(gene)
        Total=len(CombinedProteins)
        TempFreqList=[]
        TempAadict=aminoAcidDictionary(CombinedProteins)
        for amino in TempAadict:
            TempAadict[amino]=TempAadict[amino]/Total
        for label in labels:
            if label in TempAadict:
                TempFreqList.append(TempAadict[label])
            else:
                TempFreqList.append(0)
        FreqList.append(TempFreqList)
    return FreqList
    
def createChart(xLabels, freqList, freqLabels, edgeList=None):
    import numpy as np
    import matplotlib.pyplot as plt
    x = np.arange(len(xLabels))
    width = 0.6/len(freqList)
    offset=-0.3+(width/2)
    fig, ax = plt.subplots()
        
        
    for i in range(len(freqList)):
        AminoMeans=freqList[i]
        rects=ax.bar(x+offset+i*width,AminoMeans, width)
        
        
    ax.set_ylabel("Frequency")
    ax.set_title("Protein Sequencing")
    ax.set_xticks(x+width/2)
    ax.set_xticklabels(xLabels)
    plt.show()
    return

def makeEdgeList(labels, biggestDiffs):
    ColorList=[]
    DiffAminos=[]
    for list in biggestDiffs:
        DiffAminos.append(list[0])
    for label in labels:
        if label in DiffAminos:
            ColorList.append("black")
        else:
            ColorList.append("white")
    return ColorList

def runFullProgram():
    humanGene=synthesizeProteins("human_p53.txt")
    elephantGene=synthesizeProteins("elephant_p53.txt")
    commonalities=commonProteins(humanGene,elephantGene)
    differences=findAminoAcidDifferences(humanGene,elephantGene)
    print(displayTextResults(commonalities, differences))
    genes=[humanGene,elephantGene]
    labels=makeAminoAcidLabels(genes)
    data=setupChartData(labels,genes)
    edgeColor=makeEdgeList(labels,differences)
    freqLabels=["Human Gene","Elephant Gene"]
    chart=createChart(labels,data,freqLabels,edgeList=edgeColor)
    return 

#### WEEK 3 TESTS ####

def testMakeAminoAcidLabels():
    print("Testing makeAminoAcidLabels()...", end="")
    plist1 = [ [ "Start", "Pro", "Val", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "His", "Stop" ] ]
    plist2 = [ [ "Start", "Cys", "Cys", "Tyr", "Stop" ], ["Start", "Glu", "Asp", "Stop" ],
               [ "Start", "His", "Stop" ], [ "Start", "Stop" ], [ "Start", "Met", "Leu", "Stop" ] ]
    plist3 = [ [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Lys", "Stop" ], 
               [ "Start", "Asn", "Asn", "Asn", "Asn", "Stop" ] ]
    geneList = [ plist1, plist2, plist3 ]
    assert(makeAminoAcidLabels(geneList) == [ "Asn", "Asp", "Cys", "Glu", "His", 
                "Leu", "Lys", "Met", "Phe", "Pro", "Start", "Stop", "Tyr", "Val" ])
    print("... done!")

def testSetupChartData():
    print("Testing setupChartData()...", end="")
    plist1 = [ [ "Start", "Pro", "Val", "Stop" ], [ "Start", "Phe", "Stop" ],
            [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "His", "Stop" ] ]
    plist2 = [ [ "Start", "Cys", "Cys", "Tyr", "Stop" ], ["Start", "Glu", "Asp", "Stop" ],
            [ "Start", "His", "Stop" ], [ "Start", "Stop" ], [ "Start", "Met", "Leu", "Stop" ] ]
    plist3 = [ [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Lys", "Stop" ], 
               [ "Start", "Asn", "Asn", "Asn", "Asn", "Stop" ] ]
    geneList = [ plist1, plist2, plist3 ]
    labels = makeAminoAcidLabels(geneList)
    result = setupChartData(labels, geneList)
    assert(len(result) == 3 and len(result[0]) == 14)
    assert(result[0][0] == 0 and result[1][0] == 0 and result[2][0] == 0.2)
    print("... done!")

def testCreateChart():
    print("Testing createChart()...", end="")
    plist1 = [ [ "Start", "Pro", "Val", "Stop" ], [ "Start", "Phe", "Stop" ],
            [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "His", "Stop" ] ]
    plist2 = [ [ "Start", "Cys", "Cys", "Tyr", "Stop" ], ["Start", "Glu", "Asp", "Stop" ],
            [ "Start", "His", "Stop" ], [ "Start", "Stop" ], [ "Start", "Met", "Leu", "Stop" ] ]
    plist3 = [ [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Lys", "Stop" ], 
               [ "Start", "Asn", "Asn", "Asn", "Asn", "Stop" ] ]
    geneList = [ plist1, plist2, plist3 ]
    labels = makeAminoAcidLabels(geneList)
    freqList = setupChartData(labels, geneList)
    freqLabels = ["Ex1", "Ex2", "Ex3"]
    createChart(labels, freqList, freqLabels)
    print("... check your chart!")

def testMakeEdgeList():
    print("Testing makeEdgeList()...", end="")
    plist1 = [ [ "Start", "Pro", "Val", "Stop" ], [ "Start", "Phe", "Stop" ],
            [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "His", "Stop" ] ]
    plist3 = [ [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Phe", "Stop" ],
               [ "Start", "Asp", "Glu", "Stop" ], [ "Start", "Lys", "Stop" ], 
               [ "Start", "Asn", "Asn", "Asn", "Asn", "Stop" ] ]
    geneList = [ plist1, plist3 ]
    labels = makeAminoAcidLabels(geneList)
    biggestDiffs = findAminoAcidDifferences(plist1, plist3)
    result = makeEdgeList(labels, biggestDiffs)
    assert(result == ['white', 'black', 'black', 'white', 'white', 'black', 
                      'white', 'white', 'white', 'white'])
    print("... done!")

def week3Tests():
    testMakeAminoAcidLabels()
    testSetupChartData()
    testCreateChart()
    testMakeEdgeList()

week3Tests()
runFullProgram()