import os
import sample
from extractNclrData import extractNclr
from estimateOmega import estimateAllLoci
from estimateGeneTrees import estimateGeneTrees


def autorun(dataDir, resultRoot, taxa, ctlPrototype, codon="123", codonParition=False):
    dataDir = os.path.abspath(dataDir)
    resultRoot = os.path.abspath(resultRoot)
    ctlPrototype = os.path.abspath(ctlPrototype)
    outBase = ".".join([taxon[0] for taxon in taxa])
    print outBase
    alnDir = os.path.join (resultRoot, "I_Alignments", outBase)
    if codon == "123":
        omegaDir = os.path.join(resultRoot, "II_Omega", outBase)
    elif codon == "12":
        omegaDir = os.path.join(resultRoot, "II_Omega", outBase+"12")
    elif codon == "3":
        omegaDir = os.path.join(resultRoot, "II_Omega", outBase+"3")
    else:
        raise IOError("Do not know how to assign omega directory to codon "+codon)
    if codon != "123" and not os.path.isdir(omegaDir):
        raise IOError("The omega information of 123 codon is need")
    treeDir = os.path.join(resultRoot, "III_Gene_Trees", outBase)

    taxaCommonList = [commonName for commonName, abbrName in taxa]
    taxaAbbrList = [abbrName for commonName, abbrName in taxa]

    print "Extracting alignments..."
    extractNclr(dataDir, alnDir, taxaAbbrList, taxaCommonList, codon=codon)
    if codon == "123":
        print "Calculating omega..."
        estimateAllLoci(alnDir, ctlPrototype, omegaDir)
    print "Inferring gene trees..."
    estimateGeneTrees(alnDir, omegaDir, treeDir, codon=codon, codonParition=codonParition)
    print "GaoDing!"


def run_Primate_Scandentia_Glires(dataDir, resultRoot, ctlPrototype):
    primateList = sample.psgInfo["Primate"]
    scandentiaList = sample.psgInfo["Scandentia"]
    gliresList = sample.psgInfo["Glires"]
    outgroupList = sample.psgInfo["Outgroup"]

    for primate in primateList:
        for scandentia in scandentiaList:
            for glires in gliresList:
                for outgroup in outgroupList:
                    taxa = [primate, scandentia, glires, outgroup]
                    autorun(dataDir, resultRoot, taxa, ctlPrototype)

def run_Primate_Scandentia_Glires_12(dataDir, resultRoot, ctlPrototype, codon="12"):
    primateList = sample.psgInfo["Primate"]
    scandentiaList = sample.psgInfo["Scandentia"]
    gliresList = sample.psgInfo["Glires"]
    outgroupList = sample.psgInfo["Outgroup"]

    for primate in primateList:
        for scandentia in scandentiaList:
            for glires in gliresList:
                for outgroup in outgroupList:
                    taxa = [primate, scandentia, glires, outgroup]
                    autorun(dataDir, resultRoot, taxa, ctlPrototype)

def run_Primate_Scandentia_Glires_3(dataDir, resultRoot, ctlPrototype, codon="3"):
    primateList = sample.psgInfo["Primate"]
    scandentiaList = sample.psgInfo["Scandentia"]
    gliresList = sample.psgInfo["Glires"]
    outgroupList = sample.psgInfo["Outgroup"]

    for primate in primateList:
        for scandentia in scandentiaList:
            for glires in gliresList:
                for outgroup in outgroupList:
                    taxa = [primate, scandentia, glires, outgroup]
                    autorun(dataDir, resultRoot, taxa, ctlPrototype)

def run_Boreoeutheria_Afrotheria_Xenarthra(dataDir, resultRoot, ctlPrototype):
    boreoeutheriaList = sample.baxInfo["Boreoeutheria"]
    afrotheriaList = sample.baxInfo["Afrotheria"]
    xenarthraList = sample.baxInfo["Xenarthra"]
    outgroupList = sample.baxInfo["Outgroup"]

    for boreoeutheria in boreoeutheriaList:
        for afrotheria in afrotheriaList:
            for xenarthra in xenarthraList:
                for outgroup in outgroupList:
                    taxa = [boreoeutheria, afrotheria, xenarthra, outgroup]
                    autorun(dataDir, resultRoot, taxa, ctlPrototype)

def run_Perissodactyla_Chiroptera_Carnivora(dataDir, resultRoot, ctlPrototype):
    perissodatylaList = sample.pccInfo["Perissodactyla"]
    chiropteraList = sample.pccInfo["Chiroptera"]
    carnivoraList = sample.pccInfo["Carnivora"]
    outgroupList = sample.pccInfo["Outgroup"]

    for perissodatyla in perissodatylaList:
        for chiroptera in chiropteraList:
            for carnivora in carnivoraList:
                for outgroup in outgroupList:
                    taxa = [perissodatyla, chiroptera, carnivora, outgroup]
                    autorun(dataDir, resultRoot, taxa, ctlPrototype)

def run_Perissodactyla_Chiroptera_Cetartiodactyla(dataDir, resultRoot, ctlPrototype):
    perissodatylaList = sample.pcc2Info["Perissodactyla"]
    chiropteraList = sample.pcc2Info["Chiroptera"]
    cetartiodactylaList = sample.pcc2Info["Cetartiodactyla"]
    outgroupList = sample.pcc2Info["Outgroup"]

    for perissodatyla in perissodatylaList:
        for chiroptera in chiropteraList:
            for cetartiodactyla in cetartiodactylaList:
                for outgroup in outgroupList:
                    taxa = [perissodatyla, chiroptera, cetartiodactyla, outgroup]
                    autorun(dataDir, resultRoot, taxa, ctlPrototype)


if __name__ == '__main__':
    dataDir = r"/home/hc/Neutral_Phylogeny/Mammal_CDS"
    resultRoot = os.path.join("..", "..", "Results", "Empirical_results", "Primate_Scandentia_Glires")
    ctlPrototype = os.path.join("..", "codeml.ctl")
    run_Primate_Scandentia_Glires(dataDir, resultRoot, ctlPrototype)
