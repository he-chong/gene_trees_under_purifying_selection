import platform
import numpy as np
import os, re, shutil, glob, sys, pprint
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment


def extractNclr(alnDir, outRoot, taxaList, newNameList=None, key_function=None, codon="123"):
    if not os.path.isdir(alnDir):
        raise IOError(alnDir+"directory doesn\'t exist.")
    if os.path.isdir(outRoot):
        return False

    if not newNameList:
        newNameList = taxaList
    count = 0
    for alnFile in glob.glob(os.path.join(alnDir, "*.fasta"))+glob.glob(os.path.join(alnDir, "*.fas")):
        alnBase = os.path.basename(alnFile)
        outDir = os.path.join(outRoot, alnBase.strip(".fasta"))
        if "MT-" not in alnBase:
            aln = AlignIO.read(alnFile, "fasta")
            if key_function:
                alnDict = SeqIO.to_dict(aln,key_function=key_function)
            else:
                alnDict = SeqIO.to_dict(aln,key_function=lambda r:r.id)
            outRecords = []
            # print alnDict.keys()
            if not newNameList:
                for taxon in taxaList:
                    if taxon in alnDict:
                        if codon == "123":
                            record = alnDict[taxon]
                        elif codon == "12":
                            record = alnDict[taxon]
                            record.seq = Seq(''.join(np.array(record.seq)[np.where((np.arange(len(record.seq))+1)%3)]))
                        outRecords.append(record)
            else:
                for taxon, newName in zip(taxaList, newNameList):
                    if taxon in alnDict:
                        if codon == "123":
                            record = alnDict[taxon]
                            record.id = newName
                            record.description = ""
                        elif codon == "12":
                            record = alnDict[taxon]
                            record.id = newName
                            record.description = ""
                            record.seq = Seq(''.join(np.array(record.seq)[np.where((np.arange(len(record.seq))+1)%3)]))
                        elif codon == "3":
                            record = alnDict[taxon]
                            record.id = newName
                            record.description = ""
                            iList = []
                            for i in np.arange(len(record.seq)):
                                if (i + 1) % 3 == 0:
                                    iList.append(i)
                            record.seq = Seq(''.join(np.array(record.seq)[iList]))
                        outRecords.append(record)
            if len(outRecords) == len(taxaList):
                outAln = MultipleSeqAlignment(outRecords)
                outAlnFile = os.path.join(outDir, "alignment.fasta")
                if not os.path.isdir(outDir):
                    os.makedirs(outDir)
                AlignIO.write(outAln, outAlnFile, "fasta")
                count += 1
    return count



            
                    



