import os, shutil, glob, subprocess, multiprocessing, time
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment


def prepareCodeMLAlignment(alnFile, outDir, id_function = lambda i:i.split("|")[-1]):
    if not os.path.isfile(alnFile):
        raise IOError(alnFile+" doesn\'t exist")

    aln = AlignIO.read(alnFile, "fasta")
    outFile = os.path.join(outDir, "for_codeml.nuc")
    records = list(aln)

    for i in range(len(records)):
        triplet = records[i].seq[-3:]
        if "-" not in triplet:            
            if triplet.translate() == "*":
                records[i].seq = records[i].seq[:-3]+Seq("???")
        for site in range(-6, -aln.get_alignment_length(), -3):
            triplet = records[i].seq[site:site+3]
            if "-" not in triplet:
                if triplet.translate() == "*":            
                    records[i].seq = records[i].seq[:site]+Seq("???")+records[i].seq[site+3:]


    newAln = MultipleSeqAlignment(records)

    with open(outFile,"w") as outHandle:
        outHandle.write(" "+str(len(aln))+" "+str(aln.get_alignment_length())+" GC\n")
        for record in aln:
            ID = id_function(record.id)
            outHandle.write(ID+(" "*(50-len(ID)))+str(record.seq+"\n"))


def prepareCodeMLCtl(ctlPrototype, out):
    shutil.copy(ctlPrototype, out)


def estimateOmega(alnFile, result, ctlPrototype):
    if not os.path.isfile(alnFile):
        raise IOError(alnFile+" doesn\'t exist")
    if not os.path.isdir(result):
        os.makedirs(result)

    prepareCodeMLAlignment(alnFile, result)
    prepareCodeMLCtl(ctlPrototype, result)

    os.chdir(result)
    p = subprocess.Popen(["codeml"], stdout=subprocess.PIPE)
    p.communicate()


def estimateAllLoci(alnDir, ctlPrototype, outDir):
    if not os.path.isdir(alnDir):
        raise IOError(alnDir+" doesn\'t exist")
    if os.path.isdir(outDir):
        return False

    pool = multiprocessing.Pool(processes=30)
    for cds in os.listdir(alnDir):
        alnFile = os.path.join(alnDir, cds, "alignment.fasta")
        result = os.path.join(outDir, cds)
        pool.apply_async(estimateOmega, (alnFile, result, ctlPrototype,))
        # estimateOmega(alnFile, ctlPrototype, result)
    pool.close()
    pool.join()





    
