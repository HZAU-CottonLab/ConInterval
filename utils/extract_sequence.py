'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-15 14:33:32
LastEditors: zpliu
LastEditTime: 2021-07-15 16:39:40
@param: 
'''
import pandas as pd 
import pysam
import os 
import argparse


if __name__=="__main__":
        
    parser=argparse.ArgumentParser()
    parser.add_argument('-g',help='geneom sequence file *.fa')
    parser.add_argument('-b',help='request interval')
    parser.add_argument('-o',help='output file')
    args=parser.parse_args()
    # genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
    genomeFile = args.g 
    genomeObject = pysam.FastaFile(genomeFile) 

    ############################
    #sequence coordinate
    ############################
    QTLsite=pd.read_csv(args.b,header=None,index_col=None,sep="\t")
    out=[]
    #*index start with 0
    # Ghir_A01	164476	164577	QTL1*+	Ghir_D01	1	1162194	Ghir_D01G000260-Ghir_D01G000230*+
    for item in QTLsite.values:
        QTLid=item[3].strip("\*+")
        chromsome=item[0]
        start=item[1]-1
        end=item[2]-1
        pairedRegion=[chromsome,start,end]
        geneID=">"+QTLid+"-"+chromsome+":"+str(start+1)+"-"+str(end+1)+"\n"
        sequence=genomeObject.fetch(start=pairedRegion[1],end=pairedRegion[2],region=pairedRegion[0])
        out.append(geneID+sequence)
    with open(args.o,'w') as File:
        for line in out:
            File.write(line+"\n")
