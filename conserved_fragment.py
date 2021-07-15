'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-15 10:15:06
LastEditors: zpliu
LastEditTime: 2021-07-15 16:47:30
@param: 
'''
import pandas as pd
import pybedtools
import sys 
print("read file...")
QTL_align_region = pd.read_csv(
    sys.argv[1], header=None, index_col=None, sep="\t")

QTL_blastData = pd.read_csv(
    sys.argv[2], header=None, index_col=None, sep="\t")
print("complemt!")
out = []
for QTLitem in QTL_align_region.values:
    QTLchrom, QTLstart, QTLend, QTLId, targetChrom, targetstart, targetEnd, targeHomoeologene = QTLitem
    #! format the request QTL id
    requestQTLid = QTLId.split("*")[0]+"-" + \
        QTLchrom+":"+str(QTLstart)+"-"+str(QTLend)
    QTL_blast_item = QTL_blastData.loc[(QTL_blastData[0] == requestQTLid) & (
        QTL_blastData[1] == targetChrom)]
    # print(requestQTLid)
    # print(QTL_blast_item)
    if QTL_blast_item.empty:
        #! QTL without conserved fragment
        out.append(
            (QTLchrom, QTLstart, QTLend, QTLId,
             "-", "-", "-", targeHomoeologene, "-")
        )
    else:
        BlastRegion_str = ''
        for blastregion in QTL_blast_item.values:
            #! align region
            '''
            #* chromosome
            #* start
            #* end 
            #* bitscore
            '''
            if blastregion[4] <= blastregion[5]:
                BlastRegion_str += blastregion[1]+"\t" + str(blastregion[4])+"\t"+str(
                    blastregion[5])+"\t"+str(blastregion[7])+"\n"
            else:
                BlastRegion_str += blastregion[1]+"\t" + str(blastregion[5])+"\t"+str(
                    blastregion[4])+"\t"+str(blastregion[7])+"\n"
        #! intersectBed with targe reqion
        targetRegion_interval = pybedtools.BedTool(
            targetChrom+"\t"+str(targetstart)+"\t"+str(targetEnd),
            from_string=True
        )
        Blast_interval = pybedtools.BedTool(
            BlastRegion_str.strip("\n"), from_string=True)
        interSectBedOut = targetRegion_interval.intersect(
            Blast_interval, wb=True)
        '''
        Interval Interval bitscore
        '''
        if not interSectBedOut:
            #! not conserved fragment in collineral
            out.append(
                (QTLchrom, QTLstart, QTLend, QTLId,
                 "-", "-", "-", targeHomoeologene, "-")
            )
        else:
            #! intersect with collineral region
            for item in interSectBedOut:
                print(item)
                out.append((QTLchrom, QTLstart, QTLend, QTLId, str(item[0]),
                            str(item[1]), str(item[2]), targeHomoeologene, str(item[6])))

with open(sys.argv[3], 'w') as File:
    for item in out:
        File.write("\t".join(str(i) for i in item)+"\n")
