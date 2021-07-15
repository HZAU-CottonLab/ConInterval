## ConInterval


> Scan conserved intervals based on the collinearity of genes.
> such as QTLs, DHS peaks.


### 1. Obtain the collinearity gene Id 

```bash
 python ./utils/collinear_interval.py  -r ./testData/request_QTL.txt -g ./testData/homoeolog_gene_coordinate.txt  -o  ./testData/request_flank_homoeolog_gene.txt
```

### 2. Obtain the align Interval
    need the chromosomes size file 
    1. `-rf` flank length of request Interval.
    2. `-tf` flank length of target Interval, default 1M.

```bash
python ./utils/align_coordinate.py  -c ./testData/chromsome_len.txt  -g ./testData/gene_coordiante.txt -f ./testData/request_flank_homoeolog_gene.txt  -rf 100 -o ./testData/align_coordinate.txt 
```
### 3. extract sequence of request Interval

```bash
python ./utils/extract_sequence.py  -g Ghirsutum_genome.fasta -b align_coordinate.txt  -o align_coordinate.fa 
```

### 4. Blast the request sequence to genomes

```bash
module load BLAST+/2.10.1 
blastn -query  align_coordinate.fa  -db Blast/Ghir_HZAU_V1.1  -qcov_hsp_perc 80  -outfmt '6 qseqid sseqid qstart qend sstart send evalue bitscore length nident mismatch gaps qcovs qcovhsp' -strand both  -evalue 1e-5 -out  align_coordinate.blast "
```

### 5. filter the Blast out 

```bash
    python conserved_fragment.py ./testData/align_coordinate.txt  align_coordinate.blast align_conserved.txt  
```
