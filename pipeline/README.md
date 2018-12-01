# Pipeline to analyze test data

1. Download genome and gene annotation

    + genome as Fasta
```
    mkdir -p genome
    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz genome/mm10chromFa.tar.gz
    mkdir -p genome/mm10chromFa; tar -xvf genome/mm10chromFa.tar.gz -C genome/mm10chromFa;
```
    + genome chrom sizes as csv file
```
    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes genome/
```

    + gene annotation as GTF file
```
    Export GTF of Gencode Basic VM18 gene track from UCSC Table Browser
    saved as genome/wgEncodeGencodeBasicVM18.gtf.gz
```

2. Download test GRID-seq data
```
    mkdir -p data
    wget http://fugenome.ucsd.edu/gridseq/datasets/gridseq.test10M.raw.fq.gz
    saved in data/gridseq.test10M.raw.fq.gz
    cd data; ln -s gridseq.test10M.raw.fq.gz test.fq.gz
```

3. Run Snakemake
```
    cd pipeline; mkdir -p bin; cp ../gridtools/GridTools.py bin/
    source activate grid
    snakemake -s Snakefile.py -np # Check
```

