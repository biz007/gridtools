# Pipeline to analyze test data
# Version 1.1

1. Download the *GridTools* and make it available
```
    git clone https://github.com/biz007/gridtools.git
    cp gridtools/gridtools/GridTools.py ~/bin/GridTools.py
    cd gridtools/pipeline/
```

2. Download genome and gene annotation
```
    mkdir genome && cd genome
    wget -v ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/GRCm38.primary_assembly.genome.fa.gz
    wget -v ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/gencode.vM19.annotation.gtf.gz
    cd ..
```

3. Download the GRID-seq test data
```
    mkdir data && cd data
    wget -v http://fugenome.ucsd.edu/gridseq/datasets/gridseq.test10M.raw.fq.gz
    ln -s gridseq.test10M.raw.fq.gz test.fq.gz
    cd ..
```

4. Setup *conda* environment and install required packages
```
    conda create -n grid python=3.6
    source activate grid
    conda install -c bioconda bwa samtools cutadapt numpy scipy h5py pandas pysam
    conda install snakemake
```

5. Run the pipeline
```
    source activate grid
    snakemake -s Snakefile.py
```

