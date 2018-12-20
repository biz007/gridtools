
# Setup Conda environment
```
conda create -n grid python=3.6
source activate grid
conda install -c bioconda bwa==0.7.17 samtool==1.9 cutadapt==1.18
conda install numpy scipy h5py pandas==0.23.4 pysam==0.15
conda install snakemake
```
