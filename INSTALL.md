
# Setup Conda environment
```
conda create -n grid python=3.6.4
source activate grid
pip install numpy scipy pandas>=0.23.4 pysam==0.13 h5py snakemake 
pip install cutadapt
conda install bwa
conda install htslib==1.6 # compatible with pysam-0.13
conda install samtools
```
