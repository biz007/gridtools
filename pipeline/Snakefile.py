"""
GRID-seq workflow in Nat. Protocol
"""
import os

# Setting of genome
base = os.getcwd()
genomeDir = os.path.join(base, "genome") # Genome and annotation
dataDir = os.path.join(base, "data")     # Experimental raw data
cfgDir = os.path.join(base, "gridcfg")  # GRID library configuration

gtffile = os.path.join(genomeDir, "gencode.vM19.annotation.gtf.gz") # Gene gtf
genomefile = os.path.join(genomeDir, "GRCm38.primary_assembly.genome.fa.gz") # Genome fa

# Setting per project
outDir = os.path.join(base, "results") # Results directory 
gridIdx = os.path.join(outDir, "index/grid/GRIDv1")
genomeIdx = os.path.join(outDir, "index/genome/mm10")
samples = ["test", ]


rule all:
    input:
        expand(os.path.join(outDir, "{sample}.trimed.fq.gz"), sample=samples),
        gridIdx + ".bwt",
        expand(os.path.join(outDir, "mapped/{sample}.GRIDv1.bam"), sample=samples),
        genomeIdx + ".bwt",
        expand(os.path.join(outDir, "{sample}.mate.fq.gz"), sample=samples),
        expand(os.path.join(outDir, "mapped/{sample}.mate.mrk.bam"), sample=samples),
        os.path.join(outDir, "mm10.gtf.gz"),
        expand(os.path.join(outDir, "hdf/{sample}.bak.h5"), sample=samples),
        expand(os.path.join(outDir, "qcStats/{sample}.stats.counts.txt"), sample=samples),
        expand(os.path.join(outDir, "RNA/{sample}.gene_exprs.txt"), sample=samples),
        expand(os.path.join(outDir, "DNA/{sample}.dna.txt.gz"), sample=samples),
        expand(os.path.join(outDir, "matrix/{sample}.matrix.gz"), sample=samples)
        

rule trim:
    """Step 153: Trim empty bases and the adapter sequences"""
    input:
        fq = os.path.join(dataDir, "{sample}.fq.gz")
    output:
        fq = os.path.join(outDir, "{sample}.trimed.fq.gz")
    threads: os.cpu_count()
    shell:
        """
        cutadapt -j {threads} -l 86 --max-n 5 -o {output.fq} {input.fq}
        """

rule index_grid:
    """Step 154: Construct the BWA index file of the linker sequence"""
    input:
        fa = os.path.join(cfgDir, "GRIDv1.fa")
    output:
        idx = gridIdx + ".bwt"
    params:
        idx = gridIdx
    shell:
        """
        bwa index -p {params.idx} {input.fa}
        """

rule align_linker:
    """Step 155: """
    input:
        fq = os.path.join(outDir, "{sample}.trimed.fq.gz"),
        linkeridx = gridIdx + ".bwt"
    output:
        bam = os.path.join(outDir, "mapped/{sample}.GRIDv1.bam"),
        bai = os.path.join(outDir, "mapped/{sample}.GRIDv1.bam.bai")
    params:
        linkeridx = gridIdx
    threads: os.cpu_count()
    shell:
        """
        bwa mem -t {threads} -k 5 -L 4 -B 2 -O 4 {params.linkeridx} {input.fq} | \\
        samtools view -ub | samtools sort -@ {threads} -m 4G -o {output.bam}
        samtools index -@ {threads} {output.bam}
        """    

rule separate_RNA_DNA:
    """Step 156: """
    input:
        bam = os.path.join(outDir, "mapped/{sample}.GRIDv1.bam"),
        bai = os.path.join(outDir, "mapped/{sample}.GRIDv1.bam.bai")
    output:
        h5 = os.path.join(outDir, "hdf/{sample}.h5"),
        mate = os.path.join(outDir, "{sample}.mate.fq.gz")
    params:
        sample = "{sample}"
    shell:
        """
        GridTools.py matefq -l 19 -n {params.sample} -o {output.h5} {input.bam} | gzip -c > {output.mate}
        """

rule index_genome:
    """Step 157a"""
    input:
        fa = genomefile
    output:
        idx = genomeIdx + ".bwt"
    params:
        idx = genomeIdx
    shell:
        """
        bwa index -p {params.idx} {input.fa}
        """

rule align_genome:
    """Step 157b: Align read pairs to the genome, filter unambiguously mapped reads"""
    input:
        fq = os.path.join(outDir, "{sample}.mate.fq.gz")
    output:
        rbam = os.path.join(outDir, "mapped/{sample}.mate.srt.bam"),
        mbam = os.path.join(outDir, "mapped/{sample}.mate.mrk.bam")
    params:
        idx = genomeIdx,
    threads: os.cpu_count()
    shell:
        """
        bwa mem -t {threads} -k 17 -w 1 -T 1 -p {params.idx} {input.fq} | \\
        samtools sort -@ {threads} -n -m 1G | \\
        samtools fixmate -@ {threads} -pm - - | \\
        samtools view -@ {threads} -bf 0x1 | \\
        samtools sort -@ {threads} -m 1G | \\
        samtools markdup -@ {threads} - {output.rbam}
        samtools index -@ {threads} {output.rbam}

        samtools sort -@ {threads} -n -m 1G -o {output.mbam} {output.rbam}
        """

rule index_gtf:
    """Step 158a:"""
    input:
        gtf = gtffile
    output:
        gtf = os.path.join(outDir, "mm10.gtf.gz")
    shell:
        """
        zcat {input.gtf} | grep -v ^"#" | sort -k1,1 -k4,4n | bgzip > {output.gtf}
        tabix -p gff {output.gtf}
        """

rule QC_alignment:
    """Step 158"""
    input:
        gtf = os.path.join(outDir, "mm10.gtf.gz"),
        h5 = os.path.join(outDir, "hdf/{sample}.h5"),
        mbam = os.path.join(outDir, "mapped/{sample}.mate.mrk.bam")
    output:
        h5 = os.path.join(outDir, "hdf/{sample}.bak.h5") # update the hdf5 file
    params:
        bink = 10, # kilobase size of bin
        winm = 10  # number of moving windows
    shell:
        """
        cp {input.h5} {output.h5}
        GridTools.py evaluate -k {params.bink} -m {params.winm} -g {input.gtf} -o {input.h5} {input.mbam}
        """

rule QC_stats:
    """Step 159"""
    input:
        h5 = os.path.join(outDir, "hdf/{sample}.h5"),
    output:
        base = os.path.join(outDir, "qcStats/{sample}.stats.bases.txt"),
        count = os.path.join(outDir, "qcStats/{sample}.stats.counts.txt"),
        length = os.path.join(outDir, "qcStats/{sample}.stats.lengths.txt"),
        resolution = os.path.join(outDir, "qcStats/{sample}.stats.resolution.txt")
    params:
        prefix = os.path.join(outDir, "qcStats/{sample}.stats")
    shell:
        """
        GridTools.py stats -bclr -p {params.prefix} {input.h5}
        """

rule quantify_RNA:
    """Step 160"""
    input:
        h5 = os.path.join(outDir, "hdf/{sample}.h5"),
    output:
        exprs = os.path.join(outDir, "RNA/{sample}.gene_exprs.txt"), 
        scope = os.path.join(outDir, "RNA/{sample}.gene_scope.txt")
    shell:
        """
        GridTools.py RNA -e {output.exprs} -s {output.scope} {input.h5}
        """

rule quantify_DNA:
    """Step 161"""
    input:
        h5 = os.path.join(outDir, "hdf/{sample}.h5"),
    output:
        dna = os.path.join(outDir, "DNA/{sample}.dna.txt.gz"), 
    shell:
        """
        GridTools.py DNA {input.h5} | gzip -c > {output.dna} 
        """

rule quantify_matrix:
    """Step 162"""
    input:
        h5 = os.path.join(outDir, "hdf/{sample}.h5"),
    output:
        mtx = os.path.join(outDir, "matrix/{sample}.matrix.gz")
    params:
        grpk = 100, # cutoff of reads per kilobase in the genebody
        drpk = 10   # cutoff of reads per kilobase in the max bin
    shell:
        """
        GridTools.py matrix -k {params.grpk} -x {params.drpk} {input.h5} | gzip -c > {output.mtx} 
        """
