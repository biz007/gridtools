"""
GRID-seq workflow in Nat. Protocol
"""
import os

# Setting of genome
# cd ~; git clone gridtools;
base = os.path.join(os.environ["HOME"], "gridtools/pipeline")
genomeDir = os.path.join(base, "genome") # Genome and annotation
dataDir = os.path.join(base, "data")     # Experimental raw data
girdDir = os.path.join(base, "gridcfg")  # GRID library configuration
gtf = os.path.join(genomeDir, "wgEncodeGencodeBasicVM18.gtf.gz") # Gene

# Setting per project
outDir = os.path.join(base, "results") # Results directory 
gridIdx = os.path.join(outDir, "grid/GRIDv1")
genomeIdx = os.path.join(outDir, "genome/mm10")
samples = ["test", ]


rule all:
    input:
        expand(os.path.join(outDir, "trim/{sample}.trimed.fq.gz"), sample=samples),
        gridIdx + ".bwt",
        expand(os.path.join(outDir, "linkerAlign/{sample}.bam"),  sample=samples),
        genomeIdx + ".bwt",
        expand(os.path.join(outDir, "matefq/{sample}.h5"), sample=samples),
        expand(os.path.join(outDir, "matefq/{sample}.mate.fq.gz"), sample=samples),
        expand(os.path.join(outDir, "mapped/{sample}.mate.mrk.bam"), sample=samples),
        os.path.join(outDir, "gene/mm10.gtf.gz"),
        expand(os.path.join(outDir, "qcAlignment/{sample}.h5"), sample=samples),
        expand(os.path.join(outDir, "qcStats/{sample}.stats.bases.txt"), sample=samples),
        expand(os.path.join(outDir, "RNA/{sample}.gene_expr.txt"), sample=samples),
        expand(os.path.join(outDir, "RNA/{sample}.gene_scope.txt"), sample=samples),
        expand(os.path.join(outDir, "DNA/{sample}.dna.txt.gz"), sample=samples),
        expand(os.path.join(outDir, "matrix/{sample}.matrix.gz"), sample=samples),
      

rule trim:
    """Step 153: Trim empty bases and the adapter sequences"""
    input:
        fq = os.path.join(dataDir, "{sample}.fq.gz"),
    output:
        fq = os.path.join(outDir, "trim/{sample}.trimed.fq.gz"),
    shell:
        """
        cutadapt -j 11 -a AGATCGGA --trim-n -o {output.fq} {input.fq}
        """

rule index_grid:
    """Step 154: Construct the BWA index file of the linker sequence"""
    input:
        fa = os.path.join(girdDir, "GRIDv1.fa")
    output:
        idx = gridIdx + ".bwt"
    params:
        outbase = gridIdx
    shell:
        """
        bwa index -p {params.outbase} {input.fa}
        """

rule align_linker:
    """Step 155: """
    input:
        fq = os.path.join(outDir, "trim/{sample}.trimed.fq.gz"),
        linkeridx = gridIdx + ".bwt",
    output:
        bam = os.path.join(outDir, "linkerAlign/{sample}.bam"),
        bai = os.path.join(outDir, "linkerAlign/{sample}.bam.bai"),
    params:
        linkeridxbase = gridIdx
    threads: 12
    shell:
        """
        bwa mem -t {threads} -k 5 -L 4 -B 2 -O 4 -o {output.bam}.sam {params.linkeridxbase} {input.fq} 
        samtools view -ub -o {output.bam}.bam {output.bam}.sam
        samtools sort -@ {threads} -m 4G -o {output.bam} {output.bam}.bam 
        samtools index -@ {threads} {output.bam}
        """    

rule separate_RNA_DNA:
    """Step 156: """
    input:
        bam = os.path.join(outDir, "linkerAlign/{sample}.bam"),
        bai = os.path.join(outDir, "linkerAlign/{sample}.bam.bai"),
    output:
        h5 = os.path.join(outDir, "matefq/{sample}.h5"),
        mate = os.path.join(outDir, "matefq/{sample}.mate.fq.gz"), 
    shell:
        """
        ./bin/GridTools.py matefq -l 19 -n {wildcards.sample} -o {output.h5} {input.bam} | gzip -c > {output.mate}
        """

rule index_genome:
    """Step 157a"""
    input:
        faDir = os.path.join(genomeDir, "mm10chromFa"),
    output:
        idx = genomeIdx + ".bwt"
    params:
        resdir = os.path.join(outDir, "genome"),
        outbase = genomeIdx
    shell:
        """
        cat {input.faDir}/*.fa > {params.resdir}/genome.fa
        bwa index -p {params.outbase} {params.resdir}/genome.fa
        """

rule align_genome:
    """Step 157b: Align read pairs to the genome, filter unambiguously mapped reads"""
    input:
        fq = os.path.join(outDir, "matefq/{sample}.mate.fq.gz"), 
        idx = genomeIdx + ".bwt",
    output:
        srt_bam = os.path.join(outDir, "mapped/{sample}.mate.srt.bam"),
    params:
        idxbase = genomeIdx,
    threads: 12
    shell:
        """
        bwa mem -o {output.srt_bam}.sam -t {threads} -k 7 -U 0 -T 0 -ap {params.idxbase} {input.fq} 
        samtools view -@ {threads} -u -b -q 1 -o {output.srt_bam}.bam {output.srt_bam}.sam
        samtools sort -n -@ {threads} -m 1G -o {output.srt_bam} {output.srt_bam}.bam
        """

rule fixmate:
    """Step 157c: mark duplicated read mates"""
    input:
        srt_bam = os.path.join(outDir, "mapped/{sample}.mate.srt.bam"),
    output:
        fix_bam = os.path.join(outDir, "mapped/{sample}.mate.fix.bam"),
        mrk_bam = os.path.join(outDir, "mapped/{sample}.mate.mrk.bam"),
    threads: 12
    shell:
        """
        samtools fixmate -@ {threads} -pm -O BAM {input.srt_bam} {output.fix_bam}.bam
        samtools sort -@ {threads} -m 1G -o {output.fix_bam} {output.fix_bam}.bam
        samtools markdup -@ {threads} {output.fix_bam} {output.mrk_bam}
        """

rule index_gtf:
    """Step 158a:"""
    input:
        gtf = gtf
    output:
        gtf = os.path.join(outDir, "gene/mm10.gtf.gz"),
    shell:
        """
        gunzip -c {input.gtf} | grep -v "^#" | sort -k1,1 -k4,4n | bgzip > {output.gtf}
        tabix -p gff {output.gtf}
        """

rule QC_alignment:
    """Step 158"""
    input:
        gtf = os.path.join(outDir, "gene/mm10.gtf.gz"),
        mrk_bam = os.path.join(outDir, "mapped/{sample}.mate.mrk.bam"),
        h5 = os.path.join(outDir, "matefq/{sample}.h5"),
    output:
        h5 = os.path.join(outDir, "qcAlignment/{sample}.h5"), # update the hdf5 file
    params:
        odir = os.path.join(outDir, "qcAlignment"),
    shell:
        """
        mkdir -p {params.odir}
        cp {input.h5} {output.h5}
        ./bin/GridTools.py evaluate -k 1 -m 10 -g {input.gtf} -o {output.h5} {input.mrk_bam}
        """

rule QC_stats:
    """Step 159"""
    input:
        h5 = os.path.join(outDir, "qcAlignment/{sample}.h5"),
    output:
        base = os.path.join(outDir, "qcStats/{sample}.stats.bases.txt"),
        count = os.path.join(outDir, "qcStats/{sample}.stats.counts.txt"),
        length = os.path.join(outDir, "qcStats/{sample}.stats.lengths.txt"),
        qual = os.path.join(outDir, "qcStats/{sample}.stats.quals.txt"),
        resolution = os.path.join(outDir, "qcStats/{sample}.stats.resolution.txt"),
    param:
        prefix = os.path.join(outDir, "qcStats/{sample}.stats"),
    shell:
        """
        ./bin/GridTools.py stats -clbqr -p {params.prefix} {input.h5}
        """

rule quantify_RNA:
    """Step 160"""
    input:
        h5 = os.path.join(outDir, "qcAlignment/{sample}.h5"),
    output:
        expr = os.path.join(outDir, "RNA/{sample}.gene_expr.txt"), 
        scope = os.path.join(outDir, "RNA/{sample}.gene_scope.txt"),  
    shell:
        """
        ./bin/GridTools.py RNA -e {output.expr} -s {output.scope} {input.h5}
        """

rule quantify_DNA:
    """Step 161"""
    input:
        h5 = os.path.join(outDir, "qcAlignment/{sample}.h5"),
    output:
        dna = os.path.join(outDir, "DNA/{sample}.dna.txt.gz"), 
    shell:
        """
        ./bin/GridTools.py DNA {input.h5} | gzip -c > {output.dna} 
        """

rule quantify_matrix:
    """Step 162"""
    input:
        h5 = os.path.join(outDir, "qcAlignment/{sample}.h5"),
    output:
        mat = os.path.join(outDir, "matrix/{sample}.matrix.gz"), 
    shell:
        """
        ./bin/GridTools.py matrix -k 100 -x 10 {input.h5} | gzip -c > {output.mat} 
        """
