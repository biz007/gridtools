#!/bin/bash

idx_grid="GRIDv1"

gx="GRCm38"
idx_mm38="$genome/idx_bwa/$gx"
gsize="$genome/${gx}.chrom.sizes"
gbed="$genome/${gx}.gencode.main.bed.gz"
enhr="mm10.enhancers.mESC.bed"


##1. pre-processing
if [[ ! -e ${idx_grid}.bwt ]]; then
	bwa index -p $idx_grid "${idx_grid}.fa"
fi

for fq in gridseq.test*.raw.fq.gz; do 
	fqx=${fq/fq/trimmed.fq}

if [[ ! -e $fqx ]]; then
	echo $fqx
	touch $fqx

	cutadapt -j 11 -a AGATCGGA --trim-n -o $fqx $fq
fi

	lbam=${fq/fq.gz/GRIDv1.bam}
if [[ ! -e $lbam ]]; then
	echo $lbam
	touch $lbam

	bwa mem -t 11 -k 5 -L 4 -B 2 -O 4 $idx_grid $fqx | \
	samtools view -ub - | samtools sort -@ 11 -m 1G -o $lbam -
	samtools index -@ 11 $lbam
fi


##2. mapping to the genome
	prf=$(basename $lbam); prf=${prf/.raw.*bam}

	h5f="${prf}.h5"
	fqm="${prf}.mate.fq.gz"
	mid=${prf/gridseq./}

	echo $fqm
	python GridTools.py matefq -n $mid -o $h5f $lbam | gzip -c > $fqm


	bamx="${prf}.mate.${gx}.srt.bam"
	bamy=${bamx/srt/fix}
	bamz=${bamx/srt/mrk}

	echo $bamz

	bwa mem -t 31 -k 7 -U 0 -T 0 -ap $idx_mm38 $fqm | \
	samtools view -ubq 1 - | samtools sort -n -@ 11 -m 1G -o $bamx -
	samtools fixmate -@ 31 -pm $bamx - | samtools sort -@ 11 -m 1G -o $bamy -
	samtools markdup -@ 31 $bamy $bamz 
	samtools index -@ 11 $bamz


	bam=$bamz
	dbin=${bam/mrk.bam/dnabin.txt.gz}
	gexpr=${dbin/dnabin.txt.gz/gene_exprs.txt}
	gscop=${dbin/dnabin.txt.gz/gene_scope.txt}
	gmat=${dbin/dnabin/matrix}
	gnet=${dbin/dnabin/net}
	stat=${dbin/dnabin.txt.gz/stats}


	gridtools evaluate -k 1 -m 10 -g data/GRCm38.gencode.gtf.gz -o $h5f $bam
	gridtools stats -clbqr -p $stat $h5f

	gridtools DNA $h5f | gzip -c > $dbin
	gridtools RNA -e $gexpr -s $gscop $h5f

	gridtools matrix $h5f | gzip -c > $gmat
	gridtools model -e $elmt $h5f | gzip -c > $gnet

	
	echo "GRID-seq pipeline done."
done
