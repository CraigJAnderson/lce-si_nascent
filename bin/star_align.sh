#!/bin/bash
#star_align.sh
SAMPLE=$1
FASTQ=$2
STAR_INDEX=$3
$WRKDIR=$4
#alignment module : STAR align and process trimmed reads
cd ${WRKDIR}/data

##the following is how I would align paired RNAseq data. It can be complicated to preempt the experimental setup, so the starting point for this pipeline is the last line here
STAR --runThreadN 6 --runMode alignReads --genomeDir ${STAR_INDEX} --chimOutType WithinBAM --readFilesCommand zcat --readFilesIn ${FASTQ}/${SAMPLE}.p1.fq.gz ${FASTQ}/${SAMPLE}.p2.fq.gz --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --outFileNamePrefix ${SAMPLE}.raw.
#change name
ln -s ${SAMPLE}.raw.Aligned.sortedByCoord.out.bam ${SAMPLE}.raw.bam
samtools sort -o ${SAMPLE}.sort.bam -@ 6 -O bam ${SAMPLE}.raw.Aligned.sortedByCoord.out.bam
picard MarkDuplicates I=${SAMPLE}.sort.bam O=${SAMPLE}.splice.merge.bam M=${SAMPLE}.mark_dup_metrics REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=${WRKDIR}/tmp

