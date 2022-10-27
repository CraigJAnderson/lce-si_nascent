#!/bin/bash
#tpm_module.sh
##example requirements: bsub -M 20000 -R "rusage[mem=20000]" -W 12:00
WRKDIR=$1
SAMPLE=$2
STRAIN=$3
cd ${WRKDIR}/analysis
#now count reads across each of the genes- sort first
for y in exons_template exons_nontemplate introns_template introns_nontemplate ; do samtools sort ${SAMPLE}_${y}.mate_sorted.bam > ${SAMPLE}_${y}.mate_sorted2.bam ; samtools index ${SAMPLE}_${y}.mate_sorted2.bam ; done
#now use bedtools to count reads across features.
for y in exons_template exons_nontemplate introns_template introns_nontemplate ; do while read line ; do grep -w -F ${line} ${WRKDIR}/src/${STRAIN}/${STRAIN}.gtf.sorted.tsv  | awk '{ if ($10 == "gene") print $0}' | bedtools multicov -bams ${SAMPLE}_${y}.mate_sorted2.bam -p -D -q 10 -bed stdin | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$11}' ; done < ${WRKDIR}/src/${STRAIN}/geneid_NOVERLAP.csv >> ../output/${SAMPLE}_${y}.bed.counts ; done
cd ${WRKDIR}/output
for y in exons_template exons_nontemplate ; do ${WRKDIR}/bin/tpm_calc2.sh ${SAMPLE} ${y} ${WRKDIR}/src/${STRAIN}/parts/exon_sizes.txt ; done
for y in introns_template introns_nontemplate ; do ${WRKDIR}/bin/tpm_calc2.sh ${SAMPLE} ${y} ${WRKDIR}/src/${STRAIN}/parts/intron_sizes.txt ; done
##some genes overlap so extensively with others that they have no data. I can't make a fair comparison if the search space is zero in any set, therefore they should be removed. If search space is more than zero for both introns and exons, then they are considered zero.

