#/!/bin/sh
###cluster_bin_count.sh
##sample name is SAMP
##path to merged, STAR aligned reads files to look at as PATH e.g. /hps/nobackup2/flicek/user/mst/lce-si/analysis/staralign
SAMP=$1
WRKDIR=$2
PICARD=$3
STRAIN=$4

samtools view -h -f 83 -b ${WRKDIR}/data/${SAMP}.splice.merge.bam | samtools view -h -F 1284 | samtools sort > ${SAMP}.83.bam
samtools view -h -f 163 -b ${WRKDIR}/data/${SAMP}.splice.merge.bam | samtools view -h -F 1284 | samtools sort > ${SAMP}.163.bam
samtools merge ${SAMP}.83_163.bam ${SAMP}.163.bam ${SAMP}.83.bam && samtools sort ${SAMP}.83_163.bam > ${SAMP}.83_163.sorted.bam && samtools index ${SAMP}.83_163.sorted.bam

samtools view -h -f 99 -b ${WRKDIR}/data/${SAMP}.splice.merge.bam | samtools view -h -F 1284 | samtools sort > ${SAMP}.99.bam
samtools view -h -f 147 -b ${WRKDIR}/data/${SAMP}.splice.merge.bam | samtools view -h -F 1284 | samtools sort > ${SAMP}.147.bam
samtools merge ${SAMP}.99_147.bam ${SAMP}.147.bam ${SAMP}.99.bam && samtools sort ${SAMP}.99_147.bam > ${SAMP}.99_147.sorted.bam && samtools index ${SAMP}.99_147.sorted.bam

#take specific orientations of gene features and pull reads associated with them
awk '{if ($5 == "+") print $0}' ${WRKDIR}/src/${STRAIN}_all.exons.can.bed | bedtools intersect -split -f 1.0 -b stdin -abam ${SAMP}.83_163.sorted.bam > ${SAMP}_exons.83_163.tmpnontemplate.bam
awk '{if ($5 == "-") print $0}' ${WRKDIR}/src/${STRAIN}_all.exons.can.bed | bedtools intersect -split -f 1.0 -b stdin -abam ${SAMP}.83_163.sorted.bam > ${SAMP}_exons.83_163.tmptemplate.bam
awk '{if ($5 == "+") print $0}' ${WRKDIR}/src/${STRAIN}_all.exons.can.bed | bedtools intersect -split -f 1.0 -b stdin -abam ${SAMP}.99_147.sorted.bam > ${SAMP}_exons.99_147.tmptemplate.bam
awk '{if ($5 == "-") print $0}' ${WRKDIR}/src/${STRAIN}_all.exons.can.bed | bedtools intersect -split -f 1.0 -b stdin -abam ${SAMP}.99_147.sorted.bam > ${SAMP}_exons.99_147.tmpnontemplate.bam

samtools view ${SAMP}_exons.83_163.tmpnontemplate.bam | awk '{print $1}' | sort | uniq -c | awk -F" " '{if ($1 == 2) print $2}' > ${SAMP}_exons.83_163.tmpnontemplate.tag_names.txt
samtools view ${SAMP}_exons.83_163.tmptemplate.bam | awk '{print $1}' | sort | uniq -c | awk -F" " '{if ($1 == 2) print $2}' > ${SAMP}_exons.83_163.tmptemplate.tag_names.txt
samtools view ${SAMP}_exons.99_147.tmptemplate.bam | awk '{print $1}' | sort | uniq -c | awk -F" " '{if ($1 == 2) print $2}' > ${SAMP}_exons.99_147.tmptemplate.tag_names.txt
samtools view ${SAMP}_exons.99_147.tmpnontemplate.bam | awk '{print $1}' | sort | uniq -c | awk -F" " '{if ($1 == 2) print $2}' > ${SAMP}_exons.99_147.tmpnontemplate.tag_names.txt

java -Xmx50G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_exons.83_163.tmpnontemplate.bam FILTER=includeReadList RLF=${SAMP}_exons.83_163.tmpnontemplate.tag_names.txt O=${SAMP}_exons.83_163.nontemplate.bam MAX_RECORDS_IN_RAM=50000
java -Xmx50G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_exons.83_163.tmptemplate.bam FILTER=includeReadList RLF=${SAMP}_exons.83_163.tmptemplate.tag_names.txt O=${SAMP}_exons.83_163.template.bam MAX_RECORDS_IN_RAM=50000
java -Xmx50G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_exons.99_147.tmptemplate.bam FILTER=includeReadList RLF=${SAMP}_exons.99_147.tmptemplate.tag_names.txt O=${SAMP}_exons.99_147.template.bam MAX_RECORDS_IN_RAM=50000
java -Xmx50G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_exons.99_147.tmpnontemplate.bam FILTER=includeReadList RLF=${SAMP}_exons.99_147.tmpnontemplate.tag_names.txt O=${SAMP}_exons.99_147.nontemplate.bam MAX_RECORDS_IN_RAM=50000

##merge respective reads groups
samtools merge -sorted ${SAMP}_exons_nontemplate.bam ${SAMP}_exons.83_163.nontemplate.bam ${SAMP}_exons.99_147.nontemplate.bam ; samtools sort ${SAMP}_exons_nontemplate.bam > ${SAMP}_exons_nontemplate.sorted.bam ; samtools index ${SAMP}_exons_nontemplate.sorted.bam
samtools merge -sorted ${SAMP}_exons_template.bam ${SAMP}_exons.83_163.template.bam ${SAMP}_exons.99_147.template.bam ; samtools sort ${SAMP}_exons_template.bam > ${SAMP}_exons_template.sorted.bam ; samtools index ${SAMP}_exons_template.sorted.bam

##merge these to get read names. Because of how I've grouped exonic reads, any read pair mapping exclusively to exonic features is removed based upon this list. The remaining reads are intronic.
samtools merge -sorted ${SAMP}_all_exons.sorted.bam ${SAMP}_exons_nontemplate.sorted.bam ${SAMP}_exons_template.sorted.bam 
samtools view ${SAMP}_all_exons.sorted.bam | awk '{print $1}' > ${SAMP}_all_exons.tag_names.txt

#start getting intronic
awk '{ if (($9 == "protein_coding") && ($10 =="gene") && ($6 == "+")) print $1"\t"$2"\t"$3"\t"$7"\t"$6}' ${WRKDIR}/src/${STRAIN}.gtf.sorted.tsv  | bedtools intersect -b stdin -abam ${SAMP}.83_163.sorted.bam > ${SAMP}_introns.83_163.tmpnontemplate.bam
awk '{ if (($9 == "protein_coding") && ($10 =="gene") && ($6 == "-")) print $1"\t"$2"\t"$3"\t"$7"\t"$6}' ${WRKDIR}/src/${STRAIN}.gtf.sorted.tsv  | bedtools intersect -b stdin -abam ${SAMP}.83_163.sorted.bam > ${SAMP}_introns.83_163.tmptemplate.bam
awk '{ if (($9 == "protein_coding") && ($10 =="gene") && ($6 == "+")) print $1"\t"$2"\t"$3"\t"$7"\t"$6}' ${WRKDIR}/src/${STRAIN}.gtf.sorted.tsv  | bedtools intersect -b stdin -abam ${SAMP}.99_147.sorted.bam > ${SAMP}_introns.99_147.tmptemplate.bam
awk '{ if (($9 == "protein_coding") && ($10 =="gene") && ($6 == "-")) print $1"\t"$2"\t"$3"\t"$7"\t"$6}' ${WRKDIR}/src/${STRAIN}.gtf.sorted.tsv  | bedtools intersect -b stdin -abam ${SAMP}.99_147.sorted.bam > ${SAMP}_introns.99_147.tmpnontemplate.bam

java -Xmx24G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_introns.83_163.tmpnontemplate.bam FILTER=excludeReadList RLF=${SAMP}_all_exons.tag_names.txt O=${SAMP}_introns.83_163.tmp2nontemplate.bam MAX_RECORDS_IN_RAM=50000
java -Xmx24G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_introns.83_163.tmptemplate.bam FILTER=excludeReadList RLF=${SAMP}_all_exons.tag_names.txt O=${SAMP}_introns.83_163.tmp2template.bam MAX_RECORDS_IN_RAM=50000
java -Xmx24G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_introns.99_147.tmptemplate.bam FILTER=excludeReadList RLF=${SAMP}_all_exons.tag_names.txt O=${SAMP}_introns.99_147.tmp2template.bam MAX_RECORDS_IN_RAM=50000
java -Xmx24G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_introns.99_147.tmpnontemplate.bam FILTER=excludeReadList RLF=${SAMP}_all_exons.tag_names.txt O=${SAMP}_introns.99_147.tmp2nontemplate.bam MAX_RECORDS_IN_RAM=50000

samtools view ${SAMP}_introns.83_163.tmp2nontemplate.bam | awk '{print $1}' | sort | uniq -c | awk -F" " '{if ($1 == 2) print $2}' > ${SAMP}_introns.83_163.tmp2nontemplate.tag_names.txt
samtools view ${SAMP}_introns.83_163.tmp2template.bam | awk '{print $1}' | sort | uniq -c | awk -F" " '{if ($1 == 2) print $2}' > ${SAMP}_introns.83_163.tmp2template.tag_names.txt
samtools view ${SAMP}_introns.99_147.tmp2template.bam | awk '{print $1}' | sort | uniq -c | awk -F" " '{if ($1 == 2) print $2}' > ${SAMP}_introns.99_147.tmp2template.tag_names.txt
samtools view ${SAMP}_introns.99_147.tmp2nontemplate.bam | awk '{print $1}' | sort | uniq -c | awk -F" " '{if ($1 == 2) print $2}' > ${SAMP}_introns.99_147.tmp2nontemplate.tag_names.txt

java -Xmx24G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_introns.83_163.tmp2nontemplate.bam FILTER=includeReadList RLF=${SAMP}_introns.83_163.tmp2nontemplate.tag_names.txt O=${SAMP}_introns.83_163.nontemplate.bam MAX_RECORDS_IN_RAM=50000
java -Xmx24G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_introns.83_163.tmp2template.bam FILTER=includeReadList RLF=${SAMP}_introns.83_163.tmp2template.tag_names.txt O=${SAMP}_introns.83_163.template.bam MAX_RECORDS_IN_RAM=50000
java -Xmx24G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_introns.99_147.tmp2template.bam FILTER=includeReadList RLF=${SAMP}_introns.99_147.tmp2template.tag_names.txt O=${SAMP}_introns.99_147.template.bam MAX_RECORDS_IN_RAM=50000
java -Xmx24G -jar ${PICARD}/picard.jar FilterSamReads I=${SAMP}_introns.99_147.tmp2nontemplate.bam FILTER=includeReadList RLF=${SAMP}_introns.99_147.tmp2nontemplate.tag_names.txt O=${SAMP}_introns.99_147.nontemplate.bam MAX_RECORDS_IN_RAM=50000

#merge respective reads groups
samtools merge -sorted ${SAMP}_introns_nontemplate.bam ${SAMP}_introns.83_163.nontemplate.bam ${SAMP}_introns.99_147.nontemplate.bam ; samtools sort ${SAMP}_introns_nontemplate.bam > ${SAMP}_introns_nontemplate.sorted.bam ; samtools index ${SAMP}_introns_nontemplate.sorted.bam

samtools merge -sorted ${SAMP}_introns_template.bam ${SAMP}_introns.83_163.template.bam ${SAMP}_introns.99_147.template.bam ; samtools sort ${SAMP}_introns_template.bam > ${SAMP}_introns_template.sorted.bam ; samtools index ${SAMP}_introns_template.sorted.bam

rm ${SAMP}_introns_template.bam
rm ${SAMP}_introns.83_163.template.bam
rm ${SAMP}_introns.99_147.template.bam

rm ${SAMP}_introns_nontemplate.bam
rm ${SAMP}_introns.83_163.nontemplate.bam
rm ${SAMP}_introns.99_147.nontemplate.bam

rm ${SAMP}_exons_template.bam
rm ${SAMP}_exons.83_163.template.bam
rm ${SAMP}_exons.99_147.template.bam

rm ${SAMP}_exons_nontemplate.bam
rm ${SAMP}_exons.83_163.nontemplate.bam
rm ${SAMP}_exons.99_147.nontemplate.bam

rm ${SAMP}.163.bam
rm ${SAMP}.83.bam
rm ${SAMP}.83_163.bam
rm ${SAMP}.83_163.sorted.bam

rm ${SAMP}.99_147.bam
rm ${SAMP}.147.bam
rm ${SAMP}.99.bam
rm ${SAMP}.99_147.sorted.bam

rm ${SAMP}.83_163.sorted.bam.bai
rm ${SAMP}.99_147.sorted.bam.bai

rm ${SAMP}_introns.83_163.tmpnontemplate.bam
rm ${SAMP}_introns.83_163.tmptemplate.bam
rm ${SAMP}_introns.99_147.tmptemplate.bam
rm ${SAMP}_introns.99_147.tmpnontemplate.bam

