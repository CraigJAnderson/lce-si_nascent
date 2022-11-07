#!/bin/bash
#mean_tpm_module.sh
WRKDIR=$1
NUM_SAMPLE=$2
STRAIN=$3
SAMPLE_LIST=$4
cd ${WRKDIR}/output
SEQ_NUM_SAMPLE=$(seq 8 8 $(( 8 * ${NUM_SAMPLE} )))
#the following will put together a list of all the tpm scores for all genes across all samples and mean, by which the output is sorted from highest mean expression to lowest.
paste $(cat ${SAMPLE_LIST} | sed 's/$/_exons_nontemplate.tpm/g' | xargs) | cut -f1-4,$(eval echo $SEQ_NUM_SAMPLE |  tr " " ",") | grep -v -F -w -f ${WRKDIR}/src/${STRAIN}/parts/gene_list_nosize.txt | awk '{ s = 0; for (i = 5; i <= NF; i++) s += $i; print $1,$2,$3,$4,s/(NF-4) }' > ${STRAIN}.mean_exons_nontemplate.tpm
paste $(cat ${SAMPLE_LIST} | sed 's/$/_introns_nontemplate.tpm/g'| xargs) | cut -f1-4,$(eval echo $SEQ_NUM_SAMPLE |  tr " " ",") | grep -v -F -w -f ${WRKDIR}/src/${STRAIN}/parts/gene_list_nosize.txt | awk '{ s = 0; for (i = 5; i <= NF; i++) s += $i; print $1,$2,$3,$4,s/(NF-4) }' > ${STRAIN}.mean_introns_nontemplate.tpm
paste $(cat ${SAMPLE_LIST} | sed 's/$/_exons_template.tpm/g' | xargs) | cut -f1-4,$(eval echo $SEQ_NUM_SAMPLE |  tr " " ",") | grep -v -F -w -f ${WRKDIR}/src/${STRAIN}/parts/gene_list_nosize.txt | awk '{ s = 0; for (i = 5; i <= NF; i++) s += $i; print $1,$2,$3,$4,s/(NF-4) }' > ${STRAIN}.mean_exons_template.tpm
paste $(cat ${SAMPLE_LIST} | sed 's/$/_introns_template.tpm/g' | xargs) | cut -f1-4,$(eval echo $SEQ_NUM_SAMPLE |  tr " " ",") | grep -v -F -w -f ${WRKDIR}/src/${STRAIN}/parts/gene_list_nosize.txt | awk '{ s = 0; for (i = 5; i <= NF; i++) s += $i; print $1,$2,$3,$4,s/(NF-4) }' > ${STRAIN}.mean_introns_template.tpm
