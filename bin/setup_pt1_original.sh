#!/bin/bash
##setup.sh
##this bit generates the lengths of non-overlapping gene features to normalise read counts against.
WRKDIR=$1
EXONS_ANNO=$2
ANNOTATION_PTH=$3
ANNOTATION=$4
CHROM_SZ=$5
BIOMART_ANNO=$6
STRAIN=$7

cd ${WRKDIR}
cd src
cd ${STRAIN}

#see making exon boundaries readme
awk '{ print $1"\t"$2"\t"$3"\t"$4";"$5"\t"$6"\t-\t"$7"\t"$8"\texon"}' ${EXONS_ANNO} > exons.bed

##need to generate non-overlapping features for counting over, first need gtf with no overlaps- get parts of features with greater than 1 coverage
bedtools genomecov -i ${ANNOTATION_PTH}/${ANNOTATION} -g ${CHROM_SZ} -bg | awk '{ if ($4 >1) print $0}' > OVERLAPS_${ANNOTATION}
#now subtract any regions where $4 > 1
bedtools subtract -a ${ANNOTATION_PTH}/${ANNOTATION} -b OVERLAPS_${ANNOTATION} > NOOVERLAPS_${ANNOTATION}

#generate exon and intron feature lengths for each gene id, to normalise against. run in parts of 100 because it'll take a long time.
#index based upon number of genes in gtf
sed 's/;/\t/g' NOOVERLAPS_${ANNOTATION} | sed 's/GeneID=//g' | sed 's/Name=//g' | awk '{print $9}' | uniq > geneid_NOVERLAP.csv
mkdir parts
NUM=$(wc -l geneid_NOVERLAP.csv |cut -f1 -d " ")
NUM1=$(( ((${NUM})/100)*100 ))
NUM2=$(( $NUM - ${NUM1}))
NUM3=$(($(( $NUM1 / 100))))
paste <(echo $(seq 100 100 ${NUM1}) ${NUM}| sed 's/ /\n/g') <(echo $(for x in $(seq 1 1 ${NUM3}) ; do echo 100 ; done | xargs) ${NUM2} | sed 's/ /\n/g') > ./parts/parts.txt
for x in $(seq 1 1 $(wc -l ./parts/parts.txt | cut -f1 -d " ")) ; do bits=$(sed -n ${x}p ./parts/parts.txt) ; line=(${bits}) ; head -n ${line[0]} geneid_NOVERLAP.csv | tail -n ${line[1]} > ./parts/gene_subset.${x}.gff ; done

#see making exon boundaries readme
#fix biomart and get gene names
cat ${BIOMART_ANNO} | awk '{ if ( NF == 10 ) {print $0} else { if (NF == 9) print $1"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}' | cut -f1 | sort | uniq > canonical_transcript_names.txt
grep -w -F -f canonical_transcript_names.txt exons.bed | sed 's/;/\t/g' | awk '{ print $0}' > exonsFULL.bed

