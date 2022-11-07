#!/bin/bash
##setup.sh
##this bit generates the lengths of non-overlapping gene features to normalise read counts against.
WRKDIR=$1
EXONS_ANNO=$2
ANNOTATION=$3
CHROM_SZ=$4
BIOMART_ANNO=$5
STRAIN=$6
ANNOTATION_PTH=$7

cd ${WRKDIR}/src/${STRAIN}

#see making exon boundaries readme
#sed 's/\t1\t/\t+\t/g' ${EXONS_ANNO} | sed 's/\t-1\t/\t-\t/g' | awk '{ print $1"\t"$2"\t"$3"\t"$4";"$5";1\t"$6"\t-\t"$7"\t"$8"\texon"}' > exons.bed

sed 's/\t1\t/\t+\t/g' ${EXONS_ANNO} | sed 's/\t-1\t/\t-\t/g' | awk '{ print $1"\t"$2"\t"$3"\t"$4";"$5"\t"$6"\t"$7"\t"$8"\texon"}' > exons.bed

##need to generate non-overlapping features for counting over, first need gtf with no overlaps- get parts of features with greater than 1 coverage
bedtools genomecov -i ${ANNOTATION_PTH} -g ${CHROM_SZ} -bg | awk '{ if ($4 >1) print $0}' > OVERLAPS_${ANNOTATION}
#now subtract any regions where $4 > 1
bedtools subtract -a ${ANNOTATION_PTH} -b OVERLAPS_${ANNOTATION} > NOOVERLAPS_${ANNOTATION}

#generate exon and intron feature lengths for each gene id, to normalise against. run in parts of 100 because it'll take a long time.
#index based upon number of genes in gtf
sed 's/;/\t/g' NOOVERLAPS_${ANNOTATION} | sed 's/GeneID=//g' | sed 's/Name=//g' | awk '{print $9}' | uniq > geneid_NOVERLAP.csv

if [ ! -d "parts" ]; then
  mkdir parts
fi

#the next few lines making parts.txt split the gene list into 100 parts. The old version wasn't snakemake friendly, so this is more stable and will always be 100. As such, the minimum number of genes to include must be 100.
split -n l/100 geneid_NOVERLAP.csv -d ./parts/parts. --additional-suffix=".txt" --numeric-suffixes=1
for x in {001..100} ; do VAR=$(wc -l ./parts/parts.${x}.txt | cut -d" " -f1 ) ; echo $VAR >> ./parts/tmp_parts.txt ; rm ./parts/parts.${x}.txt ; done
awk '{$1=c+=$1}1' ./parts/tmp_parts.txt | paste - ./parts/tmp_parts.txt > ./parts/parts.txt

for x in {001..100} ; do bits=$(sed -n ${x}p ./parts/parts.txt) ; line=(${bits}) ; head -n ${line[0]} geneid_NOVERLAP.csv | tail -n ${line[1]} > ./parts/gene_subset.${x}.gff ; done

#see making exon boundaries readme
#fix biomart and get gene names
cat ${BIOMART_ANNO} | awk '{ if ( NF == 10 ) {print $0} else { if (NF == 9) print $1"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}' | cut -f1 | sort | uniq > canonical_transcript_names.txt
grep -w -F -f canonical_transcript_names.txt exons.bed | sed 's/;/\t/g' | awk '{ print $0}' > exonsFULL.bed

