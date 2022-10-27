#!/bin/bash
##setup.sh
##this bit generates the lengths of non-overlapping gene features to normalise read counts against.
WRKDIR=$1
EXONS_ANNO=$2
ANNOTATION=$3
CHROM_SZ=$4
BIOMART_ANNO=$5
STRAIN=$5

cd ${WRKDIR}
mkdir src
mkdir log
mkdir output
mkdir analysis

cd src
mkdir ${STRAIN}
cd ${STRAIN}

#see making exon boundaries readme
awk '{ print $1"\t"$2"\t"$3"\t"$4";"$5"\t"$6"\t-\t"$7"\t"$8"\texon"}' ${EXONS_ANNO} > exons.bed

##need to generate non-overlapping features for counting over, first need gtf with no overlaps- get parts of features with greater than 1 coverage
bedtools genomecov -i ${ANNOTATION} -g ${CHROM_SZ} -bg | awk '{ if ($4 >1) print $0}' > OVERLAPS_${ANNOTATION}
#now subtract any regions where $4 > 1
bedtools subtract -a ${ANNOTATION} -b OVERLAPS_${ANNOTATION} > NOOVERLAPS_${ANNOTATION}

#generate exon and intron feature lengths for each gene id, to normalise against. run in parts of 100 because it'll take a long time.
#index based upon number of genes in gtf
sed 's/;/\t/g' NOOVERLAPS_${ANNOTATION} | sed 's/GeneID=//g' | sed 's/Name=//g' | awk '{print $9}' | uniq > geneid_NOVERLAP.csv
mkdir parts
NUM=$(wc -l geneid_NOVERLAP.csv |cut -f1 -d " ")
NUM1=$(( ((${NUM})/100)*100 ))
NUM2=$(( $NUM - ${NUM1}))
NUM3=$(($(( $NUM1 / 100))))
paste <(echo $(seq 100 100 ${NUM1}) ${NUM}| sed 's/ /\n/g') <(echo $(for x in $(seq 1 1 ${NUM3}) ; do echo 100 ; done | xargs) ${NUM2} | sed 's/ /\n/g') > ./parts/parts.txt
for x in $(seq 1 1 ${NUM}) ; do bits=$(sed -n ${x}p ./parts/parts.txt) ; line=($bits) ; head -n ${line[0]} geneid_NOVERLAP.csv | tail -n ${line[1]} > ./${strain}_parts/gene_subset.${x}.gff ; done

#see making exon boundaries readme
#fix biomart and get gene names
cat {BIOMART_ANNO} | awk '{ if ( NF == 10 ) {print $0} else { if (NF == 9) print $1"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}' | cut -f1 | sort | uniq > canonical_transcript_names.txt
grep -w -F -f canonical_transcript_names.txt exons.bed | sed 's/;/\t/g' | awk '{ print $0}' > exonsFULL.bed

#while read line ; do bedtools intersect -a NOOVERLAPS_{ANNOTATION} -b exonsFULL.bed -wb | cut -f10- | X=${line} awk '{if ($4 == ENVIRON["X"]) print $0}' | sort -k1,1 -k2,2V | uniq | bedtools merge -i stdin | X=${line} awk '{sum += ($3-$2)} END {if (sum != "") {print ENVIRON["X"]"\t"sum} else {print ENVIRON["X"]"\t0"}}' ; done < ./parts/gene_subset.${x}.gff >> ./parts/exon_sizes_NO_OVERLAPS.${x}.txt

for x in $(seq 1 1 ${NUM}) ; do while read line ; do bedtools intersect -a NOOVERLAPS_{ANNOTATION} -b exonsFULL.bed -wb | cut -f10- | X=${line} awk '{if ($4 == ENVIRON["X"]) print $0}' | sort -k1,1 -k2,2V | uniq | bedtools merge -i stdin | X=${line} awk '{sum += ($3-$2)} END {if (sum != "") {print ENVIRON["X"]"\t"sum} else {print ENVIRON["X"]"\t0"}}' ; done < ./parts/gene_subset.${x}.gff >> ./parts/exon_sizes_NO_OVERLAPS.${x}.txt ; done 

##for non-overlapping introns
for x in $(seq 1 1 ${NUM}) ; do while read line ; do bedtools subtract -a NOOVERLAPS_{ANNOTATION} -b exonsFULL.bed | sed 's/;/\t/g' | sed 's/GeneID=//g' | sed 's/Name=//g' | X=${line} awk '{if ($9 == ENVIRON["X"]) print $1"\t"$4"\t"$5"\t"$7"\t"$9"\t"$10}' | sort -k1,1 -k2,2V | uniq | bedtools merge -i stdin | X=${line} awk '{if (($3-$2) > 0) sum += ($3-($2-1))} END {if (sum != "") {print ENVIRON["X"]"\t"sum} else {print ENVIRON["X"]"\t0"}}' ; done < ./parts/gene_subset.${x}.gff >> ./parts/intron_sizes_NO_OVERLAPS.${x}.txt ; done

#put together
cat $(for x in $(seq 1 1 ${NUM}) ; do echo ./parts/intron_sizes_NO_OVERLAPS.${x}.txt ; done | xargs) > ./parts/intron_sizes.txt
cat $(for x in $(seq 1 1 ${NUM}) ; do echo ./parts/exon_sizes_NO_OVERLAPS.${x}.txt ; done | xargs) > ./parts/exon_sizes.txt

