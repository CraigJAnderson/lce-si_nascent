#!/bin/bash
##setup.sh
##this bit generates the lengths of non-overlapping gene features to normalise read counts against.
WRKDIR=$1
ANNOTATION=$2
STRAIN=$3
GENE_LIST_PART_NUM=$4

cd ${WRKDIR}/src/${STRAIN}

while read line
 do bedtools intersect -a NOOVERLAPS_${ANNOTATION} -b exonsFULL.bed -wb | cut -f10- | X=${line} awk '{if ($4 == ENVIRON["X"]) print $0}' | sort -k1,1 -k2,2V | uniq | bedtools merge -i stdin | X=${line} awk '{sum += ($3-$2)} END {if (sum != "") {print ENVIRON["X"]"\t"sum} else {print ENVIRON["X"]"\t0"}}'
done < ./parts/gene_subset.${GENE_LIST_PART_NUM}.gff >> ./parts/exon_sizes_NO_OVERLAPS.${GENE_LIST_PART_NUM}.txt

##for non-overlapping introns
while read line
 do bedtools subtract -a NOOVERLAPS_${ANNOTATION} -b exonsFULL.bed | sed 's/;/\t/g' | sed 's/GeneID=//g' | sed 's/Name=//g' | X=${line} awk '{if ($9 == ENVIRON["X"]) print $1"\t"$4"\t"$5"\t"$7"\t"$9"\t"$10}' | sort -k1,1 -k2,2V | uniq | bedtools merge -i stdin | X=${line} awk '{if (($3-$2) > 0) sum += ($3-($2-1))} END {if (sum != "") {print ENVIRON["X"]"\t"sum} else {print ENVIRON["X"]"\t0"}}'
done < ./parts/gene_subset.${GENE_LIST_PART_NUM}.gff >> ./parts/intron_sizes_NO_OVERLAPS.${GENE_LIST_PART_NUM}.txt
