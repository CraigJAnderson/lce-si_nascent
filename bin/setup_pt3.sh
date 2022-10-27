#!/bin/bash
##setup.sh
##this bit generates the lengths of non-overlapping gene features to normalise read counts against.
WRKDIR=$1
STRAIN=$2

cd ${WRKDIR}/src/${STRAIN}
#put together
cat $(for x in {001..100} ; do echo ./parts/intron_sizes_NO_OVERLAPS.${x}.txt ; done | xargs) > ./parts/intron_sizes.txt
cat $(for x in {001..100} ; do echo ./parts/exon_sizes_NO_OVERLAPS.${x}.txt ; done | xargs) > ./parts/exon_sizes.txt

paste parts/exon_sizes.txt parts/intron_sizes.txt | awk '{if (($2 == 0) || ($4 == 0)) print $1}' > parts/gene_list_nosize.txt
