#/!/bin/sh
NOD=$1
READ_TYPE=$2
FEAT_SIZES=$3

#calculate RPKn for each gene, where RPKn = READ_COUNTSn / (GENE_LENGTHn/1000), where n is a single gene
while read line ; do bits=($line) ; LENGTH=$(X=${bits[3]} awk '{if ($1 == ENVIRON["X"]) print $0}' ${FEAT_SIZES}) ; echo $line $LENGTH | awk -F" " '{if (($7 > 0) && ($9 > 0)) {printf("%s %s:%s-%s %s %s %s %s %.6f\n", $4,$1,$2,$3,$5,$6,$7,$9,($7/($9/1000)))} else {printf("%s %s:%s-%s %s %s %s %s 0\n", $4,$1,$2,$3,$5,$6,$7,$9)}}' | sed 's/ /\t/g' >> ${NOD}_${READ_TYPE}.rpk ; done < ${NOD}_${READ_TYPE}.bed.counts

#calculate scaling factor as scaling_factor = sumRPK/1000000 (from whole gene analysis)
SCALING_FACTOR=$(awk '{SUM+=$7}END{printf("%.6f\n", SUM/1000000)}' ${NOD}_${READ_TYPE}.rpk)
echo ${NOD} ${READ_TYPE} ${SCALING_FACTOR} > ${NOD}_${READ_TYPE}.scaling_factor

#TPMn = RPKn/scaling_factor
X=$SCALING_FACTOR awk '{ if ($7 == 0) {printf("%s %.6f\n", $0,$7)} else {printf("%s %.6f\n", $0,$7/ENVIRON["X"])}}' ${NOD}_${READ_TYPE}.rpk | sed 's/ /\t/g' > ${NOD}_${READ_TYPE}.tpm
rm ${NOD}_${READ_TYPE}.rpk
