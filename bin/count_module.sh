#!/bin/bash
##count_module.sh
# orientate reads with mates and count
#bsub -M 55000 -R "rusage[mem=55000]" -W 45:00  ### original /hps/nobackup2/flicek/user/cander21/strandAsym/premRNA_abundance/
WRKDIR=$1
PICARD_LOC=$2
SAMPLE=$3
STRAIN=$4

cd ${WRKDIR}/analysis

if [ ! -d "tmp" ]; then
  mkdir tmp
fi

${WRKDIR}/bin/separator5.sh ${SAMPLE} ${WRKDIR} ${PICARD_LOC} ${STRAIN}
for x in introns_template introns_nontemplate exons_template exons_nontemplate
 do java -Xmx24G -jar ${PICARD_LOC}/picard.jar FixMateInformation I=${SAMPLE}_${x}.sorted.bam O=${SAMPLE}_${x}.mate_sorted.bam TMP_DIR=tmp SO=queryname
done


#general counts with htseq if needed to compare. Prefer my bedtools version as more nuanced control and I know exactly what I'm getting.
#for x in introns_template introns_nontemplate exons_template exons_nontemplate
# do htseq-count ${SAMPLE}_${x}.mate_sorted.bam ${ANNOTATION_PTH}/${ANNOTATION} --idattr=GeneID -t gene --stranded no -f bam --nonunique all > ../output/${SAMPLE}_${x}.counts
#done

