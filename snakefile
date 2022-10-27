import os
import pandas as pd

configfile : "config.yaml"

os.chdir(config["WRKDIR"])

if not os.path.exists("output"):
 os.mkdir("output")

if not os.path.exists("log"):
 os.mkdir("log")

if not os.path.exists("analysis"):
 os.mkdir("analysis")

if not os.path.exists("src/"+config["STRAIN"]):
 os.makedirs("src/"+config["STRAIN"])

sample_names = pd.read_table(config["WRKDIR"]+"/src/"+config["STRAIN"]+"/"+config["SAMPLE_LIST"],sep=" ",header=None)
sample_names.columns = ['nod']
SAMPLES= list(sample_names.nod)
num_samples = len(sample_names)

##this variable is a list of part numbers used to speedup the setup_pt2
GENE_LIST_PARTS=["%.3d" % i for i in range(1,101)]

rule all:
 input:
  h=config["WRKDIR"]+"/src/"+config["STRAIN"]+"/exonsFULL.bed",
  exon_sizes=expand(config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/exon_sizes_NO_OVERLAPS.{gene_list_part_num}.txt", gene_list_part_num=GENE_LIST_PARTS),
  intron_sizes=expand(config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/intron_sizes_NO_OVERLAPS.{gene_list_part_num}.txt", gene_list_part_num=GENE_LIST_PARTS),
  out=config["WRKDIR"]+"/output/"+config["STRAIN"]+".mean_introns_template.tpm",
  g=expand(config["WRKDIR"]+"/data/{sample}.splice.merge.bam", sample=SAMPLES),
  i=config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/intron_sizes.txt",
  j=config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/exon_sizes.txt",
  k=expand(config["WRKDIR"]+"/analysis/{sample}_introns_nontemplate.mate_sorted.bam", sample=SAMPLES),
  l=expand(config["WRKDIR"]+"/analysis/{sample}_exons_nontemplate.mate_sorted.bam", sample=SAMPLES),
  m=expand(config["WRKDIR"]+"/analysis/{sample}_introns_template.mate_sorted.bam", sample=SAMPLES),
  n=expand(config["WRKDIR"]+"/analysis/{sample}_exons_template.mate_sorted.bam", sample=SAMPLES),

ruleorder: setup_pt1 > setup_pt2 > setup_pt3 > count_module > tpm_module > mean_tpm_module

rule setup_pt1:
 input:
  exons_anno = config["WRKDIR"]+"/src/"+config["STRAIN"]+"/"+config["EXONS_ANNO"],
  tmp = config["WRKDIR"]+"/src/"+config["STRAIN"]+"/"+config["ANNOTATION"],
  chrom_sz = config["WRKDIR"]+"/src/"+config["STRAIN"]+"/"+config["CHROM_SZ"],
  biomart_anno = config["WRKDIR"]+"/src/"+config["STRAIN"]+"/"+config["BIOMART_ANNO"],
  annotation_pth = config["WRKDIR"]+"/src/"+config["STRAIN"]+"/"+config["ANNOTATION"],
 params:
  annotation = config["ANNOTATION"],
  wrkdir = config["WRKDIR"],
  strain = config["STRAIN"]
 output:
  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/exonsFULL.bed",
  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/parts.txt",
  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/OVERLAPS_"+config["ANNOTATION"]
 shell:
  """ {params.wrkdir}/bin/setup_pt1.sh {params.wrkdir} {input.exons_anno} {params.annotation} {input.chrom_sz} {input.biomart_anno} {params.strain} {input.annotation_pth} """

rule setup_pt2:
 input:
  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/OVERLAPS_"+config["ANNOTATION"],
 params:
  wrkdir = config["WRKDIR"],
  annotation = config["ANNOTATION"],
  strain = config["STRAIN"]
 resources:
  mem_mb=4000,
  resources='"rusage[mem=4000]"',
  time="0:20"
 output:
  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/exon_sizes_NO_OVERLAPS.{gene_list_part_num}.txt",
  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/intron_sizes_NO_OVERLAPS.{gene_list_part_num}.txt"
 shell:
  """ {params.wrkdir}/bin/setup_pt2.sh {params.wrkdir} {params.annotation} {params.strain} {wildcards.gene_list_part_num} """

rule setup_pt3:
 input:
  expand(config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/exon_sizes_NO_OVERLAPS.{gene_list_part_num}.txt", gene_list_part_num=GENE_LIST_PARTS),
  expand(config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/intron_sizes_NO_OVERLAPS.{gene_list_part_num}.txt", gene_list_part_num=GENE_LIST_PARTS)
 params:
  wrkdir = config["WRKDIR"],
  strain = config["STRAIN"]
 output:
  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/intron_sizes.txt",
  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/exon_sizes.txt"
 shell:
  """ {params.wrkdir}/bin/setup_pt3.sh {params.wrkdir} {params.strain} """

#rule align_module:
# input: 
#  config["WRKDIR"]+"/data/{sample}.p1.fq.gz", 
#  config["WRKDIR"]+"/data/{sample}.p2.fq.gz",
#  config["WRKDIR"]+"/data/{sample}.raw.Aligned.sortedByCoord.out.bam",
#  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/intron_sizes.txt",
#  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/exon_sizes.txt"
# params:
#  fastq_loc = config["WRKDIR"]+"/data",
#  star_index = config["WRKDIR"]+"/src/"+config["STRAIN"]+"/"+config["STAR_INDEX"],
#  wrkdir = config["WRKDIR"]
# output: 
#  config["WRKDIR"]+"/data/{sample}.splice.merge.bam"
# shell:
#  """ {params.wrkdir}/bin/star_align.sh {wildcards.sample} {params.fastq_loc} {params.star_index} {params.wrkdir} """

rule count_module: ##calling java directly in script- set at 24G, needs to be high but change accordingly.
 input:
  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/intron_sizes.txt",
  config["WRKDIR"]+"/src/"+config["STRAIN"]+"/parts/exon_sizes.txt",
  config["WRKDIR"]+"/data/{sample}.splice.merge.bam",
 params:
  wrkdir = config["WRKDIR"],
  picard_loc = config["PICARD_LOC"], 
  strain = config["STRAIN"]
 resources:
  mem_mb=24000,
  resources='"rusage[mem=24000]"',
  time="12:00"
 output:
  config["WRKDIR"]+"/analysis/{sample}_introns_nontemplate.mate_sorted.bam",
  config["WRKDIR"]+"/analysis/{sample}_exons_nontemplate.mate_sorted.bam",
  config["WRKDIR"]+"/analysis/{sample}_introns_template.mate_sorted.bam",
  config["WRKDIR"]+"/analysis/{sample}_exons_template.mate_sorted.bam"
 shell:
  """ {params.wrkdir}/bin/count_module.sh {params.wrkdir} {params.picard_loc} {wildcards.sample} {params.strain}"""

rule tpm_module:
 input: 
  config["WRKDIR"]+"/analysis/{sample}_introns_nontemplate.mate_sorted.bam",
  config["WRKDIR"]+"/analysis/{sample}_exons_nontemplate.mate_sorted.bam",
  config["WRKDIR"]+"/analysis/{sample}_introns_template.mate_sorted.bam",
  config["WRKDIR"]+"/analysis/{sample}_exons_template.mate_sorted.bam" 
 params:
  wrkdir = config["WRKDIR"],
  strain = config["STRAIN"],
 output: 
  config["WRKDIR"]+"/output/{sample}_introns_nontemplate.bed.counts",
  config["WRKDIR"]+"/output/{sample}_introns_template.bed.counts"
 shell:
  """ {params.wrkdir}/bin/tpm_module.sh {params.wrkdir} {wildcards.sample} {params.strain} """

rule mean_tpm_module:
 input:
  expand(config["WRKDIR"]+"/output/{sample}_introns_nontemplate.bed.counts", sample=SAMPLES),
  expand(config["WRKDIR"]+"/output/{sample}_introns_template.bed.counts", sample=SAMPLES)
 params:
  wrkdir = config["WRKDIR"],
  strain = config["STRAIN"],
  sample_list = config["WRKDIR"]+"/src/"+config["STRAIN"]+"/"+config["SAMPLE_LIST"],
  NUM_SAM = num_samples
 output: 
  config["WRKDIR"]+"/output/"+config["STRAIN"]+".mean_exons_nontemplate.tpm",
  config["WRKDIR"]+"/output/"+config["STRAIN"]+".mean_introns_nontemplate.tpm",
  config["WRKDIR"]+"/output/"+config["STRAIN"]+".mean_exons_template.tpm",
  config["WRKDIR"]+"/output/"+config["STRAIN"]+".mean_introns_template.tpm"
 shell:
  """ {params.wrkdir}/bin/mean_tpm_module.sh {params.wrkdir} {params.NUM_SAM} {params.strain} {params.sample_list}""" ##num_samples is generated from sample list
