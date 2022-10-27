Measuring transcription with measures of mRNA give insight into the steady state of transcripts in a cell, though measuring pre-mRNA we can gain better insight into the relative frequency of RNAPII passing across gene bodies. This measure of nascent RNA therefore allows us to more accurately quantify where repair via TC-NER is more likely to be triggered. There are lots of protocols for accurate measures of nascent transcription, though in bulk RNAseq exist pre-mRNA reads, normally filtered out of analyses becasue they span introns and exons. We select for such reads as a measure of nscent transcription.

this pipeline will align RNAseq data accross genome annotation. I subsequently determine non-overlapping annotations with which to count reads across in order to calculate TPM for nascent expression. We then calculate mutation rate for those same feature regions and compare mutation rate with expression.

Please cite our work if you use this pipeline: https://www.biorxiv.org/content/10.1101/2022.06.10.495644v1.full

The assumption is that you've downloaded the github repo and create a folder in src that is the name of the strain of mouse (or species, whatever) you're analysing, that you then populate with gene models, exons lengths, etc. You'll then populate the data folder with your processed BAMs.

I've created the annotation data for c3h, with my own sample names. Decompress this to see what there is before you get started. You'll need to build your own STAR index and remember that I've focussed on autosomes and the X, and have created these so there is no "chr" in the chromosome names. 

Output is mean TPM for the samples given, though TPM for each sample is output, so you can process the results as you see fit.

##############################

Making exon boundaries readme

Generate exon boundaries- get them from biomart most easily- see the mini_exons.png image in .pics for the recipe.

Here's a link to an ensembl 91 version of bl6

http://dec2017.archive.ensembl.org/biomart/martview/0d95712e13769fb29a86fd5c28c6a789?VIRTUALSCHEMANAME=default&ATTRIBUTES=mc57bl6nj_gene_ensembl.default.structure.chromosome_name|mc57bl6nj_gene_ensembl.default.structure.exon_chrom_start|mc57bl6nj_gene_ensembl.default.structure.exon_chrom_end|mc57bl6nj_gene_ensembl.default.structure.ensembl_transcript_id|mc57bl6nj_gene_ensembl.default.structure.ensembl_gene_id|mc57bl6nj_gene_ensembl.default.structure.strand|mc57bl6nj_gene_ensembl.default.structure.external_gene_name|mc57bl6nj_gene_ensembl.default.structure.gene_biotype&FILTERS=&VISIBLEPANEL=resultspanel

#save to output, examplified here: head ${EXONS_ANNO}
#1       256705  259441  MGP_C3HHeJ_G0015578;MGP_C3HHeJ_T0018487 1       -       Xkr4    protein_coding  exon
#I'm aware this is a c3h gene

#you will likely have to remove a header and add the exons column at the end- just use sed.

#####################################

Making non-overlapping exons readme

example biomart input- see recipe in .pics folder under mini_gene_anno.png.

#head ${BIOMART_ANNO}
#Gene stable ID	Gene name	Protein stable ID	Chromosome/scaffold nameGenomic coding start	Genomic coding end	CDS start	CDS end	CDS Length	Strand
#MGP_C3HHeJ_G0029249	Zcchc4	MGP_C3HHeJ_P0071329	5	51821250	51821373	1	124	1536	1
#I'm aware this is a c3h gene
##URL as:
http://dec2017.archive.ensembl.org/biomart/martview/ebfed90c012aa08f8b672710b402c3ce?VIRTUALSCHEMANAME=default&ATTRIBUTES=mc57bl6nj_gene_ensembl.default.structure.ensembl_gene_id|mc57bl6nj_gene_ensembl.default.structure.external_gene_name|mc57bl6nj_gene_ensembl.default.structure.ensembl_peptide_id|mc57bl6nj_gene_ensembl.default.structure.chromosome_name|mc57bl6nj_gene_ensembl.default.structure.genomic_coding_start|mc57bl6nj_gene_ensembl.default.structure.genomic_coding_end|mc57bl6nj_gene_ensembl.default.structure.cds_start|mc57bl6nj_gene_ensembl.default.structure.cds_end|mc57bl6nj_gene_ensembl.default.structure.cds_length|mc57bl6nj_gene_ensembl.default.structure.strand&FILTERS=&VISIBLEPANEL=resultspanel

##remove header, save as ${BIOMART_ANNO}


##############################

Make chromosome sizes readme

Two columns, tab delimited. First is chromosome name, second is nucleotide length of chromosome. Save as ${CHROM_SZ}

#############################

Make list of protein coding genes readme:

Example for hg38

wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz

zcat ${WRKDIR}/src/hg19/Homo_sapiens.GRCh38.108.gtf.gz | grep -v '^#' | perl -w ${WRKDIR}/bin/gtfToTsv.pl | awk '{ if (($10 == "gene") && ($9 == "protein_coding")) print $1"\tcurated\tgene\t"$2"\t"$3"\t.\t"$6"\t"$5"\tGeneID="$4";Name="$7}' > ${WRKDIR}/src/hg19/hg38_protein_coding_whole_genes.gff


#############################

Make STAR index folder readme:

STAR --runThreadN 12 --runMode genomeGenerate --sjdbGTFfile bl6.gtf --genomeDir bl6.starIndex --outFileNamePrefix bl6.star. --genomeFastaFiles bl6.fa


