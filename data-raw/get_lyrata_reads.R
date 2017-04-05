library(SRAdb)

sqlfile <- "~/ncbi/sradb/SRAmetadb.sqlite"
if(!file.exists("~/ncbi/sradb/SRAmetadb.sqlite")) sqlfile <<- getSRAdbFile(destdir = "~/ncbi/sradb")

file.info(sqlfile)
# sra_tables <- dbListTables(sra_con)

sra_con <- dbConnect(SQLite(),sqlfile)

rs <- getSRAinfo(c("SRP058527"), sra_con, sraType = "sra" )
sra_ids <- listSRAfile("SRP058527", sra_con)

check_sra <- function(sra_id)
{
    sra_file <- paste0("~/ncbi/sradb/SRP058527/", sra_id, ".sra")
    return(file.exists(sra_file))
}

if(!all(sapply(sra_ids$run, check_sra))) {
    getSRAfile('SRP058527', sra_con, fileType='sra', destDir = "~/ncbi/sradb/SRP058527")
}

sra_paths <- list.files("~/ncbi/sradb/SRP058527", full.names = TRUE)

dump_fastq <- function(path)
{
    system(paste("~/ncbi/sratoolkit.2.8.2-1-mac64/bin/fastq-dump",
                 "-O ~/ncbi/sradb/SRP058527_FASTQ/",
                 path))
}

sapply(sra_paths, dump_fastq)

### Metadata downloaded manually from https://www.ncbi.nlm.nih.gov/Traces/study/

### Lyrata transcriptome and proteome downloaded from
### ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Arabidopsis_lyrata/all_assembly_versions/GCF_000004255.1_v.1.0/GCF_000004255.1_v.1.0_rna.fna.gz
### ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Arabidopsis_lyrata/all_assembly_versions/GCF_000004255.1_v.1.0/GCF_000004255.1_v.1.0_protein.faa.gz

### Arabidopsis transcriptome downloaded from
### ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Arabidopsis_thaliana/latest_assembly_versions/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_protein.faa.gz
### ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Arabidopsis_thaliana/latest_assembly_versions/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_rna.fna.gz

### index generation for lyrata -----------------------------------------------------------------

# #!bin/bash
#
# ## make index
# ./Salmon-0.7.2_linux_x86_64/bin/salmon\
# index -t GCF_000004255.1_v.1.0_rna.fna\
# -i al_transcripts_index\
# --type quasi -k 31

### index generation for thaliana -----------------------------------------------------------------

# #!bin/bash
#
# ## make index
# ./Salmon-0.7.2_linux_x86_64/bin/salmon\
# index -t GCF_000001735.3_TAIR10_rna.fna\
# -i at_transcripts_index\
# --type quasi -k 31

### Mapping for lyrata -----------------------------------------------------------------

# #!bin/bash
#
#
# for f in lyrata/*.fastq; do
#
# base=${f%.fastq}
# echo $base
# base=${base#lyrata/}
#     echo $base
#
#     ./Salmon-0.7.2_linux_x86_64/bin/salmon quant \
#     -i al_transcripts_index -l U \
#     -r $f -o $base\
#     -p 8
#     done

### Mapping for thaliana -----------------------------------------------------------------

# #!bin/bash
#
#
# for f in thaliana/*.fastq; do
#
# base=${f%.fastq}
# echo $base
# base=${base#thaliana/}
#     echo $base
#
#     ./Salmon-0.7.2_linux_x86_64/bin/salmon quant \
#     -i at_transcripts_index -l U \
#     -r $f -o $base\
#     -p 8
#     done

### Make orthogroups -----------------------------------------------------------------

# #!/bin/bash
#
# python OrthoFinder-master/orthofinder.py -f lyrata/orthogroups -t 3
#
# python OrthoFinder-master/orthofinder.py -b ~/lyrata/orthogroups/Results_Apr03/WorkingDirectory
