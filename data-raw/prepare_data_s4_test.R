### Prepare Rdata for testing S4 ogset in "test_s4_structure.R"
### Mostly code lines from vignette

# devtools::load_all()

# Load edata -----------------------------------------------------------------

thaliana <- make_TPM_df("inst/exdata/thaliana")
lyrata <- make_TPM_df("inst/exdata/lyrata")


# adjust colnames ---------------------------------------------------------

sample_info <-  read.table("inst/exdata/SRP058527_metadata.txt",
                           header = TRUE,
                           sep = "\t")
reduce_info <- function(info) {
    if(info == "heat stress 0 h") "0h"
    else if(info == "heat stress 6 h 37\302\260C") "heat_6h"
    else "recovery"
}
sample_info$etime <- vapply(sample_info$exposure_time_s,
                            reduce_info,
                            character(1))
sample_info$etime <- paste(sample_info$etime,
                           c(1, 2),
                           sep = "_")

switch_name <- function(name) {
    if(name %in% sample_info$Run_s) {
        sample_info[sample_info$Run_s == name, "etime", drop = TRUE]
    } else name
}
colnames(thaliana) <- vapply(colnames(thaliana), switch_name, character(1))
colnames(thaliana) <- paste("thal",
                            colnames(thaliana),
                            sep = "_")
str(thaliana)

colnames(lyrata) <- vapply(colnames(lyrata), switch_name, character(1))
colnames(lyrata) <- paste("lyr",
                          colnames(lyrata),
                          sep = "_")
str(lyrata)




# load ogroups ------------------------------------------------------------

ogroups <- parse_orthogroups("inst/exdata/ogroups_sample.csv")


# protein ID to tx ID -----------------------------------------------------

thaliana_id_path <- "inst/exdata/thaliana_features_sample.txt"
thaliana_ids <- read.table(file = thaliana_id_path,
                           sep = "\t", header = FALSE,
                           stringsAsFactors = FALSE,
                           quote = "")

## Not sure why the header is commented with "#"
colnames(thaliana_ids) <-  scan(file = thaliana_id_path,
                                what = "",
                                nlines = 1)[-1]


lyrata_id_path <- "inst/exdata/lyrata_features_sample.txt"
lyrata_ids <- read.table(file = lyrata_id_path,
                         sep = "\t", header = FALSE,
                         stringsAsFactors = FALSE,
                         quote = "",
                         skip = 1,
                         comment.char = "") # One gene name contains "#"

## Not sure why the header is commented with "#"
colnames(lyrata_ids) <-  scan(file = lyrata_id_path,
                              what = "",
                              nlines = 1)[-1]


## thaliana conversion
ogroups <- switch_ids(ogroups = ogroups,
                      ids_table = thaliana_ids,
                      px_id = "product_accession",
                      tx_id = "related_accession",
                      mc.cores = 4)

## lyrata conversion
ogroups <- switch_ids(ogroups = ogroups,
                      ids_table = lyrata_ids,
                      px_id = "product_accession",
                      tx_id = "related_accession",
                      mc.cores = 4)


# collapse orthogroups ----------------------------------------------------

thaliana_og <- collapse_orthologs(eset = thaliana,
                                  ogroups = ogroups,
                                  mc.cores = 4)
lyrata_og <- collapse_orthologs(eset = lyrata,
                                ogroups = ogroups,
                                mc.cores = 4)

og_eset <- merge(thaliana_og, lyrata_og,
                 by = "row.names", all = TRUE)
rownames(og_eset) <- og_eset$Row.names; og_eset$Row.names <- NULL

og_NOMATCH <- og_eset[!complete.cases(og_eset), ]
dim(og_NOMATCH)

og_eset <- og_eset[complete.cases(og_eset), ]
dim(og_eset)


# fit model ---------------------------------------------------------------

og_eset_log <- log(og_eset + 1)
coldata <- data.frame(spc = sapply(colnames(og_eset_log), function(i) strsplit(i, split = "_")[[1]][1]),
                      treat = sapply(colnames(og_eset_log), function(i) strsplit(i, split = "_")[[1]][2]),
                      stringsAsFactors = FALSE)
og_fit <- add_fit(dat = og_eset_log,
                  coldata = coldata,
                  design = ~ spc + treat + spc:treat)

save(list = c("coldata",
              "lyrata",
              "og_eset",
              "og_fit",
              "sample_info",
              "thaliana",
              "ogroups"),
     file = "tests/testthat/tos4.Rdata")
