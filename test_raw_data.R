suppressMessages(
    suppressWarnings({
        library(GEOquery)
        library(ROCR)
        library(dplyr)
        library(biomaRt)
    })
)

get.annotation <- function(organism = NULL) {
    if (organism %in% c("mmu", "mmusculus_gene_ensembl")) {
        org_ensembl <- "mmusculus_gene_ensembl"
    } else if (organism %in% c("hsa", "hsapiens_gene_ensembl")) {
        org_ensembl <- "hsapiens_gene_ensembl"
    } else {
        stop("Organism needs to be passed as either mmu or hsa.")
    }

    biomart_dataset <- biomaRt::useMart(biomart = "ensembl",
                         dataset = org_ensembl,
                         host = "https://dec2021.archive.ensembl.org")#"https://www.ensembl.org")
    attributes_BM <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version",
                                                "ensembl_gene_id", "external_gene_name", "description",
                                                "transcript_biotype", "entrezgene_id", "illumina_humanwg_6_v3"),
                                    mart = biomart_dataset)
    attributes_BM <- dplyr::rename(attributes_BM, 
                                    target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id,
                                    ext_gene = external_gene_name, entrez_id = entrezgene_id,
                                    illumina = illumina_humanwg_6_v3)
    attributes_BM <- dplyr::select(attributes_BM, 
                                    c('ext_gene', 'entrez_id', 'illumina'))
    
    return(unique(attributes_BM))
}

Sys.setenv(VROOM_CONNECTION_SIZE = "500000")

gse <- "GSE33814"
annot_gpl <- TRUE
gpl <- "GPL6884"
samples <- "1XXX1X111111111XX1XXXXXX000XX0XXX00000X00X00"

gset <- GEOquery::getGEO(gse, GSEMatrix =TRUE, AnnotGPL=annot_gpl)
if (length(gset) > 1) idx <- grep(gpl, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- samples
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#head(gset)

ilmn_filter <- c("ILMN_1766405", "ILMN_1725193", "ILMN_1678841", "ILMN_1651498", 
"ILMN_1686116", "ILMN_1669523", "ILMN_1728445", "ILMN_1779448", "ILMN_1798926", 
"ILMN_2083469", "ILMN_1790689", "ILMN_1764714", "ILMN_1778357")

gene_names_filter <- c("FOS", "IGFBP1", "IGFBP2", "THBS1", "IRS2", "SOCS2", 
"UBD", "GADD45G", "DNMT3L", "GOLM1", "ANGPTL8", "EFHD1", "CRISPLD2")

expression_data <- as.data.frame(exprs(gset))
expression_data$illumina <- row.names(expression_data)

annotation_data <- get.annotation(organism="hsa")

exp_data_att <- expression_data %>%
                dplyr::left_join(annotation_data, 
                            by=c("illumina"))

exp_data_filtered <- dplyr::filter(exp_data_att, ext_gene %in% gene_names_filter) %>%
                        dplyr::distinct(ext_gene, .keep_all=TRUE)

row.names(exp_data_filtered) <- exp_data_filtered$ext_gene
exp_data_filtered <- dplyr::select(exp_data_filtered, -c(illumina, entrez_id, ext_gene)) %>%
                        t() %>%
                        as.data.frame()

exp_data_filtered$labels <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
print(head(exp_data_filtered, 3))
# quit()
for(col in gene_names_filter){
    if(col %in% names(exp_data_filtered)){
        pred <- ROCR::prediction(exp_data_filtered[[col]], exp_data_filtered[["labels"]])

        #perf <- ROCR::performance(pred, "tpr", "fpr")

        auc <- ROCR::performance(pred, "auc")
        print(paste(auc@y.name, col))
        print(auc@y.values)
    } else{
        print(paste("No column for", col))
    }
    
}
    
#plot(perf)
