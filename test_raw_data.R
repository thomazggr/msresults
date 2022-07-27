suppressMessages(
    suppressWarnings({
        library(GEOquery)
        library(ROCR)
        library(dplyr)
        library(biomaRt)
        library(stringr)
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
                         host = "https://www.ensembl.org")#"https://dec2021.archive.ensembl.org")
    attributes_BM <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version",
                                                "ensembl_gene_id", "external_gene_name", "description",
                                                "transcript_biotype", "entrezgene_id", "illumina_humanwg_6_v3",
                                                "affy_hg_u133_plus_2"),
                                    mart = biomart_dataset)
    attributes_BM <- dplyr::rename(attributes_BM, 
                                    target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id,
                                    ext_gene = external_gene_name, entrez_id = entrezgene_id,
                                    illumina = affy_hg_u133_plus_2)
    attributes_BM <- dplyr::select(attributes_BM, 
                                    c('ext_gene', 'entrez_id', 'illumina')) %>%
                        dplyr::mutate(illumina = stringr::str_replace(illumina, "_s_at", ""),
                                      entrez_id = as.character(entrez_id))
    
    return(unique(attributes_BM))
}

Sys.setenv(VROOM_CONNECTION_SIZE = "500000")

get_expression <- function(gse, annot_gpl, gpl, samples, annotation_data, gene_names_filter, labels){
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

    expression_data <- as.data.frame(exprs(gset))
    expression_data$illumina <- row.names(expression_data)
    expression_data <- expression_data %>%
        dplyr::mutate(illumina = stringr::str_replace(illumina, "_at", ""))
    expression_data$entrez_id <- expression_data$illumina
    
    exp_data_att <- expression_data %>%
                    dplyr::left_join(annotation_data, 
                                by=c("entrez_id"))

    exp_data_filtered <- dplyr::filter(exp_data_att, ext_gene %in% gene_names_filter) %>%
                            dplyr::distinct(ext_gene, .keep_all=TRUE)

    row.names(exp_data_filtered) <- exp_data_filtered$ext_gene
    print(head(exp_data_filtered))
    exp_data_filtered <- dplyr::select(exp_data_filtered, -c(illumina.x, illumina.y, entrez_id, ext_gene)) %>%
                            t() %>%
                            as.data.frame()
    
    exp_data_filtered$labels <- labels

    return(exp_data_filtered)
}

#d <- read.csv("~/Projects/msresults/automated_results/GSE37031/GSE37031.tsv", sep="\t")

#d %>% dplyr::filter(SPOT_ID %in% gene_entrez)

gene_entrez <- c("51280",
                 "3485",
                 "10537",
                 "10912",
                 "7057",
                 "2353",
                 "3484",
                 "80303",
                 "8835",
                 "8660",
                 "83716",
                 "55908",
                 "29947", "labels")

gene_names_filter <- c("FOS", "IGFBP1", "IGFBP2", "THBS1", "IRS2", "SOCS2", 
"UBD", "GADD45G", "DNMT3L", "GOLM1", "ANGPTL8", "EFHD1", "CRISPLD2", "labels")

#gse33814 <- read.csv("~/Projects/msresults/automated_results/GSE33814/GSE33814.tsv", sep="\t")

#gse33814 %>% dplyr::filter(Gene.ID %in% gene_entrez) %>% dplyr::select(c("Gene.symbol","Gene.ID", "Gene.title", "logFC"))

super_df <- data.frame(matrix(ncol = 14, nrow = 0))
colnames(super_df) <- gene_names_filter

datasets <- list(
    # list("GSE33814", TRUE, "GPL6884", "1XXX1X111111111XX1XXXXXX000XX0XXX00000X00X00", 
    #     c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
    #         TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
    #         FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
    # list("GSE89632", FALSE, "GPL14951", "X1XXXXXX11X11XX111X1X111XXXX1X11111X10000000000000000X000X00000", 
    #     c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
    #         TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
    #         FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
    list("GSE37031", FALSE, "GPL14877", "111111110000000",
        c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
)

ann <- get.annotation(organism="hsa")

for(d in datasets){
    # print(d)
    # quit()
    data <- get_expression(d[[1]], d[[2]], d[[3]], d[[4]], ann, gene_names_filter, d[[5]])
    super_df <- dplyr::bind_rows(super_df, data)
}

print(head(super_df, 3))
#quit()
for(col in gene_names_filter){
    if(col %in% names(super_df)){
        curr_df <- na.omit(dplyr::select(super_df, c(col, "labels")))
        if(nrow(curr_df) == 0){
            print("")
            print(paste("No rows for", col, "after `na.omit`"))
            print("")
        } else{
            pred <- ROCR::prediction(curr_df[[col]], curr_df[["labels"]])

            # perf <- ROCR::performance(pred, "tpr", "fpr")

            auc <- ROCR::performance(pred, "auc")
            
            print("")
            print(paste(auc@y.name, col))
            print(auc@y.values)
            print("")
            # ROCR::plot(perf)
        }
        
        
    } else{
        print("")
        print(paste("No column for", col))
        print("")
    }
    
}
    
#plot(perf)
