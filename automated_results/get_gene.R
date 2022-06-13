library("biomaRt")
library("dplyr")

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
                                                "transcript_biotype", "entrezgene_id"),
                                    mart = biomart_dataset)
    attributes_BM <- dplyr::rename(attributes_BM, 
                                    target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id,
                                    ext_gene = external_gene_name, entrez_id = entrezgene_id)
    attributes_BM <- dplyr::select(attributes_BM, 
                                    c('target_id', 'ens_gene', 'ext_gene', 'entrez_id'))
    
    return(attributes_BM)
}

list_genes <- c("51280", "3485", "10537", "10912", "7057", "2353", "3484", "80303", "8835", "8660", "83716", "55908", "29947", "54463", "54112", "3949", "4942", "2195", "387763", "90634", "3303", "9235", "56892", "11015", "388115", "9023", "55132", "1026", "5646", "159371", "8857", "6347", "3960", "64393", "54825", "5376", "6876", "771", "11075", "57016", "5284", "151126", "312", "2296", "89870", "2494", "1013", "51330", "64072", "3576", "84803", "22865", "1917", "6364", "150094", "6385", "120224", "9066", "2326", "80833", "8870", "90141", "9314", "3726", "64651", "2354", "4609", "467", "1827", "1316", "1052", "4783", "22822", "28984", "3725", "5033", "57678", "8343", "11221", "84962", "23764", "79782", "7538", "4199", "9540", "1958", "3491", "4929", "158056", "9510", "3557", "3164")

ann <- get.annotation(organism="hsa")
ss <- subset(ann, entrez_id %in% list_genes)
print(unique(ss$ext_gene))