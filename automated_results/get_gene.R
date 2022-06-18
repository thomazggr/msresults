suppressMessages(
    suppressWarnings({
        library("biomaRt")
        library("dplyr")
        library("docopt")
    })
)

create.args <- function() {
    "Prints gene names based on entrez gene ids 

    Usage:
      get_gene.R [-g --geneids] [-o --organism]

    Options:
        -h --help     Show this screen.

        -g <gids> --geneids <gids> List of entrez gene ids separated by comma

        -o <org> --organism <org> Organism name to be used in biomaRt, either mmu or hsa.
    " -> doc

    opt <- docopt::docopt(doc)

    return(opt)
}

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

args <- create.args()
list_genes <- strsplit(args$geneids, split=",")
list_genes <- as.data.frame(list_genes)
names(list_genes) <- c("entrez_id")
list_genes <- transform(list_genes, entrez_id = as.integer(entrez_id))

ann <- get.annotation(organism="hsa")
subset <- list_genes %>%
            dplyr::left_join(ann, 
                            by=c("entrez_id")) %>%
            dplyr::select(c("entrez_id", "ext_gene"))

print(unique(subset))