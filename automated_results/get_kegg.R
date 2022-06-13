library("clusterProfiler")
library("ggplot2")
library("dplyr")


pathway.enrichment <- function(gene_list = NULL, organism = NULL, path = NULL, kegg_file = NULL, go_file = NULL) {
    # - - - - - - - - - - - - -
    # Enrich with KEGG database
    # - - - - - - - - - - - - -
    if (is.null(gene_list)) { stop("No gene list has been passed.") }
    if (is.null(path)) { stop("No path has been passed.") }
    if (is.null(organism)) { stop("No organism passed. Needs to be either mmu or hsa.") }

    kegg_enrich <- clusterProfiler::enrichKEGG(gene = gene_list,
                                organism = organism,
                                pvalueCutoff = 0.05)

    print("Generate KEGG pathway results")
    # Get table of results
    kegg_table <- as.data.frame(kegg_enrich) %>% 
                    dplyr::arrange(desc(-log10(pvalue)))
    kegg_table <- head(kegg_table, n = 20)
            

    # KEGG plot
    kegg_bar <- ggplot2::ggplot(kegg_table, 
                                ggplot2::aes(x = reorder(Description, 
                                                         -log10(pvalue)), 
                                             y = -log10(pvalue))) +
      ggplot2::geom_bar(stat = 'identity', fill = "#6B2525") +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::labs(title = "KEGG Enrichment Pathways", 
                    x = "KEGG Terms") +
      ggplot2::coord_flip()

    kegg_full <- as.data.frame(kegg_enrich) %>% 
      dplyr::arrange(desc(-log10(pvalue)))
    
    write.table(kegg_full, 
                file.path(path, paste0(kegg_file, "__kegg_full_table.txt")), 
                sep = "\t", 
                quote = F, 
                row.names = F)
    
    ggplot2::ggsave(file.path(path, paste0(kegg_file, "__kegg_bar.png")), 
           plot = kegg_bar, 
           device = "png")
    
    print("Saved KEGG full table and bar plot")

    # - - - - - - - - - - - - -
    # Enrich with GO
    # - - - - - - - - - - - - -
    print("Generate GO barplot")
    if (is.null(organism)) {
        stop("No organism passed. Needs to be either mmu or hsa.")
    } else if (organism == "mmu") {
        go_enrich <- clusterProfiler::enrichGO(gene = gene_list,
                              OrgDb = "org.Mm.eg.db",
                              ont = "BP",
                              pvalueCutoff = 0.05)
    } else if (organism == "hsa") {
        go_enrich <- clusterProfiler::enrichGO(gene = gene_list,
                              OrgDb = "org.Hs.eg.db",
                              ont = "BP",
                              pvalueCutoff = 0.05)
    }

    # Get table of results
    go_table <- as.data.frame(go_enrich) %>%
                dplyr::arrange(desc(-log10(pvalue)))
    go_table <- head(go_table, n = 20)

    # Plot results
    go_bar <- ggplot2::ggplot(go_table, 
                              ggplot2::aes(x = reorder(Description, 
                                                       -log10(pvalue)), 
                                           y = -log10(pvalue))) +
      ggplot2::geom_bar(stat = 'identity', fill = "#157296") +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::labs(title = "GO Biological Process", x = "Go Terms") +
      ggplot2::coord_flip()

    go_full <- as.data.frame(go_enrich) %>% 
      dplyr::arrange(desc(-log10(pvalue)))
    
    write.table(go_full, 
                file.path(path, paste0(go_file, "__go_full_table.txt")), 
                sep = "\t", 
                quote = F, 
                row.names = F)
    
    ggplot2::ggsave(file.path(path, paste0(go_file, "__go_bar.png")), 
           plot = go_bar, 
           device = "png")
    
    print("Saved GO full table and bar plot")
}

print("padj005_logfc1 genes >> minimum 2/3 datasets")
#genes <- c("GOLM1", "EFHD1", "IGFBP2", "IRS2", "UBD", "GADD45G", "DNMT3L", "IGFBP1", "FOS", "THBS1", "SOCS2", "CRISPLD2", "AKR1B10", "C15orf52", "ANXA13", "ZNF385B", "FOXC1", "CDH15", "PIGR", "FCGBP", "CCL20", "LGALS4", "PRSS3", "KDELR3", "TNFRSF12A", "CCL2", "CA12", "ZMAT3", "EEF1A2", "TAGLN", "SLITRK3", "NR5A2", "LARP1B", "CDH23", "CDKN1A", "PMP22", "STMN2", "TRIM15", "CH25H", "LDLR", "HSPA1A", "ANGPTL8", "GPR88", "IL32", "FAT1", "C11orf96", "N4BP2L1", "OAT", "NR4A2", "ATF3", "PHLDA1", "NFIL3", "EGR1", "LRRC31", "KLF6", "KLF4", "IER3", "CEBPD", "ZFP36", "SDC4", "RCAN1", "SYT7", "MYC", "FMO1", "DUSP10", "TP53I3", "TMEM45B", "ADAMTS1", "NR4A1", "JUNB", "GPAM", "IL1RN", "P4HA1", "MAMDC4", "FOSB", "JUN", "ME1", "APOL3", "MAFF", "SIK1")
genes <- c("51280", "3485", "10537", "10912", "7057", "2353", "3484", "80303", "8835", "8660", "83716", "55908", "29947", "54463", "54112", "3949", "4942", "2195", "387763", "90634", "3303", "9235", "56892", "11015", "388115", "9023", "55132", "1026", "5646", "159371", "8857", "6347", "3960", "64393", "54825", "5376", "6876", "771", "11075", "57016", "5284", "151126", "312", "2296", "89870", "2494", "1013", "51330", "64072", "3576", "84803", "22865", "1917", "6364", "150094", "6385", "120224", "9066", "2326", "80833", "8870", "90141", "9314", "3726", "64651", "2354", "4609", "467", "1827", "1316", "1052", "4783", "22822", "28984", "3725", "5033", "57678", "8343", "11221", "84962", "23764", "79782", "7538", "4199", "9540", "1958", "3491", "4929", "158056", "9510", "3557", "3164")
pathway.enrichment(gene_list = genes, organism = "hsa", path = "~/Documents/Projects/msresults/automated_results/padj005_logfc1", kegg_file = "results", go_file = "results")


print("p005_logfc1 genes >> minimum 3/5 datasets")
#genes2 <- c("IGFBP2", "GADD45G", "IGFBP1", "FOS", "SOCS2", "NR4A2", "IGF1", "CYP1A1", "ENPP2", "GPR88", "GPAT3", "ACKR3", "GOLM1", "EFHD1", "IRS2", "UBD", "CCL2", "DNMT3L", "CRYAA", "THBS1", "CRISPLD2", "TMEM154", "ENO3", "FADS1", "FADS2", "TRHDE", "MYC", "NR4A1", "P4HA1", "FOSB", "ME1", "SIK1", "RTP3", "PHLDA1", "IP6K3", "CEBPD", "SDC4", "TMEM45B", "ADAMTS1", "IL1RN", "STMN2", "CCL20", "PEG10")
genes2 <- c("3485", "10912", "2353", "3484", "8835", "84803", "4929", "51280", "10537", "6347", "7057", "80303", "8660", "83716", "55908", "29947", "1409", "1543", "3479", "57007", "54112", "5168", "11075", "6364", "150094", "64651", "2354", "4609", "3992", "5033", "4199", "9415", "29953", "3164", "90523", "6385", "120224", "1052", "22822", "28984", "11067", "83597", "9510", "117283", "3557", "201799", "2027", "23089")
pathway.enrichment(gene_list = genes2, organism = "hsa", path = "~/Documents/Projects/msresults/automated_results/p005_logfc1", kegg_file = "results", go_file = "results")


print("p001_logfc1 genes >> minimum 3/5 datasets")
#genes3 <- c("IGFBP2", "GADD45G", "SOCS2", "GPR88", "GPAT3", "GOLM1", "EFHD1", "IRS2", "UBD", "DNMT3L", "IGFBP1", "FOS", "THBS1", "CRISPLD2", "TMEM154", "FADS2", "TRHDE", "P4HA1", "FOSB", "ME1", "RTP3", "PHLDA1", "IP6K3", "CEBPD", "SDC4", "IL1RN", "STMN2", "CCL20", "PEG10")
genes3 <- c("3485", "10912", "8835", "84803", "51280", "10537", "7057", "2353", "3484", "80303", "8660", "83716", "55908", "29947", "54112", "11075", "6364", "64651", "2354", "5033", "4199", "9415", "29953", "6385", "1052", "22822", "83597", "117283", "3557", "201799", "23089")
pathway.enrichment(gene_list = genes3, organism = "hsa", path = "~/Documents/Projects/msresults/automated_results/p001_logfc1", kegg_file = "results", go_file = "results")
