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
genes <- c("GOLM1", "EFHD1", "IGFBP2", "IRS2", "UBD", "GADD45G", "DNMT3L", "IGFBP1", "FOS", "THBS1", "SOCS2", "CRISPLD2", "AKR1B10", "C15orf52", "ANXA13", "ZNF385B", "FOXC1", "CDH15", "PIGR", "FCGBP", "CCL20", "LGALS4", "PRSS3", "KDELR3", "TNFRSF12A", "CCL2", "CA12", "ZMAT3", "EEF1A2", "TAGLN", "SLITRK3", "NR5A2", "LARP1B", "CDH23", "CDKN1A", "PMP22", "STMN2", "TRIM15", "CH25H", "LDLR", "HSPA1A", "ANGPTL8", "GPR88", "IL32", "FAT1", "C11orf96", "N4BP2L1", "OAT", "NR4A2", "ATF3", "PHLDA1", "NFIL3", "EGR1", "LRRC31", "KLF6", "KLF4", "IER3", "CEBPD", "ZFP36", "SDC4", "RCAN1", "SYT7", "MYC", "FMO1", "DUSP10", "TP53I3", "TMEM45B", "ADAMTS1", "NR4A1", "JUNB", "GPAM", "IL1RN", "P4HA1", "MAMDC4", "FOSB", "JUN", "ME1", "APOL3", "MAFF", "SIK1")
genes <- c("120224", "8857", "5284", "1026", "3557", "11075", "90634", "8660", "55132", "79782", "3960", "387763", "7538", "29947", "11221", "9510", "4609", "4199", "80303", "151126", "5646", "3726", "22822", "8870", "3484", "158056", "51330", "3485", "312", "2353", "9023", "2494", "83716", "55908", "1917", "89870", "150094", "11015", "1013", "467", "6876", "2296", "4942", "1052", "6385", "54112", "7057", "8835", "4783", "64393", "10537", "22865", "64072", "6347", "3949", "80833", "771", "9235", "1958", "3725", "51280", "6364", "2326", "9314", "9066", "10912", "4929", "1827", "3164", "1316", "57016", "9540", "5033", "57678", "23764", "5376", "2195", "2354", "3303")
pathway.enrichment(gene_list = genes, organism = "hsa", path = "~/Documents/Projects/msresults/automated_results/padj005_logfc1", kegg_file = "results", go_file = "results")


print("p005_logfc1 genes >> minimum 3/5 datasets")
genes2 <- c("IGFBP2", "GADD45G", "IGFBP1", "FOS", "SOCS2", "NR4A2", "IGF1", "CYP1A1", "ENPP2", "GPR88", "GPAT3", "ACKR3", "GOLM1", "EFHD1", "IRS2", "UBD", "CCL2", "DNMT3L", "CRYAA", "THBS1", "CRISPLD2", "TMEM154", "ENO3", "FADS1", "FADS2", "TRHDE", "MYC", "NR4A1", "P4HA1", "FOSB", "ME1", "SIK1", "RTP3", "PHLDA1", "IP6K3", "CEBPD", "SDC4", "TMEM45B", "ADAMTS1", "IL1RN", "STMN2", "CCL20", "PEG10")
genes2 <- c("120224", "3557", "83597", "1052", "6385", "54112", "11075", "8660", "29953", "7057", "2027", "8835", "9415", "29947", "9510", "4609", "10537", "23089", "4199", "5168", "6347", "80303", "84803", "22822", "3992", "117283", "3484", "51280", "6364", "3485", "3479", "2353", "201799", "10912", "4929", "3164", "5033", "83716", "1543", "57007", "2354", "150094", "1409")
pathway.enrichment(gene_list = genes2, organism = "hsa", path = "~/Documents/Projects/msresults/automated_results/p005_logfc1", kegg_file = "results", go_file = "results")


print("p001_logfc1 genes >> minimum 3/5 datasets")
genes3 <- c("IGFBP2", "GADD45G", "SOCS2", "GPR88", "GPAT3", "GOLM1", "EFHD1", "IRS2", "UBD", "DNMT3L", "IGFBP1", "FOS", "THBS1", "CRISPLD2", "TMEM154", "FADS2", "TRHDE", "P4HA1", "FOSB", "ME1", "RTP3", "PHLDA1", "IP6K3", "CEBPD", "SDC4", "IL1RN", "STMN2", "CCL20", "PEG10")
genes3 <- c("3557", "83597", "1052", "6385", "54112", "11075", "8660", "29953", "7057", "8835", "9415", "29947", "10537", "23089", "4199", "80303", "84803", "22822", "117283", "3484", "51280", "6364", "3485", "2353", "201799", "10912", "5033", "83716", "2354")
pathway.enrichment(gene_list = genes3, organism = "hsa", path = "~/Documents/Projects/msresults/automated_results/p001_logfc1", kegg_file = "results", go_file = "results")
