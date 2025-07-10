#Gene density heatmap

library(tidyr)
library(biocLite)
library(GenomicRanges)
library(ape)
library(dplyr)
library(purrr)
library(stringr)
library(tools)
setwd("C:/Users/jpaul/OneDrive - Eastern Washington University/Desktop/ThesisData/R_Stats_11_2024/gff3s/")
wd <- "C:/Users/jpaul/OneDrive - Eastern Washington University/Desktop/ThesisData/R_Stats_11_2024/gff3s/"

#load trait data
traits <- read.csv("../Traits.csv")
#keep only rarity
rarity <- traits %>%
  select(tree.name, reproductive.mode)

#create a list of the files to be processed
files <- list.files(wd)

#create a list of the species to be processed by removing the file endings
species <- (as.list(file_path_sans_ext(files)))
species <- keep(species, ~ !str_ends(., "_rm") & !str_ends(., ".scaffolds.fa"))
#species <- c("Usnea_strigosa", "Usnea_subfusca")

all_windows_list <- list()
for (x in species) {
    # Load results files and format data
    genes_file_name <- paste0(x, ".gff3")
    reps_file_name <- paste0(x, "_rm.gff3")
    genes <- read.gff(genes_file_name)
    reps <- read.gff(reps_file_name)
    #extract ID # from reps into new column
    reps$ID <- sub(".*ID=([0-9]+);.*", "\\1", reps$attributes)
    
    #load out.df
    out_df <- read.delim(paste0(x, ".scaffolds.fa.out"), header = FALSE, skip = 3, sep = "")
    colnames(out_df) <- c(
      "SW_score", "perc_div", "perc_del", "perc_ins",
      "sequence", "query_begin", "query_end", "query_left", "strand",
      "repeat", "class_family",
      "repeat_begin", "repeat_end", "repeat_left",
      "ID"
    )
    
    #filter erroneous rows
    out_df_filtered <- out_df[out_df$SW_score != "*", ]
    
    #merge .out and reps
    merged <- merge(reps, out_df_filtered[, c("ID", "class_family")], by = "ID", all.x = TRUE)
    
    #Replace 'dispersed_repeat' with the actual class name
    merged$type <- ifelse(!is.na(merged$class_family), merged$class_family, merged$type)
    
    reps <- merged[, c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")]
    
    #filter and remove simple repeats
    genes <- genes[genes$type == "gene", ]
    reps <- reps[reps$type != "Simple_repeat" & reps$type != "Low_complexity", ]
    
    genes_gr <- GRanges(
      seqnames = genes$seqid,
      ranges = IRanges(genes$start, genes$end),
      strand = genes$strand,
      type = genes$type
    )
    
    reps_gr <- GRanges(
      seqnames = reps$seqid,
      ranges = IRanges(reps$start, reps$end),
      strand = reps$strand,
      type = reps$type
    )
    
    #load chrom sizes
    chromsizes <- read.table(file = paste0("../Chromsizes/", x, ".scaffolds_chromsizes.tsv"), sep = '\t', header = FALSE)
    
    seqinfo <- Seqinfo(
      seqnames = chromsizes$V1,
      seqlengths = chromsizes$V2
    )
    
    #create windows
    
    windows <- tileGenome(seqinfo, tilewidth=100000, cut.last.tile.in.chrom=T)
    
    windows$gene_density <- countOverlaps(windows, genes_gr)
    
    windows$mge_density <- countOverlaps(windows, reps_gr)
    
    assign(paste0(x, "_windows"), windows)
    
    windows_df <- as.data.frame(windows)
    
    windows_df$tree.name <-x
    
    merged_windows_df <- merge(windows_df, rarity)
    
    assign(paste0(x, "_windowsdf"), merged_windows_df)
    
    all_windows_list[[x]] <- merged_windows_df
}

# Combine all species
combined_windows_df <- bind_rows(all_windows_list)

custom_labels <- c(
  "short" = "Short",
  "medium" = "Medium",
  "long" = "Long"
)

custom_labels <- c(
  "A" = "Asexual",
  "S" = "Sexual"
)

ggplot(combined_windows_df, aes(x = gene_density, y = mge_density)) +
  geom_bin2d() +
  scale_fill_gradient(low = "lavender", high = "darkslateblue") +
  facet_wrap(~gen_length, labeller = labeller(gen_length = custom_labels)) +
  labs(x = "Gene count per 100kb window",
       y = "MGE count per 100kb window") +
  theme_minimal()

ggplot(combined_windows_df, aes(x = gene_density, y = mge_density)) +
  geom_bin2d() +
  scale_fill_gradient(low = "lavender", high = "darkslateblue") +
  facet_wrap(~reproductive.mode, labeller = labeller(reproductive.mode = custom_labels)) +
  labs(x = "Gene count per 100kb window",
       y = "MGE count per 100kb window") +
  theme_minimal()


ggsave(filename="../presentation_figs/gene_densities_repro.png", plot=last_plot(), dpi=1200, width=6, height=5, unit="in", bg = "transparent")
