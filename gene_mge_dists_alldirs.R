###Nearest Neighbours MGE analysis###

#Julianna Paulsen #Jan_15_2025

#Load libraries
library(ape)
library(dplyr)
library(GenomicRanges)
library(stringr)
library(tidyr)
library(tools)
library(purrr)
################################################
#################################################

setwd("C:/Users/jpaul/OneDrive - Eastern Washington University/Desktop/ThesisData/R_Stats_11_2024/gff3s/")
wd <- "C:/Users/jpaul/OneDrive - Eastern Washington University/Desktop/ThesisData/R_Stats_11_2024/gff3s/"

#read genome stats data
genome_stats <- read.csv("../genome_stats.csv")

#create a list of the files to be processed
files <- list.files(wd)

#create a list of the species to be processed by removing the file endings
species <- (as.list(file_path_sans_ext(files)))
species <- keep(species, ~ !str_ends(., "_rm") & !str_ends(., ".scaffolds.fa"))
#species <- c("Acarospora_socialis")
#create a named list 
mge_types <- list(
  retro = c("LTR", "LINE", "SINE", "PLE"),
  dna = c("DNA", "RC")
)

#initialize empty list to store results later
all_results <- list()

nearest_within_contig <- function(repeats_gr, genes_gr) {
  # Split repeats and genes by contig
  repeats_split <- split(repeats_gr, seqnames(repeats_gr))
  genes_split <- split(genes_gr, seqnames(genes_gr))
  
  # Initialize a vector to hold indices of nearest genes in genes_gr
  nearest_indices <- integer(length(repeats_gr))
  nearest_indices[] <- NA
  
  # Loop through each contig
  for (ctg in names(repeats_split)) {
    if (ctg %in% names(genes_split)) {
      reps_ctg <- repeats_split[[ctg]]
      genes_ctg <- genes_split[[ctg]]
      
      # Find nearest genes regardless of direction
      idx <- nearest(reps_ctg, genes_ctg)
      
      ctg_repeats_idx <- which(seqnames(repeats_gr) == ctg)
      
      # Only update those repeats where idx is not NA
      non_na_idx <- which(!is.na(idx))
      if (length(non_na_idx) > 0) {
        matched_genes <- match(genes_ctg[idx[non_na_idx]], genes_gr)
        
        # Only use matches that are not NA
        valid_matches <- which(!is.na(matched_genes))
        if (length(valid_matches) > 0) {
          nearest_indices[ctg_repeats_idx[non_na_idx[valid_matches]]] <- matched_genes[valid_matches]
        }
      }
    }
  }
  
  return(nearest_indices)
}


for (x in species) {
  tryCatch({
    # Load results files and format data
    genes_file_name <- paste0(x, ".gff3")
    reps_file_name <- paste0(x, "_rm.gff3")
    genes <- read.gff(genes_file_name)
    reps <- read.gff(reps_file_name)
    
    # Extract ID # from reps into new column
    reps$ID <- sub(".*ID=([0-9]+);.*", "\\1", reps$attributes)
    
    # Load .out file
    out_df <- read.delim(paste0(x, ".scaffolds.fa.out"), header = FALSE, skip = 3, sep = "")
    colnames(out_df) <- c(
      "SW_score", "perc_div", "perc_del", "perc_ins",
      "sequence", "query_begin", "query_end", "query_left", "strand",
      "repeat", "class_family",
      "repeat_begin", "repeat_end", "repeat_left",
      "ID"
    )
    
    # Filter erroneous rows
    out_df_filtered <- out_df[out_df$SW_score != "*", ]
    
    # Merge .out and reps
    merged <- merge(reps, out_df_filtered[, c("ID", "class_family")], by = "ID", all.x = TRUE)
    
    # Replace 'dispersed_repeat' with the actual class name
    merged$type <- ifelse(!is.na(merged$class_family), merged$class_family, merged$type)
    
    reps <- merged[, c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")]
    
    # Loop through each MGE class (retro, dna)
    for (z in names(mge_types)) {
      
      values_to_filter <- mge_types[[z]]
      
      reps_subset <- reps %>%
        filter(str_detect(type, paste(values_to_filter, collapse = "|")))
      
      # Filter genes and repeats based on strand and type
      genes_filtered <- genes[genes$type == "gene", ]
      reps_filtered <- reps_subset[!reps_subset$type %in% c("Simple_repeat", "Low_complexity"), ]
      
      # Remove mitochondrial contigs
      genes2 <- genes_filtered %>%
        filter(str_detect(seqid, "^(contig|Scaffold|SCAF|Contig|scaffold|bac_rub)"))
      
      reps2 <- reps_filtered %>%
        filter(str_detect(seqid, "^(contig|Scaffold|SCAF|Contig|scaffold|bac_rub)"))
      
      # Create GenomicRanges objects
      genes_gr <- GRanges(
        seqnames = genes2$seqid,
        ranges = IRanges(genes2$start, genes2$end),
        strand = genes2$strand,
        type = genes2$type
      )
      
      reps_gr <- GRanges(
        seqnames = reps2$seqid,
        ranges = IRanges(reps2$start, reps2$end),
        strand = reps2$strand,
        type = reps2$type
      )
      
      # Resize both genes and repeats to include their midpoints
      gene_mids <- resize(genes_gr, width = 1, fix = "center")
      repeat_mids <- resize(reps_gr, width = 1, fix = "center")
      
      
      # Find nearest genes
      nearest_genes <- nearest_within_contig(repeat_mids, gene_mids)
      
      
      valid_indices <- !is.na(nearest_genes)
      filtered_repeats <- repeat_mids[valid_indices]
      filtered_genes <- gene_mids[nearest_genes[valid_indices]]
      
      
      if (length(filtered_repeats) > 0 && length(filtered_genes) > 0) {
        distances <- distance(filtered_repeats, filtered_genes)
        
        result_df <- data.frame(
          repetitive_element_chr = as.character(seqnames(filtered_repeats)),
          repetitive_element_start = start(filtered_repeats),
          repetitive_element_end = end(filtered_repeats),
          nearest_gene_chr = as.character(seqnames(filtered_genes)),
          nearest_gene_start = start(filtered_genes),
          distance = distances
        )
        
        mean_result_per_sp <- result_df %>%
          summarize(
            mean_distance_to_gene = mean(distance, na.rm = TRUE)
          ) %>%
          mutate(tree.name = x, .before = 1) %>%
          left_join(genome_stats %>% select(tree.name, Total_Size, Num_Genes), by = "tree.name") %>%
          mutate(
            gene_density = Num_Genes / Total_Size,
            adjusted_distance_value = mean_distance_to_gene * gene_density
          )
        
        mean_result <- mean_result_per_sp %>%
          select(tree.name, adjusted_distance_value)
        
      } else {
        # Assign NA when no valid repeat-gene pairs found
        mean_result <- data.frame(
          tree.name = x,
          adjusted_distance_value = NA
        )
      }
      
      # Store result
      if (!x %in% names(all_results)) all_results[[x]] <- list()
      all_results[[x]][[z]] <- mean_result
      
      # (Optional) Also assign to env
      assign(paste0(x, "_", z, "_meanresult"), mean_result)
      
    }
    
  }, error = function(e) {
    message(paste("Error processing species", x, ":", e$message))
  })
}

# Flatten the nested list into a single data frame
results_df <- bind_rows(
  lapply(names(all_results), function(species) {
    lapply(names(all_results[[species]]), function(mge_type) {
      df <- all_results[[species]][[mge_type]]
      df$mge_type <- mge_type
      return(df)
    }) %>% bind_rows()
  })
)

results_wide <- results_df %>%
  pivot_wider(names_from = mge_type, values_from = adjusted_distance_value)

# Export to CSV
write.csv(results_wide, "../Spatial_dists_MGE_analyses/MGE_nearest_neighbor_alldirections_gene.csv", row.names = FALSE)
