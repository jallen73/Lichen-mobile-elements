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

#function to split data by contig and calculate nearest gene within contigs only

# Replace nearest_within_contig function with gene-to-MGE version
nearest_downstream_mge <- function(genes_gr, repeats_gr, direction = "follow") {
  genes_split <- split(genes_gr, seqnames(genes_gr))
  repeats_split <- split(repeats_gr, seqnames(repeats_gr))
  
  nearest_indices <- integer(length(genes_gr))
  nearest_indices[] <- NA
  
  for (ctg in names(genes_split)) {
    if (ctg %in% names(repeats_split)) {
      genes_ctg <- genes_split[[ctg]]
      reps_ctg <- repeats_split[[ctg]]
      
      if (direction == "follow") {
        idx <- follow(genes_ctg, reps_ctg)
      } else {
        idx <- precede(genes_ctg, reps_ctg)
      }
      
      ctg_genes_idx <- which(seqnames(genes_gr) == ctg)
      non_na_idx <- which(!is.na(idx))
      if (length(non_na_idx) > 0) {
        nearest_indices[ctg_genes_idx[non_na_idx]] <- match(reps_ctg[idx[non_na_idx]], repeats_gr)
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
    
    ########################################################
    for (z in names(mge_types)) {
      
      values_to_filter <- mge_types[[z]]
      
      reps_subset <- reps %>%
        filter(str_detect(type, paste(values_to_filter, collapse = "|")))
      
      strand_results <- list()
      for (y in c("pos", "neg")) {
        if (y == "pos") {
          genes_filtered <- genes[genes$strand == "+" & genes$type == "mRNA", ]
          reps_filtered <- reps_subset[reps_subset$strand == "+" & reps_subset$type != "Simple_repeat" & reps_subset$type != "Low_complexity", ]
        } else {
          genes_filtered <- genes[genes$strand == "-" & genes$type == "gene", ]
          reps_filtered <- reps_subset[reps_subset$strand == "-" & reps_subset$type != "Simple_repeat" & reps_subset$type != "Low_complexity", ]
        }
        
        # Remove mitochondrial contigs
        genes2 <- genes_filtered %>%
          filter(str_detect(seqid, "^contig") | str_detect(seqid, "^SCAF") | str_detect(seqid, "^bac_rub") | str_detect(seqid, "^Contig") | str_detect(seqid, "^scaffold") | str_detect(seqid, "^Scaffold"))
        
        reps2 <- reps_filtered %>%
          filter(str_detect(seqid, "^contig") | str_detect(seqid, "^SCAF") | str_detect(seqid, "^bac_rub") | str_detect(seqid, "^Contig") | str_detect(seqid, "^scaffold") | str_detect(seqid, "^Scaffold"))
        
        # GenomicRanges
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
        
        # Resize gene GRanges to transcription end (for +) or start (for -)
        gene_points <- resize(genes_gr, width = 1, fix = ifelse(y == "pos", "end", "start"))
        
        # Resize MGE GRanges to transcription start (for +) or end (for -)
        mge_points <- resize(reps_gr, width = 1, fix = ifelse(y == "pos", "start", "end"))
        
        # Find nearest MGE downstream of gene
        nearest_mges <- nearest_downstream_mge(gene_points, mge_points, direction = ifelse(y == "pos", "follow", "precede"))
        
        # Filter out NAs
        valid_indices <- !is.na(nearest_mges)
        filtered_genes <- gene_points[valid_indices]
        
        if (length(filtered_genes) != 0) {
          filtered_mges <- mge_points[nearest_mges[valid_indices]]
          
          distances <- distance(filtered_genes, filtered_mges)
          
          result_df <- data.frame(
            gene_chr = as.character(seqnames(filtered_genes)),
            gene_start = start(filtered_genes),
            mge_chr = as.character(seqnames(filtered_mges)),
            mge_start = start(filtered_mges),
            distance = distances
          )
          
          strand_results[[y]] <- result_df
          
          assign(paste0(x, "_", z, "_", y, "_distances"), result_df, envir = .GlobalEnv)
          
        } else {
          all_results[[x]][[z]] <- data.frame(
            tree.name = x,
            adjusted_distance_value = NA
          )
        }
      }
      
      
      #combine stranded data into a single df
      combined_result_df <- bind_rows(strand_results)
      
      assign(paste0(x, "_", z, "_combined_distances"), combined_result_df)
      
      # Use dplyr to find the mean distance to gene across the whole genome
      mean_result_per_sp <- combined_result_df %>%
        summarize(
          mean_distance_to_gene = mean(distance, na.rm = TRUE)
        ) %>%
        mutate(tree.name = x, .before = 1) %>%
        left_join(genome_stats %>% select(tree.name, Total_Size, Num_Genes), by = "tree.name") %>% # add genome stats
        mutate(
          gene_density = Num_Genes / Total_Size  # calculate "gene density" value
        ) %>%
        mutate(
          adjusted_distance_value = mean_distance_to_gene * gene_density  # multiply by mean distance to adjust for 
        )
      
      # Reformat the result to exclude intermediate columns
      mean_result <- mean_result_per_sp %>%
        select(tree.name, adjusted_distance_value)
      
      if (!x %in% names(all_results)) all_results[[x]] <- list()
      if (!z %in% names(all_results[[x]])) all_results[[x]][[z]] <- list()
      all_results[[x]][[z]] <- mean_result
      
      # Name df after species name and strand type
      assign(paste0(x, "_", "_meanresult"), mean_result)
      
    }
  }, error = function(e) {
    message(paste("Error processing species", x, ":", e$message))
    # Continue with the next iteration
  })
}

# Initialize empty data frame

combined_results <- data.frame()


# Initialize empty list to collect rows
all_df_list <- list()

# Loop through each species
for (species in names(all_results)) {
  tryCatch({
    retro_df <- all_results[[species]][["retro"]]
    dna_df   <- all_results[[species]][["dna"]]
    
    # Safely extract the adjusted values
    retro_val <- if (!is.null(retro_df)) retro_df$adjusted_distance_value else NA
    dna_val   <- if (!is.null(dna_df)) dna_df$adjusted_distance_value else NA
    
    # Create one row per species
    row <- data.frame(
      tree.name = species,
      retro_adjusted_value = retro_val,
      dna_adjusted_value = dna_val
    )
    
    all_df_list[[species]] <- row
  }, error = function(e) {
    message(paste("Error processing species", species, ":", e$message))
  })
}
# Combine into one data frame
final_df <- bind_rows(all_df_list)

# Save the combined results for positive strand to a CSV file
write.csv(final_df, file = "../Spatial_dists_MGE_analyses/MGE_nearest_neighbor_unstranded_class_mrna.csv", row.names = FALSE)

