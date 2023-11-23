raw_to_short <- list(
  'Apilactobacillus_apinorum' = 'A. apinorum',
  'Apilactobacillus_kunkeei' = 'A. kunkeei',
  'Bartonella_apis' = 'Bart. apis',
  'Bifidobacterium_asteroides' = 'B. asteroides',
  'Bifidobacterium_coryneforme_indicum' = 'B. cor./indicum',
  'Bombella_apis' = 'Bom. apis',
  'Bombella_sp' = 'B. sp.',
  'Bombilactobacillus_mellifer' = 'B. mellifer', 
  'Bombilactobacillus_mellis' = 'B. mellis',
  'Commensalibacter_sp' = 'C. sp.',
  'Frischella_perrara' = 'F. perrara',
  'Gilliamella_apicola' = 'G. apicola',
  'Gilliamella_apis' = 'G. apis',
  'Gilliamella_sp' = 'G. sp.',
  'Lactobacillus_apis' = 'L. apis',
  'Lactobacillus_helsingborgensis' = 'L. helsingborgensis',
  'Lactobacillus_kimbladii' = 'L. kimbladii',
  'Lactobacillus_kullabergensis' = 'L. kullabergensis',
  'Lactobacillus_melliventris' = 'L. melliventris',
  'Serratia_marcescens' = 'S. marcescens',
  'Snodgrassella_alvi' = 'S. alvi'
)

cleanup_species_names <- function(names_vec, shorten=FALSE) {
  
  if (shorten) {
    for (raw in names(raw_to_short)) {
      if (raw %in% names_vec) {
        names_vec[which(names_vec == raw)] <- raw_to_short[[raw]]
      }
    }
    
    if ('Bifidobacterium_coryneforme' %in% names_vec) {
      names_vec[which(names_vec == 'Bifidobacterium_coryneforme')] <- raw_to_short[['Bifidobacterium_coryneforme_indicum']]
    }
    
  } else {
    
    names_vec <- gsub('_', ' ', names_vec)
    if (length(grep('Bifidobacterium coryneforme', names_vec)) > 0) {
      names_vec[grep('Bifidobacterium coryneforme', names_vec)] <- 'Bifidobacterium cor./indicum'
    }
    
    names_vec <- gsub('sp$', 'sp.', names_vec)
    
  }
  
  return(names_vec)
  
}


split_multi_category_rows <- function(in_df, category_col, delimiter = ",", num_cores = 1) {
  
  multi_category_row_i <- grep(delimiter, in_df[, category_col])
  
  in_df_nonmulti <- in_df[-multi_category_row_i, , drop = FALSE]
  in_df_multi <- in_df[multi_category_row_i, , drop = FALSE]
  
  in_df_multi_merged_raw <- parallel::mclapply(1:nrow(in_df_multi),
                                               function(i) {
                                                 row_subset <- in_df_multi[i, , drop = FALSE]
                                                 category_split <- base::strsplit(x = row_subset[, category_col], split = ",")[[1]]
                                                 split_row <- row_subset[rep(1, length(category_split)), , drop = FALSE]
                                                 split_row[, category_col] <- category_split
                                                 return(split_row)
                                               },
                                               mc.cores = num_cores)
  
  in_df_multi_merged <- do.call(rbind, in_df_multi_merged_raw)
  
  return(rbind(in_df_nonmulti, in_df_multi_merged))
}

species_plot_colours <- c(
  'pink1', # 'Apilactobacillus apinorum'
  'pink2', # 'Apilactobacillus kunkeei'
  'yellow2', # 'Bartonella apis'
  '#63bfaf', # 'Bifidobacterium asteroides'
  '#62d07a', # 'Bifidobacterium cor./indicum'
  'cyan', # 'Bombella apis'
  'cyan2', # 'Bombella sp.'
  '#d56c24', # 'Bombilactobacillus mellifer'
  '#bf4d81', # 'Bombilactobacillus mellis'
  '#d2cb3b', # 'Commensalibacter sp.'          
  'brown', #'Frischella perrara'
  'royalblue1', # 'Gilliamella apicola'
  'skyblue1', # 'Gilliamella apis'
  'mediumblue', # 'Gilliamella sp.'
  'grey40', # 'Lactobacillus apis'            
  'grey50', #'Lactobacillus helsingborgensis'
  'grey60', # 'Lactobacillus kimbladii'
  'grey70', # 'Lactobacillus kullabergensis'
  'grey80', # 'Lactobacillus melliventris'
  'black', # 'Serratia marcescens'           
  '#9f4ad4' # 'Snodgrassella alvi'
)

names(species_plot_colours) <- c('Apilactobacillus_apinorum',
                                 'Apilactobacillus_kunkeei',
                                 'Bartonella_apis',
                                 'Bifidobacterium_asteroides',
                                 'Bifidobacterium_coryneforme_indicum',
                                 'Bombella_apis',
                                 'Bombella_sp',
                                 'Bombilactobacillus_mellifer',
                                 'Bombilactobacillus_mellis',
                                 'Commensalibacter_sp',
                                 'Frischella_perrara',
                                 'Gilliamella_apicola',
                                 'Gilliamella_apis',
                                 'Gilliamella_sp',
                                 'Lactobacillus_apis',
                                 'Lactobacillus_helsingborgensis',
                                 'Lactobacillus_kimbladii',
                                 'Lactobacillus_kullabergensis',
                                 'Lactobacillus_melliventris',
                                 'Serratia_marcescens',
                                 'Snodgrassella_alvi')
