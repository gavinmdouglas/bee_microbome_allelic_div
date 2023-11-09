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