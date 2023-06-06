# Function to add pseudocount to all samples and apply centred log-ratio (CLR) transformation.
# Note that this functions assumes that columns are samples and rows are features.

pseudocount_and_clr <- function(in_df, pseudocount) {
  return(as.data.frame(apply(in_df + pseudocount,
                             2,
                             function(x){log(x) - mean(log(x))})))
}