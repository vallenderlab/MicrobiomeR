## Initialize a list for storing heat tree plots
get_heat_tree_plots <- function(obj, database = NULL, load_flag = NULL, filter_flag = NULL, save_flag = NULL, ambiguous_filter = NULL,
                                kingdom_filter = NULL, rank_list = NULL) {
  if (is.null(rank_list)) {
    rank_list <- ranks
  }
  htrees <- list()
  # Create a metacoder object from a phyloseq/metacoder/RData file
  metacoder_object <- object_handler(obj)
  # Create a list of heat_tree plots for saving
  for (rank in rank_list) {
    rank_level <- rank_index[[rank]]
    smb_data <- metacoder_object %>% filter_obs(data = c("statistical_data", "taxa_proportions", "taxa_abundance", "stats_tax_data"), n_supertaxa < rank_level, drop_taxa = TRUE)
    title <- sprintf("Bacterial Abundance (%s Level)", rank)
    default_heat_tree_parameters <- get_default_parameters(smb_data = smb_data, title = title, func_names = c("heat_tree"))
    # Filter by Taxonomy Rank and then create a heat tree.
    htrees[[rank]] <- do.call(what = "heat_tree", args = c(default_heat_tree_parameters))
    # Made plot title centered
    htrees[[rank]] <- htrees[[rank]] + theme(
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 24, family = "Arial")
    ) + ggtitle(title)
    htrees[[rank]]
  }
  htrees[["metacoder_object"]] <- smb_data
  return(htrees)
}
