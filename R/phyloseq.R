#' @title Get Phyloseq Object
#'
#' @description Import or Create a phyloseq object
#'
#' @param b_file DESCRIPTION.
#' @param t_file DESCRIPTION.
#' @param md_file DESCRIPTION.
#' @param parse_func DESCRIPTION.
#' @param rdata_file DESCRIPTION.
#'
#' @return A phyloseq object
#' @export
get_phyloseq_obj <- function(b_file=NULL, t_file=NULL, md_file=NULL, parse_func=NULL, rdata_file=NULL) {
  if (!is.null(rdata_file)) {
    load(rdata_file)
    phyloseq_obj <- get("phyloseq_obj")
    return(phyloseq_obj)
  } else {
    # Load the biom file and the Full tree file into a phyloseq object
    phyloseq_biom <- import_biom(BIOMfilename = b_file, parseFunction = parse_func, treefilename = t_file)
    # Write to a tree file that subsets our data
    if (!is.null(t_file)) {
      phyloseq_tree <- phy_tree(phyloseq_biom)
      # If not rooted, then root the tree using the longest edge(s) as an ougroup
      if (!ape::is.rooted(phyloseq_tree)) {
        phyloseq_tree <- root_by_longest_edge(phyloseq_tree)
      }
      # mkdir("output")
      # write.tree(biom_tree, "output/OTU_table.tre")
    }
    # Extract metadata from the phyloseq object:
    phyloseq_metadata <- import_qiime_sample_data(md_file)

    # Create a phyloseq object with all of our data
    phyloseq_obj <- merge_phyloseq(phyloseq_biom, phyloseq_metadata)
    if (!is.null(t_file)) {
      phy_tree(phyloseq_obj) <- root_by_longest_edge(phy_tree(phyloseq_obj))
    }
    return(phyloseq_obj)
  }
}


#' @title Pick New Outgroup
#'
#' @description The tree rooting functions inspired by recent comments in https://github.com/joey711/phyloseq/issues/597
#'
#' @param unrooted_tree DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
pick_new_outgroup <- function(unrooted_tree) {
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(unrooted_tree$edge),
      data.table(length = unrooted_tree$edge.length)
    )[1:Ntip(unrooted_tree)] %>%
    cbind(data.table(id = unrooted_tree$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}


#' @title Root by Longest Edge
#'
#' @description FUNCTION_DESCRIPTION
#'
#' @param unrooted_tree DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
root_by_longest_edge <- function(unrooted_tree) {
  new.outgroup <- pick_new_outgroup(unrooted_tree)
  rootedTree <- ape::root(unrooted_tree, outgroup = new.outgroup, resolve.root = TRUE)
  return(rootedTree)
}

#' @title Parse SILVA Taxonomy
#'
#' @description This function parses SILVA data.
#'
#' @param char.vec DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
parse_taxonomy_silva <- function(char.vec) {
  # Use default to assign names to elements in case problem with greengenes prefix
  char.vec <- parse_taxonomy_default(char.vec)
  # Check for unassigned taxa
  if (char.vec["Rank1"] == "Unassigned") {
    char.vec <- c(
      Rank1 = "D_0__Unassigned", Rank2 = "D_1__Unassigned", Rank3 = "D_2__Unassigned", Rank4 = "D_3__Unassigned",
      Rank5 = "D_4__Unassigned", Rank6 = "D_5__Unassigned", Rank7 = "D_6__Unassigned"
    )
  }
  # Define the meaning of each prefix according to GreenGenes taxonomy
  Tranks <- c(D_0 = "Kingdom", D_1 = "Phylum", D_2 = "Class", D_3 = "Order", D_4 = "Family", D_5 = "Genus", D_6 = "Species")
  # Check for prefix using regexp, warn if there were none. trim indices, ti
  ti <- grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
  if (length(ti) == 0L) {
    warning(
      "No silva prefixes were found. \n",
      "Consider using parse_taxonomy_delfault() instead if true for all OTUs. \n",
      "Dummy ranks may be included among taxonomic ranks now."
    )
    # Will want to return without further modifying char.vec
    taxvec <- char.vec
    # Replace names of taxvec according to prefix, if any present...
  } else {
    # Format character vectors for Ambiguous taxa
    if (length(ti) < 7) {
      for (key in names(char.vec)) {
        if (char.vec[key] == "Ambiguous_taxa") {
          tax_no <- (as.numeric(substr(key, 5, 5)) - 1)
          char.vec[key] <- sprintf("D_%s__Ambiguous_taxa", tax_no)
        }
      }
      # Reset the trimmed indicies if Ambiguous taxa
      ti <- grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
    }
    # Remove prefix using sub-"" regexp, call result taxvec
    taxvec <- gsub("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", "", char.vec)
    # Define the ranks that will be replaced
    repranks <- Tranks[substr(char.vec[ti], 1, 3)]
    # Replace, being sure to avoid prefixes notK present in Tranks
    names(taxvec)[ti[!is.na(repranks)]] <- repranks[!is.na(repranks)]
  }
  return(taxvec)
}

#' @title Parse Taxonomy Greengenes2
#'
#' @description  A custom Greengenes parsing function
#'
#' @param char.vec DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
parse_taxonomy_greengenes2 <- function(char.vec) {
  # Use default to assign names to elements in case problem with greengenes prefix
  char.vec <- parse_taxonomy_default(char.vec)
  # Check for unassigned taxa
  if (char.vec["Rank1"] == "Unassigned") {
    char.vec <- c(
      Rank1 = "k__Unassigned", Rank2 = "p__Unassigned", Rank3 = "c__Unassigned", Rank4 = "o__Unassigned",
      Rank5 = "f__Unassigned", Rank6 = "g__Unassigned", Rank7 = "s__Unassigned"
    )
  }
  # Define the meaning of each prefix according to GreenGenes taxonomy
  Tranks <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
  # Check for prefix using regexp, warn if there were none. trim indices, ti
  ti <- grep("[[:alpha:]]{1}\\_\\_", char.vec)
  if (length(ti) == 0L) {
    warning(
      "No greengenes prefixes were found. \n",
      "Consider using parse_taxonomy_default() instead if true for all OTUs. \n",
      "Dummy ranks may be included among taxonomic ranks now."
    )
    # Will want to return without further modifying char.vec
    taxvec <- char.vec
    # Replace names of taxvec according to prefix, if any present...
  } else {
    # Remove prefix using sub-"" regexp, call result taxvec
    taxvec <- gsub("[[:alpha:]]{1}\\_\\_", "", char.vec)
    # Define the ranks that will be replaced
    repranks <- Tranks[substr(char.vec[ti], 1, 1)]
    # Replace, being sure to avoid prefixes not present in Tranks
    names(taxvec)[ti[!is.na(repranks)]] <- repranks[!is.na(repranks)]
  }
  return(taxvec)
}

#' @title OTU taxa to dataframe
#'
#' @description A function for manipulating a phyloseq object into a list containing the following tibbles:
#' Tax table, OTU table, Metadata table, Stressed/Control Samples, OTU names, Merged Tax/OTU table
#'
#' @param phyloseq_obj DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
otutax_to_dataframe <- function(phyloseq_obj) {
  tax_otu <- list()
  tax_otu$tax <- data.frame(tax_table(phyloseq_obj))
  tax_otu$otu <- data.frame(otu_table(phyloseq_obj))
  tax_otu$sam <- data.frame(sample_data(phyloseq_obj))
  tax_otu$stressed_samples <- rownames(tax_otu$sam[tax_otu$sam$TreatmentGroup == "Stressed", ])
  tax_otu$control_samples <- rownames(tax_otu$sam[tax_otu$sam$TreatmentGroup == "Control", ])
  tax_otu$sam <- as.tibble(tax_otu$sam)
  tax_otu$otu_names <- rownames(tax_otu$tax)
  tax_otu$tax <- as.tibble(mutate(tax_otu$tax, OTU = rownames(tax_otu$tax)))
  tax_otu$otu <- as.tibble(mutate(tax_otu$otu, OTU = rownames(tax_otu$otu)))
  tax_otu$data <- as.tibble(left_join(x = tax_otu$otu, y = tax_otu$tax, by = "OTU"))
  return(tax_otu)
}

#' @title Remove Ambiguous Taxa
#'
#' @description A function for removing named ambiguous taxonomy at specific rank(s)
#'
#' @param phyloseq_obj DESCRIPTION.
#' @param ranks DESCRIPTION.
#' @param ambiguous_names DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
remove_ambiguous_taxa <- function(phyloseq_obj, ranks, ambiguous_names) {
  tax_otu <- otutax_to_dataframe(phyloseq_obj)
  pdat <- tax_otu$data
  for (r in ranks) {
    if (r == "Kingdom") {
      pdat <- dplyr::filter(pdat, !Kingdom %in% ambiguous_names)
    } else if (r == "Phylum") {
      pdat <- dplyr::filter(pdat, !Phylum %in% ambiguous_names)
    } else if (r == "Class") {
      pdat <- dplyr::filter(pdat, !Class %in% ambiguous_names)
    } else if (r == "Order") {
      pdat <- dplyr::filter(pdat, !Order %in% ambiguous_names)
    } else if (r == "Famlily") {
      pdat <- dplyr::filter(pdat, !Family %in% ambiguous_names)
    } else if (r == "Genus") {
      pdat <- dplyr::filter(pdat, !Genus %in% ambiguous_names)
    } else if (r == "Species") {
      pdat <- dplyr::filter(pdat, !Speices %in% ambiguous_names)
    }
  }

  tax <- dplyr::select(pdat, colnames(tax_otu$tax)) %>%
    as.data.frame() %>%
    column_to_rownames("OTU") %>%
    as.matrix()

  otu <- dplyr::select(pdat, colnames(tax_otu$otu)) %>%
    as.data.frame() %>%
    column_to_rownames("OTU") %>%
    as.matrix()
  # Merge into phyloseq data;  Include tree if applicable
  if (inherits(phyloseq_obj, "taxonomyTable")) {
    return(tax_table(tax))
  } else {
    err_ <- try(phy_tree(phyloseq_obj), silent = TRUE)
    if (stringr::str_detect(err_, "Error")) {
      phyloseq_obj <- merge_phyloseq(tax_table(tax), otu_table(otu, taxa_are_rows = TRUE), sample_data(phyloseq_obj))
    } else {
      phyloseq_obj <- merge_phyloseq(tax_table(tax), phy_tree(phyloseq_obj), otu_table(otu, taxa_are_rows = TRUE), sample_data(phyloseq_obj))
    }
    return(phyloseq_obj)
  }
}

#' @title Preprocess Phyloseq
#'
#' @description The preprocessing function
#'
#' @param phyloseq_object DESCRIPTION.
#' @param process_list DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
preprocess_phyloseq <- function(phyloseq_object, process_list = NULL) {
  pp_p <- phyloseq_object
  if (is.null(process_list)) {
    params <- get_default_parameters(smb_data = phyloseq_object, func_names = c("preprocess_phyloseq"))
    process_list <- params$process_list
  }
  pr_l <- process_list

  ## REMOVAL
  # Remove empty samples
  pp_p <- prune_taxa(taxa_sums(pp_p) > 0, pp_p)

  ## FILTER
  for (proc in names(process_list)) {
    # print(proc)
    # print(process_list[[proc]])
    if (proc == "master_thresh") {
      # Filter anything at a certain threshhold
      # phyloseq vignette - master_thresh = 1e-5
      pp_p <- phyloseq::filter_taxa(pp_p, function(x) {
        mean(x) > pr_l[[proc]]
      }, TRUE)
    }
    if (proc == "taxon_filter" & length(pr_l[[proc]]) == 3) {
      # Taxonomic Prevalence filtering
      # Filter OTUs that do not appear more than a certian amount of times in a certain percentage of samples
      # at the specified agglomerated rank
      # phyloseq vignette - taxa_thresh = c(5, 0.5)
      taxon_filter <- pr_l[[proc]]

      for (tf in taxon_filter) {
        tf_rank <- tf$rank
        tf_min_abund <- taxon_filter$min_a
        tf_req_samp_perc <- taxon_filter$r_s_p
        tf_glom <- tax_glom(pp_p, tf_rank, NArm = FALSE)
        tf_unfiltered <- genefilter_sample(tf_glom, filterfun_sample(function(x) x > tf_min_abund),
                                           A = tf_req_samp_perc * nsamples(tf_glom)
        )
        # tf_filtered <- names(tf_unfiltered[tf_unfiltered==TRUE])

        p_filter <- prune_taxa(tf_unfiltered, tf_glom)

        # Formating for standard evaluation
        r_filter <- get_taxa_unique(p_filter, taxonomic.rank = tf_rank)
        if (tf_rank == "Kingdom") {
          pp_p <- subset_taxa(pp_p, Kingdom %in% r_filter)
        } else if (tf_rank == "Phylum") {
          pp_p <- subset_taxa(pp_p, Phylum %in% r_filter)
        } else if (tf_rank == "Class") {
          pp_p <- subset_taxa(pp_p, Class %in% r_filter)
        } else if (tf_rank == "Order") {
          pp_p <- subset_taxa(pp_p, Order %in% r_filter)
        } else if (tf_rank == "Famlily") {
          pp_p <- subset_taxa(pp_p, Family %in% r_filter)
        } else if (tf_rank == "Genus") {
          pp_p <- subset_taxa(pp_p, Genus %in% r_filter)
        } else if (tf_rank == "Species") {
          pp_p <- subset_taxa(pp_p, Speices %in% r_filter)
        }
      }
    }
    if (proc == "prevelance_filter" & length(pr_l[[proc]]) == 2) {
      # Prevelance filtering
      # Filter OTUs that do not appear more than a certian amount of times in a certain percentage of samples.
      # phyloseq vignette - taxa_thresh = c(5, 0.5)
      prevalence_filter <- pr_l[[proc]]
      min_abund <- prevalence_filter[1]
      required_sample_percentage <- prevalence_filter[2]
      tx <- genefilter_sample(pp_p, filterfun_sample(function(x) x > min_abund),
                              A = required_sample_percentage * nsamples(pp_p)
      )
      pp_p <- prune_taxa(tx, pp_p)
    }
    if (proc == "glom_rank") {
      pp_p <- tax_glom(pp_p, taxrank = pr_l[[proc]])
    }
    if (proc == "ambiguous") {
      ambiguous <- pr_l[[proc]]
      amb_ranks <- ambiguous$amb_ranks
      amb_items <- ambiguous$amb_items

      pp_p <- remove_ambiguous_taxa(pp_p, amb_ranks, amb_items)
    }
    ## Standardize
    if (proc == "coeff_of_variation") {
      # Standardize abundances to the median sequencing depth
      total <- median(sample_sums(pp_p))
      standf <- function(x, t=total) round(t * (x / sum(x)))
      pp_trans <- transform_sample_counts(pp_p, standf)
      # Filter the taxa using the coefficient of variation
      # phyloseq vignette - coeff_of_variation = 0.3
      pp_filter <- phyloseq::filter_taxa(pp_trans, function(x) sd(x) / mean(x) > pr_l[[proc]])
      pp_p <- prune_taxa(pp_filter, pp_p)
    }

    ## TRANSFORM
    if (proc == "transform" & class(pr_l[[proc]]) == "function") {
      # Transform to get relative abundance
      # Note:  always do this before removing taxa that contribute to rel abund
      transform <- pr_l[[proc]]
      pp_p <- transform_sample_counts(pp_p, transform)
    }
    if (proc == "merge_samples") {
      sample_group <- pr_l[[proc]][["group"]]
      func <- pr_l[[proc]][["func"]]
      pp_p <- merge_samples(pp_p, group = sample_group, fun = func)
    }
  }
  ## REMOVAL
  # Remove empty samples
  pp_p <- prune_taxa(taxa_sums(pp_p) > 0, pp_p)
  return(pp_p)
}

#' @title Get Distance Methods
#'
#' @description The function for manipulating phyloseqs native distanceMethodList
#'
#' @param tree_methods DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
get_distance_methods <- function(tree_methods=TRUE) {
  ### Create list of distance methods
  dist_methods <- unlist(distanceMethodList)
  if (!tree_methods) {
    ### Remove the two distance-methods that require a tree, and the generic custom method that requires
    ### user-defined distance arguments.
    dist_methods <- dist_methods[-(1:3)]
  }
  ### Remove the user-defined distance
  dist_methods <- dist_methods[-which(dist_methods == "ANY")]
  return(dist_methods)
}
