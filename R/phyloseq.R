#' @title Get Phyloseq Object
#' @description Create a phyloseq object using a biom file, phylogenetic tree file, and a
#'  metadata file. Alternatively, an .Rdata file can be used.  For now this funciton takes
#'  output from Qiime.  The data provided by this package was processed by the NIH's Nephele
#'  pipeline.
#' @param biom_file A file in the BIOM format. Default: NULL
#' @param tree_file A Phylogenetic tree file. Default: NULL
#' @param metadata_file Sample metadata in tab delimited format. Default: NULL
#' @param parse_func The parse function used to parse taxonomy strings from
#'  Greengenes or SILVA database. Default: NULL
#' @param rdata_file A .Rdata file.  Default: NULL
#' @param path A path for the output tree file.  Default: NULL
#' @return A phyloseq object, and a phylogenetic tree file if one does not already exist.
#' @pretty_print TRUE
#' @details This function heavily relys on the phyloseq package to import data into R.
#'  It also requires you to use an absolute path to your data files or for your data files
#'  to be in the working directory.
#' @examples
#' \dontrun{
#'
#' > biom_file <- "input_data/silva_OTU.biom"
#' > md_file <- "input_data/nephele_metadata.txt"
#' > phy_obj <- get_phyloseq_obj(biom_file=biom_file, metadata_file=md_file)
#' }
#' @export
#' @family Data Importers
#' @rdname get_phyloseq_obj
#' @seealso
#'  \code{\link[phyloseq]{import_biom}}, \code{\link[phyloseq:phy_tree-methods]{phyloseq::phy_tree()}}, \code{\link[phyloseq]{import_qiime_sample_data}}, \code{\link[phyloseq]{merge_phyloseq}}
#'  \code{\link[ape]{root}}
#'  \code{\link[MicrobiomeR]{root_by_longest_edge}}

#' @importFrom phyloseq import_biom phy_tree import_qiime_sample_data merge_phyloseq
#' @importFrom ape is.rooted write.tree
get_phyloseq_obj <- function(biom_file = NULL, tree_file = NULL, metadata_file = NULL, parse_func = NULL, rdata_file = NULL, path = NULL) {
  if (!is.null(rdata_file)) {
    load(rdata_file)
    phyloseq_obj <- get("phyloseq_obj")
    return(phyloseq_obj)
  } else {
    # Load the biom file and the Full tree file into a phyloseq object
    phyloseq_biom <- phyloseq::import_biom(BIOMfilename = biom_file, parseFunction = parse_func, treefilename = tree_file)
    # Write to a tree file that subsets our data
    if (!is.null(tree_file)) {
      phyloseq_tree <- phyloseq::phy_tree(phyloseq_biom)
      # If not rooted, then root the tree using the longest edge(s) as an ougroup
      if (!ape::is.rooted(phyloseq_tree)) {
        ape_tree <- MicrobiomeR::root_by_longest_edge(phyloseq_tree)
      }
      tree_path <- mkdir(dirname = "output", path = path)
      ape::write.tree(ape_tree, sprintf("%s/rooted_%s", tree_path, tree_file))
    }
    # Extract metadata from the phyloseq object:
    phyloseq_metadata <- phyloseq::import_qiime_sample_data(metadata_file)

    # Create a phyloseq object with all of our data
    phyloseq_obj <- phyloseq::merge_phyloseq(phyloseq_biom, phyloseq_metadata)
    if (!is.null(tree_file)) {
      phy_tree(phyloseq_obj) <- MicrobiomeR::root_by_longest_edge(phy_tree(phyloseq_obj))
    }
    return(phyloseq_obj)
  }
}



#' @title Pick Outgroup for Tree
#' @description Pick an outgroup for rooting a phylogenetic tree.
#' @param unrooted_tree An unrooted tree object.
#' @return A new tree with the longest bracnh as the outgroup.
#' @pretty_print TRUE
#' @details This function preprocess a phylogenetic tree for rooting by the longest edge.
#' Please see this issue for more details \url{https://github.com/joey711/phyloseq/issues/597}.
#' @export
#' @family Phylogenetic Tree Manipulators
#' @rdname pick_new_outgroup
#' @seealso
#'  \code{\link[data.table]{data.table-package}}
#'
#'  \code{\link[ape]{summary.phylo}}
#' @importFrom data.table data.table
#' @importFrom ape Ntip
pick_new_outgroup <- function(unrooted_tree) {
  # Set up data
  tree_DT <- data.table::data.table(unrooted_tree$edge)
  tree_DT_len <- data.table::data.table(length = unrooted_tree$edge.length)
  tree_DT_id <- data.table::data.table(id = unrooted_tree$tip.label)
  tree_tips <- ape::Ntip(unrooted_tree)
  # tablify parts of tree that we need.
  treeDT <- cbind(tree_DT,tree_DT_len)[1:tree_tips] %>% cbind(tree_DT_id)
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}



#' @title Root Phylogenetic Tree
#' @description This function roots a phylogenetic tree object by it's longest edge.
#' @param unrooted_tree An unrooted phylogenetic tree object.
#' @return A rooted phlogenetic tree object.
#' @pretty_print TRUE
#' @details Please see this issue for more details \url{https://github.com/joey711/phyloseq/issues/597}.
#' @export
#' @family Phylogenetic Tree Manipulators
#' @rdname root_by_longest_edge
#' @seealso
#'  \code{\link[ape]{root}}
#' @importFrom ape root
root_by_longest_edge <- function(unrooted_tree) {
  new.outgroup <- pick_new_outgroup(unrooted_tree)
  rootedTree <- ape::root(unrooted_tree, outgroup = new.outgroup, resolve.root = TRUE)
  return(rootedTree)
}



#' @title Parse elements of a taxonomy vector
#' @description These are provided as both example and default functions for
#' parsing a character vector of taxonomic rank information for a single taxa.
#' As default functions, these are intended for cases where the data adheres to
#' the naming convention used by greengenes
#' the naming convention used by greengenes and silva.
#' (\url{http://greengenes.lbl.gov/cgi-bin/nph-index.cgi})
#' or where the convention is unknown, respectively.
#' To work, these functions -- and any similar custom function you may want to
#' create and use -- must take as input a single character vector of taxonomic
#' ranks for a single OTU, and return a \strong{named} character vector that has
#' been modified appropriately (according to known naming conventions,
#' desired length limits, etc.
#' The length (number of elements) of the output named vector does \strong{not}
#' need to be equal to the input, which is useful for the cases where the
#' source data files have extra meaningless elements that should probably be
#' removed, like the ubiquitous
#' ``Root'' element often found in greengenes/QIIME taxonomy labels.
#' In the case of \code{parse_taxonomy_default}, no naming convention is assumed and
#' so dummy rank names are added to the vector.
#' More usefully if your taxonomy data is based on greengenes, the
#' \code{parse_taxonomy_greengenes} function clips the first 3 characters that
#' identify the rank, and uses these to name the corresponding element according
#' to the appropriate taxonomic rank name used by greengenes
#' (e.g. \code{"p__"} at the beginning of an element means that element is
#' the name of the phylum to which this OTU belongs).
#' If you taxonomy data is based on SILVA, the \code{parse_taxonomy_silva_128} function
#' clips the first 5 characters that identify rank, and uses these to name the
#' corresponding element according to the appropriate taxonomic rank name used
#' by SILVA (e.g. \code{"D_1__"} at the beginning of an element means that element
#' is the name of the phylum to which this OTU belongs.
#' Alternatively you can create your own function to parse this data.
#' Most importantly, the expectations for these functions described above
#' make them compatible to use during data import,
#' specifcally the \code{\link{import_biom}} function, but
#' it is a flexible structure that will be implemented soon for all phyloseq
#' import functions that deal with taxonomy (e.g. \code{\link{import_qiime}}).
#' @param char.vec (Required). A single character vector of taxonomic
#'  ranks for a single OTU, unprocessed (ugly).
#' @return A character vector in which each element is a different
#'  taxonomic rank of the same OTU, and each element name is the name of
#'  the rank level. For example, an element might be \code{"Firmicutes"}
#'  and named \code{"phylum"}.
#'  These parsed, named versions of the taxonomic vector should
#'  reflect embedded information, naming conventions,
#'  desired length limits, etc; or in the case of \code{\link{parse_taxonomy_default}},
#'  not modified at all and given dummy rank names to each element.
#' @pretty_print TRUE
#' @details This function is currently under PR review by phyloseq in a well supported
#' pull request: \url{https://github.com/joey711/phyloseq/pull/854}.  If you use this function,
#' then please comment on the GitHub PR to encourage merging this feature.
#' @examples \dontrun{
#'
#'  > taxvec1 = c("Root", "k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae")
#'  > parse_taxonomy_default(taxvec1)
#'  > parse_taxonomy_greengenes(taxvec1)
#'  > taxvec2 = c("Root;k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Staphylococcaceae")
#'  > parse_taxonomy_qiime(taxvec2)
#'  > taxvec3 = c("D_0__Bacteria", "D_1__Firmicutes", "D_2__Bacilli", "D_3__Staphylococcaceae")
#'  > parse_taxonomy_silva_128(taxvec3)
#'  }
#' @export
#' @usage parse_taxonomy_default(char.vec)
#' parse_taxonomy_greengenes(char.vec)
#' parse_taxonomy_silva_128(char.vec)
#' parse_taxonomy_qiime(char.vec)
#' @rdname parse_taxonomy_silva_128
#' @seealso
#'  \code{\link[phyloseq]{parse_taxonomy_default}}
#'
#'  \code{\link[phyloseq]{parse_taxonomy_greengenes}}
#'
#'  \code{\link[phyloseq]{parse_taxonomy_qiime}}
#'
#'  \code{\link[phyloseq]{import_biom}}
#'
#'  \code{\link[phyloseq]{import_qiime}}
parse_taxonomy_silva_128 <- function(char.vec) {
  # Use default to assign names to elements in case problem with greengenes prefix
  char.vec <- phyloseq::parse_taxonomy_default(char.vec)
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

#'  @title Parse Taxonomy Greengenes2
#'
#'  @description  A custom Greengenes parsing function
#'
#'  @param char.vec DESCRIPTION.
#'
#'  @return RETURN_DESCRIPTION
#'  @export
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




#' @title Convert Phyloseq Objects to Tibbles
#' @description This function converts a phyloseq object to a tibble, while
#' taking into account the treatment groups found in the sample metadata.
#' @return A list of tibbles containing all of the phyloseq data in a tidy format.
#' @pretty_print TRUE
#' @details MicrobiomeR was developed for processing Qiime data delivered by
#' the NIH's Nephele pipeline (\url{https://nephele.niaid.nih.gov/}).   MicrobiomeR
#' was also written using the tidyverse in order to make the workflow reproducible.
#' This function was written with that data in mind.  In the future this function
#' can be expanded to accept other data models.  The returned list will contain
#' a tax tibble, an otu tibble, a sample data tibble, a full data tibble, a
#' vector of otu names, and a vectors containing the names of treatment group
#' samples.
#'
#' @examples
#' \dontrun{
#'
#' > phy_obj <- get_phyloseq_object(...)
#' > phy_tibble <- phyloseq_to_tibble(phy_obj, treatment_groups = list("Stressed", "Control"))
#' > phy_tibble
#' $tax
#' # A tibble: 3 x 8
#'   Kingdom Phylum       Class         Order            Family           Genus           Species          OTU
#'   <fct>   <fct>        <fct>         <fct>            <fct>            <fct>           <fct>            <chr>
#' 1 Archaea Euryarchaeo~ Methanobacte~ Methanobacteria~ Methanobacteria~ Methanobreviba~ uncultured arch~ JX522720.1.1265
#' 2 Archaea Euryarchaeo~ Methanobacte~ Methanobacteria~ Methanobacteria~ Methanobreviba~ uncultured arch~ ADJS01014867.75~
#' 3 Archaea Euryarchaeo~ Methanobacte~ Methanobacteria~ Methanobacteria~ Methanobreviba~ Ambiguous_taxa   AB906006.1.1368
#'
#' $otu
#' # A tibble: 3 x 49
#'   Sample_1 Sample_2 Sample_3 Sample_6 Sample_9 Sample_10 Sample_13 Sample_14 Sample_17 Sample_21 Sample_22 Sample_24
#'      <dbl>    <dbl>    <dbl>    <dbl>    <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
#' 1  0.             0  0.             0  0.       0          8.92e-6         0 0.0000154 0.0000268 0.0000851 0
#' 2  0.             0  1.99e-5        0  0.       0          2.68e-5         0 0.0000154 0.0000625 0.0000638 0.000121
#' 3  8.14e-5        0  1.15e-3        0  1.11e-5  0.000137   1.78e-5         0 0.0000154 0.0000179 0.000106  0.0000404
#' #   Sample_56 <dbl>, Sample_58 <dbl>, Sample_8 <dbl>, Sample_4 <dbl>, Sample_52 <dbl>, Sample_36 <dbl>, OTU <chr>
#'
#' $sam
#' # A tibble: 48 x 8
#'    X.SampleID BarcodeSequence LinkerPrimerSequ~ ForwardFastqFile  ReverseFastqFile TreatmentGroup SampleName Description
#'  * <fct>      <lgl>           <lgl>             <fct>             <fct>            <fct>               <int> <fct>
#'  1 Sample_1   NA              NA                33749_S1_L001_R1~ 33749_S1_L001_R~ Stressed            33749 Fecal
#'  2 Sample_2   NA              NA                33739_S2_L001_R1~ 33739_S2_L001_R~ Stressed            33739 Fecal
#'  3 Sample_3   NA              NA                33737_S3_L001_R1~ 33737_S3_L001_R~ Stressed            33737 Fecal
#' # ... with 38 more rows
#'
#' $stressed_samples
#'  [1] "Sample_1"  "Sample_2"  "Sample_3"  "Sample_6"  "Sample_9"  "Sample_10" "Sample_13" "Sample_14" "Sample_17"
#' [10] "Sample_21" "Sample_22" "Sample_24" "Sample_25" "Sample_26" "Sample_30" "Sample_59" "Sample_7"  "Sample_27"
#' [19] "Sample_5"  "Sample_20" "Sample_29" "Sample_11" "Sample_12" "Sample_15" "Sample_16" "Sample_18" "Sample_19"
#' [28] "Sample_23" "Sample_28" "Sample_31" "Sample_8"  "Sample_4"
#'
#' $control_samples
#'  [1] "Sample_50" "Sample_60" "Sample_61" "Sample_45" "Sample_57" "Sample_35" "Sample_40" "Sample_41" "Sample_46"
#' [10] "Sample_47" "Sample_51" "Sample_55" "Sample_56" "Sample_58" "Sample_52" "Sample_36"
#'
#' $otu_names
#' [1] "JX522720.1.1265"      "ADJS01014867.75.1537" "AB906006.1.1368"
#'
#' $data
#' # A tibble: 3 x 56
#'   Sample_1 Sample_2 Sample_3 Sample_6 Sample_9 Sample_10 Sample_13 Sample_14 Sample_17 Sample_21 Sample_22 Sample_24
#'      <dbl>    <dbl>    <dbl>    <dbl>    <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
#' 1  0.             0  0.             0  0.       0          8.92e-6         0 0.0000154 0.0000268 0.0000851 0
#' 2  0.             0  1.99e-5        0  0.       0          2.68e-5         0 0.0000154 0.0000625 0.0000638 0.000121
#' 3  8.14e-5        0  1.15e-3        0  1.11e-5  0.000137   1.78e-5         0 0.0000154 0.0000179 0.000106  0.0000404
#' # ... with 44 more variables: Sample_25 <dbl>, Sample_26 <dbl>, Sample_30 <dbl>, Sample_50 <dbl>, Sample_59 <dbl>,
#' #   Sample_56 <dbl>, Sample_58 <dbl>, Sample_8 <dbl>, Sample_4 <dbl>, Sample_52 <dbl>, Sample_36 <dbl>, OTU <chr>,
#' #   Kingdom <fct>, Phylum <fct>, Class <fct>, Order <fct>, Family <fct>, Genus <fct>, Species <fct>
#'
#'  }
#' @export
#' @family Data Manipulators
#' @inherit phyloseq_to_dataframe params
#' @rdname phyloseq_to_tibble
#' @seealso
#'  \code{\link[phyloseq:tax_table-methods]{phyloseq::tax_table()}}, \code{\link[phyloseq:otu_table-methods]{phyloseq::otu_table()}}, \code{\link[phyloseq:sample_data-methods]{phyloseq::sample_data()}}
#'
#'  \code{\link[stringr]{case}}, \code{\link[tibble]{as_tibble}}, \code{\link[dplyr]{join}}
#' @importFrom phyloseq tax_table otu_table sample_data
#' @importFrom stringr str_to_lower
#' @importFrom tibble as.tibble
#' @importFrom dplyr left_join
phyloseq_to_tibble <- function(phyloseq_obj, treatment_groups) {
  df_tax_otu <- phyloseq_to_dataframe(phyloseq_obj = phyloseq_obj, treatment_groups = treatment_groups)
  df_tax_otu$sam <- tibble::as.tibble(df_tax_otu$sam)
  df_tax_otu$tax <- tibble::as.tibble(df_tax_otu$tax)
  df_tax_otu$otu <- tibble::as.tibble(df_tax_otu$otu)
  df_tax_otu$data <- tibble::as.tibble(df_tax_otu$data)
  return(df_tax_otu)
}





#' @title Convert Phyloseq Object to Statistical Data Frame
#' @description This function takes a phyloseq object and converts it to a list of
#' dataframes followed by analyzing the treatment groups and adding statistical data to
#' the same list.
#' @return A list of dataframes containing all of the phyloseq data with additional statistical data.
#' @pretty_print TRUE
#' @examples
#' \dontrun{
#' > phy_obj <- get_phyloseq_object(...)
#' > stats_df <- phyloseq_to_stats_dataframe(phy_obj, treatment_groups = list("Stressed", "Control"))
#' > stats_df$stats
#' TODO: Add a stats dataframe here
#' }
#' @export
#' @family Data Manipulators
#' @inherit phyloseq_to_tibble params seealso
#' @rdname phyloseq_to_stats_dataframe
phyloseq_to_stats_dataframe <- function(phyloseq_obj, treatment_groups) {
  df_tax_otu <- phyloseq_to_dataframe(phyloseq_obj = phyloseq_obj, treatment_groups = treatment_groups)
  rownames(df_tax_otu$data) <- df_tax_otu$data$OTU
  # TODO: Remove Stress/Control from syntax.  Make dynamic for other users.
  stressed_samples <- df_tax_otu$stressed_samples
  control_samples <- df_tax_otu$control_samples
  temp <- df_tax_otu$data[, !(colnames(df_tax_otu$data) %in% c("OTU"))][c(stressed_samples, control_samples)]
  wp_values <- apply(as.matrix(temp), 1, function(x) wilcox.test(x[stressed_samples], x[control_samples])$p.value)
  mean_stressed <- apply(as.matrix(temp), 1, function(x) mean(x[stressed_samples]))
  mean_control <- apply(as.matrix(temp), 1, function(x) mean(x[control_samples]))
  log2_mean_ratio <- apply(as.matrix(temp), 1, function(x) {
    log_ratio <- log2(mean(x[stressed_samples]) / mean(x[control_samples]))
    if (is.nan(log_ratio)) {
      log_ratio <- 0
    }
    if (is.infinite(log_ratio)) {
      log_ratio <- 0
    }
    log_ratio
  })
  stats <- data.frame(wilcox_p_value = wp_values, mean_stressed = mean_stressed, mean_control = mean_control, log2_mean_ratio = log2_mean_ratio)
  stats$OTU <- names(wp_values)
  df_tax_otu$stats <- stats
  df_tax_otu$data <- left_join(df_tax_otu$data, stats, by = "OTU")
  return(df_tax_otu)
}


#' @title Convert Phyloseq Object to Statistical Excel File
#' @description This function takes a phyloseq object and converts it to a list of
#' dataframes followed by analyzing the treatment groups and adding statistical data to
#' the same list.  It then creates an excel file from this object.
#' @return An Excel file containing all of the phyloseq data with additional statistical data.
#' @param phyloseq_obj A phyloseq object.
#' @param treatment_groups A list of the treatment groups in the sample metadata.
#' @param rank The agglomerated rank of the phyloseq object
#' @param file_path Absolute path of the Excel file that will be created.
#' @return No objects are returned, but an Excel file is created.
#' @pretty_print TRUE
#' @export
#' @family Data Manipulators
#' @rdname phyloseq_to_excel
#' @seealso
#'  \code{\link[dplyr]{select_all}}
#'  \code{\link[openxlsx]{createWorkbook}},\code{\link[openxlsx]{addWorksheet}},\code{\link[openxlsx]{writeDataTable}},\code{\link[openxlsx]{saveWorkbook}}
#' @importFrom dplyr select_if
#' @importFrom openxlsx createWorkbook addWorksheet writeDataTable saveWorkbook
phyloseq_to_excel <- function(phyloseq_obj, treatment_groups, rank, file_path) {
  df_tax_otu <- phyloseq_to_stats_dataframe(phyloseq_obj = phyloseq_obj, treatment_groups = treatment_groups)
  stressed_samples <- df_tax_otu$stressed_samples
  control_samples <- df_tax_otu$control_samples
  # TODO:  Remove Control/Stress from syntax.  Make dynamic for other users.
  main_df <- df_tax_otu$data[49:60][, c("OTU", rank, "mean_stressed", "mean_control", "wilcox_p_value", "log2_mean_ratio", unlist(ranks[ranks != rank]))]
  main_df <- main_df %>% dplyr::select_if(~sum(!is.na(.)) > 0)
  wb <- openxlsx::createWorkbook(file_path)
  openxlsx::addWorksheet(wb, "Stats_Taxonomy")
  openxlsx::addWorksheet(wb, "Stressed")
  openxlsx::addWorksheet(wb, "Control")
  openxlsx::addWorksheet(wb, "Master")
  openxlsx::writeDataTable(wb, 1, main_df)
  openxlsx::writeDataTable(wb, 2, df_tax_otu$data[, c("OTU", stressed_samples)])
  openxlsx::writeDataTable(wb, 3, df_tax_otu$data[, c("OTU", control_samples)])
  openxlsx::writeDataTable(wb, 4, df_tax_otu$data)
  openxlsx::saveWorkbook(wb, file = file_path, overwrite = TRUE)
}


#' @title Convert Phyloseq Object to Dataframe
#' @description This function takes a phyloseq object and converts it to a list of
#' dataframes for further processing.
#' @param phyloseq_obj A phyloseq object.
#' @param treatment_groups A list of the treatment groups in the sample metadata.
#' @return A list of dataframes containing all of the phyloseq data
#' @pretty_print TRUE
#' @export
#' @family Data Manipulators
#' @rdname phyloseq_to_dataframe
#' @seealso
#'  \code{\link[phyloseq:tax_table-methods]{phyloseq::tax_table()}},\code{\link[phyloseq:otu_table-methods]{phyloseq::otu_table()}},\code{\link[phyloseq:sample_data-methods]{phyloseq::sample_data}}
#'  \code{\link[dplyr]{mutate}},\code{\link[dplyr]{join}}
#'  \code{\link[stringr]{case}}
#' @importFrom phyloseq tax_table otu_table sample_data
#' @importFrom dplyr mutate left_join
#' @importFrom stringr str_to_lower
phyloseq_to_dataframe <- function(phyloseq_obj, treatment_groups) {
  df_tax_otu <- list()
  p_tax <- phyloseq::tax_table(phyloseq_obj)
  p_otu <- phyloseq::otu_table(phyloseq_obj)
  p_sam <- phyloseq::sample_data(phyloseq_obj)
  df_tax_otu$tax <- data.frame(p_tax)
  df_tax_otu$otu_names <- rownames(df_tax_otu$tax)
  df_tax_otu$tax <- dplyr::mutate(df_tax_otu$tax, OTU = rownames(df_tax_otu$tax))
  df_tax_otu$otu <- data.frame(p_otu)
  df_tax_otu$otu <- dplyr::mutate(df_tax_otu$otu, OTU = rownames(df_tax_otu$otu))
  df_tax_otu$sam <- data.frame(p_sam)
  df_tax_otu$treat_groups <- treatment_groups
  df_tax_otu$tg_index <- c()
  for (tgroup in treatment_groups) {
    tgindex <- stringr::str_to_lower(tgroup)
    tgindex <- stringr::str_replace(tgindex, " ", "_")
    tgindex <- sprintf("%s_samples", tgindex)
    df_tax_otu[[tgindex]] <- rownames(df_tax_otu$sam[df_tax_otu$sam$TreatmentGroup == tgroup, ])
    df_tax_otu$tg_index <- c(df_tax_otu$tg_index, tgindex)
  }
  df_tax_otu$data <- dplyr::left_join(x = df_tax_otu$otu, y = df_tax_otu$tax, by = "OTU")
  return(df_tax_otu)
}

#' @title Remove Ambiguous Taxonomy
#' @description This function removes ambiguous taxonomy names from specified ranks
#' in a phyloseq object.
#' @param phyloseq_obj A phyloseq object to filter.
#' @param ranks A list of ranks to iterate over.
#' @param ambiguous_names a list of ambiguous names to remove.
#' @return Returns a phyloseq object without the ambiguous names in the specified ranks.
#' MicrobiomeR was also written using the tidyverse in order to make the workflow reproducible.
#' @pretty_print TRUE
#' @details This
#' @examples
#' \dontrun{
#'
#' > phy_obj <- get_phyloseq_object(...)
#' > phy_obj <- remove_ambiguous_taxa(phy_obj, list("Phylum", "Class", "Order"),
#' list("Unkown" , "Ambiguous", "uncultured"))
#'
#' }
#' @export
#' @family Filters
#' @rdname remove_ambiguous_taxa
#' @seealso
#'  \code{\link[phyloseq:tax_table-methods]{phyloseq::tax_table}}, \code{\link[phyloseq:phy_tree-methods]{phyloseq::phy_tree}}, \code{\link[phyloseq:otu_table-methods]{phyloseq::otu_table}}, \code{\link[phyloseq:sample_data-methods]{phyloseq::sample_data()}}, \code{\link[phyloseq]{merge_phyloseq}}
#'
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{select}}, \code{\link[tibble]{rownames}}, \code{\link[stringr]{str_detect}}
#' @importFrom dplyr filter select
#' @importFrom tibble column_to_rownames
#' @importFrom phyloseq tax_table phy_tree otu_table sample_data merge_phyloseq phy_tree
#' @importFrom stringr str_detect
remove_ambiguous_taxa <- function(phyloseq_obj, ranks, ambiguous_names) {
  df_tax_otu <- phyloseq_to_tibble(phyloseq_obj)
  phyloseq_data <- df_tax_otu$data
  # Loop through ranks and remove ambiguous names from various ranks.
  for (r in ranks) {
    if (r == "Kingdom") {
      phyloseq_data <- dplyr::filter(phyloseq_data, !Kingdom %in% ambiguous_names)
    } else if (r == "Phylum") {
      phyloseq_data <- dplyr::filter(phyloseq_data, !Phylum %in% ambiguous_names)
    } else if (r == "Class") {
      phyloseq_data <- dplyr::filter(phyloseq_data, !Class %in% ambiguous_names)
    } else if (r == "Order") {
      phyloseq_data <- dplyr::filter(phyloseq_data, !Order %in% ambiguous_names)
    } else if (r == "Famlily") {
      phyloseq_data <- dplyr::filter(phyloseq_data, !Family %in% ambiguous_names)
    } else if (r == "Genus") {
      phyloseq_data <- dplyr::filter(phyloseq_data, !Genus %in% ambiguous_names)
    } else if (r == "Species") {
      phyloseq_data <- dplyr::filter(phyloseq_data, !Speices %in% ambiguous_names)
    }
  }

  tax <- dplyr::select(phyloseq_data, colnames(df_tax_otu$tax)) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("OTU") %>%
    as.matrix()

  otu <- dplyr::select(phyloseq_data, colnames(df_tax_otu$otu)) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("OTU") %>%
    as.matrix()
  # Merge into phyloseq data;  Include tree if applicable
  tax <- phyloseq::tax_table(tax)
  if (inherits(phyloseq_obj, "taxonomyTable")) {
    return(tax)
  } else {
    err_ <- try(phy_tree(phyloseq_obj), silent = TRUE)
    otu <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
    sam <- phyloseq::sample_data(phyloseq_obj)
    if (stringr::str_detect(err_, "Error")) {
      phyloseq_obj <- phyloseq::merge_phyloseq(tax, otu, sam)
    } else {
      ptree <- phyloseq::phy_tree(phyloseq_obj)
      phyloseq_obj <- phyloseq::merge_phyloseq(tax, ptree, otu, sam)
    }
    return(phyloseq_obj)
  }
}



#' @title Preprocessing for Phyloseq Objects
#' @description This function uses techniques from the Bioconductor Workflow in order to
#' preprocess phyloseq data for downstream analysis.
#' @param phyloseq_object A phyloseq object
#' @param process_list This parameter is used to control the way the phyloseq object is
#' processed.  It can be one of three values:
#' \describe{
#'   \item{NULL}{Will utilize the default processing strategy or the dot parameters (...).}
#'   \item{file}{The absolute path to a YAML config file.}
#'   \item{list}{Equivalent to putting dot parameters in a list (list(...)).}
#' }
#'
#' @param ... The dot parameters can be any combination of the default processing keywords.
#' See the "Keywords for Processing" section below for more details
#' @return Returns a phyloseq object that has undergone the specified processing strategies.
#' @pretty_print TRUE
#' @details DETAILS
#' @examples
#' \dontrun{
#'
#' > phy_obj <- get_phyloseq_object(...)
#' > phy_obj
#' phyloseq-class experiment-level object
#' otu_table()   OTU Table:         [ 1955 taxa and 48 samples ]
#' sample_data() Sample Data:       [ 48 samples by 8 sample variables ]
#' tax_table()   Taxonomy Table:    [ 1955 taxa by 7 taxonomic ranks ]
#' phy_tree()    Phylogenetic Tree: [ 1955 tips and 1954 internal nodes ]
#' > pp_phy_obj <- preprocess_phyloseq(phy_obj, master_thresh = 1e-5,
#'                                     taxon_filter = list("Phylum"= list("min_a"=5, "r_s_p"=0.5),
#'                                                         "Class"=list("min_a"=3, "r_s_p"=0.3)),
#'                                     prevalence_filter = list("min_a"=5, "r_s_p"=0.5), glom_rank = NULL,
#'                                     ambiguous=list(amb_ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
#'                                     amb_items = c(NA, "", "uncharacterized", "uncultured", "Unassigned", "Ambiguous", "Ambiguous_taxa")),
#'                                     coeff_of_variation = 0.55, trans_function = function(x){x / sum(x)}, merge_samp = NULL)
#' > pp_phy_obj
#'
#' }
#' @export
#' @family Filters
#' @rdname preprocess_phyloseq
#' @seealso
#'  \code{\link[phyloseq:prune_taxa-methods]{phyloseq::prune_taxa()}}, \code{\link[phyloseq]{taxa_sums}}, \code{\link[phyloseq]{filter_taxa}}, \code{\link[phyloseq]{tax_glom}}, \code{\link[phyloseq:nsamples-methods]{phyloseq::nsamples()}},
#'  \code{\link[phyloseq]{filterfun_sample}}, \code{\link[phyloseq:genefilter_sample-methods]{phyloseq::genefilter_sample()}}, \code{\link[phyloseq]{get_taxa_unique}}, \code{\link[phyloseq:transformcounts]{phyloseq::transform_sample_counts()}},
#'  \code{\link[phyloseq:merge_samples-methods]{phyloseq::merge_samples()}}, \code{\link[phyloseq:subset_taxa-methods]{phyloseq::subset_taxa()}}
#'
#'  \code{\link[yaml]{yaml.load}}\cr\cr
#'  Bioconductor Workflow - \url{https://f1000research.com/articles/5-1492/v2}\cr
#'  Phyloseq Website - \url{https://joey711.github.io/phyloseq/index.html}
#' @section Keywords for Pre-processing:
#' \describe{
#'   \item{master_thresh}{DEFAULT: 1e-5.  Filters any taxa that do not meat the mean threshold.}
#'   \item{taxon_filter}{DEFAULT: list("Phylum"= list("min_a"=5, "r_s_p"=0.5),
#' "Class"=list("min_a"=3, "r_s_p"=0.3)).  Filters OTUs that do not appear more than a certian
#' amount of times in a certain percentage of samples at the specified agglomerated rank.}
#'   \item{prevelance_filter}{DEFAULT: list("min_a"=5, "r_s_p"=0.5).  Filters OTUs that do not
#'   appear more than a certian amount of times in a certain percentage of samples.}
#'   \item{glom_rank}{DEFAULT: NULL.  Agglomerates the data at the specified rank.}
#'   \item{ambiguous}{DEFAULT: list(amb_ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
#' amb_items = c(NA, "", "uncharacterized", "uncultured", "Unassigned", "Ambiguous", "Ambiguous_taxa")).
#' Removes OTUs that are labeled with the specified ambiguous items.  This is done for each specified rank.}
#'   \item{coeff_of_variation}{DEFAULT: 0.55.  Standardizes abundances to the median sequencing depth}
#'   \item{trans_function}{DEFAULT:  function(x){x / sum(x)}.  Transforms the abundance values to relative
#'   abundance values.}
#'   \item{merge_samples}{DEFAULT:  NULL.}
#' }
#' @importFrom yaml yaml.load_file
#' @importFrom phyloseq taxa_sums prune_taxa filter_taxa tax_glom nsamples filterfun_sample genefilter_sample get_taxa_unique transform_sample_counts merge_samples
preprocess_phyloseq <- function(phyloseq_object, process_list = NULL, ...) {
  dotparam <- list(...)
  # Set up processing parameters
  if (length(dotparam) != 0) {
    process_list <- dotparam
    } else if (is.null(process_list)) {
      warning("Warning:  Using default parameters in extdata/process_list.yml")
      process_list <- yaml::yaml.load_file("extdata/process_list.yml")
      } else if (is.character(process_list) & (file.exists(process_list))) {
        process_list <- yaml::yaml.load_file(process_list)
        }

  ## REMOVAL
  # Remove empty samples
  t_sums <- phyloseq::taxa_sums(phyloseq_object)
  processed_phy_obj <- phyloseq::prune_taxa(t_sums > 0, phyloseq_object)

  ## FILTER
  for (proc in names(process_list)) {
    if (proc == "master_thresh" & !is.null(process_list[[proc]])) {
      # Filter anything at a certain threshhold
      # phyloseq vignette - master_thresh = 1e-5
      processed_phy_obj <- phyloseq::filter_taxa(processed_phy_obj, function(x) {
        mean(x) > process_list[[proc]]
      }, TRUE)
    }
    if (proc == "taxon_filter" & !is.null(process_list[[proc]])) {
      taxon_filter <- process_list[[proc]]
      # Taxonomic Prevalence filtering
      # Filter OTUs that do not appear more than a certian amount of times in a certain percentage of samples
      # at the specified agglomerated rank
      # phyloseq vignette - taxa_thresh = c(5, 0.5)
      counter <- 0
      for (tf_rank in names(taxon_filter)) {
          if (length(taxon_filter[[tf_rank]]) == 2) {
              # Set up data for filtering
              tf_min_abund <- taxon_filter[[tf_rank]][["min_a"]]
              tf_req_samp_perc <- taxon_filter[[tf_rank]][["r_s_p"]]
              tf_glom <- phyloseq::tax_glom(processed_phy_obj, tf_rank, NArm = FALSE)
              n_samp <- phyloseq::nsamples(tf_glom)
              filter_fun <- phyloseq::filterfun_sample(function(x) x > tf_min_abund)

              # Filter
              tf_unfiltered <- phyloseq::genefilter_sample(tf_glom, filter_fun,
                A = tf_req_samp_perc * n_samp)
              p_filter <- phyloseq::prune_taxa(tf_unfiltered, tf_glom)

              # Formating for standard evaluation
              r_filter <- phyloseq::get_taxa_unique(p_filter, taxonomic.rank = tf_rank)
              if (tf_rank == "Kingdom") {
                processed_phy_obj <- phyloseq::subset_taxa(processed_phy_obj, Kingdom %in% r_filter)
              } else if (tf_rank == "Phylum") {
                processed_phy_obj <- phyloseq::subset_taxa(processed_phy_obj, Phylum %in% r_filter)
              } else if (tf_rank == "Class") {
                processed_phy_obj <- phyloseq::subset_taxa(processed_phy_obj, Class %in% r_filter)
              } else if (tf_rank == "Order") {
                processed_phy_obj <- phyloseq::subset_taxa(processed_phy_obj, Order %in% r_filter)
              } else if (tf_rank == "Famlily") {
                processed_phy_obj <- phyloseq::subset_taxa(processed_phy_obj, Family %in% r_filter)
              } else if (tf_rank == "Genus") {
                processed_phy_obj <- phyloseq::subset_taxa(processed_phy_obj, Genus %in% r_filter)
              } else if (tf_rank == "Species") {
                processed_phy_obj <- phyloseq::subset_taxa(processed_phy_obj, Speices %in% r_filter)
              }
          }
        }
      }
    if (proc == "prevelance_filter" & length(process_list[[proc]]) == 2) {
      prevalence_filter <- process_list[[proc]]
      # Prevelance filtering
      # Filter OTUs that do not appear more than a certian amount of times in a certain percentage of samples.
      # phyloseq vignette - taxa_thresh = c(5, 0.5)

      # Set up data for filtering
      min_abund <- prevalence_filter$min_a
      required_sample_percentage <- prevalence_filter$r_s_p
      n_samp <- phyloseq::nsamples(processed_phy_obj)
      filter_fun <- phyloseq::filterfun_sample(function(x) x > min_abund)

      # Filter
      tx <- genefilter_sample(processed_phy_obj, filter_fun,
        A = required_sample_percentage * n_samp)
      processed_phy_obj <- phyloseq::prune_taxa(tx, processed_phy_obj)
      }
    if (proc == "glom_rank" & !is.null(process_list[[proc]])) {
      processed_phy_obj <- phyloseq::tax_glom(processed_phy_obj, taxrank = process_list[[proc]])
    }
    if (proc == "ambiguous" & !is.null(process_list[[proc]])) {
      ambiguous <- process_list[[proc]]
      amb_ranks <- ambiguous$amb_ranks
      amb_items <- ambiguous$amb_items
      processed_phy_obj <- remove_ambiguous_taxa(processed_phy_obj, amb_ranks, amb_items)
    }
    ## Standardize
    if (proc == "coeff_of_variation" & !is.null(process_list[[proc]])) {
      # Standardize abundances to the median sequencing depth
      total <- median(sample_sums(processed_phy_obj))
      standf <- function(x, t = total) round(t * (x / sum(x)))
      pp_trans <- phyloseq::transform_sample_counts(processed_phy_obj, standf)
      # Filter the taxa using the coefficient of variation
      # phyloseq vignette - coeff_of_variation = 0.3
      pp_filter <- phyloseq::filter_taxa(pp_trans, function(x) sd(x) / mean(x) > process_list[[proc]])
      processed_phy_obj <- phyloseq::prune_taxa(pp_filter, processed_phy_obj)
    }

    ## TRANSFORM
    if (proc == "trans_function" & class(process_list[[proc]]) == "function") {
      # Transform to get relative abundance
      # Note:  always do this before removing taxa that contribute to rel abund
      trans_function <- process_list[[proc]]
      processed_phy_obj <- phyloseq::transform_sample_counts(processed_phy_obj, trans_function)
    }
    if (proc == "merge_samp" & !is.null(process_list[[proc]])) {
      sample_group <- process_list[[proc]][["group"]]
      func <- process_list[[proc]][["func"]]
      processed_phy_obj <- phyloseq::merge_samples(processed_phy_obj, group = sample_group, fun = func)
    }
  }
  ## REMOVAL
  # Remove empty samples
  processed_phy_obj <- phyloseq::prune_taxa(taxa_sums(processed_phy_obj) > 0, processed_phy_obj)
  return(processed_phy_obj)
}

#'  @title Get Distance Methods
#'
#'  @description The function for manipulating phyloseqs native distanceMethodList
#'
#'  @param tree_methods DESCRIPTION.
#'
#'  @return RETURN_DESCRIPTION
#'  @export
get_distance_methods <- function(tree_methods = TRUE) {
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
