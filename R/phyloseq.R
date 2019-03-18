#' @title Get Phyloseq Object
#' @description Create a phyloseq object using a biom file, phylogenetic tree file, and a
#'  metadata file. Alternatively, an .Rdata file can be used.  For now this function takes
#'  output from Qiime.  The data provided by this package was processed by the NIH's Nephele
#'  pipeline.
#' @param biom_file A file in the \href{http://biom-format.org/}{BIOM format}.  Default: NULL
#' @param tree_file A phylogenetic tree file.  This function will generate a rooted
#' tree file if your phylogenetic tree isn't already rooted.  Default: NULL
#' @param metadata_file Sample metadata in tab delimited format.  It's very important that you make
#'  sure that you have a \strong{sample_id} and \strong{TreatmentGroup} column in your metadata
#'  file.  You might otherwise run into problems.  Default: NULL
#' @param treatment_group The column number or name in the metadata file that contains the
#' treatment group names.  Default: NULL
#' @param parse_func The parse function used to parse taxonomy strings from
#'  Greengenes or SILVA database.  Default: NULL
#' @param rdata_file A .Rdata file.  Default: NULL
#' @param save_rooted_tree A logical that determines weather or not the rooted tree
#' is saved.  Default: FALSE
#' @param recursive_save If the directory doesn't exists create the parent directories
#' that don't exist as well.  Default: FALSE
#' @return A phyloseq object, and a phylogenetic tree file if one does not already exist.
#' @pretty_print TRUE
#' @details This function heavily relies on the phyloseq package to import data into R.
#'  It also requires you to use an absolute or relative path to your data files or for your data files
#'  to be in the working directory.
#'
#'  It's very important that you make sure that you have a \strong{sample_id} and \strong{TreatmentGroup}
#'  column in your metadata file.
#' @examples
#' \dontrun{
#' if(interactive()){
#' library(MicrobiomeR)
#' # Get the datafiles from the package
#' biom_file <- system.file("extdata", "silva_OTU.biom", package = "MicrobiomeR")
#' tree_file <- system.file("extdata", "silva.tre", package = "MicrobiomeR")
#' md_file <- system.file("extdata", "nephele_metadata2.txt", package = "MicrobiomeR")
#' parse_func <- parse_parse_taxonomy_silva_128
#' }
#' # Create a phyloseq object from the data files.
#' phy_obj <- create_phyloseq(biom_file = biom_file,
#'                            tree_file = tree_file,
#'                            metadata_file = md_file,
#'                            parse_func = parse_func)
#' }
#' @export
#' @family Data Importers
#' @rdname create_phyloseq
#' @seealso
#'  \code{\link[phyloseq]{import_biom}}, \code{\link[phyloseq:sample_data-methods]{phyloseq::sample_data()}}, \code{\link[phyloseq:phy_tree-methods]{phyloseq::phy_tree()}}, \code{\link[phyloseq]{import_qiime_sample_data}}, \code{\link[phyloseq]{merge_phyloseq}}
#'
#'  \code{\link[ape]{root}}, \code{\link[ape:root]{is.rooted}}, \code{\link[ape]{read.tree}}, \code{\link[ape]{write.tree}}
#'
#'  \code{\link[MicrobiomeR]{root_by_longest_edge}}
#' @importFrom phyloseq import_biom sample_data phy_tree import_qiime_sample_data merge_phyloseq
#' @importFrom ape is.rooted write.tree read.tree
#' @importFrom tools file_path_as_absolute
create_phyloseq <- function(biom_file = NULL, tree_file = NULL, metadata_file = NULL,
                             treatment_group = NULL, parse_func = NULL, rdata_file = NULL,
                             save_rooted_tree = FALSE, recursive_save = FALSE) {
  if (!is.null(rdata_file)) {
    load(rdata_file)
    phyloseq_object <- get("phyloseq_object")
    return(phyloseq_object)
  } else {
    # Load the biom file data amd qiime metadata separately, and then merge them into a phyloseq object
    phyloseq_biom <- phyloseq::import_biom(BIOMfilename = biom_file, parseFunction = parse_func, treefilename = tree_file)
    phyloseq_metadata <- phyloseq::import_qiime_sample_data(metadata_file)
    phyloseq_object <- phyloseq::merge_phyloseq(phyloseq_biom, phyloseq_metadata)

    # Create the "X.TreatmentGroup" variable for the metadata
    # This step gives other functions a standard treatment group variable to work with
    if (!is.null(treatment_group)){
      if (is.numeric(treatment_group) || is.character(treatment_group)){
        phyloseq::sample_data(phyloseq_object)[["TreatmentGroup"]] <- phyloseq::sample_data(phyloseq_object)[[treatment_group]]
      } else {
        warning("The treatment_group parameter must be numeric or a character string.")
        warning("The data might not be appropriate for other MicrobiomeR functions.  Please try again.")
        stop()
      }
    }
    phyloseq_object <- root_phyloseq_tree(phyloseq_object = phyloseq_object, tree_path = tree_file,
                                          save_rooted_tree = save_rooted_tree,
                                          recursive = recursive_save)
    return(phyloseq_object)
  }
}

#' @title Root phyloseq tree.
#' @description A function for rooting and saving a phylogenetic tree.
#' @param phyloseq_object A phyloseq object that contains a phylogenetic tree.
#' @param tree_path The path to the existing or desired phylogenetic tree file.
#' @param save_rooted_tree A logical that determines weather or not the rooted tree
#' is saved.
#' @param recursive If the directory doesn't exist, create the parent directories
#' that don't exist as well, Default: TRUE
#' @return Returns a phyloseq object with a rooted tree.
#' @pretty_print TRUE
#' @details This function is a helper function to get a proper phyloseq object for
#' downstream analysis.  Some analyses require a rooted tree.  The function saves
#' the rooted tree in the phyloseq object.  It can also save the tree as a file
#' if desired.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @family Phylogenetic Tree Manipulators
#' @rdname root_phyloseq_tree
#' @seealso
#'  \code{\link[ape]{root}}, \code{\link[ape]{read.tree}}, \code{\link[ape]{write.tree}}
#'
#'  \code{\link[phyloseq:phy_tree-methods]{phy_tree}}
#'
#'  \code{\link[tools]{fileutils}}
#' @importFrom ape is.rooted read.tree write.tree
#' @importFrom phyloseq phy_tree
#' @importFrom tools file_path_as_absolute file_ext
root_phyloseq_tree <- function(phyloseq_object, tree_path, save_rooted_tree, recursive = TRUE) {
  # Write to a tree file that subsets our data and update the phyloseq object
  if (!is.null(tree_path) && !ape::is.rooted(ape::read.tree(tree_path))) {
    # Root the tree
    phyloseq_tree <- phyloseq::phy_tree(phyloseq_object)
    ape_tree <- root_by_longest_edge(phyloseq_tree)
    # Write the tree file
    file_isdir <- file.info(tree_path)[1, "isdir"]
    tp <- dirname(tools::file_path_as_absolute(tree_path))
    if (file_isdir == TRUE) { # Is existing directory
      tf <- "rooted_tree.tre"
    } else if (file_isdir == FALSE) { # Is existing file
      tf <- basename(tree_path)
    } else if (is.na(file_isdir)) { # Does not exist
      if (tools::file_ext(tree_path) == "") { # Is it a directory string?
        if (!dir.create(tree_path, recursive = recursive) == TRUE) {
          stop("Unable to create directory.  Please check the path to your tree file.")
        }
      }
    } else { # It's a file string.
        tf <- basename(tree_path)
        if (!dir.exists(tp)) {
          if (!dir.create(tree_path, recursive = recursive) == TRUE) {
            stop("Unable to create directory.  Please check the path to your tree file.")
          }
        }
      }
    new_tree_path <- sprintf("%s/rooted_%s", tp, tf)
    if (!file.exists(new_tree_path)) {
      if (save_rooted_tree == TRUE) {
        ape::write.tree(ape_tree, new_tree_path)
        }
    }
    # Update and return the phyloseq object
    phyloseq::phy_tree(phyloseq_object) <- ape_tree
    return(phyloseq_object)
  } else {
    return(phyloseq_object)
  }
}

#' @title Pick Outgroup for Tree
#' @description Pick an outgroup for rooting a phylogenetic tree.
#' @param unrooted_tree An unrooted tree object.
#' @return A new tree with the longest branch as the outgroup.
#' @pretty_print TRUE
#' @details This function preprocess a phylogenetic tree for rooting by the longest edge.
#' Please see this issue for more details \url{https://github.com/joey711/phyloseq/issues/597}.
#' @export
#' @family Phylogenetic Tree Manipulators
#' @rdname pick_new_outgroup
#' @seealso
#'  \code{\link[data.table]{data.table}}
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
  treeDT <- cbind(tree_DT, tree_DT_len)[1:tree_tips] %>% cbind(tree_DT_id)
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
#'  \code{\link[MicrobiomeR]{pick_new_outgroup}}
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
#' specifically the \code{\link{import_biom}} function, but
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
#'  > taxvec1 = c("Root", "k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales",
#'  "f__Staphylococcaceae")
#'  > parse_taxonomy_default(taxvec1)
#'  > parse_taxonomy_greengenes(taxvec1)
#'  > taxvec2 = c("Root;k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Staphylococcaceae")
#'  > parse_taxonomy_qiime(taxvec2)
#'  > taxvec3 = c("D_0__Bacteria", "D_1__Firmicutes", "D_2__Bacilli", "D_3__Staphylococcaceae")
#'  > parse_taxonomy_silva_128(taxvec3)
#'  }
#' @export
#' @rdname parse_taxonomy_silva_128
#' @seealso
#'  \code{\link[phyloseq:parseTaxonomy-functions]{parse_taxonomy_default}}
#'
#'  \code{\link[phyloseq:parseTaxonomy-functions]{parse_taxonomy_greengenes}}
#'
#'  \code{\link[phyloseq:parseTaxonomy-functions]{parse_taxonomy_qiime}}
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
