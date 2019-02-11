#' Phyloseq Data with 2 Treatment Groups
#'
#' Data from the Nephele-16S-Qiime-Silva microbiome pipeline
#' in the form of a phyloseq object.
#'
#' @docType data
#' @usage data(phyloseq_silva_2)
#' @format An object of class \code{"phyloseq"}
#' @keywords datasets
#'
#' @references Zhang et al. (2019)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/}{PubMed})
#'
#' @source \href{https://nephele.niaid.nih.gov/user_guide_pipes/#amplicon_pipes}{Nephele-16S Qiime Pipeline}
"phyloseq_silva_2"

#' Phyloseq Data with 3 Treatment Groups
#'
#' Data from the Nephele-16S-Qiime-Silva microbiome pipeline
#' in the form of a phyloseq object.
#'
#' @docType data
#' @usage data(phyloseq_silva_3)
#' @format An object of class \code{"phyloseq"}
#' @keywords datasets
#'
#' @references Zhang et al. (2019)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/}{PubMed})
#'
#' @source \href{https://nephele.niaid.nih.gov/user_guide_pipes/#amplicon_pipes}{Nephele-16S Qiime Pipeline}
"phyloseq_silva_3"

#' "raw_format" Data with 2 Treatment Groups
#'
#' Data from the Nephele-16S-Qiime-Silva microbiome pipeline
#' in the form of a taxmap object.  This taxmap object has been
#' formatted into the MicrobiomeR "raw_format".
#'
#' @docType data
#' @usage data(raw_silva_2)
#' @format An object of class \code{"Taxmap"} with customized observation tables:
#' \describe{
#' \item{otu_abundance}{An abundance table with otu_id, taxon_id, and Sample_* columns.}
#' \item{otu_annotations}{A taxonomy annotations table with otu_id, taxon_id, and multiple rank columns.}
#' \item{sample_data}{A metadata table with SampleID, BarcodeSequence, LinkerPrimerSequence,
#' ForwardFastqFile, ReverseFastqFile, TreatmentGroup, SampleName, Description}
#' \item{phy_tree}{A phylogenetic tree.}
#' }
#' @keywords datasets
#'
#' @references Zhang et al. (2019)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/}{PubMed})
#'
#' @source \href{https://nephele.niaid.nih.gov/user_guide_pipes/#amplicon_pipes}{Nephele-16S Qiime Pipeline}
"raw_silva_2"

#' "raw_format" Data with 3 Treatment Groups
#'
#' Data from the Nephele-16S-Qiime-Silva microbiome pipeline
#' in the form of a taxmap object.  This taxmap object has been
#' formatted into the MicrobiomeR "raw_format".
#'
#' @docType data
#' @usage data(raw_silva_3)
#' @format An object of class \code{"Taxmap"} with customized observation tables:
#' \describe{
#' \item{otu_abundance}{An abundance table with otu_id, taxon_id, and Sample_* columns.}
#' \item{otu_annotations}{A taxonomy annotations table with otu_id, taxon_id, and multiple rank columns.}
#' \item{sample_data}{A metadata table with SampleID, BarcodeSequence, LinkerPrimerSequence,
#' ForwardFastqFile, ReverseFastqFile, TreatmentGroup, SampleName, Description}
#' \item{phy_tree}{A phylogenetic tree.}
#' }
#' @keywords datasets
#'
#' @references Zhang et al. (2019)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/}{PubMed})
#'
#' @source \href{https://nephele.niaid.nih.gov/user_guide_pipes/#amplicon_pipes}{Nephele-16S Qiime Pipeline}
"raw_silva_3"
