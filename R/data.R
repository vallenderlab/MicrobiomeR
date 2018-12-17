#' Stress Microbiome Phylseq Data (SILVA)
#'
#' Data from the Nephele-16S-Qiime-Silva microbiome pipeline
#' in the form of a phyloseq object.
#'
#' @docType data
#' @usage data(phyloseq_silva)
#' @format An object of class \code{"phyloseq"}
#' @keywords datasets
#'
#' @references Zhang et al. (2019)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/}{PubMed})
#'
#' @source \href{https://nephele.niaid.nih.gov/user_guide_pipes/#amplicon_pipes}{Nephele-16S Qiime Pipeline}
"phyloseq_silva"

#' Stress Microbiome Phylseq Data (Greengenes)
#'
#' Data from the Nephele-16S-Qiime-Greengenes microbiome pipeline
#' in the form of a phyloseq object.
#'
#' @docType data
#' @usage data(phyloseq_greengenes)
#' @format An object of class \code{"phyloseq"}
#' @keywords datasets
#'
#' @references Zhang et al. (2019)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/}{PubMed})
#'
#' @source \href{https://nephele.niaid.nih.gov/user_guide_pipes/#amplicon_pipes}{Nephele-16S Qiime Pipeline}
"phyloseq_greengenes"

#' Stress Microbiome "raw_format" Data (Silva)
#'
#' Data from the Nephele-16S-Qiime-Silva microbiome pipeline
#' in the form of a taxmap object.  This taxmap object has been
#' formatted into the MicrobiomeR "raw_format".
#'
#' @docType data
#' @usage data(raw_silva)
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
"raw_silva"

#' Stress Microbiome "raw_format" Data (Greengenes)
#'
#' Data from the Nephele-16S-Qiime-Greengenes microbiome pipeline
#' in the form of a taxmap object.  This taxmap object has been
#' formatted into the MicrobiomeR "raw_format".
#'
#' @docType data
#' @usage data(raw_greengenes)
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
"raw_greengenes"

#' Stress Microbiome "basic_format" Data (Silva)
#'
#' Data from the Nephele-16S-Qiime-Silva microbiome pipeline
#' in the form of a taxmap object.  This taxmap object has been
#' formatted into the MicrobiomeR "basic_format".
#'
#' @docType data
#' @usage data(basic_silva)
#' @format An object of class \code{"Taxmap"} with customized observation tables:
#' \describe{
#' \item{otu_abundance}{An abundance table with otu_id, taxon_id, and Sample_* columns.}
#' \item{otu_annotations}{A taxonomy annotations table with otu_id, taxon_id, and multiple rank columns.}
#' \item{otu_proportions}{A proportion table derived from the otu_abundance table.}
#' \item{sample_data}{A metadata table with SampleID, BarcodeSequence, LinkerPrimerSequence,
#' ForwardFastqFile, ReverseFastqFile, TreatmentGroup, SampleName, Description}
#' \item{phy_tree}{A phylogenetic tree.}
#' \item{taxa_abundance}{A taxonomy table interpolated from the otu_abundance table with taxon_id,
#' and Sample_* columns.}
#' \item{taxa_proportions}{A taxonomy table derived from the taxa_abundance table.}
#' }
#' @keywords datasets
#'
#' @references Zhang et al. (2019)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/}{PubMed})
#'
#' @source \href{https://nephele.niaid.nih.gov/user_guide_pipes/#amplicon_pipes}{Nephele-16S Qiime Pipeline}
"basic_silva"

#' Stress Microbiome "basic_format" Data (Greengenes)
#'
#' Data from the Nephele-16S-Qiime-Greengenes microbiome pipeline
#' in the form of a taxmap object.  This taxmap object has been
#' formatted into the MicrobiomeR "basic_format".
#'
#' @docType data
#' @usage data(basic_greengenes)
#' @format An object of class \code{"Taxmap"} with customized observation tables:
#' \describe{
#' \item{otu_abundance}{An abundance table with otu_id, taxon_id, and Sample_* columns.}
#' \item{otu_annotations}{A taxonomy annotations table with otu_id, taxon_id, and multiple rank columns.}
#' \item{otu_proportions}{A proportion table derived from the otu_abundance table.}
#' \item{sample_data}{A metadata table with SampleID, BarcodeSequence, LinkerPrimerSequence,
#' ForwardFastqFile, ReverseFastqFile, TreatmentGroup, SampleName, Description}
#' \item{phy_tree}{A phylogenetic tree.}
#' \item{taxa_abundance}{A taxonomy table interpolated from the otu_abundance table with taxon_id,
#' and Sample_* columns.}
#' \item{taxa_proportions}{A taxonomy table derived from the taxa_abundance table.}
#' }
#' @keywords datasets
#'
#' @references Zhang et al. (2019)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/}{PubMed})
#'
#' @source \href{https://nephele.niaid.nih.gov/user_guide_pipes/#amplicon_pipes}{Nephele-16S Qiime Pipeline}
"basic_greengenes"

#' Stress Microbiome "analyzed_format" Data (Silva)
#'
#' Data from the Nephele-16S-Qiime-Silva microbiome pipeline
#' in the form of a taxmap object.  This taxmap object has been
#' formatted into the MicrobiomeR "analyzed_format".
#'
#' @docType data
#' @usage data(analyzed_silva)
#' @format An object of class \code{"Taxmap"} with customized observation tables:
#' \describe{
#' \item{otu_abundance}{An abundance table with otu_id, taxon_id, and Sample_* columns.}
#' \item{otu_annotations}{A taxonomy annotations table with otu_id, taxon_id, and multiple rank columns.}
#' \item{otu_proportions}{A proportion table derived from the otu_abundance table.}
#' \item{sample_data}{A metadata table with SampleID, BarcodeSequence, LinkerPrimerSequence,
#' ForwardFastqFile, ReverseFastqFile, TreatmentGroup, SampleName, Description}
#' \item{phy_tree}{A phylogenetic tree.}
#' \item{taxa_abundance}{A taxonomy table interpolated from the otu_abundance table with taxon_id,
#' and Sample_* columns.}
#' \item{taxa_proportions}{A taxonomy table derived from the taxa_abundance table.}
#' }
#' \item{statistical_data}{A data table that contains statistical analysis comparing treatment groups.}
#' \item{stats_tax_data}{A data table that combines taxa and statistical datasets for plotting.}
#' @keywords datasets
#'
#' @references Zhang et al. (2019)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/}{PubMed})
#'
#' @source \href{https://nephele.niaid.nih.gov/user_guide_pipes/#amplicon_pipes}{Nephele-16S Qiime Pipeline}
"analyzed_silva"

#' Stress Microbiome "analyzed_format" Data (Greengenes)
#'
#' Data from the Nephele-16S-Qiime-Greengenes microbiome pipeline
#' in the form of a taxmap object.  This taxmap object has been
#' formatted into the MicrobiomeR "analyzed_format".
#'
#' @docType data
#' @usage data(analyzed_greengenes)
#' @format An object of class \code{"Taxmap"} with customized observation tables:
#' \describe{
#' \item{otu_abundance}{An abundance table with otu_id, taxon_id, and Sample_* columns.}
#' \item{otu_annotations}{A taxonomy annotations table with otu_id, taxon_id, and multiple rank columns.}
#' \item{otu_proportions}{A proportion table derived from the otu_abundance table.}
#' \item{sample_data}{A metadata table with SampleID, BarcodeSequence, LinkerPrimerSequence,
#' ForwardFastqFile, ReverseFastqFile, TreatmentGroup, SampleName, Description}
#' \item{phy_tree}{A phylogenetic tree.}
#' \item{taxa_abundance}{A taxonomy table interpolated from the otu_abundance table with taxon_id,
#' and Sample_* columns.}
#' \item{taxa_proportions}{A taxonomy table derived from the taxa_abundance table.}
#' }
#' \item{statistical_data}{A data table that contains statistical analysis comparing treatment groups.}
#' \item{stats_tax_data}{A data table that combines taxa and statistical datasets for plotting.}
#' @keywords datasets
#'
#' @references Zhang et al. (2019)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/}{PubMed})
#'
#' @source \href{https://nephele.niaid.nih.gov/user_guide_pipes/#amplicon_pipes}{Nephele-16S Qiime Pipeline}
"analyzed_greengenes"
