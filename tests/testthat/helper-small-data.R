small_phyloseq_fixture <- local({
  cached_fixture <- NULL

  function(samples_per_group = 4L, taxa_limit = 100L) {
    if (is.null(cached_fixture)) {
      full_fixture <- MicrobiomeR::phyloseq_silva_2
      sample_metadata <- data.frame(phyloseq::sample_data(full_fixture))
      sample_ids_by_group <- split(
        rownames(sample_metadata),
        sample_metadata$TreatmentGroup
      )
      selected_sample_ids <- unlist(
        lapply(sample_ids_by_group, function(sample_ids) {
          head(sort(sample_ids), samples_per_group)
        }),
        use.names = FALSE
      )

      sample_fixture <- phyloseq::prune_samples(selected_sample_ids, full_fixture)
      ranked_taxa <- names(sort(phyloseq::taxa_sums(sample_fixture), decreasing = TRUE))

      # Retain abundant taxa so the fixture still exercises every analysis rank.
      cached_fixture <- phyloseq::prune_taxa(
        head(ranked_taxa, taxa_limit),
        sample_fixture
      )
    }

    cached_fixture
  }
})

small_taxmap_fixture <- local({
  cached_formats <- NULL

  function(format = c(
    "phyloseq_format",
    "raw_format",
    "basic_format",
    "analyzed_format"
  )) {
    format <- match.arg(format)

    if (is.null(cached_formats)) {
      phyloseq_format <- MicrobiomeR::create_taxmap(small_phyloseq_fixture())
      raw_format <- MicrobiomeR::as_raw_format(phyloseq_format$clone())
      basic_format <- MicrobiomeR::as_basic_format(raw_format$clone())
      analyzed_format <- MicrobiomeR::as_analyzed_format(basic_format$clone())

      cached_formats <- list(
        phyloseq_format = phyloseq_format,
        raw_format = raw_format,
        basic_format = basic_format,
        analyzed_format = analyzed_format
      )
    }

    # Taxmap uses reference semantics, so each test receives an isolated clone.
    cached_formats[[format]]$clone()
  }
})
