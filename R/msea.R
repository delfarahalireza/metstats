#' Pathway Enrichment Analysis
#'
#' @description To perform Pathway Enrichment Analysis with any omics data type
#'
#' @param data_input a dataframe, must contain an id column and
#' a ranking metric column
#' if an optional 'feature_name' column exists in the input data,
#'' id_column' will be replaced by 'feature_name' in the final enrichment table
#'
#' @param id_column name of the column to be used as id
#'
#' @param feature_sets_file path to the feature sets file with a
#' standard GSEA tab-delimited text file format, see the link below:
#' https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
#' for metabolite enrichment analysis with KEGG's pathways
#' this file is embedded within the package:
#' "/metstats/inst/extdata/kegg_msea.txt"
#'
#' @param ranking_metric column name from input data
#' to be used as the metric for enrichment analysis
#'
#' @param metric_type type of the metric being used
#' for enrichment analysis 'abs' or 'non-abs'
#'
#' @param pval pvalue threshold to filter the enrichment results
#'
#' @returns a list that contains
#' 1) a tbl of enrichment statistics
#' 2) an enrichment dot plot
#' 3) a list of mountain plots
#'
#' @export
pathway_enrichment <- function(data_input,
                               id_column,
                               feature_sets_file,
                               ranking_metric,
                               metric_type = "non-abs",
                               pval = 1) {
  checkmate::assertString(id_column)
  checkmate::assertString(ranking_metric)
  checkmate::assertString(metric_type)

  ## check if the input data is a dataframe and has 'feature_name',
  ## 'id_column', and 'ranking_metric' columns
  if (any(class(data_input) == "data.frame")) {
    if (!all(c(id_column, ranking_metric) %in% colnames(data_input))) {
      stop("The id_column or ranking_metric is missing in the input data")
    }
  } else {
    stop("The input data must be a dataframe")
  }

  ## check if the feature_sets_file is a .txt file
  if (!grepl(".txt", feature_sets_file)) {
    stop("feature_sets_file must be a tab-delimited .txt file")
  }

  if (!metric_type %in% c("abs", "non-abs")) {
    stop("metric_type paramether must be either 'asb' or 'non-abs'")
  }

  ## remove features with missing id
  data_input <- data_input %>%
    dplyr::filter(!(is.na(data_input[id_column]) | data_input[id_column] == ""))

  ## create a ranked list from the input data and sort
  ranked_list <- data_input %>%
    dplyr::mutate(ranking_metric := !!dplyr::sym(ranking_metric)) %>%
    dplyr::select(id_column, ranking_metric)

  names(ranked_list)[names(ranked_list) == ranking_metric] <- "ranking_metric"
  ranked_list_mut <- ranked_list$ranking_metric
  names(ranked_list_mut) <- ranked_list[[id_column]]
  ranked_list_mut <- sort(ranked_list_mut, decreasing = TRUE)
  ranked_list_mut <- ranked_list_mut[!duplicated(names(ranked_list_mut))]

  ## check if there are duplicate features in the data
  if (any(duplicated(names(ranked_list_mut)))) {
    warning("There are duplicates in the list of features")
    ranked_list_mut <- ranked_list_mut[!duplicated(names(ranked_list_mut))]
  }

  if (!all(order(ranked_list_mut, decreasing = TRUE) == 1:length(ranked_list_mut))) {
    ranked_list_mut <- sort(ranked_list_mut, decreasing = TRUE)
  }

  ## generate a list of features sets from the .txt feature_sets_file
  enrichment_sets <- fgsea::gmtPathways(feature_sets_file)

  ## run enrichment analysis via fgsea package
  enrich_res <- fgsea::fgsea(
    pathways = enrichment_sets,
    stats = ranked_list_mut,
    nPermSimple = 100000
  ) %>%
    as.data.frame() %>%
    dplyr::filter(padj < !!pval)

  ## sort the enrichment results with NES (normalized enrichment score)
  enrich_res <- enrich_res %>% dplyr::arrange(dplyr::desc(NES))

  ## add a column to indicate direction of enrichment
  enrich_res$Enrichment <- ifelse(enrich_res$NES > 0,
    "Up-regulated",
    "Down-regulated"
  )

  ## update leadingEdge column of the enrichment table
  ## to show feature_name instead of id
  ## this would be useful for metabolite enrichment analysis
  ## where ids are not intuitive
  if ("feature_name" %in% colnames(data_input)) {
    tbl_temp <- enrich_res %>%
      tidyr::unnest(cols = leadingEdge) %>%
      dplyr::left_join(
        data_input %>%
          dplyr::select(feature_name, id_column) %>%
          unique(),
        by = c("leadingEdge" = id_column)
      ) %>%
      dplyr::mutate(
        name_or_id =
          ifelse(is.na(feature_name),
            leadingEdge,
            feature_name
          )
      ) %>%
      dplyr::select(-leadingEdge, -feature_name) %>%
      dplyr::rename(leadingEdge = name_or_id)

    enrichment_table <- dplyr::left_join(
      enrich_res %>%
        dplyr::select(-leadingEdge),
      stats::aggregate(
        leadingEdge ~ pathway,
        tbl_temp, paste
      ),
      by = "pathway"
    )
  } else {
    enrichment_table <- enrich_res
  }

  ## only keep positive enrichments if the metric being used is absolute
  if (metric_type == "abs") {
    enrichment_table <- enrichment_table[enrichment_table$NES > 0, ]
  }

  ## add size of pathway to the enrichment table
  enrichment_table <- enrichment_table %>%
    dplyr::left_join(
      tibble::tibble(pathway = names(enrichment_sets)) %>%
        rowwise() %>%
        mutate(setsize = length(enrichment_sets[[pathway]][enrichment_sets[[pathway]] != ""])),
      by = "pathway"
    ) %>%
    dplyr::relocate(setsize, .after = size)

  ## select top 10 up and top 10 down sets for plotting
  ## or select top 20 sets if the ranking metric being used is absolute
  if (metric_type == "abs") {
    enrich_sub <- enrich_res[enrich_res$NES > 0, ]
    enrich_sub <- utils::head(enrich_sub, n = 20)
  } else {
    enrich_sub <- rbind(
      utils::head(enrich_res, n = 10),
      utils::tail(enrich_res, n = 10)
    )
  }

  enrich_sub <- enrich_sub %>% dplyr::arrange(NES)
  enrich_sub$pathway <- factor(enrich_sub$pathway,
    levels = unique(enrich_sub$pathway)
  )

  plot_theme <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1.5),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black", size = 0.65),
    axis.text.x = ggplot2::element_text(size = 8, angle = 0, colour = "black", face = "bold"),
    axis.text.y = ggplot2::element_text(siz = 8, colour = "black", face = "bold"),
    axis.ticks.x = ggplot2::element_line(colour = "black"),
    axis.ticks.y = ggplot2::element_line(colour = "black"),
    axis.title.y = ggplot2::element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x = ggplot2::element_text(size = 12, face = "bold", colour = "black"),
    legend.text = ggplot2::element_text(size = 12, face = "bold"),
    plot.title = ggplot2::element_text(size = 12, face = "bold"),
    strip.text.x = ggplot2::element_text(size = 12, face = "bold")
  )

  ## plot a dot plot of all enriched sets
  if (metric_type == "abs") {
    enrichment_plot <- ggplot2::ggplot(
      enrich_sub,
      ggplot2::aes(
        x = NES,
        y = pathway,
        size = -log10(pval)
      )
    ) +
      ggplot2::geom_point(ggplot2::aes(colour = NES, size = -log10(pval))) +
      ggplot2::scale_color_gradient2(
        low = "white",
        high = "red", space = "Lab",
        limits = c(0, plyr::round_any(max(abs(enrich_sub$NES)), 0.5, f = ceiling))
      ) +
      ggplot2::labs(y = "Pathway", x = "Normalized Enrichment Score", title = ranking_metric) +
      plot_theme
  } else {
    enrichment_plot <- ggplot2::ggplot(
      enrich_sub,
      ggplot2::aes(x = NES, y = pathway, size = -log10(pval))
    ) +
      ggplot2::geom_point(ggplot2::aes(colour = NES, size = -log10(pval))) +
      ggplot2::scale_color_gradient2(
        low = "blue", mid = "white",
        high = "red", space = "Lab",
        limits = c(
          -plyr::round_any(max(abs(enrich_sub$NES)), 0.5, f = ceiling),
          plyr::round_any(max(abs(enrich_sub$NES)), 0.5, f = ceiling)
        )
      ) +
      ggplot2::labs(y = "Pathway", x = "Normalized Enrichment Score", title = ranking_metric) +
      plot_theme
  }

  # generate mountain plots of enriched sets
  mountain_plot <- list()
  mountain_plot <- lapply(enrich_sub$pathway, function(pathway) {
    p <- fgsea::plotEnrichment(enrichment_sets[[pathway]], ranked_list_mut) +
      ggplot2::ggtitle(pathway, ranking_metric)
    return(p)
  })

  names(mountain_plot) <- enrich_sub$pathway

  output <- list(
    "enrichment_table" = enrichment_table,
    "enrichment_plot" = enrichment_plot,
    "mountain_plot" = mountain_plot
  )

  return(output)
}
