library(shiny)
library(bslib)
library(plotly)

# ============================================================
# FUNCTIONAL ENRICHMENT: DATA AND HELPER FUNCTIONS
# ============================================================

candidate_choices <- c(
  "FLRT3", "MDK", "MMP11", "SALL4", "TGFB2", "VTN"
)

read_csv_safe <- function(filename) {
  if (file.exists(filename)) {
    tryCatch(
      read.csv(
        filename,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }
}


# Remove species labels appended to annotation names in plot labels only.
# The original annotation data and CSV files are not modified.
clean_annotation_labels <- function(x) {
  x <- as.character(x)
  x <- gsub(
    "\\s*\\(?\\s*Homo sapiens\\s*\\)?\\s*",
    " ",
    x,
    ignore.case = TRUE
  )
  x <- gsub(
    "\\s*\\(?\\s*humans?\\s*\\)?\\s*",
    " ",
    x,
    ignore.case = TRUE
  )
  trimws(gsub("\\s+", " ", x))
}


# Prepare the supplied Reactome gene-pathway records for the pathway plots.
# The individual-gene file provides the gene-to-pathway relationships;
# the hierarchy file supplies the Reactome hierarchy used by the sunburst.
make_reactome_plot_data <- function(individual_data, hierarchy_data, candidates) {

  required_individual <- c(
    "gene", "pathway_id", "pathway_name",
    "pathway_hierarchy", "top_level_pathway"
  )

  if (is.null(individual_data) ||
      !all(required_individual %in% names(individual_data))) {
    return(data.frame())
  }

  plot_data <- individual_data[
    !is.na(individual_data$gene) &
      trimws(as.character(individual_data$gene)) %in% candidates,
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) == 0) {
    return(data.frame())
  }

  # Use the separately supplied hierarchy file as the hierarchy source.
  if (!is.null(hierarchy_data) &&
      all(c("pathway_id", "pathway_hierarchy", "top_level_pathway") %in%
          names(hierarchy_data))) {

    hierarchy_index <- match(
      as.character(plot_data$pathway_id),
      as.character(hierarchy_data$pathway_id)
    )

    matched_hierarchy <- !is.na(hierarchy_index)

    plot_data$pathway_hierarchy[matched_hierarchy] <-
      as.character(hierarchy_data$pathway_hierarchy[hierarchy_index[matched_hierarchy]])

    plot_data$top_level_pathway[matched_hierarchy] <-
      as.character(hierarchy_data$top_level_pathway[hierarchy_index[matched_hierarchy]])
  }

  result <- data.frame(
    candidate = trimws(as.character(plot_data$gene)),
    annotation = clean_annotation_labels(plot_data$pathway_name),
    annotation_id = trimws(as.character(plot_data$pathway_id)),
    pathway_hierarchy = as.character(plot_data$pathway_hierarchy),
    top_level_pathway = as.character(plot_data$top_level_pathway),
    stringsAsFactors = FALSE
  )

  result <- unique(result)
  result[
    !is.na(result$candidate) &
      !is.na(result$annotation) &
      !is.na(result$annotation_id) &
      nzchar(result$candidate) &
      nzchar(result$annotation) &
      nzchar(result$annotation_id),
    ,
    drop = FALSE
  ]
}

GO_data <- read_csv_safe("GO_analysis.csv")
KEGG_data <- read_csv_safe("KEGG_analysis.csv")
focused_data <- read_csv_safe("focused_unique_features.csv")
Reactome_individual_data <- read_csv_safe("Reactome_individual_gene_pathways.csv")
Reactome_hierarchy_data <- read_csv_safe("Reactome_hierarchy_data.csv")

reactome_top_level_choices <- if (
  !is.null(Reactome_individual_data) &&
    "top_level_pathway" %in% names(Reactome_individual_data)
) {
  available_top_levels <- sort(unique(na.omit(trimws(as.character(
    Reactome_individual_data$top_level_pathway
  )))))
  c("All top-level pathways" = "__ALL__", stats::setNames(
    available_top_levels,
    available_top_levels
  ))
} else {
  c("All top-level pathways" = "__ALL__")
}


# ============================================================
# DIRECT / INDIRECT REGENERATION RELEVANCE CLASSIFICATION
# These are curated interpretation labels for the supplied
# annotation records, not statistical or causal evidence.
# ============================================================

# GO terms are classified at candidate + GO-ID level because the
# same broad GO term can have different relevance in different genes.
direct_go_pairs <- c(
  "FLRT3::GO:0098742",
  "FLRT3::GO:0005925",

  "MDK::GO:0007162",
  "MDK::GO:0009611",
  "MDK::GO:0010718",
  "MDK::GO:0010838",
  "MDK::GO:0030335",
  "MDK::GO:0032330",
  "MDK::GO:0042246",
  "MDK::GO:0045785",
  "MDK::GO:0046850",
  "MDK::GO:0051781",
  "MDK::GO:0061036",
  "MDK::GO:1900026",
  "MDK::GO:1904036",

  "MMP11::GO:0022617",
  "MMP11::GO:0030198",
  "MMP11::GO:0030199",
  "MMP11::GO:0030574",
  "MMP11::GO:0071711",
  "MMP11::GO:0031012",

  "SALL4::GO:0030326",
  "SALL4::GO:0035019",

  "TGFB2::GO:0001501",
  "TGFB2::GO:0001837",
  "TGFB2::GO:0008284",
  "TGFB2::GO:0008285",
  "TGFB2::GO:0009611",
  "TGFB2::GO:0010634",
  "TGFB2::GO:0010718",
  "TGFB2::GO:0016477",
  "TGFB2::GO:0030199",
  "TGFB2::GO:0030307",
  "TGFB2::GO:0030308",
  "TGFB2::GO:0030326",
  "TGFB2::GO:0033630",
  "TGFB2::GO:0042060",
  "TGFB2::GO:0042127",
  "TGFB2::GO:0045216",
  "TGFB2::GO:0045787",
  "TGFB2::GO:0048103",
  "TGFB2::GO:0050680",
  "TGFB2::GO:0051781",
  "TGFB2::GO:0060391",
  "TGFB2::GO:0062023",

  "VTN::GO:0007155",
  "VTN::GO:0007160",
  "VTN::GO:0010811",
  "VTN::GO:0016477",
  "VTN::GO:0030155",
  "VTN::GO:0030198",
  "VTN::GO:0033627",
  "VTN::GO:0061302",
  "VTN::GO:0090303",
  "VTN::GO:0097421",
  "VTN::GO:0005604",
  "VTN::GO:0062023",
  "VTN::GO:0071062",
  "VTN::GO:0098637",
  "VTN::GO:0005178",
  "VTN::GO:0005201",
  "VTN::GO:0005518",
  "VTN::GO:0050840"
)

# KEGG pathways with a close mechanistic connection to tissue repair,
# growth, cell behaviour or extracellular-matrix interactions.
direct_kegg_ids <- c(
  "hsa04010", # MAPK signalling
  "hsa04110", # Cell cycle
  "hsa04350", # TGF-beta signalling
  "hsa04380", # Osteoclast differentiation / tissue remodelling
  "hsa04390", # Hippo signalling
  "hsa04510", # Focal adhesion
  "hsa04512", # ECM-receptor interaction
  "hsa04518", # Integrin signalling
  "hsa04151"  # PI3K-Akt signalling
)

classify_go_relevance <- function(data) {
  if (is.null(data) || nrow(data) == 0 ||
      !all(c("candidate", "annotation_id") %in% names(data))) {
    return(data.frame())
  }

  data$regeneration_relevance <- "Indirect"
  relevance_key <- paste(
    trimws(as.character(data$candidate)),
    trimws(as.character(data$annotation_id)),
    sep = "::"
  )
  data$regeneration_relevance[relevance_key %in% direct_go_pairs] <- "Direct"
  data
}

classify_kegg_relevance <- function(data) {
  if (is.null(data) || nrow(data) == 0 ||
      !"annotation_id" %in% names(data)) {
    return(data.frame())
  }

  data$regeneration_relevance <- ifelse(
    trimws(as.character(data$annotation_id)) %in% direct_kegg_ids,
    "Direct",
    "Indirect"
  )
  data
}

classify_focused_pathway_relevance <- function(data) {
  if (is.null(data) || nrow(data) == 0 ||
      !all(c("annotation_type", "annotation_id") %in% names(data))) {
    return(data.frame())
  }

  is_kegg <- tolower(trimws(as.character(data$annotation_type))) == "kegg"
  data$regeneration_relevance <- "Indirect"
  data$regeneration_relevance[
    is_kegg & trimws(as.character(data$annotation_id)) %in% direct_kegg_ids
  ] <- "Direct"
  data
}


get_selected_candidates <- function(data, selected) {
  if (is.null(data) || nrow(data) == 0 || is.null(selected)) {
    return(data.frame())
  }

  if (!"candidate" %in% names(data)) {
    return(data.frame())
  }

  data <- data[
    !is.na(data$candidate) &
      data$candidate %in% selected,
    ,
    drop = FALSE
  ]

  data
}


filter_focused_annotations <- function(data, selected) {
  if (is.null(data) || nrow(data) == 0 || is.null(selected)) {
    return(data.frame())
  }

  if (!"genes" %in% names(data)) {
    return(data.frame())
  }

  data$genes <- as.character(data$genes)

  keep <- vapply(
    data$genes,
    function(gene_string) {
      if (is.na(gene_string) || !nzchar(gene_string)) {
        return(FALSE)
      }

      genes <- trimws(
        unlist(
          strsplit(gene_string, "[,;|]")
        )
      )

      any(genes %in% selected)
    },
    logical(1)
  )

  data[keep, , drop = FALSE]
}


# Convert focused pathway records into candidate-annotation pairs.
# Each selected candidate gets its own row and network node.

make_focused_candidate_data <- function(data, selected) {

  if (is.null(data) || nrow(data) == 0) {
    return(data.frame())
  }

  if (!all(c("genes", "annotation") %in% names(data))) {
    return(data.frame())
  }

  result <- lapply(seq_len(nrow(data)), function(i) {

    gene_string <- as.character(data$genes[i])

    if (is.na(gene_string) || !nzchar(gene_string)) {
      return(NULL)
    }

    genes <- trimws(
      unlist(
        strsplit(gene_string, "[,;|]")
      )
    )

    matching_genes <- intersect(
      genes,
      selected
    )

    if (length(matching_genes) == 0) {
      return(NULL)
    }

    data.frame(
      candidate = matching_genes,
      annotation = rep(
        as.character(data$annotation[i]),
        length(matching_genes)
      ),
      stringsAsFactors = FALSE
    )
  })

  result <- result[
    !vapply(result, is.null, logical(1))
  ]

  if (length(result) == 0) {
    return(data.frame(
      candidate = character(),
      annotation = character()
    ))
  }

  result <- do.call(rbind, result)

  result <- unique(result)

  rownames(result) <- NULL

  result
}


# ============================================================
# GENERAL ANNOTATION BAR PLOT
# Used by GO and Pathways sections.
# ============================================================

make_annotation_bar <- function(data, selected, label_column = "annotation") {

  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }

  if (!all(c("candidate", label_column) %in% names(data))) {
    return(NULL)
  }

  plot_data <- data[
    !is.na(data[[label_column]]) &
      nzchar(trimws(as.character(data[[label_column]]))),
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  counts <- aggregate(
    candidate ~ annotation,
    data = data.frame(
      annotation = clean_annotation_labels(
        plot_data[[label_column]]
      ),
      candidate = as.character(plot_data$candidate)
    ),
    FUN = length
  )

  names(counts)[2] <- "annotation_count"

  counts <- counts[
    order(counts$annotation_count),
    ,
    drop = FALSE
  ]

  counts$annotation <- factor(
    counts$annotation,
    levels = counts$annotation
  )

  plot_ly(
    counts,
    x = ~annotation_count,
    y = ~annotation,
    type = "bar",
    orientation = "h",
    marker = list(color = "#B98791"),
    hovertemplate = paste(
      "<b>%{y}</b>",
      "<br>Number of annotation records: %{x}",
      "<extra></extra>"
    )
  ) %>%
    layout(
      xaxis = list(
        title = "Number of annotation records",
        rangemode = "tozero"
      ),
      yaxis = list(
        title = "",
        automargin = TRUE
      ),
      margin = list(
        l = 180,
        r = 30,
        t = 30,
        b = 55
      ),
      paper_bgcolor = "#FFFDFC",
      plot_bgcolor = "#FFFDFC"
    )
}


# ============================================================
# KEGG BAR PLOT
# Default plot displayed when the KEGG section opens.
# ============================================================

make_kegg_bar <- function(data) {

  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }

  if (!all(c("candidate", "annotation") %in% names(data))) {
    return(NULL)
  }

  plot_data <- data[
    !is.na(data$annotation) &
      nzchar(trimws(as.character(data$annotation))),
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  plot_data$annotation <- clean_annotation_labels(
    plot_data$annotation
  )

  kegg_counts <- as.data.frame(
    table(plot_data$annotation),
    stringsAsFactors = FALSE
  )

  names(kegg_counts) <- c(
    "annotation",
    "annotation_count"
  )

  kegg_counts <- kegg_counts[
    kegg_counts$annotation_count > 0,
    ,
    drop = FALSE
  ]

  if (nrow(kegg_counts) == 0) {
    return(NULL)
  }

  kegg_counts <- kegg_counts[
    order(kegg_counts$annotation_count),
    ,
    drop = FALSE
  ]

  kegg_counts$annotation <- factor(
    kegg_counts$annotation,
    levels = kegg_counts$annotation
  )

  plot_ly(
    kegg_counts,
    x = ~annotation_count,
    y = ~annotation,
    type = "bar",
    orientation = "h",
    marker = list(
      color = "#B98791",
      line = list(
        color = "#8F568B",
        width = 0.5
      )
    ),
    hovertemplate = paste(
      "<b>%{y}</b>",
      "<br>Number of KEGG annotation records: %{x}",
      "<extra></extra>"
    )
  ) %>%
    layout(
      title = list(
        text = "KEGG Pathway Annotation Records",
        x = 0.02,
        font = list(
          size = 16,
          color = "#242020"
        )
      ),
      xaxis = list(
        title = "Number of annotation records",
        rangemode = "tozero"
      ),
      yaxis = list(
        title = "",
        automargin = TRUE
      ),
      margin = list(
        l = 190,
        r = 30,
        t = 65,
        b = 60
      ),
      paper_bgcolor = "#FFFDFC",
      plot_bgcolor = "#FFFDFC",
      font = list(color = "#242020")
    )
}


# ============================================================
# ANNOTATION NETWORK PLOT
# Bipartite Sankey network:
# Candidate genes -> Annotation terms
# ============================================================

make_annotation_network <- function(data) {

  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }

  if (!all(c("candidate", "annotation") %in% names(data))) {
    return(NULL)
  }

  plot_data <- data.frame(
    candidate = trimws(as.character(data$candidate)),
    annotation = clean_annotation_labels(data$annotation),
    stringsAsFactors = FALSE
  )

  plot_data <- plot_data[
    !is.na(plot_data$candidate) &
      !is.na(plot_data$annotation) &
      nzchar(plot_data$candidate) &
      nzchar(plot_data$annotation),
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  edge_data <- aggregate(
    list(value = rep(1, nrow(plot_data))),
    by = list(
      candidate = plot_data$candidate,
      annotation = plot_data$annotation
    ),
    FUN = sum
  )

  if (nrow(edge_data) == 0) {
    return(NULL)
  }

  candidates <- unique(edge_data$candidate)
  annotations <- unique(edge_data$annotation)

  candidate_ids <- paste0(
    "candidate_",
    seq_along(candidates)
  )

  annotation_ids <- paste0(
    "annotation_",
    seq_along(annotations)
  )

  nodes <- data.frame(
    node_id = c(
      candidate_ids,
      annotation_ids
    ),
    label = c(
      candidates,
      annotations
    ),
    node_type = c(
      rep("candidate", length(candidates)),
      rep("annotation", length(annotations))
    ),
    stringsAsFactors = FALSE
  )

  edge_data$source <- match(
    edge_data$candidate,
    candidates
  ) - 1L

  edge_data$target <- match(
    edge_data$annotation,
    annotations
  ) - 1L + length(candidates)

  if (
    anyNA(edge_data$source) ||
    anyNA(edge_data$target)
  ) {
    return(NULL)
  }

  plot_ly(
    type = "sankey",
    orientation = "h",
    arrangement = "snap",

    node = list(
      label = nodes$label,
      pad = 18,
      thickness = 20,
      line = list(
        color = "#B98791",
        width = 0.7
      ),
      color = ifelse(
        nodes$node_type == "candidate",
        "#DCE8DC",
        "#E8D5D2"
      ),
      hovertemplate = paste(
        "<b>%{label}</b>",
        "<extra></extra>"
      )
    ),

    link = list(
      source = as.integer(edge_data$source),
      target = as.integer(edge_data$target),
      value = as.numeric(edge_data$value),
      color = "rgba(185,135,145,0.40)",
      hovertemplate = paste(
        "Connection strength: %{value}",
        "<extra></extra>"
      )
    )
  ) %>%
    layout(
      font = list(
        size = 11,
        color = "#242020"
      ),
      paper_bgcolor = "#FFFDFC",
      plot_bgcolor = "#FFFDFC",
      margin = list(
        l = 15,
        r = 15,
        t = 25,
        b = 25
      )
    )
}


# ============================================================
# GO GENE TREEMAP
# GO category -> annotation term for the selected candidate gene.
# Each unique GO annotation is one leaf; tile area represents one
# annotation record (the supplied file has no term-level scores).
# ============================================================

make_go_gene_treemap <- function(data) {

  required_columns <- c(
    "candidate",
    "category",
    "annotation_id",
    "annotation"
  )

  if (
    is.null(data) ||
      nrow(data) == 0 ||
      !all(required_columns %in% names(data))
  ) {
    return(NULL)
  }

  plot_data <- data.frame(
    candidate = trimws(as.character(data$candidate)),
    category = toupper(trimws(as.character(data$category))),
    annotation_id = trimws(as.character(data$annotation_id)),
    annotation = clean_annotation_labels(data$annotation),
    stringsAsFactors = FALSE
  )

  plot_data <- plot_data[
    !is.na(plot_data$candidate) &
      !is.na(plot_data$category) &
      !is.na(plot_data$annotation_id) &
      !is.na(plot_data$annotation) &
      nzchar(plot_data$candidate) &
      nzchar(plot_data$category) &
      nzchar(plot_data$annotation_id) &
      nzchar(plot_data$annotation),
    ,
    drop = FALSE
  ]

  # A GO term should appear only once for a candidate.
  plot_data <- unique(plot_data)

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  category_labels <- c(
    BP = "Biological Process",
    CC = "Cellular Component",
    MF = "Molecular Function"
  )

  category_colors <- c(
    BP = "#B98791",
    CC = "#8F568B",
    MF = "#D6A66B"
  )

  plot_data <- plot_data[
    plot_data$category %in% names(category_labels),
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  plot_data$category_label <- unname(
    category_labels[plot_data$category]
  )

  plot_data$tile_color <- unname(
    category_colors[plot_data$category]
  )

  selected_gene <- unique(plot_data$candidate)
  selected_gene <- paste(selected_gene, collapse = ", ")

  category_order <- names(category_labels)[
    names(category_labels) %in% unique(plot_data$category)
  ]

  parent_rows <- lapply(category_order, function(category_code) {

    category_data <- plot_data[
      plot_data$category == category_code,
      ,
      drop = FALSE
    ]

    data.frame(
      id = paste0("category_", category_code),
      label = unname(category_labels[category_code]),
      parent = "",
      value = nrow(category_data),
      color = unname(category_colors[category_code]),
      hover = paste0(
        "<b>", unname(category_labels[category_code]), "</b>",
        "<br>Candidate gene: ", selected_gene,
        "<br>Unique GO annotations: ", nrow(category_data)
      ),
      stringsAsFactors = FALSE
    )
  })

  parent_data <- do.call(rbind, parent_rows)

  leaf_data <- data.frame(
    id = paste0(
      "term_",
      make.unique(
        paste(plot_data$candidate, plot_data$annotation_id, sep = "_")
      )
    ),
    label = vapply(
      plot_data$annotation,
      function(term) {
        term <- trimws(as.character(term))
        if (is.na(term) || !nzchar(term)) return("")

        words <- strsplit(term, "\\s+")[[1]]
        lines <- character(0)
        current_line <- ""

        # Wrap only between words so each word remains intact.
        # A line is kept to approximately 15 characters for compact tiles.
        for (word in words) {
          candidate <- if (nzchar(current_line)) {
            paste(current_line, word)
          } else {
            word
          }

          if (nchar(candidate, type = "width") > 15 && nzchar(current_line)) {
            lines <- c(lines, current_line)
            current_line <- word
          } else {
            current_line <- candidate
          }
        }

        lines <- c(lines, current_line)
        paste(lines, collapse = "<br>")
      },
      character(1)
    ),
    parent = paste0("category_", plot_data$category),
    value = 1,
    color = plot_data$tile_color,
    hover = paste0(
      "<b>", plot_data$annotation, "</b>",
      "<br>Candidate gene: ", plot_data$candidate,
      "<br>GO category: ", plot_data$category_label,
      "<br>GO ID: ", plot_data$annotation_id
    ),
    stringsAsFactors = FALSE
  )

  treemap_data <- rbind(
    parent_data,
    leaf_data
  )

  plot_ly(
    type = "treemap",
    ids = treemap_data$id,
    labels = treemap_data$label,
    parents = treemap_data$parent,
    values = treemap_data$value,
    branchvalues = "total",
    textinfo = "label",
    textfont = list(
      size = 18,
      color = "#FFFFFF",
      family = "Arial"
    ),
    insidetextorientation = "horizontal",
    tiling = list(
      packing = "squarify",
      pad = 4
    ),
    marker = list(
      colors = treemap_data$color,
      line = list(
        color = "#FFFDFC",
        width = 1.5
      )
    ),
    customdata = treemap_data$hover,
    hovertemplate = "%{customdata}<extra></extra>"
  ) %>%
    layout(
      title = list(
        text = paste("GO Annotation Treemap:", selected_gene),
        x = 0.02,
        xanchor = "left",
        font = list(
          size = 20,
          color = "#242020"
        )
      ),
      uniformtext = list(
        minsize = 13,
        mode = "show"
      ),
      margin = list(
        l = 10,
        r = 10,
        t = 60,
        b = 10
      ),
      paper_bgcolor = "#FFFDFC",
      plot_bgcolor = "#FFFDFC"
    )
}


# ============================================================
# ANNOTATION HEATMAP
# Retained for GO and Pathways sections.
# ============================================================

make_annotation_heatmap <- function(data) {

  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }

  if (!all(c("candidate", "annotation") %in% names(data))) {
    return(NULL)
  }

  plot_data <- data.frame(
    candidate = trimws(as.character(data$candidate)),
    annotation = clean_annotation_labels(data$annotation),
    stringsAsFactors = FALSE
  )

  plot_data <- plot_data[
    !is.na(plot_data$candidate) &
      !is.na(plot_data$annotation) &
      nzchar(plot_data$candidate) &
      nzchar(plot_data$annotation),
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  count_table <- as.data.frame(
    table(
      candidate = plot_data$candidate,
      annotation = plot_data$annotation
    ),
    stringsAsFactors = FALSE
  )

  count_table <- count_table[
    count_table$Freq > 0,
    ,
    drop = FALSE
  ]

  if (nrow(count_table) == 0) {
    return(NULL)
  }

  heatmap_data <- xtabs(
    Freq ~ candidate + annotation,
    data = count_table
  )

  plot_ly(
    x = colnames(heatmap_data),
    y = rownames(heatmap_data),
    z = unname(heatmap_data),
    type = "heatmap",
    colors = c(
      "#F7F0E8",
      "#B98791",
      "#8F568B"
    ),
    hovertemplate = paste(
      "Candidate: %{y}",
      "<br>Annotation: %{x}",
      "<br>Records: %{z}",
      "<extra></extra>"
    )
  ) %>%
    layout(
      xaxis = list(
        title = "",
        tickangle = -45,
        automargin = TRUE
      ),
      yaxis = list(
        title = "",
        automargin = TRUE
      ),
      margin = list(
        l = 80,
        r = 20,
        t = 20,
        b = 130
      ),
      paper_bgcolor = "#FFFDFC",
      plot_bgcolor = "#FFFDFC"
    )
}


# ============================================================
# PATHWAY CHORD DIAGRAM
# Candidate genes are connected to their associated pathway terms.
# Link width represents the number of annotation records.
# ============================================================

make_pathway_chord <- function(data) {

  required_columns <- c("candidate", "annotation", "annotation_id")
  if (is.null(data) || nrow(data) == 0 ||
      !all(required_columns %in% names(data))) {
    return(NULL)
  }

  plot_data <- unique(data.frame(
    candidate = trimws(as.character(data$candidate)),
    annotation = clean_annotation_labels(data$annotation),
    annotation_id = trimws(as.character(data$annotation_id)),
    stringsAsFactors = FALSE
  ))

  plot_data <- plot_data[
    !is.na(plot_data$candidate) &
      !is.na(plot_data$annotation) &
      !is.na(plot_data$annotation_id) &
      nzchar(plot_data$candidate) &
      nzchar(plot_data$annotation_id),
    ,
    drop = FALSE
  ]
  if (nrow(plot_data) == 0) return(NULL)

  edge_data <- unique(data.frame(
    candidate = plot_data$candidate,
    pathway_id = plot_data$annotation_id,
    pathway_name = plot_data$annotation,
    candidate_node = paste0("gene::", plot_data$candidate),
    pathway_node = paste0("pathway::", plot_data$annotation_id),
    stringsAsFactors = FALSE
  ))

  candidate_nodes <- unique(data.frame(
    node_id = paste0("gene::", edge_data$candidate),
    label = edge_data$candidate,
    node_type = "Candidate gene",
    detail = edge_data$candidate,
    stringsAsFactors = FALSE
  ))

  pathway_nodes <- unique(data.frame(
    node_id = paste0("pathway::", edge_data$pathway_id),
    label = edge_data$pathway_name,
    node_type = "Reactome pathway",
    detail = paste0(edge_data$pathway_name, " (", edge_data$pathway_id, ")"),
    stringsAsFactors = FALSE
  ))

  nodes <- rbind(candidate_nodes, pathway_nodes)
  n_nodes <- nrow(nodes)
  if (n_nodes < 2) return(NULL)

  angles <- seq(0, 2 * pi, length.out = n_nodes + 1)[-(n_nodes + 1)]
  node_x <- cos(angles)
  node_y <- sin(angles)
  names(node_x) <- nodes$node_id
  names(node_y) <- nodes$node_id

  pathway_degree <- table(edge_data$pathway_node)
  p <- plot_ly()

  for (i in seq_len(nrow(edge_data))) {
    x0 <- unname(node_x[edge_data$candidate_node[i]])
    y0 <- unname(node_y[edge_data$candidate_node[i]])
    x1 <- unname(node_x[edge_data$pathway_node[i]])
    y1 <- unname(node_y[edge_data$pathway_node[i]])

    t <- seq(0, 1, length.out = 40)
    curve_x <- (1 - t)^2 * x0 + 2 * (1 - t) * t * 0 + t^2 * x1
    curve_y <- (1 - t)^2 * y0 + 2 * (1 - t) * t * 0 + t^2 * y1

    shared_gene_count <- as.numeric(pathway_degree[edge_data$pathway_node[i]])
    link_width <- 1 + 1.2 * sqrt(shared_gene_count)

    p <- add_trace(
      p,
      x = curve_x,
      y = curve_y,
      type = "scatter",
      mode = "lines",
      line = list(color = "rgba(143,86,139,0.38)", width = link_width),
      text = paste0(
        edge_data$candidate[i], " → ",
        edge_data$pathway_name[i], " (", edge_data$pathway_id[i], ")"
      ),
      hovertemplate = "%{text}<extra></extra>",
      showlegend = FALSE
    )
  }

  node_hover <- paste0("<b>", nodes$detail, "</b><br>", nodes$node_type)

  p <- add_trace(
    p,
    x = node_x,
    y = node_y,
    type = "scatter",
    mode = "markers+text",
    text = ifelse(nodes$node_type == "Candidate gene", nodes$label, ""),
    textposition = "top center",
    textfont = list(size = 11, color = "#242020"),
    marker = list(
      size = ifelse(nodes$node_type == "Candidate gene", 15, 8),
      color = ifelse(nodes$node_type == "Candidate gene", "#DCE8DC", "#B98791"),
      line = list(color = "#59494A", width = 0.8)
    ),
    hovertext = node_hover,
    hovertemplate = "%{hovertext}<extra></extra>",
    showlegend = FALSE
  )

  p %>% layout(
    xaxis = list(visible = FALSE, range = c(-1.4, 1.4), fixedrange = TRUE),
    yaxis = list(
      visible = FALSE, range = c(-1.4, 1.4), fixedrange = TRUE,
      scaleanchor = "x", scaleratio = 1
    ),
    margin = list(l = 25, r = 25, t = 25, b = 25),
    paper_bgcolor = "#FFFDFC",
    plot_bgcolor = "#FFFDFC",
    showlegend = FALSE,
    annotations = list(list(
      x = 0, y = -1.31, xref = "x", yref = "y",
      text = "Candidate genes are linked to their Reactome pathways; hover over nodes or links for details.",
      showarrow = FALSE,
      font = list(size = 10, color = "#59494A")
    ))
  )
}


# ============================================================
# PATHWAY SUNBURST
# Reactome hierarchy -> pathway -> associated candidate gene.
# Segment size represents the number of gene-pathway associations.
# ============================================================

make_pathway_sunburst <- function(data) {

  required_columns <- c(
    "candidate", "annotation", "annotation_id",
    "pathway_hierarchy", "top_level_pathway"
  )
  if (is.null(data) || nrow(data) == 0 ||
      !all(required_columns %in% names(data))) {
    return(NULL)
  }

  plot_data <- unique(data.frame(
    candidate = trimws(as.character(data$candidate)),
    annotation = clean_annotation_labels(data$annotation),
    annotation_id = trimws(as.character(data$annotation_id)),
    pathway_hierarchy = as.character(data$pathway_hierarchy),
    top_level_pathway = trimws(as.character(data$top_level_pathway)),
    stringsAsFactors = FALSE
  ))

  plot_data <- plot_data[
    !is.na(plot_data$candidate) &
      !is.na(plot_data$annotation) &
      !is.na(plot_data$annotation_id) &
      !is.na(plot_data$top_level_pathway) &
      nzchar(plot_data$candidate) &
      nzchar(plot_data$annotation_id) &
      nzchar(plot_data$top_level_pathway),
    ,
    drop = FALSE
  ]
  if (nrow(plot_data) == 0) return(NULL)

  node_records <- list()
  record_index <- 0L

  add_sunburst_node <- function(id, label, parent) {
    record_index <<- record_index + 1L
    node_records[[record_index]] <<- data.frame(
      id = id, label = label, parent = parent, value = 1,
      stringsAsFactors = FALSE
    )
  }

  root_id <- "reactome_root"

  for (i in seq_len(nrow(plot_data))) {
    top_level <- plot_data$top_level_pathway[i]
    hierarchy <- plot_data$pathway_hierarchy[i]
    pathway_name <- plot_data$annotation[i]
    pathway_id <- plot_data$annotation_id[i]
    candidate <- plot_data$candidate[i]

    hierarchy_parts <- if (!is.na(hierarchy) && nzchar(trimws(hierarchy))) {
      trimws(strsplit(hierarchy, "\\s*>\\s*")[[1]])
    } else {
      character(0)
    }

    if (length(hierarchy_parts) > 0 &&
        identical(hierarchy_parts[1], top_level)) {
      hierarchy_parts <- hierarchy_parts[-1]
    }

    if (length(hierarchy_parts) > 0 &&
        identical(tail(hierarchy_parts, 1), pathway_name)) {
      hierarchy_parts <- head(hierarchy_parts, -1)
    }

    parent_id <- root_id
    top_id <- paste0("top::", top_level)
    add_sunburst_node(top_id, top_level, parent_id)
    parent_id <- top_id

    if (length(hierarchy_parts) > 0) {
      for (part in hierarchy_parts) {
        if (!is.na(part) && nzchar(part)) {
          part_id <- paste0(parent_id, "::", part)
          add_sunburst_node(part_id, part, parent_id)
          parent_id <- part_id
        }
      }
    }

    pathway_node_id <- paste0(parent_id, "::pathway::", pathway_id)
    add_sunburst_node(pathway_node_id, pathway_name, parent_id)

    gene_node_id <- paste0(pathway_node_id, "::gene::", candidate)
    add_sunburst_node(gene_node_id, candidate, pathway_node_id)
  }

  records <- do.call(rbind, node_records)
  node_values <- aggregate(
    value ~ id + label + parent,
    data = records,
    FUN = sum
  )

  sunburst_data <- rbind(
    data.frame(
      id = root_id,
      label = "Reactome pathway hierarchy",
      parent = "",
      value = nrow(plot_data),
      stringsAsFactors = FALSE
    ),
    node_values
  )

  plot_ly(
    data = sunburst_data,
    type = "sunburst",
    ids = ~id,
    labels = ~label,
    parents = ~parent,
    values = ~value,
    branchvalues = "total",
    insidetextorientation = "radial",
    hovertemplate = paste(
      "<b>%{label}</b>",
      "<br>Candidate-pathway associations: %{value}",
      "<extra></extra>"
    ),
    marker = list(line = list(color = "#FFFDFC", width = 1))
  ) %>%
    layout(
      margin = list(l = 10, r = 10, t = 20, b = 10),
      paper_bgcolor = "#FFFDFC",
      plot_bgcolor = "#FFFDFC"
    )
}


# ============================================================
# USER INTERFACE
# ============================================================

ui <- page_sidebar(
  title = div(
    style = "display:flex;align-items:center;gap:10px;",
    tags$img(
      src = "axolotl.png",
      style = "width:48px;height:48px;object-fit:contain;border-radius:50%;"
    ),
    div(
      span(
        "ReForm",
        style = "font-weight:850;letter-spacing:1px;color:#242020;font-size:25px;"
      ),
      br(),
      span(
        "The science of becoming whole again.",
        style = "font-size:12px;color:#59494A;"
      )
    )
  ),

  theme = bs_theme(
    version = 5,
    bg = "#F7F0E8",
    fg = "#242020",
    primary = "#C28F99",
    secondary = "#DCE8DC",
    base_font = font_google("DM Sans"),
    heading_font = font_google("DM Sans")
  ),

  sidebar = sidebar(
    width = 263,
    bg = "#E8D5D2",
    fg = "#242020",

    div(
      style = "padding:10px 4px 20px;",
      p(
        "AXOLOTL LIMB REGENERATION",
        style = paste(
          "color:#59494A;font-size:10px;",
          "letter-spacing:2px;font-weight:800;"
        )
      ),
      p(
        "Explore the molecular story behind regeneration.",
        style = "color:#242020;font-size:13px;"
      )
    ),

    navset_pill_list(
      id = "page",
      widths = c(12, 12),

      nav_panel("Home", value = "home", icon = icon("house")),
      nav_panel("About the Project", value = "about", icon = icon("book-open")),
      nav_panel("Differential Expression", value = "de", icon = icon("chart-bar")),
      nav_panel("Gene Mapping", value = "genes", icon = icon("dna")),
      nav_panel("Functional Enrichment", value = "enrichment", icon = icon("project-diagram")),
      nav_panel("Results & Interpretation", value = "deg", icon = icon("search"))
    ),

    br(),

    p(
      "ReForm | RNA-seq exploration",
      style = "font-size:11px;color:#59494A;"
    )
  ),

  tags$head(
    tags$style(HTML("
      body {
        background: #F7F0E8;
        color: #242020;
      }

      html,
      body {
        width: 100%;
        max-width: 100%;
        overflow-x: hidden;
      }

      .bslib-page-sidebar,
      .bslib-sidebar-layout {
        width: 100%;
        max-width: 100%;
        min-width: 0;
      }

      .bslib-sidebar-layout > .main,
      .bslib-sidebar-layout > .main > .container-fluid,
      .bslib-sidebar-layout > .main > .container {
        width: 100%;
        max-width: 100%;
        min-width: 0;
      }

      .main {
        min-width: 0;
        overflow-x: hidden;
      }

      .main img,
      .main svg,
      .main canvas {
        max-width: 100%;
      }

      #metadata_table,
      #upregulated_table,
      #downregulated_table,
      #mapping_results_table {
        display: block;
        width: 100%;
        max-width: 100%;
        overflow-x: auto;
        -webkit-overflow-scrolling: touch;
      }

      #metadata_table table,
      #upregulated_table table,
      #downregulated_table table,
      #mapping_results_table table {
        max-width: 100%;
      }

      .results-navigation-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 18px;
        width: 100%;
      }

      .results-navigation-grid .comparison-option-button {
        width: 100%;
        min-height: 150px;
        text-align: left;
        white-space: normal;
      }

      .results-navigation-grid .comparison-option-button h4,
      .results-navigation-grid .comparison-option-button p {
        white-space: normal;
      }

      .results-plot-selector {
        display: flex;
        flex-wrap: wrap;
        gap: 10px;
        margin: 18px 0;
      }

      @media (max-width: 700px) {
        .results-navigation-grid {
          grid-template-columns: 1fr;
        }
      }

      .mapping-results-scroll {
        width: 100%;
        max-width: 100%;
        max-height: 65vh;
        overflow: auto;
        -webkit-overflow-scrolling: touch;
      }

      #mapping_results_table {
        display: block;
        width: max-content;
        min-width: 100%;
        max-width: none;
        overflow: visible;
      }

      #mapping_results_table table {
        width: max-content;
        min-width: 100%;
        max-width: none;
      }

      @media (max-width: 768px) {
        .main {
          padding-left: 12px;
          padding-right: 12px;
        }

        .hero {
          padding: 24px 18px;
        }

        .workflow-layout {
          grid-template-columns: minmax(0, 1fr);
        }
      }

      .bslib-page-sidebar,
      .bslib-sidebar-layout,
      .main {
        background-color: #F7F0E8;
        color: #242020;
      }

      .bslib-page-sidebar > .navbar > .container-fluid,
      .navbar,
      .navbar-default {
        background-color: #E8D5D2;
      }

      .hero {
        background: linear-gradient(120deg, #E7C5C5, #F7F0E8);
        color: #242020;
        padding: 36px;
        border: 1px solid #E7C5C5;
        border-radius: 22px;
        margin-bottom: 28px;
      }

      .hero h1 {
        color: #242020;
        font-size: clamp(30px, 4vw, 48px);
        font-weight: 850;
        letter-spacing: -1.5px;
      }

      .hero p {
        color: #242020;
        line-height: 1.8;
        max-width: 720px;
      }

      .eyebrow {
        color: #59494A;
        font-size: 11px;
        font-weight: 800;
        letter-spacing: 2px;
      }

      .section-title {
        font-weight: 800;
        margin-top: 34px;
        margin-bottom: 8px;
        color: #242020;
      }

      .info-card,
      .stage-card,
      .candidate,
      .workflow-card {
        background: #FFFDFC;
        color: #242020;
        border: 1px solid #E8D5D2;
        border-radius: 16px;
        padding: 20px;
        height: 100%;
        transition: transform .2s ease, box-shadow .2s ease;
      }

      .info-card:hover,
      .stage-card:hover,
      .candidate:hover,
      .workflow-card:hover {
        transform: translateY(-4px);
        box-shadow: 0 10px 25px rgba(89, 73, 74, .12);
      }

      .card-symbol {
        font-size: 26px;
        color: #B98791;
        margin-bottom: 12px;
      }

      .info-card h4,
      .stage-card h4,
      .workflow-card h4,
      .candidate strong {
        font-weight: 800;
        color: #242020;
      }

      .info-card p,
      .stage-card p,
      .workflow-card p,
      .candidate p {
        color: #393232;
        font-size: 13px;
        line-height: 1.7;
      }

      .stage-number {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        width: 34px;
        height: 34px;
        border-radius: 50%;
        background: #DCE8DC;
        color: #242020;
        font-weight: 800;
        margin-bottom: 12px;
      }

      .focus-panel {
        background: #F0E3DF;
        color: #242020;
        border-left: 5px solid #B98791;
        border-radius: 14px;
        padding: 25px;
        margin-top: 12px;
      }

      .focus-panel h3,
      .focus-panel h4 {
        color: #242020;
        font-weight: 800;
      }

      .focus-panel p {
        color: #393232;
        line-height: 1.8;
      }

      .footer-band {
        background: #E8D5D2;
        color: #242020;
        padding: 18px 22px;
        border-radius: 14px;
        margin-top: 28px;
        font-size: 13px;
      }

      .small-label {
        font-size: 10px;
        letter-spacing: 1.5px;
        font-weight: 800;
        color: #59494A;
      }

      #metadata_table table,
      table.table {
        border-collapse: collapse;
        border: 1px solid #E8D5D2 !important;
        background-color: #FFFDFC !important;
        color: #242020 !important;
      }

      #metadata_table th,
      table.table th {
        background-color: #E8D5D2 !important;
        color: #242020 !important;
        border: 1px solid #D6BDBB !important;
        padding: 12px;
      }

      #metadata_table td,
      table.table td {
        background-color: #FFFDFC !important;
        color: #242020 !important;
        border: 1px solid #E8D5D2 !important;
        padding: 12px;
      }

      .table-striped > tbody > tr:nth-of-type(odd) > * {
        background-color: #F7F0E8 !important;
        color: #242020 !important;
      }

      .btn-outline-success {
        color: #59494A;
        border-color: #B98791;
      }

      .btn-outline-success:hover {
        background-color: #B98791;
        color: #242020;
        border-color: #B98791;
      }

      .btn-light {
        background-color: #FFFDFC;
        color: #242020;
        border-color: #E8D5D2;
      }

      .sidebar .nav-link {
        color: #242020;
      }

      .sidebar .nav-link.active {
        background-color: #C28F99;
        color: #242020;
      }

      .comparison-option-button {
        display: block;
        width: 100%;
        text-align: left;
        cursor: pointer;
        color: #242020;
        border: 1px solid #E8D5D2;
        background: #FFFDFC;
        white-space: normal;
      }

      .comparison-option-button:hover,
      .comparison-option-button:focus-visible {
        background: #F0E3DF;
        color: #242020;
        border-color: #B98791;
        box-shadow: 0 8px 20px rgba(89, 73, 74, .12);
        transform: translateY(-3px);
      }

      .comparison-option-button h4 {
        font-weight: 800;
        color: #242020;
        margin-top: 4px;
      }

      .comparison-option-button p {
        color: #393232;
        font-size: 13px;
        line-height: 1.7;
        margin-bottom: 0;
      }

      .comparison-summary {
        background: #F0E3DF;
        border-left: 4px solid #B98791;
        border-radius: 10px;
        padding: 12px 16px;
        margin-top: 12px;
        color: #242020;
        overflow-wrap: anywhere;
      }

      .comparison-summary p {
        margin-bottom: 6px;
        line-height: 1.6;
      }

      .comparison-summary p:last-child {
        margin-bottom: 0;
      }

      .comparison-table-bottom-scroll {
        position: sticky;
        bottom: 0;
        z-index: 5;
        width: 100%;
        max-width: 100%;
        height: 18px;
        overflow-x: auto;
        overflow-y: hidden;
        background: #fff;
        border-top: 1px solid #e7e1df;
        -webkit-overflow-scrolling: touch;
      }

      .comparison-table-bottom-scroll-inner {
        height: 1px;
      }

      .comparison-table-scroll {
        width: 100%;
        max-width: 100%;
        max-height: 70vh;
        overflow-x: hidden;
        overflow-y: auto;
        -webkit-overflow-scrolling: touch;
      }

      #comparison_degs_table {
        display: block;
        width: 100%;
        max-width: 100%;
        overflow: auto;
      }

      #comparison_degs_table table {
        min-width: 1050px;
        width: max-content;
        max-width: none;
      }

      .workflow-intro {
        color: #59494A;
        max-width: 850px;
        line-height: 1.7;
        margin-bottom: 22px;
      }

      .workflow-layout {
        display: grid;
        grid-template-columns: minmax(0, 1.7fr) minmax(230px, 0.85fr);
        gap: 22px;
        align-items: start;
        max-width: 1250px;
        margin: 0 auto;
      }

      .workflow-diagram {
        position: relative;
        width: 100%;
        min-width: 0;
        aspect-ratio: 1000 / 620;
        min-height: 400px;
        overflow: visible;
      }

      .workflow-connectors {
        position: absolute;
        inset: 0;
        width: 100%;
        height: 100%;
        z-index: 1;
        pointer-events: none;
        overflow: visible;
      }

      .workflow-node {
        position: absolute;
        z-index: 2;
        width: 30%;
        min-height: 58px;
        display: flex;
        align-items: stretch;
      }

      .workflow-node .btn {
        width: 100%;
        min-height: 58px;
        white-space: normal;
        border-radius: 13px;
        border: 1px solid #B98791;
        background: #D77ACB;
        color: #242020;
        font-weight: 750;
        font-size: clamp(11px, 1.3vw, 16px);
        line-height: 1.35;
        padding: 10px 12px;
        box-shadow: 0 2px 5px rgba(89, 73, 74, .08);
        transition:
          background-color .18s ease,
          border-color .18s ease,
          transform .18s ease,
          box-shadow .18s ease;
      }

      .workflow-node .btn:hover,
      .workflow-node .btn:focus-visible {
        background: #DCE8DC;
        color: #242020;
        border-color: #91A991;
        transform: translateY(-2px);
        box-shadow: 0 7px 15px rgba(89, 73, 74, .16);
      }

      .workflow-node .btn.workflow-node-selected {
        background: #D77ACB;
        color: #242020;
        border: 2px solid #7D5A83;
        box-shadow: 0 0 0 2px rgba(125, 90, 131, .12);
      }

      .workflow-node .btn.workflow-node-selected:hover,
      .workflow-node .btn.workflow-node-selected:focus-visible {
        background: #DCE8DC;
        color: #242020;
        border-color: #91A991;
      }

      .workflow-prep {
        left: 35%;
        top: 2%;
      }

      .workflow-deseq {
        left: 35%;
        top: 19%;
      }

      .workflow-significant {
        left: 5%;
        top: 35.5%;
      }

      .workflow-go {
        left: 65%;
        top: 35.5%;
      }

      .workflow-human {
        left: 35%;
        top: 50%;
      }

      .workflow-enrichment {
        left: 35%;
        top: 67%;
      }

      .workflow-visualization {
        left: 35%;
        top: 84%;
      }

      .workflow-detail {
        background: #FFFDFC;
        border: 1px solid #E8D5D2;
        border-left: 5px solid #B98791;
        border-radius: 16px;
        padding: 22px;
        min-width: 0;
        min-height: 170px;
        overflow-wrap: anywhere;
      }

      .workflow-detail-label {
        color: #8F568B;
        font-size: 10px;
        font-weight: 800;
        letter-spacing: 1.6px;
        text-transform: uppercase;
        margin-bottom: 10px;
      }

      .workflow-detail h3 {
        margin-top: 0;
        margin-bottom: 12px;
        color: #242020;
        font-size: 21px;
        font-weight: 800;
      }

      .workflow-detail p {
        color: #393232;
        line-height: 1.75;
        font-size: 14px;
        margin-bottom: 12px;
      }

      .workflow-detail .workflow-why {
        padding-top: 12px;
        margin-top: 12px;
        border-top: 1px solid #E8D5D2;
      }

      .workflow-detail .workflow-why strong {
        color: #59494A;
      }

      .workflow-hint {
        font-size: 12px;
        color: #71696A;
        line-height: 1.6;
        margin-top: 14px;
      }

      @media (max-width: 900px) {
        .workflow-layout {
          grid-template-columns: minmax(0, 1fr);
        }

        .workflow-diagram {
          min-height: 0;
          aspect-ratio: 1000 / 620;
        }

        .workflow-node .btn {
          min-height: 48px;
          padding: 8px;
          font-size: 12px;
        }

        .workflow-detail {
          min-height: 0;
        }
      }

      @media (max-width: 550px) {
        .workflow-diagram {
          display: flex;
          flex-direction: column;
          gap: 12px;
          aspect-ratio: auto;
          min-height: 0;
        }

        .workflow-connectors {
          display: none;
        }

        .workflow-node {
          position: relative;
          left: auto;
          top: auto;
          width: 100%;
          min-height: 0;
        }

        .workflow-node .btn {
          min-height: 48px;
          font-size: 14px;
        }

        .workflow-node:not(:last-child)::after {
          content: '↓';
          position: absolute;
          left: 50%;
          top: 100%;
          transform: translateX(-50%);
          height: 12px;
          color: #8B777B;
          font-size: 16px;
          line-height: 12px;
          z-index: 2;
        }
      }

      @media (prefers-reduced-motion: reduce) {
        .info-card,
        .stage-card,
        .candidate,
        .workflow-card,
        .workflow-node .btn {
          transition: none;
        }
      }

      .de-column {
        min-width: 0;
      }

      .de-column .card-symbol {
        font-size: 26px;
        color: #B98791;
        margin-bottom: 12px;
      }

      .de-column h4 {
        font-weight: 800;
        color: #242020;
      }

      .de-column p {
        color: #393232;
        font-size: 13px;
        line-height: 1.7;
      }

      .de-column .btn-primary {
        background-color: #C28F99;
        border-color: #C28F99;
        color: #242020;
      }

      .de-column .btn-primary:hover,
      .de-column .btn-primary:focus {
        background-color: #B98791;
        border-color: #B98791;
        color: #242020;
      }

      .gene-mapping-section {
        background: #FFFDFC;
        color: #242020;
        border: 1px solid #E8D5D2;
        border-radius: 16px;
        padding: 22px;
        min-height: 190px;
        height: 100%;
      }

      .gene-mapping-section h3 {
        color: #242020;
        font-weight: 800;
        margin-top: 0;
        margin-bottom: 12px;
      }

      .gene-mapping-section p {
        color: #393232;
        font-size: 13px;
        line-height: 1.7;
      }

      .gene-mapping-section .card-symbol {
        font-size: 26px;
        color: #B98791;
        margin-bottom: 12px;
      }

      .mapping-result-stat-card {
        padding: 12px 16px;
        min-height: 0;
        height: auto;
      }

      .mapping-result-stat-card h4 {
        margin: 0 0 4px;
        font-size: 14px;
        font-weight: 750;
      }

      .mapping-result-stat-card .gene-mapping-number {
        font-size: 25px;
        line-height: 1.15;
        margin: 0;
      }

      .gene-mapping-button {
        width: 100%;
        text-align: left;
        background: #FFFDFC;
        color: #242020;
        border: 1px solid #E8D5D2;
        border-radius: 16px;
        padding: 22px;
        min-height: 180px;
        transition: transform .2s ease, box-shadow .2s ease;
      }

      .gene-mapping-button:hover,
      .gene-mapping-button:focus {
        background: #F0E3DF;
        color: #242020;
        border-color: #B98791;
        transform: translateY(-3px);
        box-shadow: 0 8px 20px rgba(89, 73, 74, .12);
      }

      .gene-mapping-button h3 {
        font-weight: 800;
        margin: 10px 0;
        color: #242020;
      }

      .gene-mapping-button p {
        color: #393232;
        font-size: 13px;
        line-height: 1.7;
        margin-bottom: 0;
      }

      .gene-mapping-button .card-symbol {
        font-size: 26px;
        color: #B98791;
      }

      .gene-mapping-number {
        font-size: 34px;
        font-weight: 850;
        color: #8F568B;
        margin: 5px 0;
      }

      .gene-mapping-flow {
        display: flex;
        align-items: center;
        justify-content: center;
        gap: 18px;
        flex-wrap: wrap;
        margin-top: 28px;
      }

      .gene-mapping-flow-step {
        background: #F0E3DF;
        border: 1px solid #E8D5D2;
        border-radius: 14px;
        padding: 22px;
        min-width: 150px;
        text-align: center;
        font-weight: 750;
      }

      .gene-mapping-flow-step strong {
        display: block;
        font-size: 30px;
        color: #8F568B;
        margin-bottom: 6px;
      }

      .gene-mapping-flow-arrow {
        font-size: 26px;
        font-weight: 800;
        color: #8F568B;
      }

      @media (max-width: 550px) {
        .gene-mapping-flow {
          flex-direction: column;
        }

        .gene-mapping-flow-arrow {
          transform: rotate(90deg);
        }
      }

      /* FUNCTIONAL ENRICHMENT */

      .enrichment-panel {
        background: #FFFDFC;
        border: 1px solid #E8D5D2;
        border-radius: 16px;
        padding: 22px;
        min-width: 0;
        height: 100%;
      }

      .enrichment-panel h3 {
        font-size: 19px;
        font-weight: 800;
        color: #242020;
        margin-top: 0;
        margin-bottom: 14px;
      }

      .enrichment-panel .form-group {
        margin-bottom: 14px;
      }

      .enrichment-plot {
        width: 100%;
        min-width: 0;
        overflow-x: auto;
      }

      .enrichment-note {
        color: #59494A;
        font-size: 12px;
        line-height: 1.6;
      }

      .enrichment-layout {
        display: grid;
        grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: 20px;
        align-items: stretch;
        width: 100%;
      }

      .enrichment-option-button {
        width: 100%;
        height: 100%;
        min-height: 190px;
        padding: 25px;
        text-align: left;
        white-space: normal;
        background: #FFFDFC;
        color: #242020;
        border: 1px solid #E8D5D2;
        border-radius: 16px;
        transition: transform .2s ease, box-shadow .2s ease,
                    background-color .2s ease;
      }

      .enrichment-option-button:hover,
      .enrichment-option-button:focus-visible {
        background: #F0E3DF;
        color: #242020;
        border-color: #B98791;
        transform: translateY(-4px);
        box-shadow: 0 10px 25px rgba(89, 73, 74, .12);
      }

      .enrichment-option-button .card-symbol {
        font-size: 30px;
        color: #B98791;
        margin-bottom: 16px;
      }

      .enrichment-option-button h3 {
        font-size: 20px;
        font-weight: 800;
        color: #242020;
        margin: 0 0 12px;
      }

      .enrichment-option-button p {
        color: #59494A;
        font-size: 13px;
        line-height: 1.7;
        margin: 0;
      }

      .enrichment-detail {
        width: 100%;
        max-width: 1200px;
        margin: 0 auto;
      }

      .enrichment-detail .enrichment-panel {
        width: 100%;
      }

      @media (max-width: 1050px) {
        .enrichment-layout {
          grid-template-columns: repeat(2, minmax(0, 1fr));
        }
      }

      @media (max-width: 600px) {
        .enrichment-layout {
          grid-template-columns: minmax(0, 1fr);
        }

        .enrichment-option-button {
          min-height: 150px;
        }
      }

      /* ReForm animated axolotl cursor */
      body.reform-axolotl-cursor-enabled,
      body.reform-axolotl-cursor-enabled * {
        cursor: none !important;
      }

      body.reform-axolotl-cursor-enabled input,
      body.reform-axolotl-cursor-enabled textarea,
      body.reform-axolotl-cursor-enabled [contenteditable='true'] {
        cursor: text !important;
      }

      #reform-axolotl-cursor {
        position: fixed;
        left: 0;
        top: 0;
        width: 62px;
        height: 51px;
        z-index: 999999;
        pointer-events: none;
        display: none;
        transform: translate3d(0, 0, 0);
        will-change: transform;
      }

      #reform-axolotl-cursor img {
        display: block;
        width: 62px;
        height: 51px;
        object-fit: contain;
        border: 0 !important;
        border-radius: 0 !important;
        box-shadow: none !important;
        transform: none !important;
        user-select: none;
        -webkit-user-drag: none;
      }

      .reform-axolotl-bubble {
        position: fixed;
        z-index: 999998;
        width: var(--bubble-size, 10px);
        height: var(--bubble-size, 10px);
        border: 1.6px solid rgba(174, 196, 232, .96);
        border-top-color: rgba(255, 255, 255, 1);
        border-right-color: rgba(242, 190, 222, .98);
        border-bottom-color: rgba(255, 255, 255, .98);
        border-radius: 50%;
        background: radial-gradient(circle at 30% 25%,
          rgba(255,255,255,.92) 0%,
          rgba(231,241,255,.38) 28%,
          rgba(203,222,248,.20) 58%,
          rgba(242,203,229,.16) 100%);
        box-shadow:
          inset 1px 1px 2px rgba(255,255,255,.95),
          inset -1px -1px 2px rgba(187,190,224,.22),
          0 0 7px rgba(207,195,235,.42),
          0 0 2px rgba(255,255,255,.9);
        opacity: .95;
        pointer-events: none;
        animation: reformBubbleFloat 1250ms ease-out forwards;
      }

      @keyframes reformBubbleFloat {
        0% { opacity: .95; transform: translate3d(0, 0, 0) scale(.65); }
        35% { opacity: .88; }
        70% { opacity: .58; }
        100% { opacity: 0; transform: translate3d(var(--bubble-drift, 8px), -34px, 0) scale(1.2); }
      }

      @media (hover: none), (pointer: coarse) {
        #reform-axolotl-cursor,
        .reform-axolotl-bubble {
          display: none !important;
        }
      }
    ")),
    tags$script(HTML("
      (function () {
        function initialiseReformAxolotlCursor() {
          if (window.matchMedia && window.matchMedia('(hover: none), (pointer: coarse)').matches) return;
          if (document.getElementById('reform-axolotl-cursor')) return;

          var assets = {
            baby: 'axolotl_cursor_transparent.png',
            sleeping: 'Sleeping Kawaii Pink Axolotl with Zzz.png',
            waving: 'Kawaii Axolotl Waving Animation.gif'
          };

          var cursor = document.createElement('div');
          cursor.id = 'reform-axolotl-cursor';
          document.body.classList.add('reform-axolotl-cursor-enabled');
          var character = document.createElement('img');
          character.alt = '';
          character.draggable = false;
          character.src = assets.baby;
          cursor.appendChild(character);
          document.body.appendChild(cursor);

          var currentState = 'baby';
          var idleTimer = null;
          var lastBubbleAt = 0;
          var lastX = 0;
          var lastY = 0;
          var hasMoved = false;

          function setState(state) {
            if (currentState === state) return;
            currentState = state;
            character.src = assets[state];
          }

          function addBubble(x, y) {
            var bubble = document.createElement('span');
            var size = 7 + Math.random() * 8;
            bubble.className = 'reform-axolotl-bubble';
            bubble.style.left = (x + (Math.random() * 28 - 14)) + 'px';
            bubble.style.top = (y + 22 + Math.random() * 24) + 'px';
            bubble.style.setProperty('--bubble-size', size + 'px');
            bubble.style.setProperty('--bubble-drift', (Math.random() * 18 - 9) + 'px');
            document.body.appendChild(bubble);
            window.setTimeout(function () { bubble.remove(); }, 1050);
          }

          document.addEventListener('mousemove', function (event) {
            var dx = event.clientX - lastX;
            var dy = event.clientY - lastY;
            var distance = Math.sqrt(dx * dx + dy * dy);
            lastX = event.clientX;
            lastY = event.clientY;

            cursor.style.display = 'block';
            cursor.style.transform = 'translate3d(' + event.clientX + 'px, ' + event.clientY + 'px, 0)';

            if (distance > 1.5) {
              hasMoved = true;
              setState('waving');
              var now = Date.now();
              if (now - lastBubbleAt > 75) {
                addBubble(event.clientX + 18, event.clientY + 18);
                lastBubbleAt = now;
              }
              if (idleTimer) window.clearTimeout(idleTimer);
              idleTimer = window.setTimeout(function () {
                setState('sleeping');
              }, 1250);
            } else if (!hasMoved) {
              setState('baby');
            }
          }, { passive: true });

          document.addEventListener('mouseleave', function () {
            cursor.style.display = 'none';
          });

          document.addEventListener('mouseenter', function () {
            if (hasMoved) cursor.style.display = 'block';
          });

          document.addEventListener('mousedown', function () {
            if (currentState !== 'waving') setState('baby');
          });
        }

        if (document.readyState === 'loading') {
          document.addEventListener('DOMContentLoaded', initialiseReformAxolotlCursor);
        } else {
          initialiseReformAxolotlCursor();
        }
      })();
    "))
  ),

  uiOutput("home_content")
)


# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # Supplied standalone-analysis results. These are displayed directly; PCA and
  # sample correlations are NOT recalculated inside the Shiny app.
  pca_coordinates <- data.frame(
    Sample = c(
    "Intact1_Diaphysis_L183055_Track-224733",
    "Intact2_Diaphysis_L183057_Track-224735",
    "Intact3_Diaphysis_L183059_Track-224737",
    "Intact1_Epiphysis_L183054_Track-224732",
    "Intact2_Epiphysis_L183056_Track-224734",
    "Intact3_Epiphysis_L183058_Track-224736",
    "3dpa_Diaphysis_1_L183063_Track-224741",
    "3dpa_Diaphysis_2_L183064_Track-224742",
    "3dpa_Diaphysis_3_L183065_Track-224743",
    "5dpa_Diaphysis_1_L183069_Track-224603",
    "5dpa_Diaphysis_2_L183070_Track-224604",
    "5dpa_Diaphysis_3_L183071_Track-224605",
    "3dpa_Epiphysis_1_L183060_Track-224738",
    "3dpa_Epiphysis_2_L183061_Track-224739",
    "3dpa_Epiphysis_3_L183062_Track-224740",
    "5dpa_Epiphysis_1_L183066_Track-224744",
    "5dpa_Epiphysis_2_L183067_Track-224745",
    "5dpa_Epiphysis_3_L183068_Track-224746"
    ),
    PC1 = c(-84.2956, -93.5339, -101.99, -70.3531, -75.9809, -89.7532, 35.82848, 7.828733, 34.58973, 36.88921, 59.45602, 48.40662, 53.27107, 25.25894, 28.29276, 62.16395, 59.70276, 64.21878),
    PC2 = c(54.9105, -12.3178, -71.981, 98.75017, 20.93139, -76.0265, 60.3129, 12.32226, 9.430685, -22.4912, -16.1496, -6.29167, 34.31231, -4.89871, -20.9124, -14.6464, -14.2721, -30.9827),
    PC3 = c(-9.36127, -5.75847, 10.26299, -27.8149, -3.61632, 10.65437, 60.42118, 65.4367, 15.37675, -3.30182, 3.025266, 4.629348, -31.4154, -35.6714, -23.507, 1.322076, -17.7829, -12.8991),
    stringsAsFactors = FALSE
  )

  # Recover the experimental annotations from each sample's established label.
  pca_coordinates$Condition <- ifelse(grepl("^Intact", pca_coordinates$Sample), "Intact",
    sub("^([^_]+)_.*$", "\\1", pca_coordinates$Sample))
  pca_coordinates$Region <- ifelse(grepl("Diaphysis", pca_coordinates$Sample), "Diaphysis", "Epiphysis")
  pca_coordinates$Replicate <- ifelse(
    pca_coordinates$Condition == "Intact",
    sub("^Intact([0-9]+)_.*$", "\\1", pca_coordinates$Sample),
    sub("^[^_]+_[^_]+_([0-9]+)_.*$", "\\1", pca_coordinates$Sample)
  )
  pca_coordinates$Group <- paste(pca_coordinates$Condition, pca_coordinates$Region, sep = "_")
  pca_coordinates$Group <- factor(pca_coordinates$Group, levels = c(
    "Intact_Diaphysis", "3dpa_Diaphysis", "5dpa_Diaphysis",
    "Intact_Epiphysis", "3dpa_Epiphysis", "5dpa_Epiphysis"
  ))

  correlation_values <- rbind(
    c(1, 0.973126, 0.946713, 0.983, 0.97885, 0.938322, 0.931715, 0.944854, 0.939339, 0.930341, 0.914655, 0.925342, 0.925262, 0.935081, 0.931667, 0.912561, 0.919073, 0.906691),
    c(0.973126, 1, 0.974324, 0.947421, 0.986533, 0.964473, 0.911901, 0.940597, 0.93681, 0.937383, 0.920125, 0.923863, 0.91119, 0.93683, 0.93704, 0.914963, 0.920773, 0.914077),
    c(0.946713, 0.974324, 1, 0.909871, 0.956243, 0.989073, 0.889405, 0.92458, 0.918135, 0.924153, 0.906358, 0.90844, 0.891361, 0.920606, 0.928547, 0.904645, 0.906126, 0.903755),
    c(0.983, 0.947421, 0.909871, 1, 0.965628, 0.906343, 0.930468, 0.927353, 0.929135, 0.914533, 0.904625, 0.915515, 0.928854, 0.923307, 0.917243, 0.904006, 0.908788, 0.893866),
    c(0.97885, 0.986533, 0.956243, 0.965628, 1, 0.955848, 0.933703, 0.949468, 0.94726, 0.942489, 0.928379, 0.935674, 0.92989, 0.943003, 0.941784, 0.925507, 0.930881, 0.921811),
    c(0.938322, 0.964473, 0.989073, 0.906343, 0.955848, 1, 0.894362, 0.923736, 0.920279, 0.925664, 0.910303, 0.913902, 0.898845, 0.9221, 0.932611, 0.912632, 0.910817, 0.910308),
    c(0.931715, 0.911901, 0.889405, 0.930468, 0.933703, 0.894362, 1, 0.962057, 0.966157, 0.947675, 0.954131, 0.954236, 0.962858, 0.94419, 0.946418, 0.95523, 0.950514, 0.942617),
    c(0.944854, 0.940597, 0.92458, 0.927353, 0.949468, 0.923736, 0.962057, 1, 0.974235, 0.963086, 0.956198, 0.960402, 0.943694, 0.953089, 0.956874, 0.956539, 0.955047, 0.950128),
    c(0.939339, 0.93681, 0.918135, 0.929135, 0.94726, 0.920279, 0.966157, 0.974235, 1, 0.977416, 0.973131, 0.973039, 0.968995, 0.972059, 0.973418, 0.972303, 0.975434, 0.968214),
    c(0.930341, 0.937383, 0.924153, 0.914533, 0.942489, 0.925664, 0.947675, 0.963086, 0.977416, 1, 0.97777, 0.978985, 0.959578, 0.96686, 0.972873, 0.972995, 0.980815, 0.975087),
    c(0.914655, 0.920125, 0.906358, 0.904625, 0.928379, 0.910303, 0.954131, 0.956198, 0.973131, 0.97777, 1, 0.977246, 0.963695, 0.96229, 0.965324, 0.977042, 0.978556, 0.9801),
    c(0.925342, 0.923863, 0.90844, 0.915515, 0.935674, 0.913902, 0.954236, 0.960402, 0.973039, 0.978985, 0.977246, 1, 0.960878, 0.961674, 0.964396, 0.973847, 0.97873, 0.973537),
    c(0.925262, 0.91119, 0.891361, 0.928854, 0.92989, 0.898845, 0.962858, 0.943694, 0.968995, 0.959578, 0.963695, 0.960878, 1, 0.964089, 0.968139, 0.969404, 0.968542, 0.959023),
    c(0.935081, 0.93683, 0.920606, 0.923307, 0.943003, 0.9221, 0.94419, 0.953089, 0.972059, 0.96686, 0.96229, 0.961674, 0.964089, 1, 0.971795, 0.963289, 0.967979, 0.965469),
    c(0.931667, 0.93704, 0.928547, 0.917243, 0.941784, 0.932611, 0.946418, 0.956874, 0.973418, 0.972873, 0.965324, 0.964396, 0.968139, 0.971795, 1, 0.969971, 0.973138, 0.966654),
    c(0.912561, 0.914963, 0.904645, 0.904006, 0.925507, 0.912632, 0.95523, 0.956539, 0.972303, 0.972995, 0.977042, 0.973847, 0.969404, 0.963289, 0.969971, 1, 0.977451, 0.975386),
    c(0.919073, 0.920773, 0.906126, 0.908788, 0.930881, 0.910817, 0.950514, 0.955047, 0.975434, 0.980815, 0.978556, 0.97873, 0.968542, 0.967979, 0.973138, 0.977451, 1, 0.980374),
    c(0.906691, 0.914077, 0.903755, 0.893866, 0.921811, 0.910308, 0.942617, 0.950128, 0.968214, 0.975087, 0.9801, 0.973537, 0.959023, 0.965469, 0.966654, 0.975386, 0.980374, 1)
  )
  sample_correlation_matrix <- correlation_values
  rownames(sample_correlation_matrix) <- pca_coordinates$Sample
  colnames(sample_correlation_matrix) <- pca_coordinates$Sample

  if (nrow(pca_coordinates) != 18L || anyNA(pca_coordinates) ||
      !isTRUE(all.equal(sample_correlation_matrix, t(sample_correlation_matrix), tolerance = 1e-8)) ||
      any(abs(diag(sample_correlation_matrix) - 1) > 1e-8)) {
    stop("The supplied PCA coordinates or sample-correlation matrix failed validation.")
  }

  # LOAD SAMPLE METADATA

  metadata <- read.csv(
    "Metadata1.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  # LOAD THE NCBI GENE-LEVEL COUNT MATRIX FOR THE EXISTING VOLCANO-PLOT ANALYSIS
  # Sample columns are matched to metadata by their sample labels, not row order.

  counts_file <- "GSE334432_Processed-Data-Gene-Counts.txt.gz"
  sample_analysis <- NULL
  sample_analysis_error <- NULL

  if (!file.exists(counts_file)) {

    sample_analysis_error <- paste(
      "The NCBI gene-count file was not found:",
      counts_file,
      "Place it in the same folder as app.R."
    )

  } else if (!requireNamespace("DESeq2", quietly = TRUE)) {

    sample_analysis_error <- paste(
      "The DESeq2 R package is required for the volcano-plot re-analysis.",
      "Install it with:",
      "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager');",
      "BiocManager::install('DESeq2')"
    )

  } else if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
             !requireNamespace("ggplot2", quietly = TRUE)) {

    sample_analysis_error <- paste(
      "Required plotting/assay packages are missing. Install them with:",
      "BiocManager::install('SummarizedExperiment');",
      "install.packages('ggplot2')"
    )

  } else {

    analysis_attempt <- tryCatch({

      if (!is.data.frame(metadata) || nrow(metadata) != 18) {
        stop("Metadata1.csv must contain exactly 18 sample rows.")
      }

      # The supplied metadata has a trailing space in the region column name.
      # Identify the column after trimming names without changing the shared metadata.
      metadata_column_names <- trimws(names(metadata))
      region_column_index <- which(metadata_column_names == "region")
      if (length(region_column_index) != 1) {
        stop("Metadata1.csv must contain exactly one 'region' column.")
      }
      region_values <- trimws(as.character(metadata[[region_column_index]]))

      required_metadata_columns <- c("sample_id", "condition", "replicate")
      missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
      if (length(missing_metadata_columns) > 0) {
        stop(paste(
          "Metadata1.csv is missing required column(s):",
          paste(missing_metadata_columns, collapse = ", ")
        ))
      }

      metadata_for_deseq <- data.frame(
        sample_id = trimws(as.character(metadata$sample_id)),
        condition = trimws(as.character(metadata$condition)),
        region = region_values,
        replicate = trimws(as.character(metadata$replicate)),
        stringsAsFactors = FALSE
      )

      if (anyNA(metadata_for_deseq) ||
          any(!nzchar(metadata_for_deseq$sample_id)) ||
          anyDuplicated(metadata_for_deseq$sample_id)) {
        stop("Metadata1.csv contains missing or duplicated sample identifiers/annotations.")
      }

      allowed_conditions <- c("Intact", "3dpa", "5dpa")
      allowed_regions <- c("Diaphysis", "Epiphysis")
      if (any(!metadata_for_deseq$condition %in% allowed_conditions) ||
          any(!metadata_for_deseq$region %in% allowed_regions)) {
        stop("Metadata1.csv contains an unexpected condition or tissue-region label.")
      }

      # Read the tab-delimited NCBI file. Use an actual tab delimiter (not the
      # two-character string backslash+t).
      count_table <- utils::read.delim(
        counts_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        quote = "",
        comment.char = "",
        row.names = NULL
      )

      if (ncol(count_table) != 19) {
        stop(paste(
          "Expected one gene-ID column and 18 sample columns in the NCBI count file;",
          "found", ncol(count_table) - 1, "sample columns."
        ))
      }

      observed_sample_labels <- sub(
        "_L[0-9]+_Track-.*$", "", names(count_table)[-1]
      )

      expected_sample_labels <- vapply(seq_len(nrow(metadata_for_deseq)), function(i) {
        condition_value <- metadata_for_deseq$condition[[i]]
        region_value <- metadata_for_deseq$region[[i]]
        replicate_value <- metadata_for_deseq$replicate[[i]]

        if (identical(condition_value, "Intact")) {
          paste0("Intact", replicate_value, "_", region_value)
        } else {
          paste0(condition_value, "_", region_value, "_", replicate_value)
        }
      }, character(1))

      if (anyDuplicated(observed_sample_labels)) {
        stop("The NCBI count matrix contains duplicated sample labels.")
      }

      sample_column_order <- match(expected_sample_labels, observed_sample_labels)
      if (anyNA(sample_column_order)) {
        missing_samples <- expected_sample_labels[is.na(sample_column_order)]
        stop(paste(
          "The following metadata samples were not found in the NCBI count matrix:",
          paste(missing_samples, collapse = ", ")
        ))
      }

      # Reorder count columns to exactly match Metadata1.csv.
      count_table <- count_table[, c(1, sample_column_order + 1), drop = FALSE]

      gene_ids <- trimws(as.character(count_table[[1]]))
      if (anyNA(gene_ids) || any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
        stop("The NCBI count matrix contains missing, blank, or duplicated gene IDs.")
      }

      numeric_counts <- lapply(count_table[-1], function(x) {
        suppressWarnings(as.numeric(as.character(x)))
      })

      if (any(vapply(numeric_counts, function(x) anyNA(x), logical(1)))) {
        stop("One or more NCBI sample columns contain non-numeric or missing count values.")
      }

      count_matrix <- do.call(cbind, numeric_counts)
      storage.mode(count_matrix) <- "double"

      if (any(!is.finite(count_matrix)) || any(count_matrix < 0) ||
          any(count_matrix != floor(count_matrix))) {
        stop("The NCBI gene-count matrix must contain finite, non-negative integer counts.")
      }

      rownames(count_matrix) <- gene_ids
      colnames(count_matrix) <- metadata_for_deseq$sample_id
      storage.mode(count_matrix) <- "integer"

      metadata_for_deseq$group <- factor(
        paste(metadata_for_deseq$condition, metadata_for_deseq$region, sep = "_"),
        levels = c(
          "Intact_Diaphysis", "3dpa_Diaphysis", "5dpa_Diaphysis",
          "Intact_Epiphysis", "3dpa_Epiphysis", "5dpa_Epiphysis"
        )
      )
      rownames(metadata_for_deseq) <- metadata_for_deseq$sample_id

      if (anyNA(metadata_for_deseq$group)) {
        stop("One or more samples could not be assigned to the six experimental groups.")
      }

      group_sizes <- table(metadata_for_deseq$group)
      if (length(group_sizes) != 6 || any(group_sizes != 3)) {
        stop("Expected exactly three replicates in each of the six experimental groups.")
      }

      # Retain genes with at least 10 total reads across the 18 samples.
      keep_genes <- rowSums(count_matrix) >= 10
      count_matrix <- count_matrix[keep_genes, , drop = FALSE]
      if (nrow(count_matrix) < 2) {
        stop("Fewer than two genes remain after the minimum-count filter.")
      }

      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = metadata_for_deseq,
        design = ~ group
      )
      dds <- DESeq2::DESeq(dds, quiet = TRUE)

      vsd <- DESeq2::vst(
        dds,
        blind = TRUE,
        nsub = min(1000L, nrow(dds))
      )

      list(ok = TRUE, value = list(
        dds = dds,
        vsd = vsd,
        metadata = metadata_for_deseq,
        gene_count = nrow(count_matrix),
        sample_count = ncol(count_matrix)
      ))

    }, error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    })

    if (isTRUE(analysis_attempt$ok)) {
      sample_analysis <- analysis_attempt$value
      sample_analysis_error <- NULL
    } else {
      sample_analysis <- NULL
      sample_analysis_error <- analysis_attempt$message
      message("DESeq2 sample-analysis preparation failed: ", sample_analysis_error)
    }
  }


  # LOAD SIGNIFICANT DESeq2 RESULTS FOR THE 3 dpa vs 5 dpa COMPARISON

  significant_degs_file <- "Significant_DEGs.csv"
  significant_degs <- read_csv_safe(significant_degs_file)

  comparison_required_columns <- c(
    "gene_id",
    "comparison",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "regulation"
  )

  comparison_data_error <- NULL

  if (is.null(significant_degs)) {
    comparison_data_error <- paste(
      "Significant DESeq2 results file not found or could not be read:",
      significant_degs_file,
      "Place it in the same folder as app.R."
    )
  } else {
    missing_comparison_columns <- setdiff(
      comparison_required_columns,
      names(significant_degs)
    )

    if (length(missing_comparison_columns) > 0) {
      comparison_data_error <- paste(
        "The significant DEGs CSV is missing required column(s):",
        paste(missing_comparison_columns, collapse = ", ")
      )
    } else {
      significant_degs$comparison <- trimws(
        as.character(significant_degs$comparison)
      )

      significant_degs$gene_id <- as.character(significant_degs$gene_id)
      significant_degs$regulation <- as.character(significant_degs$regulation)

      numeric_comparison_columns <- c(
        "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
      )

      for (column_name in numeric_comparison_columns) {
        significant_degs[[column_name]] <- suppressWarnings(
          as.numeric(significant_degs[[column_name]])
        )
      }
    }
  }


  # LOAD UPREGULATED GENES

  upregulated_file <- "Top_Upregulated.csv"

  upregulated_genes <- if (file.exists(upregulated_file)) {
    read.csv(
      upregulated_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    NULL
  }

  if (!is.null(upregulated_genes)) {
    upregulated_genes <- head(upregulated_genes, 10)
  }


  # LOAD DOWNREGULATED GENES

  downregulated_file <- "Top_Downregulated.csv"

  downregulated_genes <- if (file.exists(downregulated_file)) {
    read.csv(
      downregulated_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    NULL
  }

  if (!is.null(downregulated_genes)) {
    downregulated_genes <- head(downregulated_genes, 10)
  }


  # LOAD AXOLOTL CANDIDATE GENES AND HUMAN GENE COMPARISON RESULTS

  candidate_genes_file <- "Axolotl_39_Candidate_Genes.csv"
  human_mapping_file <- "Axolotl_Human_Gene_Comparison.csv"

  candidate_genes <- if (file.exists(candidate_genes_file)) {
    read.csv(
      candidate_genes_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    NULL
  }

  human_mapping <- if (file.exists(human_mapping_file)) {
    read.csv(
      human_mapping_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    NULL
  }

  mapping_data <- NULL
  mapping_data_error <- NULL

  if (is.null(candidate_genes)) {
    mapping_data_error <- paste(
      "Candidate-gene file not found:", candidate_genes_file,
      "Place it in the same folder as app.R."
    )
  } else if (is.null(human_mapping)) {
    mapping_data_error <- paste(
      "Human gene comparison file not found:", human_mapping_file,
      "Place it in the same folder as app.R."
    )
  } else {

    required_candidate_columns <- c(
      "candidate_number",
      "axolotl_annotation"
    )

    required_mapping_columns <- c(
      "candidate",
      "axolotl_gene_id",
      "description",
      "ncbi_gene_id",
      "human_gene_symbol",
      "human_gene_description",
      "search_status"
    )

    missing_candidate_columns <- setdiff(
      required_candidate_columns,
      names(candidate_genes)
    )

    missing_mapping_columns <- setdiff(
      required_mapping_columns,
      names(human_mapping)
    )

    if (length(missing_candidate_columns) > 0) {

      mapping_data_error <- paste(
        "The candidate-gene CSV is missing required column(s):",
        paste(missing_candidate_columns, collapse = ", ")
      )

    } else if (length(missing_mapping_columns) > 0) {

      mapping_data_error <- paste(
        "The human gene comparison CSV is missing required column(s):",
        paste(missing_mapping_columns, collapse = ", ")
      )

    } else {

      candidate_genes$axolotl_annotation <- trimws(
        as.character(candidate_genes$axolotl_annotation)
      )

      human_mapping$candidate <- trimws(
        as.character(human_mapping$candidate)
      )

      mapping_data <- merge(
        candidate_genes,
        human_mapping,
        by.x = "axolotl_annotation",
        by.y = "candidate",
        all.x = TRUE,
        sort = FALSE
      )

      mapping_data <- mapping_data[
        order(mapping_data$candidate_number),
      ]

      mapping_record_counts <- table(
        human_mapping$candidate
      )

      mapping_data$mapping_record_count <- as.integer(
        mapping_record_counts[
          mapping_data$axolotl_annotation
        ]
      )

      mapping_data$mapping_record_count[
        is.na(mapping_data$mapping_record_count)
      ] <- 0L

      mapping_data$mapping_match_type <- ifelse(
        mapping_data$mapping_record_count > 1,
        "Multiple human records",
        ifelse(
          mapping_data$mapping_record_count == 1,
          "Single human record",
          "No human record"
        )
      )

      mapping_data$mapping_status <- ifelse(
        !is.na(mapping_data$ncbi_gene_id) &
          nzchar(trimws(as.character(mapping_data$ncbi_gene_id))),
        "Human gene record found",
        "No human gene record found"
      )

      mapping_data <- mapping_data[, c(
        "candidate_number",
        "axolotl_annotation",
        "axolotl_gene_id",
        "description",
        "ncbi_gene_id",
        "human_gene_symbol",
        "human_gene_description",
        "search_status",
        "mapping_status",
        "mapping_record_count",
        "mapping_match_type"
      )]

      rownames(mapping_data) <- NULL
    }
  }


  # PCA and correlation are supplied standalone-analysis results. The volcano
  # plot below retains its existing DESeq2/count-matrix workflow.
  output$sample_analysis_status <- renderUI({
    div(
      class = "focus-panel",
      style = "margin: 12px 0;",
      tags$strong("PCA and sample-correlation results loaded"),
      p("The supplied 18-sample PCA coordinates and Pearson correlation matrix are displayed directly. The volcano plot remains a separate DESeq2 re-analysis.", style = "margin:6px 0 0;")
    )
  })


  # SUPPLIED PCA COORDINATES AND SAMPLE-CORRELATION HEATMAP
  # These precomputed results are displayed directly with interactive Plotly
  # outputs. They do not depend on the DESeq2/count-matrix volcano workflow.


  output$sample_pca_plot <- plotly::renderPlotly({
    p <- plotly::plot_ly(
      data = pca_coordinates,
      x = ~PC1, y = ~PC2,
      type = "scatter", mode = "markers",
      color = ~Group,
      colors = c("#B98278", "#8D9F91", "#C49A56", "#7C91B5", "#9A82AD", "#5E9C9A"),
      marker = list(size = 11, opacity = 0.88, line = list(width = 0.7, color = "white")),
      text = ~paste0(
        "Sample: ", Sample,
        "<br>Group: ", Group,
        "<br>Condition: ", Condition,
        "<br>Region: ", Region,
        "<br>Replicate: ", Replicate,
        "<br>PC1: ", round(PC1, 4),
        "<br>PC2: ", round(PC2, 4),
        "<br>PC3: ", round(PC3, 4)
      ),
      hoverinfo = "text"
    )
    p <- plotly::layout(
      p,
      title = list(text = "PCA of the 18 samples", x = 0.02),
      xaxis = list(title = "PC1", zeroline = FALSE),
      yaxis = list(title = "PC2", zeroline = FALSE),
      legend = list(title = list(text = "Experimental group")),
      hovermode = "closest"
    )
    p
  })

  output$sample_correlation_plot <- plotly::renderPlotly({
    sample_labels <- colnames(sample_correlation_matrix)
    plotly::plot_ly(
      x = sample_labels,
      y = sample_labels,
      z = unname(sample_correlation_matrix),
      type = "heatmap",
      colors = c("#3B6F9C", "#FFFFFF", "#B64B55"),
      zmin = -1, zmax = 1, zmid = 0,
      colorbar = list(title = "Pearson r"),
      hovertemplate = paste(
        "Sample 1: %{x}<br>",
        "Sample 2: %{y}<br>",
        "Pearson correlation: %{z:.4f}<extra></extra>"
      )
    ) |>
      plotly::layout(
        title = list(text = "Pairwise sample correlation", x = 0.02),
        xaxis = list(title = "", tickangle = -55, automargin = TRUE),
        yaxis = list(title = "", autorange = "reversed", automargin = TRUE),
        margin = list(b = 150, l = 190, r = 35, t = 60)
      )
  })

  # ==========================================================
  # COMPLETE VOLCANO PLOTS FOR ALL SIX DESeq2 COMPARISONS
  # Nonsignificant genes are retained in each contrast result.
  # ==========================================================

  output$deseq_volcano_plot <- plotly::renderPlotly({

    req(input$volcano_comparison)

    volcano_file_lookup <- c(
      "3dpa_Diaphysis_vs_Intact_Diaphysis" = "3dpa_Diaphysis_vs_Intact_Diaphysis_DESeq2_Results.csv",
      "3dpa_Epiphysis_vs_Intact_Epiphysis" = "3dpa_Epiphysis_vs_Intact_Epiphysis_DESeq2_Results.csv",
      "5dpa_Diaphysis_vs_Intact_Diaphysis" = "5dpa_Diaphysis_vs_Intact_Diaphysis_DESeq2_Results.csv",
      "5dpa_Epiphysis_vs_Intact_Epiphysis" = "5dpa_Epiphysis_vs_Intact_Epiphysis_DESeq2_Results.csv",
      "5dpa_vs_3dpa_Diaphysis" = "5dpa_vs_3dpa_Diaphysis_DESeq2_Results.csv",
      "5dpa_vs_3dpa_Epiphysis" = "5dpa_vs_3dpa_Epiphysis_DESeq2_Results.csv"
    )

    selected_file <- unname(volcano_file_lookup[input$volcano_comparison])

    shiny::validate(
      shiny::need(
        length(selected_file) == 1 && !is.na(selected_file),
        "Select one of the six supported DESeq2 comparisons."
      ),
      shiny::need(
        file.exists(selected_file),
        paste0(
          "The DESeq2 results file was not found: ", selected_file,
          ". Place all six comparison CSV files in the same folder as app.R."
        )
      )
    )

    volcano_data <- read_csv_safe(selected_file)

    required_volcano_columns <- c(
      "gene_id", "log2FoldChange", "padj", "regulation"
    )

    shiny::validate(
      shiny::need(
        !is.null(volcano_data),
        paste0("Could not read the selected DESeq2 results file: ", selected_file)
      ),
      shiny::need(
        all(required_volcano_columns %in% names(volcano_data)),
        paste0(
          "The selected DESeq2 results file is missing required columns: ",
          paste(
            setdiff(required_volcano_columns, names(volcano_data)),
            collapse = ", "
          )
        )
      )
    )

    volcano_data$log2FoldChange <- suppressWarnings(
      as.numeric(volcano_data$log2FoldChange)
    )
    volcano_data$padj <- suppressWarnings(as.numeric(volcano_data$padj))

    if ("neglog10_padj" %in% names(volcano_data)) {
      volcano_data$neglog10_padj <- suppressWarnings(
        as.numeric(volcano_data$neglog10_padj)
      )
    } else {
      volcano_data$neglog10_padj <- ifelse(
        is.na(volcano_data$padj),
        0,
        -log10(pmax(volcano_data$padj, .Machine$double.xmin))
      )
    }

    volcano_data$regulation <- trimws(as.character(volcano_data$regulation))
    volcano_data$regulation[
      is.na(volcano_data$regulation) |
        !volcano_data$regulation %in%
          c("Not significant", "Upregulated", "Downregulated")
    ] <- "Not significant"

    volcano_data <- volcano_data[
      is.finite(volcano_data$log2FoldChange) &
        is.finite(volcano_data$neglog10_padj),
      ,
      drop = FALSE
    ]

    shiny::validate(
      shiny::need(
        nrow(volcano_data) > 0,
        "No valid gene-level results are available for this comparison."
      )
    )

    volcano_data$regulation <- factor(
      volcano_data$regulation,
      levels = c("Not significant", "Upregulated", "Downregulated")
    )

    hover_columns <- intersect(
      c("baseMean", "pvalue", "padj"),
      names(volcano_data)
    )

    volcano_data$hover_text <- paste0(
      "<b>", volcano_data$gene_id, "</b>",
      "<br>Expression status: ", as.character(volcano_data$regulation),
      "<br>log2 fold change: ", signif(volcano_data$log2FoldChange, 4),
      "<br>-log10 adjusted p-value: ",
      signif(volcano_data$neglog10_padj, 4)
    )

    if ("baseMean" %in% hover_columns) {
      volcano_data$hover_text <- paste0(
        volcano_data$hover_text,
        "<br>Base mean: ", signif(as.numeric(volcano_data$baseMean), 4)
      )
    }
    if ("pvalue" %in% hover_columns) {
      volcano_data$hover_text <- paste0(
        volcano_data$hover_text,
        "<br>p-value: ", signif(as.numeric(volcano_data$pvalue), 4)
      )
    }
    if ("padj" %in% hover_columns) {
      volcano_data$hover_text <- paste0(
        volcano_data$hover_text,
        "<br>Adjusted p-value: ",
        ifelse(
          is.na(volcano_data$padj),
          "NA",
          format(
            signif(volcano_data$padj, 4),
            scientific = TRUE,
            trim = TRUE
          )
        )
      )
    }

    volcano_plot <- plotly::plot_ly(
      data = volcano_data,
      x = ~log2FoldChange,
      y = ~neglog10_padj,
      type = "scatter",
      mode = "markers",
      color = ~regulation,
      colors = c(
        "Not significant" = "#B8B0B0",
        "Upregulated" = "#C45A62",
        "Downregulated" = "#4C83A8"
      ),
      text = ~hover_text,
      hovertemplate = "%{text}<extra></extra>",
      marker = list(
        size = 6,
        opacity = 0.72,
        line = list(width = 0)
      )
    )

    volcano_plot <- plotly::layout(
      volcano_plot,
      title = list(
        text = gsub("_", " ", input$volcano_comparison),
        x = 0.02,
        xanchor = "left"
      ),
      xaxis = list(
        title = "log2 fold change",
        zeroline = FALSE
      ),
      yaxis = list(
        title = "-log10 adjusted p-value",
        rangemode = "tozero",
        zeroline = FALSE
      ),
      legend = list(title = list(text = "Expression status")),
      hovermode = "closest",
      shapes = list(
        list(
          type = "line",
          x0 = -1, x1 = -1, y0 = 0, y1 = 1,
          xref = "x", yref = "paper",
          line = list(color = "#6F6262", dash = "dash", width = 1)
        ),
        list(
          type = "line",
          x0 = 1, x1 = 1, y0 = 0, y1 = 1,
          xref = "x", yref = "paper",
          line = list(color = "#6F6262", dash = "dash", width = 1)
        ),
        list(
          type = "line",
          x0 = 0, x1 = 1, y0 = -log10(0.05), y1 = -log10(0.05),
          xref = "paper", yref = "y",
          line = list(color = "#6F6262", dash = "dash", width = 1)
        )
      ),
      margin = list(l = 75, r = 35, t = 65, b = 65),
      paper_bgcolor = "#FFFDFC",
      plot_bgcolor = "#FFFDFC",
      font = list(color = "#242020")
    )

    volcano_plot
  })

  # TABLES CORRESPONDING TO THE SUPPLIED PCA/CORRELATION RESULTS AND VOLCANO PLOT

  output$sample_pca_table <- renderTable({
    pca_coordinates[, c("Sample", "Group", "Condition", "Region", "Replicate", "PC1", "PC2", "PC3")]
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s", digits = 4,
  rownames = FALSE, width = "100%")


  output$sample_correlation_table <- renderTable({
    correlation_table <- data.frame(
      Sample = rownames(sample_correlation_matrix),
      round(sample_correlation_matrix, 4),
      check.names = FALSE,
      row.names = NULL
    )
    names(correlation_table)[-1] <- colnames(sample_correlation_matrix)
    correlation_table
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s",
  rownames = FALSE, width = "100%")


  output$deseq_volcano_table <- renderTable({
    shiny::validate(
      shiny::need(
        !is.null(sample_analysis),
        if (is.null(sample_analysis_error)) {
          "The full DESeq2 analysis could not be prepared."
        } else {
          sample_analysis_error
        }
      )
    )
    req(input$volcano_comparison)

    contrast_lookup <- list(
      "3dpa_Diaphysis_vs_Intact_Diaphysis" = c("group", "3dpa_Diaphysis", "Intact_Diaphysis"),
      "3dpa_Epiphysis_vs_Intact_Epiphysis" = c("group", "3dpa_Epiphysis", "Intact_Epiphysis"),
      "5dpa_Diaphysis_vs_Intact_Diaphysis" = c("group", "5dpa_Diaphysis", "Intact_Diaphysis"),
      "5dpa_Epiphysis_vs_Intact_Epiphysis" = c("group", "5dpa_Epiphysis", "Intact_Epiphysis"),
      "5dpa_vs_3dpa_Diaphysis" = c("group", "5dpa_Diaphysis", "3dpa_Diaphysis"),
      "5dpa_vs_3dpa_Epiphysis" = c("group", "5dpa_Epiphysis", "3dpa_Epiphysis")
    )
    selected_contrast <- contrast_lookup[[input$volcano_comparison]]
    shiny::validate(shiny::need(!is.null(selected_contrast), "Select one of the six supported DESeq2 comparisons."))

    result_table <- as.data.frame(DESeq2::results(
      sample_analysis$dds, contrast = selected_contrast
    ))
    result_table$gene_id <- rownames(result_table)
    result_table <- result_table[!is.na(result_table$log2FoldChange), , drop = FALSE]
    shiny::validate(shiny::need(nrow(result_table) > 0, "No testable genes were returned for this comparison."))

    result_table$regulation <- ifelse(
      !is.na(result_table$padj) & result_table$padj < 0.05 & result_table$log2FoldChange >= 1,
      "Upregulated",
      ifelse(
        !is.na(result_table$padj) & result_table$padj < 0.05 & result_table$log2FoldChange <= -1,
        "Downregulated", "Not significant"
      )
    )
    result_table <- result_table[, c(
      "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "regulation"
    )]
    names(result_table) <- c(
      "Gene ID", "Base mean", "Log2 fold change", "LFC standard error",
      "Wald statistic", "p-value", "Adjusted p-value", "Regulation"
    )
    result_table
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s",
  rownames = FALSE, digits = 4, width = "100%")


  # DISPLAY METADATA TABLE

  # ==========================================================
  # UPSET PLOT: OVERLAP OF SIGNIFICANT DEGs ACROSS SIX CONTRASTS
  # Uses the supplied Significant_DEGs.csv records and the same
  # thresholds as the user's regulation summary:
  # adjusted p-value < 0.05 and |log2FoldChange| >= 1.
  # ==========================================================

  upset_comparison_order <- c(
    "3dpa_Diaphysis_vs_Intact_Diaphysis",
    "3dpa_Epiphysis_vs_Intact_Epiphysis",
    "5dpa_Diaphysis_vs_Intact_Diaphysis",
    "5dpa_Epiphysis_vs_Intact_Epiphysis",
    "5dpa_vs_3dpa_Diaphysis",
    "5dpa_vs_3dpa_Epiphysis"
  )

  upset_comparison_labels <- c(
    "3 dpa Diaphysis vs Intact Diaphysis",
    "3 dpa Epiphysis vs Intact Epiphysis",
    "5 dpa Diaphysis vs Intact Diaphysis",
    "5 dpa Epiphysis vs Intact Epiphysis",
    "5 dpa vs 3 dpa Diaphysis",
    "5 dpa vs 3 dpa Epiphysis"
  )

  upset_data <- reactive({

    shiny::validate(
      shiny::need(
        is.null(comparison_data_error) && !is.null(significant_degs),
        if (is.null(comparison_data_error)) {
          "The supplied Significant_DEGs.csv file is not available."
        } else {
          comparison_data_error
        }
      ),
      shiny::need(
        all(upset_comparison_order %in% significant_degs$comparison),
        "The Significant_DEGs.csv file does not contain all six expected comparisons."
      )
    )

    upset_required <- significant_degs[
      significant_degs$comparison %in% upset_comparison_order &
        !is.na(significant_degs$gene_id) &
        nzchar(trimws(as.character(significant_degs$gene_id))) &
        !is.na(significant_degs$padj) &
        significant_degs$padj < 0.05 &
        !is.na(significant_degs$log2FoldChange) &
        abs(significant_degs$log2FoldChange) >= 1,
      ,
      drop = FALSE
    ]

    upset_required$gene_id <- trimws(as.character(upset_required$gene_id))

    upset_gene_sets <- stats::setNames(
      lapply(upset_comparison_order, function(comparison_name) {
        unique(upset_required$gene_id[
          upset_required$comparison == comparison_name
        ])
      }),
      upset_comparison_order
    )

    shiny::validate(
      shiny::need(
        all(lengths(upset_gene_sets) > 0),
        "At least one of the six comparisons has no genes meeting the DEG thresholds."
      )
    )

    all_upset_genes <- sort(unique(unlist(upset_gene_sets, use.names = FALSE)))

    membership_matrix <- vapply(
      upset_gene_sets,
      function(gene_set) all_upset_genes %in% gene_set,
      logical(length(all_upset_genes))
    )
    colnames(membership_matrix) <- upset_comparison_order

    membership_pattern <- apply(
      membership_matrix,
      1,
      function(membership) paste0(as.integer(membership), collapse = "")
    )

    genes_by_pattern <- split(all_upset_genes, membership_pattern)
    observed_patterns <- names(genes_by_pattern)

    intersection_labels <- vapply(
      observed_patterns,
      function(pattern) {
        included_sets <- which(strsplit(pattern, "", fixed = TRUE)[[1]] == "1")
        paste(upset_comparison_labels[included_sets], collapse = " + ")
      },
      character(1)
    )

    intersection_summary <- data.frame(
      pattern = observed_patterns,
      gene_count = as.integer(lengths(genes_by_pattern)[observed_patterns]),
      intersection_label = intersection_labels,
      stringsAsFactors = FALSE
    )

    intersection_summary <- intersection_summary[
      order(-intersection_summary$gene_count, intersection_summary$intersection_label),
      ,
      drop = FALSE
    ]
    rownames(intersection_summary) <- NULL

    set_size_summary <- data.frame(
      Comparison = upset_comparison_labels,
      Significant_DEGs = as.integer(lengths(upset_gene_sets)[upset_comparison_order]),
      stringsAsFactors = FALSE
    )

    list(
      filtered_results = upset_required,
      gene_sets = upset_gene_sets,
      genes_by_pattern = genes_by_pattern,
      intersection_summary = intersection_summary,
      top_intersections = utils::head(intersection_summary, 20L),
      set_size_summary = set_size_summary,
      total_unique_genes = length(all_upset_genes)
    )
  })


  output$upset_plot <- plotly::renderPlotly({

    upset_result <- upset_data()
    top_intersections <- upset_result$top_intersections
    number_of_intersections <- nrow(top_intersections)

    shiny::validate(
      shiny::need(number_of_intersections > 0, "No observed gene intersections are available to plot.")
    )

    intersection_positions <- seq_len(number_of_intersections)
    selected_pattern <- input$upset_intersection
    bar_colors <- if (is.null(selected_pattern) || !nzchar(selected_pattern)) {
      rep("#B98791", number_of_intersections)
    } else {
      ifelse(top_intersections$pattern == selected_pattern, "#6F3947", "#B98791")
    }

    # A stable Plotly source and per-bar key make each intersection clickable.
    upset_plot <- plotly::plot_ly(source = "upset_plot_source")

    upset_plot <- plotly::add_bars(
      upset_plot,
      x = intersection_positions,
      y = top_intersections$gene_count,
      marker = list(color = bar_colors),
      text = paste0(
        "<b>", top_intersections$gene_count, " genes</b><br>",
        top_intersections$intersection_label
      ),
      hoverinfo = "text",
      name = "Intersection size",
      showlegend = FALSE
    )

    # Draw a vertical connector through the active sets in each intersection.
    connector_x <- numeric(0)
    connector_y <- numeric(0)

    for (intersection_index in intersection_positions) {
      active_sets <- which(
        strsplit(
          top_intersections$pattern[intersection_index],
          "",
          fixed = TRUE
        )[[1]] == "1"
      )

      if (length(active_sets) > 1) {
        connector_x <- c(connector_x, intersection_index, intersection_index, NA)
        connector_y <- c(connector_y, min(active_sets), max(active_sets), NA)
      }
    }

    if (length(connector_x) > 0) {
      upset_plot <- plotly::add_trace(
        upset_plot,
        x = connector_x,
        y = connector_y,
        type = "scatter",
        mode = "lines",
        yaxis = "y2",
        line = list(color = "#59494A", width = 2),
        hoverinfo = "skip",
        showlegend = FALSE
      )
    }

    matrix_points <- expand.grid(
      intersection_index = intersection_positions,
      set_index = seq_along(upset_comparison_order),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    matrix_points$pattern <- top_intersections$pattern[
      matrix_points$intersection_index
    ]
    matrix_points$present <- vapply(
      seq_len(nrow(matrix_points)),
      function(i) {
        substr(
          matrix_points$pattern[i],
          matrix_points$set_index[i],
          matrix_points$set_index[i]
        ) == "1"
      },
      logical(1)
    )
    matrix_points$set_label <- upset_comparison_labels[matrix_points$set_index]
    matrix_points$intersection_label <- top_intersections$intersection_label[
      matrix_points$intersection_index
    ]

    inactive_points <- matrix_points[!matrix_points$present, , drop = FALSE]
    active_points <- matrix_points[matrix_points$present, , drop = FALSE]

    upset_plot <- plotly::add_trace(
      upset_plot,
      data = inactive_points,
      x = ~intersection_index,
      y = ~set_index,
      type = "scatter",
      mode = "markers",
      yaxis = "y2",
      marker = list(color = "#E8D5D2", size = 8),
      hoverinfo = "skip",
      showlegend = FALSE
    )

    upset_plot <- plotly::add_trace(
      upset_plot,
      data = active_points,
      x = ~intersection_index,
      y = ~set_index,
      type = "scatter",
      mode = "markers",
      yaxis = "y2",
      marker = list(color = "#59494A", size = 10),
      key = ~paste0("DOT::", pattern),
      text = ~paste0(
        "<b>", set_label, "</b><br>",
        "Intersection: ", intersection_label,
        "<br>Genes in this intersection: ",
        vapply(pattern, function(p) {
          length(upset_result$genes_by_pattern[[p]])
        }, integer(1)),
        "<br>",
        ifelse(
          vapply(pattern, function(p) sum(strsplit(p, "", fixed = TRUE)[[1]] == "1"), integer(1)) > 1,
          "Shared across comparisons — click to inspect genes",
          "Unique to this comparison — click has no effect"
        )
      ),
      hoverinfo = "text",
      showlegend = FALSE
    )

    plotly::layout(
      upset_plot,
      title = list(
        text = "Top 20 observed intersections across six comparisons",
        x = 0.02
      ),
      xaxis = list(
        title = "Intersection rank (hover over bars for comparison membership)",
        tickmode = "array",
        tickvals = intersection_positions,
        ticktext = as.character(intersection_positions),
        range = c(0.5, number_of_intersections + 0.5),
        showgrid = FALSE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = "Number of genes",
        domain = c(0.48, 1),
        rangemode = "tozero",
        showgrid = TRUE,
        gridcolor = "#EEE7E5",
        zeroline = FALSE
      ),
      yaxis2 = list(
        title = "",
        domain = c(0, 0.36),
        anchor = "x",
        tickmode = "array",
        tickvals = seq_along(upset_comparison_order),
        ticktext = upset_comparison_labels,
        autorange = "reversed",
        showgrid = FALSE,
        zeroline = FALSE
      ),
      margin = list(l = 270, r = 25, t = 75, b = 75),
      showlegend = FALSE,
      hovermode = "closest",
      clickmode = "event+select",
      paper_bgcolor = "#FFFDFC",
      plot_bgcolor = "#FFFDFC"
    ) |>
      plotly::config(displaylogo = FALSE, responsive = TRUE)
  })


  # Only active membership-matrix dots are clickable.
  # A dot in an intersection shared by 2+ comparisons selects that shared
  # intersection and updates the gene table. A dot in a single-comparison
  # intersection is a unique-gene set and intentionally does nothing.
  observeEvent(
    plotly::event_data("plotly_click", source = "upset_plot_source"),
    {
      clicked_point <- plotly::event_data(
        "plotly_click",
        source = "upset_plot_source"
      )

      if (!is.null(clicked_point$key) && length(clicked_point$key) > 0) {
        clicked_key <- as.character(clicked_point$key[[1]])

        if (startsWith(clicked_key, "DOT::")) {
          clicked_pattern <- sub("^DOT::", "", clicked_key)
          number_of_comparisons <- sum(
            strsplit(clicked_pattern, "", fixed = TRUE)[[1]] == "1"
          )

          if (number_of_comparisons > 1) {
            shiny::updateSelectInput(
              session,
              "upset_intersection",
              selected = clicked_pattern
            )
          }
        }
      }
    },
    ignoreInit = TRUE
  )


  output$upset_set_sizes_table <- renderTable({
    upset_data()$set_size_summary
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s",
  rownames = FALSE, width = "100%")


  output$upset_intersection_selector <- renderUI({
    upset_result <- upset_data()
    intersection_summary <- upset_result$intersection_summary

    choice_labels <- paste0(
      "Intersection ", seq_len(nrow(intersection_summary)),
      " (", intersection_summary$gene_count, " genes): ",
      intersection_summary$intersection_label
    )
    choices <- stats::setNames(intersection_summary$pattern, choice_labels)

    selectInput(
      "upset_intersection",
      "Select an intersection to inspect its genes",
      choices = choices,
      selected = intersection_summary$pattern[1],
      width = "100%"
    )
  })


  output$upset_intersection_table <- renderTable({
    upset_result <- upset_data()
    req(input$upset_intersection)

    selected_genes <- upset_result$genes_by_pattern[[input$upset_intersection]]
    shiny::validate(
      shiny::need(length(selected_genes) > 0, "No genes are available for this intersection.")
    )

    utils::head(
      data.frame(Gene_ID = selected_genes, stringsAsFactors = FALSE),
      100L
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s",
  rownames = FALSE, width = "100%")


  output$download_upset_genes <- downloadHandler(
    filename = function() {
      paste0("UpSet_intersection_", gsub("[^A-Za-z0-9]+", "_", input$upset_intersection), ".csv")
    },
    content = function(file) {
      upset_result <- upset_data()
      req(input$upset_intersection)
      selected_genes <- upset_result$genes_by_pattern[[input$upset_intersection]]
      selected_results <- upset_result$filtered_results[
        upset_result$filtered_results$gene_id %in% selected_genes,
        ,
        drop = FALSE
      ]
      selected_results <- selected_results[
        order(
          selected_results$gene_id,
          match(selected_results$comparison, upset_comparison_order)
        ),
        ,
        drop = FALSE
      ]
      write.csv(selected_results, file, row.names = FALSE)
    }
  )


  output$metadata_table <- renderTable({
    req(nrow(metadata) > 0)
    metadata
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")


  # DISPLAY THE 3 dpa vs 5 dpa SIGNIFICANT DEGs SUMMARY

  output$comparison_summary <- renderUI({

    if (!is.null(comparison_data_error)) {
      return(
        div(
          class = "alert alert-warning",
          comparison_data_error
        )
      )
    }

    req(input$comparison_region)

    comparison_name <- paste0(
      "5dpa_vs_3dpa_",
      input$comparison_region
    )

    comparison_rows <- significant_degs[
      significant_degs$comparison == comparison_name,
      ,
      drop = FALSE
    ]

    if (nrow(comparison_rows) == 0) {
      return(
        div(
          class = "alert alert-warning",
          paste("No rows were found for comparison:", comparison_name)
        )
      )
    }

    up_count <- sum(comparison_rows$regulation == "Upregulated", na.rm = TRUE)
    down_count <- sum(comparison_rows$regulation == "Downregulated", na.rm = TRUE)

    div(
      class = "comparison-summary",
      p(
        strong("Comparison: "),
        paste("5 dpa versus 3 dpa —", input$comparison_region)
      ),
      p(
        strong("Significant genes: "),
        format(nrow(comparison_rows), big.mark = ","),
        "  |  ",
        strong("Higher at 5 dpa: "),
        format(up_count, big.mark = ","),
        "  |  ",
        strong("Higher at 3 dpa: "),
        format(down_count, big.mark = ",")
      )
    )
  })


  # DISPLAY THE SIGNIFICANT DEGs FOR THE SELECTED 3 dpa vs 5 dpa REGION

  output$comparison_degs_table <- renderTable({

    shiny::validate(
      shiny::need(
        is.null(comparison_data_error),
        if (is.null(comparison_data_error)) "" else comparison_data_error
      )
    )

    req(input$comparison_region)

    comparison_name <- paste0(
      "5dpa_vs_3dpa_",
      input$comparison_region
    )

    comparison_rows <- significant_degs[
      significant_degs$comparison == comparison_name,
      ,
      drop = FALSE
    ]

    shiny::validate(
      shiny::need(
        nrow(comparison_rows) > 0,
        paste("No significant DEGs were found for", comparison_name)
      )
    )

    comparison_rows <- comparison_rows[
      order(comparison_rows$padj, na.last = TRUE),
      ,
      drop = FALSE
    ]

    display_columns <- c(
      "gene_id",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj",
      "regulation"
    )

    comparison_rows <- comparison_rows[, display_columns, drop = FALSE]

    # Keep small p-values readable instead of rounding them to 0.0000.
    comparison_rows$pvalue <- ifelse(
      is.na(comparison_rows$pvalue),
      NA_character_,
      formatC(comparison_rows$pvalue, format = "e", digits = 3)
    )

    comparison_rows$padj <- ifelse(
      is.na(comparison_rows$padj),
      NA_character_,
      formatC(comparison_rows$padj, format = "e", digits = 3)
    )

    names(comparison_rows) <- c(
      "Gene ID",
      "Base mean",
      "log2 fold change",
      "LFC standard error",
      "Wald statistic",
      "p-value",
      "Adjusted p-value",
      "Regulation"
    )

    rownames(comparison_rows) <- NULL
    comparison_rows

  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s",
  rownames = FALSE, digits = 4, width = "100%")


  # DISPLAY THE INTERACTIVE 3 dpa / 5 dpa LOG2 FOLD-CHANGE PROFILE
  # Values are taken from the existing Significant_DEGs.csv contrast results.

  output$expression_line_plot <- renderPlotly({

    shiny::validate(
      shiny::need(
        is.null(comparison_data_error),
        if (is.null(comparison_data_error)) "" else comparison_data_error
      )
    )

    req(input$expression_region)

    # If neither dropdown has a selection, show all six candidate genes.
    # Once a gene is selected in either dropdown, display the selected gene(s).
    selected_genes <- unique(c(
      input$expression_gene1,
      input$expression_gene2
    ))
    selected_genes <- selected_genes[
      !is.na(selected_genes) & nzchar(selected_genes)
    ]

    if (length(selected_genes) == 0) {
      selected_genes <- candidate_choices
    }

    tissue <- input$expression_region

    # Extract the gene symbol from IDs formatted as 'MSTRG.x | SYMBOL'.
    plot_source <- significant_degs
    plot_source$gene_symbol <- trimws(
      sub("^.*\\|\\s*", "", as.character(plot_source$gene_id))
    )

    get_contrast_lfc <- function(gene, timepoint) {
      comparison_name <- paste0(
        timepoint, "_", tissue, "_vs_Intact_", tissue
      )

      matching_rows <- plot_source[
        plot_source$gene_symbol == gene &
          plot_source$comparison == comparison_name,
        ,
        drop = FALSE
      ]

      if (nrow(matching_rows) == 0) {
        return(NA_real_)
      }

      matching_rows$log2FoldChange[[1]]
    }

    plot_object <- plot_ly()
    has_plottable_values <- FALSE

    for (gene in selected_genes) {

      value_3dpa <- get_contrast_lfc(gene, "3dpa")
      value_5dpa <- get_contrast_lfc(gene, "5dpa")

      if (!is.na(value_3dpa) || !is.na(value_5dpa)) {
        has_plottable_values <- TRUE
      }

      plot_object <- add_trace(
        plot_object,
        x = c("3 dpa", "5 dpa"),
        y = c(value_3dpa, value_5dpa),
        type = "scatter",
        mode = "lines+markers",
        name = gene,
        connectgaps = FALSE,
        marker = list(size = 9),
        line = list(width = 2),
        hovertemplate = paste0(
          "Gene: ", gene,
          "<br>Timepoint: %{x}",
          "<br>log2 fold change: %{y:.3f}<extra></extra>"
        )
      )
    }

    shiny::validate(
      shiny::need(
        has_plottable_values,
        "No significant-DEG values are available for the selected genes, timepoints, and tissue region."
      )
    )

    layout(
      plot_object,
      xaxis = list(
        title = "Timepoint",
        categoryorder = "array",
        categoryarray = c("3 dpa", "5 dpa")
      ),
      yaxis = list(
        title = "log2 fold change vs Intact",
        zeroline = TRUE,
        zerolinewidth = 1
      ),
      legend = list(title = list(text = "Candidate gene")),
      hovermode = "closest",
      margin = list(l = 70, r = 30, b = 65, t = 25)
    )
  })


  # DISPLAY TOP 10 UPREGULATED GENES TABLE

  output$upregulated_table <- renderTable({

    shiny::validate(
      shiny::need(
        !is.null(upregulated_genes),
        paste0(
          "Upregulated genes CSV not found. Place '",
          upregulated_file,
          "' in the same folder as app.R, or update upregulated_file to your CSV filename."
        )
      )
    )

    shiny::validate(
      shiny::need(
        nrow(upregulated_genes) > 0,
        "The upregulated genes CSV does not contain any rows."
      )
    )

    upregulated_genes

  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")


  # DISPLAY TOP 10 DOWNREGULATED GENES TABLE

  output$downregulated_table <- renderTable({

    shiny::validate(
      shiny::need(
        !is.null(downregulated_genes),
        paste0(
          "Downregulated genes CSV not found. Place '",
          downregulated_file,
          "' in the same folder as app.R, or update downregulated_file to your CSV filename."
        )
      )
    )

    shiny::validate(
      shiny::need(
        nrow(downregulated_genes) > 0,
        "The downregulated genes CSV does not contain any rows."
      )
    )

    downregulated_genes

  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")


  # GENE MAPPING SUMMARY VALUES

  output$mapping_input_count <- renderText({
    req(!is.null(mapping_data))
    length(unique(mapping_data$axolotl_annotation))
  })

  output$mapping_human_record_count <- renderText({
    req(!is.null(mapping_data))
    nrow(human_mapping)
  })


  # FILTER THE GENE MAPPING RESULTS

  filtered_mapping_data <- reactive({

    req(!is.null(mapping_data))

    results <- mapping_data

    search_text <- if (is.null(input$mapping_search)) {
      ""
    } else {
      trimws(input$mapping_search)
    }

    if (nzchar(search_text)) {

      searchable_columns <- c(
        "axolotl_annotation",
        "axolotl_gene_id",
        "description",
        "ncbi_gene_id",
        "human_gene_symbol",
        "human_gene_description"
      )

      search_match <- Reduce(
        `|`,
        lapply(searchable_columns, function(column_name) {

          grepl(
            search_text,
            ifelse(
              is.na(results[[column_name]]),
              "",
              as.character(results[[column_name]])
            ),
            ignore.case = TRUE,
            fixed = TRUE
          )
        })
      )

      search_match[is.na(search_match)] <- FALSE

      results <- results[
        search_match,
        ,
        drop = FALSE
      ]
    }

    selected_status <- if (is.null(input$mapping_status_filter)) {
      "All statuses"
    } else {
      input$mapping_status_filter
    }

    if (!identical(selected_status, "All statuses")) {

      results <- results[
        results$mapping_status == selected_status,
        ,
        drop = FALSE
      ]
    }

    selected_match_type <- if (is.null(input$mapping_match_filter)) {
      "All match types"
    } else {
      input$mapping_match_filter
    }

    if (!identical(selected_match_type, "All match types")) {

      results <- results[
        results$mapping_match_type == selected_match_type,
        ,
        drop = FALSE
      ]
    }

    results
  })


  # DISPLAY THE NUMBER OF FILTERED RESULTS

  output$mapping_filtered_count <- renderText({

    req(!is.null(mapping_data))

    paste(
      "Showing",
      nrow(filtered_mapping_data()),
      "of",
      nrow(mapping_data),
      "human gene records"
    )
  })


  # DISPLAY FILTERED GENE MAPPING RESULTS

  output$mapping_results_table <- renderTable({

    req(!is.null(mapping_data))

    filtered_mapping_data()

  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")


  # DOWNLOAD FILTERED GENE MAPPING RESULTS

  output$download_filtered_mapping <- downloadHandler(

    filename = function() {
      paste0(
        "Axolotl_Human_Gene_Mapping_Filtered_",
        Sys.Date(),
        ".csv"
      )
    },

    content = function(file) {
      write.csv(
        filtered_mapping_data(),
        file,
        row.names = FALSE,
        na = ""
      )
    }
  )


  # DOWNLOAD COMPLETE GENE MAPPING RESULTS

  output$download_all_mapping <- downloadHandler(

    filename = function() {
      paste0(
        "Axolotl_Human_Gene_Mapping_All_",
        Sys.Date(),
        ".csv"
      )
    },

    content = function(file) {

      req(!is.null(mapping_data))

      write.csv(
        mapping_data,
        file,
        row.names = FALSE,
        na = ""
      )
    }
  )


  # ==========================================================
  # FUNCTIONAL ENRICHMENT PLOT OUTPUTS
  # ==========================================================


  # ----------------------------------------------------------
  # KEGG PLOT
  # Only candidates with retrieved KEGG pathway mappings
  # are included in the visualisations.
  # ----------------------------------------------------------

  output$kegg_plot <- renderPlotly({

    req(input$kegg_candidates)

    shiny::validate(
      shiny::need(
        !is.null(KEGG_data),
        "KEGG_analysis.csv was not found or could not be read."
      ),

      shiny::need(
        all(c(
          "candidate",
          "annotation_id",
          "annotation"
        ) %in% names(KEGG_data)),
        "The KEGG file must contain candidate, annotation_id and annotation columns."
      ),

      shiny::need(
        length(input$kegg_candidates) > 0,
        "Select a candidate gene."
      )
    )

    # Only candidates with retrieved KEGG pathway mappings
    # are eligible for the KEGG visualisations.

    kegg_mapped_candidates <- c(
      "TGFB2",
      "VTN"
    )

    selected_kegg_candidates <- intersect(
      input$kegg_candidates,
      kegg_mapped_candidates
    )

    shiny::validate(
      shiny::need(
        length(selected_kegg_candidates) > 0,
        paste(
          "No KEGG pathway mappings are available for the selected candidate.",
          "Only TGFB2 and VTN have retrieved KEGG pathway associations."
        )
      )
    )

    plot_data <- get_selected_candidates(
      KEGG_data,
      selected_kegg_candidates
    )

    # Remove missing or blank pathway annotations.

    plot_data <- plot_data[
      !is.na(plot_data$annotation) &
        nzchar(trimws(as.character(plot_data$annotation))),
      ,
      drop = FALSE
    ]

    plot_data <- classify_kegg_relevance(plot_data)
    selected_relevance <- if (is.null(input$kegg_relevance)) "Direct" else input$kegg_relevance
    plot_data <- plot_data[
      plot_data$regeneration_relevance == selected_relevance,
      ,
      drop = FALSE
    ]

    shiny::validate(
      shiny::need(
        nrow(plot_data) > 0,
        paste(
          "No", tolower(selected_relevance),
          "regeneration-relevance KEGG pathways were found for the selected candidate."
        )
      )
    )

    selected_graph <- input$kegg_graph

    if (
      is.null(selected_graph) ||
      !nzchar(selected_graph)
    ) {
      selected_graph <- "Bar plot"
    }

    if (selected_graph == "Bar plot") {

      make_kegg_bar(plot_data)

    } else if (selected_graph == "Network plot") {

      make_annotation_network(plot_data)

    } else {

      make_kegg_bar(plot_data)
    }
  })


  # ----------------------------------------------------------
  # GO PLOT
  # ----------------------------------------------------------

  output$go_plot <- renderPlotly({

    req(input$go_candidates, input$go_graph)

    shiny::validate(
      shiny::need(
        !is.null(GO_data),
        "GO_analysis.csv was not found or could not be read."
      ),

      shiny::need(
        all(c(
          "candidate",
          "annotation_id",
          "annotation"
        ) %in% names(GO_data)),
        "The GO file must contain candidate, annotation_id and annotation columns."
      ),

      shiny::need(
        length(input$go_candidates) > 0,
        "Select a candidate gene."
      )
    )

    plot_data <- get_selected_candidates(
      GO_data,
      input$go_candidates
    )

    plot_data <- classify_go_relevance(plot_data)
    selected_relevance <- if (is.null(input$go_relevance)) "Direct" else input$go_relevance
    plot_data <- plot_data[
      plot_data$regeneration_relevance == selected_relevance,
      ,
      drop = FALSE
    ]

    shiny::validate(
      shiny::need(
        nrow(plot_data) > 0,
        paste(
          "No", tolower(selected_relevance),
          "regeneration-relevance GO annotations were found for the selected candidate."
        )
      )
    )

    if (input$go_graph == "Gene treemap") {

      make_go_gene_treemap(plot_data)

    } else {

      make_annotation_bar(
        plot_data,
        input$go_candidates,
        "annotation"
      )
    }
  })


  # ----------------------------------------------------------
  # FOCUSED PATHWAYS PLOT
  # ----------------------------------------------------------

  output$pathway_plot <- renderPlotly({

    req(input$pathway_graph)

    shiny::validate(
      shiny::need(
        !is.null(Reactome_individual_data),
        "Reactome_individual_gene_pathways.csv was not found or could not be read."
      ),
      shiny::need(
        all(c(
          "gene",
          "pathway_id",
          "pathway_name",
          "pathway_hierarchy",
          "top_level_pathway"
        ) %in% names(Reactome_individual_data)),
        "The Reactome individual pathway file is missing one or more required columns."
      )
    )

    pathway_plot_data <- make_reactome_plot_data(
      Reactome_individual_data,
      Reactome_hierarchy_data,
      candidate_choices
    )

    shiny::validate(
      shiny::need(
        nrow(pathway_plot_data) > 0,
        "No Reactome gene-to-pathway associations were found for the six candidate genes."
      )
    )

    if (input$pathway_graph == "Chord diagram") {

      make_pathway_chord(pathway_plot_data)

    } else if (input$pathway_graph == "Sunburst plot") {

      make_pathway_sunburst(pathway_plot_data)

    } else {

      make_annotation_heatmap(pathway_plot_data)
    }
  })



  # ----------------------------------------------------------
  # INDIVIDUAL REACTOME PATHWAYS TABLE
  # ----------------------------------------------------------

  output$individual_pathway_table <- renderTable({

    req(input$individual_pathway_candidate, input$individual_pathway_top_level)

    shiny::validate(
      shiny::need(
        !is.null(Reactome_individual_data),
        "Reactome_individual_gene_pathways.csv was not found or could not be read."
      ),
      shiny::need(
        all(c(
          "gene",
          "uniprot_id",
          "pathway_id",
          "pathway_name",
          "pathway_species",
          "pathway_type",
          "has_diagram",
          "is_inferred",
          "pathway_hierarchy",
          "top_level_pathway",
          "top_level_pathway_id"
        ) %in% names(Reactome_individual_data)),
        "The Reactome pathway file is missing one or more required columns."
      )
    )

    individual_data <- Reactome_individual_data[
      !is.na(Reactome_individual_data$gene) &
        trimws(as.character(Reactome_individual_data$gene)) ==
          input$individual_pathway_candidate,
      ,
      drop = FALSE
    ]

    if (input$individual_pathway_top_level != "__ALL__") {
      individual_data <- individual_data[
        !is.na(individual_data$top_level_pathway) &
          trimws(as.character(individual_data$top_level_pathway)) ==
            input$individual_pathway_top_level,
        ,
        drop = FALSE
      ]
    }

    if (nrow(individual_data) == 0) {
      selected_group <- if (
        input$individual_pathway_top_level == "__ALL__"
      ) {
        "any top-level pathway"
      } else {
        input$individual_pathway_top_level
      }

      return(data.frame(
        Message = paste(
          "No Reactome pathway records were found for",
          input$individual_pathway_candidate,
          "under",
          selected_group,
          "."
        ),
        check.names = FALSE
      ))
    }

    individual_data <- individual_data[
      order(
        tolower(as.character(individual_data$top_level_pathway)),
        tolower(as.character(individual_data$pathway_hierarchy)),
        tolower(as.character(individual_data$pathway_name))
      ),
      ,
      drop = FALSE
    ]

    data.frame(
      `Pathway ID` = as.character(individual_data$pathway_id),
      `Pathway name` = clean_annotation_labels(individual_data$pathway_name),
      `Pathway hierarchy` = as.character(individual_data$pathway_hierarchy),
      `Top-level pathway` = as.character(individual_data$top_level_pathway),
      `Top-level pathway ID` = as.character(individual_data$top_level_pathway_id),
      `UniProt accession` = as.character(individual_data$uniprot_id),
      Species = as.character(individual_data$pathway_species),
      `Pathway type` = as.character(individual_data$pathway_type),
      `Has diagram` = as.character(individual_data$has_diagram),
      Inferred = as.character(individual_data$is_inferred),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  },
  striped = TRUE,
  bordered = TRUE,
  hover = TRUE,
  spacing = "s",
  width = "100%"
  )


  # ==========================================================
  # TRACK CURRENT PAGE
  # ==========================================================

  current_page <- reactiveVal("home")

  workflow_selected <- reactiveVal("prep")

  gene_mapping_page <- reactiveVal("overview")

  enrichment_selected <- reactiveVal("overview")

  # Results & Interpretation sub-navigation
  results_section <- reactiveVal("overview")


  # ==========================================================
  # MAIN PAGE CONTENT
  # ==========================================================

  output$home_content <- renderUI({

    page <- current_page()


    # ABOUT THE PROJECT

    if (page == "about") {

      return(
        tagList(

          h2("About the Project"),

          p("Explore the dataset and the analysis workflow."),

          layout_columns(

            div(
              class = "workflow-card",

              div(
                class = "card-symbol",
                icon("table")
              ),

              h4("Metadata"),

              p("View the sample information and experimental groups."),

              actionButton(
                "open_metadata",
                "Open Metadata",
                icon = icon("arrow-right"),
                class = "btn-primary"
              )
            ),

            div(
              class = "workflow-card",

              div(
                class = "card-symbol",
                icon("diagram-project")
              ),

              h4("Complete Workflow"),

              p("View the stages of the RNA-seq analysis."),

              actionButton(
                "open_workflow",
                "View Workflow",
                icon = icon("arrow-right"),
                class = "btn-primary"
              )
            ),

            col_widths = c(6, 6)
          )
        )
      )
    }


    # METADATA PAGE

    if (page == "metadata") {

      return(
        tagList(

          actionButton(
            "back_about_metadata",
            "Back to About the Project",
            icon = icon("arrow-left"),
            class = "btn-light"
          ),

          br(),
          br(),

          h2("Sample Metadata"),

          tableOutput("metadata_table")
        )
      )
    }


    # COMPLETE WORKFLOW PAGE

    if (page == "workflow") {

      selected <- workflow_selected()

      workflow_stages <- list(

        prep = list(
          title = "RNA-seq preparation",
          what = paste(
            "The sample metadata and gene-count data provide the inputs",
            "for downstream RNA-seq analysis."
          ),
          why = paste(
            "These inputs help connect gene-expression",
            "to the correct samples and experimental groups."
          )
        ),

        deseq = list(
          title = "DESeq2 (Differential expression)",
          what = paste(
            "DESeq2 analyses gene-count data to identify differences",
            "in gene expression between the groups being compared."
          ),
          why = paste(
            "It helps identify genes whose expression changes",
            "between the selected experimental conditions."
          )
        ),

        significant = list(
          title = "Significant DESeq2 results",
          what = paste(
            "The DESeq2 results are reviewed to identify genes that",
            "meet the significance criteria used in the analysis."
          ),
          why = paste(
            "Filtering helps focus subsequent investigation on genes",
            "with evidence of differential expression."
          )
        ),

        go = list(
          title = "Gene Ontology (GO)",
          what = paste(
            "Gene Ontology analysis examines whether genes are associated",
            "with biological processes, molecular functions, or cellular",
            "components more often than expected."
          ),
          why = paste(
            "It helps interpret gene lists in terms of biological",
            "functions and processes. It helps to identify top",
            "upregulated and downregulated genes."
          )
        ),

        human = list(
          title = "Human Gene Mapping",
          what = paste(
            "This stage maps axolotl gene identifiers to corresponding",
            "human gene identifiers or symbols where a suitable mapping",
            "is available."
          ),
          why = paste(
            "Mapped identifiers can make it easier to use human-focused",
            "annotation resources. Mappings need to be checked rather",
            "than assumed to be one-to-one."
          )
        ),

        enrichment = list(
          title = "Functional Enrichment",
          what = paste(
            "Functional enrichment analyses examine whether a gene list",
            "is associated with particular biological functions or pathways."
          ),
          why = paste(
            "It helps summarise gene-level results into broader biological",
            "themes. Statistical significance should be checked in the",
            "enrichment results before drawing conclusions."
          )
        ),

        visualization = list(
          title = "Visualisation",
          what = paste(
            "Visual summaries present selected results in plots and",
            "other graphical formats."
          ),
          why = paste(
            "Plots can make expression patterns and analysis results",
            "easier to inspect and compare."
          )
        )
      )


      return(
        tagList(

          actionButton(
            "back_about_workflow",
            "Back to About the Project",
            icon = icon("arrow-left"),
            class = "btn-light"
          ),

          br(),
          br(),

          h2("Complete RNA-seq Workflow"),

          div(
            class = "workflow-layout",

            div(
              class = "workflow-diagram",

              tags$svg(
                class = "workflow-connectors",
                viewBox = "0 0 1000 620",
                preserveAspectRatio = "none",

                tags$defs(
                  tags$marker(
                    id = "workflow-arrowhead",
                    markerWidth = "7",
                    markerHeight = "7",
                    refX = "6",
                    refY = "3",
                    orient = "auto",
                    markerUnits = "userSpaceOnUse",

                    tags$path(
                      d = "M0,0 L6,3 L0,6 Z",
                      fill = "#71696A"
                    )
                  )
                ),

                tags$g(
                  fill = "none",
                  stroke = "#71696A",
                  `stroke-width` = "4",
                  `stroke-linejoin` = "round",
                  `stroke-linecap` = "round",

                  tags$path(
                    d = "M500 70 V116",
                    `marker-end` = "url(#workflow-arrowhead)"
                  ),

                  tags$path(
                    d = "M350 147 H200 V220",
                    `marker-end` = "url(#workflow-arrowhead)"
                  ),

                  tags$path(
                    d = "M650 147 H800 V220",
                    `marker-end` = "url(#workflow-arrowhead)"
                  ),

                  tags$path(
                    d = "M800 278 V340 H650",
                    `marker-end` = "url(#workflow-arrowhead)"
                  ),

                  tags$path(
                    d = "M500 368 V414",
                    `marker-end` = "url(#workflow-arrowhead)"
                  ),

                  tags$path(
                    d = "M500 473 V519",
                    `marker-end` = "url(#workflow-arrowhead)"
                  )
                )
              ),

              div(
                class = "workflow-node workflow-prep",

                actionButton(
                  "workflow_stage_prep",
                  "RNA-seq preparation",
                  class = paste(
                    "workflow-node-btn btn",
                    if (selected == "prep")
                      "workflow-node-selected"
                    else ""
                  )
                )
              ),

              div(
                class = "workflow-node workflow-deseq",

                actionButton(
                  "workflow_stage_deseq",
                  "DESeq2",
                  class = paste(
                    "workflow-node-btn btn",
                    if (selected == "deseq")
                      "workflow-node-selected"
                    else ""
                  )
                )
              ),

              div(
                class = "workflow-node workflow-significant",

                actionButton(
                  "workflow_stage_significant",
                  "Significant DESeq2",
                  class = paste(
                    "workflow-node-btn btn",
                    if (selected == "significant")
                      "workflow-node-selected"
                    else ""
                  )
                )
              ),

              div(
                class = "workflow-node workflow-go",

                actionButton(
                  "workflow_stage_go",
                  "Gene Ontology (GO)",
                  class = paste(
                    "workflow-node-btn btn",
                    if (selected == "go")
                      "workflow-node-selected"
                    else ""
                  )
                )
              ),

              div(
                class = "workflow-node workflow-human",

                actionButton(
                  "workflow_stage_human",
                  "Human Gene Mapping",
                  class = paste(
                    "workflow-node-btn btn",
                    if (selected == "human")
                      "workflow-node-selected"
                    else ""
                  )
                )
              ),

              div(
                class = "workflow-node workflow-enrichment",

                actionButton(
                  "workflow_stage_enrichment",
                  "Functional Enrichment",
                  class = paste(
                    "workflow-node-btn btn",
                    if (selected == "enrichment")
                      "workflow-node-selected"
                    else ""
                  )
                )
              ),

              div(
                class = "workflow-node workflow-visualization",

                actionButton(
                  "workflow_stage_visualization",
                  "Visualisation",
                  class = paste(
                    "workflow-node-btn btn",
                    if (selected == "visualization")
                      "workflow-node-selected"
                    else ""
                  )
                )
              )
            ),

            div(
              class = "workflow-detail",

              div(
                class = "workflow-detail-label",
                "Selected stage"
              ),

              h3(workflow_stages[[selected]]$title),

              p(workflow_stages[[selected]]$what),

              div(
                class = "workflow-why",

                p(
                  strong("Why this stage matters"),
                  br(),
                  workflow_stages[[selected]]$why
                )
              )
            )
          )
        )
      )
    }


    # DIFFERENTIAL EXPRESSION PAGE

    if (page == "de") {

      return(
        tagList(

          h2(
            "Differential Expression",
            style = paste(
              "font-weight:800;",
              "color:#242020;",
              "margin-bottom:22px;"
            )
          ),

          layout_columns(

            div(
              class = "workflow-card de-column",

              div(
                class = "card-symbol",
                icon("table")
              ),

              h4("Comparison Table"),

              p("View comparisons between experimental groups."),

              actionButton(
                "de_comparison",
                "View Comparisons",
                icon = icon("arrow-right"),
                class = "btn-primary"
              )
            ),

            div(
              class = "workflow-card de-column",

              div(
                class = "card-symbol",
                icon("dna")
              ),

              h4("Regulated Genes"),

              p("Explore upregulated and downregulated genes."),

              actionButton(
                "de_regulated",
                "Explore Genes",
                icon = icon("arrow-right"),
                class = "btn-primary"
              )
            ),

            col_widths = c(6, 6)
          )
        )
      )
    }


    # DESeq2 COMPARISON LANDING PAGE

    if (page == "comparison") {

      return(
        tagList(

          actionButton(
            "back_de_comparison",
            "Back to Differential Expression",
            icon = icon("arrow-left"),
            class = "btn-light"
          ),

          br(),
          br(),

          h2(
            "DESeq2 Comparison",
            style = paste(
              "font-weight:800;",
              "color:#242020;",
              "margin-bottom:10px;"
            )
          ),

          p(
            "Choose how you want to explore the 3 dpa versus 5 dpa results.",
            style = "color:#59494A;line-height:1.7;margin-bottom:22px;"
          ),

          layout_columns(

            actionButton(
              "open_comparison_table",
              label = tagList(
                div(
                  class = "card-symbol",
                  icon("table")
                ),
                h4("Table"),
                p("View significant differential expression results for 5 dpa versus 3 dpa.")
              ),
              class = "workflow-card de-column comparison-option-button"
            ),

            actionButton(
              "open_expression_line_plot",
              label = tagList(
                div(
                  class = "card-symbol",
                  icon("chart-line")
                ),
                h4("Expression Line Plot"),
                p("Open the separate expression line plot page.")
              ),
              class = "workflow-card de-column comparison-option-button"
            ),

            col_widths = c(6, 6)
          )
        )
      )
    }


    # 3 dpa vs 5 dpa COMPARISON TABLE PAGE

    if (page == "comparison_table") {

      return(
        tagList(

          actionButton(
            "back_comparison_options",
            "Back to DESeq2 Comparison",
            icon = icon("arrow-left"),
            class = "btn-light"
          ),

          br(),
          br(),

          h2(
            "3 vs 5 dpa — Significant DEGs",
            style = paste(
              "font-weight:800;",
              "color:#242020;",
              "margin-bottom:10px;"
            )
          ),

          p(
            paste(
              "Positive log2FoldChange values indicate higher expression",
              "at 5 dpa; negative values indicate higher expression at 3 dpa."
            ),
            style = "color:#59494A;line-height:1.7;margin-bottom:18px;"
          ),

          div(
            class = "workflow-card",
            selectInput(
              "comparison_region",
              "Select tissue region",
              choices = c(
                "Diaphysis" = "Diaphysis",
                "Epiphysis" = "Epiphysis"
              ),
              selected = "Diaphysis",
              multiple = FALSE
            ),

            uiOutput("comparison_summary")
          ),

          br(),

          div(
            class = "workflow-card",
            h4("Significant differential expression results"),
            p(
              "Rows are ordered by adjusted p-value (padj), from smallest to largest.",
              style = "font-size:13px;color:#59494A;"
            ),
            div(
              id = "comparison_table_scroll",
              class = "comparison-table-scroll",
              tableOutput("comparison_degs_table"),
              div(
                id = "comparison_table_bottom_scroll",
                class = "comparison-table-bottom-scroll",
                div(
                  id = "comparison_table_bottom_scroll_inner",
                  class = "comparison-table-bottom-scroll-inner"
                )
              )
            ),
            tags$script(HTML("
              (function() {
                function setupComparisonTableScroll() {
                  var bottom = document.getElementById('comparison_table_bottom_scroll');
                  var bottomInner = document.getElementById('comparison_table_bottom_scroll_inner');
                  var tableArea = document.getElementById('comparison_table_scroll');
                  var tableOutput = document.getElementById('comparison_degs_table');
                  if (!bottom || !bottomInner || !tableArea || !tableOutput) return;
                  if (bottom.dataset.syncReady === 'true') return;
                  bottom.dataset.syncReady = 'true';

                  var updateWidth = function() {
                    var table = tableOutput.querySelector('table');
                    if (!table) return;
                    var width = Math.max(table.scrollWidth, table.offsetWidth, tableArea.clientWidth);
                    bottomInner.style.width = width + 'px';
                    bottom.scrollLeft = tableArea.scrollLeft;
                  };

                  var syncing = false;
                  bottom.addEventListener('scroll', function() {
                    if (syncing) return;
                    syncing = true;
                    tableArea.scrollLeft = bottom.scrollLeft;
                    syncing = false;
                  });
                  tableArea.addEventListener('scroll', function() {
                    if (syncing) return;
                    syncing = true;
                    bottom.scrollLeft = tableArea.scrollLeft;
                    syncing = false;
                  });

                  if (window.ResizeObserver) {
                    var resizeObserver = new ResizeObserver(updateWidth);
                    resizeObserver.observe(tableOutput);
                  }
                  if (window.MutationObserver) {
                    var mutationObserver = new MutationObserver(updateWidth);
                    mutationObserver.observe(tableOutput, { childList: true, subtree: true });
                  }
                  updateWidth();
                  setTimeout(updateWidth, 100);
                }
                if (document.readyState === 'loading') {
                  document.addEventListener('DOMContentLoaded', setupComparisonTableScroll);
                } else {
                  setupComparisonTableScroll();
                }
                document.addEventListener('shiny:connected', setupComparisonTableScroll);
                document.addEventListener('shiny:value', function(event) {
                  if (event.target && event.target.id === 'comparison_degs_table') {
                    setupComparisonTableScroll();
                  }
                });
              })();
            "))
          )
        )
      )
    }


    # EXPRESSION LINE PLOT PAGE

    if (page == "expression_line_plot") {

      return(
        tagList(

          actionButton(
            "back_comparison_options",
            "Back to DESeq2 Comparison",
            icon = icon("arrow-left"),
            class = "btn-light"
          ),

          br(),
          br(),

          h2(
            "Expression Line Plot",
            style = paste(
              "font-weight:800;",
              "color:#242020;",
              "margin-bottom:10px;"
            )
          ),

          p(
            "The plot initially displays all six candidate genes. Use Gene 1 and Gene 2 to select up to two genes for a focused comparison; each dropdown contains all six candidate genes.",
            style = "color:#59494A;line-height:1.7;margin-bottom:20px;"
          ),

          div(
            class = "workflow-card",
            h4("Select genes and tissue"),

            div(
              style = paste(
                "display:grid;",
                "grid-template-columns:repeat(2,minmax(0,1fr));",
                "gap:16px;",
                "width:100%;",
                "align-items:start;"
              ),
              div(
                h5("Gene 1"),
                selectInput(
                  "expression_gene1",
                  label = NULL,
                  choices = c("Select a gene" = "", candidate_choices),
                  selected = "",
                  width = "100%"
                )
              ),
              div(
                h5("Gene 2"),
                selectInput(
                  "expression_gene2",
                  label = NULL,
                  choices = c("Select a gene" = "", candidate_choices),
                  selected = "",
                  width = "100%"
                )
              )
            ),

            selectInput(
              "expression_region",
              "Select tissue region",
              choices = c(
                "Diaphysis" = "Diaphysis",
                "Epiphysis" = "Epiphysis"
              ),
              selected = "Diaphysis",
              multiple = FALSE
            )
          ),

          br(),

          div(
            class = "workflow-card",
            h4("Candidate genes: log2 fold-change profile"),
            plotlyOutput("expression_line_plot", height = "500px")
          )
        )
      )
    }


    # REGULATED GENES PAGE

    if (page == "regulated") {

      return(
        tagList(

          actionButton(
            "back_de_regulated",
            "Back to Differential Expression",
            icon = icon("arrow-left"),
            class = "btn-light"
          ),

          br(),
          br(),

          h2(
            "Regulated Genes",
            style = paste(
              "font-weight:800;",
              "color:#242020;",
              "margin-bottom:22px;"
            )
          ),

          layout_columns(

            div(
              class = "workflow-card de-column",

              div(
                class = "card-symbol",
                icon("arrow-trend-up")
              ),

              h4("Upregulated Genes"),

              p("Explore genes with increased expression."),

              actionButton(
                "view_upregulated",
                "View the genes",
                icon = icon("arrow-right"),
                class = "btn-primary"
              )
            ),

            div(
              class = "workflow-card de-column",

              div(
                class = "card-symbol",
                icon("arrow-trend-down")
              ),

              h4("Downregulated Genes"),

              p("Explore genes with decreased expression."),

              actionButton(
                "view_downregulated",
                "View the genes",
                icon = icon("arrow-right"),
                class = "btn-primary"
              )
            ),

            col_widths = c(6, 6)
          )
        )
      )
    }


    # TOP 10 UPREGULATED GENES TABLE PAGE

    if (page == "upregulated") {

      return(
        tagList(

          actionButton(
            "back_regulated_up",
            "Back to Regulated Genes",
            icon = icon("arrow-left"),
            class = "btn-light"
          ),

          br(),
          br(),

          h2(
            "Top 10 Upregulated Genes",
            style = paste(
              "font-weight:800;",
              "color:#242020;",
              "margin-bottom:22px;"
            )
          ),

          div(
            class = "workflow-card",
            tableOutput("upregulated_table")
          )
        )
      )
    }


    # TOP 10 DOWNREGULATED GENES TABLE PAGE

    if (page == "downregulated") {

      return(
        tagList(

          actionButton(
            "back_regulated_down",
            "Back to Regulated Genes",
            icon = icon("arrow-left"),
            class = "btn-light"
          ),

          br(),
          br(),

          h2(
            "Top 10 Downregulated Genes",
            style = paste(
              "font-weight:800;",
              "color:#242020;",
              "margin-bottom:22px;"
            )
          ),

          div(
            class = "workflow-card",
            tableOutput("downregulated_table")
          )
        )
      )
    }


    # ========================================================
    # FUNCTIONAL ENRICHMENT PAGE
    # ========================================================

    if (page == "enrichment") {

      selected_enrichment <- enrichment_selected()


      # FUNCTIONAL ENRICHMENT LANDING PAGE

      if (selected_enrichment == "overview") {

        return(
          tagList(

            h2(
              "Functional Enrichment",
              style = paste(
                "font-weight:800;",
                "color:#242020;",
                "margin-bottom:24px;"
              )
            ),

            div(
              class = "enrichment-layout",

              actionButton(
                "open_kegg",

                label = tagList(
                  div(
                    class = "card-symbol",
                    icon("diagram-project")
                  ),
                  h3("KEGG"),
                  p(
                    "Explore KEGG pathway annotations associated with the candidate genes."
                  )
                ),

                class = "enrichment-option-button"
              ),

              actionButton(
                "open_go",

                label = tagList(
                  div(
                    class = "card-symbol",
                    icon("sitemap")
                  ),
                  h3("Gene Ontology (GO)"),
                  p(
                    "Explore biological process, cellular component and molecular function annotations."
                  )
                ),

                class = "enrichment-option-button"
              ),

              actionButton(
                "open_pathways",

                label = tagList(
                  div(
                    class = "card-symbol",
                    icon("route")
                  ),
                  h3("Pathways"),
                  p(
                    "Explore the focused KEGG pathway records associated with the selected candidates."
                  )
                ),

                class = "enrichment-option-button"
              )
            )
          )
        )
      }


      # ------------------------------------------------------
      # KEGG SECTION
      # ------------------------------------------------------

      if (selected_enrichment == "kegg") {

        return(
          tagList(

            actionButton(
              "back_enrichment",
              "Back to Functional Enrichment",
              icon = icon("arrow-left"),
              class = "btn-light"
            ),

            br(),
            br(),

            h2(
              icon("diagram-project"),
              " KEGG",
              style = paste(
                "font-weight:800;",
                "color:#242020;",
                "margin-bottom:18px;"
              )
            ),

            div(
              class = "enrichment-detail",

              div(
                class = "enrichment-panel",

                p(
                  "Explore KEGG pathway annotations associated with the candidate genes.",
                  class = "enrichment-note"
                ),

                div(
                  style = paste(
                    "background:#F0E3DF;",
                    "border-left:4px solid #B98791;",
                    "border-radius:10px;",
                    "padding:15px 18px;",
                    "margin-bottom:20px;"
                  ),

                  h4(
                    "KEGG mapping note",
                    style = paste(
                      "font-size:15px;",
                      "font-weight:800;",
                      "color:#242020;",
                      "margin-top:0;"
                    )
                  ),

                  p(
                    paste(
                      "FLRT3, MDK, MMP11 and SALL4 were excluded from",
                      "the KEGG pathway visualisations because no direct",
                      "pathway associations were retrieved for their",
                      "mapped human orthologs in the KEGG database."
                    ),
                    style = "font-size:13px;line-height:1.7;color:#393232;"
                  ),

                  p(
                    paste(
                      "Although these four candidates had KEGG Orthology",
                      "(KO) annotations, no corresponding pathway mappings",
                      "were identified. A KO annotation does not necessarily",
                      "indicate an association with a specific KEGG pathway."
                    ),
                    style = "font-size:13px;line-height:1.7;color:#393232;"
                  ),

                  p(
                    paste(
                      "Therefore, only TGFB2 and VTN, which had retrieved",
                      "KEGG pathway associations, are included in the",
                      "visualisations. This does not imply that the other",
                      "four genes have no biological pathway involvement."
                    ),
                    style = "font-size:13px;line-height:1.7;color:#393232;margin-bottom:0;"
                  )
                ),

                layout_columns(

                  selectInput(
                    "kegg_candidates",
                    "Select candidate gene",
                    choices = c(
                      "TGFB2",
                      "VTN"
                    ),
                    selected = "TGFB2",
                    multiple = FALSE
                  ),

                  selectInput(
                    "kegg_graph",
                    "Select graph",
                    choices = c(
                      "Bar plot",
                      "Network plot"
                    ),
                    selected = "Bar plot"
                  ),

                  selectInput(
                    "kegg_relevance",
                    "Regeneration relevance",
                    choices = c("Direct", "Indirect"),
                    selected = "Direct",
                    multiple = FALSE
                  ),

                  col_widths = c(4, 4, 4)
                ),

                div(
                  class = "enrichment-plot",

                  plotlyOutput(
                    "kegg_plot",
                    height = "450px"
                  )
                )
              ),

              div(
                class = "focus-panel",

                h4("Interpretation note"),

                p(
                  paste(
                    "These visualisations display KEGG pathway annotation",
                    "records associated with the selected candidate gene.",
                    "They are not statistical enrichment significance plots:",
                    "the supplied annotation files do not contain p-values",
                    "or adjusted p-values."
                  )
                ),

                p(
                  paste(
                    "The KEGG annotations are based on mapped human",
                    "orthologs and should therefore be interpreted as",
                    "human-ortholog-based pathway annotations rather",
                    "than direct axolotl pathway annotations."
                  )
                )
              )
            )
          )
        )
      }


      # ------------------------------------------------------
      # GENE ONTOLOGY SECTION
      # ------------------------------------------------------

      if (selected_enrichment == "go") {

        return(
          tagList(

            actionButton(
              "back_enrichment",
              "Back to Functional Enrichment",
              icon = icon("arrow-left"),
              class = "btn-light"
            ),

            br(),
            br(),

            h2(
              icon("sitemap"),
              " Gene Ontology (GO)",
              style = paste(
                "font-weight:800;",
                "color:#242020;",
                "margin-bottom:18px;"
              )
            ),

            div(
              class = "enrichment-detail",

              div(
                class = "enrichment-panel",

                p(
                  "Explore biological process, cellular component and molecular function annotations.",
                  class = "enrichment-note"
                ),

                layout_columns(

                  selectInput(
                    "go_candidates",
                    "Select candidate gene",
                    choices = candidate_choices,
                    selected = candidate_choices[1],
                    multiple = FALSE
                  ),

                  selectInput(
                    "go_graph",
                    "Select graph",
                    choices = c(
                      "Bar plot",
                      "Gene treemap"
                    ),
                    selected = "Bar plot"
                  ),

                  selectInput(
                    "go_relevance",
                    "Regeneration relevance",
                    choices = c("Direct", "Indirect"),
                    selected = "Direct",
                    multiple = FALSE
                  ),

                  col_widths = c(4, 4, 4)
                ),

                div(
                  class = "enrichment-plot",

                  plotlyOutput(
                    "go_plot",
                    height = "600px"
                  )
                )
              ),

              div(
                class = "focus-panel",

                h4("Interpretation note"),

                p(
                  paste(
                    "These visualisations display annotation records associated",
                    "with the selected candidate gene. They are not statistical",
                    "enrichment significance plots: the supplied annotation files",
                    "do not contain p-values or adjusted p-values."
                  )
                )
              )
            )
          )
        )
      }


      # ------------------------------------------------------
      # FOCUSED PATHWAYS SECTION
      # ------------------------------------------------------

      if (selected_enrichment == "pathways") {

        return(
          tagList(

            actionButton(
              "back_enrichment",
              "Back to Functional Enrichment",
              icon = icon("arrow-left"),
              class = "btn-light"
            ),

            br(),
            br(),

            h2(
              icon("route"),
              " Pathways",
              style = paste(
                "font-weight:800;",
                "color:#242020;",
                "margin-bottom:18px;"
              )
            ),

            actionButton(
              "open_individual_pathways",
              "View Individual Pathways",
              icon = icon("list"),
              class = "btn-primary",
              style = "margin-bottom:18px;"
            ),

            div(
              class = "enrichment-detail",

              div(
                class = "enrichment-panel",

                p(
                  "Explore Reactome pathway associations across all six candidate genes.",
                  class = "enrichment-note"
                ),

                layout_columns(

                  selectInput(
                    "pathway_graph",
                    "Select graph",
                    choices = c("Heatmap", "Chord diagram", "Sunburst plot"),
                    selected = "Heatmap"
                  ),

                  col_widths = c(12)
                ),

                div(
                  class = "enrichment-plot",

                  plotlyOutput(
                    "pathway_plot",
                    height = "650px"
                  )
                )
              ),

              div(
                class = "focus-panel",

                h4("Interpretation note"),

                p(
                  paste(
                    "These visualisations display annotation records associated",
                    "with the six candidate genes. They are not statistical",
                    "enrichment significance plots: the supplied annotation files",
                    "do not contain p-values or adjusted p-values."
                  )
                ),

                p(
                  paste(
                    "The heatmap and chord diagram use the Reactome",
                    "gene-to-pathway associations. The sunburst follows the",
                    "Reactome pathway hierarchy and shows associated candidate genes."
                  )
                ),

              )
            )
          )
        )
      }


      # ------------------------------------------------------
      # INDIVIDUAL PATHWAYS TABLE
      # ------------------------------------------------------

      if (selected_enrichment == "individual_pathways") {

        return(
          tagList(

            actionButton(
              "back_to_pathways",
              "Back to Pathways",
              icon = icon("arrow-left"),
              class = "btn-light"
            ),

            br(),
            br(),

            h2(
              icon("list"),
              " Individual Pathways",
              style = paste(
                "font-weight:800;",
                "color:#242020;",
                "margin-bottom:18px;"
              )
            ),

            div(
              class = "enrichment-detail",

              div(
                class = "enrichment-panel",

                p(
                  "View the Reactome pathway records associated with an individual candidate gene and filter them by their top-level Reactome pathway.",
                  class = "enrichment-note"
                ),

                layout_columns(

                  selectInput(
                    "individual_pathway_candidate",
                    "Select candidate gene",
                    choices = candidate_choices,
                    selected = candidate_choices[1],
                    multiple = FALSE
                  ),

                  selectInput(
                    "individual_pathway_top_level",
                    "Top-level Reactome pathway",
                    choices = reactome_top_level_choices,
                    selected = "__ALL__",
                    multiple = FALSE
                  ),

                  col_widths = c(6, 6)
                ),

                div(
                  class = "enrichment-plot",
                  tableOutput("individual_pathway_table")
                )
              ),

              div(
                class = "focus-panel",

                h4("Interpretation note"),

                p(
                  paste(
                    "This table displays the gene-to-pathway associations",
                    "provided in Reactome_individual_gene_pathways.csv.",
                    "A candidate with no matching Reactome record in this",
                    "file will have no pathway rows displayed."
                  )
                ),

                p(
                  paste(
                    "The hierarchy shows the pathway's position within",
                    "Reactome's pathway organization. Top-level pathway",
                    "groups are organizational categories; they are not",
                    "Direct or Indirect regeneration-relevance labels, nor",
                    "statistical or causal evidence."
                  )
                )
              )
            )
          )
        )
      }
    }


    # ========================================================
    # GENE MAPPING PAGE
    # ========================================================

    if (page == "genes") {

      selected_gene_page <- gene_mapping_page()


      if (selected_gene_page == "overview") {

        return(
          tagList(

            h2(
              "Gene Mapping",
              style = paste(
                "font-weight:800;",
                "color:#242020;",
                "margin-bottom:22px;"
              )
            ),

            layout_columns(

              div(
                actionButton(
                  "open_gene_mapping_overview",

                  label = tagList(
                    div(
                      class = "card-symbol",
                      icon("chart-simple")
                    ),
                    h3("Overview"),
                    p(
                      "View the gene-mapping summary and the workflow."
                    ),
                    icon("arrow-right")
                  ),

                  class = "gene-mapping-button"
                )
              ),

              div(
                actionButton(
                  "open_gene_mapping_results",

                  label = tagList(
                    div(
                      class = "card-symbol",
                      icon("table")
                    ),
                    h3("Mapping Results"),
                    p(
                      "View the gene-level output from the human gene mapping step."
                    ),
                    icon("arrow-right")
                  ),

                  class = "gene-mapping-button"
                )
              ),

              col_widths = c(6, 6)
            )
          )
        )
      }


      # OVERVIEW PAGE: SUMMARY AND WORKFLOW TOGETHER

      if (selected_gene_page == "overview_details") {

        return(
          tagList(

            actionButton(
              "back_gene_mapping",
              "Back to Gene Mapping",
              icon = icon("arrow-left"),
              class = "btn-light"
            ),

            br(),
            br(),

            h2(
              "Overview",
              style = paste(
                "font-weight:800;",
                "color:#242020;",
                "margin-bottom:22px;"
              )
            ),

            layout_columns(

              div(
                class = "gene-mapping-section",

                div(
                  class = "card-symbol",
                  icon("file-import")
                ),

                h3("Input genes"),

                div(
                  class = "gene-mapping-number",
                  textOutput("mapping_input_count", inline = TRUE)
                ),

                p("Unique, non-common genes used as input for gene mapping.")
              ),

              div(
                class = "gene-mapping-section",

                div(
                  class = "card-symbol",
                  icon("dna")
                ),

                h3("Mapped human genes"),

                div(
                  class = "gene-mapping-number",
                  textOutput("mapping_human_record_count", inline = TRUE)
                ),

                p("Human genes obtained after the gene-mapping step.")
              ),

              col_widths = c(6, 6)
            ),

            br(),

            div(
              class = "gene-mapping-section",

              h3("Gene Mapping Workflow"),

              p(
                paste(
                  "The candidate-gene list and human gene comparison file are used to display the Axolotl-to-human gene records below."
                )
              ),

              div(
                class = "gene-mapping-flow",

                div(
                  class = "gene-mapping-flow-step",
                  strong("39"),
                  "Input genes"
                ),

                div(
                  class = "gene-mapping-flow-arrow",
                  icon("arrow-right")
                ),

                div(
                  class = "gene-mapping-flow-step",
                  strong("Mapping"),
                  "Human gene mapping"
                ),

                div(
                  class = "gene-mapping-flow-arrow",
                  icon("arrow-right")
                ),

                div(
                  class = "gene-mapping-flow-step",
                  strong("42"),
                  "Mapped human genes"
                )
              )
            )
          )
        )
      }


      # MAPPING RESULTS PAGE

      if (selected_gene_page == "mapping_results") {

        return(
          tagList(

            actionButton(
              "back_gene_mapping_overview",
              "Back to Gene Mapping",
              icon = icon("arrow-left"),
              class = "btn-light"
            ),

            br(),
            br(),

            h2(
              "Mapping Results",
              style = paste(
                "font-weight:800;",
                "color:#242020;",
                "margin-bottom:22px;"
              )
            ),

            p(
              "Explore the correspondence between the 39 Axolotl candidate genes and the human gene records returned by the comparison step.",
              style = "color:#59494A;line-height:1.7;"
            ),

            if (!is.null(mapping_data)) {

              tagList(

                layout_columns(

                  div(
                    class = "gene-mapping-section mapping-result-stat-card",

                    h4("Axolotl candidates"),

                    div(
                      class = "gene-mapping-number",
                      textOutput("mapping_input_count", inline = TRUE)
                    )
                  ),

                  div(
                    class = "gene-mapping-section mapping-result-stat-card",

                    h4("Human gene records"),

                    div(
                      class = "gene-mapping-number",
                      textOutput("mapping_human_record_count", inline = TRUE)
                    )
                  ),

                  col_widths = c(6, 6)
                ),

                br(),

                div(
                  class = "gene-mapping-section",

                  h4("Find and filter gene records"),

                  textInput(
                    "mapping_search",
                    "Search by Axolotl gene, gene ID, human gene symbol or description",
                    placeholder = "e.g. IER5L, 138..., MDK"
                  ),

                  layout_columns(

                    selectInput(
                      "mapping_status_filter",
                      "Human gene record status",
                      choices = c(
                        "All statuses",
                        "Human gene record found",
                        "No human gene record found"
                      ),
                      selected = "All statuses"
                    ),

                    selectInput(
                      "mapping_match_filter",
                      "Number of human records per Axolotl candidate",
                      choices = c(
                        "All match types",
                        "Single human record",
                        "Multiple human records",
                        "No human record"
                      ),
                      selected = "All match types"
                    ),

                    col_widths = c(6, 6)
                  ),

                  textOutput("mapping_filtered_count"),

                  div(
                    class = "mapping-results-scroll",
                    tableOutput("mapping_results_table")
                  ),

                  br(),

                  div(
                    style = "display:flex;gap:10px;flex-wrap:wrap;",

                    downloadButton(
                      "download_filtered_mapping",
                      "Download filtered results",
                      class = "btn-primary"
                    ),

                    downloadButton(
                      "download_all_mapping",
                      "Download all mapping results",
                      class = "btn-light"
                    )
                  ),

                  br(),

                  p(
                    "Interpretation note: 'Human gene record found' indicates that a human NCBI gene record was found. The supplied file marks orthology as unverified, so these matches should not be interpreted as confirmed orthologues.",
                    style = "color:#59494A;font-size:12px;line-height:1.7;"
                  )
                )
              )

            } else {

              div(
                class = "focus-panel",

                h4("Mapping files not available"),

                p(mapping_data_error)
              )
            }
          )
        )
      }
    }


    # ========================================================
    # RESULTS & INTERPRETATION
    # ========================================================

    if (page == "deg") {

      selected_results <- results_section()

      if (selected_results == "overview") {
        return(
          tagList(
            h2(
              "Results & Interpretation",
              style = "font-weight:800;color:#242020;margin-bottom:10px;"
            ),
            p(
              "Explore the differential-expression findings, visual summaries, gene overlaps and biological interpretation from the axolotl limb-regeneration RNA-seq analysis.",
              style = "color:#59494A;line-height:1.7;margin-bottom:22px;"
            ),
            div(
              class = "results-navigation-grid",
              actionButton(
                "results_open_overview",
                label = tagList(
                  div(class = "card-symbol", icon("chart-bar")),
                  h4("Results Overview"),
                  p("Review the overall findings and key regeneration-associated genes.")
                ),
                class = "workflow-card de-column comparison-option-button"
              ),
              actionButton(
                "results_open_deseq2",
                label = tagList(
                  div(class = "card-symbol", icon("chart-line")),
                  h4("DESeq2 Plots"),
                  p("Explore PCA, sample correlation and volcano plots.")
                ),
                class = "workflow-card de-column comparison-option-button"
              ),
              actionButton(
                "results_open_upset",
                label = tagList(
                  div(class = "card-symbol", icon("diagram-project")),
                  h4("UpSet Plot"),
                  p("Explore shared and unique differentially expressed genes.")
                ),
                class = "workflow-card de-column comparison-option-button"
              ),
              actionButton(
                "results_open_interpretation",
                label = tagList(
                  div(class = "card-symbol", icon("book-open")),
                  h4("Interpretation"),
                  p("Examine the biological relevance of axolotl–human gene matches.")
                ),
                class = "workflow-card de-column comparison-option-button"
              )
            )
          )
        )
      }

      if (selected_results == "overview_details") {
        # Build the regional 3 dpa vs 5 dpa summary from the same
        # Significant_DEGs.csv data used by the Comparison Table page.
        overview_regions <- c("Diaphysis", "Epiphysis")
        overview_comparison_rows <- lapply(overview_regions, function(region_name) {
          if (!is.null(comparison_data_error) || is.null(significant_degs)) {
            return(data.frame(
              region = region_name,
              total = NA_integer_,
              up = NA_integer_,
              down = NA_integer_,
              stringsAsFactors = FALSE
            ))
          }

          comparison_name <- paste0("5dpa_vs_3dpa_", region_name)
          region_rows <- significant_degs[
            significant_degs$comparison == comparison_name,
            ,
            drop = FALSE
          ]

          data.frame(
            region = region_name,
            total = nrow(region_rows),
            up = sum(region_rows$regulation == "Upregulated", na.rm = TRUE),
            down = sum(region_rows$regulation == "Downregulated", na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        })
        overview_comparison_summary <- do.call(rbind, overview_comparison_rows)

        overview_count_text <- function(value) {
          if (is.na(value)) "Not available" else format(value, big.mark = ",")
        }

        overview_table_body <- lapply(seq_len(nrow(overview_comparison_summary)), function(i) {
          row <- overview_comparison_summary[i, , drop = FALSE]
          tags$tr(
            tags$td(row$region),
            tags$td(overview_count_text(row$total)),
            tags$td(overview_count_text(row$up)),
            tags$td(overview_count_text(row$down))
          )
        })

        overview_annotation_count <- function(data, gene_name) {
          if (is.null(data) || !is.data.frame(data) ||
              !"candidate" %in% names(data)) {
            return(NA_integer_)
          }
          sum(trimws(as.character(data$candidate)) == gene_name, na.rm = TRUE)
        }

        overview_gene_rows <- lapply(candidate_choices, function(gene_name) {
          go_count <- overview_annotation_count(GO_data, gene_name)
          kegg_count <- overview_annotation_count(KEGG_data, gene_name)
          tags$tr(
            tags$td(tags$strong(gene_name)),
            tags$td(if (is.na(go_count)) "Not available" else format(go_count, big.mark = ",")),
            tags$td(if (is.na(kegg_count)) "Not available" else format(kegg_count, big.mark = ","))
          )
        })

        return(
          tagList(
            actionButton("results_back_home", "Back to Results & Interpretation", icon = icon("arrow-left"), class = "btn-light"),
            br(), br(),
            h2("Results Overview", style = "font-weight:800;color:#242020;margin-bottom:8px;"),
            p(
              "An integrated overview of the RNA-seq study, the regional differential-expression comparison, and the candidate-gene annotations being examined in the context of axolotl limb regeneration.",
              style = "color:#59494A;line-height:1.7;margin-bottom:22px;"
            ),

            # 1. OBJECTIVE AND STUDY DESIGN
            div(
              class = "workflow-card",
              h3("1. Objective of the RNA-seq Analysis", style = "font-weight:750;color:#302728;margin-bottom:10px;"),
              p(
                "The objective is to identify genes associated with axolotl limb regeneration and examine whether any of these genes are shared with humans. The broader research question is whether regeneration-associated genes in axolotl can provide useful molecular insights for studying tissue repair and regeneration in humans.",
                style = "color:#59494A;line-height:1.75;"
              ),
              div(
                style = "display:flex;flex-wrap:wrap;gap:10px;margin-top:14px;",
                div(style = "flex:1;min-width:135px;background:#F8F3F1;border-radius:12px;padding:13px;", strong("18"), br(), span("RNA-seq samples", style = "color:#6B5A5A;font-size:12px;")),
                div(style = "flex:1;min-width:135px;background:#F8F3F1;border-radius:12px;padding:13px;", strong("6"), br(), span("Experimental groups", style = "color:#6B5A5A;font-size:12px;")),
                div(style = "flex:1;min-width:135px;background:#F8F3F1;border-radius:12px;padding:13px;", strong("3"), br(), span("Biological replicates per group", style = "color:#6B5A5A;font-size:12px;")),
                div(style = "flex:1;min-width:135px;background:#F8F3F1;border-radius:12px;padding:13px;", strong("GSE334432"), br(), span("NCBI GEO dataset", style = "color:#6B5A5A;font-size:12px;"))
              ),
              p(
                "The six groups represent Intact, 3 days post-amputation (3 DPA), and 5 days post-amputation (5 DPA), examined in the diaphysis and epiphysis regions.",
                style = "color:#59494A;line-height:1.7;margin-top:14px;margin-bottom:0;"
              )
            ),
            br(),

            # 2. OVERALL DIFFERENTIAL-EXPRESSION RESULTS
            div(
              class = "workflow-card",
              h3("2. Overall Differential Expression Results", style = "font-weight:750;color:#302728;margin-bottom:10px;"),
              p(
                "The current significant-DEG comparison summarises the change from 3 DPA to 5 DPA separately for the diaphysis and epiphysis. A positive log2 fold change indicates higher expression at 5 DPA, while a negative value indicates higher expression at 3 DPA.",
                style = "color:#59494A;line-height:1.7;"
              ),
              if (!is.null(comparison_data_error)) {
                div(class = "alert alert-warning", comparison_data_error)
              } else {
                tags$table(
                  class = "table table-bordered table-striped",
                  style = "margin-top:12px;margin-bottom:8px;background:#FFFFFF;",
                  tags$thead(
                    tags$tr(
                      tags$th("Tissue region"),
                      tags$th("Significant DEGs"),
                      tags$th("Higher at 5 DPA"),
                      tags$th("Higher at 3 DPA")
                    )
                  ),
                  do.call(tags$tbody, overview_table_body)
                )
              },
              
            ),
            br(),

            # 3. REGIONAL COMPARISON
            div(
              class = "workflow-card",
              h3("3. Differences Between Diaphysis and Epiphysis", style = "font-weight:750;color:#302728;margin-bottom:10px;"),
              p(
                "The diaphysis and epiphysis are analysed as separate tissue regions. Comparing their 3 DPA versus 5 DPA results helps distinguish expression changes observed in both regions from changes that may be more region-specific.",
                style = "color:#59494A;line-height:1.75;"
              ),
              div(
                style = "display:flex;flex-wrap:wrap;gap:12px;margin-top:12px;",
                div(style = "flex:1;min-width:220px;border-left:4px solid #B98278;background:#FBF7F5;border-radius:8px;padding:13px 15px;", h4("Diaphysis", style = "font-weight:700;margin-top:0;"), p("Review the 3 DPA to 5 DPA expression changes in the shaft region of the regenerating limb.", style = "color:#59494A;line-height:1.6;margin-bottom:0;")),
                div(style = "flex:1;min-width:220px;border-left:4px solid #8D9F91;background:#F6F8F5;border-radius:8px;padding:13px 15px;", h4("Epiphysis", style = "font-weight:700;margin-top:0;"), p("Review the corresponding expression changes in the region near the limb ends.", style = "color:#59494A;line-height:1.6;margin-bottom:0;"))
              ),
              
            ),
            br(),

            # 4. FUNCTIONAL ANNOTATION AND CANDIDATE GENES
            div(
              class = "workflow-card",
              h3("4. Functional Annotation and Pathway Associations", style = "font-weight:750;color:#302728;margin-bottom:10px;"),
              p(
                "Functional annotation is used to explore the biological context of the candidate genes through Gene Ontology (GO), KEGG, and Reactome resources. The focused pathway review includes processes relevant to cell behaviour, tissue remodelling, and extracellular-matrix interactions.",
                style = "color:#59494A;line-height:1.7;"
              ),
              tags$table(
                class = "table table-bordered table-striped",
                style = "margin-top:12px;margin-bottom:8px;background:#FFFFFF;",
                tags$thead(tags$tr(tags$th("Candidate gene"), tags$th("GO records"), tags$th("KEGG records"))),
                do.call(tags$tbody, overview_gene_rows)
              ),
             
              
            ),
            br(),

            # 5. OVERALL CONCLUSION
            div(
              class = "focus-panel",
              h3("5. Overall Conclusion", style = "font-weight:750;color:#302728;margin-bottom:10px;"),
              p(
                "Together, the differential-expression comparisons and functional annotations provide a framework for investigating the molecular changes that accompany axolotl limb regeneration. The cross-species comparison then asks whether any of the regeneration-associated genes are also represented in the human gene set.",
                style = "color:#59494A;line-height:1.75;"
              ),
              p(
                "The central research direction is to identify shared genes and assess whether their known biological associations may help guide future studies of human tissue regeneration. Shared gene matches are a starting point for investigation, not evidence that the same regenerative process occurs in humans.",
                style = "color:#59494A;line-height:1.75;margin-bottom:0;"
              )
            )
          )
        )
      }

      if (selected_results == "deseq2") {
        return(
          tagList(
            actionButton("results_back_home", "Back to Results & Interpretation", icon = icon("arrow-left"), class = "btn-light"),
            br(), br(),
            h2("DESeq2 Plots", style = "font-weight:800;color:#242020;"),
            p("Select a plot to view its visualisation and corresponding results table.", style = "color:#59494A;line-height:1.7;"),
            div(
              class = "workflow-card",
              selectInput(
                "deseq_plot_type",
                "Select plot",
                choices = c(
                  "PCA" = "pca",
                  "Sample Correlation" = "correlation",
                  "Volcano Plot" = "volcano"
                ),
                selected = "pca",
                multiple = FALSE,
                width = "100%"
              )
            ),
            uiOutput("sample_analysis_status"),
            conditionalPanel(
              condition = "input.deseq_plot_type == 'pca'",
              div(
                class = "focus-panel",
                h4("Principal Component Analysis (PCA)"),
                p("PCA displays the supplied PC1 and PC2 coordinates for all 18 samples. Hover over a point to inspect its sample ID, experimental group and PC1–PC3 scores.", style = "color:#59494A;line-height:1.7;"),
                plotlyOutput("sample_pca_plot", height = "500px"),
                h4("PCA Coordinates Table", style = "margin-top:20px;"),
                tableOutput("sample_pca_table")
              )
            ),
            conditionalPanel(
              condition = "input.deseq_plot_type == 'correlation'",
              div(
                class = "focus-panel",
                h4("Sample Correlation Heatmap"),
                p("The supplied Pearson correlation matrix for all 18 samples. Hover over any cell to inspect the sample pair and correlation coefficient.", style = "color:#59494A;line-height:1.7;"),
                plotlyOutput("sample_correlation_plot", height = "620px"),
                h4("Sample Correlation Matrix", style = "margin-top:20px;"),
                tableOutput("sample_correlation_table")
              )
            ),
            conditionalPanel(
              condition = "input.deseq_plot_type == 'volcano'",
              div(
                class = "focus-panel",
                h4("Volcano Plot"),
                p("Select a comparison to explore the supplied DESeq2 results for all tested genes, including non-significant genes. Hover over points for gene-level statistics.", style = "color:#59494A;line-height:1.7;"),
                selectInput(
                  "volcano_comparison",
                  "Select comparison",
                  choices = c(
                    "3 dpa Diaphysis vs Intact Diaphysis" = "3dpa_Diaphysis_vs_Intact_Diaphysis",
                    "3 dpa Epiphysis vs Intact Epiphysis" = "3dpa_Epiphysis_vs_Intact_Epiphysis",
                    "5 dpa Diaphysis vs Intact Diaphysis" = "5dpa_Diaphysis_vs_Intact_Diaphysis",
                    "5 dpa Epiphysis vs Intact Epiphysis" = "5dpa_Epiphysis_vs_Intact_Epiphysis",
                    "5 dpa vs 3 dpa Diaphysis" = "5dpa_vs_3dpa_Diaphysis",
                    "5 dpa vs 3 dpa Epiphysis" = "5dpa_vs_3dpa_Epiphysis"
                  ),
                  selected = "3dpa_Diaphysis_vs_Intact_Diaphysis"
                ),
                plotlyOutput("deseq_volcano_plot", height = "520px"),
                h4("Differential Expression Results Table", style = "margin-top:20px;"),
                tableOutput("deseq_volcano_table"),
                p("The interactive volcano plot is displayed from the supplied DESeq2 comparison-result CSV. Dashed lines mark log2 fold change = ±1 and adjusted p-value = 0.05.", style = "font-size:12px;color:#6A5E5E;line-height:1.6;margin-top:10px;")
              )
            )
          )
        )
      }

      if (selected_results == "upset") {
        return(
          tagList(
            actionButton("results_back_home", "Back to Results & Interpretation", icon = icon("arrow-left"), class = "btn-light"),
            br(), br(),
            h2("UpSet Plot", style = "font-weight:800;color:#242020;"),
            p(
              "Compare shared and comparison-specific significant DEGs across the six regeneration contrasts. Genes are matched using the exact gene_id in Significant_DEGs.csv.",
              style = "color:#59494A;line-height:1.7;"
            ),
            div(
              class = "focus-panel",
              p(
                "Inclusion criteria: adjusted p-value < 0.05 and absolute log2 fold change >= 1. The plot displays the 20 largest observed intersections. Click an active membership-matrix dot to inspect genes only when that intersection is shared by two or more comparisons. Dots belonging to unique, single-comparison intersections do nothing. Hover over bars and dots for details; use the selector to inspect any observed intersection.",
                style = "color:#59494A;line-height:1.7;"
              ),
              plotlyOutput("upset_plot", height = "760px"),
              h4("Significant DEG set sizes", style = "margin-top:22px;"),
              tableOutput("upset_set_sizes_table"),
              tags$hr(),
              h4("Inspect genes in an intersection", style = "margin-top:20px;"),
              uiOutput("upset_intersection_selector"),
              div(
                style = "display:flex;justify-content:flex-end;margin:8px 0 12px;",
                downloadButton(
                  "download_upset_genes",
                  "Download all genes in this intersection",
                  class = "btn-primary"
                )
              ),
              p(
                "The table previews the first 100 gene identifiers. The download contains all genes in the selected intersection, with their comparison-specific DESeq2 results.",
                style = "font-size:12px;color:#6A5E5E;line-height:1.6;"
              ),
              div(
                style = "max-height:420px;overflow-y:auto;",
                tableOutput("upset_intersection_table")
              )
            )
          )
        )
      }

      if (selected_results == "interpretation") {

        # ======================================================
        # INTERPRETATION: 20 DATA-LINKED QUESTIONS
        # All calculations below use the supplied DEG and
        # candidate annotation files already loaded by the app.
        # ======================================================

        interpretation_comparisons <- c(
          "3dpa_Diaphysis_vs_Intact_Diaphysis",
          "3dpa_Epiphysis_vs_Intact_Epiphysis",
          "5dpa_Diaphysis_vs_Intact_Diaphysis",
          "5dpa_Epiphysis_vs_Intact_Epiphysis",
          "5dpa_vs_3dpa_Diaphysis",
          "5dpa_vs_3dpa_Epiphysis"
        )

        interpretation_labels <- c(
          "3 dpa Diaphysis vs Intact Diaphysis",
          "3 dpa Epiphysis vs Intact Epiphysis",
          "5 dpa Diaphysis vs Intact Diaphysis",
          "5 dpa Epiphysis vs Intact Epiphysis",
          "5 dpa vs 3 dpa Diaphysis",
          "5 dpa vs 3 dpa Epiphysis"
        )
        names(interpretation_labels) <- interpretation_comparisons

        interpretation_degs <- data.frame()
        if (!is.null(significant_degs) &&
            all(c("gene_id", "comparison", "log2FoldChange", "padj") %in% names(significant_degs))) {
          interpretation_degs <- significant_degs[
            significant_degs$comparison %in% interpretation_comparisons &
              !is.na(significant_degs$padj) & significant_degs$padj < 0.05 &
              !is.na(significant_degs$log2FoldChange) &
              abs(significant_degs$log2FoldChange) >= 1,
            , drop = FALSE
          ]
        }

        fmt_num <- function(x, digits = 3) {
          ifelse(is.na(x), "—", formatC(x, format = "f", digits = digits))
        }

        fmt_padj <- function(x) {
          ifelse(is.na(x), "—", formatC(x, format = "e", digits = 2))
        }

        make_interpretation_table <- function(x) {
          if (is.null(x) || nrow(x) == 0) {
            return(tags$p("No matching records are available in the supplied files.", class = "interpretation-note"))
          }
          header <- tags$thead(tags$tr(lapply(names(x), function(nm) tags$th(nm))))
          body <- tags$tbody(lapply(seq_len(nrow(x)), function(i) {
            tags$tr(lapply(x[i, , drop = FALSE], function(value) {
              value <- as.character(value[[1]])
              if (is.na(value) || !nzchar(value)) value <- "—"
              tags$td(value)
            }))
          }))
          tags$div(
            class = "interpretation-table-wrap",
            tags$table(class = "table table-sm table-striped table-hover interpretation-table", header, body)
          )
        }

        make_question <- function(question, answer) {
          tags$details(
            class = "interpretation-question",
            tags$summary(question),
            tags$div(class = "interpretation-answer", answer)
          )
        }

        make_group <- function(heading, questions) {
          tags$section(
            class = "interpretation-group",
            tags$h3(heading),
            tags$div(class = "interpretation-question-list", questions)
          )
        }

        comparison_counts <- do.call(rbind, lapply(interpretation_comparisons, function(comp) {
          x <- interpretation_degs[interpretation_degs$comparison == comp, , drop = FALSE]
          data.frame(
            Comparison = interpretation_labels[[comp]],
            `Significant DEGs` = nrow(x),
            Upregulated = sum(x$log2FoldChange >= 1, na.rm = TRUE),
            Downregulated = sum(x$log2FoldChange <= -1, na.rm = TRUE),
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
        }))

        stage_counts <- comparison_counts[match(
          c("3 dpa Diaphysis vs Intact Diaphysis", "3 dpa Epiphysis vs Intact Epiphysis",
            "5 dpa Diaphysis vs Intact Diaphysis", "5 dpa Epiphysis vs Intact Epiphysis"),
          comparison_counts$Comparison
        ), , drop = FALSE]
        stage_counts$Stage <- c("3 dpa", "3 dpa", "5 dpa", "5 dpa")
        stage_totals <- aggregate(
          cbind(`Significant DEGs`, Upregulated, Downregulated) ~ Stage,
          data = stage_counts,
          FUN = sum
        )
        stage_totals <- stage_totals[, c("Stage", "Significant DEGs", "Upregulated", "Downregulated"), drop = FALSE]

        top_gene_table <- function(direction = c("up", "down")) {
          direction <- match.arg(direction)
          pieces <- lapply(interpretation_comparisons, function(comp) {
            x <- interpretation_degs[interpretation_degs$comparison == comp, , drop = FALSE]
            if (direction == "up") {
              x <- x[x$log2FoldChange >= 1, , drop = FALSE]
              x <- x[order(-x$log2FoldChange, x$padj), , drop = FALSE]
            } else {
              x <- x[x$log2FoldChange <= -1, , drop = FALSE]
              x <- x[order(x$log2FoldChange, x$padj), , drop = FALSE]
            }
            x <- head(x, 10)
            if (nrow(x) == 0) return(NULL)
            data.frame(
              Comparison = interpretation_labels[[comp]],
              Rank = seq_len(nrow(x)),
              Gene = as.character(x$gene_id),
              `log2 fold change` = fmt_num(x$log2FoldChange, 3),
              `Adjusted p-value` = fmt_padj(x$padj),
              check.names = FALSE,
              stringsAsFactors = FALSE
            )
          })
          pieces <- Filter(Negate(is.null), pieces)
          if (length(pieces) == 0) return(data.frame())
          do.call(rbind, pieces)
        }

        top_up_table <- top_gene_table("up")
        top_down_table <- top_gene_table("down")

        # GO and KEGG annotations in focused_unique_features.csv are
        # candidate-gene annotation records, not statistical enrichment
        # results from the complete DEG lists.
        focused_annotations <- if (!is.null(focused_data) &&
          all(c("annotation_type", "category", "annotation", "genes") %in% names(focused_data))) {
          focused_data
        } else {
          data.frame()
        }
        go_annotations <- if (nrow(focused_annotations) > 0) {
          focused_annotations[toupper(focused_annotations$annotation_type) == "GO", , drop = FALSE]
        } else data.frame()
        kegg_annotations <- if (nrow(focused_annotations) > 0) {
          focused_annotations[toupper(focused_annotations$annotation_type) == "KEGG", , drop = FALSE]
        } else data.frame()

        reactome_records <- if (!is.null(Reactome_individual_data) &&
          all(c("gene", "pathway_id", "pathway_name", "top_level_pathway") %in% names(Reactome_individual_data))) {
          unique(Reactome_individual_data[, c("gene", "pathway_id", "pathway_name", "top_level_pathway"), drop = FALSE])
        } else data.frame()

        reactome_summary <- data.frame()
        reactome_level_summary <- data.frame()
        if (nrow(reactome_records) > 0) {
          reactome_summary <- aggregate(
            gene ~ pathway_id + pathway_name + top_level_pathway,
            data = reactome_records,
            FUN = function(x) paste(sort(unique(as.character(x))), collapse = "; ")
          )
          names(reactome_summary) <- c("Pathway ID", "Reactome pathway", "Top-level category", "Associated candidates")
          reactome_level_summary <- as.data.frame(table(reactome_records$top_level_pathway), stringsAsFactors = FALSE)
          names(reactome_level_summary) <- c("Top-level category", "Pathway records")
          reactome_level_summary <- reactome_level_summary[order(-reactome_level_summary$`Pathway records`), , drop = FALSE]
        }

        candidate_expression <- data.frame()
        candidate_expression_source <- data.frame()
        if (!is.null(significant_degs) &&
            all(c("gene_id", "comparison", "log2FoldChange", "padj") %in% names(significant_degs))) {
          candidate_expression_source <- significant_degs[
            significant_degs$comparison %in% interpretation_comparisons &
              !is.na(significant_degs$padj) & significant_degs$padj < 0.05 &
              !is.na(significant_degs$log2FoldChange),
            , drop = FALSE
          ]
        }
        if (nrow(candidate_expression_source) > 0) {
          candidate_expression_rows <- lapply(candidate_choices, function(gene_name) {
            gene_pattern <- paste0("\\|\\s*", gene_name, "\\s*$")
            x <- candidate_expression_source[grepl(gene_pattern, trimws(as.character(candidate_expression_source$gene_id))), , drop = FALSE]
            if (nrow(x) == 0) return(NULL)
            data.frame(
              Candidate = gene_name,
              Comparison = unname(interpretation_labels[x$comparison]),
              `log2 fold change` = fmt_num(x$log2FoldChange, 3),
              `Adjusted p-value` = fmt_padj(x$padj),
              Direction = ifelse(x$log2FoldChange > 0, "Higher in first-named group", "Lower in first-named group"),
              check.names = FALSE,
              stringsAsFactors = FALSE
            )
          })
          candidate_expression_rows <- Filter(Negate(is.null), candidate_expression_rows)
          if (length(candidate_expression_rows) > 0) candidate_expression <- do.call(rbind, candidate_expression_rows)
        }

        candidate_function_table <- data.frame()
        if (nrow(focused_annotations) > 0) {
          candidate_function_table <- do.call(rbind, lapply(candidate_choices, function(gene_name) {
            go_terms <- if (nrow(go_annotations) > 0) {
              as.character(go_annotations$annotation[grepl(paste0("(^|;\\s*)", gene_name, "($|;)") , as.character(go_annotations$genes))])
            } else character(0)
            kegg_terms <- if (nrow(kegg_annotations) > 0) {
              as.character(kegg_annotations$annotation[grepl(paste0("(^|;\\s*)", gene_name, "($|;)") , as.character(kegg_annotations$genes))])
            } else character(0)
            reactome_terms <- if (nrow(reactome_records) > 0) {
              as.character(reactome_records$pathway_name[reactome_records$gene == gene_name])
            } else character(0)
            data.frame(
              Candidate = gene_name,
              `GO annotations` = if (length(go_terms)) paste(unique(go_terms), collapse = "; ") else "No GO record in focused file",
              `KEGG annotations` = if (length(kegg_terms)) paste(unique(kegg_terms), collapse = "; ") else "No KEGG record in focused file",
              `Reactome associations` = if (length(reactome_terms)) paste(unique(reactome_terms), collapse = "; ") else "No Reactome record",
              check.names = FALSE,
              stringsAsFactors = FALSE
            )
          }))
        }

        human_candidate_table <- data.frame()
        if (!is.null(human_mapping) &&
            all(c("candidate", "human_gene_symbol", "human_gene_description", "search_status") %in% names(human_mapping))) {
          human_candidate_table <- do.call(rbind, lapply(candidate_choices, function(gene_name) {
            x <- human_mapping[
              trimws(as.character(human_mapping$candidate)) == gene_name &
                trimws(as.character(human_mapping$human_gene_symbol)) == gene_name,
              , drop = FALSE
            ]
            if (nrow(x) == 0) {
              return(data.frame(Candidate = gene_name, `Human gene record` = "No exact-symbol match in supplied mapping", `Evidence status` = "Not established", check.names = FALSE))
            }
            data.frame(
              Candidate = gene_name,
              `Human gene record` = paste0(as.character(x$human_gene_symbol[1]), " — ", as.character(x$human_gene_description[1])),
              `Evidence status` = as.character(x$search_status[1]),
              check.names = FALSE,
              stringsAsFactors = FALSE
            )
          }))
        }

        go_bp <- if (nrow(go_annotations) > 0 && "category" %in% names(go_annotations)) {
          go_annotations[toupper(as.character(go_annotations$category)) == "BP", , drop = FALSE]
        } else data.frame()
        go_cc <- if (nrow(go_annotations) > 0 && "category" %in% names(go_annotations)) {
          go_annotations[toupper(as.character(go_annotations$category)) == "CC", , drop = FALSE]
        } else data.frame()
        go_mf <- if (nrow(go_annotations) > 0 && "category" %in% names(go_annotations)) {
          go_annotations[toupper(as.character(go_annotations$category)) == "MF", , drop = FALSE]
        } else data.frame()

        q1 <- make_question(
          "1. How does gene expression change between 3 dpa and 5 dpa in the diaphysis?",
          tagList(
            p("For the direct 5 dpa vs 3 dpa Diaphysis comparison, positive log2 fold changes indicate higher expression at 5 dpa; negative values indicate higher expression at 3 dpa."),
            make_interpretation_table(comparison_counts[comparison_counts$Comparison == interpretation_labels[["5dpa_vs_3dpa_Diaphysis"]], , drop = FALSE]),
            p("With the analysis cut-offs of adjusted p-value < 0.05 and |log2 fold change| >= 1, more genes are lower at 5 dpa than at 3 dpa in this tissue. This describes the direction of the detected changes; it does not by itself identify the biological processes causing them.", class = "interpretation-note")
          )
        )

        q2 <- make_question(
          "2. How does gene expression change between 3 dpa and 5 dpa in the epiphysis?",
          tagList(
            p("For the direct 5 dpa vs 3 dpa Epiphysis comparison, positive log2 fold changes indicate higher expression at 5 dpa and negative values indicate higher expression at 3 dpa."),
            make_interpretation_table(comparison_counts[comparison_counts$Comparison == interpretation_labels[["5dpa_vs_3dpa_Epiphysis"]], , drop = FALSE]),
            p("The downregulated 5 dpa group is larger than the upregulated 5 dpa group under the stated cut offs. This indicates a broad shift in transcript abundance between the sampled stages, not necessarily a simple increase or decrease in regeneration activity.", class = "interpretation-note")
          )
        )

        q3 <- make_question(
          "3. Which regeneration stage shows more differentially expressed genes?",
          tagList(
            p("Comparing each stage against its matching intact tissue, the combined number of thresholded DEG records is higher at 5 dpa than at 3 dpa. The tissue-specific counts are shown separately because the response differs between regions."),
            make_interpretation_table(stage_counts),
            make_interpretation_table(stage_totals),
            p("Across both tissues, 3 dpa has 10,145 DEG records and 5 dpa has 10,908, a difference of 763 records. However, diaphysis has slightly more DEGs at 3 dpa, whereas epiphysis has more at 5 dpa. These are comparison-specific DEG records, not necessarily unique genes across tissues.", class = "interpretation-note")
          )
        )

        q4 <- make_question(
          "4. What do these changes suggest about the progression of limb regeneration?",
          tagList(
            p("The direct stage to stage contrasts identify 958 DEGs in the diaphysis and 820 in the epiphysis using adjusted p-value < 0.05 and (log2 fold change) >= 1."),
            make_interpretation_table(comparison_counts[comparison_counts$Comparison %in% c(interpretation_labels[["5dpa_vs_3dpa_Diaphysis"]], interpretation_labels[["5dpa_vs_3dpa_Epiphysis"]]), , drop = FALSE]),
            p("In both tissues, more genes have negative than positive log2 fold changes in the 5 dpa vs 3 dpa contrast. Thus, among genes meeting the thresholds, more are expressed lower at 5 dpa than at 3 dpa. The data support stage associated transcriptional.", class = "interpretation-note")
          )
        )

        q5 <- make_question(
          "5. What are the top 10 upregulated genes in each comparison?",
          tagList(
            p("The table lists the ten largest positive log2 fold changes in each contrast among genes meeting adjusted p-value < 0.05 and (log2 fold change) >= 1. The comparison label defines the reference direction."),
            p("Many identifiers are transcript assembly IDs (MSTRG) or LOC identifiers. Their rank reflects fold change magnitude, not proven biological importance.", class = "interpretation-note")
          )
        )

        q6 <- make_question(
          "6. What are the top 10 downregulated genes in each comparison?",
          tagList(
            p("The table lists the ten most negative log2 fold changes in each contrast among genes meeting adjusted p-value < 0.05 and (log2 fold change) >= 1."),
            p("For the 5 dpa vs 3 dpa contrasts, negative values mean lower expression at 5 dpa. For the stage versus intact contrasts, negative values mean lower expression in the named post-amputation tissue relative to its intact reference.", class = "interpretation-note")
          )
        )

        q7 <- make_question(
          "7. What are the potential biological functions of these highly differentially expressed genes?",
          tagList(
            p("The supplied top-gene tables contain a mixture of named gene symbols, LOC identifiers and StringTie-predicted transcript identifiers"),
            p("A defensible functional interpretation requires gene annotation (for example, verified gene symbols, protein descriptions, orthology support or GO/Pathway records) for each top DEG. The current focused annotation file covers the six selected candidate genes, not every top DEG. Therefore, this app does not invent functions for uncharacterised MSTRG/LOC entries.", class = "interpretation-note")
          )
        )

        q8 <- make_question(
          "8. Are any of the top genes associated with known processes involved in tissue regeneration?",
          tagList(
            p("The six selected candidate genes have supplied functional annotations linked to tissue repair, cell migration, extracellular-matrix remodelling, wound healing and stem-cell or proliferation related processes. Their individual expression patterns are summarised in the candidate-gene section below."),
            
          )
        )

        q9 <- make_question(
          "9. What are the major biological processes represented by the GO annotations?",
          tagList(
            p(paste0("The focused candidate annotation file contains ", nrow(go_bp), " Biological Process (BP) GO records.")),
            make_interpretation_table(if (nrow(go_bp) > 0) go_bp[, intersect(c("annotation_id", "annotation", "genes"), names(go_bp)), drop = FALSE] else go_bp),
            p("The listed terms include tissue regeneration, positive regulation of cell migration, extracellular-matrix disassembly, collagen catabolism, stem-cell population maintenance, regulation of cell proliferation and wound healing. These are candidate-level functional annotations, not statistically enriched GO terms from the complete DEG set.", class = "interpretation-note")
          )
        )

        q10 <- make_question(
          "10. Which cellular components and molecular functions are represented?",
          tagList(
            p(paste0("The supplied focused GO records include ", nrow(go_cc), " Cellular Component and ", nrow(go_mf), " Molecular Function.")),
            p("The available Cellular Components record is focal adhesion (FLRT3). No MF records are found. This should not be interpreted as evidence that molecular functions are absent; it means they are not represented in this supplied focused file.", class = "interpretation-note")
          )
        )

        q11 <- make_question(
          "11. What do these GO annotations suggest about the molecular mechanisms involved in regeneration?",
          tagList(
            p("Taken together, the candidate annotations point to several biological themes: cell migration and tissue regeneration (MDK), extracellular-matrix and collagen breakdown (MMP11), stem-cell population maintenance (SALL4), cell-population proliferation and wound healing (TGFB2), and wound-healing or focal adhesion related annotations (VTN and FLRT3)."),
            p("These annotations provide hypotheses about candidate gene functions.", class = "interpretation-note")
          )
        )

        q12 <- make_question(
          "12. Which major pathways are associated with the candidate regeneration genes?",
          tagList(
            p(paste0("The Reactome individual-gene file contains ", nrow(reactome_records), " gene to pathway records across ", if (nrow(reactome_records) > 0) length(unique(reactome_records$pathway_id)) else 0, " distinct pathways.")),
            
            p("The annotation file also lists these KEGG associations: MAPK signalling and TGF-beta signalling (TGFB2), and ECM-receptor interaction, focal adhesion and PI3K-Akt signalling (VTN). These are gene-pathway annotation matches rather than pathway-enrichment statistics.", class = "interpretation-note"),
            make_interpretation_table(if (nrow(kegg_annotations) > 0) kegg_annotations[, intersect(c("annotation_id", "annotation", "genes"), names(kegg_annotations)), drop = FALSE] else kegg_annotations)
          )
        )

        q13 <- make_question(
          "13. What are the potential roles of these pathways in tissue repair and regeneration?",
          tagList(
            p("The supplied Reactome records connect the candidates to extracellular-matrix organisation and degradation, receptor-mediated signal transduction, developmental biology and stem-cell transcriptional programmes."),
            p("Examples include collagen degradation and activation of matrix metalloproteinases (MMP11); TGF-beta receptor signalling and SMAD activation (TGFB2); ALK-related signalling (MDK); FGFR and ROBO signalling (FLRT3); and pluripotency-associated transcriptional pathways (SALL4). These pathway associations suggest plausible mechanisms to investigate, but do not demonstrate pathway activation in the experimental tissue.", class = "interpretation-note")
          )
        )

        q14 <- make_question(
          "14. Which pathways are associated with cell proliferation, signalling, extracellular-matrix remodelling and tissue development?",
          tagList(
            p("The Reactome hierarchy groups have candidate associations into broad biological categories. The table counts distinct pathway records per top-level category."),
            make_interpretation_table(reactome_level_summary),
            p("Extracellular-matrix organisation includes collagen degradation, matrix metalloproteinase activation, proteoglycans and integrin-related interactions. Signal transduction includes ALK, FGFR and TGF-beta-related signalling. Developmental Biology includes ROBO signalling and stem-cell or pluripotency-associated records. GO/KEGG annotations additionally include cell-population proliferation, wound healing, MAPK, PI3K-Akt and focal adhesion.", class = "interpretation-note")
          )
        )

        q15 <- make_question(
          "15. Are there differences in enriched pathways between 3 dpa and 5 dpa or between the two tissues?",
          tagList(
            p("The supplied Reactome and focused GO and KEGG files are candidate gene annotation datasets. They do not contain a comparison, timepoint or tissue column that links each pathway to a specific 3 dpa and 5 dpa or diaphysis and epiphysis contrast."),
            
          )
        )

        q16 <- make_question(
          "16. What are the expression patterns of FLRT3, MDK, MMP11, SALL4, TGFB2 and VTN?",
          tagList(
            
            p("The candidate genes show distinct patterns: FLRT3 and TGFB2 have negative fold changes in the regeneration vs intact comparisons; MDK, MMP11, SALL4 and VTN have positive fold changes in most or all of those comparisons. The direct 5 dpa vs 3 dpa contrasts provide the stage specific directionn.", class = "interpretation-note")
          )
        )

        q17 <- make_question(
          "17. What are the known biological functions of these candidate genes?",
          tagList(
            make_interpretation_table(candidate_function_table),
            p("The functions shown are based only on the supplied focused GO, KEGG and Reactome annotations. They should be treated as annotation supported functional context, not as experimentally validated roles in axolotl limb regeneration.", class = "interpretation-note")
          )
        )

        q18 <- make_question(
          "18. How might these candidate genes contribute to tissue regeneration?",
          tagList(
            p("The annotation records suggest complementary hypotheses: MDK may be relevant to cell migration and tissue regeneration signalling; MMP11 to extracellular-matrix and collagen remodelling; SALL4 to stem-cell population maintenance and developmental transcription; TGFB2 to TGF-beta signalling, proliferation regulation and wound healing; VTN to wound healing and cell-matrix interactions; and FLRT3 to focal adhesion and developmental signalling."),
            
          )
        )

        q19 <- make_question(
          "19. What do the combined differential-expression and pathway results reveal about axolotl limb regeneration?",
          tagList(
            p("The DEG comparisons show extensive transcriptional differences between intact tissue and the 3 dpa and 5 dpa samples, as well as hundreds of genes changing directly between 3 dpa and 5 dpa in each tissue. The direction and magnitude of these changes differ between diaphysis and epiphysis."),
            p("The six selected candidates connect the expression results to annotations involving cell migration, tissue repair, extracellular-matrix remodelling, wound healing, proliferation-related regulation and developmental signalling. Together, these findings support the view that regeneration is associated with coordinated changes across multiple biological systems.", class = "interpretation-note")
          )
        )

        q20 <- make_question(
          "20. Which findings could provide insights into potential mechanisms of human tissue regeneration?",
          tagList(
            p("The six candidate genes have matching human gene symbol records in the supplied comparison file. The file explicitly labels these records as 'Human gene record found; orthology unverified'. A matching symbol or gene record must not be presented as a confirmed one-to-one orthologue."),
            make_interpretation_table(human_candidate_table),
            p("The most useful translational leads are therefore hypotheses around conserved biological themes such as extracellular-matrix remodelling, TGF-beta, MAPK, PI3K-Akt signalling, wound healing, cell migration and stem-cell regulation. The current analysis does not show that these genes will induce regeneration in humans. Orthology verification, comparison with human tissue-repair datasets and experimental validation are required.", class = "interpretation-note")
          )
        )

        return(
          tagList(
            actionButton("results_back_home", "Back to Results & Interpretation", icon = icon("arrow-left"), class = "btn-light"),
            br(), br(),
            h2("Interpretation", style = "font-weight:800;color:#242020;"),
            p("Explore the biological interpretation of the axolotl RNA-seq results. Expand a question to view its answer, supporting results and interpretation.", class = "workflow-intro"),
            
            tags$style(HTML("\n              .interpretation-group { margin: 24px 0 30px; }\n              .interpretation-group h3 { color:#3B3031; font-size:20px; font-weight:800; padding:0 0 10px; border-bottom:1px solid #E7DAD6; margin-bottom:12px; }\n              .interpretation-question { background:#FFFCFA; border:1px solid #E8DCD7; border-radius:12px; margin:10px 0; overflow:hidden; }\n              .interpretation-question summary { cursor:pointer; list-style:none; padding:15px 48px 15px 17px; position:relative; color:#382E2F; font-weight:700; line-height:1.5; }\n              .interpretation-question summary::-webkit-details-marker { display:none; }\n              .interpretation-question summary:after { content:'\\25BC'; position:absolute; right:18px; top:15px; color:#9A6D73; font-size:13px; transition:transform .15s ease; }\n              .interpretation-question[open] summary:after { transform:rotate(180deg); }\n              .interpretation-question[open] summary { border-bottom:1px solid #EFE5E1; background:#FBF4F1; }\n              .interpretation-answer { padding:16px 18px 18px; color:#514748; line-height:1.7; }\n              .interpretation-answer p { margin:0 0 12px; }\n              .interpretation-note { background:#F7F0EC; border-left:3px solid #C28F99; padding:10px 13px; border-radius:4px; margin-top:12px !important; }\n              .interpretation-table-wrap { width:100%; overflow-x:auto; margin:12px 0 16px; }\n              .interpretation-table { min-width:620px; margin-bottom:0; font-size:13px; }\n              .interpretation-table th { white-space:nowrap; color:#493B3D; background:#F4EAE6; }\n              .interpretation-table td { vertical-align:top; }\n            ")),
            make_group("1. Comparison of 3 dpa and 5 dpa", tagList(q1, q2, q3, q4)),
            make_group("2. Top differentially expressed genes", tagList(q5, q6, q7, q8)),
            make_group("3. GO functional interpretation", tagList(q9, q10, q11)),
            make_group("4. KEGG and Reactome pathway interpretation", tagList(q12, q13, q14, q15)),
            make_group("5. Candidate gene interpretation", tagList(q16, q17, q18)),
            make_group("6. Overall biological interpretation", tagList(q19, q20))
          )
        )
      }
    }


    # PLACEHOLDER FOR OTHER SECTIONS

    if (page != "home") {

      return(
        tagList(

          actionButton(
            "back_home",
            "Back to Home",
            icon = icon("arrow-left"),
            class = "btn-light"
          ),

          br(),
          br(),

          div(
            class = "focus-panel",

            h3("This section is coming next"),

            p(
              paste(
                "The homepage is ready. The data tables,",
                "differential-expression results and plots",
                "will be connected in the next stage."
              )
            )
          )
        )
      )
    }


    # HOME PAGE

    tagList(

      div(
        class = "hero",

        p(
          class = "eyebrow",
          "AXOLOTL LIMB REGENERATION · RNA-seq"
        ),

        h1("From amputation to regeneration"),

        br(),

        p(
          paste(
            "The axolotl (Ambystoma mexicanum) can regenerate lost limbs,",
            "making it a powerful model to study the molecular and cellular",
            "mechanisms of regeneration."
          )
        )
      ),

      # AXOLOTL LIMB REGENERATION STAGES
      tags$style(HTML("
        @import url('https://fonts.googleapis.com/css2?family=Caveat:wght@500;600;700&display=swap');

        .home-regeneration-stages {
          display: grid;
          grid-template-columns: repeat(3, minmax(0, 1fr));
          gap: 20px;
          width: 100%;
          margin: 18px 0 26px;
          align-items: start;
        }
        .home-regeneration-stage {
          position: relative;
          background: #FFFDFC;
          border: 1px solid #E8D5D2;
          border-radius: 16px;
          overflow: hidden;
          box-shadow: 0 3px 12px rgba(89, 73, 74, .06);
          transition: border-color .25s ease, box-shadow .25s ease;
        }
        .home-regeneration-stage summary {
          display: block;
          cursor: pointer;
          list-style: none;
          padding: 12px;
          color: #242020;
          font-weight: 800;
        }
        .home-regeneration-stage summary::-webkit-details-marker {
          display: none;
        }
        .home-regeneration-stage summary:focus-visible {
          outline: 3px solid #C28F99;
          outline-offset: -3px;
        }
        .home-regeneration-stage summary img {
          display: block;
          width: 100%;
          height: 220px;
          object-fit: cover;
          border-radius: 11px;
          margin-bottom: 12px;
        }
        .home-regeneration-stage summary img {
          transition: filter .3s ease;
        }
        .home-regeneration-stage:hover summary img,
        .home-regeneration-stage:focus-within summary img,
        .home-regeneration-stage[open] summary img {
          filter: brightness(.48);
        }
        .home-regeneration-stage summary .stage-heading {
          display: flex;
          align-items: center;
          gap: 10px;
          font-size: 16px;
        }
        .home-regeneration-stage .stage-number {
          flex: 0 0 30px;
          width: 30px;
          height: 30px;
          margin: 0;
          font-size: 13px;
        }
        .home-regeneration-stage .stage-description {
          position: absolute;
          top: 12px;
          left: 12px;
          right: 12px;
          height: 220px;
          box-sizing: border-box;
          display: flex;
          align-items: center;
          justify-content: center;
          padding: 20px;
          margin: 0;
          overflow-y: auto;
          border-radius: 11px;
          background: rgba(255, 253, 252, .94);
          color: #59494A;
          font-size: 13px;
          line-height: 1.7;
          text-align: left;
          opacity: 0;
          visibility: hidden;
          pointer-events: none;
          transition: opacity .3s ease, visibility .3s ease;
          z-index: 2;
        }
        .home-regeneration-stage:hover .stage-description,
        .home-regeneration-stage:focus-within .stage-description,
        .home-regeneration-stage[open] .stage-description {
          opacity: 1;
          visibility: visible;
        }
        .home-regeneration-stage[open] {
          border-color: #C28F99;
          box-shadow: 0 7px 20px rgba(89, 73, 74, .10);
        }
        @media (max-width: 850px) {
          .home-regeneration-stages {
            grid-template-columns: 1fr;
          }
          .home-regeneration-stage summary img {
            height: auto;
            max-height: 340px;
            object-fit: contain;
            background: #F7F0E8;
          }
        }

        .fun-facts-box {
          background: #DCEEF6;
          border: 1px solid #D2E7F0;
          border-radius: 16px;
          padding: 18px;
          margin: 8px 0 24px;
          color: #294B5D;
        }
        .fun-facts-title {
          font-family: 'Caveat', 'Comic Sans MS', cursive;
          font-size: 38px;
          font-weight: 600;
          line-height: 1.1;
          text-align: center;
          color: #315F78;
          margin: 0 0 3px;
        }
        .fun-facts-subtitle {
          text-align: center;
          color: #416D82;
          font-size: 11px;
          margin: 0 0 12px;
        }
        .fun-facts-grid {
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 9px;
        }
        .fun-fact-item {
          background: rgba(255,255,255,.76);
          border: 1px solid #DCEBF1;
          border-radius: 10px;
          padding: 11px;
          min-width: 0;
        }
        .fun-fact-number {
          display: block;
          color: #3F7188;
          font-size: 10px;
          font-weight: 800;
          letter-spacing: 1px;
          margin-bottom: 3px;
        }
        .fun-fact-item h4 {
          color: #294F63;
          font-size: 14px;
          font-weight: 800;
          margin: 0 0 4px;
        }
        .fun-fact-item p {
          color: #344F5D;
          font-size: 12px;
          line-height: 1.5;
          margin: 0;
        }
        .fun-facts-footer {
          color: #66899A;
          text-align: center;
          font-size: 11px;
          font-style: italic;
          margin: 12px 0 0;
        }
        @media (max-width: 650px) {
          .fun-facts-box { padding: 14px; }
          .fun-facts-grid { grid-template-columns: 1fr; }
          .fun-facts-title { font-size: 35px; }
        }
      ")),

      h3(
        class = "section-title",
        "Stages of Axolotl Limb Regeneration"
      ),

      p(
        "Tap on the image to learn what happens at each stage.",
        style = "color:#59494A;line-height:1.7;margin-bottom:14px;"
      ),

      div(
        class = "home-regeneration-stages",

        tags$details(
          class = "home-regeneration-stage",
          tags$summary(
            tags$img(
              src = "axoltl_wound.png",
              alt = "Axolotl limb during the wound-healing stage of regeneration"
            ),
            div(
              class = "stage-heading",
              span(class = "stage-number", "1"),
              span("Wound healing")
            )
          ),
          div(
            class = "stage-description",
            "Following amputation, the wound surface closes and a specialised wound epidermis forms over the injured limb. This establishes the local environment needed for regeneration."
          )
        ),

        tags$details(
          class = "home-regeneration-stage",
          tags$summary(
            tags$img(
              src = "axoltl_blastoma.png",
              alt = "Axolotl limb showing blastema formation during regeneration"
            ),
            div(
              class = "stage-heading",
              span(class = "stage-number", "2"),
              span("Blastema formation")
            )
          ),
          div(
            class = "stage-description",
            "Cells accumulate beneath the wound epidermis and form a blastema, a group of proliferating, regeneration-competent cells that contributes to rebuilding the missing limb."
          )
        ),

        tags$details(
          class = "home-regeneration-stage",
          tags$summary(
            tags$img(
              src = "axoltl_redifferentiation.png",
              alt = "Axolotl limb during redifferentiation and limb outgrowth"
            ),
            div(
              class = "stage-heading",
              span(class = "stage-number", "3"),
              span("Redifferentiation")
            )
          ),
          div(
            class = "stage-description",
            "Blastema cells differentiate into specialised tissues, which organise and grow to reconstruct the structures of the missing limb, including the digits."
          )
        )
      ),

      # AXOLOTL FUN FACTS
      div(
        class = "fun-facts-box",

        h2(
          class = "fun-facts-title",
          "Fun Facts"
        ),

        p(
          class = "fun-facts-subtitle",
          "A little creature with extraordinary abilities!"
        ),

        div(
          class = "fun-facts-grid",

          div(
            class = "fun-fact-item",
            span(class = "fun-fact-number", "01"),
            h4("Forever young"),
            p("Axolotls retain juvenile features, including their external gills, even after reaching sexual maturity.")
          ),

          div(
            class = "fun-fact-item",
            span(class = "fun-fact-number", "02"),
            h4("Super regenerators"),
            p("They can regenerate limbs, tails, spinal cord tissue and parts of organs such as the heart and brain.")
          ),

          div(
            class = "fun-fact-item",
            span(class = "fun-fact-number", "03"),
            h4("Three ways to breathe"),
            p("They exchange gases through their feathery external gills, lungs and skin.")
          ),

          div(
            class = "fun-fact-item",
            span(class = "fun-fact-number", "04"),
            h4("Remarkable tissue repair"),
            p("Axolotls can regenerate many injured tissues with little or no typical permanent scarring.")
          ),

          div(
            class = "fun-fact-item",
            span(class = "fun-fact-number", "05"),
            h4("Metamorphosis is possible"),
            p("Although they usually remain aquatic, axolotls can undergo metamorphosis under certain hormonal conditions.")
          ),

          div(
            class = "fun-fact-item",
            span(class = "fun-fact-number", "06"),
            h4("Tiny predators"),
            p("They are carnivores that feed on worms, insect larvae, small crustaceans and other aquatic animals.")
          ),

          div(
            class = "fun-fact-item",
            span(class = "fun-fact-number", "07"),
            h4("A giant genome"),
            p("Their genome contains about 32 billion DNA base pairs — roughly ten times the size of the human genome.")
          ),

          div(
            class = "fun-fact-item",
            span(class = "fun-fact-number", "08"),
            h4("A species at risk"),
            p("Wild axolotls are critically endangered, with habitat loss, pollution and invasive species threatening their survival.")
          )
        ),

        p(
          class = "fun-facts-footer",
          "Small amphibian. Big surprises!"
        )
      ),

      h3(
        class = "section-title",
        "The focus of this project"
      ),

      div(
        class = "focus-panel",

        h4("What changes during regeneration?"),

        p(
          paste(
            "This project investigates gene-expression",
            "changes across the 3 and 5 days",
            "post-amputation stages, with the aim of",
            "identifying molecular patterns and candidate",
            "genes associated with regeneration."
          )
        )
      ),

      div(
        class = "footer-band",

        icon("leaf"),

        " Explore the data. Investigate the genes. Understand regeneration."
      )
    )
  })


  # ==========================================================
  # NAVIGATION BUTTONS
  # ==========================================================

  observeEvent(input$open_metadata, {
    current_page("metadata")
  })


  observeEvent(input$back_about_metadata, {
    current_page("about")
  })


  observeEvent(input$open_workflow, {
    current_page("workflow")
  })


  observeEvent(input$back_about_workflow, {
    current_page("about")
  })


  # WORKFLOW STAGE SELECTION

  observeEvent(input$workflow_stage_prep, {
    workflow_selected("prep")
  })

  observeEvent(input$workflow_stage_deseq, {
    workflow_selected("deseq")
  })

  observeEvent(input$workflow_stage_significant, {
    workflow_selected("significant")
  })

  observeEvent(input$workflow_stage_go, {
    workflow_selected("go")
  })

  observeEvent(input$workflow_stage_human, {
    workflow_selected("human")
  })

  observeEvent(input$workflow_stage_enrichment, {
    workflow_selected("enrichment")
  })

  observeEvent(input$workflow_stage_visualization, {
    workflow_selected("visualization")
  })


  # RESULTS & INTERPRETATION NAVIGATION

  observeEvent(input$results_open_overview, {
    results_section("overview_details")
  })

  observeEvent(input$results_open_deseq2, {
    results_section("deseq2")
  })

  observeEvent(input$results_open_upset, {
    results_section("upset")
  })

  observeEvent(input$results_open_interpretation, {
    results_section("interpretation")
  })

  observeEvent(input$results_back_home, {
    results_section("overview")
  })


  # SIDEBAR NAVIGATION

  observeEvent(input$page, {

    req(input$page)

    current_page(input$page)

    if (input$page == "genes") {
      gene_mapping_page("overview")
    }

    if (input$page == "enrichment") {
      enrichment_selected("overview")
    }

  }, ignoreInit = TRUE)


  # ==========================================================
  # FUNCTIONAL ENRICHMENT NAVIGATION
  # ==========================================================

  observeEvent(input$open_kegg, {
    enrichment_selected("kegg")
  })


  observeEvent(input$open_go, {
    enrichment_selected("go")
  })


  observeEvent(input$open_pathways, {
    enrichment_selected("pathways")
  })


  observeEvent(input$open_individual_pathways, {
    enrichment_selected("individual_pathways")
  })


  observeEvent(input$back_to_pathways, {
    enrichment_selected("pathways")
  })


  observeEvent(input$back_enrichment, {
    enrichment_selected("overview")
  })


  # GENE MAPPING NAVIGATION

  observeEvent(input$open_gene_mapping_overview, {
    gene_mapping_page("overview_details")
  })


  observeEvent(input$open_gene_mapping_results, {
    gene_mapping_page("mapping_results")
  })


  observeEvent(input$back_gene_mapping, {
    gene_mapping_page("overview")
  })


  observeEvent(input$back_gene_mapping_overview, {
    gene_mapping_page("overview")
  })


  # OPEN THE DESeq2 COMPARISON LANDING PAGE

  observeEvent(input$de_comparison, {
    current_page("comparison")
  })


  # OPEN THE COMPARISON TABLE PAGE

  observeEvent(input$open_comparison_table, {
    current_page("comparison_table")
  })


  # OPEN THE EXPRESSION LINE PLOT PAGE

  observeEvent(input$open_expression_line_plot, {
    current_page("expression_line_plot")
  })


  # RETURN TO THE DESeq2 COMPARISON OPTIONS

  observeEvent(input$back_comparison_options, {
    current_page("comparison")
  })


  # RETURN TO DIFFERENTIAL EXPRESSION

  observeEvent(input$back_de_comparison, {
    current_page("de")
  })


  # OPEN REGULATED GENES PAGE

  observeEvent(input$de_regulated, {
    current_page("regulated")
  })


  # RETURN TO DIFFERENTIAL EXPRESSION

  observeEvent(input$back_de_regulated, {

    current_page("de")

    bslib::nav_select(
      session,
      id = "page",
      selected = "de"
    )
  })


  # OPEN THE UPREGULATED GENES TABLE

  observeEvent(input$view_upregulated, {
    current_page("upregulated")
  })


  # RETURN TO THE REGULATED GENES PAGE

  observeEvent(input$back_regulated_up, {
    current_page("regulated")
  })


  # OPEN THE DOWNREGULATED GENES TABLE

  observeEvent(input$view_downregulated, {
    current_page("downregulated")
  })


  # RETURN TO THE REGULATED GENES PAGE

  observeEvent(input$back_regulated_down, {
    current_page("regulated")
  })


  # RETURN TO HOME FROM PLACEHOLDER PAGES

  observeEvent(input$back_home, {

    current_page("home")

    bslib::nav_select(
      session,
      id = "page",
      selected = "home"
    )
  })


  # EXISTING NAVIGATION BUTTONS

  observeEvent(input$explore_data, {

    bslib::nav_select(
      session,
      id = "page",
      selected = "de"
    )
  })


  observeEvent(input$go_de, {

    bslib::nav_select(
      session,
      id = "page",
      selected = "de"
    )
  })
}


shinyApp(ui = ui, server = server)
