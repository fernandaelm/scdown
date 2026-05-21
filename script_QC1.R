# ref: https://doi.org/10.1038/s41467-025-63752-0

setwd("/scratch/16974917")
getwd() # conferir se ta no diretório certo

rm(list = ls()) 

set.seed(123)

########################
### loading library  ###
########################

library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(SingleCellExperiment)
library(scDblFinder)
library(tidyverse)
library(chromo) 

##################
### Functions  ###
##################

markersFunction <- function(seurat_obj = scdown) {

  markers_df_filt <- markers_df %>%
    filter(marker %in% rownames(seurat_obj))

  features_named <- setNames(markers_df_filt$marker, markers_df_filt$feature_label)

  auxmark <- markers_df_filt %>%
    group_by(celltype) %>%
    summarise(n_markers = n(), .groups = "drop") %>%
    mutate(celltype = factor(celltype, levels = unique(markers_df_filt$celltype))) %>%
    arrange(celltype) %>%
    mutate(
      y    = length(unique(seurat_obj@meta.data$seurat_clusters)) * 1.05,
      yend = y,
      x    = 1,
      xend = n_markers,
      celltype = case_when(
        celltype == "NPC"          ~ "Neural Progenitor Cell",
        celltype == "Exc"          ~ "Excitatory Neuron",
        celltype == "Inh"          ~ "Inhibitory Neuron",
        celltype == "Inh Striatum" ~ "Inhibitory Neuron (Striatum)",
        celltype == "Astro P"      ~ "Astrocyte Progenitor",
        celltype == "OPC"          ~ "Oligodendrocyte Progenitor",
        celltype == "Oligo"        ~ "Oligodendrocyte",
        celltype == "Micro"        ~ "Microglia",
        celltype == "Vasc"         ~ "Vascular Cell",
        TRUE ~ as.character(celltype)
      )
    ) %>%
    as.data.frame()

  for (i in 2:nrow(auxmark)) {
    auxmark[i, "x"]    <- auxmark[i - 1, "xend"] + 1
    auxmark[i, "xend"] <- auxmark[i, "x"] + auxmark[i, "n_markers"] - 1
  }

  dot_plot <- DotPlot(
    object         = seurat_obj,
    features       = features_named,
    cols           = c("blue", "red"),
    cluster.idents = TRUE,
    scale          = TRUE
  ) +
    scale_x_discrete(labels = names(features_named)) +
    theme(
      plot.background   = element_rect(fill = "white"),
      panel.background  = element_rect(fill = "white"),
      panel.border      = element_blank(),
      axis.line         = element_line(color = "black"),
      legend.background = element_rect(fill = "white"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
      axis.text.y = element_text(size = 17, face = "bold")
    ) +
    labs(x = NULL, y = NULL)

  for (i in 1:nrow(auxmark)) {
    dot_plot <- dot_plot +
      annotate(
        "segment",
        x    = auxmark[i, "x"],
        y    = auxmark[i, "y"],
        xend = auxmark[i, "xend"],
        yend = auxmark[i, "yend"]
      ) +
      annotate(
        "text",
        x        = (auxmark[i, "x"] + auxmark[i, "xend"]) / 2,
        y        = auxmark[i, "y"] + 0.15,
        label    = auxmark[i, "celltype"],
        color    = "black",
        fontface = "bold",
        size     = 4
      ) +
      coord_cartesian(clip = "off")
  }

  return(dot_plot)
}

barfunc <- function(
    df,
    x_axis,
    y_axis,
    fill_groups  = FALSE,
    annot_loc    = "center",
    fill_col     = FALSE,
    title_name   = NULL,
    annot_text   = y_axis,
    size_annot        = 4,
    size_xaxis_text   = 15,
    size_xaxis_title  = 15,
    size_yaxis_text   = 15,
    size_yaxis_title  = 15,
    rotate_x     = 0,
    x_title      = FALSE,
    y_title      = FALSE,
    extra_annot  = NULL
) {
  bar_plot <- ggplot(df, aes(x = !!sym(x_axis), y = !!sym(y_axis))) +
    geom_bar(stat = "identity") +
    labs(title = title_name, x = x_axis, y = y_axis) +
    theme_minimal() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line    = element_line(color = "black", linewidth = 0.5),
      axis.title.x = element_text(face = "bold", size = size_xaxis_title),
      axis.title.y = element_text(face = "bold", size = size_yaxis_title),
      axis.text.y  = element_text(face = "bold", size = size_yaxis_text),
      axis.text.x  = element_text(
        face  = "bold",
        angle = rotate_x,
        hjust = ifelse(rotate_x != 0, 1, 0.5),
        size  = size_xaxis_text
      ),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5)
    )

  if (size_yaxis_text == 0) {
    bar_plot <- bar_plot + theme(axis.ticks.y = element_blank())
  }

  if (!isFALSE(fill_groups)) {
    bar_plot <- bar_plot +
      aes(fill = !!sym(fill_groups)) +
      labs(fill = fill_groups)
  }

  if (is.list(fill_col)) {
    bar_plot <- bar_plot + scale_fill_manual(values = fill_col)
  } else if (!isFALSE(fill_col)) {
    bar_plot <- bar_plot +
      geom_bar(stat = "identity", fill = fill_col) +
      theme(legend.position = "none")
  }

  if (!is.null(extra_annot)) {
    bar_plot <- bar_plot +
      annotate(
        geom  = "text",
        x     = extra_annot[[x_axis]],
        y     = extra_annot$y,
        label = extra_annot$annotation_text,
        color = extra_annot$color,
        size  = 3, angle = 90, fontface = "bold"
      )
  } else if (annot_loc == "center") {
    bar_plot <- bar_plot +
      geom_text(aes(label = !!sym(annot_text)),
                position = position_stack(vjust = 0.5), size = size_annot)
  } else if (annot_loc == "repel") {
    bar_plot <- bar_plot +
      geom_text_repel(aes(label = !!sym(annot_text)),
                      direction = "y", size = size_annot)
  } else if (annot_loc == "top") {
    y_max <- if (!isFALSE(fill_groups)) {
      max(df %>%
            group_by(!!sym(x_axis)) %>%
            summarise(total = sum(!!sym(y_axis)), .groups = "drop") %>%
            pull(total))
    } else {
      max(df[[y_axis]])
    }
    bar_plot <- bar_plot +
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 1.1 * y_max)
      ) +
      geom_text(aes(label = !!sym(annot_text)), vjust = -0.2, size = size_annot)
  }

  if (!isFALSE(x_title)) bar_plot <- bar_plot + xlab(x_title)
  if (!isFALSE(y_title)) bar_plot <- bar_plot + ylab(y_title)

  return(bar_plot)
}

############################
### creating directories ###
############################

if (!file.exists("results_thyreoid.v1")) {
  dir.create("results_thyreoid.v1")
  dir.create("results_thyreoid.v1/exploratory_analysis")
}

############################
###       counts         ###
############################

runs <- c( 
  "GSM8591182", "GSM8591180", "GSM8591184", "GSM8591183", "GSM8591181",  # Down
  "GSM8591178", "GSM8591179", "GSM8591177", "GSM8591176", "GSM8591175"   # Euploid
)

# raw counts
counts <- list()
for (run in runs) {
  cat("Lendo:", run, "\n")
  path <- paste0("data/counts/", run, "/")
  counts[[run]] <- ReadMtx(
    mtx      = paste0(path, "matrix.mtx.gz"),
    cells    = paste0(path, "barcodes.tsv.gz"),
    features = paste0(path, "features.tsv.gz"),
    feature.column = 2,
    mtx.transpose = FALSE
  )
}

# ver se os 10 arquivos carregaram
names(counts)
length(counts)

for (run in names(counts)) {
  cat(run, "->", dim(counts[[run]]), "\n")
}
 

# separate Seurat objects
seurats <- list()
for (run in names(counts)) {
  seurats[[run]] <- CreateSeuratObject(counts = counts[[run]], project = run)
}

# clear counts
rm(counts)
gc()

scdown <- merge(
  x = seurats$GSM8591182,
  y = seurats[names(seurats) != "GSM8591182"],
  add.cell.ids = runs,
  project = "allRuns"
)

# clear seurats
rm(seurats)
gc()

# joining layers
scdown <- JoinLayers(scdown)


############################
###      meta data       ###
############################

meta_raw <- read.csv("data/meta.csv")

meta <- meta_raw %>%
  dplyr::select(Sample.Name, sex, treatment, Developmental_stage) %>%
  dplyr::rename(
    donor               = Sample.Name,
    condition           = treatment,
    developmental_stage = Developmental_stage
  ) %>%
  dplyr::mutate(
    developmental_stage = gsub(" ", "", developmental_stage),
    condition = dplyr::case_when(
      condition == "Euploid"       ~ "CT",
      condition == "Down Syndrome" ~ "DS"
    )
  ) %>%
  dplyr::mutate_all(factor)

scdown@meta.data <- scdown@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(meta, by = c("orig.ident" = "donor")) %>%
  column_to_rownames("cell")

scdown[["percent.mt"]] <- PercentageFeatureSet(scdown, pattern = "^MT-")


#############################
###   genes of interest   ###
#############################

int_genes_df <- read.csv("data/genes_down.csv")
int_genes <- int_genes_df$genesofinterest


#cell type markers
markers_df <- read.csv("data/cellmarkers1.csv") %>%
  group_by(marker) %>%
  mutate(
    occurrence    = row_number(),
    feature_label = ifelse(n() > 1, paste0(marker, ".", occurrence), marker)
  ) %>%
  ungroup() %>%
  as.data.frame()

cellmarkers1 <- markers_df$marker

#genes_fig1 <- c('C1R', 'C1S', 'C4A', 'C4B', 'C5', 'C8G', 'CD46', 'CD59', 'CFI', 'FCN2','SLC1A2', 'SLC1A3')


############################
###   filtering cells    ###
############################

#scdown <- subset(scdown, subset = nCount_RNA > 200 & nCount_RNA < 60000) nao usaram no artigo

scdown <- subset(scdown, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

scdown <- subset(scdown, subset = percent.mt < 5)

# doublet detection
sce <- as.SingleCellExperiment(scdown)
sce <- scDblFinder(sce, samples = "orig.ident") 
scdown$doublet <- sce$scDblFinder.class
scdown <- subset(scdown, subset = doublet == "singlet")


##########################
### filtering features ###
##########################

scdown <- scdown[!grepl("MT-", rownames(scdown)), ]

chrObj <- chromoInitiate(
  data.frame(Symbol = rownames(scdown), log2fc = 0, pval = 1),
  gene_col = "Symbol",
  fc_col   = "log2fc",
  p_col    = "pval"
)
all_features <- chrObj@data %>%
  dplyr::select(-log2fc, -pval, -DEG, -entrezgene_id)

missing_int <- int_genes[!int_genes %in% all_features$Symbol]
if (length(missing_int) > 0) {
  warning("Os seguintes genes de interesse foram removidos pelo filtro chromo: ",
          paste(missing_int, collapse = ", "))
}

scdown <- scdown[rownames(scdown) %in% all_features$Symbol, ]

nonzero_genes <- rowSums(scdown[["RNA"]]$counts) > 0
scdown <- scdown[nonzero_genes, ]

all_features <- all_features[all_features$Symbol %in% rownames(scdown), ]


##################
### Processing ###
##################

scdown <- NormalizeData(scdown)

gc() # clears garbage (next step requires a lot of memory)

scdown <- FindVariableFeatures(scdown) # default 2000 features

top10 <- head(VariableFeatures(scdown), 10) # Highlight the 10 most highly variable genes

scdown <- ScaleData(scdown, features = rownames(scdown))

scdown <- RunPCA(scdown, features = VariableFeatures(object = scdown))

ElbowPlot(scdown, ndims = 30)

scdown <- RunUMAP(scdown, dims = 1:10) #definir direito dps de ver o elbow plot

scdown <- FindNeighbors(scdown, dims = 1:10, k.param = 30)

for (res in c(0.2, 0.5, 0.8)) {
  scdown <- FindClusters(scdown, resolution = res)
}

Idents(scdown) <- "RNA_snn_res.0.2"


#################################
### Batch effect verification ###
#################################

for (i in c("seurat_clusters", "orig.ident", "sex", "condition")) {
  clusters <- DimPlot(
    scdown,
    reduction = "umap",
    group.by  = i,
    pt.size   = 0.1,
    label     = ifelse(i == "seurat_clusters", TRUE, FALSE),
    raster    = FALSE
  ) +
    ggtitle("") +
    labs(color = i) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    coord_fixed() +
    theme(
      legend.text  = element_text(size = 26),
      legend.title = element_text(size = 28),
      axis.text    = element_text(size = 16),
      axis.title   = element_text(size = 22)
    )

  ggsave(paste0("results_thyreoid.v1/exploratory_analysis/umap_", i, ".png"),
         plot = clusters, width = 10, height = 10)
}

for (res in c("RNA_snn_res.0.2", "RNA_snn_res.0.5", "RNA_snn_res.0.8")) {
  p <- DimPlot(
    scdown,
    reduction = "umap",
    group.by  = res,
    pt.size   = 0.1,
    label     = TRUE,
    raster    = FALSE
  ) +
    ggtitle(res) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    coord_fixed() +
    theme(
      legend.text  = element_text(size = 26),
      legend.title = element_text(size = 28),
      axis.text    = element_text(size = 16),
      axis.title   = element_text(size = 22)
    )

  ggsave(paste0("results_thyreoid.v1/exploratory_analysis/umap_", res, ".png"),
         plot = p, width = 10, height = 10)
}

#########################
### Bias verification ###
#########################

aux <- scdown@meta.data %>%
  dplyr::select(condition, sex) %>%
  group_by(condition, sex) %>%
  summarise(num_cells = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(pct = num_cells / sum(num_cells) * 100) %>%
  ungroup()

bar_plot <- barfunc(aux, x_axis = "condition", y_axis = "pct", fill_groups = "sex",
                    annot_loc = "center", annot_text = "num_cells")
ggsave("results_thyreoid.v1/exploratory_analysis/pct_condition_sex.png",
       plot = bar_plot, width = 3, height = 3)


#######################
### Finding markers ###
#######################

DefaultAssay(scdown)

Idents(scdown) <- "seurat_clusters"

for (res in c("RNA_snn_res.0.2", "RNA_snn_res.0.5", "RNA_snn_res.0.8")) {
  Idents(scdown) <- res
  dp <- markersFunction()
  ggsave(paste0("results_thyreoid.v1/exploratory_analysis/canonical_markers_", res, ".png"),
         plot = dp, height = 5, width = 20)
}

Idents(scdown) <- "RNA_snn_res.0.2"

###################
### saving data ###
###################

saveRDS(scdown, file = "results_thyreoid.v1/scdown_processed.rds")
