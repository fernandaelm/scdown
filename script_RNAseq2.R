# ref: https://doi.org/10.1038/s41591-026-04211-1


setwd("C:/Users/User/Documents/fernanda")
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


##################
### functions  ###
##################

markersFunction <- function(){

  celltypes_markers <- markers_df %>%
    filter(marker %in% rownames(scdown)) %>%
    pull(marker)

  markers_df <- markers_df %>%
    filter(marker %in% rownames(scdown))

  auxmark <- markers_df %>%
    group_by(celltype) %>%
    summarise(n_markers = n()) %>%
    mutate(celltype = factor(celltype, levels = unique(markers_df$celltype))) %>%
    arrange(celltype) %>%
    mutate(
      y    = length(unique(scdown@meta.data$seurat_clusters)) * 1.05,
      yend = y,
      x    = 1,
      xend = n_markers,
      colour   = "#000000",
      celltype = case_when(
        celltype == "RG/Ast"       ~ "Radial Glia / Astrocytes",
        celltype == "AST"          ~ "Astrocytes",
        celltype == "NEU"          ~ "Neuron",
        celltype == "EXC"          ~ "Excitatory Neuron",  
        celltype == "INH"          ~ "Inhibitory Neuron",                   
        celltype == "VASC"         ~ "Vascular Cell",
        celltype == "OPC"          ~ "Oligodendrocyte Progenitor Cell",               
        celltype == "IPC"          ~ "Intermediate Progenitor Cell",                     
        celltype == "Vasc"         ~ "Vascular Cell"                  
      )
    ) %>%
    as.data.frame()

  for (i in 2:nrow(auxmark)) {
    auxmark[i, "x"]    <- auxmark[i - 1, "xend"] + 1
    auxmark[i, "xend"] <- auxmark[i, "x"] + auxmark[i, "n_markers"] - 1
  }

  # canonical markers dot plot
  dot_plot <- DotPlot(
    object         = scdown,
    features       = celltypes_markers,
    cols           = c("blue", "red"),
    cluster.idents = T,
    scale          = T
  ) +
    theme(
      plot.background   = element_rect(fill = "white"),
      panel.background  = element_rect(fill = "white"),
      panel.border      = element_blank(),
      axis.line         = element_line(color = "black"),
      legend.background = element_rect(fill = "white"),
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
      axis.text.y       = element_text(size = 17, face = "bold")
    )

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
        x      = (auxmark[i, "x"] + auxmark[i, "xend"]) / 2,
        y      = auxmark[i, "y"] + 0.15,
        label  = auxmark[i, "celltype"],
        color  = "black",
        fontface = "bold",
        size   = 4
      ) +
      coord_cartesian(clip = "off")
  }

  return(dot_plot)
}

barfunc <- function(
    df,
    x_axis,
    y_axis,
    fill_groups  = F,
    annot_loc    = "center",
    fill_col     = F,
    title_name   = NULL,
    annot_text   = y_axis,
    size_annot   = 4,
    size_xaxis_text  = 15,
    size_xaxis_title = 15,
    size_yaxis_text  = 15,
    size_yaxis_title = 15,
    rotate_x     = 0,
    x_title      = F,
    y_title      = F,
    extra_annot  = NULL
) {
  bar_plot <- ggplot(df, aes(x = !!sym(x_axis), y = !!sym(y_axis))) +
    geom_bar(stat = "identity") +
    labs(title = title_name, x = x_axis, y = y_axis) +
    theme_minimal() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      panel.border      = element_blank(),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(color = "black", linewidth = 0.5),
      axis.title.x      = element_text(face = "bold", size = size_xaxis_title),
      axis.title.y      = element_text(face = "bold", size = size_yaxis_title),
      axis.text.y       = element_text(face = "bold", size = size_yaxis_text),
      axis.text.x       = element_text(face = "bold", angle = rotate_x,
                                       hjust = ifelse(rotate_x != 0, 1, 0.5),
                                       size = size_xaxis_text),
      axis.ticks.x      = element_blank(),
      axis.ticks.y      = element_line(color = "black", linewidth = 0.5)
    )

  if (size_yaxis_text == 0) {
    bar_plot <- bar_plot + theme(axis.ticks.y = element_blank())
  }

  if (fill_groups != F) {
    bar_plot <- bar_plot +
      aes(fill = !!sym(fill_groups)) +
      labs(fill = fill_groups)
  }

  if (is.list(fill_col)) {
    bar_plot <- bar_plot + scale_fill_manual(values = fill_col)
  } else if (fill_col != F) {
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
    bar_plot <- bar_plot +
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 1.1 * max(ifelse(
          fill_groups != F,
          df %>% group_by(!!sym(fill_groups)) %>% pull(!!sym(y_axis)),
          df[[y_axis]]
        )))
      ) +
      geom_text(aes(label = !!sym(annot_text)), vjust = -0.2, size = size_annot)
  }

  if (x_title != F) bar_plot <- bar_plot + xlab(x_title)
  if (y_title != F) bar_plot <- bar_plot + ylab(y_title)

  return(bar_plot)
}

############################
### creating directories ###
############################

if (!file.exists("results_thyreoid.v2")) {
  dir.create("results_thyreoid.v2")
  dir.create("results_thyreoid.v2/exploratory_analysis")
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
  path <- paste0("data/counts/", run, "/")
  counts[[run]] <- ReadMtx(
    mtx = paste0(path, "matrix.mtx"),
    cells = paste0(path, "barcodes.tsv"),
    features = paste0(path, "features.tsv")
  )
}

# ver se os 10 arquivos carregaram
names(counts)
length(counts)

for (run in names(counts)) {
  cat(run, "->", dim(counts[[run]]), "\n")
}

# resultado
#GSM8591182 -> 36601 1007348 
#GSM8591180 -> 36601 1749432 
#GSM8591184 -> 36601 2322946 
#GSM8591183 -> 36601 1590680 
#GSM8591181 -> 36601 1705696 
#GSM8591178 -> 36601 1830185 
#GSM8591179 -> 36601 1669894 
#GSM8591177 -> 36601 1686505 
#GSM8591176 -> 36601 1659270 
#GSM8591175 -> 36601 831017 

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

# mudei essa parte pq o select nao tava usando o pacote do dplyr
meta_raw <- read.csv("data/meta.csv")

meta <- meta_raw %>%
  dplyr::select(Run, sex, treatment, Developmental_stage) %>%
  dplyr::rename(
    condition = treatment,
    donor = Run,
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

# including the sample meta data into the merged cell meta data (mudei tbm)
scdown@meta.data <- scdown@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(meta, by = c("orig.ident" = "donor")) %>%
  column_to_rownames("cell")


# creating mitochondrial reads percentage
scdown[["percent.mt"]] <- PercentageFeatureSet(scdown, pattern = "^MT-")


#############################
###   genes of interest   ###
#############################

int_genes_df <- read.csv("data/genes_down.csv")
int_genes    <- int_genes_df$genesofinterest
comp_genes   <- character(0)


#cell type markers
markers_df <- read.csv("data/cellmarkers2.csv")
markers_df <- markers_df[!duplicated(markers_df$marker), ] # keep first occurrence per marker
cellmarkers2 <- markers_df$marker

#genes_fig1 <- c('C1R', 'C1S', 'C4A', 'C4B', 'C5', 'C8G', 'CD46', 'CD59', 'CFI', 'FCN2','SLC1A2', 'SLC1A3')


############################
###   filtering cells    ###
############################

scdown <- subset(scdown, subset = nCount_RNA > 500 & nCount_RNA < 30000)

#scdown <- subset(scdown, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) # was defined separately

scdown <- subset(scdown, subset = percent.mt < 2)



# doublet detection
sce <- as.SingleCellExperiment(scdown)
sce <- scDblFinder(sce, samples = "orig.ident") # mudei
scdown$doublet <- sce$scDblFinder.class
scdown <- subset(scdown, subset = doublet == "singlet")


# clear memory
rm(sce)
gc()

# save rds
saveRDS(scdown, file = "results_thyreoid.v2/scdown_singlets.rds")

dim(scdown) # ver quantos genes e cels tem

# scdown <- readRDS("results_thyreoid.v2/scdown_singlets.rds")