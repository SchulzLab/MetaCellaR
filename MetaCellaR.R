### To summerize single cells based on k-medoids ###
### The clustering is done per cell type ###
### parameter k, if not provided by the user, is defined based on the number of potential groups that have at least 30 cells inside (k <- ncol(CT_cluster) / 30)

set.seed(0)
library(Seurat)
library(cluster)
########################
summary_method <- "kmed_means" #kmed
iter_flag <- F
csv_flag <- F
#args <- commandArgs(trailingOnly= T)
#file_name <- args[1] # Path to where the Seurat object is stored
#RNA_count_expression <- args[2]
#celltype_expression <- args[3]
########################
## Parse arguments
argsexpr <- commandArgs(trailingOnly= T)
defined_args <- c("-file", "-RNA", "-celltype", "-output", "-k")
arg_tokens <- unlist(strsplit(argsexpr, split= " "))
file_hit <- which(arg_tokens == defined_args[1])
RNA_hit <- which(arg_tokens == defined_args[2])
celltype_hit <- which(arg_tokens == defined_args[3])
output_hit <- which(arg_tokens == defined_args[4])
k_hit <- which(arg_tokens == defined_args[5])
if(length(file_hit)){
  if(length(RNA_hit) && length(celltype_hit)){
    file_name <- arg_tokens[file_hit + 1]
    RNA_count_expression <- arg_tokens[RNA_hit + 1]
    celltype_expression <- arg_tokens[celltype_hit + 1]
  }
  else{# CSV file
    if(!length(celltype_hit)){
      stop("For the csv input file the -celltype options must be present! Please provide a two column file with headers, linking the cell names (first column) to cell types (second columns).")
    }
    print(paste("Reading", arg_tokens[file_hit + 1], "..."))
    csv_flag <- T
    csv_data <- read.csv(arg_tokens[file_hit + 1], row.names= 1)
    csv_cells <- read.table(arg_tokens[celltype_hit + 1], header= T)
    print("Done!")
  }
}else{
  stop("Argument -file is missing. Please provide the path to your input file using the -file option followed by the path to your file!")
}
if(!length(output_hit)){
  output_file <- file_name
}else{
  output_file <- arg_tokens[output_hit + 1]
  ifelse(!dir.exists(output_file), dir.create(output_file), FALSE)
}
########################
if(!csv_flag){
  Rds_name <- load(file_name)
  Sub <- get(Rds_name) ##To make it generalizable for Suerat object stored with any name
  print("Done loading Rds!")
  RNAcounts <- eval(parse(text= paste0(Rds_name, "@", RNA_count_expression))) #Sub@assays$RNA@counts
  celltypes <- eval(parse(text= paste0(Rds_name, "@", celltype_expression)))
  cell_types <- unique(as.character(celltypes))
}else{
  cell_types <- unique(csv_cells[, 2])
}
clusters <- list()
all_mediods <- list()

for(ct in cell_types){
  print(ct)
  if(!csv_flag){
    CT_cluster <- RNAcounts[, rownames(eval(parse(text= paste0(Rds_name, "@meta.data")))[celltypes == ct, ])]
  }else{
    CT_cluster <- csv_data[, csv_cells[csv_cells[, 2] == ct, 1]]
  }

  if(length(k_hit)){
    k <- as.integer(arg_tokens[k_hit + 1])
  }else{
    k <- floor(ncol(CT_cluster) / 30)
    print(paste("Setting k to", k))
  }
  if(summary_method == "kmed" || summary_method == "kmed_means"){
    if(length(k) >0 && k > 3){
      class(CT_cluster)
      if(iter_flag){
        for(iter in seq(10)){
          clara_res <- clara(t(as.matrix(CT_cluster)), k, metric = "euclidean", stand = FALSE, samples= 30, pamLike = FALSE)
          dim(clara_res$medoids)
          dim(CT_cluster)

          clusters[[ct]] <- clara_res
          all_mediods[[ct]] <- rbind(all_mediods[[ct]], clara_res$medoids)
        }
      }
      else{
        clusters[[ct]] <- clara(t(as.matrix(CT_cluster)), k, metric = "euclidean", stand = FALSE, samples= 30, pamLike = FALSE)
      }
      print(paste("Done clustering", ct))
    }
  }else if(summary_method == "kmeans"){
    if(length(k) >0 && k > 3){
      if(class(CT_cluster) == "numeric"){
        next;
      }
      clusters[[ct]] <- kmeans(t(as.matrix(CT_cluster)), k)
    }
  }
  else{
    error("Undefined method of summarization. Please pick either kmed or kmeans!")
  }
}
#save(clusters, file= "clusters_heart_debug.Rdata")

#load("clusters_heart_debug.Rdata")

mat <- NULL
for(i in seq(length(clusters))){
  if(summary_method == "kmed"){
    mat <- cbind(mat, t(clusters[[i]]$medoids))
  }else if(summary_method == "kmeans"){
    mat <- cbind(mat, t(clusters[[i]]$centers))
  }
  else{##kmed_means
    temp_cl <- NULL
    for(clst in unique(clusters[[i]]$clustering)){
      if(class(clusters[[i]]$data[clusters[[i]]$cluster == clst, ]) == "numeric"){
        temp_cl <- rbind(temp_cl, clusters[[i]]$data[clusters[[i]]$cluster == clst, ])
      }
      else{
        temp_cl <- rbind(temp_cl, apply(clusters[[i]]$data[clusters[[i]]$cluster == clst, ], 2, FUN= mean))
      }
    }
    mat <- cbind(mat, t(temp_cl))
  }
}


mc_names <- NULL;
for(i in seq(length(clusters)))
  if(summary_method == "kmed" || summary_method == "kmed_means"){
    mc_names <- c(mc_names, paste(names(clusters)[i], seq(nrow(clusters[[i]]$medoids)), sep= "_"))
  }else if(summary_method == "kmeans"){
    mc_names <- c(mc_names, paste(names(clusters)[i], seq(nrow(clusters[[i]]$centers)), sep= "_"))
  }
colnames(mat) <- mc_names

write.csv(mat, paste0(output_file, "_cellSummarized_", summary_method, ".csv"))
save(clusters, mc_names, file= paste0(output_file, "_", summary_method, "_clustered.RData"))
print(paste("Done writing results to the output files:", paste0(output_file, "_cellSummarized_", summary_method, ".csv"), "and", paste0(output_file, "_", summary_method, "_clustered.RData")))
Rtnse_plot <- F
if(Rtnse_plot){
  library(Rtsne)
  Rtsne_whole_res <- Rtsne(as.matrix(RNAcounts), check_duplicates= F)
  pdf(paste0(output_file, "_", summary_method, "_Rtsne.pdf"))
  plot(Rtsne_whole_res$Y)
  dev.off()
}
