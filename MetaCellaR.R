### To summerize single cells based on k-medoids ###
### The clustering is done per cell type ###
### parameter k is defined based on the number of potential groups that have at least 30 cells inside (k <- ncol(CT_cluster) / 30)

set.seed(0)
library(cluster)
library(knn.covertree)
library(irlba)
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
defined_args <- c("-file", "-RNA", "-celltype", "-output", "-k", "-assay", "-umap")
arg_tokens <- unlist(strsplit(argsexpr, split= " "))
file_hit <- which(arg_tokens == defined_args[1])
RNA_hit <- which(arg_tokens == defined_args[2])
celltype_hit <- which(arg_tokens == defined_args[3])
output_hit <- which(arg_tokens == defined_args[4])
k_hit <- which(arg_tokens == defined_args[5])
assay_hit <- which(arg_tokens == defined_args[6])
umap_hit <- which(arg_tokens == defined_args[7])
if(length(umap_hit)){
	umap_flag <- as.logical(arg_tokens[umap_hit + 1])
}else{
	umap_flag <- F
}
if(length(file_hit)){
  if(length(RNA_hit) && length(celltype_hit)){
    file_name <- arg_tokens[file_hit + 1]
    RNA_count_expression <- arg_tokens[RNA_hit + 1]
    celltype_expression <- arg_tokens[celltype_hit + 1]
    assay_expression <- arg_tokens[assay_hit + 1]
  }
  else{# CSV file
    if(!length(celltype_hit)){
      stop("For the csv input file the -celltype options must be present! Please provide a two column file with headers, linking the cell names (first column) to cell types (second columns).")
    }
    print(paste("Reading", arg_tokens[file_hit + 1], "..."))
    csv_flag <- T
    csv_data <- read.csv(arg_tokens[file_hit + 1], row.names= 1)
    csv_cells <- read.table(arg_tokens[celltype_hit + 1], header= T)
		general_data <- csv_data
    print("Done!")
  }
}else{
  stop("Argument -file is missing. Please provide the path to your input file using the -file option followed by the path to your file!")
}
if(!length(output_hit)){
  output_file <- getwd()
}else{
  output_file <- arg_tokens[output_hit + 1]
}
########################
if(!csv_flag){
  library(Seurat)
  if(strsplit(file_name, "\\.")[[1]][2] == "rds"){
    print("Reading the Seurat object")
    Sub <- readRDS(file_name)
    Rds_name <- "Sub"
  }else{
    print("Loading the Seurat object")
    Rds_name <- load(file_name)
    Sub <- get(Rds_name) ##To make it generalizable for Suerat object stored with any name
  }
  print("Done loading the Seurat object")
  print("Done loading Rds!")
  RNAcounts <- eval(parse(text= paste0(Rds_name, "@", RNA_count_expression))) #Sub@assays$RNA@counts
	RNAcounts <- as.matrix(RNAcounts)
  celltypes <- eval(parse(text= paste0(Rds_name, "@", celltype_expression)))
  cell_types <- unique(as.character(celltypes))
  if(length(assay_hit)){
	  assays <- eval(parse(text= paste0(Rds_name, "@", assay_expression)))

	  rna_hits <- which(tolower(assays) == "scrna-seq" | tolower(assays) == "scrna" | tolower(assays) == "scrna" | tolower(assays) == "rna")
	  ATACcounts <- RNAcounts[, -rna_hits]
	  RNAcounts <- RNAcounts[, rna_hits]
		general_data <- RNAcounts
	  RNA_celltypes <- celltypes[rna_hits]
	  ATAC_celltypes <- celltypes[-rna_hits]
		RNA_umap <- Sub[["umap"]]@cell.embeddings[rna_hits, ]
		ATAC_umap <- Sub[["umap"]]@cell.embeddings[-rna_hits, ]
print("ATAC_umap")
print(head(ATAC_umap))
	  celltypes <- RNA_celltypes
	  cell_types <- unique(celltypes)
  }
}else{
  cell_types <- unique(csv_cells[, 2])
}
if(umap_flag){
	print("Computing UMAP")
	ptm <- proc.time()
	if(csv_flag){
		## Florian's way ##
		#pca_res <- prcomp_irlba(csv_data, n= 30, scale.= T)
		#umap_res <- umap::umap(t(csv_data), n_components= 20)
		#umap_res <- umap::umap(pca_res$rotation, n_components= 20, n_threads= 10)
		##Seurat way###
    ptm_umap <- proc.time(); 
    umap_layout <- uwot::umap(t(as.matrix(csv_data)), pca= 30, pca_center= T, n_components= 20)
		print(paste("UMAP computation time:", proc.time() - ptm_umap))
		celltypes <- csv_cells[, 2]
	}else{
		## Florian's way ##
		#print("Running PCA")
		#pca_res <- prcomp_irlba(RNAcounts, n=30, scale.= T)
		#print("Done running PCA")
		#print(paste("PCA computation time:", (proc.time() - ptm)))
		#umap_res <- umap::umap(t(as.matrix(RNAcounts)), n_components= 20)
		#ptm2 <- proc.time()
		#print("Starting UMAP calculation...")
		#umap_res <- umap::umap(pca_res$rotation, n_components= 20, n_threads= 10)
		#print("UMAP calculation done!")
		#print(paste("UMAP computation time:", (proc.time() - ptm)))
		##Seurat way###
		ptm_umap <- proc.time();
		umap_layout <- uwot::umap(t(as.matrix(RNAcounts)), pca= 30, pca_center= T, n_components= 20)
		print(paste("UMAP computation time:", proc.time() - ptm_umap))
	}
}
clusters <- list()
all_mediods <- list()
RNA_metacell_umap <- NULL
cell2metacell_info <- NULL
cluster_data <- list()
print(dim(RNAcounts))

for(ct in cell_types){
  print(ct)
  if(!csv_flag){
	  if(!length(assay_hit)){
	    CT_cluster <- RNAcounts[, rownames(eval(parse(text= paste0(Rds_name, "@meta.data")))[celltypes == ct, ])]
		}
    CT_cluster <- RNAcounts[, celltypes == ct]
  }else{
    CT_cluster <- csv_data[, csv_cells[csv_cells[, 2] == ct, 1]]
  }
	if(umap_flag){
		CT_cluster <- t(umap_layout[celltypes == ct, ]) # I need to transform it bcuz the data gets transformed in the clara call
	}
	original_CT_cluster <- CT_cluster
  print(c("dim(CT_cluster):", dim(CT_cluster)))
  if(length(k_hit)){
    k <- as.integer(arg_tokens[k_hit + 1])
  }else{
    k <- floor(ncol(CT_cluster) / 30)
    print(paste("Setting k to", k))
  }
  if(summary_method == "kmed" || summary_method == "kmed_means"){
    print(paste("k=", k, "ncol(CT_cluster)=", ncol(CT_cluster)))
    if(length(k) >0 && k > 3 && k < ncol(CT_cluster)){
      class(CT_cluster)
			cluster_data[[ct]] <- t(as.matrix(original_CT_cluster))
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
			RNA_metacell_umap_ct <- NULL
			for(i in unique(clusters[[ct]]$clustering)){
				data_subset <- RNA_umap[which(clusters[[ct]]$clustering == i), ];
				if(is.null(dim(data_subset))){
					RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, data_subset);
				}else{
					RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, colMeans(data_subset))
				}
			}
			rownames(RNA_metacell_umap_ct) <- paste(ct, seq(nrow(RNA_metacell_umap_ct)), sep= "_")
			cell2metacell_info <- c(cell2metacell_info, paste(ct, clusters[[ct]]$clustering, sep= "_"))
      print(paste("Done clustering", ct))
			RNA_metacell_umap <- rbind(RNA_metacell_umap, RNA_metacell_umap_ct)
    }
  }else if(summary_method == "kmeans"){
    if(length(k) > 0 && k > 3){
      print(dim(CT_cluster))
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
save(cluster_data, clusters, file= paste0(output_file, "_clusters_heart_debug.Rdata"))

#load("clusters_heart_debug.Rdata")

print("Done clustering!")
mat <- NULL
mat_sum <- NULL
for(i in seq(length(clusters))){
  if(summary_method == "kmed"){
    mat <- cbind(mat, t(clusters[[i]]$medoids))
  }else if(summary_method == "kmeans"){
    mat <- cbind(mat, t(clusters[[i]]$centers))
  }
  else{##kmed_means
    temp_cl <- NULL
    temp_cl_sum <- NULL
		if(umap_flag){
			for(clst in unique(clusters[[i]]$clustering)){
				if(class(cluster_data[[i]][clusters[[i]]$cluster == clst, ]) == "numeric"){
					temp_cl <- rbind(temp_cl, cluster_data[[i]][clusters[[i]]$cluster == clst, ])
					temp_cl_sum <- rbind(temp_cl_sum, cluster_data[[i]][clusters[[i]]$cluster == clst, ])
				}else{
					temp_cl <- rbind(temp_cl, apply(cluster_data[[i]][clusters[[i]]$cluster == clst, ], 2, FUN= mean))
					temp_cl_sum <- rbind(temp_cl_sum, apply(cluster_data[[i]][clusters[[i]]$cluster == clst, ], 2, FUN= sum))
				}
			}
		}else{
    	for(clst in unique(clusters[[i]]$clustering)){
      	if(class(clusters[[i]]$data[clusters[[i]]$cluster == clst, ]) == "numeric"){
        	temp_cl <- rbind(temp_cl, clusters[[i]]$data[clusters[[i]]$cluster == clst, ])
        	temp_cl_sum <- rbind(temp_cl_sum, clusters[[i]]$data[clusters[[i]]$cluster == clst, ])
      	}
      	else{
        	temp_cl <- rbind(temp_cl, apply(clusters[[i]]$data[clusters[[i]]$cluster == clst, ], 2, FUN= mean))
        	temp_cl_sum <- rbind(temp_cl_sum, apply(clusters[[i]]$data[clusters[[i]]$cluster == clst, ], 2, FUN= sum))
      	}
    	}
		}
		print(dim(temp_cl))
		print(dim(temp_cl_sum))
    mat <- cbind(mat, t(temp_cl))
    mat_sum <- cbind(mat_sum, t(temp_cl_sum))
  }
}

print("done making mat")

mc_names <- NULL;
for(i in seq(length(clusters)))
  if(summary_method == "kmed" || summary_method == "kmed_means"){
    mc_names <- c(mc_names, paste(names(clusters)[i], seq(nrow(clusters[[i]]$medoids)), sep= "_"))
  }else if(summary_method == "kmeans"){
    mc_names <- c(mc_names, paste(names(clusters)[i], seq(nrow(clusters[[i]]$centers)), sep= "_"))
  }
colnames(mat) <- mc_names
colnames(mat_sum) <- mc_names

print("creating directory!")
dir.create(output_file)
if(length(assay_hit)){
	kk <- 5L
	knn_res <- class::knn(train= RNA_metacell_umap, test= ATAC_umap, cl= rownames(RNA_metacell_umap), k= kk)
	atac2metacell_info <- data.frame(barcode= rownames(ATAC_umap), metacell= knn_res)
	write.csv(atac2metacell_info, paste0(output_file, "/ATAC_cell2metacell_info_", summary_method, ".csv"), row.names= F)
}
rna2metacell_info <- data.frame(barcode= rownames(RNA_umap), metacell= cell2metacell_info)
write.csv(rna2metacell_info, paste0(output_file, "/RNA_cell2metacell_info_", summary_method, ".csv"), row.names= F)
write.csv(mat, paste0(output_file, "/cellSummarized_", summary_method, ".csv"))
write.csv(mat_sum, paste0(output_file, "/cellSummarized_", summary_method, "_sum.csv"))
save(atac2metacell_info, ATACcounts, clusters, RNA_metacell_umap, ATAC_umap, mc_names, file= paste0(output_file, "/", summary_method, "_clustered.RData"))

uniq_mc <- unique(atac2metacell_info$metacell)
atac_metacell <- NULL;
for(i in seq(length(uniq_mc))){
	hits <- atac2metacell_info$barcode[which(atac2metacell_info$metacell == uniq_mc[i])];
	if(length(hits) > 1){
		atac_metacell <- cbind(atac_metacell, rowMeans(as.matrix(ATACcounts)[, hits]))
	}else{
		atac_metacell <- cbind(atac_metacell, as.matrix(ATACcounts)[, hits])
	}
}
colnames(atac_metacell) <- uniq_mc
write.csv(atac_metacell, paste0(output_file, "/cellSummarized_ATAC_", summary_method, ".csv"))

final_umap_res <- umap::umap(t(mat))

celltypes <- sapply(colnames(mat), function(i) strsplit(i, "_")[[1]][1])
df <- data.frame(UMAP1= final_umap_res$layout[, 1], UMAP2= final_umap_res$layout[, 2], celltype= celltypes)
library(ggplot2)
pdf(paste0(output_file, "/umap_", summary_method, ".pdf"))
ggplot(df, aes(x= UMAP1, y= UMAP2)) + geom_point(aes(color= celltype)) + theme_classic() + geom_text(aes(label= celltype),hjust=0, vjust=0, size= 3, check_overlap = T)
dev.off()
print("Done!")
Rtnse_plot <- F
if(Rtnse_plot){
  library(Rtsne)
  Rtsne_whole_res <- Rtsne(as.matrix(RNAcounts), check_duplicates= F)
  pdf(paste0(output_file, "/tSNE_", summary_method, "_Rtsne.pdf"))
  plot(Rtsne_whole_res$Y)
  dev.off()
}
