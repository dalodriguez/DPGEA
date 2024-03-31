#' DPGEA: differential pseudospatial gradient expression analysis
#'
#' @param cds a cds file where pseudotime has already been calculated
#' @param condition the column of "cds@clusters@listData that contains the condition names
#' @param threshold the p value threshold, default is 0.01
#' @return A list of dataframes including correlations, p-values and the final classification
#' @import Seurat Hmisc monocle3
#' @export

DPGEA <- function(cds, condition, threshold= 0.01) {

        pseudotime <- data.frame(pseudotime = cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]],
                                 cell_type = cds@clusters@listData[["UMAP"]][["clusters"]],
                                 condition = cds@colData@listData[[condition]],
                                 cell = rownames(as.data.frame(cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]])) )

        pseudotime <- pseudotime[order(pseudotime$pseudotime),] ; rownames(pseudotime) <- NULL

        DPGEA <- list()

        for (i in 1:length(unique(pseudotime$condition))){
              print(paste0("Starting condition ",i,": ",unique(pseudotime$condition)[i]))
                  pseu <- pseudotime[pseudotime$condition == pseudotime$condition[i],]
                  pseu <- pseu[order(pseu$pseudotime),]
              message(paste0("...Generating binary matrix"))
                  for (t in 1:length(unique(pseudotime$cell_type))){
                    for (j in 0:2){
                      if (t+j < length(unique(pseudotime$cell_type))+1){
                        pseu[,(ncol(pseu)+1)] = 0
                        pseu[,ncol(pseu)][pseu$cell_type %in% c(unique(pseudotime$cell_type)[t:(t+j)])] = 1
                        colnames(pseu)[ncol(pseu)] = do.call(paste, c(as.list(paste(unique(pseudotime$cell_type)[t:(t+j)])), sep = "_"))
                      }
                    }
                  }
                  matrix <- exprs(cds)
                  matrix <- matrix[,pseu$cell]
                  matrix <- t(apply(matrix,1,function(x){x}))
                  pseubis <- t(as.matrix(pseu)[,5:ncol(pseu)])
                  colnames(pseubis) <- pseu$cell
                  matrix <- rbind(pseubis,matrix)

              message(paste0("...Calculating correlations"))
                 cors <- rcorr(t(matrix))

              message(paste0("...Adjusting p-values"))
                  cors_r <- t(cors$r[c(1:nrow(pseubis)),])
                  cors_r <- as.data.frame(cors_r[c((nrow(pseubis)+1):nrow(cors_r)),])
  
                  DPGEA[["correlations"]][[unique(pseudotime$condition)[i]]] <- cors_r
    
                  cors_p <- t(cors$P[c(1:nrow(pseubis)),])
                  cors_p <- as.data.frame(cors_p[c((nrow(pseubis)+1):nrow(cors_p)),])
                  cors_p <- apply(cors_p, 2, function(x) p.adjust(x, method = "fdr", n = length(x)))
                  cors_p <- as.data.frame(cors_p)
    
                  DPGEA[["p_values"]][[unique(pseudotime$condition)[i]]] <- cors_p

              message(paste0("...Determining only significant correlations"))
                  only_sig <- list()
                  for (k in 1:ncol(cors_p)){
                    crr_p <- cors_p[cors_p[k] < threshold,]
                    crr_r <- cors_r[rownames(cors_r) %in% rownames(crr_p),]
                    crr_r <- crr_r[crr_r[k] > 0,]
                    crr_r <- crr_r[order(-crr_r[,k]),]
                    only_sig <- c(only_sig, list(crr_r))
                  }
    
                  names(only_sig) <- colnames(cors_r)
                  sig_genes <- c()
                  for (s in 1:length(only_sig)){
                    genes <- as.vector(rownames(only_sig[[s]]))
                    sig_genes <- c(sig_genes,genes)}
    
                  sig_genes <- unique(sig_genes)
                  cors_r_only_sig <- cors_r[rownames(cors_r) %in% sig_genes,]

              message(paste0("...Generating classification"))
                  pseu_binary <- pseu[,5:ncol(pseu)]
                  prod_matrix <- (matrix[(nrow(pseubis)+1):nrow(matrix),] > 0) %*% as.matrix(pseu_binary)
    
                  prod_matrix  <- prod_matrix / rowSums(matrix[(nrow(pseubis)+1):nrow(matrix),] > 0)
                  prod_matrix <- as.matrix(prod_matrix)  * as.matrix(cors_r)
                  prod_matrix[is.na(prod_matrix)] = 0
                  prod_matrix <- prod_matrix[rownames(prod_matrix) %in%rownames(cors_r_only_sig),]
    
                  prod_matrix_class <- prod_matrix
                  prod_matrix_class[prod_matrix_class<0.2] = 0
                  prod_matrix_class = prod_matrix_class[rowSums(prod_matrix_class) > 0,]
                  classification <- which(prod_matrix_class == rowMax(prod_matrix_class),arr.ind = T)
                  classification <- as.data.frame(classification)
                  classification$classification <- colnames(prod_matrix)[classification$col]
                  classification$row <- NULL
                  classification$genes <- rownames(classification)
                  rownames(classification) <- NULL
                  classification$condition <- unique(pseudotime$condition)[i]
                  classification$col <- NULL
    
                  DPGEA[["classification"]][[unique(pseudotime$condition)[i]]] <- classification
              message("DONE!")
              
        }
        return(DPGEA)
        }

