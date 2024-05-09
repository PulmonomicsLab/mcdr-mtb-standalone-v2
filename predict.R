args=commandArgs(trailingOnly=TRUE)
# print(length(args))
if (length(args) != 3 )
  stop("Invalid number of arguments to Rscript.")

script_path <- args[1]
input_file <- args[2]
output_path <- args[3]
# print(script_path)
# print(input_file)
# print(output_path)

convertRI <- function(diff){
  if(diff < 0.1)
    return(1)

  return(ceiling(round(diff*10, digits = 0)/2))
}

library(reshape2, quietly=TRUE)
library(caret, quietly=TRUE)

A_col <- unlist(read.table(paste0(script_path, "A_col.txt")))

In <- read.table(input_file, comment.char = "", sep = "\t", header = TRUE)

if(length(unique(In$SAMPLE))>1)
  stop("More than 1 samples are found in the VCF file.")

Ind <- In[,c("POS", "REF", "ALT", "QUAL", "SAMPLE", "AO", "DP")]
Ind$AO[Ind$AO=="."] <- 0
Ind$DP[Ind$DP=="."] <- 0
Ind <- transform(Ind, DP = as.numeric(DP))
Ind <- transform(Ind, AO = as.numeric(AO))
Ind$MUT <- paste(Ind$REF, Ind$POS, Ind$ALT, sep = "")
Ind$VAL <- Ind$AO/Ind$DP
Ind$VAL[Ind$VAL == "NaN"] <- 0.0000
Ind$VAL <- round(Ind$VAL, digits = 4)
Ind_l <- dcast(Ind[,c("SAMPLE", "MUT", "VAL")], SAMPLE~MUT, value.var = "VAL")

data_2 <- data.frame(matrix(ncol = length(A_col), nrow = nrow(Ind_l)))
colnames(data_2) <- A_col

for(i in colnames(data_2)){
  if(i %in% colnames(Ind_l)){
    data_2[1,i] <- Ind_l[1,i]
  }else{
    data_2[1,i] <- 0.00
  }
}


#### Prediction with full MLP model ####

mlp <- readRDS(paste0(script_path, "MLP_model.rds"))

pred_prob_mlp_ind <- predict(mlp, data_2, type="prob")
pred_mlp_ind <- predict(mlp, data_2)
# print(colnames(pred_prob_mlp_ind))
# print(pred_prob_mlp_ind)

prob <- n <- pred_prob_mlp_ind[1, c("S", "M", "P", "X")]

MLP_out <- prob

  MLP_out$sum <- rowSums(prob)
  MLP_out$S <- MLP_out$S / MLP_out$sum
  MLP_out$M <- MLP_out$M / MLP_out$sum
  MLP_out$P <- MLP_out$P / MLP_out$sum
  MLP_out$X <- MLP_out$X / MLP_out$sum
  m <- MLP_out[1, c("S", "M", "P", "X")]

# prob <- unname(unlist(pred_prob_mlp_ind))
diff <- max(m)-max(m[m!=max(m)]) #highest - second highest
RI <- convertRI(diff)

# out <- pred_prob
out <- MLP_out[, c("S", "M", "P", "X")]
colnames(out)[1:4] <- c("Susceptible", "MDR", "Pre-XDR", "XDR")
out[,c("Susceptible", "MDR", "Pre-XDR", "XDR")] <- round(out[,c("Susceptible", "MDR", "Pre-XDR", "XDR")], digits = 4)
if(pred_mlp_ind == "X") {
  out$Class <- "XDR"
} else if(pred_mlp_ind == "P") {
  out$Class <- "Pre-XDR"
} else if(pred_mlp_ind == "M") {
  out$Class <- "MDR"
} else {
  out$Class <- "Susceptible"
}

out$Sample <- Ind_l$SAMPLE
out$RI <- RI
out <- out[c(6, 1:4, 5, 7)]
cat("\n")
cat("Prediction result with full model\n")
cat("\n")
print(out)
cat("\n")

write.table(out, file = paste0(output_path, "prediction.tsv"), row.names = FALSE)

# cat("\n")
# cat("\tSAMPLE =", Ind_l$SAMPLE, "\n")
# cat("\tPredicted class =", suppressWarnings(names(sort(pred_prob_svm_R_ind, decreasing=TRUE)))[1], "\n")
# cat("\tProbability =", max(prob), "\n")
# cat("\tRelability Index (RI) =", RI, "\n\n")


#### SHAP execution with 37 feature model ####

suppressPackageStartupMessages(library(DALEX, quietly = TRUE))

explainer_f <- readRDS(paste0(script_path, "explainer_mlp_37.rds"))
shap_MLP_f <- predict_parts_shap(explainer_f, new_observation = data_2, B = 25)
# print(colnames(shap_MLP_f))

pred_label <- paste0("MLP.", pred_mlp_ind)
shap_MLP_class_f <- shap_MLP_f[shap_MLP_f$label == pred_label & shap_MLP_f$variable_value > 0 & shap_MLP_f$variable_name != "CATEGORY",]
shap_MLP_class_f_agg <- aggregate(shap_MLP_class_f$contribution, by=list(shap_MLP_class_f$variable), FUN=mean)
colnames(shap_MLP_class_f_agg) <- c("variable", "contribution")
shap_MLP_class_f_agg <- shap_MLP_class_f_agg[order(abs(shap_MLP_class_f_agg$contribution), decreasing = TRUE),]

shap_MLP_f_plot <- plot(shap_MLP_f[shap_MLP_f$label == pred_label & shap_MLP_f$variable_value > 0 & shap_MLP_f$variable_name != "CATEGORY",], max_features = 38)

cat("\n")
cat("SHAP result with 37-feature model - ")
cat(Ind_l$SAMPLE)
cat("\n\n")
print(shap_MLP_class_f_agg)
cat("\n")

write.table(shap_MLP_class_f_agg, file = paste0(output_path, "shap_result_37_features_", Ind_l$SAMPLE, ".tsv"), row.names = FALSE)
ggsave(file = paste0(output_path, "shap_plot_37_features_", Ind_l$SAMPLE, ".svg"), plot = shap_MLP_f_plot)


#### SHAP execution with 100 feature model ####

explainer_f <- readRDS(paste0(script_path, "explainer_mlp_100.rds"))
shap_MLP_f <- predict_parts_shap(explainer_f, new_observation = data_2, B = 25)
# print(colnames(shap_MLP_f))

pred_label <- paste0("MLP.", pred_mlp_ind)
shap_MLP_class_f <- shap_MLP_f[shap_MLP_f$label == pred_label & shap_MLP_f$variable_value > 0 & shap_MLP_f$variable_name != "CATEGORY",]
shap_MLP_class_f_agg <- aggregate(shap_MLP_class_f$contribution, by=list(shap_MLP_class_f$variable), FUN=mean)
colnames(shap_MLP_class_f_agg) <- c("variable", "contribution")
shap_MLP_class_f_agg <- shap_MLP_class_f_agg[order(abs(shap_MLP_class_f_agg$contribution), decreasing = TRUE),]

shap_MLP_f_plot <- plot(shap_MLP_f[shap_MLP_f$label == pred_label & shap_MLP_f$variable_value > 0 & shap_MLP_f$variable_name != "CATEGORY",], max_features = 101)

cat("\n")
cat("SHAP result with 100-feature model - ")
cat(Ind_l$SAMPLE)
cat("\n\n")
print(shap_MLP_class_f_agg)
cat("\n")

write.table(shap_MLP_class_f_agg, file = paste0(output_path, "shap_result_100_features_", Ind_l$SAMPLE, ".tsv"), row.names = FALSE)
ggsave(file = paste0(output_path, "shap_plot_100_features_", Ind_l$SAMPLE, ".svg"), plot = shap_MLP_f_plot)
