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

Input <- read.table(input_file, sep = "\t", comment.char = "", header = TRUE, stringsAsFactors=FALSE)
Ind_batch <- Input[,c("POS", "REF", "ALT", "QUAL", "SAMPLE", "AO.1", "DP.1")]
Ind_batch$AO.1[Ind_batch$AO.1=="."] <- 0
Ind_batch$DP.1[Ind_batch$DP.1=="."] <- 0
Ind_batch <- transform(Ind_batch, DP.1 = as.numeric(DP.1))
Ind_batch <- transform(Ind_batch, AO.1 = as.numeric(AO.1))
Ind_batch$MUT <- paste(Ind_batch$REF, Ind_batch$POS, Ind_batch$ALT, sep = "")
Ind_batch$VAL <- Ind_batch$AO.1/Ind_batch$DP.1
Ind_batch$VAL[Ind_batch$VAL == "NaN"] <- 0.0000
Ind_batch$VAL <- round(Ind_batch$VAL, digits = 4)
Ind_batch_l <- dcast(Ind_batch[,c("SAMPLE", "MUT", "VAL")], SAMPLE~MUT, value.var = "VAL")

A_col <- unlist(read.table(paste0(script_path, "A_col.txt")))
mlp <- readRDS(paste0(script_path, "MLP_model.rds"))

data_2_bat <- data.frame(matrix(ncol = length(A_col), nrow = nrow(Ind_batch_l)))
colnames(data_2_bat) <- A_col

for(j in 1:nrow(data_2_bat)){
  for(i in colnames(data_2_bat)){
    if(i %in% colnames(Ind_batch_l)){
      data_2_bat[j,i] <- Ind_batch_l[j,i]
    }else{
      data_2_bat[j,i] <- 0.00
    }
  }
}

pred_prob_bat <- predict(mlp, data_2_bat, type="prob")
pred_bat <- predict(mlp, data_2_bat)
# print(colnames(pred_prob_bat))
# print(pred_prob_bat)

MLP_out <- pred_prob_bat
for(i in 1:nrow(pred_prob_bat)){
  row_sum <- rowSums(pred_prob_bat[i, c("S", "M", "P", "X")])
  MLP_out[i, "S"] <- MLP_out[i, "S"] / row_sum
  MLP_out[i, "M"] <- MLP_out[i, "M"] / row_sum
  MLP_out[i, "P"] <- MLP_out[i, "P"] / row_sum
  MLP_out[i, "X"] <- MLP_out[i, "X"] / row_sum
}
# print(MLP_out)

# pred_bat <- unlist(pred_prob_svm_R_ind_bat)
output <- data.frame(matrix(ncol = 7, nrow = nrow(pred_prob_bat)))
colnames(output) <- c("SAMPLE", "MDR", "Susceptible", "Pre-XDR", "XDR", "CLASS", "RI")
for(i in 1:nrow(pred_prob_bat)){
  n <- which.max(pred_prob_bat[i,])
  output[i,"SAMPLE"] <- Ind_batch_l$SAMPLE[i]
  output[i, "MDR"] <- MLP_out[i, "M"]
  output[i, "Pre-XDR"] <- MLP_out[i, "P"]
  output[i, "Susceptible"] <- MLP_out[i, "S"]
  output[i, "XDR"] <- MLP_out[i, "X"]
#   output[i, 2:5] <- pred_prob_bat[i, 1:4]
  output[i, 2:5] <- round(output[i, 2:5], digits = 4)
#   output[i,"CLASS"] <- colnames(pred_prob_bat[i,])[n]
  if(colnames(MLP_out[i,])[n] == "M"){
    output[i,"CLASS"] <- "MDR"
  } else if(colnames(MLP_out[i,])[n] == "S"){
    output[i,"CLASS"] <- "Susceptible"
  } else if(colnames(MLP_out[i,])[n] == "P"){
    output[i,"CLASS"] <- "Pre-XDR"
  } else{
    output[i,"CLASS"] <- "XDR"
  }
#   output[i,3] <- max(pred_prob_bat[i,])- max(pred_prob_bat[i,-n])
#   output[i,"RI"] <- round((max(pred_prob_bat[i,])- max(pred_prob_bat[i,-n]))*10, digits = 0)
  output[i,"RI"] <- convertRI((max(MLP_out[i,])- max(MLP_out[i,-n])))
}
# colnames(output) <- c("SAMPLE", "PREDICTION", "DIFF", "RI")

cat("\n")
cat("Prediction result with full model\n")
cat("\n")
print(output)
cat("\n")

write.table(output, file = paste0(output_path, "prediction.tsv"), row.names = FALSE)

suppressPackageStartupMessages(library(DALEX, quietly = TRUE))

explainer_f <- readRDS(paste0(script_path, "explainer_mlp_37.rds"))
for(i in 1:length(pred_bat)){
  shap_MLP_f <- predict_parts_shap(explainer_f, new_observation = data_2_bat[i,], B = 25)
#   print(colnames(shap_MLP_f))

  pred_label <- paste0("MLP.", pred_bat[i])
  shap_MLP_class_f <- shap_MLP_f[shap_MLP_f$label == pred_label & shap_MLP_f$variable_value > 0 & shap_MLP_f$variable_name != "CATEGORY", ]
  shap_MLP_class_f_agg <- aggregate(shap_MLP_class_f$contribution, by=list(shap_MLP_class_f$variable), FUN=mean)
  colnames(shap_MLP_class_f_agg) <- c("variable", "contribution")
  shap_MLP_class_f_agg <- shap_MLP_class_f_agg[order(abs(shap_MLP_class_f_agg$contribution), decreasing = TRUE),]

  shap_MLP_f_plot <- plot(shap_MLP_f[shap_MLP_f$label == pred_label & shap_MLP_f$variable_value > 0 & shap_MLP_f$variable_name != "CATEGORY",], max_features = 38)

  cat("\n")
  cat("SHAP result with 37-feature model - ")
  cat(Ind_batch_l$SAMPLE[i])
  cat("\n\n")
  print(shap_MLP_class_f_agg)
  cat("\n")

  write.table(shap_MLP_class_f_agg, file = paste0(output_path, "shap_result_37_features_", Ind_batch_l$SAMPLE[i], ".tsv"), row.names = FALSE)
  ggsave(file = paste0(output_path, "shap_plot_37_features_", Ind_batch_l$SAMPLE[i], ".svg"), plot = shap_MLP_f_plot)
}

explainer_f <- readRDS(paste0(script_path, "explainer_mlp_100.rds"))
for(i in 1:length(pred_bat)){
  shap_MLP_f <- predict_parts_shap(explainer_f, new_observation = data_2_bat[i, ], B = 25)
#   print(colnames(shap_MLP_f))

  pred_label <- paste0("MLP.", pred_bat[i])
  shap_MLP_class_f <- shap_MLP_f[shap_MLP_f$label == pred_label & shap_MLP_f$variable_value > 0 & shap_MLP_f$variable_name != "CATEGORY",]
  shap_MLP_class_f_agg <- aggregate(shap_MLP_class_f$contribution, by=list(shap_MLP_class_f$variable), FUN=mean)
  colnames(shap_MLP_class_f_agg) <- c("variable", "contribution")
  shap_MLP_class_f_agg <- shap_MLP_class_f_agg[order(abs(shap_MLP_class_f_agg$contribution), decreasing = TRUE),]

  shap_MLP_f_plot <- plot(shap_MLP_f[shap_MLP_f$label == pred_label & shap_MLP_f$variable_value > 0 & shap_MLP_f$variable_name != "CATEGORY",], max_features = 101)

  cat("\n")
  cat("SHAP result with 100-feature model - ")
  cat(Ind_batch_l$SAMPLE[i])
  cat("\n\n")
  print(shap_MLP_class_f_agg)
  cat("\n")

  write.table(shap_MLP_class_f_agg, file = paste0(output_path, "shap_result_100_features_", Ind_batch_l$SAMPLE[i], ".tsv"), row.names = FALSE)
  ggsave(file = paste0(output_path, "shap_plot_100_features_", Ind_batch_l$SAMPLE[i], ".svg"), plot = shap_MLP_f_plot)
}
