library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
source("refer.ML.R")

Train_data <- read.table("data.train.txt", header = T, sep = "\t", check.names=F, row.names=1, stringsAsFactors=F)
Train_expr=Train_data[,1:(ncol(Train_data)-1),drop=F]
Train_class=Train_data[,ncol(Train_data),drop=F]

Test_data <- read.table("data.test.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors = F)
Test_expr=Test_data[,1:(ncol(Test_data)-1),drop=F]
Test_class=Test_data[,ncol(Test_data),drop=F]
Test_class$Cohort=gsub("(.*)\\_(.*)\\_(.*)", "\\1", row.names(Test_class))
Test_class=Test_class[,c("Cohort", "Type")]


comgene <- intersect(colnames(Train_expr), colnames(Test_expr))
Train_expr <- as.matrix(Train_expr[,comgene])
Test_expr <- as.matrix(Test_expr[,comgene])
Train_set = scaleData(data=Train_expr, centerFlags=T, scaleFlags=T) 
names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort))
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)

methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)
methods=methodRT$Model
methods <- gsub("-| ", "", methods)



classVar = "Type" 
min.selected.var = 3
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))


preTrain.var <- list() 
set.seed(seed = 123)      
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method,              #机器学习方法
                                 Train_set = Train_set,        #训练组的基因表达数据
                                 Train_label = Train_class,    #训练组的分类数据
                                 mode = "Variable",            #选择运行模式(筛选变量)
                                 classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(Train_set)



model <- list()            
set.seed(123)       
Train_set_bk = Train_set

for (method in methods) {
  cat(match(method, methods), ":", method, "\n")
  method_name = method
  method <- strsplit(method, "\\+")[[1]]
  
  if (length(method) == 1) method <- c("simple", method)
  

  if (is.null(preTrain.var[[method[1]]]) || length(preTrain.var[[method[1]]]) == 0) {
    cat("No variables selected for method:", method_name, " - Skipping\n")
    next 
  }
  
  Variable = preTrain.var[[method[1]]]
  
  
  if (!all(Variable %in% colnames(Train_set_bk))) {
    cat("Variable(s) not found in Train_set_bk for method:", method_name, " - Skipping\n")
    next 
  }
  
  Train_set = Train_set_bk[, Variable]
  Train_label = Train_class
  
  model[[method_name]] <- RunML(method = method[2],           
                                Train_set = Train_set,       
                                Train_label = Train_label,  
                                mode = "Model",             
                                classVar = classVar)
  
  if (length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    cat("Selected variables less than threshold for method:", method_name, " - Skipping\n")
    model[[method_name]] <- NULL
    next
  }
}
Train_set = Train_set_bk; rm(Train_set_bk)



#保存所有机器学习模型的结果
saveRDS(model, "model.MLmodel.rds")


#根据基因表达量计算每个样本的分类得分
model <- readRDS("model.MLmodel.rds")
methodsValid <- names(model)                     
common_cols <- intersect(colnames(Train_set), colnames(Test_set))

Train_set_common <- Train_set[, common_cols]
Test_set_common <- Test_set[, common_cols]


combined_data <- rbind.data.frame(Train_set_common, Test_set_common)


RS_list <- list()
for (method in methodsValid) {
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], new_data = combined_data)
}

riskTab=as.data.frame(t(do.call(rbind, RS_list)))
riskTab=cbind(id=row.names(riskTab), riskTab)
write.table(riskTab, "model.riskMatrix.txt", sep="\t", row.names=F, quote=F)


common_cols <- intersect(colnames(Train_set), colnames(Test_set))


Train_set_common <- Train_set[, common_cols]
Test_set_common <- Test_set[, common_cols]


combined_data <- rbind.data.frame(Train_set_common, Test_set_common)


Class_list <- list()
for (method in methodsValid) {
  Class_list[[method]] <- PredictClass(fit = model[[method]], new_data = combined_data)
}


if (length(Class_list) > 0) {
  Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
} else {
  stop("Class_list is empty, no predictions to combine.")
}


if (exists("Class_mat")) {
  classTab <- cbind(id = row.names(Class_mat), Class_mat)
  write.table(classTab, "model.classMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  stop("Class_mat does not exist, cannot proceed with saving results.")
}

classTab=cbind(id=row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep="\t", row.names=F, quote=F)


fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}
fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file="model.genes.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#计算每个模型的AUC值
AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],     
                                Test_set = Test_set,    
                                Test_label = Test_class,    
                                Train_set = Train_set,      
                                Train_label = Train_class, 
                                Train_name = "Train",      
                                cohortVar = "Cohort",     
                                classVar = classVar)      
}
AUC_mat <- do.call(rbind, AUC_list)
aucTab=cbind(Method=row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep="\t", row.names=F, quote=F)



AUC_mat <- read.table("model.AUCmatrix.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)


avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = T)
AUC_mat <- AUC_mat[names(avg_AUC),]

fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))


CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")
names(CohortCol) <- colnames(AUC_mat)


cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = AUC_mat,      
                    avg_Cindex = avg_AUC,       
                    CohortCol = CohortCol,     
                    barCol = "steelblue",      
                    cellwidth = cellwidth, cellheight = cellheight,   
                    cluster_columns = F, cluster_rows = F)      


pdf(file="model.AUCheatmap.pdf", width=cellwidth * ncol(AUC_mat) + 6, height=cellheight * nrow(AUC_mat) * 0.45)
draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()




library(pROC)

rsFile="model.riskMatrix.txt"   
method="plsRglm+RF"              


riskRT=read.table(rsFile, header=T, sep="\t", check.names=F, row.names=1)

CohortID=gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID=gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort=CohortID


for(Cohort in unique(riskRT$Cohort)){

  rt=riskRT[riskRT$Cohort==Cohort,]
  y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
  y=ifelse(y=="Control", 0, 1)
  

  roc1=roc(y, as.numeric(rt[,method]))     
  ci1=ci.auc(roc1, method="bootstrap")     
  ciVec=as.numeric(ci1)
  pdf(file=paste0("ROC.", Cohort, ".pdf"), width=5, height=4.75)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=Cohort)
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}


library(pROC)

# 设置文件路径和方法
rsFile = "model.riskMatrix.txt"  
method = "plsRglm+RF"         

riskRT = read.table(rsFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

CohortID = gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID = gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort = CohortID

for (Cohort in unique(riskRT$Cohort)) {

  rt = riskRT[riskRT$Cohort == Cohort,]
  y = gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
  y = ifelse(y == "Control", 0, 1)
  

  roc1 = roc(y, as.numeric(rt[, method]))  
  ci1 = ci.auc(roc1, method = "bootstrap")
  ciVec = as.numeric(ci1)
  

  pdf(file = paste0("ROC.", Cohort, ".pdf"), width = 5, height = 4.75)
  plot(roc1, print.auc = TRUE, col = "#1f77b4", legacy.axes = TRUE, main = Cohort) # 设置曲线颜色
  

  polygon(c(1, roc1$specificities, 0), c(0, roc1$sensitivities, 0), col = rgb(31/255, 119/255, 180/255, 0.2), border = NA) 
  
  text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "#1f77b4") # 添加置信区间文字
  dev.off()  
}




library(rms)
library(rmda)

inputFile="data.train.txt"       
geneFile="model_gene-1.txt"    


data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)  
row.names(data)=gsub("-", "_", row.names(data))

geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]  

data=t(data) 
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))  
rt=cbind(as.data.frame(data), Type=group) 
paste(colnames(data), collapse="+") 

ddist=datadist(rt)
options(datadist="ddist")

lrmModel=lrm(Type~ RPL36AL+RPS27A+KLF6+HELZ+LCOR+CHD9, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.0001,0.1,0.3,0.5,0.7,0.9,0.99),
              lp=F, funlabel="Risk of Disease")  
pdf("B细胞Nomo.pdf", width=8, height=6)  
plot(nomo)  
dev.off()


cali=calibrate(lrmModel, method="boot", B=1000)  
pdf("B细胞Calibration.pdf", width=5.5, height=5.5)  
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F) 
dev.off()


rt$Type=ifelse(rt$Type=="Control", 0, 1)  
table(rt$Type)


dc=decision_curve(Type ~  RPL36AL+RPS27A+KLF6+HELZ+LCOR+CHD9, data=rt,
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)  
pdf(file="B细胞DCA.pdf", width=5.5, height=5.5)  
plot_decision_curve(dc,
                    curve.names="Model",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="red",
                    confidence.intervals=FALSE,
                    standardize=FALSE)  
dev.off()

######单基因ROC######

library(glmnet)
library(pROC)

# 设置文件路径
expFile="merge.normalize.txt"   
geneFile="model_gene.txt"


rt = read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

y = gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y = ifelse(y == "Control", 0, 1)

geneRT = read.table(geneFile, header=F, sep="\t", check.names=F)

bioCol = rainbow(nrow(geneRT), s=0.9, v=0.9)   
aucText = c()
k = 0
for (x in as.vector(geneRT[,1])) {
  k = k + 1
  roc1 = roc(y, as.numeric(rt[x,]))  
  if (k == 1) {
    pdf(file="ROC.genes.pdf", width=5, height=4.75)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    aucText = c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  } else {
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    aucText = c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  }
}

legend("bottomright", aucText, lwd=2, bty="n", col=bioCol[1:(ncol(rt)-1)])
dev.off()


rt = rt[as.vector(geneRT[,1]),]
rt = as.data.frame(t(rt))
logit = glm(y ~ ., family=binomial(link='logit'), data=rt)
pred = predict(logit, newx=rt)    

roc1 = roc(y, as.numeric(pred))      
ci1 = ci.auc(roc1, method="bootstrap")     
ciVec = as.numeric(ci1)
pdf(file="ROC.model.pdf", width=5, height=4.75)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Model")

polygon(c(roc1$specificities, 1), c(roc1$sensitivities, 0), col=rgb(1, 0, 0, 0.2), border=NA)
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.3f", ciVec[1]), "-", sprintf("%.3f", ciVec[3])), col="red")
dev.off()


##########SHAP############

library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)
library(RColorBrewer)
library(pheatmap)
library(limma)


rt=read.table("merge.normalize.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


gene=read.table("gene.txt", header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]),rownames(data))
geneExp=data[sameGene,]


out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="LASSO-Gexp.txt",sep="\t",quote=F,col.names=F)


out <- cbind(Gene = rownames(geneExp), geneExp)

# 保存为 CSV 文件
write.csv(out, file = "geneexp.csv", row.names = FALSE, quote = FALSE)

work_dir <- "D:\\AD自测数据\\所有数据合并\\1\\B_cell\\SHAP"
result_dir <- "results"  # 结果输出文件夹
random_seed <- 12345
num_colors <- 20  # pastel色系最多9色，超出自动渐变扩展
train_ratio <- 0.7
cv_folds <- 5
roc_file <- file.path(result_dir, "model_roc_curve.pdf")
barplot_file <- file.path(result_dir, "shap_importance_barplot.pdf")
bee_file <- file.path(result_dir, "shap_importance_beeswarm.pdf")
dependence_file <- file.path(result_dir, "shap_dependence_all.pdf")
waterfall_file <- file.path(result_dir, "shap_waterfall_sample1.pdf")
force_file <- file.path(result_dir, "shap_force_sample1.pdf")
heatmap_file <- file.path(result_dir, "shap_value_heatmap.pdf")
top20_file <- file.path(result_dir, "shap_top_genes.csv")
data_file <- "geneexp.csv"
pdf_width <- 6
pdf_height <- 6
progress_style <- 3

ml_methods <- data.frame(
  ModelName = c("RF", "SVM", "XGB", "GBM", "KNN"),
  MethodID = c("rf", "svmRadial", "xgbTree", "gbm", "knn")
)

cat("步骤1/8: 设置工作目录和随机种子...\n")
if (!dir.exists(work_dir)) stop("工作目录不存在，请检查路径！")
setwd(work_dir)
set.seed(random_seed)


if (!dir.exists(result_dir)) dir.create(result_dir)

# ===================== 步骤2：配色 =====================
cat("步骤2/8: 生成可爱风格配色方案...\n")
if (num_colors > 9) {
  color_palette <- colorRampPalette(brewer.pal(9, "Pastel1"))(num_colors)
} else {
  color_palette <- brewer.pal(num_colors, "Pastel1")
}

# ===================== 步骤3：数据读取与预处理 =====================
cat("步骤3/8: 读取并预处理表达数据...\n")
if (!file.exists(data_file)) stop("数据文件不存在！")
raw_data <- read.csv(data_file, header=TRUE, check.names=FALSE, row.names=1)
if (nrow(raw_data) == 0) stop("数据文件为空！")
transposed_data <- t(raw_data)
sample_ids <- rownames(transposed_data)
if (!all(grepl("(_Control$|_Treat$)", sample_ids))) stop("样本名需以_Control或_Treat结尾！")
group_labels <- ifelse(grepl("_Control$", sample_ids), "Control",
                       ifelse(grepl("_Treat$", sample_ids), "Treatment", NA))
if (any(is.na(group_labels))) stop("分组标签生成失败！")
expr_data <- as.data.frame(transposed_data)
expr_data$Group <- as.factor(group_labels)

# ===================== 步骤4：训练集/测试集划分 =====================
cat("步骤4/8: 划分训练集和测试集...\n")
if (length(unique(expr_data$Group)) != 2) stop("分组数不是2，无法二分类！")
split_idx <- createDataPartition(y=expr_data$Group, p=train_ratio, list=FALSE)
train_set <- expr_data[split_idx, ]
test_set <- expr_data[-split_idx, ]
test_labels <- test_set$Group
test_features <- test_set[, -ncol(test_set)]

# ===================== 步骤5：模型训练与评估 =====================
cat("步骤5/8: 多模型训练与ROC评估...\n")
auc_results <- c()
model_auc_map <- list()
roc_colors <- color_palette[1:nrow(ml_methods)]
roc_ci_list <- list()
pb <- txtProgressBar(min=0, max=nrow(ml_methods), style=progress_style)

pdf(file=roc_file, width=pdf_width, height=pdf_height)
par(mar=c(5, 5, 4, 2) + 0.1)
for (i in 1:nrow(ml_methods)) {
  setTxtProgressBar(pb, i)
  mdl_name <- ml_methods$ModelName[i]
  mdl_id <- ml_methods$MethodID[i]
  if (nrow(train_set) < 10) stop("训练集样本太少！")
  if (mdl_id == "svmRadial") {
    mdl <- train(Group ~ ., data=train_set, method=mdl_id, prob.model=TRUE,
                 trControl=trainControl(method="repeatedcv", number=cv_folds, savePredictions=TRUE))
  } else {
    mdl <- train(Group ~ ., data=train_set, method=mdl_id,
                 trControl=trainControl(method="repeatedcv", number=cv_folds, savePredictions=TRUE))
  }
  pred_prob <- predict(mdl, newdata=test_features, type="prob")
  if (!"Treatment" %in% colnames(pred_prob)) stop("预测概率列名不含Treatment！")
  roc_obj <- roc(as.numeric(test_labels)-1, as.numeric(pred_prob[,"Treatment"]))
  auc_val <- as.numeric(roc_obj$auc)
  auc_ci <- ci.auc(roc_obj, conf.level=0.95)
  roc_ci_list[[mdl_id]] <- auc_ci
  auc_results <- c(auc_results,
                   sprintf("%s: %.03f [%.03f, %.03f]", mdl_name, auc_val, auc_ci[1], auc_ci[3]))
  model_auc_map[[mdl_id]] <- auc_val
  if (i == 1) {
    plot(roc_obj, print.auc=FALSE, legacy.axes=TRUE, main="", col=roc_colors[i], lwd=3)
  } else {
    plot(roc_obj, print.auc=FALSE, legacy.axes=TRUE, main="", col=roc_colors[i], lwd=3, add=TRUE)
  }
}
close(pb)
legend('bottomright', auc_results, col=roc_colors, lwd=3, bty='n', cex=0.9)
dev.off()

# ===================== 步骤6：最优模型选择与SHAP计算 =====================
cat("步骤6/8: 选择最优模型并计算SHAP值...\n")
auc_vec <- unlist(model_auc_map)
if (length(auc_vec) == 0) stop("AUC结果为空！")
best_method <- names(which.max(auc_vec))
cat("AUC最高模型为：", best_method, "\n")
final_model <- train(Group ~ ., data=train_set, method=best_method,
                     trControl=trainControl(method="repeatedcv", number=cv_folds, savePredictions=TRUE))
shap_fit <- kernelshap(
  object = final_model,
  X = train_set[, -ncol(train_set)],
  pred_fun = function(obj, newdata) {
    pred <- predict(obj, newdata, type="prob")
    if (!"Treatment" %in% colnames(pred)) stop("SHAP预测概率列名不含Treatment！")
    pred[,"Treatment"]
  }
)

# ===================== 步骤7：SHAP可视化 =====================
cat("步骤7/8: SHAP可视化...\n")
shap_obj <- shapviz(shap_fit, X_pred = train_set[, -ncol(train_set)], X = train_set[, -ncol(train_set)], interactions=TRUE)
shap_importance <- sort(colMeans(abs(shap_obj$S)), decreasing=TRUE)
top_features <- names(shap_importance)

# 统一主题
custom_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right")

options(shapviz.colors = color_palette)

# SHAP重要性条形图
pdf(file=barplot_file, width=pdf_width, height=pdf_height)
print(sv_importance(shap_obj, kind="bar", show_numbers=TRUE) + custom_theme)
dev.off()

# SHAP重要性蜂群图
pdf(file=bee_file, width=pdf_width, height=pdf_height)
print(sv_importance(shap_obj, kind="bee", show_numbers=TRUE) + custom_theme)
dev.off()

# 所有Top特征依赖图
pdf(file=dependence_file, width=12, height=pdf_height)
print(sv_dependence(shap_obj, v=top_features) + custom_theme)
dev.off()

# 瀑布图（第1个样本）
pdf(file=waterfall_file, width=pdf_width, height=pdf_height)
print(sv_waterfall(shap_obj, row_id=1) + custom_theme)
dev.off()

# 力图（第1个样本）
pdf(file=force_file, width=pdf_width, height=pdf_height)
print(sv_force(shap_obj, row_id=1) + custom_theme)
dev.off()

# 针对每个Top基因，分别绘制依赖图、密度图、散点图
cat("为每个Top基因绘制单独的依赖图、密度图、散点图...\n")
for (gene in top_features) {
  safe_gene <- tolower(gsub("[^A-Za-z0-9_]", "_", gene))
  # 依赖图
  pdf(file = file.path(result_dir, paste0("shap_dependence_", safe_gene, ".pdf")), width = pdf_width, height = pdf_height)
  print(sv_dependence(shap_obj, v = gene) + custom_theme + ggtitle(gene))
  dev.off()
  # 密度图
  shap_df <- data.frame(
    SHAP = shap_obj$S[, gene],
    Group = train_set$Group
  )
  pdf(file = file.path(result_dir, paste0("shap_density_", safe_gene, ".pdf")), width=pdf_width, height=pdf_height)
  print(ggplot(shap_df, aes(x=SHAP, fill=Group)) +
          geom_density(alpha=0.5) +
          custom_theme +
          ggtitle(paste("SHAP Density for", gene)))
  dev.off()
  # SHAP-表达量散点图
  shap_df <- data.frame(
    SHAP = shap_obj$S[, gene],
    Expr = train_set[, gene],
    Group = train_set$Group
  )
  pdf(file = file.path(result_dir, paste0("shap_scatter_", safe_gene, ".pdf")), width=pdf_width, height=pdf_height)
  print(ggplot(shap_df, aes(x=Expr, y=SHAP, color=Group)) +
          geom_point(alpha=0.7) +
          custom_theme +
          ggtitle(paste("SHAP vs Expression for", gene)))
  dev.off()
}

# Top2基因交互探索散点图
gene1 <- top_features[1]
gene2 <- top_features[2]
safe_gene1 <- tolower(gsub("[^A-Za-z0-9_]", "_", gene1))
safe_gene2 <- tolower(gsub("[^A-Za-z0-9_]", "_", gene2))
shap_df <- data.frame(
  SHAP = shap_obj$S[, gene1],
  Expr1 = train_set[, gene1],
  Expr2 = train_set[, gene2],
  Group = train_set$Group
)
pdf(file.path(result_dir, paste0("shap_interaction_scatter_", safe_gene1, "_vs_", safe_gene2, ".pdf")), width=pdf_width, height=pdf_height)
print(
  ggplot(shap_df, aes(x=Expr1, y=Expr2, color=SHAP)) +
    geom_point(size=2, alpha=0.8) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
    labs(title=paste("Interaction-like plot:", gene1, "vs", gene2),
         x=gene1, y=gene2, color="SHAP") +
    custom_theme
)
dev.off()

# SHAP聚类热图
pdf(heatmap_file, width=pdf_width, height=pdf_height)
pheatmap(shap_obj$S,
         cluster_rows=TRUE, cluster_cols=TRUE,
         show_rownames=FALSE, show_colnames=TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))
dev.off()

# ===================== 步骤8：输出基因 =====================
cat("步骤8/8: 输出基因...\n")
top20_genes <- head(shap_importance, 20)
write.csv(data.frame(Gene=names(top20_genes), Importance=as.numeric(top20_genes)),
          top20_file, row.names=FALSE)

# SHAP特征累计贡献曲线
cum_contrib <- cumsum(abs(shap_importance)) / sum(abs(shap_importance))
pdf(file.path(result_dir, "shap_cumulative_contribution_curve.pdf"), width = pdf_width, height = pdf_height)
plot(seq_along(cum_contrib), cum_contrib, type = "s", lwd = 2,
     xlab = "Number of top features", ylab = "Cumulative SHAP Contribution",
     main = "Cumulative SHAP Contribution Curve")
abline(h = 0.8, col = "red", lty = 2)  # 80%参考线
dev.off()

#分组SHAP重要性条形图
group_levels <- levels(train_set$Group)
shap_group <- sapply(group_levels, function(grp)
  colMeans(abs(shap_obj$S[train_set$Group == grp, , drop = FALSE]))
)

#分组SHAP重要性条形图
shap_group <- as.data.frame(shap_group)
shap_group$Gene <- rownames(shap_group)
library(reshape2)
shap_group_long <- reshape2::melt(shap_group, id.vars = "Gene", variable.name = "Group", value.name = "SHAP_Importance")

pdf(file.path(result_dir, "shap_group_importance_barplot.pdf"), width = max(8, pdf_width*1.5), height = pdf_height)
print(
  ggplot(shap_group_long, aes(x = reorder(Gene, SHAP_Importance), y = SHAP_Importance, fill = Group)) +
    geom_bar(stat='identity', position = "dodge") +
    coord_flip() +
    labs(title="Group-wise SHAP Feature Importance", x="Gene", y="Mean(|SHAP|)") +
    custom_theme
)
dev.off()
feature_importance_df <- data.frame(
  Gene = names(shap_importance),
  MeanAbsSHAP = as.numeric(shap_importance)
)
write.csv(feature_importance_df, file.path(result_dir, "shap_feature_importance.csv"), row.names = FALSE)
group_levels <- levels(train_set$Group)
shap_group <- sapply(group_levels, function(grp)
  colMeans(abs(shap_obj$S[train_set$Group == grp, , drop = FALSE]))
)
shap_group_df <- data.frame(Gene = rownames(shap_group), shap_group)
write.csv(shap_group_df, file.path(result_dir, "shap_groupwise_importance.csv"), row.names = FALSE)

auc_table <- data.frame(
  Model = ml_methods$ModelName,
  AUC = sapply(ml_methods$MethodID, function(id) as.numeric(model_auc_map[[id]])),
  CI_lower = sapply(ml_methods$MethodID, function(id) roc_ci_list[[id]][1]),
  CI_upper = sapply(ml_methods$MethodID, function(id) roc_ci_list[[id]][3])
)
write.csv(auc_table, file.path(result_dir, "model_auc_summary.csv"), row.names = FALSE)

sample_influence <- rowSums(abs(shap_obj$S))
influential_samples <- order(sample_influence, decreasing=TRUE)[1:10]
influence_df <- data.frame(
  Sample = rownames(train_set)[influential_samples],
  InfluenceScore = sample_influence[influential_samples]
)
write.csv(influence_df, file.path(result_dir, "influential_samples.csv"), row.names=FALSE)

# 将y转为0/1
y_numeric <- ifelse(train_set$Group == "Treatment", 1, 0)
library(DALEX)
explainer <- explain(
  final_model,
  data = train_set[, -ncol(train_set)],
  y = y_numeric,
  label = "BestModel"
)
vi <- model_parts(explainer, loss_function = loss_one_minus_auc)
vi_df <- as.data.frame(vi)
pdf(file.path(result_dir, "feature_permutation_importance.pdf"), width=pdf_width*1.2, height=pdf_height)
print(
  ggplot(vi_df[vi_df$variable != "_full_model_", ], aes(x=reorder(variable, -dropout_loss), y=dropout_loss)) +
    geom_bar(stat="identity", fill="tomato") +
    coord_flip() +
    labs(title="Permutation Feature Importance (AUC Drop)", x="Feature", y="AUC Drop") +
    custom_theme
)
dev.off()
write.csv(vi_df, file.path(result_dir, "feature_permutation_importance.csv"), row.names=FALSE)



cat("全部流程已完成！\n")




######pd-ftd-ob###########
# 加载所需的包
library(glmnet)
library(pROC)

# 设置文件路径
expFile="5fd.normalize.txt"      # 表达数据文件
geneFile="model_gene.txt"         # 基因列表文件(!!!!!!!!需要去掉id)
#setwd("D:\\AD自测数据\\所有数据合并\\1\\1\\5xfad\\9")    # 设置工作目录

# 读取表达数据文件
rt = read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

# 提取样本的分组信息
y = gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y = ifelse(y == "Control", 0, 1)

# 读取基因列表文件
geneRT = read.table(geneFile, header=F, sep="\t", check.names=F)

# 循环遍历每个基因并生成ROC曲线
bioCol = rainbow(nrow(geneRT), s=0.9, v=0.9)    # 设置曲线的颜色
aucText = c()
k = 0
for (x in as.vector(geneRT[,1])) {
  k = k + 1
  # 生成ROC曲线
  roc1 = roc(y, as.numeric(rt[x,]))     # 获取ROC曲线的参数
  if (k == 1) {
    pdf(file="ROC.genes.pdf", width=5, height=4.75)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    aucText = c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  } else {
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    aucText = c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  }
}
# 添加图例，显示每条ROC曲线对应的基因和AUC值
legend("bottomright", aucText, lwd=2, bty="n", col=bioCol[1:(ncol(rt)-1)])
dev.off()

# 生成逻辑回归模型
rt = rt[as.vector(geneRT[,1]),]
rt = as.data.frame(t(rt))
logit = glm(y ~ ., family=binomial(link='logit'), data=rt)
pred = predict(logit, newx=rt)     # 获取模型的预测值

# 计算并绘制模型的ROC曲线
roc1 = roc(y, as.numeric(pred))      # 获取模型ROC曲线的参数
ci1 = ci.auc(roc1, method="bootstrap")     # 获取ROC曲线下面积的置信区间
ciVec = as.numeric(ci1)
pdf(file="ROC.model.pdf", width=5, height=4.75)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Model")
# 填充ROC曲线下的面积
polygon(c(roc1$specificities, 1), c(roc1$sensitivities, 0), col=rgb(1, 0, 0, 0.2), border=NA)
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.3f", ciVec[1]), "-", sprintf("%.3f", ciVec[3])), col="red")
dev.off()

#######UpSetR#########

library(UpSetR)
outFile="intersectGenes.txt"        #输出交集基因文件
outPic="upset.pdf"                  #输出图片
setwd("D:\\AD自测数据\\所有数据合并\\1\\upsetR")      #设置工作目录

files=dir()                          #获取目录下所有文件
files=grep("txt$",files,value=T)     #提取.txt结尾的文件
geneList=list()

#获取所有txt文件中的基因信息，保存到geneList
for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile==outFile){next}
  rt=read.table(inputFile,header=F)         #读取输入文件
  geneNames=as.vector(rt[,1])               #提取基因名称
  geneNames=gsub("^ | $","",geneNames)      #去掉基因首尾的空格
  uniqGene=unique(geneNames)                #基因取unique，唯一基因列表
  header=unlist(strsplit(inputFile,"\\.|\\-"))
  geneList[[header[1]]]=uniqGene
  uniqLength=length(uniqGene)
  print(paste(header[1],uniqLength,sep=" "))
}

#绘制UpSetͼ
upsetData=fromList(geneList)
pdf(file=outPic,onefile = FALSE,width=9,height=6)
upset(upsetData,
      nsets = length(geneList),               #展示多少个数据
      nintersects = 50,                       #展示基因集数目
      order.by = "freq",                      #按照数目排序
      show.numbers = "yes",                   #柱状图上方是否显示数值
      number.angles = 20,                     #字体角度
      point.size = 2,                         #点的大小
      matrix.color="red",                     #交集点颜色
      line.size = 0.8,                        #线条粗细
      mainbar.y.label = "Gene Intersections", 
      sets.x.label = "Set Size")
dev.off()

#保存交集文件
intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)
