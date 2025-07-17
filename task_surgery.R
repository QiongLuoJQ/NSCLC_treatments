

# 肺癌疾病数据提取【单列】
# 数据预处理与样本筛选【单列】
# 分类变量onehot编码，没有连续变量【不用标准化】
####设置task######
rm(list = ls())
load("./SEER_lungtreat/1data/data_8_lungtreat_fs.Rdata")
load("./SEER_lungtreat/1data/data_12_lungtreat_fs.Rdata")
library(tidyverse)

## 提取需要的变量
y = "Surg_prim"
y_target = "Surgery"
Feature = c("Age", "Gender", "Race_ethnicity", "Martial", "Income", "Rural", "Site",
            "Tumor_Histology", "Tumor_side", "Stage", "T", "N", "M", "Lung_nodule")
target_data = data_8_lungtreat_fs %>% 
  dplyr::select(Patient_ID, y, Feature)

## 转换成任务
library(mlr3verse)
task_f <- function(taskdata,taskname,taks_target) {
  tasks = mlr3::TaskClassif$new(taskname, 
                                backend = taskdata, 
                                target = taks_target,
                                positive = y_target)
  return(tasks)
}
## 数据划分
target_task = task_f(target_data, "target", y)
set.seed(180615)
split_task <- partition(target_task, ratio = 0.7)
target_task_train <- target_task$clone()$filter(split_task$train)
target_task_test <- target_task$clone()$filter(split_task$test)

target_data_train = target_task_train$data() %>% select(Patient_ID, y, Feature)
target_data_test = target_task_test$data() %>% select(Patient_ID, y, Feature)

## 进行onenot编码
dummy_f <- function(data_treat) {
  library(caret)
  dummies <- dummyVars(~., data = data_treat)
  dummy_data <- predict(dummies, newdata = data_treat) |> as.data.frame()
  return(dummy_data)
}
X_train_dummy = dummy_f(target_data_train %>% select(-Patient_ID, -y))
X_test_dummy = dummy_f(target_data_test %>% select(-Patient_ID, -y))
y_train = target_data_train %>% dplyr::select(y)
y_test = target_data_test %>% dplyr::select(y)

### 组合起来，形成新的task
task_train = task_f(cbind(y_train, X_train_dummy), "train", y)
task_test = task_f(cbind(y_test, X_test_dummy), "test", y)

## 外部验证数据
valid_data = data_12_lungtreat_fs %>% 
  filter(!Patient_ID %in% intersect(data_8_lungtreat_fs$Patient_ID,
                                     data_12_lungtreat_fs$Patient_ID)) %>% 
  dplyr::select(colnames(target_data)) 
X_valid_dummy = dummy_f(valid_data %>% select(-Patient_ID, -y))
y_valid = valid_data %>% dplyr::select(y)
task_valid = task_f(cbind(y_valid, X_valid_dummy), "valid", y)

## 保存现有的数据  
path0 <- "./SEER_lungtreat/3result_features/00dataset/"
if(!file.exists(path0)){
  dir.create(path0)
}
write.csv(data_8_lungtreat_fs, row.names = F,
          file = paste0(path0, "lungtreat_multilabels.csv"))
save(target_data_train, target_data_test, valid_data,
     X_train_dummy, X_test_dummy,
     y_train, y_test,
     task_train, task_test, task_valid,
     file = paste0(path0,y,"_task_list.Rdata"))
write.csv(target_data_train, row.names = F,
          file = paste0(path0, y, "_data_train.csv"))
write.csv(target_data_test, row.names = F,
          file = paste0(path0, y, "_data_test.csv"))
####table1分析######
# table 1 对比结果
## 纳入排除对比【单独】
## 对比train和test
### 设置路径，保存结果，01table
path1 <- "./SEER_lungtreat/3result_features/01table/"
if(!file.exists(path1)){
  dir.create(path1)
}
table_1_f <- function(dataset, group) {
  library(gtsummary)
  table_res = dataset |> 
    tbl_summary(
      by = group,
      type = list(
        #### 设置一些需要的指标，转换格式 
      ),
      statistic = list(all_continuous() ~ "{mean} ({sd})",
                       all_categorical() ~ "{n} ({p}%)"
      ),
      percent = "column",
      digits = list(all_continuous()~c(2,2),
                    all_categorical()~c(0,2)),
      missing = "no") |>
    add_difference(everything() ~ "smd") %>%
    modify_header(
      estimate ~ "**SMD**", # rename Difference to SMD
      all_stat_cols() ~ "**{level}**  \nN = {n}"
    ) %>%
    modify_column_hide(columns = ci) %>% # remove 95% CI
    add_p(pvalue_fun = ~style_pvalue(.x, digits = 3),simulate.p.value = TRUE) |>
    add_significance_stars() |>
    add_overall() 
  return(table_res)
}

train_test_table = rbind(target_data_train %>% mutate(group = "train"),
                         target_data_test %>% mutate(group = "test")) %>% 
  select(-Patient_ID) %>% table_1_f(group = group) %>% as.data.frame()
openxlsx::write.xlsx(train_test_table, paste0(path1,y,"_train_test_table.xlsx"))

## 对比develop和valid
dev_valid_table = rbind(target_data %>% mutate(group = "develop"),
                        valid_data %>% mutate(group = "valid")) %>% 
  select(-Patient_ID) %>% table_1_f(group = group) %>% as.data.frame()
openxlsx::write.xlsx(dev_valid_table, paste0(path1,y,"_dev_valid_table.xlsx"))

## train 内部 y对比
train_table = target_data_train %>% select(-Patient_ID) %>% 
  table_1_f(group = y) %>% as.data.frame()
openxlsx::write.xlsx(train_table, paste0(path1,y,"_train_table.xlsx"))

## develop 内部 y 对比
develop_table = target_data %>% select(-Patient_ID) %>% 
  table_1_f(group = y) %>% as.data.frame()
openxlsx::write.xlsx(develop_table, paste0(path1,y,"_develop_table.xlsx"))

####特征选择######
# rm(list = ls())
load(paste0(path0,y,"_task_list.Rdata"))
path2 <- "./SEER_lungtreat/3result_features/02featureselection/"
if(!file.exists(path2)){ dir.create(path2)}
# lasso
### 可以有两种选择方式，一种是变量筛选，一种是特征筛选
#### 如果X_data是多分类，需要对X_data进行数值转换
#### 需要注意的是，单纯的多分类变量。onehot性能会更好。
# X_data = target_data_train %>% select(-Patient_ID, -y) %>%
#   mutate(across(where(is.factor), ~ as.integer(factor(., levels = unique(.)))))
set.seed(180615)
X_data = X_train_dummy %>%
  select(-all_of(findCorrelation(cor(X_train_dummy),
                                 cutoff = 0.9, names = TRUE))) %>%
  as.data.table()

y_data = ifelse(y_train[[y]] == y_target, 1, 0)

X_data_lasso <- function(X_data, y_data) {
  library(glmnet)
  X_data <- as.matrix(X_data)
  
  #### Lasso筛选变量动态过程图
  set.seed(180615)
  la.md <- glmnet(X_data, y_data, family = "gaussian", 
                  intercept = F, alpha = 1, 
                  nlambda = 100) 
  pdf(paste0(path2, y, "lasso_lambda.pdf"), width = 6, height = 6)
  plot(la.md, xvar = "lambda", label = TRUE)
  dev.off()
  
  set.seed(180615)
  mod_cv <- cv.glmnet(x = X_data, y = y_data, 
                      family = "gaussian", 
                      nfolds = 10,
                      intercept = F, alpha = 1)
  
  pdf(paste0(path2, y, "lasso_cv.pdf"), width = 6, height = 6)
  plot(mod_cv)
  abline(v = log(c(mod_cv$lambda.min,mod_cv$lambda.1se)), lty = "dashed")
  dev.off()
  
  min_lambda <- mod_cv$lambda.min ##选择合适的namida值 0.004353376
  se_lambda <- mod_cv$lambda.1se
  ##这里并没有选择fit1的结果进行变量筛选
  # lasso_coef <- coef(mod_cv, s = "lambda.min")
  lasso_coef <- coef(mod_cv, s = "lambda.1se")
  index <- which(lasso_coef != 0)
  actCoef <- lasso_coef[index]
  lasso_var <- row.names(lasso_coef)[index]
  
  varCoef <- as.data.frame(cbind(Var = lasso_var,Coef = actCoef))
  varCoef$Coef <- as.numeric(varCoef$Coef)
  varCoef$positve = ifelse(varCoef$Coef > 0,"positive","negative")
  
  lasso_plot <- ggplot(aes(x = reorder(lasso_var,actCoef),
                           y = actCoef,fill = positve),
                       data = varCoef) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(x = "") +
    ggtitle("LASSO identified variables") +
    theme(legend.position = "")
  pdf(paste0(path2, y, "lasso_imp_1se.pdf"), width = 6, height = 6)
  print(lasso_plot)
  dev.off()
  
  lasso_plot_abs <- ggplot(aes(x = reorder(lasso_var,abs(Coef)),
                               y = abs(Coef),fill = positve),
                           data = varCoef) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(x = "") +
    ggtitle("LASSO identified variables") +
    theme(legend.position = "")
  pdf(paste0(path2, y, "lasso_imp_1se_abs.pdf"), width = 6, height = 6)
  print(lasso_plot_abs)
  dev.off()
  
  train_lasso = cbind(y_train, X_data %>% 
                        as.data.frame() %>% 
                        select(lasso_var)) %>% as.data.table()
  return(train_lasso)
}
train_lasso = X_data_lasso(X_data, y_data)
# train_lasso
save(train_lasso,
     file = paste0(path2, y, "_data_lasso.Rdata"))
write.csv(train_lasso,file = paste0(path2, y, "data_lasso.csv"))

# RFECV
library(mlr3filters)
library(mlr3learners)
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")
## 加速计算
library(future)
plan("multisession",workers = 2)
task_train_lasso = mlr3::TaskClassif$new("train_lasso", 
                                         backend = train_lasso, 
                                         target = y,
                                         positive = y_target)

fs_rfecv_lgb <- function(tasks) {
  fs_rfecv = fs("rfecv",
                n_features = 1,
                feature_number = 1)
  lrn_lgb = lrn("classif.lightgbm",predict_type = "prob")
  set.seed(180615)
  instance <- fselect(fselector = fs_rfecv,
                      task = tasks,
                      learner = lrn_lgb,
                      resampling = rsmp("cv",folds = 10),
                      measures = msr("classif.auc"),
                      terminator = trm("none"),
                      store_models = T )
  
  print(instance$result_feature_set)
  print(instance$result_y)
  return(instance)
}
fs_res = fs_rfecv_lgb(task_train_lasso)

## 提取trainlasso RFECV之后的结果
task_train_lasso$select(fs_res$result_feature_set)

## 保存结果
save(fs_res, 
     file = paste0(path2, y, "_fs_res.Rdata"))
save(task_train_lasso,
     file = paste0(path2, y, "_task_rfe.Rdata"))

#############RFECV plot
res_auc_f <- function(fsdata) {
  fs_res = as.data.table(fsdata$archive)[!is.na(iteration), ]
  # nmax = ncol(fs_res)
  # fs_res$sum = apply(fs_res %>% select(-nmax), 1, sum)
  return(fs_res)
}

rfecv_res = res_auc_f(fs_res) |> 
  select(n_features, classif.auc) |> 
  mutate(n_features = as.numeric(n_features)) |> 
  group_by(n_features) |> 
  summarise(auc = mean(classif.auc),
            sd = sd(classif.auc)) |> 
  as.data.frame()

library(ggplot2)
# 自己提取数据
library(viridisLite)
library(mlr3misc)
n_max = rfecv_res |> filter(auc == max(auc)) |> select(n_features)
rfecv_plots = rfecv_res |> 
  ggplot(aes(x = n_features, y = auc)) +
  geom_line(
    color = viridis(1, begin = 0.5),
    linewidth = 1) +
  geom_point(
    fill = viridis(1, begin = 0.5),
    shape = 21,
    size = 3,
    stroke = 0.5,
    alpha = 0.8) +
  geom_vline(
    xintercept = n_max$n_features,
    colour = viridis(1, begin = 0.33),
    linetype = 3) +
  xlab("Number of Features") +
  ylab("Mean AUC") +
  scale_x_reverse() +
  theme_minimal()

pdf(paste0(path2, y, "_rfecv.pdf"), width = 6, height = 6)
print(rfecv_plots)
dev.off()

####default benchmark######
# rm(list = ls())
path3 <- "./SEER_lungtreat/3result_features/03default_benchmark/"
if(!file.exists(path3)){ dir.create(path3)}
load(paste0("./SEER_lungtreat/3result_features/02featureselection/",y,"_task_rfe.Rdata"))

library(mlr3verse)
library(tidyverse)
lrns = list(
  lrn("classif.log_reg", predict_type = "prob", id = "LogisticRegression"), #二分类适合
  lrn("classif.AdaBoostM1", predict_type = "prob", id = "AdaBoost"),
  # lrn("classif.catboost", predict_type = "prob", id = "CatBoost"),
  lrn("classif.glmnet", predict_type = "prob", s = 0.01, id = "ElasticNet"),
  lrn("classif.kknn", predict_type = "prob", id = "Kknn"),
  lrn("classif.lightgbm", predict_type = "prob", id = "LightGBM"),
  lrn("classif.naive_bayes", predict_type = "prob", id = "Naive_Bayes"),
  # lrn("classif.nnet", predict_type = "prob", id = "ANN"),
  lrn("classif.ranger", predict_type = "prob", id = "RandomForest"),
  lrn("classif.rpart", predict_type = "prob", id = "DecisionTree"),
  # lrn("classif.svm", type = "C-classification",
  #     kernel = "radial", predict_type = "prob", id = "SVM"),
  lrn("classif.xgboost", predict_type = "prob", id = "XGBoost")
)
lrns_default = lrns

msr_classif = c("classif.acc","classif.auc",
                "classif.bacc","classif.bbrier",
                "classif.prauc","classif.fbeta",
                "classif.ce",#"classif.costs",
                "classif.fdr","classif.mcc",
                "classif.fnr","classif.fpr",
                "classif.tnr","classif.tpr",
                "classif.npv","classif.ppv",
                "classif.precision","classif.recall",
                "classif.sensitivity","classif.specificity")

design <- benchmark_grid(task_train_lasso, lrns_default, rsmp("cv", folds = 10))
# lgr::get_logger("mlr3")$set_threshold("warn")
# lgr::get_logger("bbotk")$set_threshold("warn")

set.seed(180615)
bmr <- benchmark(design,store_models = T)
bmr

benchmark_res = bmr$aggregate(msrs(msr_classif))[,c("task_id","learner_id",..msr_classif)]
benchmark_res_table = cbind(benchmark_res[,c(1,2)],benchmark_res[,-c(1:2)] %>% round(3))
benchmark_res_table
openxlsx::write.xlsx(benchmark_res_table,
                     paste0(path3,y,"_benchmark_default.xlsx"))
# 基于benchmark的测试与验证
# rm(list = ls())
load(paste0("./SEER_lungtreat/3result_features/00dataset/",y,"_task_list.Rdata"))
load(paste0("./SEER_lungtreat/3result_features/02featureselection/",y,"_task_rfe.Rdata"))
###保留原始变量的分析
# task_test_valid_f <- function(target_data,task_train_lasso) {
#   X_data_test = target_data %>% select(-Patient_ID, -y) %>%
#     mutate(across(where(is.factor), ~ as.integer(factor(., levels = unique(.)))))
#   y_data_test = target_data %>% select(y)
#   task_lasso = mlr3::TaskClassif$new("test_lasso", 
#                                      backend = cbind(y_data_test,X_data_test) %>% 
#                                        select(y, task_train_lasso$feature_names),
#                                      target = y, positive = y_target) 
#   return(task_lasso)
# }

# task_test_lasso = task_test_valid_f(target_data_test, task_train_lasso)
# task_valid_lasso = task_test_valid_f(valid_data, task_train_lasso)
### dummy变量的分析
task_test_lasso = task_test
task_valid_lasso = task_valid
task_test_lasso$select(task_train_lasso$feature_names)
task_valid_lasso$select(task_train_lasso$feature_names)

test_valid_task <- function(lrns, task_train, task_test,
                            validname, boot, n, path) {
  
  classif_total <- vector("list")
  classif_pre_all <- vector("list")
  for(k in 1:length(lrns)){
    print(k)
    classif_lrn = lrns[[k]]
    
    set.seed(180615)
    classif_lrn$train(task_train)
    print(classif_lrn)
    
    if(boot == TRUE){
      classif_res_all <- vector("list", n)
      set.seed(180615)
      task_rsmp <- rsmp("bootstrap",repeats = n)
      task_rsmp$instantiate(task_test)
      
      for(j in 1:n){
        print(j)
        test_boot = task_test$data()[task_rsmp$train_set(j),]
        test_boot_task = mlr3::TaskClassif$new("bootstrap", 
                                               backend = test_boot, 
                                               target = y,
                                               positive = y_target)
        set.seed(180615)
        classif_pre = classif_lrn$predict(test_boot_task)
        classif_res = classif_pre$score(msrs(msr_classif)) %>% round(4)
        print(classif_res)
        classif_res_all[[j]] <- classif_res
      }
      classif_df <- do.call(rbind, classif_res_all) %>% 
        as.data.frame() %>% 
        mutate(lrn = lrns[[k]]$id)
      classif_total[[lrns[[k]]$id]] = classif_df
      classif_pre_all[[lrns[[k]]$id]] = classif_pre
    } else {
      set.seed(180615)
      classif_pre = classif_lrn$predict(task_test)
      classif_res = classif_pre$score(msrs(msr_classif)) %>% round(4)
      print(classif_res)
      
      classif_total[[lrns[[k]]$id]] = c(lrn = lrns[[k]]$id, classif_res)
      classif_pre_all[[lrns[[k]]$id]] = classif_pre
    }
  }
  save(classif_pre_all, file = paste0(path, y,"_", validname ,"_predict.Rdata"))
  
  save(classif_total, file = paste0(path, y,"_", validname ,"_res.Rdata"))
  classif_total_res = do.call(rbind, classif_total) %>% as.data.frame()
  openxlsx::write.xlsx(classif_total_res,
                       paste0(path,y,"_",validname,"_res.xlsx"))
  return(classif_pre_all)
}

test_res = test_valid_task(lrns_default, task_train_lasso, task_test_lasso,
                           "test", boot = FALSE, n = 100, path3)
valid_res = test_valid_task(lrns_default, task_train_lasso, task_valid_lasso,
                            "valid", boot = FALSE, n = 100, path3)

data_train_rfe = task_train_lasso$data()
data_test_rfe = task_test$data()
write.csv(data_train_rfe, row.names = F,
          file = paste0(path0, y, "_data_train_rfe.csv"))
write.csv(data_test_rfe, row.names = F,
          file = paste0(path0, y, "_data_test_rfe.csv"))

####lrn tuning######
# 基于特征选择模型超参数调优【随机搜索】
library(mlr3tuning)
library(mlr3tuningspaces)
path4 <- "./SEER_lungtreat/3result_features/04tuning/"
if(!file.exists(path4)){ dir.create(path4)}
# 设置调优策略和参数搜索空间
lrns_tuning  = lrns[-1]

search_classif <- list(
  ps(#adaboost
    P = p_int(90, 100), #N Estimators 弱学习最大数量
    S = p_int(1, 100), # Min Samples Leaf (S) 最小样本数
    I = p_int(1, 100) # Max Depth (I)最大深度
  ),
  # ps(#catboost
  #   depth = p_int(1, 16),
  #   iterations = p_int(1, 1000),
  #   learning_rate = p_dbl(0.01, 0.3)
  # ),
  ps(#glmnet
    s = p_dbl(0, 1),
    alpha = p_dbl(0, 1)
  ),
  ps(#kknn
    k = p_int(1, 20, default = 7),
    distance = p_int(1, 5, default = 2)
  ),
  ps(#lgb
    learning_rate = p_dbl(0.01, 0.3),
    max_depth = p_int(1, 100),
    num_iterations = p_int(1,1000)
  ),
  ps(#naive_bayes
    laplace = p_dbl(0, 1),
    threshold = p_dbl(0.001, 1)
  ),
  # ps(#nnet
  #   size = p_int(1, 10),
  #   decay = p_dbl(0, 0.1),
  #   maxit = p_int(100, 1000)
  # ),
  ps(#ranger
    mtry.ratio = p_dbl(0.1, 1),
    num.trees = p_int(100, 1000),
    sample.fraction = p_dbl(0.1, 1.0, default = 0.5)
  ),
  ps(#rpart
    minsplit = p_int(1, 20),
    cp = p_dbl(0, 1),
    minbucket = p_int(1, 64)
  ),
  # ps(#svm
  #   cost = p_dbl(0.1, 10),
  #   gamma = p_dbl(0, 5)
  # ),
  ps(#xgboost
    max_depth = p_int(1, 10),
    nrounds = p_int(1, 1000),
    eta = p_dbl(0.01, 0.3)
  )
)
instance_tune = vector("list", length(lrns_tuning))

for(i in 1:length(lrns_tuning)){
  set.seed(180615)
  instance <- tune(
    tuner = tnr("grid_search", resolution = 5),
    # tuner = tnr("random_search"),
    # terminator = trm("evals", n_evals = 100),
    task = task_train_lasso,
    learner = lrns_tuning[[i]],
    resampling = rsmp("cv",folds = 5),
    measures = msr("classif.auc"),
    search_space = search_classif[[i]]
  )
  instance$result_learner_param_vals
  lrns_tuning[[i]]$param_set$values <- instance$result_learner_param_vals
  instance_tune[[i]] = instance
}
save(instance_tune, lrns_tuning, file = paste0(path4, y, "_tune.Rdata"))

####lrn valid######
# 根据tuning之后的模型，再测试验证一遍。
path5 <- "./SEER_lungtreat/3result_features/05tuning_benchmark/"
if(!file.exists(path5)){ dir.create(path5)}
# lrns_all = c(lrns[[1]], lrns_tuning)
#### 如果是跨系统调配数据，那么需要重新对这些lrn进行设置parameter
lrns_all = lrns
for(i in 2:length(lrns)){
  lrns_all[[i]]$param_set$values <- lrns_tuning[[i-1]]$param_set$values
}

design_tuning <- benchmark_grid(task_train_lasso, lrns_all, rsmp("cv", folds = 10))
# lgr::get_logger("mlr3")$set_threshold("warn")
# lgr::get_logger("bbotk")$set_threshold("warn")

set.seed(180615)
bmr_tuning <- benchmark(design_tuning,store_models = T)
bmr_tuning

benchmark_res_tuning = bmr_tuning$aggregate(msrs(msr_classif))[,c("task_id","learner_id",..msr_classif)]
benchmark_res_table_tuning = cbind(benchmark_res_tuning[,c(1,2)],
                                   benchmark_res_tuning[,-c(1:2)] %>% 
                                     round(3))
benchmark_res_table_tuning
openxlsx::write.xlsx(benchmark_res_table_tuning,
                     paste0(path5,y,"_benchmark_tuning.xlsx"))

test_res_tune = test_valid_task(lrns_all, task_train_lasso, task_test_lasso,
                                "test_tune", boot = FALSE, n = 100, path5)
valid_res_tune = test_valid_task(lrns_all, task_train_lasso, task_valid_lasso,
                                 "valid_tune", boot = FALSE, n = 100, path5)
test_res_boot = test_valid_task(lrns_all, task_train_lasso, task_test_lasso,
                                "test_boot", boot = TRUE, n = 100, path5)
valid_res_boot = test_valid_task(lrns_all, task_train_lasso, task_valid_lasso,
                                 "valid_boot", boot = TRUE, n = 100, path5)

###########################
path5 <- "./SEER_lungtreat/3result_features/05tuning_benchmark/"

boot_res_f <- function(path5, y, dataname) {
  
  load(paste0(path5, y, "_",dataname, "_boot_res.Rdata"))
  library(tidyverse)
  classif_res = do.call(rbind, classif_total) %>% as.data.frame()
  
  classif_res_long = classif_res %>% pivot_longer(cols = !c("lrn"),
                                                  names_to = "metrics",
                                                  values_to = "value")
  
  # 定义计算平均值和CI的函数
  mean_ci <- function(x){
    mean_x <- mean(x, na.rm = TRUE)
    n <- length(x)
    se <- sd(x, na.rm = TRUE) / sqrt(n)
    ci <- qt(0.975, df = n - 1) * se
    list(mean = mean_x, lower = mean_x - ci, upper = mean_x + ci)
  }
  
  # 计算每个组合的CI
  classif_all <- classif_res_long %>%
    group_by(lrn, metrics) %>% 
    summarise(
      mean = mean_ci(value)[[1]],
      lower = mean_ci(value)[[2]],
      upper = mean_ci(value)[[3]]
    ) %>%
    mutate(ci = paste0(sprintf("%0.3f (%0.3f, %0.3f)",
                               mean, lower, upper))) %>% 
    dplyr::select(lrn, metrics, ci) %>% 
    pivot_wider(id_cols = lrn,
                names_from = metrics,
                values_from = ci)
  
  openxlsx::write.xlsx(classif_all, 
                       file = paste0(path5, y,"_",dataname,"_boot_res.xlsx"))

}
boot_res_f(path5, y,"test")
boot_res_f(path5, y,"valid")

# ############Stack性能检测########
# ###额外增加，Stacking模型
# ##选择最佳的lightgbm，glmnet和xgboost作为基础学习者，lr作为meta学习者
# # 1,3,5,9
# pre_f <- function(lrn, train_task, test_task, valid_task) {
#   set.seed(180615)
#   # 定义5折交叉验证的重采样方案
#   rcv = rsmp("repeated_cv", repeats = 5, folds = 5)
#   # 使用交叉验证训练模型
#   rr = resample(train_task, lrn, rcv)
#   # rr$aggregate(msr("surv.cindex"))
#   rrp_train = data.table(id = rr$prediction()$row_ids,
#                          prob = rr$prediction()$prob[,2]) %>%
#     group_by(id) %>% summarise(prob = mean(prob))
#   
#   set.seed(180615)
#   lrn$train(train_task)
#   rrp_pre = lrn$predict(test_task)
#   rrp_test = data.table(id = rrp_pre$row_ids,
#                         prob = rrp_pre$prob[,2])
#   
#   valid_pre = lrn$predict(valid_task)
#   rrp_valid = data.table(id = valid_pre$row_ids,
#                          prob = valid_pre$prob[,2])
#   return(list(rrp_train,rrp_test,rrp_valid))
# }
# 
# lgb_pre = pre_f(lrns_all[[5]], 
#                 task_train_lasso, 
#                 task_test_lasso,
#                 task_valid_lasso)
# lgb_train = lgb_pre[[1]]
# lgb_test = lgb_pre[[2]]
# lgb_valid = lgb_pre[[3]]
# 
# glm_pre = pre_f(lrns_all[[3]], 
#                 task_train_lasso, 
#                 task_test_lasso, 
#                 task_valid_lasso)
# glm_train = glm_pre[[1]]
# glm_test = glm_pre[[2]]
# glm_valid = glm_pre[[3]]
# 
# xgb_pre = pre_f(lrns_all[[9]], 
#                 task_train_lasso, 
#                 task_test_lasso,
#                 task_valid_lasso)
# xgb_train = xgb_pre[[1]]
# xgb_test = xgb_pre[[2]]
# xgb_valid = xgb_pre[[3]]
# ##计算后面就是软投票。如果直接算分类，就是硬投票
# ##只需要三个概率就是普通stack，纳入原始特征就是特征组合。
# train_ensemble = cbind(task_train_lasso$data(),
#                        lgb = lgb_train$prob,
#                        glm = glm_train$prob,
#                        xgb = xgb_train$prob) %>% 
#   as.data.table() |> 
#   mutate(prob = (lgb+glm+xgb) / 3) |> 
#   select(y, prob)
# task_train_ensemble = mlr3::TaskClassif$new("train_ensemble", 
#                                             backend = train_ensemble,
#                                             target = y, positive = y_target)
# test_ensemble = cbind(task_test_lasso$data(),
#                       lgb = lgb_test$prob,
#                       glm = glm_test$prob,
#                       xgb = xgb_test$prob) %>% 
#   as.data.table() |> 
#   mutate(prob = (lgb+glm+xgb) / 3) |> 
#   select(y, prob)
# task_test_ensemble = mlr3::TaskClassif$new("test_ensemble", 
#                                            backend = test_ensemble,
#                                            target = y, positive = y_target)
# 
# valid_ensemble = cbind(task_valid_lasso$data(),
#                        lgb = lgb_valid$prob,
#                        glm = glm_valid$prob,
#                        xgb = xgb_valid$prob) %>% 
#   as.data.table() |> 
#   mutate(prob = (lgb+glm+xgb) / 3) |> 
#   select(y, prob)
# task_valid_ensemble = mlr3::TaskClassif$new("valid_ensemble", 
#                                             backend = valid_ensemble,
#                                             target = y, positive = y_target)
# 
# lrn_meta = lrns_all[[1]]
# set.seed((180615))
# lrn_meta$train(task_train_ensemble)
# 
# stack_pre = lrn_meta$predict(task_test_ensemble)
# classif_res = stack_pre$score(msrs(msr_classif)) %>% round(3)
# print(classif_res)
# 
# stack_pre_valid = lrn_meta$predict(task_valid_ensemble)
# classif_res_valid = stack_pre_valid$score(msrs(msr_classif)) %>% round(3)
# print(classif_res_valid)
# #####
# ####查看结果显示，并没有提升性能。所以还是作罢。


# 模型ROC曲线绘制
# 模型校准曲线绘制
# 模型DCA曲线绘制
#########
Curves_f <- function(lrns, classif_pre_all,validname, path5) {
  # 模型ROC曲线绘制
  lrns_all = lrns
  level_names = lrns_all %>% sapply(function(x) x$id)
  
  library(pROC)
  roc_list = vector("list", length(lrns_all))
  roc_res = vector("list", length(lrns_all))
  
  for(ri in 1:length(lrns_all)){
    roc_list[[ri]] = roc(classif_pre_all[[lrns_all[[ri]]$id]]$truth,
                         as.data.frame(classif_pre_all[[lrns_all[[rj]]$id]]$prob)[[y_target]])
    names(roc_list)[ri] = lrns_all[[ri]]$id
    
    roc_res[[ri]] = cbind(sensitivity = roc_list[[ri]]$sensitivities,
                          specificity = roc_list[[ri]]$specificities,
                          models = names(roc_list)[ri],
                          auc = sprintf("%.3f",roc_list[[ri]]$auc))
  }
  
  roc_res_total = do.call(rbind, roc_res) %>% 
    as.data.frame() %>% 
    mutate(Sensitivity  = as.numeric(sensitivity),
           Specificity = 1 - as.numeric(specificity),
           Model = paste0(models," (AUC = ",auc,")"))
  # roc_res_total$Model = factor(roc_res_total$Model, levels = level_names)
  
  roc_plots = roc_res_total %>% 
    ggplot(aes(x = Specificity, 
               y = Sensitivity,
               color = Model)) +
    geom_path(lwd = 0.5, alpha = 0.8)+
    geom_abline(lty = 3) +
    theme_minimal() +
    theme_classic() +
    # scale_color_lancet() +
    theme(legend.background = element_rect(fill = rgb(1,1,1,alpha = 0.001),
                                           colour = NA)) +
    scale_x_continuous(name = "1 - Specificity", limits = c(0,1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    theme(legend.position = c(0.8,0.3)) +
    theme(panel.border = element_rect(fill = NA,
                                      color = "black", 
                                      linewidth = 0.5, 
                                      linetype = "solid"))
  pdf(paste0(path5, y, "_",validname,"_tune_roc.pdf"), width = 6, height = 6)
  print(roc_plots)
  dev.off()
  
  ###
  delong_list = expand.grid(lrni = sapply(lrns_all, function(learner) learner$id),
                            lrnj = sapply(lrns_all, function(learner) learner$id)) %>% 
    filter(lrni != lrnj)
  for(ti in 1:nrow(delong_list)){
    li = delong_list$lrni[ti]
    lj = delong_list$lrnj[ti]
    delong_test = roc.test(classif_pre_all[[li]]$truth,
                           as.data.frame(classif_pre_all[[lrns_all[[li]]$id]]$prob)[[y_target]],
                           as.data.frame(classif_pre_all[[lrns_all[[lj]]$id]]$prob)[[y_target]])
    # print(delong_test)
    delong_list$auc1[ti] = sprintf("%.3f",delong_test$roc1$auc) 
    delong_list$auc2[ti] = sprintf("%.3f",delong_test$roc2$auc)
    delong_list$Statistic[ti] = sprintf("%.3f",delong_test$statistic)
    delong_list$p.value[ti] = ifelse(delong_test$p.value < 0.001,
                                     "<0.001",
                                     sprintf("%.3f",delong_test$p.value))
  }
  delong_res = delong_list %>% 
    mutate(lrni = factor(lrni, levels = lrns_all %>% sapply(function(x) x$id)),
           lrnj = factor(lrnj, levels = lrns_all %>% sapply(function(x) x$id))) %>% 
    select(lrni, lrnj, p.value) %>% 
    pivot_wider(id_cols = lrni,
                names_from = lrnj, 
                values_from = p.value) %>% 
    arrange(lrni)
  delong_stat = delong_list %>% 
    mutate(lrni = factor(lrni, levels = lrns_all %>% sapply(function(x) x$id)),
           lrnj = factor(lrnj, levels = lrns_all %>% sapply(function(x) x$id))) %>% 
    select(lrni, lrnj, Statistic) %>% 
    pivot_wider(id_cols = lrni,
                names_from = lrnj, 
                values_from = Statistic) %>% 
    arrange(lrni)
  
  openxlsx::write.xlsx(list(delong_res, delong_stat, delong_list),
                       paste0(path5,y,"_",validname,"_delong.xlsx"))
  
  ####PRAUC 
  library(PRROC)
  prroc_list = vector("list", length(lrns_all))
  prroc_res = vector("list", length(lrns_all))
  
  for(pri in 1:length(lrns_all)){
    prroci = classif_pre_all[[lrns_all[[pri]]$id]]$data
    prroc_list[[pri]] = pr.curve(scores.class0 = as.data.frame(prroci$prob)[[y_target]],
                                 weights.class0 = prroci$truth == y_target, 
                                 curve = T)
    names(prroc_list)[pri] = lrns_all[[pri]]$id
    
    prroc_res[[pri]] = cbind(precision = prroc_list[[pri]]$curve[,1],
                             recall = prroc_list[[pri]]$curve[,2],
                             models = names(prroc_list)[pri],
                             prauc = sprintf("%.3f",prroc_list[[pri]]$auc.integral))
  }
  
  prroc_res_total = do.call(rbind, prroc_res) %>% 
    as.data.frame() %>% 
    mutate(Precision  = as.numeric(precision),
           Recall = as.numeric(recall),
           Model = paste0(models," (PRAUC = ",prauc,")"))
  # roc_res_total$Model = factor(roc_res_total$Model, levels = level_names)
  
  prroc_plots = prroc_res_total %>% 
    ggplot(aes(x = Precision, 
               y = Recall,
               color = Model)) +
    geom_path(lwd = 0.5, alpha = 0.8) +
    # geom_abline(lty = 3) +
    theme_minimal() +
    theme_classic() +
    # scale_color_lancet() +
    theme(legend.background = element_rect(fill = rgb(1,1,1,alpha = 0.001),
                                           colour = NA)) +
    scale_x_continuous(name = "Precision", limits = c(0,1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    theme(legend.position = c(0.8,0.3)) +
    theme(panel.border = element_rect(fill = NA,
                                      color = "black", 
                                      linewidth = 0.5, 
                                      linetype = "solid"))
  pdf(paste0(path5, y, "_",validname,"_tune_prroc.pdf"), width = 6, height = 6)
  print(prroc_plots)
  dev.off()
  
  # 模型校准曲线绘制
  cal_list = vector("list", length(lrns_all))
  cal_res = vector("list", length(lrns_all))
  
  library(probably)
  for(cj in 1:length(lrns_all)){
    cal_list[[cj]] = cbind(id = classif_pre_all[[lrns_all[[cj]]$id]]$row_ids,
                           y = classif_pre_all[[lrns_all[[cj]]$id]]$truth == y_target,
                           x = as.data.frame(classif_pre_all[[lrns_all[[cj]]$id]]$prob)[[y_target]]) %>% 
      as.data.frame() %>% 
      mutate(outcomes = y, predictions = as.numeric(x)) %>% 
      cal_plot_breaks(outcomes, predictions, 
                      num_breaks = 10 , # 选几个点
                      include_rug = F, # 是否添加地毯线
                      include_ribbon = T, # 是否添加可信区间
                      conf_level = 0.95 # 可信区间范围 
      )
    cal_res[[cj]] = cal_list[[cj]]$data %>% 
      mutate(n = row_number(), models = lrns_all[[cj]]$id)
  }
  
  cal_res_total = do.call(rbind, cal_res) %>%  as.data.frame() 
  cal_res_slop = cal_res_total |> 
    mutate(slop = event_rate / predicted_midpoint) |> 
    select(models, n, slop) |> 
    pivot_wider(id_cols = n,
                names_from = models, 
                values_from = slop)
  
  openxlsx::write.xlsx(list(cal_res_slop, cal_res_total),
                       paste0(path5,y,"_",validname,"_cal_slop.xlsx"))
  
  cal_plots = cal_res_total %>% ggplot(aes(x = predicted_midpoint, 
                                           y = event_rate, 
                                           color = models, shape = models))+
    geom_path(lwd = 0.5, alpha = 0.8) +
    geom_abline(lty = 3) +
    theme_minimal() +
    theme_classic() +
    # scale_color_lancet() +
    theme(legend.background = element_rect(fill = rgb(1,1,1,alpha = 0.001),
                                           colour = NA)) +
    xlab("Predicted Probability")+
    ylab("Observed Risk")+
    scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    theme(legend.position = c(0.2,0.7)) +
    theme(panel.border = element_rect(fill = NA,
                                      color = "black", 
                                      linewidth = 0.5, 
                                      linetype = "solid"))
  pdf(paste0(path5, y, "_",validname,"_tune_cal.pdf"), width = 6, height = 6)
  print(cal_plots)
  dev.off()
  
  # 模型DCA曲线绘制
  #DCA决策曲线
  dca_list = vector("list", length(lrns_all))
  for(di in 1:length(lrns_all)){
    dca_list[[di]] = cbind(id = classif_pre_all[[lrns_all[[di]]$id]]$row_ids,
                           y = classif_pre_all[[lrns_all[[di]]$id]]$truth == y_target,
                           x = as.data.frame(classif_pre_all[[lrns_all[[dj]]$id]]$prob)[[y_target]],
                           models = lrns_all[[di]]$id)
  }
  dca_data = do.call(rbind, dca_list) %>% 
    as.data.frame() %>% 
    mutate(value = as.numeric(x)) %>% 
    pivot_wider(id_cols = c(id,y), names_from = models, values_from = value) %>% 
    select(-id) %>% as.data.frame()
  
  library(dcurves)
  dca_res_total = dca_data %>% dca(y == 1 ~ ., data = .)
  dca_res = dca_res_total$dca |> as.data.frame()
  dca_threshold = dca_res |> 
    group_by(label) |> 
    filter(net_benefit > 0) |> 
    filter(threshold == max(threshold))
  
  openxlsx::write.xlsx(list(dca_threshold,dca_res),
                       paste0(path5,y,"_",validname,"_dca_threshold.xlsx"))
  
  dca_plots = dca_res_total  %>% 
    plot(type = 'net_benefit', smooth = FALSE, show_ggplot_code = FALSE) +
    theme_minimal() +
    theme_classic() +
    # scale_color_lancet() +
    theme(legend.background = element_rect(fill = rgb(1,1,1,alpha = 0.001),
                                           colour = NA)) +
    # scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
    # scale_y_continuous(limits = c(-0.2,1), expand = c(0,0)) +
    theme(legend.position = c(0.8,0.7)) +
    theme(panel.border = element_rect(
      fill = NA,color = "black", linewidth = 0.5, 
      linetype = "solid"))
  pdf(paste0(path5, y, "_",validname,"_tune_dca.pdf"), width = 6, height = 6)
  print(dca_plots)
  dev.off()
}

Curves_f(lrns, test_res_tune, "test_tune", path5)
Curves_f(lrns, valid_res_tune, "valid_tune", path5)

# 模型SHAP可解释性分析【这里需要用到原始变量，python实现更合适】
