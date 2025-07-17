
rm(list = ls())
library(mlr3verse)
library(survival)
library(data.table)
library(tidyverse)

####
#train数据，生存数据
y = "Surg_prim"
path7 = "./SEER_lungtreat/3result_features/07PSM/"
if(!file.exists(path7)){ dir.create(path7)}
# 首先根据target train数据，确定patientID，并与原始数据匹配生存时间
# 之后，将target train数据进行dummy转换，结合train lasso 筛选需要的特征。
load(paste0("./SEER_lungtreat/3result_features/00dataset/",y, "_task_list.Rdata"))
load("./SEER_lungtreat/1data/data_8_lungtreat_fs.RData")

surv_data = target_data_train |> 
  left_join(data_8_lungtreat_fs |> 
              select(Patient_ID, Survival_time, Status,
                     Caus.death, Surg_type),
            by = "Patient_ID")

##########
###train lasso 顺序是一致的。所以直接复制过来。但注意排序必须是一致的。
load(paste0("./SEER_lungtreat/3result_features/04tuning/",y,"_tune.Rdata"))
###直接使用XGBoost就可以。###通过前面的分析 还是觉得 LR模型就可以了。
load(paste0("./SEER_lungtreat/3result_features/02featureselection/",y,"_task_rfe.Rdata"))
# lrn_best = lrn("classif.log_reg", predict_type = "prob", id = "LogisticRegression")
lrn_best = lrn("classif.lightgbm", predict_type = "prob", id = "LightGBM")
# lrn_best$param_set$values <- lrns_tuning[[4]]$param_set$values
##这里lightgbm估计更新了，需要注意paramater单独设置。
lrn_best$param_set$values$learning_rate = lrns_tuning[[4]]$param_set$values$learning_rate
lrn_best$param_set$values$max_depth = lrns_tuning[[4]]$param_set$values$max_depth
lrn_best$param_set$values$num_iterations = lrns_tuning[[4]]$param_set$values$num_iterations

# lrn_best = lrn("classif.xgboost", predict_type = "prob", id = "XGBoost")
# lrn_best$param_set$values <- lrns_tuning[[8]]$param_set$values

set.seed(180615)
# 定义5折交叉验证的重采样方案
rcv = rsmp("repeated_cv", repeats = 5, folds = 5)
# 使用交叉验证训练模型
rr = resample(task_train_lasso, lrn_best, rcv)
rrp = data.table(id = rr$prediction()$row_ids,
                 # truth = rr$prediction()$truth,
                 eps = rr$prediction()$prob[,1]) %>% #注意这里提取的概率。
  group_by(id) %>% summarise(eps = mean(eps))

##################
library(gtsummary)
res_origin = surv_data |> 
  select(-Patient_ID) |> 
  tbl_summary(
    by = y,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_dichotomous() ~ "{p}%")
  ) %>%
  add_difference(everything() ~ "smd") %>%
  modify_header(
    estimate ~ "**SMD**", # rename Difference to SMD
    all_stat_cols() ~ "**{level}**  \nN = {n}"
  ) %>%
  modify_column_hide(columns = ci) # remove 95% CI

openxlsx::write.xlsx(res_origin |> as.data.frame(),
                     file = paste0(path7, y, "_psm_origin.xlsx"))

surv_fs = surv_data |> select(-c(Patient_ID,Caus.death,Surg_type))

coxph(Surv(Survival_time, Status)~., data = surv_fs)

library(MatchIt)
surv_fs$eps = rrp$eps
##必须加上这一步，不然对照组和阳性组兑换了。
surv_fs$Surg_prim = factor(surv_fs$Surg_prim, levels = c("None_surgery", "Surgery"))
################
varlist_origin = c("Age","Gender", "Race_ethnicity", "Martial",
                   "Income", "Rural", "Site", "Tumor_Histology",
                   "Tumor_side", "Stage", "T", "N", "M", "Lung_nodule")
varlist_formula = paste0(varlist_origin, collapse = " + ")
psm_res_f <- function(survdatas, rate, varlist_formula) {
  set.seed(180615)
  surgery_matchit_res <- MatchIt::matchit(
    formula= as.formula(paste0(y, " ~", varlist_formula)),
    data=survdatas,
    # distance = "logit",
    distance = survdatas$eps, # 自己估计的eps
    method='nearest',
    # replace=TRUE,
    # discard='both',
    # ratio=2,
    caliper = rate, # starting typically with 0.25
    std.caliper = TRUE,
    replace = FALSE
  )
  
  # Extract matched data for analysis and rename variable "subclass" to "matched_group_id"
  surgery_post_match <- MatchIt::match.data(surgery_matchit_res) %>%
    rename(matched_group_id=subclass)
  print(dim(surgery_post_match))
  
  res_match = surgery_post_match |> 
    select(y,
           Survival_time,Status,
           Age, Gender,
           Race_ethnicity,
           Martial, Income, Rural,
           Site, 
           Tumor_Histology, Tumor_side, 
           Stage, T, N, M,
           Lung_nodule) |> 
    tbl_summary(
      by = y,
      statistic = list(all_continuous() ~ "{mean} ({sd})",
                       all_dichotomous() ~ "{p}%")
    ) %>%
    add_difference(everything() ~ "smd") %>%
    modify_header(
      estimate ~ "**SMD**", # rename Difference to SMD
      all_stat_cols() ~ "**{level}**  \nN = {n}"
    ) #%>%
    # modify_column_hide(columns = ci) # remove 95% CI
  
  # Cox Regression for marginal HR-Estimating using the cluster robust standard error
  coxMatchedclust <- 
    coxph(Surv(Survival_time, Status)~ Surg_prim2,
          robust = TRUE,
          weights = weights,
          cluster = matched_group_id,
          data = surgery_post_match
    ) %>%
    tbl_regression(exponentiate = TRUE)
  
  
  # Frailty model
  coxMatchedfrail <- 
    coxph(
      Surv(Survival_time, Status)~ Surg_prim2 + frailty(matched_group_id),
      data = surgery_post_match
    ) %>% 
    tbl_regression(exponentiate = TRUE, include = y)
  
  res_cox = list(coxMatchedclust, coxMatchedfrail) %>%
    tbl_merge(tab_spanner = c("**Marginal HR model**", "**Frailty model**"))
  return(list(res_match, res_cox,surgery_post_match))
}
psm_all = psm_res_f(surv_fs, 0.01, varlist_formula)

psm_all_match = psm_all[[1]] |> as.data.frame()
psm_all_cox = psm_all[[2]] |> as.data.frame()
surv_matchdata = psm_all[[3]]

openxlsx::write.xlsx(list(psm_all_match,psm_all_cox),
                     file = paste0(path7, y, "_psm_all.xlsx"))

km_plot_f <- function(surv_matchdata,targets, k) {
  surv_matchdata$Survival_time = surv_matchdata$Survival_time / 12
  # 从Cox模型中获取预测生存曲线
  fit_cox <- survfit(Surv(Survival_time, Status) ~ Surg_prim, 
                     data = surv_matchdata)
  summary(fit_cox)
  
  library(survminer)
  surv_plot = ggsurvplot(fit_cox, data = surv_matchdata, 
                         pval = TRUE,                 
                         conf.int = TRUE,              
                         risk.table = T,
                         risk.table.height = 0.25,
                         pval.size = 5,
                         surv.median.line = "hv",
                         xlab = "Time (years)",
                         y_lab = "Overall Survival Probability",
                         title = ifelse(k == "all", "", k),
                         legend.labs = c("no Surgery", "Surgery"),
                         legend.title = "Treatment",
                         palette = c("red", "blue")
  )
  # 绘制Cox模型的预测生存曲线
  pdf(paste0(path7, y, "_", targets, "_", k, "_survplot.pdf"), width = 6, height = 6)
  print(surv_plot)
  dev.off()
}
km_plot_f(surv_matchdata, "all","all")
km_plot_f(surv_data, "origin","all")
################
psm_res_subtype <- function(surv_fs,target_var, rates) {
  varlist_origin = c("Age","Gender", "Race_ethnicity", "Martial",
                     "Income", "Rural", "Site", "Tumor_Histology",
                     "Tumor_side", "Stage", "T", "N", "M", "Lung_nodule")
  varlist_remove = varlist_origin[-which(varlist_origin == target_var)]
  varlist_formula = paste0(varlist_remove, collapse = " + ")
  
  varlevels = levels(surv_fs[[target_var]])
  
  match_list = list()
  cox_list = list()
  
  for(k in varlevels){
    surv_fs$target_vars = surv_fs[[target_var]]
    surv_subtype = surv_fs |> filter(target_vars == k)
    psm_subtype = psm_res_f(surv_subtype, rates, varlist_formula)
    
    
    match_list[[k]] = psm_subtype[[1]] |> as.data.frame()
    cox_list[[k]] = psm_subtype[[2]] |> as.data.frame()
    surv_type_data = psm_subtype[[3]]
    km_plot_f(surv_type_data, target_var,k)
    
  }
  return(list(match_list, cox_list))
  openxlsx::write.xlsx(match_list,
                       file =paste0(path7, y, "_",target_var,"_match_list.xlsx"))
  openxlsx::write.xlsx(cox_list,
                       file =paste0(path7, y, "_",target_var,"_cox_res.xlsx"))
  
}

stage_psm = psm_res_subtype(surv_fs, target_var = "Stage", 0.01)
m_psm = psm_res_subtype(surv_fs, target_var = "M", 0.01)
n_psm = psm_res_subtype(surv_fs, target_var = "N", 0.01)
age_psm = psm_res_subtype(surv_fs, target_var = "Age", 0.01)

stage_res = do.call(rbind, stage_psm[[2]]) %>% na.omit() 
colnames(stage_res) = c("Characteristic","HR_Marginal","CI_Marginal","P_Marginal",
                        "HR_Frailty","CI_Frailty", "P_Frailty")
stage_res$hr_ci1 = paste0(stage_res$HR_Marginal, " (", stage_res$CI_Marginal, ")")
stage_res = stage_res |> rownames_to_column("Facet")

m_res = do.call(rbind, m_psm[[2]]) %>% na.omit() 
colnames(m_res) = c("Characteristic","HR_Marginal","CI_Marginal","P_Marginal",
                        "HR_Frailty","CI_Frailty", "P_Frailty")
m_res$hr_ci1 = paste0(m_res$HR_Marginal, " (", m_res$CI_Marginal, ")")
m_res = m_res |> rownames_to_column("Facet")

n_res = do.call(rbind, n_psm[[2]]) %>% na.omit()
colnames(n_res) = c("Characteristic","HR_Marginal","CI_Marginal","P_Marginal",
                        "HR_Frailty","CI_Frailty", "P_Frailty")
n_res$hr_ci1 = paste0(n_res$HR_Marginal, " (", n_res$CI_Marginal, ")")
n_res = n_res |> rownames_to_column("Facet")

openxlsx::write.xlsx(list(stage_res, m_res, n_res),
                     file = paste0(path7, y, "_allpsm_res.xlsx"))

age_res = do.call(rbind, age_psm[[2]]) %>% na.omit()
colnames(age_res) = c("Characteristic","HR_Marginal","CI_Marginal","P_Marginal",
                    "HR_Frailty","CI_Frailty", "P_Frailty")
age_res$hr_ci1 = paste0(age_res$HR_Marginal, " (", age_res$CI_Marginal, ")")
age_res = age_res |> rownames_to_column("Facet")

openxlsx::write.xlsx(age_res,
                     file = paste0(path7, y, "_age_psm_res.xlsx"))

