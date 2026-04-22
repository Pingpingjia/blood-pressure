# 生成Cox模型的代码，并绘制生存曲线，包括虚拟数据的生成
# 加载必要的库
library(survival)
library(survminer)
# 生成虚拟数据
set.seed(123) # 设置随机种子以确保结果可重复
n <- 200 # 样本量
# 生成生存时间和事件状态
survival_time <- rexp(n, rate = 0.1) # 生成生存时间，服从指数分布
event_status <- rbinom(n, 1, 0.7) # 生成事件状态，70%发生事件
# 生成一个二分类变量（如治疗组和对照组）
treatment_group <- rbinom(n, 1, 0.5) # 生成治疗组变量，50%在治疗组
# 创建数据框
data <- data.frame(survival_time, event_status, treatment_group)
# 拟合Cox比例风险模型
cox_model <- coxph(Surv(survival_time, event_status) ~ treatment_group, data = data)
# 输出模型摘要
summary(cox_model)
# 绘制生存曲线
surv_fit <- survfit(cox_model)
ggsurvplot(surv_fit, data = data, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, legend.labs = c("Control", "Treatment"),
           xlab = "Time", ylab = "Survival Probability",
           title = "Cox Proportional Hazards Model Survival Curves")