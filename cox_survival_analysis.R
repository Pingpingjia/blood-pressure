# ============================================================
# Cox比例风险模型与生存曲线分析
# 包含：虚拟数据生成、Cox模型、KM曲线、模型诊断
# ============================================================

# ---- 1. 安装/加载所需包 ----
packages <- c("survival", "survminer", "ggplot2", "dplyr", "broom")
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

set.seed(42)

# ============================================================
# 2. 生成虚拟数据
# ============================================================
n <- 500  # 样本量

# 协变量
age       <- rnorm(n, mean = 60, sd = 10)
sex       <- rbinom(n, 1, 0.5)           # 0=女, 1=男
sbp       <- rnorm(n, mean = 130, sd = 20)  # 收缩压
treatment <- rbinom(n, 1, 0.5)           # 0=对照, 1=治疗

# 线性预测子（真实系数）
lp <- 0.03 * age + 0.4 * sex + 0.015 * sbp - 0.6 * treatment

# 从Weibull分布生成生存时间
scale <- exp(-lp)
time_true <- rweibull(n, shape = 1.5, scale = scale * 10)

# 随机删失（约30%删失率）
censor_time <- runif(n, min = 1, max = 15)
time   <- pmin(time_true, censor_time)
status <- as.integer(time_true <= censor_time)  # 1=事件, 0=删失

# 构建数据框
dat <- data.frame(
  id        = 1:n,
  time      = time,
  status    = status,
  age       = age,
  sex       = factor(sex, levels = c(0, 1), labels = c("Female", "Male")),
  sbp       = sbp,
  treatment = factor(treatment, levels = c(0, 1), labels = c("Control", "Treatment"))
)

cat("--- 数据概览 ---\n")
cat(sprintf("样本量: %d | 事件数: %d | 删失率: %.1f%%\n",
            n, sum(dat$status), 100 * mean(dat$status == 0)))
cat(sprintf("随访时间: 中位 %.2f, 范围 [%.2f, %.2f]\n\n",
            median(dat$time), min(dat$time), max(dat$time)))

# ============================================================
# 3. Kaplan-Meier 生存曲线
# ============================================================

# 3a. 整体KM曲线
km_overall <- survfit(Surv(time, status) ~ 1, data = dat)

# 3b. 按治疗组分层KM曲线
km_trt <- survfit(Surv(time, status) ~ treatment, data = dat)

# Log-rank检验
logrank_test <- survdiff(Surv(time, status) ~ treatment, data = dat)
logrank_p    <- 1 - pchisq(logrank_test$chisq, df = 1)
cat(sprintf("Log-rank检验 p值: %.4f\n\n", logrank_p))

# 绘制KM曲线
p_km <- ggsurvplot(
  km_trt,
  data          = dat,
  pval          = TRUE,
  pval.method   = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.height = 0.25,
  xlab          = "时间（年）",
  ylab          = "生存概率",
  title         = "Kaplan-Meier 生存曲线（按治疗组）",
  legend.labs   = c("对照组", "治疗组"),
  palette       = c("#E7298A", "#1B9E77"),
  ggtheme       = theme_bw()
)

# 保存KM图
ggsave("km_curve.pdf", plot = print(p_km), width = 8, height = 7)
cat("KM曲线已保存至 km_curve.pdf\n")

# ============================================================
# 4. Cox 比例风险模型
# ============================================================

# 4a. 单因素Cox
univariate_vars <- c("age", "sex", "sbp", "treatment")
uni_results <- lapply(univariate_vars, function(var) {
  formula <- as.formula(paste("Surv(time, status) ~", var))
  fit     <- coxph(formula, data = dat)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    mutate(variable = var)
}) |> bind_rows()

cat("--- 单因素Cox结果 (HR, 95%CI, p) ---\n")
print(uni_results[, c("variable", "term", "estimate", "conf.low", "conf.high", "p.value")],
      digits = 3)

# 4b. 多因素Cox
cox_multi <- coxph(
  Surv(time, status) ~ age + sex + sbp + treatment,
  data = dat
)

cat("\n--- 多因素Cox模型摘要 ---\n")
print(summary(cox_multi))

# 整洁格式输出
cox_tbl <- tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE)
cat("\n--- HR及95%CI ---\n")
print(cox_tbl[, c("term", "estimate", "conf.low", "conf.high", "p.value")], digits = 3)

# ============================================================
# 5. 模型诊断
# ============================================================

# 5a. Schoenfeld残差检验（PH假设）
ph_test <- cox.zph(cox_multi)
cat("\n--- PH假设检验 (Schoenfeld残差) ---\n")
print(ph_test)

# 绘制Schoenfeld残差图
pdf("schoenfeld_residuals.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
plot(ph_test)
dev.off()
cat("Schoenfeld残差图已保存至 schoenfeld_residuals.pdf\n")

# 5b. Martingale残差（识别非线性）
dat$martingale_resid <- residuals(cox_multi, type = "martingale")

p_mart <- ggplot(dat, aes(x = age, y = martingale_resid)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Martingale残差 vs 年龄",
       x = "年龄", y = "Martingale残差") +
  theme_bw()

ggsave("martingale_residuals.pdf", plot = p_mart, width = 6, height = 5)
cat("Martingale残差图已保存至 martingale_residuals.pdf\n")

# 5c. Concordance index (C-statistic)
c_stat <- concordance(cox_multi)
cat(sprintf("\nC-statistic (区分度): %.3f (se=%.3f)\n",
            c_stat$concordance, sqrt(c_stat$var)))

# ============================================================
# 6. 调整后生存曲线（基于Cox模型）
# ============================================================

# 构建新预测数据（固定协变量在均值/参考水平，只改变treatment）
new_data <- data.frame(
  age       = rep(mean(dat$age), 2),
  sex       = factor(c("Female", "Female"), levels = levels(dat$sex)),
  sbp       = rep(mean(dat$sbp), 2),
  treatment = factor(c("Control", "Treatment"), levels = levels(dat$treatment))
)

adj_surv <- survfit(cox_multi, newdata = new_data)

p_adj <- ggsurvplot(
  adj_surv,
  data       = new_data,
  conf.int   = TRUE,
  xlab       = "时间（年）",
  ylab       = "生存概率",
  title      = "Cox模型调整后生存曲线\n（年龄=均值, 女性, SBP=均值）",
  legend.labs = c("对照组", "治疗组"),
  palette    = c("#E7298A", "#1B9E77"),
  ggtheme    = theme_bw()
)

ggsave("cox_adjusted_survival.pdf", plot = print(p_adj), width = 7, height = 6)
cat("调整后生存曲线已保存至 cox_adjusted_survival.pdf\n")

# ============================================================
# 7. Forest Plot（多因素Cox结果）
# ============================================================

p_forest <- ggplot(cox_tbl, aes(x = estimate, y = term,
                                 xmin = conf.low, xmax = conf.high)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbarh(height = 0.2, color = "steelblue") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "多因素Cox模型 Forest Plot",
       x = "风险比 HR (95% CI)",
       y = "") +
  theme_bw()

ggsave("cox_forest_plot.pdf", plot = p_forest, width = 7, height = 5)
cat("Forest Plot已保存至 cox_forest_plot.pdf\n")

cat("\n--- 所有分析完成 ---\n")
