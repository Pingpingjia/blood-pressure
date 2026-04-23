if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required. Please run install.packages('survival').")
}
library(survival)

set.seed(123)
n <- 300

covariate <- rnorm(n, mean = 0, sd = 1)
# Baseline event rate for the simulated cohort.
baseline_hazard <- 0.08
# log(1.5) corresponds to a hazard ratio of 1.5 per 1 SD increase in covariate.
beta <- log(1.5)

event_time <- rexp(n, rate = baseline_hazard * exp(beta * covariate))
# Independent censoring rate chosen to create a mix of events and censored samples.
censor_time <- rexp(n, rate = 0.05)

time <- pmin(event_time, censor_time)
status <- as.integer(event_time <= censor_time)

data <- data.frame(
  time = time,
  status = status,
  covariate = covariate
)

cox_model <- coxph(Surv(time, status) ~ covariate, data = data)
print(summary(cox_model))

newdata <- data.frame(covariate = c(-1, 1))
surv_curve <- survfit(cox_model, newdata = newdata)

plot(
  surv_curve,
  col = c("steelblue", "tomato"),
  lwd = 2,
  xlab = "Time",
  ylab = "Survival Probability",
  main = "Cox Model Survival Curves (Virtual Data)"
)
legend(
  "bottomleft",
  legend = c("Covariate = -1", "Covariate = 1"),
  col = c("steelblue", "tomato"),
  lwd = 2,
  bty = "n"
)
