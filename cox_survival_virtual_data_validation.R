library(survival)

set.seed(123)
n <- 300

covariate <- rnorm(n, mean = 0, sd = 1)
baseline_hazard <- 0.08
beta <- log(1.5)

event_time <- rexp(n, rate = baseline_hazard * exp(beta * covariate))
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
