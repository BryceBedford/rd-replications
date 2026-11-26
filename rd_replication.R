pkgs <- c("haven", "fixest", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(haven)
library(fixest)
library(ggplot2)

data_dir <- "data"
out_dir  <- "output"
dir.create(out_dir, showWarnings = FALSE)

sink(file.path(out_dir, "replication_log.txt"), split = TRUE)

save_plot <- function(p, name, w = 7, h = 4.5) {
  ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 300)
  print(p)
}

show_table <- function(models, name, ...) {
  args <- list(...)
  print(do.call(etable, c(list(models), args)))
  do.call(etable, c(list(models), args,
                    list(tex = TRUE, file = file.path(out_dir, name), replace = TRUE)))
}

lm_raw <- as.data.frame(zap_labels(read_dta(file.path(data_dir, "lmdata.dta"))))
str(lm_raw)
print(summary(lm_raw))
print(table(lm_raw$quebec, lm_raw$year))

lm_raw$census <- factor(1900 + lm_raw$year)
lm_raw$se_cell <- lm_raw$sdempc / sqrt(lm_raw$nobs)
lm_raw$ci_lo   <- lm_raw$empc - 1.96 * lm_raw$se_cell
lm_raw$ci_hi   <- lm_raw$empc + 1.96 * lm_raw$se_cell

cutoff <- 29.5

emp_plot <- function(df, title) {
  ggplot(df, aes(age, empc, colour = census, shape = census)) +
    geom_line(alpha = 0.35) +
    geom_point(size = 2.2) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    labs(x = "Age at Census day", y = "Employment rate in Census week",
         colour = "Census year", shape = "Census year", title = title) +
    theme_bw()
}

save_plot(emp_plot(subset(lm_raw, quebec == 1),
                   "Quebec: employment rate by age, single men, 1986 and 1991"),
          "q3_3_quebec_all_ages.png")

save_plot(emp_plot(subset(lm_raw, quebec == 0),
                   "Rest of Canada: employment rate by age, single men, 1986 and 1991"),
          "q3_4_roc_all_ages.png")

for (qb in c(1, 0)) for (yr in c(86, 91)) {
  s <- subset(lm_raw, quebec == qb & year == yr)
  e <- setNames(s$empc, s$age)
  cat(sprintf("%s 19%d: avg yearly change ages 20-24 = %.3f, ages 25-29 = %.3f, jump 29->30 = %.3f\n",
              ifelse(qb == 1, "Quebec", "RoC   "), yr,
              (e["24"] - e["20"]) / 4, (e["29"] - e["25"]) / 4, e["30"] - e["29"]))
}

save_plot(emp_plot(subset(lm_raw, quebec == 1 & age >= 25 & age <= 39),
                   "Quebec: employment rate by age, ages 25-39"),
          "q3_5_quebec_25_39.png")

qc <- subset(lm_raw, quebec == 1 & year == 86 & age >= 25 & age <= 39)
qc$treat <- as.numeric(qc$age >= 30)
qc$x     <- qc$age - cutoff
qc$x2    <- qc$x^2
qc$x3    <- qc$x^3
qc$tx    <- qc$treat * qc$x
qc$tx2   <- qc$treat * qc$x2
qc$w     <- 1 / (qc$se_cell^2)

cat("\nMean of empc (unweighted, nobs-weighted):",
    round(mean(qc$empc), 3), round(weighted.mean(qc$empc, qc$nobs), 3), "\n")

fml <- list(
  "Linear"           = empc ~ treat + x,
  "Quadratic"        = empc ~ treat + x + x2,
  "Cubic"            = empc ~ treat + x + x2 + x3,
  "Linear spline"    = empc ~ treat + x + tx,
  "Quadratic spline" = empc ~ treat + x + x2 + tx + tx2
)

m_ols <- lapply(fml, function(f) feols(f, data = qc, vcov = "iid"))
m_wls <- lapply(fml, function(f) feols(f, data = qc, weights = ~w, vcov = "iid"))

cat("\nTable 1, column 2 replication: OLS (no weights)\n")
show_table(m_ols, "q3_8_table1_col2_ols.tex", keep = "^treat$",
           headers = names(fml), digits = 3)
cat("\nTable 1, column 2 replication: WLS (weights = 1/V_a, as in the paper)\n")
show_table(m_wls, "q3_8_table1_col2_wls.tex", keep = "^treat$",
           headers = names(fml), digits = 3)

gof <- function(m) {
  stat <- sum(qc$w * resid(m)^2)
  df   <- nobs(m) - length(coef(m))
  c(GOF = stat, df = df, p_value = pchisq(stat, df, lower.tail = FALSE))
}
cat("\nGoodness of fit (paper: linear p = 0.52, linear spline p = 0.72)\n")
print(round(t(sapply(m_wls, gof)), 3))

spline_fit <- m_wls[["Linear spline"]]
qc$fit <- predict(spline_fit, newdata = qc)
fig3 <- ggplot(qc, aes(age, empc)) +
  geom_line(aes(y = ci_lo), linetype = "dotted") +
  geom_line(aes(y = ci_hi), linetype = "dotted") +
  geom_point(size = 2.2) +
  geom_line(data = subset(qc, age < 30),  aes(y = fit), linewidth = 0.8) +
  geom_line(data = subset(qc, age >= 30), aes(y = fit), linewidth = 0.8) +
  geom_vline(xintercept = cutoff, linetype = "dashed") +
  scale_x_continuous(limits = c(24, 40), breaks = 24:40) +
  labs(x = "Age (Census day)", y = "Employment rate",
       title = "Fig. 3 replication: employment rate in Census week, Quebec 1986",
       caption = "Dots: cell means. Dotted lines: 95% confidence bounds. Solid lines: linear spline fit (WLS).") +
  theme_bw()
save_plot(fig3, "q3_7_figure3.png")

ui <- as.data.frame(zap_labels(read_dta(file.path(data_dir, "UI_RD.dta"))))
ui <- ui[, !grepl("^_est_", names(ui))]
str(ui)

num_vars <- names(ui)[sapply(ui, is.numeric)]
sumstats <- t(sapply(ui[num_vars], function(v)
  c(N = sum(!is.na(v)), Mean = mean(v, na.rm = TRUE), SD = sd(v, na.rm = TRUE),
    Min = min(v, na.rm = TRUE), Max = max(v, na.rm = TRUE))))
print(round(sumstats, 3))
write.csv(round(sumstats, 3), file.path(out_dir, "q4_b_summary_stats.csv"))
cat("\nAverage months worked in previous 7 years (expbaseline):",
    round(mean(ui$expbaseline), 2), "\n")

cat("\nPotential duration P by age in years:\n"); print(table(ui$age, ui$P))
cat("Everyone has >= 52 months of experience:", all(ui$expbaseline >= 52), "\n")
cat("Spells with negative nonemployment duration (durnonempdumvar == 0):",
    sum(ui$durnonemp < 0), "\n")
cat("Missing realgdpgrowthf1 (all in 2009):", sum(is.na(ui$realgdpgrowthf1)), "\n")

ols_c <- list(
  "(1) P only"   = feols(durnonemp ~ P, data = ui, vcov = "hetero"),
  "(2) + age"    = feols(durnonemp ~ P + age, data = ui, vcov = "hetero"),
  "(3) + demog." = feols(durnonemp ~ P + age + female + edyrs + nonger,
                         data = ui, vcov = "hetero"),
  "(4) + job"    = feols(durnonemp ~ P + age + female + edyrs + nonger +
                           tenure + expbaseline, data = ui, vcov = "hetero"),
  "(5) + macro"  = feols(durnonemp ~ P + age + female + edyrs + nonger +
                           tenure + expbaseline + unemp + chgunemp + realgdpgrowthf1,
                         data = ui, vcov = "hetero")
)
show_table(ols_c, "q4_c_ols.tex", headers = names(ols_c), digits = 3)

p_by_bin <- aggregate(P ~ agebin, data = ui, FUN = mean)
fig_e <- ggplot(p_by_bin, aes(agebin, P)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = c(42, 44), linetype = "dashed") +
  scale_y_continuous(breaks = c(12, 18, 22)) +
  labs(x = "Age at start of UI spell (15-day bins)", y = "Potential UI duration P (months)",
       title = "Potential UI duration by age") +
  theme_bw()
save_plot(fig_e, "q4_e_P_by_age.png")

rd <- subset(ui, normage < 2)
rd$age42_normage <- rd$age42 * rd$normage
cat("\nRD sample size:", nrow(rd), "\n")
cat("Check: age42 lines up with agebin >= 42 (should be diagonal):\n")
print(table(age42 = rd$age42, right_of_42 = rd$agebin >= 42))
cat("Check: P within the RD sample by side:\n"); print(table(rd$age42, rd$P))

dens <- as.data.frame(table(rd$agebin)); names(dens) <- c("agebin", "n")
dens$agebin <- as.numeric(as.character(dens$agebin))
dens$age42  <- as.numeric(dens$agebin >= 42)
dens$r      <- dens$agebin - 42
fig_f <- ggplot(dens, aes(agebin, n)) +
  geom_point(size = 1.5) +
  geom_smooth(aes(group = age42), method = "lm", formula = y ~ x, se = FALSE, colour = "black") +
  geom_vline(xintercept = 42, linetype = "dashed") +
  labs(x = "Age at start of UI spell (15-day bins)", y = "Number of spells",
       title = "Density of UI spells around the age-42 cutoff") +
  theme_bw()
save_plot(fig_f, "q4_f_density.png")
dens_reg <- feols(n ~ age42 + r + age42:r, data = dens, vcov = "hetero")
cat("\nJump in number of spells per bin at 42 (mean count per bin =",
    round(mean(dens$n), 1), "):\n")
print(coeftable(dens_reg)["age42", ])

rd_reg <- function(y, data = rd)
  feols(as.formula(paste(y, "~ age42 + normage + age42_normage")),
        data = data, cluster = ~agebin)
tab1 <- lapply(c("edyrs", "female", "nonger", "tenure"), rd_reg)
show_table(tab1, "q4_g_table1.tex", keep = "^age42$", digits = 4,
           headers = c("Years of educ.", "Female", "Foreign citizen", "Tenure last job"))

bin_vars <- c("edyrs", "female", "nonger", "tenure", "duruib", "durnonemp")
bins <- aggregate(rd[bin_vars], by = list(agebin = rd$agebin, age42 = rd$age42), FUN = mean)
long <- do.call(rbind, lapply(bin_vars, function(v)
  data.frame(agebin = bins$agebin, age42 = bins$age42, var = v, value = bins[[v]])))

rd_bin_plot <- function(vars, labels, title) {
  df <- subset(long, var %in% vars)
  df$var <- factor(df$var, levels = vars, labels = labels)
  ggplot(df, aes(agebin, value)) +
    geom_point(size = 1.2) +
    geom_smooth(aes(group = age42), method = "lm", formula = y ~ x, se = TRUE,
                colour = "black", linewidth = 0.7) +
    geom_vline(xintercept = 42, linetype = "dashed") +
    facet_wrap(~var, scales = "free_y") +
    labs(x = "Age at start of UI spell (15-day bins)", y = NULL, title = title) +
    theme_bw()
}

save_plot(rd_bin_plot(c("edyrs", "female", "nonger", "tenure"),
                      c("Years of education", "Female", "Foreign citizen", "Tenure at last job"),
                      "Predetermined characteristics around the age-42 cutoff"),
          "q4_h_covariates.png", w = 9, h = 6)

save_plot(rd_bin_plot("duruib", "Months of UI benefit receipt",
                      "Fig. IIIa replication: UI benefit duration by age"),
          "q4_i_fig3a_duruib.png")
save_plot(rd_bin_plot("durnonemp", "Months of nonemployment",
                      "Fig. IIIb replication: nonemployment duration by age"),
          "q4_i_fig3b_durnonemp.png")

rf <- lapply(c("duruib", "durnonemp"), rd_reg)
show_table(rf, "q4_j_table2_rd.tex", keep = "^age42$", digits = 3,
           headers = c("UI duration (B)", "Nonemployment (D)"))

fs <- lm(P ~ age42 + normage + age42_normage, data = rd)
cat("\nFirst stage (effect of crossing 42 on P):", round(coef(fs)["age42"], 4),
    "| residual SD:", round(sd(resid(fs)), 6), "\n")

tsls <- function(y, endog, z, w, cl) {
  W <- cbind(1, as.matrix(w)); Z <- cbind(W, z); X <- cbind(W, endog)
  Xhat <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  b <- solve(crossprod(Xhat, X), crossprod(Xhat, y))
  e <- as.vector(y - X %*% b); n <- nrow(X); k <- ncol(X)
  A <- solve(crossprod(Xhat))
  S <- rowsum(Xhat * e, cl); G <- nrow(S)
  V <- A %*% crossprod(S) %*% A * (G / (G - 1)) * ((n - 1) / (n - k))
  c(Estimate = b[k], `Std. Error` = sqrt(V[k, k]))
}
w_ctrl <- rd[, c("normage", "age42_normage")]

rB <- tsls(rd$duruib,    rd$P, rd$age42, w_ctrl, rd$agebin)
rD <- tsls(rd$durnonemp, rd$P, rd$age42, w_ctrl, rd$agebin)
print(rbind(UI_duration = rB, Nonemployment = rD))
iv_B <- rB["Estimate"]; iv_D <- rD["Estimate"]

iv_db <- feols(durnonemp ~ normage + age42_normage | duruib ~ age42,
               data = rd, cluster = ~agebin)
cat("\nIV of nonemployment on actual UI duration (dD/dB):\n")
print(coeftable(iv_db)["fit_duruib", ])

b_fs  <- coef(fs)["age42"]
b_rfB <- coef(rf[[1]])["age42"]
b_rfD <- coef(rf[[2]])["age42"]
cat(sprintf("\nBy hand: dB/dP = %.3f / %.3f = %.4f | IV output = %.4f\n",
            b_rfB, b_fs, b_rfB / b_fs, iv_B))
cat(sprintf("By hand: dD/dP = %.3f / %.3f = %.4f | IV output = %.4f\n",
            b_rfD, b_fs, b_rfD / b_fs, iv_D))
cat(sprintf("Compare naive OLS coefficient on P (col 1 / col 5): %.4f / %.4f\n",
            coef(ols_c[[1]])["P"], coef(ols_c[[5]])["P"]))

rd_pos <- subset(rd, durnonempdumvar == 1)
rf_pos <- lapply(c("duruib", "durnonemp"), rd_reg, data = rd_pos)
cat("\nRobustness: RD estimates dropping negative-duration spells\n")
show_table(rf_pos, "q4_robust_no_negative.tex", keep = "^age42$", digits = 3,
           headers = c("UI duration (B)", "Nonemployment (D)"))

sink()
cat("Done. Files are in", normalizePath(out_dir), "\n")
