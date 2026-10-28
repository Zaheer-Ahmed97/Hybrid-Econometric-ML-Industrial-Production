# =============================================================================
#  "Do hybrid econometric-machine learning models improve forecasts of
#   industrial production under commodity price shocks? A leakage-free
#   evaluation for Pakistan"
#  INPUT : the three Excel files from the replication repository, placed in DATA_DIR
#            Industrial Production Data.xlsx, Wheat data.xlsx, Energy data.xlsx

#  RUNTIME: Sections 2-4 take about one minute. The rolling-origin forecasts
#  (Section 5) take roughly 20-25 minutes on 8 cores and 45-60 minutes on 4
#  core. #
#  Tested with R 4.3.3. Random seeds are fixed, so the numbers are reproducible.
# =============================================================================
DATA_DIR <- "C:/Users/ahmed/Desktop/IfADo Projects/Kashif Thesis new/Revised_manuscript_IJF_LaTeX_with_figures"  # <-- folder that contains the three Excel files
OUT_DIR  <- "C:/Users/ahmed/Desktop/IfADo Projects/Kashif Thesis new/Revised_manuscript_IJF_LaTeX_with_figures/results"      # <-- all tables and figures are written here
N_CORES  <- max(1, parallel::detectCores() - 1)   # cores for the rolling forecasts
B_BOOT   <- 1999           # bootstrap replications for the bounds tests 
B_MCS    <- 5000           # bootstrap replications for the model confidence set 

pkgs <- c("readxl","dplyr","tidyr","zoo","urca","tseries","lmtest","sandwich","strucchange",
          "FinTS","ARDL","forecast","gbm","randomForest","MCS","ggplot2","parallel")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Please install the missing packages first:\n  install.packages(c(",
                       paste0('"', miss, '"', collapse = ", "), "))")
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(zoo); library(urca); library(tseries)
  library(lmtest); library(sandwich); library(strucchange); library(FinTS); library(ARDL)
  library(forecast); library(gbm); library(randomForest); library(MCS); library(ggplot2); library(parallel)
})
dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "cache"), showWarnings = FALSE)

banner  <- function(x) cat("\n", strrep("=", 78), "\n ", x, "\n", strrep("=", 78), "\n", sep = "")
save_tab <- function(df, file) write.csv(df, file.path(OUT_DIR, file), row.names = FALSE)
r3 <- function(df) { df[] <- lapply(df, function(v) if (is.numeric(v)) round(v, 3) else v); df }

# =============================================================================
#  SECTION 2. DATA AND PRELIMINARY EVIDENCE
# =============================================================================
parse_mon_yr <- function(x) {
  mon_map <- c(Jan="01",Feb="02",Mar="03",Apr="04",May="05",Jun="06",
               Jul="07",Aug="08",Sep="09",Oct="10",Nov="11",Dec="12")
  x <- trimws(as.character(x)); out <- as.Date(rep(NA, length(x)))
  for (i in seq_along(x)) {
    if (is.na(x[i]) || nchar(x[i]) == 0) next
    parts <- strsplit(x[i], "-")[[1]]; if (length(parts) != 2) next
    mon <- mon_map[substr(parts[1], 1, 3)]; if (is.na(mon)) next
    out[i] <- as.Date(paste0("20", parts[2], "-", mon, "-01"))
  }
  out
}
# QIM (PBS, base 2015-16)
ip_raw  <- read_excel(file.path(DATA_DIR, "Industrial Production Data.xlsx"), sheet = "Trend Sheet (2)", skip = 5, col_names = FALSE)
ip_data <- data.frame(date = as.Date(format(as.Date(ip_raw[[1]]), "%Y-%m-01")), QIM = as.numeric(ip_raw[[2]]))
ip_data <- ip_data[!is.na(ip_data$date) & !is.na(ip_data$QIM), ]; ip_data <- ip_data[!duplicated(ip_data$date), ]
# Wheat (PBS CPI average prices, 10 kg, 2015-16 base block)
wh_raw  <- read_excel(file.path(DATA_DIR, "Wheat data.xlsx"), sheet = "Sheet1", skip = 8, col_names = FALSE)
wh_data <- data.frame(date = parse_mon_yr(wh_raw[[5]]), WHEAT = as.numeric(wh_raw[[6]]))
wh_data <- wh_data[!is.na(wh_data$date) & !is.na(wh_data$WHEAT), ]; wh_data <- wh_data[!duplicated(wh_data$date), ]
# Energy (PBS CPI average prices, 2007-08 base series)
en_raw  <- read_excel(file.path(DATA_DIR, "Energy data.xlsx"), sheet = "Base Year 2007-08", skip = 1, col_names = FALSE)
colnames(en_raw) <- c("Month","Electricity","Gas_MMBTU","Petrol","Diesel","CNG")
en_data <- data.frame(date = parse_mon_yr(en_raw$Month), Electricity = as.numeric(en_raw$Electricity),
                      Gas_MMBTU = as.numeric(en_raw$Gas_MMBTU), Petrol = as.numeric(en_raw$Petrol),
                      Diesel = as.numeric(en_raw$Diesel))
en_data <- en_data[!is.na(en_data$date), ]; en_data <- en_data[!duplicated(en_data$date), ]

merged <- ip_data %>% inner_join(wh_data, by = "date") %>%
  inner_join(en_data[, c("date","Petrol","Diesel","Gas_MMBTU","Electricity")], by = "date")
merged$FUEL   <- (merged$Petrol + merged$Diesel) / 2
merged$ENERGY <- merged$Electricity
merged$GAS    <- merged$Gas_MMBTU
merged <- merged[order(merged$date), ]
for (v in c("QIM","WHEAT","FUEL","ENERGY","GAS")) merged[[v]] <- as.numeric(na.approx(merged[[v]], na.rm = FALSE))
d  <- as.data.frame(merged[complete.cases(merged[, c("QIM","WHEAT","FUEL","ENERGY","GAS")]),
                           c("date","QIM","WHEAT","FUEL","ENERGY","GAS")])
N  <- nrow(d); y <- d$QIM; PR <- c("WHEAT","FUEL","ENERGY","GAS"); X <- as.matrix(d[, PR])
mon <- as.numeric(format(d$date, "%m"))

# ---- Table 1 --------------------------------------------------------------
banner("Table 1. Series definitions, sources and construction")
tab1 <- data.frame(
  Series = c("QIM","Wheat","Fuel","Electricity","Natural gas"),
  Source = c("PBS, Quantum Index of Manufacturing, base 2015-16",
             "PBS, CPI Prices for Selected Months, Items and Cities: wheat (10 kg), Pakistan average",
             "PBS, CPI Prices for Selected Months, Items and Cities: petrol and high-speed diesel",
             "PBS, CPI Prices for Selected Months, Items and Cities: electricity per unit",
             "PBS, CPI Prices for Selected Months, Items and Cities: natural gas per MMBTU"),
  Unit = c("Index","PKR/10 kg","PKR/litre","PKR/kWh","PKR/MMBTU"),
  Construction = c("Monthly level, not seasonally adjusted","Monthly national average","(petrol + diesel)/2",
                   "Monthly value as recorded","Monthly value as recorded"))
print(tab1, right = FALSE); save_tab(tab1, "Table01_sources.csv")
cat(sprintf("\nMerged sample: %s to %s, N = %d (no missing values after merging: %s)\n",
            format(min(d$date), "%B %Y"), format(max(d$date), "%B %Y"), N, !anyNA(d)))

# ---- Table 2 --------------------------------------------------------------
banner("Table 2. Descriptive statistics")
skw <- function(x) mean((x - mean(x))^3) / sd(x)^3
kur <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3
tab2 <- data.frame(Variable = c("Industrial production (QIM)","Wheat (PKR/10 kg)","Fuel (PKR/L)",
                                "Electricity (PKR/kWh)","Gas (PKR/MMBTU)"),
                   Mean = sapply(d[, c("QIM", PR)], mean), SD = sapply(d[, c("QIM", PR)], sd))
tab2$CV_pct <- 100 * tab2$SD / tab2$Mean
tab2$Skewness <- sapply(d[, c("QIM", PR)], skw); tab2$Excess_kurtosis <- sapply(d[, c("QIM", PR)], kur)
tab2 <- tab2[c(1, 2, 3, 4, 5), ]; rownames(tab2) <- NULL
print(r3(tab2)); save_tab(tab2, "Table02_descriptives.csv")
cat(sprintf("QIM minimum %.2f in %s, %.1f standard deviations below the mean.\n",
            min(y), format(d$date[which.min(y)], "%B %Y"), (mean(y) - min(y)) / sd(y)))

# ---- Figure 1 -------------------------------------------------------------
banner("Figure 1. Monthly series as recorded (saved as Figure1_series.pdf)")
L1 <- d %>% pivot_longer(-date) %>%
  mutate(name = factor(name, levels = c("QIM","WHEAT","FUEL","ENERGY","GAS"),
    labels = c("Quantum Index of Manufacturing (2015-16 = 100)","Wheat (PKR per 10 kg)",
               "Fuel: mean of petrol and diesel (PKR per litre)","Electricity (PKR per kWh)","Natural gas (PKR per MMBTU)")))
p1 <- ggplot(L1, aes(date, value)) + geom_line(linewidth = 0.45, colour = "#1f3b57") +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = as.Date("2024-01-01"), linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0)) + labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9) + theme(strip.background = element_rect(fill = "grey93"),
                                  panel.grid.minor = element_blank(), strip.text = element_text(hjust = 0, face = "bold"))
ggsave(file.path(OUT_DIR, "Figure1_series.pdf"), p1, width = 6.3, height = 7.6)
cat("Seasonal profile (monthly means of QIM, 2020 excluded):\n")
print(round(tapply(y[format(d$date, "%Y") != "2020"], mon[format(d$date, "%Y") != "2020"], mean), 1))

# =============================================================================
#  SECTION 3. PHASE I: LINEAR MODELS AND THE LONG-RUN QUESTION
# =============================================================================
# ---- Table 3 --------------------------------------------------------------
banner("Table 3. Unit root tests in levels and first differences")
tab3 <- do.call(rbind, lapply(c("QIM", PR), function(v) {
  x <- d[[v]]; dx <- diff(x)
  a1 <- ur.df(x, type = "trend", selectlags = "AIC"); a2 <- ur.df(dx, type = "drift", selectlags = "AIC")
  k1 <- ur.kpss(x, type = "tau", lags = "short");    k2 <- ur.kpss(dx, type = "mu", lags = "short")
  cls <- if (a1@teststat[1] < a1@cval[1, 2] && k1@teststat <= k1@cval[2]) "I(0)" else
         if (a2@teststat[1] < a2@cval[1, 2] && k2@teststat <= k2@cval[2]) "I(1)" else "ambiguous"
  data.frame(Variable = v, ADF_level = a1@teststat[1], KPSS_level = k1@teststat,
             ADF_diff = a2@teststat[1], KPSS_diff = k2@teststat, Classification = cls)
}))
print(r3(tab3)); save_tab(tab3, "Table03_unit_roots.csv")
cat("5% critical values: ADF levels -3.43, KPSS levels 0.146, ADF differences -2.88, KPSS differences 0.463\n")

# ---- ARDL selection ---------------------------------------------------------
banner("Section 3.3. ARDL lag selection (AIC, maximum order 6)")
ardl_sel   <- auto_ardl(QIM ~ WHEAT + FUEL + ENERGY + GAS, data = d, max_order = 6)
ardl_model <- ardl_sel$best_model
ord <- ardl_sel$best_order; names(ord) <- c("QIM", PR)
cat("Selected order (QIM, Wheat, Fuel, Electricity, Gas):", ord, "\n")
stopifnot(ord[1] == 1)       # the UECM code below assumes p = 1 (as selected)
q  <- ord[PR]

# UECM written out explicitly (matches the ARDL package's parameterisation)
Lg <- function(v, k) c(rep(NA, k), v[seq_len(length(v) - k)])
Df <- function(v) c(NA, diff(v))
lv <- paste0(PR, "_l1")
build <- function(y, X, extra = NULL) {
  Z <- data.frame(dy = Df(y), y1 = Lg(y, 1))
  for (k in PR) Z[[paste0(k, "_l1")]] <- if (q[k] == 0) X[, k] else Lg(X[, k], 1)
  for (k in PR) { if (q[k] == 0) next
    for (j in 0:(q[k] - 1)) Z[[paste0("d", k, "_", j)]] <- Lg(Df(X[, k]), j) }
  if (!is.null(extra)) Z <- cbind(Z, extra)
  Z
}
fit_uecm <- function(Z) lm(dy ~ ., data = Z)
stats3 <- function(fit) {
  b <- coef(fit); V <- vcov(fit)
  W <- function(nm) as.numeric(t(b[nm]) %*% solve(V[nm, nm]) %*% b[nm]) / length(nm)
  c(F_overall = W(c("y1", lv)), t_dep = unname(b["y1"] / sqrt(V["y1", "y1"])), F_exog = W(lv))
}
fast_stats <- function(dy, M) {
  qr_ <- qr(M); b <- qr.coef(qr_, dy); e <- dy - M %*% b
  s2 <- sum(e^2) / (nrow(M) - ncol(M)); XtXi <- chol2inv(qr.R(qr_)); piv <- qr_$pivot
  V <- matrix(NA, ncol(M), ncol(M)); V[piv, piv] <- XtXi * s2; dimnames(V) <- list(colnames(M), colnames(M))
  names(b) <- colnames(M)
  W <- function(nm) as.numeric(t(b[nm]) %*% solve(V[nm, nm]) %*% b[nm]) / length(nm)
  c(F_overall = W(c("y1", lv)), t_dep = unname(b["y1"] / sqrt(V["y1", "y1"])), F_exog = W(lv))
}
# McNown, Sam & Goh (2018) conditional bootstrap under each null
boot_null <- function(drop, stat, B = B_BOOT, extra = NULL) {
  Z <- build(y, X, extra); cc <- complete.cases(Z); idx <- which(cc)
  M <- cbind("(Intercept)" = 1, as.matrix(Z[cc, -1])); dy0 <- Z$dy[cc]
  keep <- setdiff(colnames(M), drop)
  br <- qr.coef(qr(M[, keep]), dy0); e <- dy0 - M[, keep] %*% br; e <- as.numeric(e - mean(e))
  has_y1 <- "y1" %in% keep
  mu  <- as.numeric(M[, setdiff(keep, "y1"), drop = FALSE] %*% br[setdiff(keep, "y1")])
  lam <- if (has_y1) br["y1"] else 0
  res <- numeric(B)
  for (b in 1:B) {
    ys <- y; es <- sample(e, length(idx), replace = TRUE)
    for (i in seq_along(idx)) { t <- idx[i]; ys[t] <- ys[t - 1] + mu[i] + lam * ys[t - 1] + es[i] }
    Ms <- M; Ms[, "y1"] <- ys[idx - 1]
    res[b] <- fast_stats(diff(ys)[idx - 1], Ms)[stat]
  }
  res
}
bootrep <- function(extra, label) {
  f <- fit_uecm(build(y, X, extra)); s <- stats3(f)
  bF <- boot_null(c("y1", lv), "F_overall", extra = extra)
  bt <- boot_null("y1", "t_dep", extra = extra)
  bE <- boot_null(lv, "F_exog", extra = extra)
  data.frame(Specification = label,
             F_ov = s[1], boot_p_F = mean(bF >= s[1]), boot_cv5_F = quantile(bF, .95),
             t_y1 = s[2], boot_p_t = mean(bt <= s[2]), boot_cv5_t = quantile(bt, .05),
             F_exog = s[3], boot_p_Fexog = mean(bE >= s[3]), boot_cv5_Fexog = quantile(bE, .95),
             lambda0 = coef(f)["y1"], row.names = NULL)
}
SD  <- sapply(2:12, function(k) as.numeric(mon == k)); colnames(SD) <- paste0("m", 2:12)
COV <- cbind(covid_apr20 = as.numeric(d$date == as.Date("2020-04-01")),
             covid_may20 = as.numeric(d$date == as.Date("2020-05-01")))

# ---- Table 4 --------------------------------------------------------------
banner("Table 4. Bounds tests with bootstrap critical values (McNown et al., 2018)")
set.seed(20261022)
tab4 <- rbind(bootrep(NULL, "Baseline ARDL"), bootrep(SD, "+ monthly dummies"),
              bootrep(cbind(SD, COV), "+ monthly and COVID dummies"))
tr_idx <- d$date <= as.Date("2023-12-01")
ardl_tr <- ardl(QIM ~ WHEAT + FUEL + ENERGY + GAS, data = d[tr_idx, ], order = ord)
bf_tr <- bounds_f_test(ardl_tr, case = 3); bt_tr <- bounds_t_test(ardl_tr, case = 3)
print(r3(tab4))
cat(sprintf("Training sample only (to Dec 2023): F = %.3f (asymptotic p = %.3f), t = %.3f (asymptotic p = %.3f)\n",
            bf_tr$statistic, bf_tr$p.value, bt_tr$statistic, bt_tr$p.value))
cat("Check: ARDL package bounds F for the baseline =", round(bounds_f_test(ardl_model, case = 3)$statistic, 4), "\n")
save_tab(tab4, "Table04_bounds_bootstrap.csv")

# ---- Table 5 --------------------------------------------------------------
banner("Table 5. Long-run coefficients under alternative controls")
lr_of <- function(f) { b <- coef(f); Vv <- vcov(f); g <- -b[lv] / b["y1"]
  se <- sapply(lv, function(k) { gr <- setNames(rep(0, length(b)), names(b))
    gr[k] <- -1 / b["y1"]; gr["y1"] <- b[k] / b["y1"]^2; sqrt(t(gr) %*% Vv %*% gr) })
  data.frame(Variable = PR, Coef = g, p = 2 * pt(-abs(g / se), df.residual(f))) }
f0 <- fit_uecm(build(y, X))
tab5 <- cbind(setNames(lr_of(f0), c("Variable","Baseline","p_baseline")),
              setNames(lr_of(fit_uecm(build(y, X, SD)))[, -1], c("Monthly","p_monthly")),
              setNames(lr_of(fit_uecm(build(y, X, cbind(SD, COV))))[, -1], c("Monthly_COVID","p_monthly_COVID")))
rownames(tab5) <- NULL; print(r3(tab5)); save_tab(tab5, "Table05_longrun.csv")
th0 <- -coef(f0)[1] / coef(f0)["y1"]
cat(sprintf("Baseline long-run intercept: %.2f.  Half-life of a deviation: %.2f months.\n",
            th0, log(0.5) / log(1 + coef(f0)["y1"])))
lr_tr <- multipliers(ardl_tr, type = "lr")
cat("Training-sample long-run coefficients (text, Section 3.3):\n"); print(r3(lr_tr[, 1:2]))

# ---- Weak exogeneity (Section 3.3) -------------------------------------------
banner("Section 3.3. Weak exogeneity of the prices")
lrb <- -coef(f0)[lv] / coef(f0)["y1"]
ect <- as.numeric(y - th0 - X %*% lrb)
we <- do.call(rbind, lapply(PR, function(k) {
  Zm <- data.frame(dx = Df(X[, k]), ect1 = Lg(ect, 1), dy1 = Lg(Df(y), 1), dy2 = Lg(Df(y), 2))
  for (kk in PR) { Zm[[paste0("d", kk, "1")]] <- Lg(Df(X[, kk]), 1); Zm[[paste0("d", kk, "2")]] <- Lg(Df(X[, kk]), 2) }
  f <- lm(dx ~ ., data = Zm); ct <- coeftest(f, vcov = NeweyWest(f, lag = 4, prewhite = FALSE))
  data.frame(Variable = k, alpha_ect = ct["ect1", 1], se = ct["ect1", 2], p = ct["ect1", 4]) }))
print(r3(we)); save_tab(we, "Section3_weak_exogeneity.csv")

# ---- Parameter stability (Section 3.3) ----------------------------------------
banner("Section 3.3. Parameter stability")
Z0 <- build(y, X); Zc <- Z0[complete.cases(Z0), ]
rc <- efp(dy ~ ., data = Zc, type = "Rec-CUSUM"); sc_rc <- sctest(rc)
fs <- suppressWarnings(Fstats(dy ~ ., data = Zc, from = 0.15)); sc_f <- sctest(fs, type = "supF")
cat(sprintf("Recursive CUSUM: S = %.2f, p = %.2f\nAndrews sup-F: %.1f, p = %.4f, break dated %s\n",
            sc_rc$statistic, sc_rc$p.value, sc_f$statistic, sc_f$p.value,
            format(d$date[complete.cases(Z0)][fs$breakpoint], "%B %Y")))

# ---- Table 6 --------------------------------------------------------------
banner("Table 6. Residual diagnostics for the ARDL model")
r0 <- residuals(ardl_model); fC <- fit_uecm(build(y, X, cbind(SD, COV))); rC <- resid(fC)
tab6a <- data.frame(Test = c("Breusch-Godfrey order 1","Breusch-Godfrey order 4","Breusch-Godfrey order 12",
                             "ARCH-LM 12 lags","Jarque-Bera","RESET (powers 2,3)"),
  Statistic = c(sapply(c(1, 4, 12), function(o) bgtest(ardl_model, order = o)$statistic),
                ArchTest(r0, lags = 12)$statistic, jarque.bera.test(r0)$statistic,
                resettest(ardl_model, power = 2:3, type = "fitted")$statistic),
  p = c(sapply(c(1, 4, 12), function(o) bgtest(ardl_model, order = o)$p.value),
        ArchTest(r0, lags = 12)$p.value, jarque.bera.test(r0)$p.value,
        resettest(ardl_model, power = 2:3, type = "fitted")$p.value),
  p_with_month_COVID = c(sapply(c(1, 4, 12), function(o) bgtest(fC, order = o)$p.value),
        ArchTest(rC, lags = 12)$p.value, jarque.bera.test(rC)$p.value,
        resettest(fC, power = 2:3, type = "fitted")$p.value))
print(r3(tab6a))
bdsrow <- function(x, lab) { b <- bds.test(x, m = 5, eps = sd(x) * c(0.5, 1, 1.5, 2))
  do.call(rbind, lapply(1:4, function(i) data.frame(Series = lab, eps = c("0.5sd","1.0sd","1.5sd","2.0sd")[i],
     m2 = b$p.value[1, i], m3 = b$p.value[2, i], m4 = b$p.value[3, i], m5 = b$p.value[4, i]))) }
gfit <- suppressWarnings(garch(rC, order = c(1, 1), trace = FALSE)); hs <- fitted(gfit)[, 1]
tab6b <- rbind(bdsrow(r0, "Baseline residuals"), bdsrow(rC, "With month + COVID dummies"),
               bdsrow((rC / hs)[!is.na(hs)], "GARCH(1,1)-standardised"))
cat("\nBDS test p-values by embedding dimension m and distance eps:\n"); print(r3(tab6b))
save_tab(tab6a, "Table06a_diagnostics.csv"); save_tab(tab6b, "Table06b_BDS.csv")

# Appendix Table B.12 is computed here and printed at the end of the script
ctab <- coeftest(f0); hac <- coeftest(f0, vcov = NeweyWest(f0, lag = 4, prewhite = FALSE))
tabB12 <- data.frame(Term = rownames(ctab), Coef = ctab[, 1], OLS_se = ctab[, 2], OLS_p = ctab[, 4],
                     HAC_se = hac[, 2], HAC_p = hac[, 4], row.names = NULL)
B12_note <- sprintf("N = %d; R2 = %.3f; adjusted R2 = %.3f", nobs(f0), summary(f0)$r.squared, summary(f0)$adj.r.squared)

# =============================================================================
#  SECTION 4. PHASE II: FORECASTING DESIGN
# =============================================================================
banner("Table 7. Forecasting methods")
lab <- c(SNaive = "Seasonal naive", ARIMA = "Seasonal ARIMA", ETS = "ETS", STL_ARIMA = "STL-ARIMA",
  ARIMAX_feasible = "ARIMAX (no-change prices)", ARDL_feasible = "ARDL (no-change prices)",
  RF = "Random forest", GBM = "GBM", H1_ARIMA_GBM = "H1: ARIMA-GBM", H2_ARIMA_RF = "H2: ARIMA-RF",
  H3_STL_GBM = "H3: STL-GBM", H4_ETS_GBM = "H4: ETS-GBM", GBM_noCom = "GBM, no prices",
  H3_STL_GBM_noCom = "H3: STL-GBM, no prices", STL_GBM_direct = "STL + GBM (no ARIMA)",
  ARIMAX_conditional = "ARIMAX (actual prices, reference only)")
print(data.frame(Code = names(lab), Method = unname(lab)), right = FALSE)
feasible <- setdiff(names(lab), "ARIMAX_conditional")

# =============================================================================
#  SECTION 5. PHASE II: FORECASTING RESULTS
# =============================================================================
H <- c(1, 3, 6, 12); W0 <- 66; origins <- W0:(N - 1)

feat <- function(t, h, LAGS, sa = NULL, com = TRUE) {      # information set at origin t - h
  yy <- if (is.null(sa)) y else sa
  v <- c(sapply(0:(LAGS - 1), function(j) yy[t - h - j]))
  if (com) for (k in PR) v <- c(v, sapply(0:(LAGS - 1), function(j) X[t - h - j, k]))
  v <- c(v, mon[t], t)
  names(v) <- c(paste0("QIM_l", 0:(LAGS - 1)),
                if (com) as.vector(sapply(PR, function(k) paste0(k, "_l", 0:(LAGS - 1)))), "month", "trend")
  v
}
fmat <- function(ts_, h, LAGS, ...) do.call(rbind, lapply(ts_, feat, h = h, LAGS = LAGS, ...))
gbm_fit_pred <- function(Xtr, ytr, Xte, depth = 3, ntmax = 500, seed = 42) {
  set.seed(seed); n <- nrow(Xtr); nv <- 12
  df <- data.frame(yv = ytr, Xtr); dte <- data.frame(Xte); colnames(dte) <- colnames(df)[-1]
  f1 <- gbm(yv ~ ., data = df[1:(n - nv), ], distribution = "gaussian", n.trees = ntmax, interaction.depth = depth,
            shrinkage = 0.05, bag.fraction = 0.8, n.minobsinnode = 5, verbose = FALSE)
  pv <- predict(f1, df[(n - nv + 1):n, ], n.trees = 1:ntmax)
  best <- max(10, which.min(colMeans((pv - ytr[(n - nv + 1):n])^2)))     # time-ordered validation
  set.seed(seed)
  f2 <- gbm(yv ~ ., data = df, distribution = "gaussian", n.trees = best, interaction.depth = depth,
            shrinkage = 0.05, bag.fraction = 0.8, n.minobsinnode = 5, verbose = FALSE)
  as.numeric(predict(f2, dte, n.trees = best))
}
rf_fit_pred <- function(Xtr, ytr, Xte, seed = 42) { set.seed(seed)
  as.numeric(predict(randomForest(Xtr, ytr, ntree = 500, mtry = floor(sqrt(ncol(Xtr))), nodesize = 3), Xte)) }
ardl_fc <- function(tr, o, hmax) {             # training-only ARDL, recursive, no-change price paths
  a <- auto_ardl(QIM ~ WHEAT + FUEL + ENERGY + GAS, data = d[tr, ], max_order = 6)$best_model
  b <- coef(a); yy <- y; XX <- X
  for (k in PR) XX[(o + 1):N, k] <- X[o, k]
  fc <- numeric(hmax)
  for (s in 1:hmax) { t <- o + s; v <- b["(Intercept)"]
    for (nm in names(b)[-1]) {
      lag <- if (grepl("^L\\(", nm)) as.numeric(sub(".*, *([0-9]+)\\)$", "\\1", nm)) else 0
      var <- sub("^L\\(([A-Z]+),.*$", "\\1", nm)
      v <- v + b[nm] * (if (var == "QIM") yy[t - lag] else XX[t - lag, var]) }
    yy[t] <- v; fc[s] <- v }
  list(fc = fc, order = paste(a$order, collapse = ","))
}
one_origin <- function(o, scheme, LAGS) {
  tr <- if (scheme == "expanding") 1:o else (o - W0 + 1):o
  st <- c(as.numeric(format(d$date[tr[1]], "%Y")), mon[tr[1]])
  yt <- ts(y[tr], start = st, frequency = 12); hmax <- min(12, N - o)
  aa <- function(z, ...) auto.arima(z, ic = "aic", stepwise = TRUE, approximation = FALSE, ...)
  snaive <- sapply(1:hmax, function(h) y[o + h - 12])
  fA <- aa(yt, seasonal = TRUE); pA <- as.numeric(forecast(fA, h = hmax)$mean)
  fE <- ets(yt);                 pE <- as.numeric(forecast(fE, h = hmax)$mean)
  sl <- stl(yt, s.window = "periodic", robust = TRUE); seas <- as.numeric(sl$time.series[, "seasonal"])
  sa <- yt - seas; fS <- aa(sa, seasonal = FALSE)
  seas_f <- sapply(1:hmax, function(h) seas[length(seas) - 12 + ((h - 1) %% 12) + 1])
  pS <- as.numeric(forecast(fS, h = hmax)$mean) + seas_f
  fX <- aa(yt, xreg = X[tr, ], seasonal = TRUE)
  pXc <- as.numeric(forecast(fX, xreg = X[o + 1:hmax, , drop = FALSE])$mean)
  pXf <- as.numeric(forecast(fX, xreg = X[rep(o, hmax), , drop = FALSE])$mean)
  ad <- ardl_fc(tr, o, hmax); pD <- ad$fc
  eA <- as.numeric(y[tr] - fitted(fA)); eE <- as.numeric(y[tr] - fitted(fE)); eS <- as.numeric(sa - fitted(fS))
  sa_full <- y; sa_full[tr] <- as.numeric(sa)
  res <- list()
  for (h in H[H <= hmax]) {
    t_tr <- if (scheme == "expanding") tr[tr - h - (LAGS - 1) >= 1] else tr[tr - h - (LAGS - 1) >= tr[1]]
    tt <- o + h; idx <- match(t_tr, tr)
    F1 <- fmat(t_tr, h, LAGS); f1 <- rbind(feat(tt, h, LAGS))
    F0 <- fmat(t_tr, h, LAGS, com = FALSE); f0_ <- rbind(feat(tt, h, LAGS, com = FALSE))
    FS <- fmat(t_tr, h, LAGS, sa = sa_full); fS_ <- rbind(feat(tt, h, LAGS, sa = sa_full))
    r <- c(SNaive = snaive[h], ARIMA = pA[h], ETS = pE[h], STL_ARIMA = pS[h],
           ARIMAX_feasible = pXf[h], ARIMAX_conditional = pXc[h], ARDL_feasible = pD[h],
           RF = rf_fit_pred(F1, y[t_tr], f1), GBM = gbm_fit_pred(F1, y[t_tr], f1),
           GBM_noCom = gbm_fit_pred(F0, y[t_tr], f0_),
           H1_ARIMA_GBM = pA[h] + gbm_fit_pred(F1, eA[idx], f1, depth = 2, ntmax = 300),
           H2_ARIMA_RF  = pA[h] + rf_fit_pred(F1, eA[idx], f1),
           H3_STL_GBM   = pS[h] + gbm_fit_pred(F1, eS[idx], f1, depth = 2, ntmax = 300),
           H3_STL_GBM_noCom = pS[h] + gbm_fit_pred(F0, eS[idx], f0_, depth = 2, ntmax = 300),
           H4_ETS_GBM   = pE[h] + gbm_fit_pred(F1, eE[idx], f1, depth = 2, ntmax = 300),
           STL_GBM_direct = seas_f[h] + gbm_fit_pred(FS, sa_full[t_tr], fS_))
    res[[length(res) + 1]] <- data.frame(origin = o, h = h, target = tt, date = d$date[tt], actual = y[tt],
                                         model = names(r), fc = as.numeric(r))
  }
  do.call(rbind, res)
}
run_scheme <- function(scheme, LAGS) {
  tag <- if (LAGS == 6) scheme else paste0(scheme, "_L", LAGS)
  cf  <- file.path(OUT_DIR, "cache", paste0("rolling_", tag, ".rds"))
  if (file.exists(cf)) { cat("  using cached forecasts:", tag, "\n"); return(readRDS(cf)) }
  cat(sprintf("  computing %s (%d origins, %d core(s)) ...\n", tag, length(origins), N_CORES)); t0 <- Sys.time()
  f <- function(o) one_origin(o, scheme, LAGS)
  if (N_CORES > 1 && .Platform$OS.type == "windows") {
    cl <- makeCluster(N_CORES); on.exit(stopCluster(cl))
    clusterEvalQ(cl, suppressPackageStartupMessages({ library(forecast); library(gbm); library(randomForest); library(ARDL) }))
    clusterExport(cl, c("one_origin","feat","fmat","gbm_fit_pred","rf_fit_pred","ardl_fc","d","y","X","mon","N",
                        "PR","H","W0"), envir = environment())
    R <- parLapply(cl, origins, f)
  } else if (N_CORES > 1) {
    R <- mclapply(origins, f, mc.cores = N_CORES, mc.preschedule = FALSE)
  } else R <- lapply(origins, f)
  R <- do.call(rbind, R); R$tag <- tag; R$err <- R$actual - R$fc
  saveRDS(R, cf); cat("  done in", format(round(Sys.time() - t0, 1)), "\n"); R
}
banner("Section 5. Rolling-origin forecasts (46 origins, Jan 2022 - Oct 2025)")
RALL <- rbind(run_scheme("expanding", 6), run_scheme("rolling", 6),
              run_scheme("expanding", 3), run_scheme("expanding", 12))

# ---- accuracy, Diebold-Mariano, MCS ----------------------------------------------
WIN <- list(A = as.Date(c("2024-01-01","2025-11-01")), B = as.Date(c("2023-01-01","2025-11-01")))
mase_scale <- mean(abs(diff(y[1:W0], lag = 12)))
acc <- RALL %>% group_by(tag, h, model) %>% group_modify(function(g, k) {
  do.call(rbind, lapply(names(WIN), function(w) { s <- g[g$date >= WIN[[w]][1] & g$date <= WIN[[w]][2], ]
    data.frame(window = w, n = nrow(s), RMSE = sqrt(mean(s$err^2)), MAE = mean(abs(s$err)),
               MAPE = 100 * mean(abs(s$err / s$actual)), MASE = mean(abs(s$err)) / mase_scale) })) }) %>% ungroup() %>%
  group_by(tag, h, window) %>% mutate(relRMSE = RMSE / RMSE[model == "SNaive"]) %>% ungroup()
dm_hln <- function(e1, e2, h) { r <- tryCatch(dm.test(e1, e2, h = h, power = 2, alternative = "two.sided",
                                                      varestimator = "bartlett"), error = function(e) NULL)
  if (is.null(r)) c(NA, NA) else c(unname(r$statistic), r$p.value) }
mcs_col <- function(sh) { nm <- colnames(sh)
  if ("MCS_M" %in% nm) "MCS_M" else if ("MCS p-Value" %in% nm) "MCS p-Value" else tail(nm, 1) }
DM <- list(); MC <- list()
for (tg in unique(RALL$tag)) for (hh in H) for (w in names(WIN)) {
  s <- RALL[RALL$tag == tg & RALL$h == hh & RALL$date >= WIN[[w]][1] & RALL$date <= WIN[[w]][2] & RALL$model %in% feasible, ]
  E <- s %>% select(date, model, err) %>% pivot_wider(names_from = model, values_from = err) %>% arrange(date)
  E <- as.data.frame(E[complete.cases(E), ])
  for (cmp in setdiff(feasible, "H3_STL_GBM")) { r <- dm_hln(E$H3_STL_GBM, E[[cmp]], hh)
    DM[[length(DM) + 1]] <- data.frame(tag = tg, h = hh, window = w, set = "H3", comparator = cmp, DM = r[1], p = r[2]) }
  for (cmp in setdiff(feasible, "SNaive")) { r <- dm_hln(E[[cmp]], E$SNaive, hh)
    DM[[length(DM) + 1]] <- data.frame(tag = tg, h = hh, window = w, set = "vsSNaive", comparator = cmp, DM = r[1], p = r[2]) }
  mcs <- MCSprocedure(as.matrix(E[, feasible]^2), alpha = 0.10, B = B_MCS, statistic = "Tmax",
                      k = max(1, hh), verbose = FALSE, seed = 123)
  sh <- mcs@show; pc <- sh[, mcs_col(sh)]
  MC[[length(MC) + 1]] <- data.frame(tag = tg, h = hh, window = w, model = rownames(sh), MCS_p = pc,
                                     in90 = pc >= 0.10, in75 = pc >= 0.25, row.names = NULL)
}
DM <- do.call(rbind, DM) %>% group_by(tag, h, window, set) %>% mutate(p_holm = p.adjust(p, "holm")) %>% ungroup()
MC <- do.call(rbind, MC)
save_tab(acc, "Section5_accuracy_all.csv"); save_tab(DM, "Section5_DM_all.csv"); save_tab(MC, "Section5_MCS_all.csv")

rmse_table <- function(tg) {
  x <- acc %>% filter(tag == tg) %>% left_join(MC %>% select(tag, h, window, model, in90, in75),
                                               by = c("tag","h","window","model")) %>%
    group_by(window, h) %>% mutate(best = model %in% feasible & RMSE == min(RMSE[model %in% feasible])) %>% ungroup()
  x$cell <- paste0(sprintf("%.2f", x$RMSE), ifelse(x$best, " (best)", ""),
                   ifelse(is.na(x$in90), "", ifelse(x$in75, " **", ifelse(x$in90, " *", ""))))
  w <- x %>% select(model, window, h, cell) %>% pivot_wider(names_from = c(window, h), values_from = cell, names_prefix = "")
  w <- w[match(names(lab), w$model), c("model","A_1","A_3","A_6","A_12","B_1","B_3","B_6","B_12")]
  w$model <- lab[w$model]; names(w)[1] <- "Method"; as.data.frame(w)
}

# ---- Table 8 --------------------------------------------------------------
banner("Table 8. Out-of-sample RMSE, expanding window  (A = Jan 2024-Nov 2025, B = Jan 2023-Nov 2025)")
tab8 <- rmse_table("expanding"); print(tab8, right = FALSE); save_tab(tab8, "Table08_RMSE_expanding.csv")
cat("** = in 75% and 90% model confidence sets; * = in 90% set only; no mark = excluded.\n")

banner("Section 5.1. In-text summaries")
rk <- acc %>% filter(model %in% feasible) %>% group_by(tag, window, h) %>%
  summarise(best = lab[model[which.min(RMSE)]], best_RMSE = min(RMSE), H3_RMSE = RMSE[model == "H3_STL_GBM"],
            H3_rank = rank(RMSE)[model == "H3_STL_GBM"], SNaive_rank = rank(RMSE)[model == "SNaive"], .groups = "drop") %>%
  left_join(MC %>% group_by(tag, window, h) %>% summarise(in_MCS90 = sum(in90), in_MCS75 = sum(in75), .groups = "drop"),
            by = c("tag","window","h"))
print(as.data.frame(r3(rk %>% filter(tag == "expanding"))))
cat(sprintf("H3 rank range (expanding): %d to %d of %d feasible methods\n",
            min(rk$H3_rank[rk$tag == "expanding"]), max(rk$H3_rank[rk$tag == "expanding"]), length(feasible)))
cat(sprintf("90%% MCS size across all schemes: %d to %d of %d\n", min(rk$in_MCS90), max(rk$in_MCS90), length(feasible)))
v <- DM %>% filter(tag == "expanding", set == "vsSNaive")
best_vs_sn <- rk %>% filter(tag == "expanding") %>% rowwise() %>%
  mutate(code = names(lab)[match(best, lab)]) %>% ungroup() %>% filter(code != "SNaive") %>%
  left_join(v, by = c("tag","h","window","code" = "comparator"))
cat("Best method vs seasonal naive (HLN-DM, expanding):\n"); print(as.data.frame(r3(best_vs_sn[, c("window","h","best","DM","p","p_holm")])))
cat("Comparisons that survive Holm adjustment (all schemes):\n")
print(as.data.frame(r3(DM %>% filter(p_holm < 0.05) %>% select(tag, h, window, set, comparator, DM, p_holm))))
o_dec <- which(d$date == as.Date("2023-12-01"))
a_dec <- auto_ardl(QIM ~ WHEAT + FUEL + ENERGY + GAS, data = d[1:o_dec, ], max_order = 6)$best_model
lr_dec <- multipliers(a_dec)
cat(sprintf("Equilibrium implied by the training-only ARDL at December 2023 prices: %.1f index points\n",
            sum(c(1, X[o_dec, ]) * lr_dec$Estimate)))

# ---- Figure 2 -------------------------------------------------------------
banner("Figure 2. Rolling-origin forecasts at h = 1 and h = 12 (saved as Figure2_forecasts.pdf)")
pi_emp <- function(tg, mdl, hh) {
  s <- RALL %>% filter(tag == tg, model == mdl, h == hh) %>% arrange(origin)
  q <- t(sapply(seq_len(nrow(s)), function(i) { past <- s$err[s$target <= s$origin[i]]
    if (length(past) >= 12) s$fc[i] + quantile(past, c(.10, .90, .025, .975)) else rep(NA, 4) }))
  cbind(s, lo80 = q[, 1], hi80 = q[, 2], lo95 = q[, 3], hi95 = q[, 4])
}
keep <- c(SNaive = "Seasonal naive", ARIMA = "Seasonal ARIMA", GBM = "GBM", H1_ARIMA_GBM = "H1: ARIMA-GBM", H3_STL_GBM = "H3: STL-GBM")
fd <- RALL %>% filter(tag == "expanding", model %in% names(keep), h %in% c(1, 12), date >= WIN$A[1]) %>%
  mutate(Model = factor(keep[model], levels = keep), Hz = factor(paste0("h = ", h), levels = c("h = 1","h = 12")))
act <- RALL %>% filter(tag == "expanding", model == "SNaive", h == 1, date >= WIN$A[1]) %>% select(date, actual)
act <- rbind(cbind(act, Hz = "h = 1"), cbind(act, Hz = "h = 12")); act$Hz <- factor(act$Hz, levels = c("h = 1","h = 12"))
band <- rbind(pi_emp("expanding", "H1_ARIMA_GBM", 1), pi_emp("expanding", "H1_ARIMA_GBM", 12)) %>%
  filter(date >= WIN$A[1], !is.na(lo80)) %>% mutate(Hz = factor(paste0("h = ", h), levels = c("h = 1","h = 12")))
p2 <- ggplot() + geom_ribbon(data = band, aes(date, ymin = lo80, ymax = hi80), fill = "#d62728", alpha = 0.12) +
  geom_line(data = fd, aes(date, fc, colour = Model, linetype = Model), linewidth = 0.55) +
  geom_line(data = act, aes(date, actual), colour = "black", linewidth = 0.9) + facet_wrap(~Hz, ncol = 1) +
  scale_colour_manual(values = c("#7f7f7f","#1f77b4","#2ca02c","#d62728","#9467bd"), name = NULL) +
  scale_linetype_manual(values = c("dotted","dashed","dotdash","solid","longdash"), name = NULL) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") + labs(x = NULL, y = "QIM (index)") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                                  strip.background = element_rect(fill = "grey93"))
ggsave(file.path(OUT_DIR, "Figure2_forecasts.pdf"), p2, width = 6.5, height = 6.2)

# ---- Figure 3 -------------------------------------------------------------
banner("Figure 3. RMSE relative to the seasonal naive (saved as Figure3_relative_RMSE.pdf)")
f3m <- c(ARIMA = "Seasonal ARIMA", ETS = "ETS", STL_ARIMA = "STL-ARIMA", ARIMAX_feasible = "ARIMAX (no-change prices)",
         GBM = "GBM", GBM_noCom = "GBM, no prices", H1_ARIMA_GBM = "H1: ARIMA-GBM", H3_STL_GBM = "H3: STL-GBM",
         H3_STL_GBM_noCom = "H3: STL-GBM, no prices")
a3 <- acc %>% filter(tag == "expanding", model %in% names(f3m)) %>%
  mutate(Model = factor(f3m[model], levels = f3m),
         Window = ifelse(window == "A", "Window A: Jan 2024 - Nov 2025", "Window B: Jan 2023 - Nov 2025"))
p3 <- ggplot(a3, aes(factor(h), relRMSE, group = Model, colour = Model)) + geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.3) + facet_wrap(~Window) + scale_y_log10() +
  labs(x = "Horizon (months)", y = "RMSE relative to seasonal naive (log scale)", colour = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                                  strip.background = element_rect(fill = "grey93")) + guides(colour = guide_legend(ncol = 3))
ggsave(file.path(OUT_DIR, "Figure3_relative_RMSE.pdf"), p3, width = 6.5, height = 4.4)
print(as.data.frame(r3(a3 %>% select(window, h, model, relRMSE) %>% pivot_wider(names_from = h, values_from = relRMSE))))

# ---- Table 9 --------------------------------------------------------------
banner("Table 9. Ablations: RMSE ratio with / without the ingredient (expanding window)")
pairs <- list(c("STL_ARIMA","ARIMA","STL decomposition, no ML (STL-ARIMA / ARIMA)"),
              c("H3_STL_GBM","H1_ARIMA_GBM","STL decomposition, GBM corrector (H3 / H1)"),
              c("H3_STL_GBM","STL_ARIMA","GBM corrector on STL-ARIMA (H3 / STL-ARIMA)"),
              c("H1_ARIMA_GBM","ARIMA","GBM corrector on ARIMA (H1 / ARIMA)"),
              c("H3_STL_GBM","H3_STL_GBM_noCom","Prices in H3 corrector (H3 / H3 no prices)"),
              c("GBM","GBM_noCom","Prices in GBM (GBM / GBM no prices)"),
              c("ARIMAX_feasible","ARIMA","Prices in ARIMA (ARIMAX / ARIMA)"))
tab9 <- do.call(rbind, lapply(pairs, function(p) acc %>% filter(tag == "expanding") %>% group_by(window, h) %>%
  summarise(Ingredient = p[3], ratio = RMSE[model == p[1]] / RMSE[model == p[2]], .groups = "drop"))) %>%
  pivot_wider(names_from = c(window, h), values_from = ratio)
tab9 <- as.data.frame(tab9[, c("Ingredient","A_1","A_3","A_6","A_12","B_1","B_3","B_6","B_12")])
print(r3(tab9), right = FALSE); save_tab(tab9, "Table09_ablations.csv")
cat("DM tests for the H3-based pairs at h = 1 (note to Table 9):\n")
print(as.data.frame(r3(DM %>% filter(tag == "expanding", set == "H3", h == 1,
                                     comparator %in% c("STL_ARIMA","H1_ARIMA_GBM","H3_STL_GBM_noCom")) %>%
                         select(window, comparator, DM, p, p_holm))))

# ---- Table 10 -------------------------------------------------------------
banner("Table 10. Robustness to the estimation scheme and the number of lags")
tab10 <- rk %>% filter(tag != "expanding") %>%
  mutate(tag = factor(tag, levels = c("rolling","expanding_L3","expanding_L12"),
                      labels = c("Rolling (66 months), 6 lags","Expanding, 3 lags","Expanding, 12 lags"))) %>%
  arrange(tag, window, h) %>% select(Scheme = tag, Window = window, h, Most_accurate = best, RMSE = best_RMSE,
                                     H3_RMSE, H3_rank, SNaive_rank, In_MCS90 = in_MCS90)
print(as.data.frame(r3(tab10)), right = FALSE); save_tab(tab10, "Table10_robustness.csv")

# ---- Table 11 -------------------------------------------------------------
banner("Table 11. From a single train/test split to the leakage-free design (RMSE, Jan 2024 - Nov 2025)")
# Row (a): the single-split design, reproduced as in the preliminary version
mrg <- d; train_n <- floor(0.8 * N); y_train <- y[1:train_n]; test_dates <- d$date[(train_n + 1):N]
cf_ <- mrg; for (v in c("QIM","FUEL","ENERGY","WHEAT","GAS")) for (k in 1:6) cf_[[paste0(v, "_lag", k)]] <- dplyr::lag(mrg[[v]], k)
cf_$month <- mon; cf_$trend <- seq_len(N); cf_ <- cf_[complete.cases(cf_), ]
fcols <- grep("_lag|month|trend", names(cf_), value = TRUE)
Xtr_ml <- as.matrix(cf_[cf_$date < min(test_dates), fcols]); ytr_ml <- cf_$QIM[cf_$date < min(test_dates)]
Xte_ml <- as.matrix(cf_[cf_$date %in% test_dates, fcols]);   yte_ml <- cf_$QIM[cf_$date %in% test_dates]
rm_ <- function(p) sqrt(mean((yte_ml - p)^2))
gbm_corr_orig <- function(resid_train, X_tr, X_te) {       # corrector as in the preliminary version
  r <- resid_train[!is.na(resid_train)]; yr <- tail(r, nrow(X_tr)); df <- as.data.frame(cbind(yr = yr, X_tr))
  set.seed(42); fit <- gbm(yr ~ ., data = df, distribution = "gaussian", n.trees = 300, interaction.depth = 2,
                           shrinkage = 0.05, bag.fraction = 0.80, cv.folds = 5, verbose = FALSE)
  predict(fit, as.data.frame(X_te), n.trees = gbm.perf(fit, method = "cv", plot.it = FALSE)) }
ar_ns <- auto.arima(y_train, seasonal = TRUE, stepwise = FALSE, approximation = FALSE, ic = "aic")   # no frequency -> non-seasonal
cat("Benchmark actually selected in the single-split design:", paste0("ARIMA(", paste(arimaorder(ar_ns), collapse = ","), ")"), "\n")
p_ns <- as.numeric(forecast(ar_ns, h = length(test_dates))$mean)
set.seed(42)
g_or <- gbm(y ~ ., data = as.data.frame(cbind(y = ytr_ml, Xtr_ml)), distribution = "gaussian", n.trees = 500,
            interaction.depth = 3, shrinkage = 0.05, bag.fraction = 0.80, cv.folds = 5, verbose = FALSE)
p_gbm <- predict(g_or, as.data.frame(Xte_ml), n.trees = gbm.perf(g_or, method = "cv", plot.it = FALSE))
yfull <- ts(y, start = c(as.numeric(format(d$date[1], "%Y")), mon[1]), frequency = 12)
stl_full <- stl(yfull, s.window = "periodic", robust = TRUE); seas_full <- as.numeric(stl_full$time.series[, "seasonal"])
sa_tr <- window(yfull - seas_full, end = time(yfull)[train_n])
st_ar <- auto.arima(sa_tr, seasonal = FALSE, stepwise = FALSE, approximation = FALSE, ic = "aic")
h3_a  <- as.numeric(forecast(st_ar, h = length(test_dates))$mean) + seas_full[(train_n + 1):N] +
         gbm_corr_orig(as.numeric(st_ar$residuals), Xtr_ml, Xte_ml)
# Row (b): seasonal ARIMA on the same fixed-origin task
yts_tr <- ts(y_train, start = c(as.numeric(format(d$date[1], "%Y")), mon[1]), frequency = 12)
sar <- auto.arima(yts_tr, ic = "aic", stepwise = FALSE, approximation = FALSE)
p_sar <- as.numeric(forecast(sar, h = length(test_dates))$mean)
# Row (c): STL estimated on the training sample only
stl_tr <- stl(yts_tr, s.window = "periodic", robust = TRUE); seas_tr <- as.numeric(stl_tr$time.series[, "seasonal"])
fa <- auto.arima(yts_tr - seas_tr, seasonal = FALSE, stepwise = FALSE, approximation = FALSE, ic = "aic")
base_c <- as.numeric(forecast(fa, h = length(test_dates))$mean) + rep(tail(seas_tr, 12), length.out = length(test_dates))
h3_c <- base_c + gbm_corr_orig(as.numeric(residuals(fa)), Xtr_ml, Xte_ml)
rowd <- acc %>% filter(tag == "expanding", window == "A", h == 1)
tab11 <- data.frame(Design = c("(a) Single split; non-seasonal ARIMA; full-sample STL; mixed horizons",
                               "(b) As (a) but with a seasonal ARIMA benchmark",
                               "(c) As (b) but STL estimated on the training sample only",
                               "(d) Rolling origin, all models re-estimated, common horizon h = 1"),
                    H3_STL_GBM = c(rm_(h3_a), rm_(h3_a), rm_(h3_c), rowd$RMSE[rowd$model == "H3_STL_GBM"]),
                    Benchmark_ARIMA = c(rm_(p_ns), rm_(p_sar), rm_(p_sar), rowd$RMSE[rowd$model == "ARIMA"]),
                    GBM = c(rm_(p_gbm), rm_(p_gbm), rm_(p_gbm), rowd$RMSE[rowd$model == "GBM"]))
print(r3(tab11), right = FALSE); save_tab(tab11, "Table11_bridge.csv")
cat(sprintf("Row (c): STL-ARIMA component alone has RMSE %.2f.  Row (a): H3 is %.0f%% below ARIMA and %.0f%% below GBM.\n",
            rm_(base_c), 100 * (1 - rm_(h3_a) / rm_(p_ns)), 100 * (1 - rm_(h3_a) / rm_(p_gbm))))

# ---- Section 5.5 ------------------------------------------------------------
banner("Section 5.5. Empirical 80% / 95% predictive intervals, h = 1, window A (expanding)")
cov <- do.call(rbind, lapply(c("H1_ARIMA_GBM","ARIMA","H3_STL_GBM","STL_ARIMA","GBM","SNaive"), function(mdl)
  pi_emp("expanding", mdl, 1) %>% filter(date >= WIN$A[1], !is.na(lo80)) %>%
    summarise(model = mdl, n = n(), coverage80 = mean(actual >= lo80 & actual <= hi80),
              coverage95 = mean(actual >= lo95 & actual <= hi95), mean_width80 = mean(hi80 - lo80))))
print(r3(cov)); save_tab(cov, "Section5_5_interval_coverage.csv")

# =============================================================================
#  APPENDIX B. ADDITIONAL TABLES
# =============================================================================
banner("Table B.12. Complete unrestricted error-correction estimates (dependent variable: change in QIM)")
print(r3(tabB12)); cat(B12_note, "\n"); save_tab(tabB12, "TableB12_UECM.csv")

banner("Table B.13. Out-of-sample RMSE, fixed 66-month rolling window")
tabB13 <- rmse_table("rolling"); print(tabB13, right = FALSE); save_tab(tabB13, "TableB13_RMSE_rolling.csv")

banner("Table B.14. MAE, MAPE and MASE one month ahead, expanding window")
tabB14 <- acc %>% filter(tag == "expanding", h == 1) %>% select(model, window, MAE, MAPE, MASE) %>%
  pivot_wider(names_from = window, values_from = c(MAE, MAPE, MASE))
tabB14 <- tabB14[match(names(lab), tabB14$model), c("model","MAE_A","MAPE_A","MASE_A","MAE_B","MAPE_B","MASE_B")]
tabB14$model <- lab[tabB14$model]; print(as.data.frame(r3(tabB14)), right = FALSE); save_tab(tabB14, "TableB14_MAE_MAPE_MASE.csv")

banner("Table B.15. Diebold-Mariano tests of H3 against each competitor (DM statistic, Holm-adjusted p)")
tabB15 <- DM %>% filter(tag == "expanding", set == "H3", h %in% c(1, 12)) %>%
  mutate(cell = sprintf("%.2f (%.3f)", DM, p_holm)) %>% select(comparator, window, h, cell) %>%
  pivot_wider(names_from = c(window, h), values_from = cell)
tabB15 <- tabB15[match(intersect(names(lab), tabB15$comparator), tabB15$comparator), c("comparator","A_1","A_12","B_1","B_12")]
tabB15$comparator <- lab[tabB15$comparator]; print(as.data.frame(tabB15), right = FALSE)
save_tab(tabB15, "TableB15_DM.csv")

banner(paste("Finished. All tables and figures are in:", normalizePath(OUT_DIR)))
sessionInfo()
