# =============================================================================
#  "How evaluation design shapes the apparent accuracy of hybrid econometric-
#   machine learning forecasts: Pakistan's industrial production under domestic
#   energy and food price shocks"
#  INPUT : the three Excel files from the replication repository, placed in DATA_DIR
#            Industrial Production Data.xlsx, Wheat data.xlsx, Energy data.xlsx
#          and "Brent data.csv", monthly Brent crude oil spot price (USD per barrel,
#          U.S. EIA), in the same folder.
#  OUTPUT: every table and figure of the manuscript, numbered as in the manuscript
#  (Tables 1-18, A1-A14; Figures 1-3, A1). Table 16 summarises the robustness
#  checks and has no separate output.
#  RUNTIME: Sections 2-4 take about two minutes. The rolling-origin forecasts
#  (Section 5) take roughly 20-25 minutes on 8 cores and 45-60 minutes on 4 cores;
#  the additional variants of Section 5.6 add about 10 minutes on 4 cores, those of
#  Section 5.7 about 8 minutes, and the additional benchmarks and order-consistent
#  price comparisons of Section 5.8 about 4 minutes. Forecasts are cached in
#  OUT_DIR/cache, so later runs are fast.
#  Tested with R 4.3.3. Random seeds are fixed, so the numbers are reproducible.
# =============================================================================
DATA_DIR <- "C:/Users/ahmed/Desktop/IfADo Projects/Kashif Thesis new/MDPI_manuscript"  # <-- folder that contains the three Excel files
OUT_DIR  <- "C:/Users/ahmed/Desktop/IfADo Projects/Kashif Thesis new/MDPI_manuscript/results"      # <-- all tables and figures are written here
N_CORES  <- max(1, parallel::detectCores() - 1)   # cores for the rolling forecasts
B_BOOT   <- 1999           # bootstrap replications for the bounds tests 
B_MCS    <- 5000           # bootstrap replications for the model confidence set 

pkgs <- c("readxl","dplyr","tidyr","zoo","urca","tseries","lmtest","sandwich","strucchange",
          "FinTS","ARDL","forecast","gbm","randomForest","MCS","ggplot2","parallel","MSwM","quantreg")
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


# ---- Table A6: price shocks and the import-cost proxy --------------
# A "shock" is operationalised in two ways: (i) the price
# innovation, i.e. the residual of an AR(1) model for the monthly log change, and
# (ii) the net price increase of Hamilton (1996, 2003): the log excess of the price
# over its maximum in the previous 12 months. As an import-cost proxy we use the
# Brent crude oil spot price (USD per barrel, U.S. EIA; file "Brent data.csv").
banner("Table A6. Price shocks: innovations and net price increases")
br_raw <- read.csv(file.path(DATA_DIR, "Brent data.csv"))
BRENT  <- br_raw$Brent_USD_per_barrel[match(format(d$date, "%Y-%m"), br_raw$Month)]
stopifnot(!anyNA(BRENT))
# Brent row of Table 2 (descriptive statistics; Brent is loaded here because it is used only in the robustness checks)
tab2b <- rbind(tab2, data.frame(Variable = "Brent crude oil (USD/barrel)", Mean = mean(BRENT), SD = sd(BRENT),
                                CV_pct = 100 * sd(BRENT) / mean(BRENT), Skewness = skw(BRENT), Excess_kurtosis = kur(BRENT)))
cat("Table 2 with the Brent row:\n"); print(r3(tab2b)); save_tab(tab2b, "Table02_descriptives.csv")
# AR(1) innovations of the log change, coefficients estimated on the rows in 'rows' only
price_innov <- function(x, rows = seq_along(x)) {
  dl <- c(NA, diff(log(x))); dl1 <- c(NA, dl[-length(dl)])
  ok <- intersect(rows, which(!is.na(dl) & !is.na(dl1)))
  b <- coef(lm(dl[ok] ~ dl1[ok])); as.numeric(dl - (b[1] + b[2] * dl1))
}
# Hamilton net price increase over the previous 12 months (shorter history at the start)
nopi <- function(x) { lx <- log(x)
  sapply(seq_along(lx), function(t) if (t == 1) 0 else max(0, lx[t] - max(lx[max(1, t - 12):(t - 1)]))) }
PX <- cbind(X, BRENT = BRENT)
tabA6 <- do.call(rbind, lapply(colnames(PX), function(k) {
  x <- PX[, k]; dl <- c(NA, diff(log(x))); u <- price_innov(x); np <- nopi(x)
  big <- which(abs(u) > 2 * sd(u, na.rm = TRUE)); im <- which.max(u); inp <- which.max(np)
  data.frame(Series = k, sd_dlog_pct = 100 * sd(dl, na.rm = TRUE), months_large_innov = length(big),
             largest_innov_pct = 100 * u[im], largest_innov_month = format(d$date[im], "%b %Y"),
             months_NOPI_pos = sum(np > 0), largest_NOPI_pct = 100 * np[inp],
             largest_NOPI_month = format(d$date[inp], "%b %Y"),
             large_innov_months = paste(format(d$date[big], "%b %Y"), collapse = "; "), row.names = NULL) }))
print(r3(tabA6[, -ncol(tabA6)]), right = FALSE)
cat("Months with innovations larger than 2 SD:\n"); for (i in seq_len(nrow(tabA6))) cat(" ", tabA6$Series[i], ":", tabA6$large_innov_months[i], "\n")
cat(sprintf("Correlation of log Brent with log fuel price: %.2f (levels), %.2f (monthly changes)\n",
            cor(log(BRENT), log(X[, "FUEL"])), cor(diff(log(BRENT)), diff(log(X[, "FUEL"])))))
save_tab(tabA6, "TableA6_price_shocks.csv")

# =============================================================================
#  SECTION 3. PHASE I: LINEAR MODELS AND THE LONG-RUN QUESTION
# =============================================================================
# ---- Table 3 --------------------------------------------------------------
banner("Table 3. Unit root tests in levels and first differences")
tab3 <- do.call(rbind, lapply(c("QIM", PR), function(v) {
  x <- d[[v]]; dx <- diff(x)
  a1 <- ur.df(x, type = "trend", lags = 9, selectlags = "AIC"); a2 <- ur.df(dx, type = "drift", lags = 9, selectlags = "AIC")   # AIC over 1-9 lags
  k1 <- ur.kpss(x, type = "tau", lags = "short");    k2 <- ur.kpss(dx, type = "mu", lags = "short")
  cls <- if (a1@teststat[1] < a1@cval[1, 2] && k1@teststat <= k1@cval[2]) "I(0)" else
         if (a2@teststat[1] < a2@cval[1, 2] && k2@teststat <= k2@cval[2]) "I(1)" else "ambiguous"
  nl <- function(a) sum(grepl("z.diff.lag", rownames(a@testreg$coefficients)))
  data.frame(Variable = v, ADF_level = a1@teststat[1], ADF_lags_level = nl(a1), KPSS_level = k1@teststat,
             ADF_diff = a2@teststat[1], ADF_lags_diff = nl(a2), KPSS_diff = k2@teststat, Classification = cls)
}))
print(r3(tab3)); save_tab(tab3, "Table03_unit_roots.csv")
cat("5% critical values: ADF levels -3.43, KPSS levels 0.146, ADF differences -2.88, KPSS differences 0.463\n")
cat("KPSS bandwidth (lags = 'short'):", ur.kpss(y, type = "tau", lags = "short")@lag, "\n")
# Sensitivity of the classification of QIM: up to 12 lags, and monthly dummies in the ADF regression
a12 <- ur.df(y, type = "trend", lags = 12, selectlags = "AIC")
cat(sprintf("QIM, ADF with AIC over 1-12 lags: %.3f (%d lags)\n", a12@teststat[1], sum(grepl("z.diff.lag", rownames(a12@testreg$coefficients)))))
adf_sd <- function(x, maxl = 9) { n <- length(x); dx <- c(NA, diff(x)); best <- NULL
  for (k in 1:maxl) { Z <- data.frame(dx = dx, x1 = c(NA, x[-n]), tt = seq_len(n), m = factor(mon))
    for (j in 1:k) Z[[paste0("l", j)]] <- c(rep(NA, j), dx[1:(n - j)])
    Z <- Z[(maxl + 2):n, ]; f <- lm(dx ~ ., data = Z)
    if (is.null(best) || AIC(f) < best$aic) best <- list(aic = AIC(f), k = k, t = coef(summary(f))["x1", 3]) }
  best }
r_sd <- adf_sd(y); cat(sprintf("QIM, ADF with trend and monthly dummies (AIC over 1-9 lags): %.3f (%d lags)\n", r_sd$t, r_sd$k))

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


# ---- Table 4, panel B, and Table 5: seasonal specification as the preferred model ----
# QIM is strongly seasonal, so the preferred ARDL includes monthly (and COVID) dummies and its
# lag order is selected by AIC WITH these dummies as fixed regressors. The UECM and the
# McNown et al. (2018) bootstrap are generalised here to p > 1 lags of QIM (Delta-y lags are
# regenerated recursively in each bootstrap sample) and to an optional deterministic trend.
banner("Table 4, panel B. Preferred seasonal specification, deterministic trend and log prices")
build_g <- function(y, X, q, p = 1, extra = NULL) {
  Z <- data.frame(dy = Df(y), y1 = Lg(y, 1))
  if (p > 1) for (j in 1:(p - 1)) Z[[paste0("dy_l", j)]] <- Lg(Df(y), j)
  for (k in PR) Z[[paste0(k, "_l1")]] <- if (q[k] == 0) X[, k] else Lg(X[, k], 1)
  for (k in PR) { if (q[k] == 0) next
    for (j in 0:(q[k] - 1)) Z[[paste0("d", k, "_", j)]] <- Lg(Df(X[, k]), j) }
  if (!is.null(extra)) Z <- cbind(Z, extra)
  Z
}
boot_null_g <- function(y, X, q, p, extra, drop, stat, B = B_BOOT) {
  Z <- build_g(y, X, q, p, extra); cc <- complete.cases(Z); idx <- which(cc)
  M <- cbind("(Intercept)" = 1, as.matrix(Z[cc, -1])); dy0 <- Z$dy[cc]
  keep <- setdiff(colnames(M), drop)
  br <- qr.coef(qr(M[, keep]), dy0); br[is.na(br)] <- 0
  e <- dy0 - M[, keep] %*% br; e <- as.numeric(e - mean(e))
  dyl <- grep("^dy_l", keep, value = TRUE); fixed <- setdiff(keep, c("y1", dyl))
  mu  <- as.numeric(M[, fixed, drop = FALSE] %*% br[fixed])
  lam <- if ("y1" %in% keep) br["y1"] else 0; phi <- br[dyl]
  res <- numeric(B)
  for (b in 1:B) {
    ys <- y; es <- sample(e, length(idx), replace = TRUE)
    for (i in seq_along(idx)) { t <- idx[i]; dyt <- mu[i] + lam * ys[t - 1] + es[i]
      if (length(phi)) for (j in seq_along(phi)) dyt <- dyt + phi[j] * (ys[t - j] - ys[t - j - 1])
      ys[t] <- ys[t - 1] + dyt }
    Ms <- M; Ms[, "y1"] <- ys[idx - 1]
    for (j in seq_along(dyl)) Ms[, paste0("dy_l", j)] <- ys[idx - j] - ys[idx - j - 1]
    res[b] <- fast_stats(diff(ys)[idx - 1], Ms)[stat]
  }
  res
}
bootrep_g <- function(y, X, q, p, extra, label) {
  f <- lm(dy ~ ., data = build_g(y, X, q, p, extra)); s <- stats3(f)
  bF <- boot_null_g(y, X, q, p, extra, c("y1", lv), "F_overall")
  bt <- boot_null_g(y, X, q, p, extra, "y1", "t_dep")
  bE <- boot_null_g(y, X, q, p, extra, lv, "F_exog")
  b <- coef(f); Vv <- vcov(f); g <- -b[lv] / b["y1"]
  se <- sapply(lv, function(k) { gr <- setNames(rep(0, length(b)), names(b)); gr[k] <- -1 / b["y1"]
    gr["y1"] <- b[k] / b["y1"]^2; sqrt(t(gr) %*% Vv %*% gr) })
  data.frame(Specification = label, ARDL_order = paste(c(p, q), collapse = ","), F_ov = s[1], boot_p_F = mean(bF >= s[1]),
             t_y1 = s[2], boot_p_t = mean(bt <= s[2]), F_exog = s[3], boot_p_Fexog = mean(bE >= s[3]),
             boot_cv5_Fexog = quantile(bE, .95), lambda0 = coef(f)["y1"],
             LR_wheat = g[1], LR_fuel = g[2], LR_elec = g[3], LR_gas = g[4],
             p_wheat = 2 * pt(-abs(g[1] / se[1]), df.residual(f)), p_fuel = 2 * pt(-abs(g[2] / se[2]), df.residual(f)),
             p_elec = 2 * pt(-abs(g[3] / se[3]), df.residual(f)), p_gas = 2 * pt(-abs(g[4] / se[4]), df.residual(f)),
             n = nobs(f), row.names = NULL)
}
sel_order <- function(dat, dummies) {
  fm <- as.formula(paste("QIM ~ WHEAT + FUEL + ENERGY + GAS |", paste(colnames(dummies), collapse = " + ")))
  o <- auto_ardl(fm, data = cbind(dat, dummies), max_order = 6)$best_order
  list(p = o[1], q = setNames(o[-1], PR))
}
set.seed(20261112)
TRD <- cbind(trend = seq_len(N)); Xlog <- log(X); dlog <- d; dlog[PR] <- log(dlog[PR])
o_lev <- sel_order(d, cbind(SD, COV)); o_log <- sel_order(dlog, cbind(SD, COV))
cat("AIC order with monthly and COVID dummies: levels", o_lev$p, o_lev$q, "; log prices", o_log$p, o_log$q, "\n")
tab4b <- rbind(
  bootrep_g(y, X, q, 1, NULL, "Check: baseline via generalised code"),
  bootrep_g(y, X, o_lev$q, o_lev$p, cbind(SD, COV), "Preferred: monthly + COVID dummies, order re-selected"),
  bootrep_g(y, X, o_lev$q, o_lev$p, cbind(SD, COV, TRD), "Preferred + linear trend (Case V)"),
  bootrep_g(y, Xlog, q, 1, NULL, "Log prices: baseline order"),
  bootrep_g(y, Xlog, q, 1, cbind(SD, COV), "Log prices: + monthly and COVID dummies"),
  bootrep_g(y, Xlog, o_log$q, o_log$p, cbind(SD, COV), "Log prices: monthly + COVID dummies, order re-selected"))
print(r3(tab4b[, c("Specification","ARDL_order","F_ov","boot_p_F","t_y1","boot_p_t","F_exog","boot_p_Fexog","lambda0")]), right = FALSE)
cat("\nLong-run coefficients (delta method; descriptive when F_exog does not reject):\n")
print(r3(tab4b[, c("Specification","LR_wheat","p_wheat","LR_fuel","p_fuel","LR_elec","p_elec","LR_gas","p_gas")]), right = FALSE)
save_tab(tab4b, "Table04_panelB_preferred_trend_logs.csv")

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

# ---- Table 6: robustness to the April 2020 observation -------------
# Alternative treatments of the COVID-19 lockdown:
# impulse dummies, a lag order re-selected with the dummies as fixed regressors,
# the post-lockdown subsample, and the stability tests without April-May 2020.
banner("Table 6. Robustness of the ARDL results to the April 2020 observation")
set.seed(20261111)
D_apr <- COV[, "covid_apr20", drop = FALSE]                  # April 2020 impulse dummy
lr_tab <- function(f, lab) { o <- lr_of(f); data.frame(Specification = lab, Variable = o$Variable,
                                                       Coef = o$Coef, p = o$p, row.names = NULL) }
stab <- function(Z, drop_rows = NULL) {      # sup-F (Andrews) and recursive CUSUM on the UECM
  Zc <- Z[complete.cases(Z), ]; dts <- d$date[complete.cases(Z)]
  if (!is.null(drop_rows)) { keep <- !(dts %in% drop_rows); Zc <- Zc[keep, ]; dts <- dts[keep] }
  Zc <- Zc[, !grepl("^covid|^m[0-9]", names(Zc)), drop = FALSE]
  rc <- sctest(efp(dy ~ ., data = Zc, type = "Rec-CUSUM"))
  fs <- suppressWarnings(Fstats(dy ~ ., data = Zc, from = 0.15)); sf <- sctest(fs, type = "supF")
  data.frame(CUSUM_S = unname(rc$statistic), CUSUM_p = rc$p.value, supF = unname(sf$statistic),
             supF_p = sf$p.value, break_date = format(dts[fs$breakpoint], "%b %Y"), n = nrow(Zc))
}
rev_b <- list(); rev_lr <- list(); rev_st <- list()
for (s in list(list("Baseline", NULL), list("+ April 2020 dummy", D_apr),
               list("+ April and May 2020 dummies", COV), list("+ monthly dummies", SD),
               list("+ monthly and COVID dummies", cbind(SD, COV)))) {
  rev_b[[length(rev_b) + 1]]   <- bootrep(s[[2]], s[[1]])
  rev_lr[[length(rev_lr) + 1]] <- lr_tab(fit_uecm(build(y, X, s[[2]])), s[[1]])
}
# April-May 2020 removed from the UECM (equivalent to two impulse dummies, since p = 1)
rev_st[[1]] <- cbind(Specification = "Baseline", stab(build(y, X)))
rev_st[[2]] <- cbind(Specification = "April-May 2020 removed",
                     stab(build(y, X), as.Date(c("2020-04-01", "2020-05-01"))))
# Post-lockdown subsample, July 2020 - November 2025 (same ARDL order). The bootstrap and
# UECM functions use y, X and d from the workspace, so these are swapped temporarily.
y_full <- y; X_full <- X; d_full <- d; SD_full <- SD
sub <- d$date >= as.Date("2020-07-01")
y <- y_full[sub]; X <- X_full[sub, ]; d <- d_full[sub, ]; SD <- SD_full[sub, ]
rev_b[[length(rev_b) + 1]]   <- bootrep(NULL, "Post-lockdown subsample (Jul 2020-Nov 2025)")
rev_lr[[length(rev_lr) + 1]] <- lr_tab(fit_uecm(build(y, X)), "Post-lockdown subsample (Jul 2020-Nov 2025)")
rev_b[[length(rev_b) + 1]]   <- bootrep(SD, "Post-lockdown subsample + monthly dummies")
rev_lr[[length(rev_lr) + 1]] <- lr_tab(fit_uecm(build(y, X, SD)), "Post-lockdown subsample + monthly dummies")
rev_st[[3]] <- cbind(Specification = "Post-lockdown subsample", stab(build(y, X)))
y <- y_full; X <- X_full; d <- d_full; SD <- SD_full; rm(y_full, X_full, d_full, SD_full, sub)
# Lag order re-selected by AIC with the COVID dummies as fixed regressors
dd_cov <- cbind(d, COV)
sel_c <- auto_ardl(QIM ~ WHEAT + FUEL + ENERGY + GAS | covid_apr20 + covid_may20, data = dd_cov, max_order = 6)
bf_c <- bounds_f_test(sel_c$best_model, case = 3); bt_c <- bounds_t_test(sel_c$best_model, case = 3)
tab6_b <- do.call(rbind, rev_b); tab6_lr <- do.call(rbind, rev_lr); tab6_st <- do.call(rbind, rev_st)
print(r3(tab6_b[, c("Specification","F_ov","boot_p_F","t_y1","boot_p_t","F_exog","boot_p_Fexog","lambda0")]))
cat("\nLong-run coefficients:\n"); print(r3(tab6_lr))
cat("\nParameter stability:\n"); print(r3(tab6_st))
cat(sprintf("\nRe-selected order with COVID dummies as fixed regressors: (%s)\nBounds F = %.3f (asymptotic p = %.3f); t = %.3f (asymptotic p = %.3f)\n",
            paste(sel_c$best_order, collapse = ","), bf_c$statistic, bf_c$p.value, bt_c$statistic, bt_c$p.value))
print(r3(multipliers(sel_c$best_model)))
cat("Note: Table 6 takes the Baseline and '+ monthly and COVID' rows from Tables 4-5; the other rows from this run.\n")
save_tab(tab6_b, "Table06a_covid_bounds.csv"); save_tab(tab6_lr, "Table06b_covid_longrun.csv")
save_tab(tab6_st, "Table06c_covid_stability.csv")

# ---- Table A7 and Figure A1: structural breaks and regimes ---------
# (a) Bai-Perron multiple-break tests, (b) a two-regime
# Markov-switching model, (c) rolling-window (60-month) ARDL estimates.
banner("Table A7. Structural breaks and regimes")
# (a) Bai-Perron: mean shifts in seasonally adjusted QIM, and breaks in the full UECM
ysa <- y - ave(y, mon) + mean(y)                        # QIM minus its calendar-month mean
bp_mean <- breakpoints(ysa ~ 1, h = 0.10)              # minimum segment: 11 months
bp_mean_n <- length(na.omit(breakpoints(bp_mean)$breakpoints))
cat("Bai-Perron, mean shifts in seasonally adjusted QIM: BIC selects", bp_mean_n, "break(s)\n")
if (bp_mean_n > 0) { print(format(d$date[breakpoints(bp_mean)$breakpoints], "%b %Y"))
  ci_m <- tryCatch(confint(bp_mean)$confint, error = function(e) NULL)
  if (!is.null(ci_m)) print(apply(ci_m, 2, function(i) format(d$date[i], "%b %Y"))) }
seg_means <- tapply(ysa, breakfactor(bp_mean), mean); print(round(seg_means, 1))
cat("sup-F test for a mean shift (15% trimming):\n"); print(sctest(Fstats(ysa ~ 1, from = 0.15), type = "supF"))
# UECM with April-May 2020 removed (equivalent to impulse dummies); minimum segment 20% (21 obs)
Zb <- build(y, X); dts_b <- d$date[complete.cases(Zb)]; Zb <- Zb[complete.cases(Zb), ]
keep_b <- !(dts_b %in% as.Date(c("2020-04-01", "2020-05-01"))); Zb <- Zb[keep_b, ]; dts_b <- dts_b[keep_b]
bp_uecm <- breakpoints(dy ~ ., data = Zb, h = 0.20)
bp_uecm_n <- length(na.omit(breakpoints(bp_uecm)$breakpoints))
cat("Bai-Perron, all UECM coefficients (April-May 2020 removed): BIC selects", bp_uecm_n, "break(s)",
    if (bp_uecm_n > 0) paste("at", paste(format(dts_b[breakpoints(bp_uecm)$breakpoints], "%b %Y"), collapse = ", ")), "\n")
print(summary(bp_uecm)$RSS["BIC", ])
# (b) Markov-switching AR(1) for seasonally adjusted QIM (Hamilton, 1989): two regimes with
# switching mean and variance and a common AR coefficient; estimated on the series as
# recorded and on an outlier-adjusted series (April-May 2020 interpolated linearly)
ia <- which(d$date %in% as.Date(c("2020-04-01", "2020-05-01")))
ysa_adj <- ysa; ysa_adj[ia] <- approx(c(ia[1] - 1, ia[2] + 1), ysa[c(ia[1] - 1, ia[2] + 1)], xout = ia)$y
ms_summary <- function(z, label) {
  set.seed(20261111)
  m <- MSwM::msmFit(lm(z ~ 1), k = 2, p = 1, sw = c(TRUE, FALSE, TRUE), control = list(parallel = FALSE))
  sp <- m@Fit@smoProb[-1, , drop = FALSE]; dts <- d$date[-1]; P <- m@transMat
  out <- do.call(rbind, lapply(1:2, function(r) { inr <- sp[, r] > 0.5; rr <- rle(inr)
    data.frame(Series = label, Regime = r, Mean = m@Coef[r, 1] / (1 - m@Coef[r, 2]), AR1 = m@Coef[r, 2], SD = m@std[r],
               Stay_prob = P[r, r], Expected_duration = 1 / (1 - P[r, r]), Months = sum(inr),
               Spells = sum(rr$values), row.names = NULL) }))
  r_small <- which.min(colSums(sp > 0.5)); inr <- sp[, r_small] > 0.5; rr <- rle(inr); e <- cumsum(rr$lengths); s0 <- e - rr$lengths + 1
  attr(out, "periods") <- paste(sapply(which(rr$values), function(i) paste(format(dts[s0[i]], "%b %Y"), format(dts[e[i]], "%b %Y"), sep = "-")), collapse = "; ")
  out }
ms1 <- ms_summary(ysa, "As recorded"); ms2 <- ms_summary(ysa_adj, "April-May 2020 interpolated")
tabA7_ms <- rbind(ms1, ms2); print(r3(tabA7_ms))
cat("Smaller regime, as recorded:", attr(ms1, "periods"), "\nSmaller regime, outlier-adjusted:", attr(ms2, "periods"), "\n")
# (c) Rolling 60-month windows of the UECM (fixed order), with and without monthly dummies;
# April and May 2020 impulse dummies are included whenever these months are in the window
Wr <- 60
roll <- do.call(rbind, lapply(Wr:N, function(e) { ix <- (e - Wr + 1):e
  cc_ix <- ix[-seq_len(max(q))]                     # rows that survive the lags of the UECM
  cv <- COV[ix, colSums(COV[cc_ix, , drop = FALSE]) > 0, drop = FALSE]
  do.call(rbind, lapply(c("no monthly dummies", "monthly dummies"), function(sp_) {
    ex <- if (sp_ == "monthly dummies") cbind(SD[ix, ], cv) else cv
    if (ncol(ex) == 0) ex <- NULL
    f <- fit_uecm(build(y[ix], X[ix, ], ex)); lr <- lr_of(f); s3 <- stats3(f)
    data.frame(end = d$date[e], spec = sp_, lambda0 = coef(f)["y1"], theta_elec = lr$Coef[lr$Variable == "ENERGY"],
               p_elec = lr$p[lr$Variable == "ENERGY"], F_exog = s3["F_exog"],
               p_Fexog_asy = pf(s3["F_exog"], 4, df.residual(f), lower.tail = FALSE), row.names = NULL) })) }))
tabA7_roll <- roll %>% group_by(spec) %>% summarise(windows = n(), lambda0_min = min(lambda0), lambda0_max = max(lambda0),
  elec_min = min(theta_elec), elec_max = max(theta_elec), elec_sig5_share = mean(p_elec < 0.05),
  Fexog_sig5_share = mean(p_Fexog_asy < 0.05), .groups = "drop")
print(as.data.frame(r3(tabA7_roll)))
save_tab(roll, "TableA7_rolling_ARDL.csv"); save_tab(tabA7_ms, "TableA7_markov_switching.csv")
pA1 <- ggplot(roll %>% select(end, spec, lambda0, theta_elec) %>% pivot_longer(c(lambda0, theta_elec)) %>%
                mutate(name = factor(name, levels = c("lambda0", "theta_elec"),
                       labels = c("Adjustment~coefficient~(lambda[0])", "'Long-run electricity coefficient'"))),
              aes(end, value, colour = spec, linetype = spec)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.6) + facet_wrap(~name, ncol = 1, scales = "free_y", labeller = label_parsed) +
  scale_colour_manual(values = c("#1f77b4", "#d62728"), name = NULL) + scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") + labs(x = "End of the 60-month estimation window", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey93"))
ggsave(file.path(OUT_DIR, "FigureA1_rolling_ARDL.pdf"), pA1, width = 6.5, height = 4.8)

# ---- Table 7 --------------------------------------------------------------
banner("Table 7. Residual diagnostics for the ARDL model")
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
save_tab(tab6a, "Table07a_diagnostics.csv"); save_tab(tab6b, "Table07b_BDS.csv")

# Appendix Table A1 is computed here and printed at the end of the script
ctab <- coeftest(f0); hac <- coeftest(f0, vcov = NeweyWest(f0, lag = 4, prewhite = FALSE))
tabB12 <- data.frame(Term = rownames(ctab), Coef = ctab[, 1], OLS_se = ctab[, 2], OLS_p = ctab[, 4],
                     HAC_se = hac[, 2], HAC_p = hac[, 4], row.names = NULL)
B12_note <- sprintf("N = %d; R2 = %.3f; adjusted R2 = %.3f", nobs(f0), summary(f0)$r.squared, summary(f0)$adj.r.squared)

# =============================================================================
#  SECTION 4. PHASE II: FORECASTING DESIGN
# =============================================================================
banner("Table 8. Forecasting methods (15 feasible methods and one infeasible reference)")
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

# ---- Table 9 --------------------------------------------------------------
banner("Table 9. Out-of-sample RMSE, expanding window  (A = Jan 2024-Nov 2025, B = Jan 2023-Nov 2025)")
tab8 <- rmse_table("expanding"); print(tab8, right = FALSE); save_tab(tab8, "Table09_RMSE_expanding.csv")
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

# ---- Table 13 -------------------------------------------------------------
banner("Table 13. Ablations: RMSE ratio with / without the ingredient (expanding window)")
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
print(r3(tab9), right = FALSE); save_tab(tab9, "Table13_ablations.csv")
cat("DM tests for the H3-based pairs at h = 1 (note to Table 13):\n")
print(as.data.frame(r3(DM %>% filter(tag == "expanding", set == "H3", h == 1,
                                     comparator %in% c("STL_ARIMA","H1_ARIMA_GBM","H3_STL_GBM_noCom")) %>%
                         select(window, comparator, DM, p, p_holm))))

# ---- Table 15 -------------------------------------------------------------
banner("Table 15. Robustness to the estimation scheme and the number of lags")
tab10 <- rk %>% filter(tag != "expanding") %>%
  mutate(tag = factor(tag, levels = c("rolling","expanding_L3","expanding_L12"),
                      labels = c("Rolling (66 months), 6 lags","Expanding, 3 lags","Expanding, 12 lags"))) %>%
  arrange(tag, window, h) %>% select(Scheme = tag, Window = window, h, Lowest_RMSE_method = best, RMSE = best_RMSE,
                                     H3_RMSE, H3_rank, SNaive_rank, In_MCS90 = in_MCS90)
print(as.data.frame(r3(tab10)), right = FALSE); save_tab(tab10, "Table15_robustness.csv")

# ---- Table 12 -------------------------------------------------------------
banner("Table 12. From a single train/test split to the leakage-free design (RMSE, Jan 2024 - Nov 2025)")
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
print(r3(tab11), right = FALSE); save_tab(tab11, "Table12_bridge.csv")
cat(sprintf("Row (c): STL-ARIMA component alone has RMSE %.2f.  Row (a): H3 is %.0f%% below ARIMA and %.0f%% below GBM.\n",
            rm_(base_c), 100 * (1 - rm_(h3_a) / rm_(p_ns)), 100 * (1 - rm_(h3_a) / rm_(p_gbm))))

# ---- Section 5.5 ------------------------------------------------------------
banner("Section 5.5. Empirical 80% / 95% predictive intervals, h = 1, window A (expanding; see also Table 18)")
cov <- do.call(rbind, lapply(c("H1_ARIMA_GBM","ARIMA","H3_STL_GBM","STL_ARIMA","GBM","SNaive"), function(mdl)
  pi_emp("expanding", mdl, 1) %>% filter(date >= WIN$A[1], !is.na(lo80)) %>%
    summarise(model = mdl, n = n(), coverage80 = mean(actual >= lo80 & actual <= hi80),
              coverage95 = mean(actual >= lo95 & actual <= hi95), mean_width80 = mean(hi80 - lo80))))
print(r3(cov)); save_tab(cov, "Section5_5_interval_coverage.csv")

# =============================================================================
#  SECTION 5.6. ADDITIONAL ANALYSES: Tables 10, 17, 18 and A5
# =============================================================================
# Additional methods, computed for the expanding window with six lags only:
#   H1h-H4h        hybrids whose corrector is trained on in-sample h-step errors
#   RF_alt, GBM_d1, GBM_fixed, H1_d1, H1_fixed   alternative hyperparameters
#   H1_noCom, RF_noCom                             further price ablations
#   ARIMAX_pricefc ARIMAX with auto.arima-forecast price paths
# The baseline methods are not recomputed; their forecasts are taken from RALL. None of the
# additional methods enters the 15 feasible methods, the model confidence sets or Tables 9-15.
rf_alt_fit_pred <- function(Xtr, ytr, Xte, seed = 42) { set.seed(seed)
  as.numeric(predict(randomForest(Xtr, ytr, ntree = 500, mtry = max(1, floor(ncol(Xtr) / 3)), nodesize = 5), Xte)) }
gbm_fixed_pred <- function(Xtr, ytr, Xte, depth = 3, M = 100, seed = 42) { set.seed(seed)
  df <- data.frame(yv = ytr, Xtr); dte <- data.frame(Xte); colnames(dte) <- colnames(df)[-1]
  f <- gbm(yv ~ ., data = df, distribution = "gaussian", n.trees = M, interaction.depth = depth,
           shrinkage = 0.05, bag.fraction = 0.8, n.minobsinnode = 5, verbose = FALSE)
  as.numeric(predict(f, dte, n.trees = M)) }
one_origin_rev <- function(o, LAGS = 6) {
  tr <- 1:o
  st <- c(as.numeric(format(d$date[tr[1]], "%Y")), mon[tr[1]])
  yt <- ts(y[tr], start = st, frequency = 12); hmax <- min(12, N - o)
  aa <- function(z, ...) auto.arima(z, ic = "aic", stepwise = TRUE, approximation = FALSE, ...)
  fA <- aa(yt, seasonal = TRUE); pA <- as.numeric(forecast(fA, h = hmax)$mean)
  fE <- ets(yt);                 pE <- as.numeric(forecast(fE, h = hmax)$mean)
  sl <- stl(yt, s.window = "periodic", robust = TRUE); seas <- as.numeric(sl$time.series[, "seasonal"])
  sa <- yt - seas; fS <- aa(sa, seasonal = FALSE)
  seas_f <- sapply(1:hmax, function(h) seas[length(seas) - 12 + ((h - 1) %% 12) + 1])
  pS <- as.numeric(forecast(fS, h = hmax)$mean) + seas_f
  fX <- aa(yt, xreg = X[tr, ], seasonal = TRUE)
  Xpf <- sapply(PR, function(k) as.numeric(forecast(auto.arima(ts(X[tr, k], start = st, frequency = 12),
                                       seasonal = FALSE, ic = "aic"), h = hmax)$mean))
  Xpf <- matrix(Xpf, nrow = hmax, dimnames = list(NULL, PR))     # forecast price paths, training data only
  pXp <- as.numeric(forecast(fX, xreg = Xpf)$mean)
  eA <- as.numeric(y[tr] - fitted(fA)); eE <- as.numeric(y[tr] - fitted(fE)); eS <- as.numeric(sa - fitted(fS))
  res <- list()
  for (h in H[H <= hmax]) {
    t_tr <- tr[tr - h - (LAGS - 1) >= 1]
    tt <- o + h; idx <- match(t_tr, tr)
    F1 <- fmat(t_tr, h, LAGS); f1 <- rbind(feat(tt, h, LAGS))
    F0 <- fmat(t_tr, h, LAGS, com = FALSE); f0_ <- rbind(feat(tt, h, LAGS, com = FALSE))
    # in-sample h-step errors: h-step forecasts of the training data with the parameters
    # estimated on the current window (identical to the one-step residuals for h = 1)
    if (h == 1) { eAh <- eA; eEh <- eE; eSh <- eS } else {
      eAh <- as.numeric(y[tr] - fitted(fA, h = h)); eEh <- as.numeric(y[tr] - fitted(fE, h = h))
      eSh <- as.numeric(sa - fitted(fS, h = h)) }
    corr_h <- function(e, learner = "gbm") { k <- !is.na(e[idx])
      if (learner == "gbm") gbm_fit_pred(F1[k, , drop = FALSE], e[idx][k], f1, depth = 2, ntmax = 300)
      else rf_fit_pred(F1[k, , drop = FALSE], e[idx][k], f1) }
    r <- c(H1h = pA[h] + corr_h(eAh), H2h = pA[h] + corr_h(eAh, "rf"),
           H3h = pS[h] + corr_h(eSh), H4h = pE[h] + corr_h(eEh),
           H1_noCom = pA[h] + gbm_fit_pred(F0, eA[idx], f0_, depth = 2, ntmax = 300),
           RF_noCom = rf_fit_pred(F0, y[t_tr], f0_),
           RF_alt = rf_alt_fit_pred(F1, y[t_tr], f1),
           GBM_d1 = gbm_fit_pred(F1, y[t_tr], f1, depth = 1, ntmax = 1000),
           GBM_fixed = gbm_fixed_pred(F1, y[t_tr], f1),
           H1_d1 = pA[h] + gbm_fit_pred(F1, eA[idx], f1, depth = 1, ntmax = 300),
           H1_fixed = pA[h] + gbm_fixed_pred(F1, eA[idx], f1, depth = 2, M = 50),
           ARIMAX_pricefc = pXp[h])
    res[[length(res) + 1]] <- data.frame(origin = o, h = h, target = tt, date = d$date[tt], actual = y[tt],
                                         model = names(r), fc = as.numeric(r), n_train = length(t_tr))
  }
  do.call(rbind, res)
}
run_rev <- function() {
  cf <- file.path(OUT_DIR, "cache", "rolling_expanding_revision.rds")
  if (file.exists(cf)) { cat("  using cached forecasts: additional variants\n"); return(readRDS(cf)) }
  cat(sprintf("  computing additional variants (%d origins, %d core(s)) ...\n", length(origins), N_CORES)); t0 <- Sys.time()
  if (N_CORES > 1 && .Platform$OS.type == "windows") {
    cl <- makeCluster(N_CORES); on.exit(stopCluster(cl))
    clusterEvalQ(cl, suppressPackageStartupMessages({ library(forecast); library(gbm); library(randomForest) }))
    clusterExport(cl, c("one_origin_rev","feat","fmat","gbm_fit_pred","rf_fit_pred","rf_alt_fit_pred","gbm_fixed_pred",
                        "d","y","X","mon","N","PR","H"), envir = environment())
    R <- parLapply(cl, origins, one_origin_rev)
  } else if (N_CORES > 1) {
    R <- mclapply(origins, one_origin_rev, mc.cores = N_CORES, mc.preschedule = FALSE)
  } else R <- lapply(origins, one_origin_rev)
  R <- do.call(rbind, R); R$tag <- "expanding"; R$err <- R$actual - R$fc
  saveRDS(R, cf); cat("  done in", format(round(Sys.time() - t0, 1)), "\n"); R
}
banner("Section 5.6. Additional analyses: rolling-origin forecasts of the additional variants (expanding, 6 lags)")
RREV <- run_rev()
REX  <- rbind(RALL[RALL$tag == "expanding", c("origin","h","target","date","actual","model","fc","tag","err")],
              RREV[, c("origin","h","target","date","actual","model","fc","tag","err")])
inwin <- function(df, w) df[df$date >= WIN[[w]][1] & df$date <= WIN[[w]][2], ]
accR <- do.call(rbind, lapply(names(WIN), function(w) inwin(REX, w) %>% group_by(model, h) %>%
  summarise(window = w, n = n(), RMSE = sqrt(mean(err^2)), .groups = "drop")))
rmse_rev <- accR %>% pivot_wider(id_cols = model, names_from = c(window, h), values_from = RMSE)
rmse_rev <- as.data.frame(rmse_rev[, c("model","A_1","A_3","A_6","A_12","B_1","B_3","B_6","B_12")])
cat("RMSE of all methods, including the additional variants:\n"); print(r3(rmse_rev), row.names = FALSE)
save_tab(rmse_rev, "Rev_RMSE_all_methods.csv")

E_of <- function(w, hh, models) { s <- inwin(REX, w) %>% filter(h == hh, model %in% models)
  list(e = s %>% select(date, model, err) %>% pivot_wider(names_from = model, values_from = err) %>% arrange(date),
       f = s %>% select(date, model, fc) %>% pivot_wider(names_from = model, values_from = fc) %>% arrange(date)) }
# Clark-West (2007) adjusted MSPE test; H1: the larger (nested) model is more accurate;
# Newey-West variance with h - 1 lags; one-sided p-value from t(n - 1)
cw_test <- function(e_small, e_large, f_small, f_large, h) {
  fd <- e_small^2 - (e_large^2 - (f_small - f_large)^2); n <- length(fd)
  m <- lm(fd ~ 1); v <- if (h > 1) NeweyWest(m, lag = h - 1, prewhite = FALSE)[1, 1] else vcov(m)[1, 1]
  tstat <- unname(coef(m)[1] / sqrt(v)); c(tstat, pt(tstat, df = n - 1, lower.tail = FALSE)) }
# RMSE ratio of p[1] to p[2] with a DM or Clark-West test, for every horizon and window
pair_table <- function(pairs, test = c("dm","cw")) {
  test <- match.arg(test)
  do.call(rbind, lapply(pairs, function(p) do.call(rbind, lapply(names(WIN), function(w) do.call(rbind, lapply(H, function(hh) {
    L <- E_of(w, hh, p[1:2]); e1 <- L$e[[p[1]]]; e2 <- L$e[[p[2]]]
    tt <- if (test == "dm") dm_hln(e1, e2, hh) else cw_test(e2, e1, L$f[[p[2]]], L$f[[p[1]]], hh)
    data.frame(comparison = p[3], window = w, h = hh, ratio = sqrt(mean(e1^2)) / sqrt(mean(e2^2)),
               stat = tt[1], p = tt[2]) }))))))
}
wide_cells <- function(tb) { tb$cell <- ifelse(is.na(tb$p), sprintf("%.2f", tb$ratio), sprintf("%.2f (%.3f)", tb$ratio, tb$p))
  as.data.frame(tb %>% select(comparison, window, h, cell) %>% pivot_wider(names_from = c(window, h), values_from = cell)) }

# ---- Table 10: lowest-RMSE method and its statistical support ----
banner("Table 10. The method with the lowest RMSE and the statistical support for its advantage")
tab10r <- do.call(rbind, lapply(names(WIN), function(w) do.call(rbind, lapply(H, function(hh) {
  a <- accR %>% filter(window == w, h == hh, model %in% feasible) %>% arrange(RMSE)
  best <- a$model[1]; E <- E_of(w, hh, feasible)$e
  cmp <- setdiff(feasible, "SNaive")
  holm <- setNames(p.adjust(sapply(cmp, function(m) dm_hln(E[[m]], E$SNaive, hh)[2]), "holm"), cmp)
  data.frame(window = w, h = hh, best = lab[best], RMSE = a$RMSE[1], runner_up = lab[a$model[2]],
             gap_pct = 100 * (a$RMSE[2] / a$RMSE[1] - 1), rel_SNaive = a$RMSE[1] / a$RMSE[a$model == "SNaive"],
             p_vs_SNaive = if (best == "SNaive") NA else dm_hln(E[[best]], E$SNaive, hh)[2],
             p_vs_SNaive_holm = if (best == "SNaive") NA else unname(holm[best]),
             p_vs_ARIMA = if (best == "ARIMA") NA else dm_hln(E[[best]], E$ARIMA, hh)[2], row.names = NULL) }))))
print(r3(tab10r), right = FALSE); save_tab(tab10r, "Table10_best_vs_benchmarks.csv")

# ---- Prices, first version (all pairs with Clark-West); superseded by Table 14 in Section 5.8 ----
banner("Clark-West tests for all price pairs (first version; superseded by Table 14 in Section 5.8)")
price_pairs <- list(c("ARIMAX_feasible","ARIMA","ARIMAX / seasonal ARIMA"),
                    c("ARIMAX_pricefc","ARIMA","ARIMAX, forecast prices / seasonal ARIMA"),
                    c("GBM","GBM_noCom","GBM / GBM, no prices"), c("RF","RF_noCom","RF / RF, no prices"),
                    c("H1_ARIMA_GBM","H1_noCom","H1 / H1, no prices"), c("H3_STL_GBM","H3_STL_GBM_noCom","H3 / H3, no prices"))
tab12 <- pair_table(price_pairs, "cw"); print(wide_cells(tab12), right = FALSE)
cat("Ratios below one:", sum(tab12$ratio < 1), "of", nrow(tab12), "| Clark-West p < 0.05:", sum(tab12$p < 0.05), "\n")
print(tab12 %>% group_by(comparison) %>% summarise(below1 = sum(ratio < 1), cw05 = sum(p < 0.05),
      min_holm_within = min(p.adjust(p, "holm")), .groups = "drop"), width = Inf)
cat("Smallest Holm-adjusted p across all", nrow(tab12), "tests:", round(min(p.adjust(tab12$p, "holm")), 3), "\n")
cat("ARIMAX with forecast price paths relative to no-change paths (RMSE ratio, DM p):\n")
print(wide_cells(pair_table(list(c("ARIMAX_pricefc","ARIMAX_feasible","forecast / no-change prices")), "dm")), right = FALSE)
save_tab(tab12, "Prices_ClarkWest_all_pairs_superseded.csv")

# ---- Table 17: horizon-specific versus one-step correctors ----
banner("Table 17. Horizon-specific / one-step residual correctors: RMSE ratio (two-sided DM p)")
tab14 <- pair_table(list(c("H1h","H1_ARIMA_GBM","H1: ARIMA-GBM"), c("H2h","H2_ARIMA_RF","H2: ARIMA-RF"),
                         c("H3h","H3_STL_GBM","H3: STL-GBM"), c("H4h","H4_ETS_GBM","H4: ETS-GBM")), "dm")
tab14 <- tab14[tab14$h > 1, ]                         # both designs coincide at h = 1
print(wide_cells(tab14), right = FALSE)
cat("Holm across all", nrow(tab14), "tests, smallest p:", round(min(p.adjust(tab14$p, "holm")), 3), "\n")
cat("Holm within each column (four hybrids), significant at 5%:\n")
print(as.data.frame(tab14 %>% group_by(window, h) %>% mutate(p_holm = p.adjust(p, "holm")) %>%
      ungroup() %>% filter(p_holm < 0.05) %>% select(comparison, window, h, ratio, p, p_holm)))
tab14b <- pair_table(list(c("H1h","ARIMA","H1h / seasonal ARIMA"), c("H2h","ARIMA","H2h / seasonal ARIMA"),
                          c("H3h","STL_ARIMA","H3h / STL-ARIMA"), c("H4h","ETS","H4h / ETS")), "dm")
cat("Horizon-specific corrector relative to its own linear stage (lower panel):\n")
print(wide_cells(tab14b[tab14b$h > 1, ]), right = FALSE)
save_tab(rbind(tab14, tab14b), "Table17_horizon_specific.csv")

# ---- Table 18: empirical predictive intervals ---------------
banner("Table 18. Empirical predictive intervals, window A: all past errors vs the last 24 errors")
pi_emp_rev <- function(mdl, hh, last = Inf, minn = 12) {
  s <- REX %>% filter(model == mdl, h == hh) %>% arrange(origin)
  q <- t(sapply(seq_len(nrow(s)), function(i) { past <- s$err[s$target <= s$origin[i]]
    if (is.finite(last)) past <- tail(past, last)
    c(length(past), if (length(past) >= minn) s$fc[i] + quantile(past, c(.10, .90, .025, .975)) else rep(NA, 4)) }))
  colnames(q) <- c("n_err","lo80","hi80","lo95","hi95"); cbind(s, q) }
kupiec <- function(hits, p) { n <- length(hits); x <- sum(hits); ph <- x / n
  ll1 <- if (x == 0 || x == n) 0 else (n - x) * log(1 - ph) + x * log(ph)
  LR <- -2 * ((n - x) * log(1 - p) + x * log(p) - ll1); pchisq(LR, 1, lower.tail = FALSE) }
tab16 <- do.call(rbind, lapply(c("H1_ARIMA_GBM","ARIMA","SNaive","GBM","H3_STL_GBM"), function(m)
  do.call(rbind, lapply(c(1, 3, 6), function(hh) do.call(rbind, lapply(c(Inf, 24), function(last) {
    s <- pi_emp_rev(m, hh, last) %>% filter(date >= WIN$A[1], date <= WIN$A[2], !is.na(lo80))
    h80 <- s$actual >= s$lo80 & s$actual <= s$hi80; h95 <- s$actual >= s$lo95 & s$actual <= s$hi95
    data.frame(model = lab[m], h = hh, errors = if (is.finite(last)) "last 24" else "all past", n = nrow(s),
               n_err = paste(min(s$n_err), max(s$n_err), sep = "-"), cov80 = mean(h80), kupiec80_p = kupiec(h80, 0.8),
               width80 = mean(s$hi80 - s$lo80), cov95 = mean(h95), kupiec95_p = kupiec(h95, 0.95),
               width95 = mean(s$hi95 - s$lo95), row.names = NULL) }))))))
print(r3(tab16), row.names = FALSE); save_tab(tab16, "Table18_intervals.csv")
cat("Standard deviation of the one-step errors before and during window A:\n")
print(as.data.frame(r3(REX %>% filter(model %in% c("H1_ARIMA_GBM","ARIMA","SNaive"), h == 1) %>%
  mutate(period = ifelse(date < WIN$A[1], "2022-2023", "2024-2025")) %>%
  group_by(model, period) %>% summarise(sd_err = sd(err), n = n(), .groups = "drop"))))

# ---- Table A5: hyperparameter sensitivity -----------
banner("Table A5. Machine-learning hyperparameters: RMSE ratio alternative / baseline (two-sided DM p)")
tabA5 <- pair_table(list(c("RF_alt","RF","RF: mtry = d/3, node size 5"), c("GBM_d1","GBM","GBM: depth 1, up to 1,000 trees"),
                         c("GBM_fixed","GBM","GBM: 100 trees, no inner validation"), c("H1_d1","H1_ARIMA_GBM","H1: corrector depth 1"),
                         c("H1_fixed","H1_ARIMA_GBM","H1: corrector 50 trees, no inner validation")), "dm")
print(wide_cells(tabA5), right = FALSE)
cat("Holm across all", nrow(tabA5), "tests, smallest p:", round(min(p.adjust(tabA5$p, "holm")), 3), "\n")
cat("Training rows for the learners (first origin / December 2023 / last origin):\n")
print(as.data.frame(RREV %>% filter(model == "RF_alt") %>% group_by(h) %>%
  summarise(first = n_train[origin == min(origin)], dec2023 = n_train[origin == which(d$date == as.Date("2023-12-01"))][1],
            last = n_train[origin == max(origin)], .groups = "drop")))
save_tab(tabA5, "TableA5_hyperparameters.csv")


# =============================================================================
#  SECTION 5.7. FURTHER ROBUSTNESS ANALYSES: Tables A6-A14
# =============================================================================
# Additional variants, expanding window, six lags (none enters the 15 feasible methods):
#   price shocks: GBM_shock, H1_shock (AR(1) innovations of the log price changes
#     instead of price levels), GBM_nopi (Hamilton net price increases), ARIMAX_shock
#     (innovations as regressors, expected future innovations = 0)
#   import costs: GBM_brent, H1_brent, ARIMAX_brent (Brent crude price, USD)
#   small sample: GBM_pars, H1_pars (9 features), GBM_nested, H1_nested
#     (nested time-series validation of depth and number of trees), GBM_sel (top-8 features
#     by in-window importance), GBM_wheat/_fuel/_elec/_gas (QIM plus a single price)
#   direct hybrids: DAR (direct seasonal AR for each h), HD (DAR + GBM corrector
#     of its in-sample h-step errors)
#   price paths: ARIMAX_drift, ARDL_pricefc, ARDL_drift (drift = average monthly
#     change over the last 12 months; pricefc = auto.arima forecast of each price)
#   real time: *_rt methods use QIM only up to o - 1 (published about 45 days
#     after the reference month) and prices up to o
#   intervals: residual-bootstrap intervals of seasonal ARIMA (stored in PI2)
gbm_cv_pred <- function(Xtr, ytr, Xte, depths = 1:3, ntmax = 500, nfold = 3, vlen = 6, seed = 42) {
  # nested time-series validation: 'nfold' expanding folds, each validated on the next 'vlen'
  # months at the end of the estimation window; depth and number of trees minimise the
  # average validation MSE; the model is then refitted to the whole window
  n <- nrow(Xtr); df <- data.frame(yv = ytr, Xtr); dte <- data.frame(Xte); colnames(dte) <- colnames(df)[-1]
  best <- c(depth = NA, M = NA, mse = Inf)
  for (dp in depths) {
    curve <- 0
    for (f in seq_len(nfold)) { end_tr <- n - vlen * (nfold - f + 1); va <- (end_tr + 1):(end_tr + vlen)
      set.seed(seed)
      m <- gbm(yv ~ ., data = df[1:end_tr, ], distribution = "gaussian", n.trees = ntmax, interaction.depth = dp,
               shrinkage = 0.05, bag.fraction = 0.8, n.minobsinnode = 5, verbose = FALSE)
      pv <- predict(m, df[va, ], n.trees = 1:ntmax); curve <- curve + colMeans((pv - ytr[va])^2) / nfold }
    curve <- unname(curve); Mb <- max(10, which.min(curve)); if (curve[Mb] < best[["mse"]]) best <- c(depth = dp, M = Mb, mse = curve[Mb]) }
  set.seed(seed)
  m <- gbm(yv ~ ., data = df, distribution = "gaussian", n.trees = best[["M"]], interaction.depth = best[["depth"]],
           shrinkage = 0.05, bag.fraction = 0.8, n.minobsinnode = 5, verbose = FALSE)
  as.numeric(predict(m, dte, n.trees = best[["M"]]))
}
gbm_sel_pred <- function(Xtr, ytr, Xte, keep = 8, seed = 42) {
  # in-window variable selection: rank features by relative influence in a GBM fitted to the
  # window, keep the 'keep' most important ones and refit with the usual validation
  set.seed(seed); df <- data.frame(yv = ytr, Xtr)
  m <- gbm(yv ~ ., data = df, distribution = "gaussian", n.trees = 300, interaction.depth = 3,
           shrinkage = 0.05, bag.fraction = 0.8, n.minobsinnode = 5, verbose = FALSE)
  ri <- summary(m, plotit = FALSE); top <- as.character(ri$var[seq_len(min(keep, nrow(ri)))])
  gbm_fit_pred(Xtr[, top, drop = FALSE], ytr, Xte[, top, drop = FALSE], seed = seed)
}
# generic feature builder: QIM lags dated t - h - qoff, lags of the columns of Z dated t - h
feat_gen <- function(t, h, Z = NULL, LAGS = 6, qoff = 0, pars = FALSE) {
  if (pars) {                                    # parsimonious set (9 features)
    v <- c(y[t - h - qoff], y[t - h - qoff - 1], y[t - 12 - (if (t - 12 >= t - h - qoff) 12 else 0)],
           if (!is.null(Z)) Z[t - h, ], mon[t], t)
    names(v) <- c("QIM_l0", "QIM_l1", "QIM_s12", if (!is.null(Z)) colnames(Z), "month", "trend"); return(v) }
  v <- sapply(0:(LAGS - 1), function(j) y[t - h - qoff - j]); nm <- paste0("QIM_l", 0:(LAGS - 1))
  if (!is.null(Z)) for (k in colnames(Z)) { v <- c(v, sapply(0:(LAGS - 1), function(j) Z[t - h - j, k]))
                                            nm <- c(nm, paste0(k, "_l", 0:(LAGS - 1))) }
  v <- c(v, mon[t], t); names(v) <- c(nm, "month", "trend"); v
}
fmat_gen <- function(ts_, h, ...) do.call(rbind, lapply(ts_, feat_gen, h = h, ...))
one_origin_rev2 <- function(o, LAGS = 6) {
  tr <- 1:o; st <- c(as.numeric(format(d$date[1], "%Y")), mon[1])
  yt <- ts(y[tr], start = st, frequency = 12); hmax <- min(12, N - o)
  aa <- function(z, ...) auto.arima(z, ic = "aic", stepwise = TRUE, approximation = FALSE, ...)
  fA <- aa(yt, seasonal = TRUE); pA <- as.numeric(forecast(fA, h = hmax)$mean)
  eA <- as.numeric(y[tr] - fitted(fA))
  set.seed(1000 + o)                                       # residual-bootstrap intervals of seasonal ARIMA
  fb <- forecast(fA, h = hmax, level = c(80, 95), bootstrap = TRUE, npaths = 2000)
  # price shocks estimated on the window only
  INN <- sapply(PR, function(k) price_innov(X[, k], rows = tr)); INN[is.na(INN)] <- 0
  NOP <- sapply(PR, function(k) nopi(X[, k]))
  # ARIMAX variants (same regressors as the baseline ARIMAX unless stated)
  fX <- aa(yt, xreg = X[tr, ], seasonal = TRUE)
  drift <- (X[o, ] - X[o - 12, ]) / 12
  Xdr <- t(sapply(1:hmax, function(s) X[o, ] + s * drift)); colnames(Xdr) <- PR
  Xpf <- matrix(sapply(PR, function(k) as.numeric(forecast(auto.arima(ts(X[tr, k], start = st, frequency = 12),
                seasonal = FALSE, ic = "aic"), h = hmax)$mean)), nrow = hmax, dimnames = list(NULL, PR))
  pXnc <- as.numeric(forecast(fX, xreg = X[rep(o, hmax), , drop = FALSE])$mean)
  pXdr <- as.numeric(forecast(fX, xreg = Xdr)$mean)
  fXs <- aa(yt, xreg = INN[tr, ], seasonal = TRUE)
  pXs <- as.numeric(forecast(fXs, xreg = matrix(0, hmax, length(PR), dimnames = list(NULL, PR)))$mean)
  fXb <- aa(yt, xreg = cbind(BRENT = BRENT[tr]), seasonal = TRUE)
  pXb <- as.numeric(forecast(fXb, xreg = cbind(BRENT = rep(BRENT[o], hmax)))$mean)
  # ARDL selected once, forecast recursively under three price paths
  a <- auto_ardl(QIM ~ WHEAT + FUEL + ENERGY + GAS, data = d[tr, ], max_order = 6)$best_model; b <- coef(a)
  ardl_path <- function(P) { yy <- y; XX <- X; XX[(o + 1):(o + hmax), ] <- P; fc <- numeric(hmax)
    for (s in 1:hmax) { t <- o + s; v <- b["(Intercept)"]
      for (nm in names(b)[-1]) { lag <- if (grepl("^L\\(", nm)) as.numeric(sub(".*, *([0-9]+)\\)$", "\\1", nm)) else 0
        var <- sub("^L\\(([A-Z]+),.*$", "\\1", nm); v <- v + b[nm] * (if (var == "QIM") yy[t - lag] else XX[t - lag, var]) }
      yy[t] <- v; fc[s] <- v }; fc }
  pDpf <- ardl_path(Xpf); pDdr <- ardl_path(Xdr)
  # pseudo-real-time: QIM observed up to o - 1, prices up to o
  trr <- 1:(o - 1); ytr <- ts(y[trr], start = st, frequency = 12)
  fAr <- aa(ytr, seasonal = TRUE); pAr <- as.numeric(forecast(fAr, h = hmax + 1)$mean); eAr <- as.numeric(y[trr] - fitted(fAr))
  pEr <- as.numeric(forecast(ets(ytr), h = hmax + 1)$mean)
  fXr <- aa(ytr, xreg = X[trr, ], seasonal = TRUE)
  pXr <- as.numeric(forecast(fXr, xreg = X[rep(o, hmax + 1), , drop = FALSE])$mean)   # month o known, then no change
  res <- list(); pis <- list()
  for (h in H[H <= hmax]) {
    tt <- o + h
    t_tr <- tr[tr - h - (LAGS - 1) >= 1]; idx <- match(t_tr, tr)
    t_s  <- tr[tr - h - (LAGS - 1) >= 3]; idx_s <- match(t_s, tr)            # shock lags need two earlier months
    t_p  <- tr[tr - h - 1 >= 1 & tr - 12 - (if (h == 12) 12 else 0) >= 1]; idx_p <- match(t_p, tr)
    F1 <- fmat(t_tr, h, LAGS); f1 <- rbind(feat(tt, h, LAGS))
    F0 <- fmat(t_tr, h, LAGS, com = FALSE); f0_ <- rbind(feat(tt, h, LAGS, com = FALSE))
    FS <- fmat_gen(t_s, h, Z = INN); fs <- rbind(feat_gen(tt, h, Z = INN))
    FN <- fmat_gen(t_tr, h, Z = NOP); fn <- rbind(feat_gen(tt, h, Z = NOP))
    BZ <- cbind(BRENT = BRENT)
    FB <- fmat_gen(t_tr, h, Z = BZ); fb_ <- rbind(feat_gen(tt, h, Z = BZ))
    FP <- fmat_gen(t_p, h, Z = X, pars = TRUE); fp <- rbind(feat_gen(tt, h, Z = X, pars = TRUE))
    single <- sapply(PR, function(k) { Zk <- X[, k, drop = FALSE]
      gbm_fit_pred(fmat_gen(t_tr, h, Z = Zk), y[t_tr], rbind(feat_gen(tt, h, Z = Zk))) })
    # direct seasonal AR for horizon h, and its hybrid with a GBM corrector (direct h-step residuals)
    dar_df <- function(ts_) data.frame(setNames(as.data.frame(matrix(sapply(0:(LAGS - 1), function(j) y[ts_ - h - j]),
                                       nrow = length(ts_))), paste0("QIM_l", 0:(LAGS - 1))),
                                       month = factor(mon[ts_], levels = 1:12), trend = ts_)
    Dtr <- dar_df(t_tr); fD <- lm(y[t_tr] ~ ., data = Dtr)
    pDAR <- as.numeric(predict(fD, newdata = dar_df(tt))); eD <- as.numeric(residuals(fD))
    # pseudo-real-time ML: training targets observed by o - 1, QIM lags shifted by one month
    t_r <- trr[trr - h - 1 - (LAGS - 1) >= 1]; idx_r <- match(t_r, trr)
    FR <- fmat_gen(t_r, h, Z = X, qoff = 1); fr <- rbind(feat_gen(tt, h, Z = X, qoff = 1))
    r <- c(GBM_shock = gbm_fit_pred(FS, y[t_s], fs), H1_shock = pA[h] + gbm_fit_pred(FS, eA[idx_s], fs, depth = 2, ntmax = 300),
           GBM_nopi = gbm_fit_pred(FN, y[t_tr], fn), ARIMAX_shock = pXs[h],
           GBM_brent = gbm_fit_pred(FB, y[t_tr], fb_), H1_brent = pA[h] + gbm_fit_pred(FB, eA[idx], fb_, depth = 2, ntmax = 300),
           ARIMAX_brent = pXb[h],
           GBM_pars = gbm_fit_pred(FP, y[t_p], fp), H1_pars = pA[h] + gbm_fit_pred(FP, eA[idx_p], fp, depth = 2, ntmax = 300),
           GBM_nested = gbm_cv_pred(F1, y[t_tr], f1),
           H1_nested = pA[h] + gbm_cv_pred(F1, eA[idx], f1, depths = 1:2, ntmax = 300),
           GBM_sel = gbm_sel_pred(F1, y[t_tr], f1),
           GBM_wheat = single[["WHEAT"]], GBM_fuel = single[["FUEL"]], GBM_elec = single[["ENERGY"]], GBM_gas = single[["GAS"]],
           DAR = pDAR, HD = pDAR + gbm_fit_pred(F1, eD, f1, depth = 2, ntmax = 300),
           ARIMAX_nochange_chk = pXnc[h], ARIMAX_drift = pXdr[h], ARDL_pricefc = pDpf[h], ARDL_drift = pDdr[h],
           SNaive_rt = if (h <= 11) y[o + h - 12] else y[o - 12], ARIMA_rt = pAr[h + 1], ETS_rt = pEr[h + 1],
           ARIMAX_rt = pXr[h + 1], GBM_rt = gbm_fit_pred(FR, y[t_r], fr),
           H1_rt = pAr[h + 1] + gbm_fit_pred(FR, eAr[idx_r], fr, depth = 2, ntmax = 300))
    res[[length(res) + 1]] <- data.frame(origin = o, h = h, target = tt, date = d$date[tt], actual = y[tt],
                                         model = names(r), fc = as.numeric(r))
    pis[[length(pis) + 1]] <- data.frame(origin = o, h = h, target = tt, date = d$date[tt], actual = y[tt],
                                         model = "ARIMA", fc = pA[h], lo80 = fb$lower[h, 1], hi80 = fb$upper[h, 1],
                                         lo95 = fb$lower[h, 2], hi95 = fb$upper[h, 2])
  }
  list(fc = do.call(rbind, res), pi = do.call(rbind, pis))
}
run_rev2 <- function() {
  cf <- file.path(OUT_DIR, "cache", "rolling_expanding_revision2.rds")
  if (file.exists(cf)) { cat("  using cached forecasts: revision-2 variants\n"); return(readRDS(cf)) }
  cat(sprintf("  computing revision-2 variants (%d origins, %d core(s)) ...\n", length(origins), N_CORES)); t0 <- Sys.time()
  if (N_CORES > 1 && .Platform$OS.type == "windows") {
    cl <- makeCluster(N_CORES); on.exit(stopCluster(cl))
    clusterEvalQ(cl, suppressPackageStartupMessages({ library(forecast); library(gbm); library(randomForest); library(ARDL) }))
    clusterExport(cl, c("one_origin_rev2","feat","fmat","feat_gen","fmat_gen","gbm_fit_pred","gbm_cv_pred","gbm_sel_pred",
                        "price_innov","nopi","d","y","X","BRENT","mon","N","PR","H"), envir = environment())
    R <- parLapply(cl, origins, one_origin_rev2)
  } else if (N_CORES > 1) {
    R <- mclapply(origins, one_origin_rev2, mc.cores = N_CORES, mc.preschedule = FALSE)
  } else R <- lapply(origins, one_origin_rev2)
  out <- list(fc = do.call(rbind, lapply(R, `[[`, "fc")), pi = do.call(rbind, lapply(R, `[[`, "pi")))
  out$fc$tag <- "expanding"; out$fc$err <- out$fc$actual - out$fc$fc
  saveRDS(out, cf); cat("  done in", format(round(Sys.time() - t0, 1)), "\n"); out
}
banner("Section 5.7. Revision-2 analyses: rolling-origin forecasts of the additional variants (expanding, 6 lags)")
RV2  <- run_rev2()
REX2 <- rbind(REX, RV2$fc[, c("origin","h","target","date","actual","model","fc","tag","err")])
chk <- merge(RV2$fc[RV2$fc$model == "ARIMAX_nochange_chk", c("origin","h","fc")],
             REX[REX$model == "ARIMAX_feasible", c("origin","h","fc")], by = c("origin","h"))
cat("Check: ARIMAX with no-change prices reproduced (max abs. difference):", signif(max(abs(chk$fc.x - chk$fc.y)), 3), "\n")

accR2 <- do.call(rbind, lapply(names(WIN), function(w) inwin(REX2, w) %>% group_by(model, h) %>%
  summarise(window = w, n = n(), RMSE = sqrt(mean(err^2)), .groups = "drop")))
rmse_rev2 <- as.data.frame(accR2 %>% pivot_wider(id_cols = model, names_from = c(window, h), values_from = RMSE))
rmse_rev2 <- rmse_rev2[, c("model","A_1","A_3","A_6","A_12","B_1","B_3","B_6","B_12")]
save_tab(rmse_rev2, "Rev2_RMSE_all_methods.csv")
E_of <- function(w, hh, models) { s <- inwin(REX2, w) %>% filter(h == hh, model %in% models)   # now on REX2
  list(e = s %>% select(date, model, err) %>% pivot_wider(names_from = model, values_from = err) %>% arrange(date),
       f = s %>% select(date, model, fc) %>% pivot_wider(names_from = model, values_from = fc) %>% arrange(date)) }
rmse_cells <- function(models) { x <- rmse_rev2[match(models, rmse_rev2$model), ]
  x$model <- ifelse(is.na(lab[x$model]), x$model, lab[x$model]); r3(x) }

# ---- Price shocks and import costs, first version; superseded by Table A8 in Section 5.8 ----
banner("Price shocks and import costs, Clark-West for all pairs (first version; superseded by Table A8 in Section 5.8)")
tabA8 <- pair_table(list(
  c("GBM","GBM_noCom","GBM: price levels (baseline)"), c("GBM_shock","GBM_noCom","GBM: price innovations"),
  c("GBM_nopi","GBM_noCom","GBM: net price increases"), c("GBM_brent","GBM_noCom","GBM: Brent crude price"),
  c("H1_ARIMA_GBM","H1_noCom","H1: price levels (baseline)"), c("H1_shock","H1_noCom","H1: price innovations"),
  c("H1_brent","H1_noCom","H1: Brent crude price"),
  c("ARIMAX_feasible","ARIMA","ARIMAX: price levels (baseline)"), c("ARIMAX_shock","ARIMA","ARIMAX: price innovations"),
  c("ARIMAX_brent","ARIMA","ARIMAX: Brent crude price")), "cw")
print(wide_cells(tabA8), right = FALSE)
cat("Ratios below one by price measure:\n"); print(as.data.frame(tabA8 %>% group_by(comparison) %>%
  summarise(below1 = sum(ratio < 1), cw05 = sum(p < 0.05), holm_within_min = min(p.adjust(p, "holm")), .groups = "drop")))
cat("Smallest Holm-adjusted p across the", sum(!grepl("baseline", tabA8$comparison)), "new tests:",
    round(min(p.adjust(tabA8$p[!grepl("baseline", tabA8$comparison)], "holm")), 3), "\n")
save_tab(tabA8, "PriceShocks_ClarkWest_all_pairs_superseded.csv")

# ---- Table A9: small-sample remedies for the machine-learning stage ----
banner("Table A9. Parsimonious features, nested validation and variable selection: RMSE ratio to the baseline (DM p)")
tabA9 <- pair_table(list(
  c("GBM_pars","GBM","GBM: 9 parsimonious features"), c("GBM_nested","GBM","GBM: nested validation of depth and trees"),
  c("GBM_sel","GBM","GBM: top-8 features selected in the window"),
  c("GBM_wheat","GBM","GBM: QIM + wheat only"), c("GBM_fuel","GBM","GBM: QIM + fuel only"),
  c("GBM_elec","GBM","GBM: QIM + electricity only"), c("GBM_gas","GBM","GBM: QIM + gas only"),
  c("H1_pars","H1_ARIMA_GBM","H1: 9 parsimonious features"), c("H1_nested","H1_ARIMA_GBM","H1: nested validation")), "dm")
print(wide_cells(tabA9), right = FALSE)
cat("Holm across all", nrow(tabA9), "tests, smallest p:", round(min(p.adjust(tabA9$p, "holm")), 3), "\n")
cat("RMSE levels:\n"); print(rmse_cells(c("GBM","GBM_noCom","GBM_pars","GBM_nested","GBM_sel","GBM_wheat","GBM_fuel","GBM_elec","GBM_gas","H1_ARIMA_GBM","H1_pars","H1_nested")), row.names = FALSE)
save_tab(tabA9, "TableA9_small_sample.csv")

# ---- Table A10: recursive versus direct hybrids ----
banner("Table A10. Recursive and direct linear stages and hybrids: RMSE")
tabA10 <- rmse_cells(c("ARIMA","H1_ARIMA_GBM","H1h","DAR","HD")); print(tabA10, row.names = FALSE)
tabA10b <- pair_table(list(c("HD","H1_ARIMA_GBM","HD / H1"), c("HD","H1h","HD / H1h"), c("DAR","ARIMA","DAR / seasonal ARIMA"),
                           c("HD","DAR","HD / DAR")), "dm")
print(wide_cells(tabA10b), right = FALSE)
save_tab(rbind(transform(tabA10b)), "TableA10_direct_hybrids.csv")

# ---- Table A11: future price paths ----
banner("Table A11. Future price paths for ARIMAX and ARDL: RMSE")
print(rmse_cells(c("ARIMAX_feasible","ARIMAX_drift","ARIMAX_pricefc","ARDL_feasible","ARDL_drift","ARDL_pricefc")), row.names = FALSE)
spread <- REX2 %>% filter(model %in% c("ARIMAX_feasible","ARIMAX_drift","ARIMAX_pricefc"), date >= WIN$B[1]) %>%
  group_by(origin, h) %>% summarise(range = max(fc) - min(fc), .groups = "drop") %>% group_by(h) %>%
  summarise(mean_range = mean(range), max_range = max(range), .groups = "drop")
cat("Price-path uncertainty: range of the three ARIMAX forecasts (index points), window B targets:\n"); print(as.data.frame(r3(spread)))
save_tab(rmse_cells(c("ARIMAX_feasible","ARIMAX_drift","ARIMAX_pricefc","ARDL_feasible","ARDL_drift","ARDL_pricefc")), "TableA11_price_paths.csv")

# ---- Table A12: pseudo-real-time evaluation ----
banner("Table A12. Pseudo-real-time forecasts (QIM known to o - 1, prices to o): RMSE and ratio to the final-data design")
rt <- c(SNaive = "SNaive_rt", ARIMA = "ARIMA_rt", ETS = "ETS_rt", ARIMAX_feasible = "ARIMAX_rt", GBM = "GBM_rt", H1_ARIMA_GBM = "H1_rt")
tabA12 <- do.call(rbind, lapply(names(rt), function(m) { a <- rmse_rev2[rmse_rev2$model == m, -1]; b <- rmse_rev2[rmse_rev2$model == rt[m], -1]
  rbind(data.frame(method = lab[m], design = "final data", a), data.frame(method = lab[m], design = "pseudo-real-time", b),
        data.frame(method = lab[m], design = "ratio", b / a)) }))
print(r3(tabA12), row.names = FALSE)
cat("Lowest RMSE among the six pseudo-real-time methods:\n")
print(sapply(names(rmse_rev2)[-1], function(cn) { x <- rmse_rev2[rmse_rev2$model %in% rt, ]; x$model[which.min(x[[cn]])] }))
save_tab(tabA12, "TableA12_pseudo_real_time.csv")

# ---- Table A13: interval methods ----
banner("Table A13. Predictive intervals by construction method, window A: coverage, width and interval score")
int_score <- function(lo, hi, yv, a) (hi - lo) + (2 / a) * pmax(lo - yv, 0) + (2 / a) * pmax(yv - hi, 0)
interval_methods <- function(mdl, hh, minn = 12, gamma = 0.05) {
  s <- REX2 %>% filter(model == mdl, h == hh) %>% arrange(origin); n <- nrow(s)
  out <- list(); alpha_t <- c(`80` = 0.2, `95` = 0.05); used <- rep(FALSE, n); hitlog <- list(`80` = rep(NA, n), `95` = rep(NA, n))
  for (i in seq_len(n)) {
    obs <- which(s$target <= s$origin[i]); past <- s$err[obs]
    # adaptive conformal inference: update alpha with every interval whose target has just been observed
    for (j in obs[!used[obs]]) { for (lv in c("80","95")) if (!is.na(hitlog[[lv]][j]))
      alpha_t[lv] <- alpha_t[lv] + gamma * ((if (lv == "80") 0.2 else 0.05) - (1 - hitlog[[lv]][j])); used[j] <- TRUE }
    if (length(past) < minn) next
    row <- function(meth, l80, h80, l95, h95) data.frame(origin = s$origin[i], date = s$date[i], actual = s$actual[i], method = meth,
                                                         lo80 = l80, hi80 = h80, lo95 = l95, hi95 = h95)
    q <- quantile(past, c(.1, .9, .025, .975))
    out[[length(out) + 1]] <- row("Empirical, all past errors", s$fc[i] + q[1], s$fc[i] + q[2], s$fc[i] + q[3], s$fc[i] + q[4])
    q24 <- quantile(tail(past, 24), c(.1, .9, .025, .975))
    out[[length(out) + 1]] <- row("Empirical, last 24 errors", s$fc[i] + q24[1], s$fc[i] + q24[2], s$fc[i] + q24[3], s$fc[i] + q24[4])
    sig <- sqrt(mean(tail(past, 12)^2))
    out[[length(out) + 1]] <- row("Gaussian, RMSE of last 12 errors", s$fc[i] - 1.2816 * sig, s$fc[i] + 1.2816 * sig, s$fc[i] - 1.96 * sig, s$fc[i] + 1.96 * sig)
    a80 <- min(max(alpha_t["80"], 0.005), 0.995); a95 <- min(max(alpha_t["95"], 0.005), 0.995)
    qa <- quantile(past, c(a80 / 2, 1 - a80 / 2, a95 / 2, 1 - a95 / 2))
    out[[length(out) + 1]] <- row("Adaptive conformal (ACI)", s$fc[i] + qa[1], s$fc[i] + qa[2], s$fc[i] + qa[3], s$fc[i] + qa[4])
    hitlog[["80"]][i] <- as.numeric(s$actual[i] >= s$fc[i] + qa[1] & s$actual[i] <= s$fc[i] + qa[2])
    hitlog[["95"]][i] <- as.numeric(s$actual[i] >= s$fc[i] + qa[3] & s$actual[i] <= s$fc[i] + qa[4])
    pf <- data.frame(act = s$actual[obs], fc = s$fc[obs])
    qr <- sapply(c(.1, .9, .025, .975), function(tau) as.numeric(predict(quantreg::rq(act ~ fc, tau = tau, data = pf),
                                                                       newdata = data.frame(fc = s$fc[i]))))
    out[[length(out) + 1]] <- row("Quantile regression", min(qr[1:2]), max(qr[1:2]), min(qr[3:4]), max(qr[3:4]))
  }
  # residual bootstrap of seasonal ARIMA (model-based; needs no error history)
  bt <- RV2$pi %>% filter(h == hh) %>% arrange(origin)
  shift <- if (mdl == "ARIMA") 0 else s$fc[match(bt$origin, s$origin)] - bt$fc
  if (mdl %in% c("ARIMA", "H1_ARIMA_GBM")) out[[length(out) + 1]] <- data.frame(origin = bt$origin, date = bt$date, actual = bt$actual,
      method = "Residual bootstrap (seasonal ARIMA)", lo80 = bt$lo80 + shift, hi80 = bt$hi80 + shift, lo95 = bt$lo95 + shift, hi95 = bt$hi95 + shift)
  do.call(rbind, out)
}
tabA13 <- do.call(rbind, lapply(c("H1_ARIMA_GBM","ARIMA","SNaive"), function(m) do.call(rbind, lapply(H, function(hh) {
  x <- interval_methods(m, hh) %>% filter(date >= WIN$A[1], date <= WIN$A[2])
  x %>% group_by(method) %>% summarise(n = n(), cov80 = mean(actual >= lo80 & actual <= hi80),
    kupiec80_p = kupiec(actual >= lo80 & actual <= hi80, 0.8), width80 = mean(hi80 - lo80),
    IS80 = mean(int_score(lo80, hi80, actual, 0.2)), cov95 = mean(actual >= lo95 & actual <= hi95),
    width95 = mean(hi95 - lo95), IS95 = mean(int_score(lo95, hi95, actual, 0.05)), .groups = "drop") %>%
    mutate(model = lab[m], h = hh, .before = 1) }))))
print(as.data.frame(r3(tabA13)), row.names = FALSE)
save_tab(tabA13, "TableA13_interval_methods.csv")

# ---- Table A14: economic loss functions ----
banner("Table A14. Asymmetric (lin-lin) loss and frequency of large errors, window A, relative to the seasonal naive")
linlin <- function(e, a_under, b_over) mean(ifelse(e > 0, a_under * e, -b_over * e))   # e = actual - forecast
tabA14 <- do.call(rbind, lapply(c("A","B"), function(w) do.call(rbind, lapply(c(1, 3, 6, 12), function(hh) {
  s <- inwin(REX2, w) %>% filter(h == hh, model %in% feasible)
  x <- s %>% group_by(model) %>% summarise(L_under2 = linlin(err, 2, 1), L_over2 = linlin(err, 1, 2),
         large5 = mean(abs(err) / actual > 0.05), RMSE = sqrt(mean(err^2)), .groups = "drop")
  sn <- x[x$model == "SNaive", ]
  x %>% mutate(window = w, h = hh, rel_under2 = L_under2 / sn$L_under2, rel_over2 = L_over2 / sn$L_over2,
               rank_RMSE = rank(RMSE), rank_under2 = rank(L_under2), rank_over2 = rank(L_over2)) }))))
tabA14$method <- lab[tabA14$model]
best_by <- tabA14 %>% group_by(window, h) %>% summarise(best_RMSE = lab[model[which.min(RMSE)]],
  best_under2 = lab[model[which.min(L_under2)]], best_over2 = lab[model[which.min(L_over2)]],
  best_large5 = lab[model[which.min(large5)]], min_large5 = min(large5), H1_large5 = large5[model == "H1_ARIMA_GBM"],
  SN_large5 = large5[model == "SNaive"], .groups = "drop")
print(as.data.frame(r3(best_by)), right = FALSE)
cat("Rank correlation between RMSE ranks and lin-lin ranks (mean over columns):",
    round(mean(tabA14 %>% group_by(window, h) %>% summarise(c1 = cor(rank_RMSE, rank_under2), .groups = "drop") %>% pull(c1)), 2),
    round(mean(tabA14 %>% group_by(window, h) %>% summarise(c2 = cor(rank_RMSE, rank_over2), .groups = "drop") %>% pull(c2)), 2), "\n")
print(as.data.frame(r3(tabA14 %>% filter(window == "A", h %in% c(1, 3)) %>%
  select(h, method, rel_under2, rel_over2, large5, rank_RMSE, rank_under2, rank_over2))), row.names = FALSE)
save_tab(tabA14, "TableA14_economic_loss.csv")


# =============================================================================
#  SECTION 5.8. ADDITIONAL BENCHMARKS AND ORDER-CONSISTENT PRICE COMPARISONS
# =============================================================================
#   Benchmarks: Theta, TBATS, level-adjusted seasonal naive, seasonal ARIMA with d = D = 0,
#     seasonal ARDL (monthly and COVID dummies, order selected with the dummies, no-change
#     prices), and the equal-weight combination of seasonal naive, seasonal ARIMA and ETS
#     (computed after the run).
#   Order-consistent price comparisons: ARIMAX fitted with the ARIMA order that auto.arima
#     selects for seasonal ARIMA (domestic prices, price innovations, Brent), and seasonal
#     ARIMA fitted with the ARIMAX order, so that ARIMAX and ARIMA differ only by the prices.
one_origin_rev3 <- function(o) {
  tr <- 1:o; st <- c(as.numeric(format(d$date[1], "%Y")), mon[1])
  yt <- ts(y[tr], start = st, frequency = 12); hmax <- min(12, N - o)
  aa <- function(z, ...) auto.arima(z, ic = "aic", stepwise = TRUE, approximation = FALSE, ...)
  fA <- aa(yt, seasonal = TRUE); fX <- aa(yt, xreg = X[tr, ], seasonal = TRUE)
  oA <- arimaorder(fA); oX <- arimaorder(fX)
  refit <- function(ord, xr = NULL, like) {          # same order and mean/drift terms as 'like'
    cf <- names(coef(like))
    do.call(Arima, list(y = yt, order = ord[1:3], seasonal = list(order = ord[4:6], period = 12), xreg = xr,
          include.mean = "intercept" %in% cf, include.drift = "drift" %in% cf, method = "ML")) }
  fc_x <- function(f, xf) as.numeric(forecast(f, xreg = xf, h = hmax)$mean)
  fXs <- tryCatch(refit(oA, X[tr, ], fA), error = function(e) NULL)       # ARIMAX with the ARIMA order
  fAx <- tryCatch(refit(oX, NULL, fX), error = function(e) NULL)          # ARIMA with the ARIMAX order
  INN <- sapply(PR, function(k) price_innov(X[, k], rows = tr)); INN[is.na(INN)] <- 0
  fSs <- tryCatch(refit(oA, INN[tr, ], fA), error = function(e) NULL)
  fBs <- tryCatch(refit(oA, cbind(BRENT = BRENT[tr]), fA), error = function(e) NULL)
  pXs <- if (is.null(fXs)) rep(NA, hmax) else fc_x(fXs, X[rep(o, hmax), , drop = FALSE])
  pAx <- if (is.null(fAx)) rep(NA, hmax) else as.numeric(forecast(fAx, h = hmax)$mean)
  pSs <- if (is.null(fSs)) rep(NA, hmax) else fc_x(fSs, matrix(0, hmax, length(PR), dimnames = list(NULL, PR)))
  pBs <- if (is.null(fBs)) rep(NA, hmax) else fc_x(fBs, cbind(BRENT = rep(BRENT[o], hmax)))
  pA0 <- as.numeric(forecast(aa(yt, seasonal = TRUE, d = 0, D = 0), h = hmax)$mean)
  pTh <- as.numeric(thetaf(yt, h = hmax)$mean)
  pTB <- as.numeric(forecast(tbats(yt, use.parallel = FALSE), h = hmax)$mean)
  lvl <- mean(y[(o - 11):o]) - mean(y[(o - 23):(o - 12)])
  # seasonal ARDL: order selected with monthly and COVID dummies as fixed regressors
  DX <- cbind(X, SD, COV); dd_ <- cbind(d, SD, COV)
  fm <- as.formula(paste("QIM ~ WHEAT + FUEL + ENERGY + GAS |", paste(c(colnames(SD), colnames(COV)), collapse = " + ")))
  a <- auto_ardl(fm, data = dd_[tr, ], max_order = 6)$best_model; b <- coef(a); b[is.na(b)] <- 0
  yy <- y; XX <- DX; XX[(o + 1):(o + hmax), PR] <- matrix(X[o, ], hmax, length(PR), byrow = TRUE); pDs <- numeric(hmax)
  for (s in 1:hmax) { t <- o + s; v <- b["(Intercept)"]
    for (nm in names(b)[-1]) { lag <- if (grepl("^L\\(", nm)) as.numeric(sub(".*, *([0-9]+)\\)$", "\\1", nm)) else 0
      var <- sub("^L\\(([A-Za-z0-9_]+),.*$", "\\1", nm); v <- v + b[nm] * (if (var == "QIM") yy[t - lag] else XX[t - lag, var]) }
    yy[t] <- v; pDs[s] <- v }
  res <- list()
  for (h in H[H <= hmax]) {
    r <- c(Theta = pTh[h], TBATS = pTB[h], SNaive_adj = y[o + h - 12] + lvl, ARIMA_d0 = pA0[h], ARDL_seas = pDs[h],
           ARIMAX_sameorder = pXs[h], ARIMA_xorder = pAx[h], ARIMAX_shock_sameorder = pSs[h], ARIMAX_brent_sameorder = pBs[h])
    res[[length(res) + 1]] <- data.frame(origin = o, h = h, target = o + h, date = d$date[o + h], actual = y[o + h],
                                         model = names(r), fc = as.numeric(r))
  }
  list(fc = do.call(rbind, res),
       orders = data.frame(origin = o, ARIMA = paste0("(", paste(oA[1:3], collapse = ","), ")(", paste(oA[4:6], collapse = ","), ")"),
                           ARIMAX = paste0("(", paste(oX[1:3], collapse = ","), ")(", paste(oX[4:6], collapse = ","), ")"),
                           d_ARIMA = oA[2], d_ARIMAX = oX[2], D_ARIMA = oA[5], D_ARIMAX = oX[5]))
}
run_rev3 <- function() {
  cf <- file.path(OUT_DIR, "cache", "rolling_expanding_benchmarks.rds")
  if (file.exists(cf)) { cat("  using cached forecasts: additional benchmarks\n"); return(readRDS(cf)) }
  cat(sprintf("  computing additional benchmarks (%d origins, %d core(s)) ...\n", length(origins), N_CORES)); t0 <- Sys.time()
  if (N_CORES > 1 && .Platform$OS.type == "windows") {
    cl <- makeCluster(N_CORES); on.exit(stopCluster(cl))
    clusterEvalQ(cl, suppressPackageStartupMessages({ library(forecast); library(ARDL) }))
    clusterExport(cl, c("one_origin_rev3","price_innov","d","y","X","BRENT","SD","COV","mon","N","PR","H"), envir = environment())
    R <- parLapply(cl, origins, one_origin_rev3)
  } else if (N_CORES > 1) {
    R <- mclapply(origins, one_origin_rev3, mc.cores = N_CORES, mc.preschedule = FALSE)
  } else R <- lapply(origins, one_origin_rev3)
  out <- list(fc = do.call(rbind, lapply(R, `[[`, "fc")), orders = do.call(rbind, lapply(R, `[[`, "orders")))
  out$fc$tag <- "expanding"; out$fc$err <- out$fc$actual - out$fc$fc
  saveRDS(out, cf); cat("  done in", format(round(Sys.time() - t0, 1)), "\n"); out
}
banner("Section 5.8. Additional benchmarks and order-consistent price comparisons (expanding, 6 lags)")
RV3 <- run_rev3()
cmb <- REX %>% filter(model %in% c("SNaive","ARIMA","ETS")) %>% group_by(origin, h, target, date, actual, tag) %>%
  summarise(fc = mean(fc), .groups = "drop") %>% mutate(model = "Comb3", err = actual - fc)
REX3 <- rbind(REX2, RV3$fc[, c("origin","h","target","date","actual","model","fc","tag","err")],
              cmb[, c("origin","h","target","date","actual","model","fc","tag","err")])
cat("Orders selected by auto.arima at the 46 origins: d = 1 for seasonal ARIMA at", sum(RV3$orders$d_ARIMA == 1),
    "origins; d = 1 for ARIMAX at", sum(RV3$orders$d_ARIMAX == 1), "origins; identical orders at",
    sum(RV3$orders$ARIMA == RV3$orders$ARIMAX), "origins.\n")
save_tab(RV3$orders, "Rev3_selected_orders.csv")

accR3 <- do.call(rbind, lapply(names(WIN), function(w) inwin(REX3, w) %>% group_by(model, h) %>%
  summarise(window = w, n = n(), RMSE = sqrt(mean(err^2)), .groups = "drop")))
rmse_rev3 <- as.data.frame(accR3 %>% pivot_wider(id_cols = model, names_from = c(window, h), values_from = RMSE))
rmse_rev3 <- rmse_rev3[, c("model","A_1","A_3","A_6","A_12","B_1","B_3","B_6","B_12")]
save_tab(rmse_rev3, "Rev3_RMSE_all_methods.csv")
E_of <- function(w, hh, models) { s <- inwin(REX3, w) %>% filter(h == hh, model %in% models)   # now on REX3
  list(e = s %>% select(date, model, err) %>% pivot_wider(names_from = model, values_from = err) %>% arrange(date),
       f = s %>% select(date, model, fc) %>% pivot_wider(names_from = model, values_from = fc) %>% arrange(date)) }

# ---- Table 11: additional benchmarks, with model confidence sets on the extended set ----
banner("Table 11. Additional benchmarks: RMSE and membership of the 90% model confidence set (expanding window)")
bench <- c(Comb3 = "Combination: seasonal naive, seasonal ARIMA, ETS", Theta = "Theta", TBATS = "TBATS",
           SNaive_adj = "Seasonal naive, level-adjusted", ARIMA_d0 = "Seasonal ARIMA, d = D = 0",
           ARDL_seas = "ARDL with monthly and COVID dummies")
ext <- c(feasible, names(bench))
MC3 <- list(); DMb <- list()
for (w in names(WIN)) for (hh in H) {
  E <- E_of(w, hh, ext)$e; E <- as.data.frame(E[complete.cases(E), ])
  mcs <- MCSprocedure(as.matrix(E[, ext]^2), alpha = 0.10, B = B_MCS, statistic = "Tmax", k = max(1, hh), verbose = FALSE, seed = 123)
  sh <- mcs@show; pc <- sh[, mcs_col(sh)]
  MC3[[length(MC3) + 1]] <- data.frame(window = w, h = hh, model = rownames(sh), MCS_p = pc, in90 = pc >= 0.10, row.names = NULL)
  best15 <- accR3 %>% filter(window == w, h == hh, model %in% feasible) %>% arrange(RMSE) %>% slice(1) %>% pull(model)
  for (b in names(bench)) { r1 <- dm_hln(E[[b]], E$SNaive, hh); r2 <- if (b == best15) c(NA, NA) else dm_hln(E[[b]], E[[best15]], hh)
    DMb[[length(DMb) + 1]] <- data.frame(window = w, h = hh, model = b, best15 = best15, p_vs_SNaive = r1[2], p_vs_best15 = r2[2]) }
}
MC3 <- do.call(rbind, MC3); DMb <- do.call(rbind, DMb)
tab10b <- rmse_rev3[match(names(bench), rmse_rev3$model), ]; tab10b$model <- bench[tab10b$model]
best_row <- sapply(names(rmse_rev3)[-1], function(cn) min(rmse_rev3[rmse_rev3$model %in% feasible, cn]))
print(r3(rbind(tab10b, data.frame(model = "Lowest RMSE among the 15 methods", t(best_row)))), row.names = FALSE)
cat("Size of the 90% MCS among", length(ext), "methods, and benchmark membership:\n")
print(as.data.frame(MC3 %>% group_by(window, h) %>% summarise(size = sum(in90), Comb3 = in90[model == "Comb3"],
  Theta = in90[model == "Theta"], TBATS = in90[model == "TBATS"], SNadj = in90[model == "SNaive_adj"],
  ARIMA_d0 = in90[model == "ARIMA_d0"], ARDL_seas = in90[model == "ARDL_seas"], .groups = "drop")))
cat("Rank of each benchmark among the 15 methods plus that benchmark:\n")
print(sapply(names(bench), function(b) sapply(names(rmse_rev3)[-1], function(cn) {
  x <- rmse_rev3[rmse_rev3$model %in% c(feasible, b), ]; rank(x[[cn]])[x$model == b] })))
cat("DM tests of the benchmarks (two-sided, unadjusted):\n"); print(r3(DMb))
save_tab(rbind(tab10b, data.frame(model = "Lowest RMSE among the 15 methods", t(best_row))), "Table11_benchmarks.csv")
save_tab(MC3, "Table11_MCS_extended.csv"); save_tab(DMb, "Table11_DM_benchmarks.csv")

# ---- Power: minimum detectable RMSE ratio of the lowest-RMSE method against the seasonal naive ----
banner("Power. Smallest RMSE ratio (method / seasonal naive) that the HLN-DM test could detect at 5%")
mde <- do.call(rbind, lapply(names(WIN), function(w) do.call(rbind, lapply(H, function(hh) {
  best15 <- accR3 %>% filter(window == w, h == hh, model %in% setdiff(feasible, "SNaive")) %>% arrange(RMSE) %>% slice(1) %>% pull(model)
  E <- E_of(w, hh, c(best15, "SNaive"))$e; n <- nrow(E)
  dbar <- mean(E[[best15]]^2 - E$SNaive^2); st <- dm_hln(E[[best15]], E$SNaive, hh)[1]
  thr <- qt(0.975, n - 1) * abs(dbar / st); mse_sn <- mean(E$SNaive^2)
  data.frame(window = w, h = hh, method = lab[best15], n = n, ratio = sqrt(mean(E[[best15]]^2) / mse_sn),
             MDE_ratio = sqrt(max(0, 1 - thr / mse_sn)), row.names = NULL) }))))
print(r3(mde)); save_tab(mde, "Power_MDE.csv")

# ---- Table 14: prices with identical ARIMA orders; DM for the tree models ----
banner("Table 14. Prices: RMSE ratio with / without prices; Clark-West for nested linear pairs, DM otherwise")
t12a <- pair_table(list(c("ARIMAX_sameorder","ARIMA","ARIMAX / seasonal ARIMA, both with the seasonal-ARIMA order"),
                        c("ARIMAX_feasible","ARIMA_xorder","ARIMAX / seasonal ARIMA, both with the ARIMAX order")), "cw")
t12b <- pair_table(list(c("ARIMAX_feasible","ARIMA","ARIMAX / seasonal ARIMA, orders selected separately"),
                        c("GBM","GBM_noCom","GBM / GBM, no prices"), c("RF","RF_noCom","RF / RF, no prices"),
                        c("H1_ARIMA_GBM","H1_noCom","H1 / H1, no prices"), c("H3_STL_GBM","H3_STL_GBM_noCom","H3 / H3, no prices")), "dm")
t12a$test <- "CW"; t12b$test <- "DM"; tab12r <- rbind(t12a, t12b)
print(wide_cells(tab12r), right = FALSE)
cat("Ratios below one:", sum(tab12r$ratio < 1), "of", nrow(tab12r), "\n")
print(as.data.frame(tab12r %>% group_by(comparison, test) %>% summarise(below1 = sum(ratio < 1), p05 = sum(p < 0.05),
      holm_table_min = NA, .groups = "drop")))
tab12r$p_holm_table <- p.adjust(tab12r$p, "holm"); cat("Smallest Holm-adjusted p across all", nrow(tab12r), "tests:", round(min(tab12r$p_holm_table), 3), "\n")
print(tab12r[tab12r$p < 0.05, c("comparison","window","h","ratio","p","p_holm_table")])
cat("Conditional ARIMAX (realised future prices) RMSE, for reference:\n"); print(r3(rmse_rev3[rmse_rev3$model == "ARIMAX_conditional", ]))
save_tab(tab12r, "Table14_prices.csv")

# ---- Table A8: shocks and import costs; same-order ARIMAX, DM for tree models ----
banner("Table A8. Price shocks and import costs: RMSE ratio to the same method without price information")
a8a <- pair_table(list(c("ARIMAX_sameorder","ARIMA","ARIMAX: price levels"), c("ARIMAX_shock_sameorder","ARIMA","ARIMAX: price innovations"),
                       c("ARIMAX_brent_sameorder","ARIMA","ARIMAX: Brent crude price")), "cw")
a8b <- pair_table(list(c("GBM","GBM_noCom","GBM: price levels"), c("GBM_shock","GBM_noCom","GBM: price innovations"),
  c("GBM_nopi","GBM_noCom","GBM: net price increases"), c("GBM_brent","GBM_noCom","GBM: Brent crude price"),
  c("H1_ARIMA_GBM","H1_noCom","H1: price levels"), c("H1_shock","H1_noCom","H1: price innovations"),
  c("H1_brent","H1_noCom","H1: Brent crude price")), "dm")
a8a$test <- "CW"; a8b$test <- "DM"; tabA8r <- rbind(a8b, a8a)
print(wide_cells(tabA8r), right = FALSE)
print(as.data.frame(tabA8r %>% group_by(comparison) %>% summarise(below1 = sum(ratio < 1), p05 = sum(p < 0.05), .groups = "drop")))
tabA8r$p_holm_table <- p.adjust(tabA8r$p, "holm"); cat("Smallest Holm-adjusted p across all", nrow(tabA8r), "tests:", round(min(tabA8r$p_holm_table), 3), "\n")
print(tabA8r[tabA8r$p < 0.05, c("comparison","window","h","ratio","p","p_holm_table")])
save_tab(tabA8r, "TableA8_price_shocks_import_costs.csv")

# =============================================================================
#  APPENDIX B. ADDITIONAL TABLES
# =============================================================================
banner("Table A1. Complete unrestricted error-correction estimates (dependent variable: change in QIM)")
print(r3(tabB12)); cat(B12_note, "\n"); save_tab(tabB12, "TableA1_UECM.csv")

banner("Table A2. Out-of-sample RMSE, fixed 66-month rolling window")
tabB13 <- rmse_table("rolling"); print(tabB13, right = FALSE); save_tab(tabB13, "TableA2_RMSE_rolling.csv")

banner("Table A3. MAE, MAPE and MASE one month ahead, expanding window")
tabB14 <- acc %>% filter(tag == "expanding", h == 1) %>% select(model, window, MAE, MAPE, MASE) %>%
  pivot_wider(names_from = window, values_from = c(MAE, MAPE, MASE))
tabB14 <- tabB14[match(names(lab), tabB14$model), c("model","MAE_A","MAPE_A","MASE_A","MAE_B","MAPE_B","MASE_B")]
tabB14$model <- lab[tabB14$model]; print(as.data.frame(r3(tabB14)), right = FALSE); save_tab(tabB14, "TableA3_MAE_MAPE_MASE.csv")

banner("Table A4. Diebold-Mariano tests of H3 against each competitor (DM statistic, Holm-adjusted p)")
tabB15 <- DM %>% filter(tag == "expanding", set == "H3", h %in% c(1, 12)) %>%
  mutate(cell = sprintf("%.2f (%.3f)", DM, p_holm)) %>% select(comparator, window, h, cell) %>%
  pivot_wider(names_from = c(window, h), values_from = cell)
tabB15 <- tabB15[match(intersect(names(lab), tabB15$comparator), tabB15$comparator), c("comparator","A_1","A_12","B_1","B_12")]
tabB15$comparator <- lab[tabB15$comparator]; print(as.data.frame(tabB15), right = FALSE)
save_tab(tabB15, "TableA4_DM.csv")

banner(paste("Finished. All tables and figures are in:", normalizePath(OUT_DIR)))
sessionInfo()
