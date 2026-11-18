# Hybrid-Econometric-ML-Industrial-Production

Data and R code for **"How evaluation design shapes the apparent accuracy of hybrid econometric–machine learning forecasts: Pakistan's industrial production under domestic energy and food price shocks"** (Yaqoob, Ahmed, Danish, Jabbar and Wascher; *Forecasting*, MDPI).

## Contents

| File | Purpose |
|---|---|
| `Industrial Production Data.xlsx` | Quantum Index of Manufacturing (QIM), PBS |
| `Wheat data.xlsx` | Wheat retail price, PBS CPI average prices |
| `Energy data.xlsx` | Petrol, diesel, electricity and natural-gas prices, PBS CPI average prices |
| `Brent data.csv` | Brent crude oil spot price, USD per barrel, monthly (U.S. EIA, series RBRTE); import-cost proxy added in the second revision |
| `R_script_to_replicate_results.R` | Reproduces every table and figure of the paper (Tables 1–15, 17 and 18, Tables A1–A14, Figures 1–3 and A1; Table 16 is a summary table without separate output) |

## Data: source and transformation of each variable

All series are monthly and published by the Pakistan Bureau of Statistics (PBS): QIM at <https://www.pbs.gov.pk/content/qim>, prices at <https://www.pbs.gov.pk/price-statistics/>.

| Variable (code) | File, sheet, columns | Unit | Transformation |
|---|---|---|---|
| QIM (`QIM`) | `Industrial Production Data.xlsx`, sheet `Trend Sheet (2)`, column 1 (date) and column 2 (QIM); first 5 rows skipped | Index, 2015–16 = 100 | Level as published, not seasonally adjusted; dates set to the 1st of the month |
| Wheat (`WHEAT`) | `Wheat data.xlsx`, `Sheet1`, right-hand block "1516" (2015–16 base): column E (month label `Mon-YY`) and column F ("Pakistan", Wheat 10 kg); first 8 rows skipped | PKR per 10 kg | National average; the import starts in August 2016 (the July 2016 value is not read) |
| Fuel (`FUEL`) | `Energy data.xlsx`, sheet `Base Year 2007-08`, columns `Petrol (Per Litre)` and `Diesel (Per Litre)`; first row skipped | PKR per litre | (petrol + diesel) / 2 |
| Electricity (`ENERGY`) | same sheet, column `Electricity (Per Unit)` | PKR per kWh | As recorded |
| Natural gas (`GAS`) | same sheet, column `Gas (MMBTU)` | PKR per MMBTU | As recorded |
| Brent (`BRENT`) | `Brent data.csv`, columns `Month` (YYYY-MM) and `Brent_USD_per_barrel`; source <https://www.eia.gov/dnav/pet/hist/LeafHandler.ashx?n=PET&s=RBRTE&f=M> | USD per barrel | Monthly average spot price; used only in the robustness checks of Appendix C |

The three files are merged by calendar month (inner join). The merged sample runs from **August 2016 to November 2025 (N = 112)** and has no missing values. All of these steps are in Section 2 of `R_script_to_replicate_results.R`.

## How to reproduce the results

1. Install R (tested with **R 4.3.3**) and the packages:
   ```r
   install.packages(c("readxl","dplyr","tidyr","zoo","urca","tseries","lmtest","sandwich",
                      "strucchange","FinTS","ARDL","forecast","gbm","randomForest","MCS","ggplot2",
                      "MSwM","quantreg"))
   ```
2. Put the three Excel files, `Brent data.csv` and the script in one folder.
3. In `R_script_to_replicate_results.R`, set `DATA_DIR` (the folder with the Excel files) and `OUT_DIR` (the output folder) at the top. Then run the whole script. It writes every table (`Table*.csv`, numbered as in the paper) and figure (`Figure*.pdf`) to `OUT_DIR`.
   * Sections 2–4 take about two minutes, including the unit root tests (AIC lag selection from 1 to 9 lags), the preferred seasonal ARDL and log-price models (Table 4, panel B) and the ARDL robustness checks of Table 6.
   * The rolling-origin forecasts take about 20–25 minutes on 8 cores. Section 5 covers 46 origins (January 2022–October 2025) and four schemes.
   * The additional variants of Section 5.6 add about 10 minutes on 4 cores. These are the horizon-specific correctors, alternative hyperparameters, forecast price paths and extra price ablations for Tables 10, 17, 18 and A5.
   * The variants of Section 5.7 add about 8 minutes on 4 cores. These are the price-shock and import-cost measures, small-sample remedies, direct hybrids, price-path scenarios, pseudo-real-time forecasts and bootstrap intervals for Tables A9–A14. The price-shock table (A6) is produced in Section 2, and the structural-break analyses (Table A7, Figure A1) in Section 3.
   * Section 5.8 adds about 4 minutes on 4 cores. It computes the six additional benchmarks (Theta, TBATS, level-adjusted seasonal naive, seasonal ARIMA with d = D = 0, ARDL with monthly and COVID dummies, and the combination of seasonal naive, seasonal ARIMA and ETS) for Table 11, the power calculation for the DM tests, and the order-consistent price comparisons for Table 14 and Table A8, in which ARIMAX and seasonal ARIMA are fitted with the same ARIMA order. Earlier versions of the price tables that compared models with separately selected orders are still written, with `_superseded` in the file name.
   * Forecasts are cached in `OUT_DIR/cache`, so a second run is fast.

Random seeds are fixed, so the numbers are reproducible. Small differences in the last digit of bootstrap p-values can arise across package versions.

## Forecasting implementation (summary)

* **Origins and windows.** At every origin *o* (the last observed month) every model is re-estimated on data up to *o* only. The windows are expanding, or rolling with a fixed length of 66 months. The horizons are *h* = 1, 3, 6, 12. Evaluation window A covers January 2024–November 2025 and window B January 2023–November 2025.
* **Methods.** The paper compares 15 feasible methods. A conditional ARIMAX with realised future prices is shown for reference only and never enters rankings, tests or model confidence sets.
* **Machine learning.** The features are six lags of QIM and of each price dated at *o*, the calendar month of the target and a trend (32 features). The direct strategy fits one model per horizon. The number of GBM trees is chosen on the last 12 months of the estimation window, and the model is then refitted to the whole window.
* **Hybrids.** Forecast = linear-model *h*-step forecast + ML forecast of the linear model's residual. In the baseline, the corrector is trained on one-step residuals; in Section 5.6 of the script, also on *h*-step in-sample errors.
* **Tests.** Diebold–Mariano with the Harvey–Leybourne–Newbold correction and Holm adjustment; model confidence set (T_max, 5,000 block-bootstrap replications); and, in the revision, Clark–West tests for nested comparisons.

## License

See `LICENSE`.
