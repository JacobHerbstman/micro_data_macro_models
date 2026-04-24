###############################################################################
# Homework 3: Migration Trends and Selection
###############################################################################

#setwd("/Users/jacobherbstman/Desktop/micro_data_macro_models/assignments/code")

rm(list = ls())

library(ipumsr)
library(data.table)
library(ggplot2)

xml_file <- "../input/User Extract usa_00012.dat.xml"
output_dir <- "../output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

latex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  escaped <- vapply(
    strsplit(x, "", useBytes = TRUE),
    function(chars) {
      paste0(
        vapply(
          chars,
          function(ch) {
            switch(
              ch,
              "\\" = "\\textbackslash{}",
              "&" = "\\&",
              "%" = "\\%",
              "$" = "\\$",
              "#" = "\\#",
              "_" = "\\_",
              "{" = "\\{",
              "}" = "\\}",
              "~" = "\\textasciitilde{}",
              "^" = "\\textasciicircum{}",
              ch
            )
          },
          character(1)
        ),
        collapse = ""
      )
    },
    character(1)
  )
  escaped
}

format_decimal <- function(x, digits = 3) {
  sprintf(paste0("%.", digits, "f"), x)
}

write_msa_table_tex <- function(x, file_name, value_col, value_label, digits = 3) {
  rows <- sprintf(
    "    %s & %s %s",
    latex_escape(x$msa_name),
    format_decimal(x[[value_col]], digits),
    "\\\\"
  )
  table_lines <- c(
    "    \\centering",
    "    \\small",
    "    \\begin{tabular}{p{0.66\\textwidth}r}",
    "    \\toprule",
    sprintf("    MSA & %s \\\\", value_label),
    "    \\midrule",
    rows,
    "    \\bottomrule",
    "    \\end{tabular}"
  )
  writeLines(table_lines, file.path(output_dir, file_name), useBytes = TRUE)
}

write_selection_regression_table_tex <- function(x, file_name, digits = 3) {
  panel_label <- fifelse(
    grepl("In-migrant", x$panel),
    "A: In-migrants",
    "B: Out-migrants"
  )
  rows <- sprintf(
    "    %s & %s & %s & %s & %s %s",
    panel_label,
    format_decimal(x$slope, digits),
    format_decimal(x$standard_error, digits),
    format_decimal(x$weighted_r_squared, digits),
    format_decimal(x$p_beta_eq_1, digits),
    "\\\\"
  )
  table_lines <- c(
    "    \\centering",
    "    \\small",
    "    \\begin{tabular}{lrrrr}",
    "    \\toprule",
    "    Panel & Slope & SE & $R^2$ & $p(\\beta = 1)$ \\\\",
    "    \\midrule",
    rows,
    "    \\bottomrule",
    "    \\end{tabular}"
  )
  writeLines(table_lines, file.path(output_dir, file_name), useBytes = TRUE)
}

###############################################################################
# Read data and keep only needed variables
###############################################################################

ddi <- read_ipums_ddi(xml_file)

state_lookup <- as.data.table(ipums_val_labels(ddi, "STATEFIP"))
setnames(state_lookup, c("val", "lbl"), c("state_code", "state_name"))

msa_lookup <- as.data.table(ipums_val_labels(ddi, "MET2013"))
setnames(msa_lookup, c("val", "lbl"), c("msa_code", "msa_name"))

dt <- read_ipums_micro(
  ddi,
  vars = c(
    "YEAR", "GQ", "STATEFIP", "MET2013", "PERWT", "AGE", "BPL", "EDUC",
    "EMPSTAT", "MIGRATE1", "MIGPLAC1", "MIGMET131"
  ),
  verbose = FALSE
)

setDT(dt)
setnames(dt, tolower(names(dt)))

integer_cols <- c(
  "year", "gq", "statefip", "met2013", "age", "bpl", "educ",
  "empstat", "migrate1", "migplac1", "migmet131"
)

for (col_name in integer_cols) {
  set(dt, j = col_name, value = as.integer(dt[[col_name]]))
}

set(dt, j = "perwt", value = as.numeric(dt[["perwt"]]))

###############################################################################
# Early checks and household sample restriction
###############################################################################

cat("\nAge and birthplace checks\n")
cat("------------------------\n")
cat("Rows in raw extract:", format(nrow(dt), big.mark = ","), "\n")
cat("Min age:", min(dt$age, na.rm = TRUE), "\n")
cat("Max age:", max(dt$age, na.rm = TRUE), "\n")
cat("Max birthplace code:", max(dt$bpl, na.rm = TRUE), "\n")

stopifnot(min(dt$age, na.rm = TRUE) == 25L)
stopifnot(max(dt$age, na.rm = TRUE) == 54L)
stopifnot(max(dt$bpl, na.rm = TRUE) <= 56L)

gq_check <- dt[, .N, by = gq][order(gq)]
cat("\nGroup quarters tabulation before household restriction\n")
cat("------------------------------------------------------\n")
print(gq_check)

dt <- dt[gq %in% c(1L, 2L)]

gq_after <- dt[, .N, by = gq][order(gq)]
cat("\nGroup quarters tabulation after household restriction\n")
cat("-----------------------------------------------------\n")
print(gq_after)

###############################################################################
# Education groups
###############################################################################

dt[, education_group := NA_character_]
dt[educ %in% 1:5, education_group := "High school dropout"]
dt[educ == 6, education_group := "High school only"]
dt[educ %in% 7:9, education_group := "Some college, no BA"]
dt[educ == 10, education_group := "BA only"]
dt[educ == 11, education_group := "Post-BA"]

educ_check <- dt[
  ,
  .N,
  by = .(educ, education_group)
][order(educ, education_group)]

cat("\nEDUC to education-group crosswalk\n")
cat("---------------------------------\n")
print(educ_check)

###############################################################################
# Geography availability checks
###############################################################################

metro_check <- dt[
  ,
  .(
    current_msa_positive = sum(met2013 > 0, na.rm = TRUE),
    prior_msa_positive = sum(migmet131 > 0, na.rm = TRUE)
  ),
  by = year
][order(year)]

cat("\nMSA availability by year\n")
cat("------------------------\n")
print(metro_check)

###############################################################################
# Interstate migration trend by education
###############################################################################

trend_dt <- dt[
  !is.na(education_group) & migrate1 %in% c(1L, 2L, 3L, 4L),
  .(
    population = sum(perwt),
    interstate_movers = sum(perwt[migrate1 == 3L])
  ),
  by = .(year, education_group)
][order(year, education_group)]

trend_dt[, interstate_migration_rate := interstate_movers / population]

fwrite(trend_dt, file.path(output_dir, "hw3_migration_trend.csv"))

trend_selected_dt <- trend_dt[
  year %in% c(2000L, 2004L, 2019L, 2020L, 2024L)
][order(year, education_group)]

fwrite(
  trend_selected_dt,
  file.path(output_dir, "hw3_migration_trend_selected_years.csv")
)

education_levels <- c(
  "High school dropout",
  "High school only",
  "Some college, no BA",
  "BA only",
  "Post-BA"
)

trend_plot <- ggplot(
  trend_dt[
    ,
    education_group := factor(education_group, levels = education_levels)
  ],
  aes(x = year, y = interstate_migration_rate, color = education_group)
) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = seq(2000, 2024, by = 2)) +
  labs(
    title = "Interstate migration rates by education group",
    x = "Year",
    y = "Interstate migration rate",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 0.1),
    expand = ggplot2::expansion(mult = c(0.02, 0.05))
  )

ggsave(
  filename = file.path(output_dir, "hw3_migration_trend.pdf"),
  plot = trend_plot,
  width = 8,
  height = 5
)

###############################################################################
# Top 75 MSAs by 2006 population
###############################################################################

top75_msas <- dt[
  year == 2006L & met2013 > 0,
  .(population_2006 = sum(perwt)),
  by = met2013
][order(-population_2006)][1:75]

setnames(top75_msas, "met2013", "msa_code")
top75_msas <- merge(top75_msas, msa_lookup, by = "msa_code", all.x = TRUE)
setorder(top75_msas, -population_2006)
setcolorder(top75_msas, c("msa_code", "msa_name", "population_2006"))

cat("\nTop-75 MSA check\n")
cat("----------------\n")
cat("Unique positive MSA codes:", uniqueN(top75_msas$msa_code), "\n")

fwrite(top75_msas, file.path(output_dir, "hw3_top75_msas.csv"))

top75_codes <- top75_msas$msa_code

###############################################################################
# MSA inflows, outflows, and flow rates
###############################################################################

flow_years <- c(2006L, 2010L, 2019L)
msa_year_grid <- CJ(year = flow_years, msa_code = top75_codes)

msa_population <- dt[
  year %in% flow_years & met2013 %in% top75_codes,
  .(population = sum(perwt)),
  by = .(year, msa_code = met2013)
]

msa_inflows <- dt[
  year %in% flow_years & met2013 %in% top75_codes &
    migrate1 %in% c(2L, 3L) &
    migmet131 > 0 & migmet131 != met2013,
  .(inflows = sum(perwt)),
  by = .(year, msa_code = met2013)
]

msa_outflows <- dt[
  year %in% flow_years & migmet131 %in% top75_codes &
    migrate1 %in% c(2L, 3L) &
    met2013 > 0 &
    met2013 != migmet131,
  .(outflows = sum(perwt)),
  by = .(year, msa_code = migmet131)
]

prior_msa_support <- dt[
  year %in% flow_years & migmet131 %in% top75_codes,
  .(
    prior_msa_support_n = .N,
    prior_msa_support_weight = sum(perwt)
  ),
  by = .(year, msa_code = migmet131)
]

msa_flows <- merge(msa_year_grid, msa_population, by = c("year", "msa_code"), all.x = TRUE)
msa_flows <- merge(msa_flows, msa_inflows, by = c("year", "msa_code"), all.x = TRUE)
msa_flows <- merge(msa_flows, msa_outflows, by = c("year", "msa_code"), all.x = TRUE)
msa_flows <- merge(msa_flows, prior_msa_support, by = c("year", "msa_code"), all.x = TRUE)
msa_flows <- merge(msa_flows, top75_msas, by = "msa_code", all.x = TRUE)

for (col_name in c("population", "inflows")) {
  set(msa_flows, which(is.na(msa_flows[[col_name]])), col_name, 0)
}

msa_flows[, prior_msa_supported := !is.na(prior_msa_support_n)]
msa_flows[prior_msa_supported == TRUE & is.na(outflows), outflows := 0]
msa_flows[
  prior_msa_supported == FALSE | is.na(prior_msa_supported),
  outflows := NA_real_
]
msa_flows[, `:=`(gross_flow_rate = NA_real_, net_flow_rate = NA_real_)]
msa_flows[
  !is.na(outflows) & population > 0,
  `:=`(
    gross_flow_rate = (inflows + outflows) / population,
    net_flow_rate = (inflows - outflows) / population
  )
]
setcolorder(
  msa_flows,
  c(
    "year", "msa_code", "msa_name", "population_2006", "population",
    "inflows", "outflows", "prior_msa_supported",
    "prior_msa_support_n", "prior_msa_support_weight",
    "gross_flow_rate", "net_flow_rate"
  )
)
setorder(msa_flows, year, -population)

fwrite(msa_flows, file.path(output_dir, "hw3_msa_flows.csv"))
fwrite(prior_msa_support, file.path(output_dir, "hw3_prior_msa_support.csv"))

cat("\nSpot check of MSA flow rows\n")
cat("---------------------------\n")
print(msa_flows[year == 2006L][1:10])

rate_check <- msa_flows[
  year == 2006L & population > 0 & !is.na(outflows),
  .(
    gross_ok = all.equal(gross_flow_rate[1], (inflows[1] + outflows[1]) / population[1]),
    net_ok = all.equal(net_flow_rate[1], (inflows[1] - outflows[1]) / population[1])
  ),
  by = .(msa_code)
][1:5]

cat("\nFlow-rate formula check\n")
cat("-----------------------\n")
print(rate_check)

cat("\nPrior-MSA support check\n")
cat("-----------------------\n")
print(msa_flows[, .(unsupported_msa_years = sum(!prior_msa_supported)), by = year])

gross_flow_top10_2006 <- msa_flows[
  year == 2006L & !is.na(gross_flow_rate)
][order(-gross_flow_rate)][1:10]
net_flow_top10_2006 <- msa_flows[
  year == 2006L & !is.na(net_flow_rate)
][order(-net_flow_rate)][1:10]

fwrite(gross_flow_top10_2006, file.path(output_dir, "hw3_gross_flow_top10_2006.csv"))
fwrite(net_flow_top10_2006, file.path(output_dir, "hw3_net_flow_top10_2006.csv"))

msa_flow_wide_input <- copy(msa_flows)
setnames(msa_flow_wide_input, "population_2006", "baseline_population_2006")
msa_flow_wide <- dcast(
  msa_flow_wide_input,
  msa_code + msa_name + baseline_population_2006 ~ year,
  value.var = c("population", "inflows", "outflows", "gross_flow_rate", "net_flow_rate")
)

msa_flow_wide[, net_flow_change_2006_2010 := net_flow_rate_2010 - net_flow_rate_2006]
msa_flow_wide[, net_flow_change_2006_2019 := net_flow_rate_2019 - net_flow_rate_2006]

net_flow_change_top10_2006_2010 <- msa_flow_wide[
  !is.na(net_flow_change_2006_2010)
][order(-net_flow_change_2006_2010)][1:10]
net_flow_change_bottom10_2006_2010 <- msa_flow_wide[
  !is.na(net_flow_change_2006_2010)
][order(net_flow_change_2006_2010)][1:10]
net_flow_change_top10_2006_2019 <- msa_flow_wide[
  !is.na(net_flow_change_2006_2019)
][order(-net_flow_change_2006_2019)][1:10]
net_flow_change_bottom10_2006_2019 <- msa_flow_wide[
  !is.na(net_flow_change_2006_2019)
][order(net_flow_change_2006_2019)][1:10]

fwrite(
  net_flow_change_top10_2006_2010,
  file.path(output_dir, "hw3_net_flow_change_top10_2006_2010.csv")
)
fwrite(
  net_flow_change_bottom10_2006_2010,
  file.path(output_dir, "hw3_net_flow_change_bottom10_2006_2010.csv")
)
fwrite(
  net_flow_change_top10_2006_2019,
  file.path(output_dir, "hw3_net_flow_change_top10_2006_2019.csv")
)
fwrite(
  net_flow_change_bottom10_2006_2019,
  file.path(output_dir, "hw3_net_flow_change_bottom10_2006_2019.csv")
)

write_msa_table_tex(
  gross_flow_top10_2006,
  "hw3_gross_flow_top10_2006_table.tex",
  "gross_flow_rate",
  "Gross flow rate"
)
write_msa_table_tex(
  net_flow_top10_2006,
  "hw3_net_flow_top10_2006_table.tex",
  "net_flow_rate",
  "Net flow rate"
)
write_msa_table_tex(
  net_flow_change_top10_2006_2010,
  "hw3_net_flow_change_top10_2006_2010_table.tex",
  "net_flow_change_2006_2010",
  "Change in net flow rate"
)
write_msa_table_tex(
  net_flow_change_bottom10_2006_2010,
  "hw3_net_flow_change_bottom10_2006_2010_table.tex",
  "net_flow_change_2006_2010",
  "Change in net flow rate"
)
write_msa_table_tex(
  net_flow_change_top10_2006_2019,
  "hw3_net_flow_change_top10_2006_2019_table.tex",
  "net_flow_change_2006_2019",
  "Change in net flow rate"
)
write_msa_table_tex(
  net_flow_change_bottom10_2006_2019,
  "hw3_net_flow_change_bottom10_2006_2019_table.tex",
  "net_flow_change_2006_2019",
  "Change in net flow rate"
)

###############################################################################
# Employment rates for residents, in-migrants, and out-migrants
###############################################################################

employment_years <- c(2005L, 2010L, 2019L)
employment_grid <- CJ(year = employment_years, msa_code = top75_codes)

resident_employment <- dt[
  year %in% employment_years & met2013 %in% top75_codes & empstat %in% c(1L, 2L, 3L),
  .(
    resident_population = sum(perwt),
    resident_employed = sum(perwt[empstat == 1L])
  ),
  by = .(year, msa_code = met2013)
]

resident_employment[, resident_employment_rate := resident_employed / resident_population]

# For the employment-selection exercise, keep movers tied to the focal MSA even
# when one side of the metro move is not fully identifiable. Exclude only people
# who are explicitly observed in the same MSA one year earlier.
inmigrant_employment <- dt[
  year %in% employment_years & met2013 %in% top75_codes &
    migrate1 %in% c(2L, 3L, 4L) &
    !(migmet131 > 0 & migmet131 == met2013) &
    empstat %in% c(1L, 2L, 3L),
  .(
    inmigrant_population = sum(perwt),
    inmigrant_employed = sum(perwt[empstat == 1L])
  ),
  by = .(year, msa_code = met2013)
]

inmigrant_employment[, inmigrant_employment_rate := inmigrant_employed / inmigrant_population]

outmigrant_employment <- dt[
  year %in% employment_years & migmet131 %in% top75_codes &
    migrate1 %in% c(2L, 3L, 4L) &
    met2013 != migmet131 &
    empstat %in% c(1L, 2L, 3L),
  .(
    outmigrant_population = sum(perwt),
    outmigrant_employed = sum(perwt[empstat == 1L])
  ),
  by = .(year, msa_code = migmet131)
]

outmigrant_employment[, outmigrant_employment_rate := outmigrant_employed / outmigrant_population]

msa_employment <- merge(employment_grid, resident_employment, by = c("year", "msa_code"), all.x = TRUE)
msa_employment <- merge(msa_employment, inmigrant_employment, by = c("year", "msa_code"), all.x = TRUE)
msa_employment <- merge(msa_employment, outmigrant_employment, by = c("year", "msa_code"), all.x = TRUE)
msa_employment <- merge(msa_employment, top75_msas, by = "msa_code", all.x = TRUE)
setcolorder(
  msa_employment,
  c(
    "year", "msa_code", "msa_name", "population_2006",
    "resident_population", "resident_employed", "resident_employment_rate",
    "inmigrant_population", "inmigrant_employed", "inmigrant_employment_rate",
    "outmigrant_population", "outmigrant_employed", "outmigrant_employment_rate"
  )
)
setorder(msa_employment, year, -population_2006)

fwrite(msa_employment, file.path(output_dir, "hw3_msa_employment_rates.csv"))

###############################################################################
# 2005 to 2010 employment-rate change scatter
###############################################################################

scatter_data <- dcast(
  msa_employment,
  msa_code + msa_name + population_2006 ~ year,
  value.var = c(
    "resident_employment_rate",
    "inmigrant_employment_rate",
    "outmigrant_employment_rate",
    "inmigrant_population",
    "outmigrant_population"
  )
)

scatter_data[, resident_change := resident_employment_rate_2010 - resident_employment_rate_2005]
scatter_data[, inmigrant_change := inmigrant_employment_rate_2010 - inmigrant_employment_rate_2005]
scatter_data[, outmigrant_change := outmigrant_employment_rate_2010 - outmigrant_employment_rate_2005]
scatter_data[, inmigrant_avg_population := (inmigrant_population_2005 + inmigrant_population_2010) / 2]
scatter_data[, outmigrant_avg_population := (outmigrant_population_2005 + outmigrant_population_2010) / 2]

cat("\nScatter-data check\n")
cat("------------------\n")
cat(
  "Rows with non-missing resident and in-migrant changes:",
  scatter_data[!is.na(resident_change) & !is.na(inmigrant_change), .N],
  "\n"
)
cat(
  "Rows with non-missing resident and out-migrant changes:",
  scatter_data[!is.na(resident_change) & !is.na(outmigrant_change), .N],
  "\n"
)

fwrite(scatter_data, file.path(output_dir, "hw3_selection_scatter_data.csv"))

panel_a_dt <- scatter_data[
  !is.na(resident_change) & !is.na(inmigrant_change),
  .(
    msa_code,
    msa_name,
    resident_change,
    outcome_change = inmigrant_change,
    avg_population = inmigrant_avg_population,
    rate_2005 = inmigrant_employment_rate_2005,
    rate_2010 = inmigrant_employment_rate_2010,
    pop_2005 = inmigrant_population_2005,
    pop_2010 = inmigrant_population_2010
  )
]

panel_b_dt <- scatter_data[
  !is.na(resident_change) & !is.na(outmigrant_change),
  .(
    msa_code,
    msa_name,
    resident_change,
    outcome_change = outmigrant_change,
    avg_population = outmigrant_avg_population,
    rate_2005 = outmigrant_employment_rate_2005,
    rate_2010 = outmigrant_employment_rate_2010,
    pop_2005 = outmigrant_population_2005,
    pop_2010 = outmigrant_population_2010
  )
]

fit_panel_a <- lm(
  outcome_change ~ resident_change,
  data = panel_a_dt,
  weights = avg_population
)

fit_panel_b <- lm(
  outcome_change ~ resident_change,
  data = panel_b_dt,
  weights = avg_population
)

panel_a_coef <- coef(summary(fit_panel_a))["resident_change", ]
panel_b_coef <- coef(summary(fit_panel_b))["resident_change", ]

panel_a_beta_one_t <- (panel_a_coef["Estimate"] - 1) / panel_a_coef["Std. Error"]
panel_b_beta_one_t <- (panel_b_coef["Estimate"] - 1) / panel_b_coef["Std. Error"]

selection_regressions <- data.table(
  panel = c(
    "Panel A: In-migrant employment rate change",
    "Panel B: Out-migrant employment rate change"
  ),
  slope = c(panel_a_coef["Estimate"], panel_b_coef["Estimate"]),
  standard_error = c(panel_a_coef["Std. Error"], panel_b_coef["Std. Error"]),
  weighted_r_squared = c(summary(fit_panel_a)$r.squared, summary(fit_panel_b)$r.squared),
  n_msas = c(nrow(panel_a_dt), nrow(panel_b_dt)),
  p_beta_eq_1 = c(
    2 * pt(abs(panel_a_beta_one_t), df = fit_panel_a$df.residual, lower.tail = FALSE),
    2 * pt(abs(panel_b_beta_one_t), df = fit_panel_b$df.residual, lower.tail = FALSE)
  )
)

fwrite(
  selection_regressions,
  file.path(output_dir, "hw3_selection_scatter_regressions.csv")
)
write_selection_regression_table_tex(
  selection_regressions,
  "hw3_selection_scatter_regressions_table.tex"
)

cat("\nWeighted selection regressions\n")
cat("------------------------------\n")
print(selection_regressions)

panel_a_dt[, panel := "Panel A: In-migrant employment rate change"]
panel_b_dt[, panel := "Panel B: Out-migrant employment rate change"]

scatter_plot_dt <- rbindlist(
  list(panel_a_dt, panel_b_dt),
  use.names = TRUE,
  fill = TRUE
)

pad_range <- function(x, mult = 0.06) {
  r <- range(x, na.rm = TRUE)
  d <- diff(r)
  r + c(-1, 1) * mult * d
}

x_limits <- pad_range(scatter_plot_dt$resident_change, mult = 0.08)
y_limits <- pad_range(scatter_plot_dt$outcome_change, mult = 0.06)

selection_plot <- ggplot(
  scatter_plot_dt,
  aes(x = resident_change, y = outcome_change)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "firebrick") +
  geom_smooth(
    aes(weight = avg_population),
    method = "lm",
    se = FALSE,
    color = "steelblue",
    linewidth = 0.8
  ) +
  geom_point(aes(size = avg_population), alpha = 0.75) +
  facet_wrap(~panel, ncol = 2, scales = "fixed") +
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  scale_x_continuous(
    labels = scales::label_percent(accuracy = 1),
    breaks = scales::breaks_width(0.025)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 5),
    breaks = scales::breaks_width(0.05)
  ) +
  scale_size_continuous(name = "Avg weighted\nmovers", range = c(1.5, 4.2)) +
  labs(
    title = "Changes in MSA employment rates, 2005 to 2010",
    subtitle = "Dashed line is the no-selection benchmark; fitted lines are weighted by mover counts",
    x = "Change in resident employment rate",
    y = "Change in mover employment rate"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(output_dir, "hw3_selection_scatter.pdf"),
  plot = selection_plot,
  width = 10,
  height = 5
)
