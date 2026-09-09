suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

make_pbc_age <- function() {
  data <- as.data.table(survival::pbcseq)
  setorder(data, id, day)
  baseline <- data[, .SD[1L], by = id]
  baseline <- baseline[, .(
    id, entry_age = age, exit_age = age + futime / 365.25,
    event = as.integer(status == 2L), cr_status = status,
    male = as.integer(sex == "m"),
    stage0 = stage, trt0 = trt, bili0 = bili
  )]
  baseline <- baseline[complete.cases(baseline)]
  baseline[, log_bili_z := as.numeric(scale(log(bili0)))]
  stopifnot(!anyNA(baseline), all(baseline$entry_age < baseline$exit_age))

  data[, bili_locf := nafill(bili, type = "locf"), by = id]
  data[, albumin_locf := nafill(albumin, type = "locf"), by = id]
  data <- merge(
    data[, .(id, day, bili_locf, albumin_locf)],
    baseline[, .(id, entry_age, exit_age, event, cr_status, male, stage0, trt0)],
    by = "id", all = FALSE
  )
  data <- data[complete.cases(data)]
  data[, visit_age := entry_age + day / 365.25]
  data <- data[visit_age < exit_age]
  data[, log_bili := log(bili_locf)]
  stopifnot(!anyDuplicated(data[, .(id, visit_age)]), !anyNA(data),
            all(data$visit_age >= data$entry_age))
  list(long = as.data.frame(data), baseline = as.data.frame(baseline))
}


make_nafld_age <- function(include_preindex_labs = FALSE) {
  baseline <- as.data.table(survival::nafld1)
  baseline <- baseline[complete.cases(baseline[, .(id, age, futime, status, male, bmi)])]
  baseline <- baseline[, .(
    id, entry_age = age, exit_age = age + futime / 365.25,
    event = as.integer(status == 1L), male, bmi
  )]
  baseline[, bmi_z := as.numeric(scale(bmi))]
  stopifnot(!anyDuplicated(baseline$id), !anyNA(baseline),
            all(baseline$entry_age < baseline$exit_age))

  labs <- as.data.table(survival::nafld2)
  if (!include_preindex_labs) labs <- labs[days >= 0]
  labs <- labs[test %in% c("hdl", "chol")]
  labs <- labs[, .(value = mean(value)), by = .(id, days, test)]
  labs <- dcast(labs, id + days ~ test, value.var = "value")
  labs <- labs[complete.cases(labs)]
  data <- merge(labs, baseline, by = "id", all = FALSE)
  data[, visit_age := entry_age + days / 365.25]
  data <- data[visit_age < exit_age]
  data[, hdl_z := as.numeric(scale(hdl))]
  stopifnot(!anyDuplicated(data[, .(id, visit_age)]), !anyNA(data),
            all(data$visit_age >= data$entry_age))
  list(long = as.data.frame(data), baseline = as.data.frame(baseline))
}
