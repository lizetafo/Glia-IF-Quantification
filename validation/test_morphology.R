# ================================================================
# Validation of the skeleton measurements in both apps
#
# Loads ONLY the analysis functions from each app.R (the app is not
# launched), runs them on images with known answers, and prints the
# results next to the true values.
#
# Run from the repository folder in RStudio:
#   source("validation/test_morphology.R")
# Optional: also test the pre-fix versions with
#   TEST_ORIGINAL <- TRUE; source("validation/test_morphology.R")
# ================================================================

suppressPackageStartupMessages({
  library(EBImage)
  library(png)
})

load_app_functions <- function(app_file) {
  env <- new.env()
  for (expr in parse(app_file)) {
    is_function <- is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))
    is_infix <- is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.name(expr[[2]]) && grepl("^%", as.character(expr[[2]]))
    if (is_function || is_infix) eval(expr, env)
  }
  env
}

read_mask <- function(path, threshold = 25) {
  a <- png::readPNG(path)
  red <- if (length(dim(a)) == 3) a[, , 1] else a
  red * 255 > threshold
}

synthetic <- file.path("example_data", "synthetic")
dense <- file.path("validation", "dense_fields")

tests <- data.frame(
  image = c(file.path(synthetic, c("1_topology_shapes.png", "2_microglia_phantom.png",
                                   "4_microglia_realistic.png", "5_crossings.png")),
            file.path(dense, sprintf("dense_field_%d.png", 1:4))),
  true_endpoints = c(20, 32, 32, 16, 178, 178, 213, 248),
  true_branchpoints = c(5, 13, 13, 4, NA, NA, NA, NA),
  true_length_px = c(1134.3, NA, NA, NA, NA, NA, NA, NA),
  fiji_endpoints = c(20, NA, NA, 16, 242, 268, 318, 423),
  fiji_branchpoints = c(5, NA, NA, 5, 244, 331, 429, 600),
  stringsAsFactors = FALSE
)

apps <- c(GFAP = "apps/GFAP_IF_Quantifier/app.R", IBA1 = "apps/IBA1_IF_Quantifier/app.R")
if (exists("TEST_ORIGINAL") && isTRUE(TEST_ORIGINAL)) {
  apps <- c(apps, GFAP_before_fix = "validation/original_apps/GFAP_app_before_fix.R",
            IBA1_before_fix = "validation/original_apps/IBA1_app_before_fix.R")
}

results <- list()
for (app in names(apps)) {
  env <- load_app_functions(apps[[app]])
  analyze <- if (exists("analyze_iba1_skeleton", env, inherits = FALSE)) env$analyze_iba1_skeleton else env$analyze_gfap_skeleton
  for (i in seq_len(nrow(tests))) {
    m <- analyze(read_mask(tests$image[i]))$metrics
    get <- function(pattern) m[[grep(pattern, names(m))[1]]]
    results[[length(results) + 1]] <- data.frame(
      app = app, image = basename(tests$image[i]),
      endpoints = get("endpoints_rough$"), true_endpoints = tests$true_endpoints[i],
      branchpoints = get("branchpoints_rough$"), true_branchpoints = tests$true_branchpoints[i],
      length_px = round(get("length_px_rough$"), 1), true_length_px = tests$true_length_px[i],
      fiji_endpoints = tests$fiji_endpoints[i], fiji_branchpoints = tests$fiji_branchpoints[i])
  }
}
results <- do.call(rbind, results)
print(results, row.names = FALSE)
utils::write.csv(results, file.path("validation", "morphology_results.csv"), row.names = FALSE)

exact <- subset(results, !is.na(true_branchpoints) & !grepl("before_fix", app))
cat(sprintf("\nKnown-answer images: endpoints exact %d/%d, branch points exact %d/%d\n",
            sum(exact$endpoints == exact$true_endpoints), nrow(exact),
            sum(exact$branchpoints == exact$true_branchpoints), nrow(exact)))
cat("Dense fields: true_endpoints are estimated free process tips (see validation/dense_fields/estimate_true_tips.py).\n")
