# ================================================================
# GFAP IF QUANTIFIER
#
# Supports:
#   1. Nuance .im3 cubes with a reference-spectrum CSV (NNLS)
#   2. Nuance .im3 cubes without a library (deterministic NMF)
#   3. Composite TIFF/JPEG/PNG images
#
# For library-free NMF, fit components on one representative image,
# visually assign GFAP/DAPI/autofluorescence, lock the settings, and
# reuse the fitted spectra for the entire batch.
#
# Batch output can include:
#   - 16-bit TIFF channels
#   - PNG/JPEG previews
#   - Binary GFAP masks
#   - Colored overlays
#   - Skeleton-based GFAP morphology
#   - Scale-bar composites
#   - Quantification CSV
#   - ZIP archive containing the results
# ================================================================


# ================================================================
# WINDOWS JAVA SETUP
#
# RBioFormats uses Java. This first honors an existing JAVA_HOME,
# then looks in common Windows JDK folders. The first candidate is
# the portable JDK location already used on the current computer.
# ================================================================

configure_java_for_r <- function() {
  current_java_home <- Sys.getenv("JAVA_HOME", unset = "")

  if (nzchar(current_java_home) &&
      file.exists(file.path(current_java_home, "bin", "java.exe"))) {
    java_bin <- file.path(current_java_home, "bin")
    Sys.setenv(
      PATH = paste(java_bin, Sys.getenv("PATH"), sep = .Platform$path.sep)
    )
    return(invisible(current_java_home))
  }

  if (!identical(.Platform$OS.type, "windows")) {
    return(invisible(current_java_home))
  }

  user_profile <- Sys.getenv("USERPROFILE", unset = "")
  candidates <- c(
    file.path(
      user_profile,
      "Downloads",
      "jdk-26_windows-x64_bin",
      "jdk-26.0.2.1"
    ),
    Sys.glob(file.path(user_profile, "Downloads", "jdk-*", "jdk-*")),
    Sys.glob("C:/Program Files/Java/jdk-*"),
    Sys.glob("C:/Program Files/Java/jre*"),
    Sys.glob("C:/Program Files/Eclipse Adoptium/jdk-*")
  )

  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[
    vapply(
      candidates,
      function(path) file.exists(file.path(path, "bin", "java.exe")),
      logical(1)
    )
  ]

  if (!length(candidates)) {
    return(invisible(""))
  }

  selected_java_home <- normalizePath(
    candidates[1],
    winslash = "/",
    mustWork = TRUE
  )

  Sys.setenv(JAVA_HOME = selected_java_home)
  Sys.setenv(
    PATH = paste(
      file.path(selected_java_home, "bin"),
      Sys.getenv("PATH"),
      sep = .Platform$path.sep
    )
  )

  invisible(selected_java_home)
}

configured_java_home <- configure_java_for_r()

if (identical(.Platform$OS.type, "windows") &&
    !nzchar(configured_java_home) &&
    !nzchar(Sys.which("java"))) {
  stop(
    paste(
      "A Windows Java JDK was not found.",
      "Set JAVA_HOME or place the portable JDK in its expected Downloads folder,",
      "then restart the app."
    ),
    call. = FALSE
  )
}


# ================================================================
# REQUIRED PACKAGES
# ================================================================

required_cran <- c(
  "shiny",
  "DT",
  "ggplot2",
  "zip",
  "rJava"
)

required_bioc <- c(
  "EBImage",
  "RBioFormats"
)


# ================================================================
# CHECK INSTALLED PACKAGES
# ================================================================

missing_cran <- required_cran[
  !vapply(
    required_cran,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

missing_bioc <- required_bioc[
  !vapply(
    required_bioc,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

missing_packages <- c(
  missing_cran,
  missing_bioc
)


# ================================================================
# STOP AND SHOW INSTALLATION INSTRUCTIONS IF PACKAGES ARE MISSING
# ================================================================

if (length(missing_packages) > 0) {
  
  stop(
    paste0(
      "Missing package(s): ",
      paste(missing_packages, collapse = ", "),
      "\n\n",
      
      "Run the following commands in the RStudio Console:\n\n",
      
      "install.packages(",
      "\"BiocManager\"",
      ")\n\n",
      
      "install.packages(",
      "c(\"shiny\", \"DT\", \"ggplot2\", \"zip\", \"rJava\")",
      ")\n\n",
      
      "BiocManager::install(",
      "c(\"EBImage\", \"RBioFormats\"), ",
      "ask = FALSE, update = FALSE",
      ")\n\n",
      
      "After installation finishes, restart R and run the app again."
    ),
    call. = FALSE
  )
}


# ================================================================
# LOAD PACKAGES
# ================================================================

suppressPackageStartupMessages({
  
  library(shiny)
  
  library(DT)
  
  library(ggplot2)
  
  library(EBImage)

  library(rJava)
  
  library(RBioFormats)
  
})


# ================================================================
# SHINY UPLOAD LIMIT
# ================================================================

options(
  shiny.maxRequestSize = 2048 * 1024^2
)


# ================================================================
# GENERAL HELPERS
# ================================================================

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

normalize01 <- function(x, q_low = 0.001, q_high = 0.995) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  limits <- as.numeric(stats::quantile(x, c(q_low, q_high), na.rm = TRUE, names = FALSE))
  if (!all(is.finite(limits)) || limits[2] <= limits[1]) return(matrix(0, nrow(x), ncol(x)))
  x <- (x - limits[1]) / (limits[2] - limits[1])
  pmin(pmax(x, 0), 1)
}

draw_gray <- function(image, main = "") {
  z <- normalize01(image)
  palette <- grDevices::gray.colors(256, start = 0, end = 1)
  raster <- matrix(palette[1 + floor(z * 255)], nrow = nrow(z), ncol = ncol(z))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, ncol(z)), ylim = c(nrow(z), 0), asp = 1)
  graphics::rasterImage(raster, 0, nrow(z), ncol(z), 0)
  graphics::title(main = main)
  graphics::box()
}

plot_gray <- function(image, main = "") {
  old <- graphics::par(mar = c(1, 1, 3, 1))
  on.exit(graphics::par(old), add = TRUE)
  draw_gray(image, main)
}

plot_rgb_matrix <- function(red, green, blue, main = "") {
  rgb_image <- grDevices::rgb(normalize01(red), normalize01(green), normalize01(blue))
  dim(rgb_image) <- dim(red)
  old <- graphics::par(mar = c(1, 1, 3, 1))
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, ncol(red)), ylim = c(nrow(red), 0), asp = 1)
  graphics::rasterImage(rgb_image, 0, nrow(red), ncol(red), 0)
  graphics::title(main = main)
  graphics::box()
}

plot_mask_overlay <- function(gfap, dapi = NULL, mask, main = "") {
  red <- normalize01(gfap)
  blue <- if (is.null(dapi)) matrix(0, nrow(gfap), ncol(gfap)) else normalize01(dapi)
  green <- matrix(0, nrow(gfap), ncol(gfap))
  positive <- as.logical(mask)
  red[positive] <- 0.15 * red[positive]
  blue[positive] <- 0.15 * blue[positive]
  green[positive] <- 1
  rgb_image <- grDevices::rgb(red, green, blue)
  dim(rgb_image) <- dim(gfap)
  old <- graphics::par(mar = c(1, 1, 3, 1))
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, ncol(gfap)), ylim = c(nrow(gfap), 0), asp = 1)
  graphics::rasterImage(rgb_image, 0, nrow(gfap), ncol(gfap), 0)
  graphics::title(main = main)
  graphics::box()
}

make_skeleton_overlay_array <- function(gfap, morphology) {
  background <- 0.30 * normalize01(gfap)
  red <- background
  green <- background
  blue <- background

  skeleton <- morphology$skeleton
  branchpoints <- morphology$branchpoint_pixels
  endpoints <- morphology$endpoint_pixels

  red[skeleton] <- 1
  green[skeleton] <- 1
  blue[skeleton] <- 0

  red[branchpoints] <- 1
  green[branchpoints] <- 0
  blue[branchpoints] <- 0

  red[endpoints] <- 0
  green[endpoints] <- 1
  blue[endpoints] <- 1

  array(c(red, green, blue), dim = c(nrow(gfap), ncol(gfap), 3L))
}

plot_skeleton_overlay <- function(gfap, morphology,
                                  main = "Skeleton: yellow; branches: red; endpoints: cyan") {
  overlay <- make_skeleton_overlay_array(gfap, morphology)
  rgb_image <- grDevices::rgb(overlay[, , 1], overlay[, , 2], overlay[, , 3])
  dim(rgb_image) <- dim(gfap)
  old <- graphics::par(mar = c(1, 1, 3, 1))
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, ncol(gfap)), ylim = c(nrow(gfap), 0), asp = 1)
  graphics::rasterImage(rgb_image, 0, nrow(gfap), ncol(gfap), 0)
  graphics::title(main = main)
  graphics::box()
}

# ================================================================
# IMAGE EXPORT HELPERS
# ================================================================

safe_file_stem <- function(filename) {
  stem <- tools::file_path_sans_ext(basename(filename))
  stem <- gsub("[^A-Za-z0-9_-]+", "_", stem)
  stem <- gsub("^_+|_+$", "", stem)
  if (!nzchar(stem)) stem <- "image"
  stem
}

image_export_scale <- function(image, probability = 0.999) {
  value <- as.numeric(stats::quantile(image, probability, na.rm = TRUE, names = FALSE))
  if (!is.finite(value) || value <= 0) value <- max(image, na.rm = TRUE)
  if (!is.finite(value) || value <= 0) value <- 1
  value
}

scale_to_unit <- function(image, maximum) {
  image <- as.matrix(image)
  image[!is.finite(image)] <- 0
  pmin(pmax(image / maximum, 0), 1)
}

write_gray_export <- function(image, path, format = "tiff", fixed_scale = NULL,
                              jpeg_quality = 95) {
  if (identical(format, "tiff")) {
    if (is.null(fixed_scale) || !is.finite(fixed_scale) || fixed_scale <= 0) {
      stop("A positive fixed scale is required for TIFF export.")
    }
    prepared <- scale_to_unit(image, fixed_scale)
    EBImage::writeImage(
      EBImage::Image(t(prepared)), path, type = "tiff",
      bits.per.sample = 16, compression = "none"
    )
  } else {
    prepared <- normalize01(image)
    EBImage::writeImage(
      EBImage::Image(t(prepared)), path, type = format,
      quality = jpeg_quality
    )
  }
  invisible(path)
}

write_binary_mask_png <- function(mask, path) {
  prepared <- matrix(as.numeric(mask > 0), nrow(mask), ncol(mask))
  EBImage::writeImage(EBImage::Image(t(prepared)), path, type = "png")
  invisible(path)
}

write_skeleton_overlay_png <- function(gfap, morphology, path) {
  overlay <- make_skeleton_overlay_array(gfap, morphology)
  overlay_xyc <- aperm(overlay, c(2, 1, 3))
  EBImage::writeImage(
    EBImage::Image(overlay_xyc, colormode = "Color"),
    path, type = "png"
  )
  invisible(path)
}

make_overlay_array <- function(gfap, dapi = NULL, mask) {
  red <- normalize01(gfap)
  blue <- if (is.null(dapi)) matrix(0, nrow(gfap), ncol(gfap)) else normalize01(dapi)
  green <- matrix(0, nrow(gfap), ncol(gfap))
  positive <- as.logical(mask)
  red[positive] <- 0.15 * red[positive]
  blue[positive] <- 0.15 * blue[positive]
  green[positive] <- 1
  array(c(red, green, blue), dim = c(nrow(gfap), ncol(gfap), 3L))
}

write_overlay_png <- function(gfap, dapi, mask, path) {
  overlay <- make_overlay_array(gfap, dapi, mask)
  overlay_xyc <- aperm(overlay, c(2, 1, 3))
  EBImage::writeImage(
    EBImage::Image(overlay_xyc, colormode = "Color"),
    path, type = "png"
  )
  invisible(path)
}

scale_bar_spec <- function(image_width, image_height, pixel_size_um, scale_bar_um) {
  if (!is.finite(pixel_size_um) || pixel_size_um <= 0) {
    return(list(valid = FALSE, status = "Not drawn: calibrated pixel size is unavailable"))
  }
  if (!is.finite(scale_bar_um) || scale_bar_um <= 0) {
    return(list(valid = FALSE, status = "Not drawn: scale-bar length must be positive"))
  }

  bar_pixels <- as.integer(round(scale_bar_um / pixel_size_um))
  margin <- max(12L, as.integer(round(min(image_width, image_height) * 0.035)))
  if (bar_pixels < 1L || bar_pixels > image_width - 2L * margin) {
    return(list(valid = FALSE, status = "Not drawn: requested scale bar does not fit the image"))
  }

  list(
    valid = TRUE,
    pixels = bar_pixels,
    margin = margin,
    status = sprintf(
      "Drawn: %g um (%d pixels at %.8g um/pixel)",
      scale_bar_um, bar_pixels, pixel_size_um
    )
  )
}

write_clean_composite_png <- function(gfap, dapi, path,
                                      pixel_size_um = NA_real_,
                                      scale_bar_um = NA_real_,
                                      scale_bar_color = "white",
                                      scale_bar_thickness = 6) {
  empty_mask <- matrix(FALSE, nrow(gfap), ncol(gfap))
  composite <- make_overlay_array(gfap, dapi, empty_mask)
  spec <- scale_bar_spec(ncol(gfap), nrow(gfap), pixel_size_um, scale_bar_um)

  if (!isTRUE(spec$valid)) {
    composite_xyc <- aperm(composite, c(2, 1, 3))
    EBImage::writeImage(
      EBImage::Image(composite_xyc, colormode = "Color"),
      path, type = "png"
    )
    return(invisible(path))
  }

  scale_bar_color <- if (identical(scale_bar_color, "black")) "black" else "white"
  scale_bar_thickness <- max(2, as.numeric(scale_bar_thickness %||% 6))
  image_height <- nrow(gfap)
  image_width <- ncol(gfap)
  raster <- matrix(
    grDevices::rgb(composite[, , 1], composite[, , 2], composite[, , 3]),
    nrow = image_height,
    ncol = image_width
  )

  grDevices::png(
    filename = path,
    width = image_width,
    height = image_height,
    units = "px",
    res = 96,
    bg = "black"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, image_width), ylim = c(image_height, 0), asp = 1)
  graphics::rasterImage(raster, 0, image_height, image_width, 0, interpolate = FALSE)

  x_right <- image_width - spec$margin
  x_left <- x_right - spec$pixels
  bar_y <- image_height - spec$margin
  label <- paste0(format(scale_bar_um, trim = TRUE, scientific = FALSE), " µm")
  label_cex <- max(0.75, min(1.35, image_width / 650))
  label_y <- bar_y - scale_bar_thickness / 2 - spec$margin * 0.55
  label_width <- graphics::strwidth(label, cex = label_cex, units = "user")
  label_height <- graphics::strheight(label, cex = label_cex, units = "user")
  center_x <- (x_left + x_right) / 2
  box_left <- max(0, min(x_left - 8, center_x - label_width / 2 - 6))
  box_right <- min(image_width, max(x_right + 8, center_x + label_width / 2 + 6))
  box_top <- max(0, label_y - label_height / 2 - 5)
  box_bottom <- min(image_height, bar_y + scale_bar_thickness / 2 + 7)
  box_color <- if (identical(scale_bar_color, "white")) {
    grDevices::adjustcolor("black", alpha.f = 0.6)
  } else {
    grDevices::adjustcolor("white", alpha.f = 0.7)
  }
  graphics::rect(box_left, box_top, box_right, box_bottom, col = box_color, border = NA)
  graphics::rect(
    x_left, bar_y - scale_bar_thickness / 2,
    x_right, bar_y + scale_bar_thickness / 2,
    col = scale_bar_color, border = NA
  )
  graphics::text(
    center_x, label_y, labels = label,
    col = scale_bar_color, cex = label_cex, font = 2
  )
  invisible(path)
}

channel_export_extension <- function(format) {
  switch(format, tiff = "tif", png = "png", jpeg = "jpg", "tif")
}

stage_uploaded_files <- function(upload_table) {
  extensions <- tolower(tools::file_ext(upload_table$name))
  display_names <- make.unique(upload_table$name)
  staged_paths <- character(nrow(upload_table))
  for (i in seq_len(nrow(upload_table))) {
    staged_paths[i] <- tempfile(pattern = "gfap_", fileext = paste0(".", extensions[i]))
    if (!file.copy(upload_table$datapath[i], staged_paths[i], overwrite = TRUE)) {
      stop("Could not stage uploaded file: ", upload_table$name[i])
    }
  }
  data.frame(
    name = display_names,
    original_name = upload_table$name,
    path = staged_paths,
    ext = extensions,
    stringsAsFactors = FALSE
  )
}

# ================================================================
# IMAGE READERS
# ================================================================

extract_im3_pixel_size_um <- function(path) {
  file_size <- file.info(path)$size
  if (!length(file_size) || !is.finite(file_size) || file_size <= 0 || file_size > .Machine$integer.max) {
    return(NA_real_)
  }

  bytes <- readBin(path, what = "raw", n = as.integer(file_size))
  record_name <- "MillimetersPerPixel"
  positions <- grepRaw(record_name, bytes, fixed = TRUE, all = TRUE)
  if (!length(positions)) return(NA_real_)

  for (position in positions) {
    if (position <= 4L) next
    name_length <- readBin(
      bytes[(position - 4L):(position - 1L)],
      what = integer(), n = 1L, size = 4L, signed = TRUE, endian = "little"
    )
    if (!identical(as.integer(name_length), nchar(record_name, type = "bytes"))) next

    after_name <- position + nchar(record_name, type = "bytes")
    if (after_name + 19L > length(bytes)) next
    record_size <- readBin(
      bytes[after_name:(after_name + 3L)],
      what = integer(), n = 1L, size = 4L, signed = TRUE, endian = "little"
    )
    record_type <- readBin(
      bytes[(after_name + 4L):(after_name + 7L)],
      what = integer(), n = 1L, size = 4L, signed = TRUE, endian = "little"
    )
    if (!is.finite(record_size) || record_size < 20L || record_type != 7L) next

    payload_start <- after_name + 8L
    value_starts <- c(payload_start + 8L, payload_start + 4L)
    for (value_start in value_starts) {
      if (value_start + 3L > length(bytes)) next
      millimeters_per_pixel <- readBin(
        bytes[value_start:(value_start + 3L)],
        what = numeric(), n = 1L, size = 4L, endian = "little"
      )
      if (is.finite(millimeters_per_pixel) &&
          millimeters_per_pixel > 0 && millimeters_per_pixel < 10) {
        return(as.numeric(millimeters_per_pixel) * 1000)
      }
    }
  }
  NA_real_
}

as_yxc_cube <- function(image) {
  arr <- as.array(image)
  d <- dim(arr)
  order_names <- tryCatch(tolower(as.character(RBioFormats::dimorder(image))),
                          error = function(e) character())
  if (length(order_names) != length(d)) {
    order_names <- c("x", "y", rep("other", max(0, length(d) - 2)))
  }
  ix <- match("x", order_names)
  iy <- match("y", order_names)
  if (is.na(ix) || is.na(iy)) stop("Bio-Formats did not report X and Y dimensions.")
  remaining <- setdiff(seq_along(d), c(ix, iy))
  permuted <- aperm(arr, c(iy, ix, remaining))
  dim(permuted) <- c(d[iy], d[ix], if (length(remaining)) prod(d[remaining]) else 1L)
  storage.mode(permuted) <- "double"
  permuted[!is.finite(permuted)] <- 0
  pmax(permuted, 0)
}

read_im3_cube <- function(path) {
  image <- RBioFormats::read.image(path, normalize = FALSE)
  cube <- as_yxc_cube(image)
  if (dim(cube)[3] < 2) stop("The IM3 file did not expose multiple spectral planes.")
  cube
}

read_composite_image <- function(path) {
  image <- EBImage::readImage(path)
  arr <- as.array(image)
  if (length(dim(arr)) == 2) arr <- array(arr, dim = c(dim(arr), 1L))
  if (length(dim(arr)) > 3) dim(arr) <- c(dim(arr)[1:2], prod(dim(arr)[-c(1, 2)]))
  arr <- aperm(arr, c(2, 1, 3))
  storage.mode(arr) <- "double"
  arr[!is.finite(arr)] <- 0
  pmax(arr, 0)
}

composite_components <- function(image) {
  n_channels <- dim(image)[3]
  output <- list()
  if (n_channels >= 1) output$Red <- image[, , 1]
  if (n_channels >= 2) output$Green <- image[, , 2]
  if (n_channels >= 3) output$Blue <- image[, , 3]
  if (n_channels == 1) output <- list(Gray = image[, , 1])
  output
}

# ================================================================
# REFERENCE CSV
#
# Expected layout: one spectral plane per row and one component per
# column. An optional first column may contain wavelength/plane labels.
# ================================================================

read_basis_csv <- function(path, channels) {
  data <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(data) < 1 || nrow(data) < 2) stop("The reference CSV is empty or malformed.")
  first_name <- tolower(names(data)[1])
  first_numeric <- suppressWarnings(as.numeric(data[[1]]))
  drop_first <- first_name %in% c("wavelength", "plane", "channel", "index") ||
    (ncol(data) > 1 && any(!is.finite(first_numeric)))
  if (drop_first) data <- data[-1]
  basis <- as.matrix(data)
  storage.mode(basis) <- "double"
  if (nrow(basis) != channels) {
    stop("Reference CSV has ", nrow(basis), " planes but the image has ", channels, ".")
  }
  basis[!is.finite(basis)] <- 0
  basis <- pmax(basis, 0)
  if (any(colSums(basis) <= 0)) stop("Every reference spectrum must contain a positive value.")
  basis <- sweep(basis, 2, sqrt(colSums(basis^2)), "/")
  colnames(basis) <- make.unique(names(data) %||% paste0("Reference ", seq_len(ncol(basis))))
  basis
}

basis_to_df <- function(basis) {
  if (is.null(basis) || !length(basis)) return(data.frame())
  data.frame(
    Plane = rep(seq_len(nrow(basis)), times = ncol(basis)),
    Magnitude = as.vector(basis),
    Component = rep(colnames(basis), each = nrow(basis)),
    stringsAsFactors = FALSE
  )
}

# ================================================================
# LIBRARY-FREE NMF
#
# The NMF spectra are fitted on a reproducible sample of nonzero
# pixels. The fitted spectra are then held fixed while all pixels are
# projected with non-negative least squares.
# ================================================================

sample_cube_pixels <- function(cube, max_pixels, seed) {
  n_pixels <- dim(cube)[1] * dim(cube)[2]
  channels <- dim(cube)[3]
  total <- numeric(n_pixels)
  for (channel in seq_len(channels)) total <- total + as.vector(cube[, , channel])
  eligible <- which(is.finite(total) & total > 0)
  if (!length(eligible)) stop("The image cube contains no positive pixels.")
  set.seed(as.integer(seed))
  if (length(eligible) > max_pixels) eligible <- sample(eligible, max_pixels)
  pixels <- matrix(0, nrow = length(eligible), ncol = channels)
  for (channel in seq_len(channels)) pixels[, channel] <- as.vector(cube[, , channel])[eligible]
  pixels[!is.finite(pixels)] <- 0
  scale_value <- as.numeric(stats::quantile(pixels, 0.999, na.rm = TRUE, names = FALSE))
  if (!is.finite(scale_value) || scale_value <= 0) scale_value <- max(pixels)
  pixels / max(scale_value, .Machine$double.eps)
}

fit_nmf_basis <- function(cube, rank = 3, max_iter = 250,
                          sample_pixels = 50000, seed = 1,
                          tolerance = 1e-5) {
  channels <- dim(cube)[3]
  if (rank > channels) stop("NMF components cannot exceed the number of spectral planes.")
  Y <- sample_cube_pixels(cube, as.integer(sample_pixels), seed)
  eps <- 1e-10
  set.seed(as.integer(seed))
  W <- matrix(stats::runif(nrow(Y) * rank, 0.1, 1), nrow(Y), rank)
  H <- matrix(stats::runif(channels * rank, 0.1, 1), channels, rank)
  converged_at <- max_iter

  for (iteration in seq_len(max_iter)) {
    old_H <- H
    W <- W * ((Y %*% H) / pmax(W %*% crossprod(H), eps))
    H <- H * ((crossprod(Y, W)) / pmax(H %*% crossprod(W), eps))
    component_scale <- pmax(sqrt(colSums(H^2)), eps)
    H <- sweep(H, 2, component_scale, "/")
    W <- sweep(W, 2, component_scale, "*")
    change <- max(abs(H - old_H)) / max(1, max(abs(old_H)))
    if (is.finite(change) && change < tolerance) {
      converged_at <- iteration
      break
    }
  }

  colnames(H) <- paste0("Component ", seq_len(rank))
  reconstruction_error <- sqrt(sum((Y - W %*% t(H))^2) / sum(Y^2))
  list(
    basis = H,
    iterations = converged_at,
    sampled_pixels = nrow(Y),
    relative_error = reconstruction_error
  )
}

# ================================================================
# FIXED-BASIS NON-NEGATIVE PROJECTION
# ================================================================

project_nnls_matrix <- function(Y, basis, max_iter = 100, tolerance = 1e-6) {
  G <- crossprod(basis)
  H <- Y %*% basis
  ridge <- diag(1e-10, ncol(basis))
  B <- pmax(H %*% solve(G + ridge), 0)
  largest_eigenvalue <- max(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(largest_eigenvalue) || largest_eigenvalue <= 0) stop("Invalid component basis.")
  step_size <- 0.95 / largest_eigenvalue
  for (iteration in seq_len(max_iter)) {
    B_new <- pmax(B - step_size * (B %*% G - H), 0)
    relative_change <- max(abs(B_new - B)) / max(1, max(abs(B)))
    B <- B_new
    if (is.finite(relative_change) && relative_change < tolerance) break
  }
  B
}

unmix_fixed_basis <- function(cube, basis, max_iter = 100, chunk_pixels = 200000) {
  height <- dim(cube)[1]
  width <- dim(cube)[2]
  channels <- dim(cube)[3]
  if (nrow(basis) != channels) stop("Basis spectra do not match the image channel count.")
  Y <- matrix(cube, ncol = channels)
  Y[!is.finite(Y)] <- 0
  Y <- pmax(Y, 0)
  B <- matrix(0, nrow(Y), ncol(basis))
  starts <- seq.int(1L, nrow(Y), by = as.integer(chunk_pixels))
  for (start in starts) {
    index <- start:min(nrow(Y), start + as.integer(chunk_pixels) - 1L)
    B[index, ] <- project_nnls_matrix(Y[index, , drop = FALSE], basis, max_iter)
  }
  output <- setNames(vector("list", ncol(B)), colnames(basis))
  for (component in seq_len(ncol(B))) {
    output[[component]] <- matrix(B[, component], nrow = height, ncol = width)
  }
  output
}

# ================================================================
# SEGMENTATION AND MEASUREMENTS
# ================================================================

background_correct <- function(image, sigma = 15) {
  image <- pmax(as.matrix(image), 0)
  if (!is.finite(sigma) || sigma <= 0) return(image)
  background <- EBImage::gblur(EBImage::Image(image), sigma = sigma)
  pmax(image - as.array(background), 0)
}

suggest_threshold <- function(image) {
  minimum <- min(image, na.rm = TRUE)
  maximum <- max(image, na.rm = TRUE)
  if (!is.finite(maximum) || maximum <= minimum) return(minimum)
  normalized <- pmin(pmax((image - minimum) / (maximum - minimum), 0), 1)
  threshold_normalized <- EBImage::otsu(EBImage::Image(normalized))
  as.numeric(minimum + threshold_normalized * (maximum - minimum))
}

remove_small_objects <- function(mask, min_pixels = 0) {
  mask <- as.matrix(mask > 0)
  if (min_pixels <= 1) return(mask)
  labels <- as.integer(as.array(EBImage::bwlabel(EBImage::Image(mask))))
  sizes <- tabulate(labels)
  keep <- which(sizes >= min_pixels)
  matrix(labels %in% keep, nrow = nrow(mask), ncol = ncol(mask))
}

make_gfap_mask <- function(gfap, threshold, min_pixels = 0) {
  remove_small_objects(gfap > threshold, min_pixels)
}

zhang_suen_skeleton <- function(mask, max_iter = 500L) {
  mask <- as.matrix(mask > 0)
  height <- nrow(mask)
  width <- ncol(mask)
  if (!height || !width || !any(mask)) return(mask)

  padded <- matrix(FALSE, nrow = height + 2L, ncol = width + 2L)
  rows <- 2L:(height + 1L)
  columns <- 2L:(width + 1L)
  padded[rows, columns] <- mask

  for (iteration in seq_len(as.integer(max_iter))) {
    deleted_this_iteration <- FALSE

    for (pass in 1:2) {
      center <- padded[rows, columns, drop = FALSE]
      north <- padded[rows - 1L, columns, drop = FALSE]
      northeast <- padded[rows - 1L, columns + 1L, drop = FALSE]
      east <- padded[rows, columns + 1L, drop = FALSE]
      southeast <- padded[rows + 1L, columns + 1L, drop = FALSE]
      south <- padded[rows + 1L, columns, drop = FALSE]
      southwest <- padded[rows + 1L, columns - 1L, drop = FALSE]
      west <- padded[rows, columns - 1L, drop = FALSE]
      northwest <- padded[rows - 1L, columns - 1L, drop = FALSE]

      neighbor_count <- north + northeast + east + southeast +
        south + southwest + west + northwest
      transitions <-
        (!north & northeast) + (!northeast & east) +
        (!east & southeast) + (!southeast & south) +
        (!south & southwest) + (!southwest & west) +
        (!west & northwest) + (!northwest & north)

      connectivity_rule <- if (pass == 1L) {
        !(north & east & south) & !(east & south & west)
      } else {
        !(north & east & west) & !(north & south & west)
      }

      remove <- center & neighbor_count >= 2L & neighbor_count <= 6L &
        transitions == 1L & connectivity_rule

      if (any(remove)) {
        center[remove] <- FALSE
        padded[rows, columns] <- center
        deleted_this_iteration <- TRUE
      }
    }

    if (!deleted_this_iteration) break
  }

  padded[rows, columns, drop = FALSE]
}

neighbor_degree8 <- function(binary) {
  binary <- as.matrix(binary > 0)
  height <- nrow(binary)
  width <- ncol(binary)
  padded <- matrix(FALSE, nrow = height + 2L, ncol = width + 2L)
  rows <- 2L:(height + 1L)
  columns <- 2L:(width + 1L)
  padded[rows, columns] <- binary

  padded[rows - 1L, columns, drop = FALSE] +
    padded[rows - 1L, columns + 1L, drop = FALSE] +
    padded[rows, columns + 1L, drop = FALSE] +
    padded[rows + 1L, columns + 1L, drop = FALSE] +
    padded[rows + 1L, columns, drop = FALSE] +
    padded[rows + 1L, columns - 1L, drop = FALSE] +
    padded[rows, columns - 1L, drop = FALSE] +
    padded[rows - 1L, columns - 1L, drop = FALSE]
}

count_components8 <- function(binary) {
  binary <- as.matrix(binary > 0)
  height <- nrow(binary)
  width <- ncol(binary)
  seeds <- which(binary)
  if (!length(seeds)) return(0L)

  visited <- matrix(FALSE, nrow = height, ncol = width)
  queue_rows <- integer(length(seeds))
  queue_columns <- integer(length(seeds))
  component_count <- 0L

  for (seed in seeds) {
    seed_row <- ((seed - 1L) %% height) + 1L
    seed_column <- ((seed - 1L) %/% height) + 1L
    if (visited[seed_row, seed_column]) next

    component_count <- component_count + 1L
    head <- 1L
    tail <- 1L
    queue_rows[1L] <- seed_row
    queue_columns[1L] <- seed_column
    visited[seed_row, seed_column] <- TRUE

    while (head <= tail) {
      current_row <- queue_rows[head]
      current_column <- queue_columns[head]
      head <- head + 1L

      for (row_offset in c(-1L, 0L, 1L)) {
        for (column_offset in c(-1L, 0L, 1L)) {
          if (row_offset == 0L && column_offset == 0L) next
          candidate_row <- current_row + row_offset
          candidate_column <- current_column + column_offset
          if (candidate_row < 1L || candidate_row > height ||
              candidate_column < 1L || candidate_column > width) next
          if (binary[candidate_row, candidate_column] &&
              !visited[candidate_row, candidate_column]) {
            tail <- tail + 1L
            queue_rows[tail] <- candidate_row
            queue_columns[tail] <- candidate_column
            visited[candidate_row, candidate_column] <- TRUE
          }
        }
      }
    }
  }

  component_count
}

# ----------------------------------------------------------------
# SKELETON CLEAN-UP AND MEASUREMENT RULES (validated 2026-11)
#
# Zhang-Suen thinning leaves "staircase" pixels on diagonal lines and
# extra pixels at bends. Each of these has 3+ neighbours and was being
# counted as a branch point, and its extra steps were added to the
# length. The functions below:
#   1. remove redundant skeleton pixels (not needed for connectivity),
#   2. count junction pixels within 2 px of each other as ONE branch point,
#   3. measure length without counting a diagonal step where an
#      orthogonal path already joins the same two pixels.
# Validated against images with known answers and Fiji AnalyzeSkeleton
# (see validation/VALIDATION.md).
# ----------------------------------------------------------------

skeleton_padded <- function(binary) {
  height <- nrow(binary)
  width <- ncol(binary)
  padded <- matrix(FALSE, nrow = height + 2L, ncol = width + 2L)
  padded[2L:(height + 1L), 2L:(width + 1L)] <- binary
  padded
}

# Yokoi 8-connectivity number for one pixel of a padded skeleton.
# 1 = removing the pixel keeps its neighbours connected.
yokoi_c8 <- function(padded, r, c) {
  # neighbours in order E, NE, N, NW, W, SW, S, SE (r = row, c = column)
  v <- c(padded[r, c + 1L], padded[r - 1L, c + 1L], padded[r - 1L, c],
         padded[r - 1L, c - 1L], padded[r, c - 1L], padded[r + 1L, c - 1L],
         padded[r + 1L, c], padded[r + 1L, c + 1L])
  x <- 1L - as.integer(v)
  x <- c(x, x[1L])
  total <- 0L
  for (k in c(1L, 3L, 5L, 7L)) {
    total <- total + x[k] - x[k] * x[k + 1L] * x[k + 2L]
  }
  total
}

# Remove skeleton pixels that are not needed to keep the skeleton connected
# (inside corners of bends and staircases). Endpoints are never removed.
thin_redundant_pixels <- function(skeleton) {
  skeleton <- as.matrix(skeleton > 0)
  if (!any(skeleton)) return(skeleton)
  padded <- skeleton_padded(skeleton)
  height <- nrow(skeleton)
  width <- ncol(skeleton)

  repeat {
    rows <- 2L:(height + 1L)
    columns <- 2L:(width + 1L)
    n <- padded[rows - 1L, columns] + padded[rows + 1L, columns] +
      padded[rows, columns - 1L] + padded[rows, columns + 1L] +
      padded[rows - 1L, columns - 1L] + padded[rows - 1L, columns + 1L] +
      padded[rows + 1L, columns - 1L] + padded[rows + 1L, columns + 1L]
    east <- padded[rows, columns + 1L]
    west <- padded[rows, columns - 1L]
    north <- padded[rows - 1L, columns]
    south <- padded[rows + 1L, columns]
    orthogonal_pair <- (east & north) | (north & west) | (west & south) | (south & east)
    candidates <- which(padded[rows, columns] & n >= 2L & orthogonal_pair, arr.ind = TRUE)
    if (!nrow(candidates)) break

    removed <- 0L
    for (k in seq_len(nrow(candidates))) {
      r <- candidates[k, 1L] + 1L
      c <- candidates[k, 2L] + 1L
      if (!padded[r, c]) next
      neighbours <- sum(padded[(r - 1L):(r + 1L), (c - 1L):(c + 1L)]) - 1L
      if (neighbours < 2L) next
      e <- padded[r, c + 1L]; w <- padded[r, c - 1L]
      nn <- padded[r - 1L, c]; s <- padded[r + 1L, c]
      if (!((e && nn) || (nn && w) || (w && s) || (s && e))) next
      if (yokoi_c8(padded, r, c) != 1L) next
      padded[r, c] <- FALSE
      removed <- removed + 1L
    }
    if (removed == 0L) break
  }
  padded[2L:(height + 1L), 2L:(width + 1L), drop = FALSE]
}

# Count groups of pixels, joining pixels up to `radius` px apart
# (chessboard distance). radius = 1 is ordinary 8-connectivity.
count_components_within <- function(binary, radius = 2L) {
  binary <- as.matrix(binary > 0)
  height <- nrow(binary)
  width <- ncol(binary)
  seeds <- which(binary)
  if (!length(seeds)) return(0L)
  visited <- matrix(FALSE, nrow = height, ncol = width)
  queue_rows <- integer(length(seeds))
  queue_columns <- integer(length(seeds))
  offsets <- expand.grid(dr = -radius:radius, dc = -radius:radius)
  offsets <- offsets[!(offsets$dr == 0L & offsets$dc == 0L), ]
  count <- 0L
  for (seed in seeds) {
    seed_row <- ((seed - 1L) %% height) + 1L
    seed_column <- ((seed - 1L) %/% height) + 1L
    if (visited[seed_row, seed_column]) next
    count <- count + 1L
    head <- 1L; tail <- 1L
    queue_rows[1L] <- seed_row; queue_columns[1L] <- seed_column
    visited[seed_row, seed_column] <- TRUE
    while (head <= tail) {
      cr <- queue_rows[head] + offsets$dr
      cc <- queue_columns[head] + offsets$dc
      head <- head + 1L
      ok <- cr >= 1L & cr <= height & cc >= 1L & cc <= width
      cr <- cr[ok]; cc <- cc[ok]
      hit <- binary[cbind(cr, cc)] & !visited[cbind(cr, cc)]
      if (any(hit)) {
        cr <- cr[hit]; cc <- cc[hit]
        visited[cbind(cr, cc)] <- TRUE
        idx <- seq_along(cr) + tail
        queue_rows[idx] <- cr; queue_columns[idx] <- cc
        tail <- tail + length(cr)
      }
    }
  }
  count
}

# Weighted skeleton length: orthogonal steps = 1, diagonal steps = sqrt(2).
# A diagonal step is not counted when an orthogonal path already joins the
# two pixels (otherwise a bend is measured twice). Isolated single pixels
# have no length.
skeleton_edge_length_pixels <- function(skeleton, degree = NULL) {
  skeleton <- as.matrix(skeleton > 0)
  height <- nrow(skeleton)
  width <- ncol(skeleton)
  if (height < 2L || width < 2L) {
    return(as.numeric(sum(skeleton[-1L] & skeleton[-length(skeleton)])))
  }
  top <- seq_len(height - 1L); bottom <- 2L:height
  left <- seq_len(width - 1L); right <- 2L:width
  orthogonal <- sum(skeleton[, left] & skeleton[, right]) +
    sum(skeleton[top, ] & skeleton[bottom, ])
  # down-right diagonal (r, c) -> (r + 1, c + 1)
  down_right <- skeleton[top, left] & skeleton[bottom, right] &
    !skeleton[top, right] & !skeleton[bottom, left]
  # down-left diagonal (r, c + 1) -> (r + 1, c)
  down_left <- skeleton[top, right] & skeleton[bottom, left] &
    !skeleton[top, left] & !skeleton[bottom, right]
  as.numeric(orthogonal + sqrt(2) * (sum(down_right) + sum(down_left)))
}

analyze_gfap_skeleton <- function(mask, pixel_size_um = NA_real_) {
  mask <- as.matrix(mask > 0)
  skeleton <- zhang_suen_skeleton(mask)
  skeleton <- thin_redundant_pixels(skeleton)
  degree <- neighbor_degree8(skeleton)
  endpoint_pixels <- skeleton & degree == 1L
  branchpoint_pixels <- skeleton & degree >= 3L

  # Skeletons clipped by an image border create artificial endpoints.
  if (nrow(mask) > 0L && ncol(mask) > 0L) {
    endpoint_pixels[c(1L, nrow(mask)), ] <- FALSE
    endpoint_pixels[, c(1L, ncol(mask))] <- FALSE
    branchpoint_pixels[c(1L, nrow(mask)), ] <- FALSE
    branchpoint_pixels[, c(1L, ncol(mask))] <- FALSE
  }

  length_pixels <- skeleton_edge_length_pixels(skeleton, degree)
  endpoint_count <- sum(endpoint_pixels)
  # junction pixels within 2 px of each other form one branch point
  branchpoint_count <- count_components_within(branchpoint_pixels, 2L)
  network_fragments <- count_components8(skeleton)

  calibrated <- is.finite(pixel_size_um) && pixel_size_um > 0
  fov_area_um2 <- if (calibrated) length(mask) * pixel_size_um^2 else NA_real_
  fov_area_mm2 <- if (calibrated) fov_area_um2 / 1e6 else NA_real_
  positive_area_um2 <- if (calibrated) sum(mask) * pixel_size_um^2 else NA_real_
  skeleton_length_um <- if (calibrated) length_pixels * pixel_size_um else NA_real_

  metrics <- data.frame(
    GFAP_skeleton_length_px_rough = length_pixels,
    GFAP_skeleton_length_um_rough = skeleton_length_um,
    GFAP_skeleton_length_um_per_mm2_rough =
      if (calibrated && fov_area_mm2 > 0) skeleton_length_um / fov_area_mm2 else NA_real_,
    GFAP_internal_branchpoints_rough = branchpoint_count,
    GFAP_internal_branchpoints_per_mm2_rough =
      if (calibrated && fov_area_mm2 > 0) branchpoint_count / fov_area_mm2 else NA_real_,
    GFAP_internal_endpoints_rough = endpoint_count,
    GFAP_internal_endpoints_per_mm2_rough =
      if (calibrated && fov_area_mm2 > 0) endpoint_count / fov_area_mm2 else NA_real_,
    GFAP_skeleton_network_fragments_rough = network_fragments,
    GFAP_area_to_skeleton_width_um_rough =
      if (calibrated && skeleton_length_um > 0) positive_area_um2 / skeleton_length_um else NA_real_
  )

  list(
    skeleton = skeleton,
    endpoint_pixels = endpoint_pixels,
    branchpoint_pixels = branchpoint_pixels,
    metrics = metrics
  )
}

count_dapi_nuclei <- function(dapi) {
  if (is.null(dapi)) return(NA_integer_)
  dapi_normalized <- normalize01(dapi)
  blurred <- EBImage::gblur(EBImage::Image(dapi_normalized), sigma = 1)
  binary <- blurred > EBImage::otsu(blurred)
  binary <- EBImage::opening(binary, EBImage::makeBrush(3, shape = "disc"))
  binary <- EBImage::fillHull(binary)
  labels <- as.array(EBImage::watershed(EBImage::distmap(binary), tolerance = 1, ext = 1))
  labels[!as.array(binary)] <- 0
  as.integer(max(labels, na.rm = TRUE))
}

quantify_gfap <- function(gfap, mask, dapi = NULL, pixel_size_um = NA_real_,
                          morphology = NULL) {
  positive <- as.logical(mask)
  positive_pixels <- sum(positive)
  total_pixels <- length(positive)
  integrated_positive <- if (positive_pixels) sum(gfap[positive], na.rm = TRUE) else 0
  integrated_full_fov <- sum(gfap, na.rm = TRUE)
  calibrated <- is.finite(pixel_size_um) && pixel_size_um > 0
  pixel_area_um2 <- if (calibrated) pixel_size_um^2 else NA_real_
  positive_area_um2 <- if (calibrated) positive_pixels * pixel_area_um2 else NA_real_
  full_fov_area_um2 <- if (calibrated) total_pixels * pixel_area_um2 else NA_real_
  signal_per_positive_um2 <- if (calibrated && positive_area_um2 > 0) {
    integrated_positive / positive_area_um2
  } else {
    NA_real_
  }
  signal_per_full_fov_um2 <- if (calibrated && full_fov_area_um2 > 0) {
    integrated_full_fov / full_fov_area_um2
  } else {
    NA_real_
  }

  intensity_and_area <- data.frame(
    GFAP_area_percent = 100 * positive_pixels / total_pixels,
    GFAP_positive_area_um2 = positive_area_um2,
    GFAP_signal_per_positive_um2_AU = signal_per_positive_um2,
    GFAP_signal_per_full_FOV_um2_AU = signal_per_full_fov_um2,
    GFAP_positive_pixels = positive_pixels,
    Total_pixels = total_pixels,
    Full_FOV_area_um2 = full_fov_area_um2,
    GFAP_mean_positive = if (positive_pixels) mean(gfap[positive], na.rm = TRUE) else NA_real_,
    GFAP_integrated_positive = integrated_positive,
    GFAP_mean_full_FOV = mean(gfap, na.rm = TRUE),
    DAPI_nuclei_rough = count_dapi_nuclei(dapi)
  )

  if (is.null(morphology)) {
    morphology <- analyze_gfap_skeleton(mask, pixel_size_um)
  }
  cbind(intensity_and_area, morphology$metrics)
}

# ================================================================
# SHINY UI
# ================================================================

ui <- navbarPage(
  title = "GFAP IF Quantifier",

  tabPanel(
    "1. Import / Unmix",
    sidebarLayout(
      sidebarPanel(
        fileInput("images", "Upload images", multiple = TRUE,
                  accept = c(".im3", ".tif", ".tiff", ".jpg", ".jpeg", ".png")),
        selectInput("current_file", "Representative / preview file", choices = NULL),
        selectInput(
          "unmix_method", "IM3 unmixing method",
          choices = c(
            "Automatic: CSV if supplied, otherwise NMF" = "auto",
            "Library-free NMF" = "nmf",
            "Reference CSV + NNLS" = "csv"
          ), selected = "auto"
        ),
        fileInput("library_csv", "Optional reference-spectrum CSV", accept = ".csv"),
        tags$hr(),
        numericInput("nmf_rank", "NMF components", value = 3, min = 2, max = 6, step = 1),
        numericInput("nmf_iter", "NMF iterations", value = 250, min = 25, max = 2000, step = 25),
        numericInput("nmf_pixels", "NMF training pixels", value = 50000, min = 1000, max = 250000, step = 5000),
        numericInput("nmf_seed", "NMF random seed", value = 1, min = 1, step = 1),
        numericInput("nnls_iter", "Projection iterations", value = 100, min = 10, max = 1000, step = 10),
        actionButton("reload", "Load / decompose", class = "btn-primary"),
        helpText("For a batch, fit NMF on one representative image, verify the components, then lock the settings.")
      ),
      mainPanel(
        h4("Image information"),
        verbatimTextOutput("file_info"),
        h4("Component spectra"),
        plotOutput("spectra_plot", height = "320px"),
        DTOutput("spectra_table")
      )
    )
  ),

  tabPanel(
    "2. Components",
    sidebarLayout(
      sidebarPanel(
        selectInput("gfap_component", "GFAP component", choices = NULL),
        selectInput("dapi_component", "DAPI component", choices = NULL),
        selectInput("af_component", "Autofluorescence component", choices = NULL),
        helpText("NMF component numbers are not biological labels. Assign them from both morphology and spectral profile.")
      ),
      mainPanel(
        h4("All available components"),
        plotOutput("component_grid", height = "650px"),
        fluidRow(
          column(4, plotOutput("gfap_plot", height = "300px")),
          column(4, plotOutput("dapi_plot", height = "300px")),
          column(4, plotOutput("af_plot", height = "300px"))
        ),
        plotOutput("composite_plot", height = "420px")
      )
    )
  ),

  tabPanel(
    "3. GFAP Segmentation",
    sidebarLayout(
      sidebarPanel(
        numericInput("background_sigma", "Background sigma (pixels)", value = 15, min = 0, step = 1),
        numericInput("gfap_threshold", "GFAP threshold", value = 0, min = 0, step = 0.001),
        actionButton("auto_threshold", "Suggest Otsu threshold"),
        numericInput("min_pixels", "Remove objects smaller than (pixels)", value = 0, min = 0, step = 1),
        tags$hr(),
        actionButton("lock", "Lock basis + analysis settings", class = "btn-success"),
        verbatimTextOutput("lock_status")
      ),
      mainPanel(
        fluidRow(
          column(6, plotOutput("corrected_plot", height = "400px")),
          column(6, plotOutput("mask_plot", height = "400px"))
        ),
        plotOutput("skeleton_plot", height = "430px"),
        h4("Measurements"),
        tableOutput("metrics")
      )
    )
  ),

  tabPanel(
    "4. Batch Analysis",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput(
          "export_items",
          "Include in batch ZIP",
          choices = c(
            "Assigned GFAP channel" = "gfap",
            "Assigned DAPI channel" = "dapi",
            "Assigned autofluorescence channel" = "af",
            "All decomposition components" = "all_components",
            "Background-corrected GFAP" = "corrected",
            "Clean GFAP-red + DAPI-blue composite (PNG)" = "composite",
            "Binary GFAP mask (PNG)" = "mask",
            "GFAP skeleton QC overlay (PNG)" = "skeleton",
            "Colored GFAP/DAPI/mask overlay (PNG)" = "overlay"
          ),
          selected = c("gfap", "dapi", "corrected", "composite", "mask", "skeleton")
        ),
        selectInput(
          "channel_export_format",
          "Unmixed channel format",
          choices = c(
            "16-bit TIFF (fixed batch scale)" = "tiff",
            "PNG preview (per-image contrast)" = "png",
            "JPEG preview (per-image contrast)" = "jpeg"
          ),
          selected = "tiff"
        ),
        conditionalPanel(
          condition = "input.channel_export_format == 'jpeg'",
          sliderInput("jpeg_quality", "JPEG quality", min = 50, max = 100, value = 95, step = 1)
        ),
        helpText("TIFF channel scaling is fixed when settings are locked. PNG/JPEG channels are display previews only."),
        checkboxInput(
          "add_scale_bar",
          "Add calibrated scale bar to clean composite",
          value = FALSE
        ),
        conditionalPanel(
          condition = "input.add_scale_bar == true",
          numericInput("scale_bar_um", "Scale-bar length (µm)", value = 100, min = 1, step = 10),
          selectInput(
            "scale_bar_color", "Scale-bar color",
            choices = c("White" = "white", "Black" = "black"), selected = "white"
          ),
          numericInput(
            "scale_bar_thickness", "Scale-bar thickness (pixels)",
            value = 6, min = 2, max = 30, step = 1
          ),
          helpText("Uses MillimetersPerPixel metadata. No bar is drawn when calibration is unavailable.")
        ),
        actionButton("run_batch", "Run batch analysis", class = "btn-primary"),
        tags$br(), tags$br(),
        downloadButton("download_results", "Download CSV"),
        tags$br(), tags$br(),
        uiOutput("download_bundle_ui")
      ),
      mainPanel(verbatimTextOutput("batch_status"), DTOutput("batch_table"))
    )
  ),

  tabPanel(
    "5. Notes",
    fluidPage(
      h3("Recommended interpretation"),
      p("Primary outcome: GFAP-positive area (%)."),
      p("Library-free NMF coefficients are not calibrated fluorescence intensities; interpret intensity outcomes cautiously."),
      p("Do not refit NMF or choose a new threshold independently for every experimental image."),
      p("Fit on a representative image, verify the biological assignments, lock the component basis and segmentation settings, then apply them consistently."),
      p("Validate the component assignment against morphology, known controls, and any available Nuance composite or vendor-unmixed image before publication."),
      p("The DAPI nuclei count is a rough quality-control estimate, not a validated stereological cell count."),
      p("Batch TIFF exports are 16-bit images scaled with fixed limits from the locked representative image. PNG/JPEG channel exports are contrast-enhanced previews, not quantitative intensity files."),
      p("An optional scale bar can be added to the clean red-GFAP/blue-DAPI composite. It is drawn only when calibrated MillimetersPerPixel metadata is present in the IM3 file."),
      p("Calibrated area is GFAP-positive mask area in square micrometers, not the area of an individual astrocyte. Individual-cell area requires a separate cell-segmentation method."),
      p("Signal-per-area values are reported as arbitrary NMF fluorescence units per square micrometer (AU/µm²), not absolute calibrated photon counts."),
      p("Skeleton rules (validated November 2026): after Zhang-Suen thinning, redundant corner and staircase pixels are removed; skeleton pixels with three or more neighbours that lie within 2 pixels of each other count as one branch point; diagonal steps already bridged by an orthogonal path are not added to length. See validation/VALIDATION.md."),
      p("Skeleton measurements are secondary field-level descriptors of the thresholded GFAP network. They include approximate skeleton length, internal branch points, internal endpoints, network fragments, and an area-to-skeleton estimate of process thickness."),
      p("Branch points and endpoints are marked as rough because touching astrocytes, projection of three-dimensional processes into two dimensions, image borders, threshold choice, and GFAP-negative fine processes prevent valid cell-by-cell reconstruction."),
      p("GFAP expression and morphology alone do not establish a complete reactive-astrocyte state; interpret them alongside disease context and, when available, additional molecular or functional markers.")
    )
  )
)

# ================================================================
# SHINY SERVER
# ================================================================

server <- function(input, output, session) {
  state <- reactiveValues(
    files = NULL,
    current = NULL,
    locked = NULL,
    results = NULL,
    export_dir = NULL,
    exported_count = 0L,
    zip_path = NULL
  )

  load_file_row <- function(row) {
    if (row$ext == "im3") {
      cube <- read_im3_cube(row$path)
      pixel_size_um <- extract_im3_pixel_size_um(row$path)
      use_csv <- identical(input$unmix_method, "csv") ||
        (identical(input$unmix_method, "auto") && !is.null(input$library_csv))

      if (identical(input$unmix_method, "csv") && is.null(input$library_csv)) {
        stop("Reference CSV mode was selected, but no CSV was uploaded.")
      }

      if (use_csv) {
        basis <- read_basis_csv(input$library_csv$datapath, dim(cube)[3])
        fit_info <- NULL
        mode <- "Reference CSV + non-negative projection"
      } else {
        fit <- fit_nmf_basis(
          cube,
          rank = input$nmf_rank,
          max_iter = input$nmf_iter,
          sample_pixels = input$nmf_pixels,
          seed = input$nmf_seed
        )
        basis <- fit$basis
        fit_info <- fit
        mode <- "Library-free NMF + non-negative projection"
      }

      components <- unmix_fixed_basis(cube, basis, max_iter = input$nnls_iter)
      list(
        type = "IM3 multispectral cube",
        filename = row$original_name,
        dimensions = dim(cube),
        mode = mode,
        basis = basis,
        fit_info = fit_info,
        pixel_size_um = pixel_size_um,
        components = components
      )
    } else {
      image <- read_composite_image(row$path)
      list(
        type = "Composite image",
        filename = row$original_name,
        dimensions = dim(image),
        mode = "Direct RGB/gray channels (no spectral decomposition)",
        basis = NULL,
        fit_info = NULL,
        pixel_size_um = NA_real_,
        components = composite_components(image)
      )
    }
  }

  do_load <- function(display_name) {
    req(state$files, display_name)
    index <- match(display_name, state$files$name)
    if (is.na(index)) return()
    row <- state$files[index, , drop = FALSE]
    withProgress(message = "Loading / decomposing image", value = 0, {
      state$current <- tryCatch({
        incProgress(0.15, detail = "Reading image")
        result <- load_file_row(row)
        incProgress(0.85, detail = "Finished")
        result
      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
        NULL
      })
    })
  }

  observeEvent(input$images, {
    req(input$images)
    state$files <- stage_uploaded_files(input$images)
    state$locked <- NULL
    state$results <- NULL
    state$export_dir <- NULL
    state$exported_count <- 0L
    state$zip_path <- NULL
    first_name <- state$files$name[1]
    updateSelectInput(session, "current_file", choices = state$files$name, selected = first_name)
    do_load(first_name)
  })

  observeEvent(input$current_file, {
    req(input$current_file, state$files)
    if (is.null(state$current) || !identical(input$current_file, state$current$filename)) {
      do_load(input$current_file)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$reload, {
    req(input$current_file, state$files)
    do_load(input$current_file)
  })

  observeEvent(state$current, {
    req(state$current)
    choices <- names(state$current$components)
    gfap_default <- if ("Red" %in% choices) "Red" else choices[1]
    dapi_default <- if ("Blue" %in% choices) "Blue" else choices[min(2, length(choices))]
    af_default <- if (length(choices) >= 3 && !"Blue" %in% choices) choices[3] else ""
    updateSelectInput(session, "gfap_component", choices = choices, selected = gfap_default)
    updateSelectInput(session, "dapi_component", choices = c(None = "", choices), selected = dapi_default)
    updateSelectInput(session, "af_component", choices = c(None = "", choices), selected = af_default)
  })

  output$file_info <- renderPrint({
    req(state$current)
    cat("File:", state$current$filename, "\n")
    cat("Type:", state$current$type, "\n")
    cat("Dimensions [Y x X x planes/channels]:", paste(state$current$dimensions, collapse = " x "), "\n")
    cat("Mode:", state$current$mode, "\n")
    cat("Components:", paste(names(state$current$components), collapse = ", "), "\n")
    if (is.finite(state$current$pixel_size_um)) {
      cat("Calibrated pixel size:", signif(state$current$pixel_size_um, 9), "µm/pixel\n")
    } else {
      cat("Calibrated pixel size: not available\n")
    }
    if (!is.null(state$current$fit_info)) {
      cat("NMF training pixels:", state$current$fit_info$sampled_pixels, "\n")
      cat("NMF iterations:", state$current$fit_info$iterations, "\n")
      cat("NMF relative reconstruction error:", signif(state$current$fit_info$relative_error, 5), "\n")
    }
  })

  output$spectra_plot <- renderPlot({
    req(state$current)
    data <- basis_to_df(state$current$basis)
    if (!nrow(data)) {
      plot.new(); text(0.5, 0.5, "Composite image: no spectral basis"); return()
    }
    print(
      ggplot(data, aes(Plane, Magnitude, color = Component)) +
        geom_line(linewidth = 0.9) + geom_point(size = 1.2) +
        theme_minimal(base_size = 12) + labs(x = "Spectral plane", y = "Normalized magnitude")
    )
  })

  output$spectra_table <- renderDT({
    req(state$current)
    datatable(basis_to_df(state$current$basis), rownames = FALSE,
              options = list(pageLength = 8, scrollX = TRUE))
  })

  selected_gfap <- reactive({
    req(state$current, input$gfap_component)
    state$current$components[[input$gfap_component]]
  })

  selected_dapi <- reactive({
    req(state$current)
    if (is.null(input$dapi_component) || !nzchar(input$dapi_component)) return(NULL)
    state$current$components[[input$dapi_component]]
  })

  selected_af <- reactive({
    req(state$current)
    if (is.null(input$af_component) || !nzchar(input$af_component)) return(NULL)
    state$current$components[[input$af_component]]
  })

  output$component_grid <- renderPlot({
    req(state$current)
    components <- state$current$components
    columns <- if (length(components) <= 3) length(components) else 3
    rows <- ceiling(length(components) / columns)
    old <- par(mfrow = c(rows, columns), mar = c(1, 1, 3, 1))
    on.exit(par(old), add = TRUE)
    for (name in names(components)) draw_gray(components[[name]], name)
  })

  output$gfap_plot <- renderPlot({
    req(selected_gfap())
    plot_gray(selected_gfap(), paste("GFAP:", input$gfap_component))
  })

  output$dapi_plot <- renderPlot({
    dapi <- selected_dapi()
    if (is.null(dapi)) { plot.new(); text(0.5, 0.5, "No DAPI selected") }
    else plot_gray(dapi, paste("DAPI:", input$dapi_component))
  })

  output$af_plot <- renderPlot({
    af <- selected_af()
    if (is.null(af)) { plot.new(); text(0.5, 0.5, "No autofluorescence selected") }
    else plot_gray(af, paste("Autofluorescence:", input$af_component))
  })

  output$composite_plot <- renderPlot({
    gfap <- selected_gfap()
    dapi <- selected_dapi()
    req(gfap)
    if (is.null(dapi)) dapi <- matrix(0, nrow(gfap), ncol(gfap))
    plot_rgb_matrix(gfap, matrix(0, nrow(gfap), ncol(gfap)), dapi,
                    "Assigned GFAP (red) + DAPI (blue)")
  })

  corrected_gfap <- reactive({
    req(selected_gfap(), input$background_sigma)
    background_correct(selected_gfap(), input$background_sigma)
  })

  update_threshold <- function() {
    threshold <- suggest_threshold(corrected_gfap())
    updateNumericInput(session, "gfap_threshold", value = signif(threshold, 6), min = 0,
                       step = max(signif(threshold / 100, 3), 0.001))
  }

  observeEvent(list(state$current, input$gfap_component, input$background_sigma), {
    req(selected_gfap())
    update_threshold()
  }, ignoreInit = TRUE)

  observeEvent(input$auto_threshold, {
    req(corrected_gfap())
    update_threshold()
  })

  current_mask <- reactive({
    req(corrected_gfap(), input$gfap_threshold)
    make_gfap_mask(corrected_gfap(), input$gfap_threshold, input$min_pixels)
  })

  output$corrected_plot <- renderPlot({
    req(corrected_gfap())
    plot_gray(corrected_gfap(), "Background-corrected GFAP")
  })

  output$mask_plot <- renderPlot({
    req(corrected_gfap(), current_mask())
    plot_mask_overlay(corrected_gfap(), selected_dapi(), current_mask(),
                      "GFAP-positive pixels = green")
  })

  current_morphology <- reactive({
    req(current_mask(), state$current)
    analyze_gfap_skeleton(
      current_mask(),
      pixel_size_um = state$current$pixel_size_um
    )
  })

  output$skeleton_plot <- renderPlot({
    req(corrected_gfap(), current_morphology())
    plot_skeleton_overlay(corrected_gfap(), current_morphology())
  })

  current_metrics <- reactive({
    req(corrected_gfap(), current_mask(), current_morphology())
    quantify_gfap(
      corrected_gfap(), current_mask(), selected_dapi(),
      pixel_size_um = state$current$pixel_size_um,
      morphology = current_morphology()
    )
  })

  output$metrics <- renderTable({
    data <- current_metrics()
    data.frame(
      Measurement = names(data),
      Value = vapply(data, function(column) format(column[[1]], digits = 6, scientific = FALSE), character(1))
    )
  })

  current_parameters <- reactive({
    req(state$current)
    component_scales <- vapply(
      state$current$components,
      image_export_scale,
      numeric(1)
    )
    corrected_scale <- image_export_scale(corrected_gfap())
    list(
      gfap_component = input$gfap_component,
      dapi_component = input$dapi_component %||% "",
      af_component = input$af_component %||% "",
      background_sigma = input$background_sigma,
      threshold = input$gfap_threshold,
      min_pixels = input$min_pixels,
      nnls_iter = input$nnls_iter,
      basis = state$current$basis,
      mode = state$current$mode,
      source_type = state$current$type,
      component_scales = component_scales,
      corrected_gfap_scale = corrected_scale
    )
  })

  observeEvent(input$lock, {
    state$locked <- current_parameters()
    showNotification("Component basis, assignments, and segmentation settings locked.", type = "message")
  })

  output$lock_status <- renderPrint({
    if (is.null(state$locked)) {
      cat("Settings are NOT locked.")
    } else {
      cat("LOCKED SETTINGS\n")
      cat("Mode:", state$locked$mode, "\n")
      cat("GFAP:", state$locked$gfap_component, "\n")
      cat("DAPI:", state$locked$dapi_component, "\n")
      cat("Autofluorescence:", state$locked$af_component, "\n")
      cat("Background sigma:", state$locked$background_sigma, "\n")
      cat("GFAP threshold:", state$locked$threshold, "\n")
      cat("Minimum pixels:", state$locked$min_pixels, "\n")
      if (!is.null(state$locked$basis)) cat("Fixed basis components:", ncol(state$locked$basis), "\n")
      cat("Fixed TIFF scaling: representative-image 99.9th percentiles\n")
    }
  })

  write_batch_exports <- function(row, components, gfap, dapi, af,
                                  corrected, mask, morphology, pixel_size_um, parameters,
                                  export_options, export_dir) {
    items <- export_options$items %||% character()
    if (!length(items)) return(invisible(character()))

    stem <- safe_file_stem(row$original_name)
    format <- export_options$channel_format
    extension <- channel_export_extension(format)
    quality <- export_options$jpeg_quality
    written <- character()

    write_channel <- function(image, label, scale) {
      if (is.null(image)) return(invisible(NULL))
      path <- file.path(export_dir, paste0(stem, "_", label, ".", extension))
      write_gray_export(
        image, path, format = format, fixed_scale = scale,
        jpeg_quality = quality
      )
      written <<- c(written, path)
      invisible(path)
    }

    if ("gfap" %in% items) {
      write_channel(gfap, "GFAP_unmixed", parameters$component_scales[[parameters$gfap_component]])
    }
    if ("dapi" %in% items && !is.null(dapi)) {
      write_channel(dapi, "DAPI_unmixed", parameters$component_scales[[parameters$dapi_component]])
    }
    if ("af" %in% items && !is.null(af)) {
      write_channel(af, "autofluorescence_unmixed", parameters$component_scales[[parameters$af_component]])
    }
    if ("all_components" %in% items) {
      for (component_name in names(components)) {
        write_channel(
          components[[component_name]],
          paste0("component_", safe_file_stem(component_name)),
          parameters$component_scales[[component_name]]
        )
      }
    }
    if ("corrected" %in% items) {
      write_channel(corrected, "GFAP_background_corrected", parameters$corrected_gfap_scale)
    }
    if ("mask" %in% items) {
      path <- file.path(export_dir, paste0(stem, "_GFAP_mask.png"))
      write_binary_mask_png(mask, path)
      written <- c(written, path)
    }
    if ("skeleton" %in% items) {
      path <- file.path(export_dir, paste0(stem, "_GFAP_skeleton_QC.png"))
      write_skeleton_overlay_png(corrected, morphology, path)
      written <- c(written, path)
    }
    if ("composite" %in% items) {
      requested_bar <- isTRUE(export_options$add_scale_bar)
      bar_spec <- scale_bar_spec(
        ncol(gfap), nrow(gfap), pixel_size_um, export_options$scale_bar_um
      )
      suffix <- if (requested_bar && isTRUE(bar_spec$valid)) {
        bar_name <- safe_file_stem(
          paste0(format(export_options$scale_bar_um, trim = TRUE, scientific = FALSE), "um")
        )
        paste0("_GFAP_red_DAPI_blue_composite_", bar_name, "_scalebar.png")
      } else {
        "_GFAP_red_DAPI_blue_composite.png"
      }
      path <- file.path(export_dir, paste0(stem, suffix))
      write_clean_composite_png(
        gfap, dapi, path,
        pixel_size_um = if (requested_bar) pixel_size_um else NA_real_,
        scale_bar_um = if (requested_bar) export_options$scale_bar_um else NA_real_,
        scale_bar_color = export_options$scale_bar_color,
        scale_bar_thickness = export_options$scale_bar_thickness
      )
      written <- c(written, path)
    }
    if ("overlay" %in% items) {
      path <- file.path(export_dir, paste0(stem, "_GFAP_overlay.png"))
      write_overlay_png(gfap, dapi, mask, path)
      written <- c(written, path)
    }

    invisible(written)
  }

  write_export_documentation <- function(parameters, export_options, export_dir) {
    component_table <- data.frame(
      Image = paste0("component: ", names(parameters$component_scales)),
      Locked_scale_max = as.numeric(parameters$component_scales),
      stringsAsFactors = FALSE
    )
    scaling_table <- rbind(
      component_table,
      data.frame(
        Image = "background-corrected GFAP",
        Locked_scale_max = parameters$corrected_gfap_scale,
        stringsAsFactors = FALSE
      )
    )
    scaling_table$TIFF_interpretation <- paste(
      "Approximate original value = stored 16-bit fraction * Locked_scale_max;",
      "values above the maximum are clipped."
    )
    utils::write.csv(
      scaling_table,
      file.path(export_dir, "image_export_scaling.csv"),
      row.names = FALSE
    )

    notes <- c(
      "GFAP IF Quantifier batch image export",
      "",
      paste("Channel format:", export_options$channel_format),
      "TIFF channels are 16-bit and use fixed scales from the locked representative image.",
      "Use image_export_scaling.csv to interpret TIFF pixel values.",
      "PNG/JPEG channel files use per-image contrast normalization and are previews only.",
      "The clean composite displays GFAP in red and DAPI in blue with no green mask.",
      paste("Scale bar requested:", isTRUE(export_options$add_scale_bar)),
      paste("Requested scale-bar length (um):", export_options$scale_bar_um),
      "Scale bars use each IM3 file's MillimetersPerPixel metadata; files without calibration receive no bar.",
      "GFAP_positive_area_um2 = GFAP-positive pixels multiplied by calibrated pixel area.",
      "GFAP_signal_per_positive_um2_AU = integrated positive GFAP signal divided by GFAP-positive area.",
      "GFAP_signal_per_full_FOV_um2_AU = total background-corrected GFAP signal divided by calibrated FOV area.",
      "Signal density is in arbitrary component units per square micrometer, not calibrated photon counts.",
      "GFAP-positive area is not individual-cell area because this pipeline does not separate astrocytes.",
      "Skeleton QC colors: yellow = skeleton, red = branch-point region, cyan = internal endpoint.",
      "Skeleton outputs are rough field-level descriptors, not branches per individual astrocyte.",
      "Skeleton length and point densities use each file's calibrated pixel size when available.",
      "Area-to-skeleton width is an approximate GFAP process-thickness index and includes somata/overlaps.",
      "Binary masks and colored overlays are lossless PNG files.",
      "Quantification in GFAP_quantification.csv used the original numeric component values,",
      "not the scaled values in exported display images."
    )
    writeLines(notes, file.path(export_dir, "README_image_exports.txt"))
  }

  create_batch_zip <- function(export_dir) {
    files <- list.files(
      export_dir,
      recursive = TRUE,
      all.files = FALSE,
      include.dirs = FALSE
    )
    if (!length(files)) stop("No batch-export files were created.")

    zip_path <- tempfile(pattern = "GFAP_batch_results_", fileext = ".zip")
    zip::zip(
      zipfile = zip_path,
      files = files,
      root = export_dir,
      mode = "mirror",
      include_directories = FALSE,
      compression_level = 6
    )

    if (!file.exists(zip_path) || is.na(file.info(zip_path)$size) || file.info(zip_path)$size <= 0) {
      stop("The ZIP archive was not created correctly.")
    }
    zip_path
  }

  process_batch_file <- function(row, parameters, export_options, export_dir) {
    pixel_size_um <- NA_real_
    if (row$ext == "im3") {
      if (is.null(parameters$basis)) stop("No locked spectral basis is available for IM3 processing.")
      cube <- read_im3_cube(row$path)
      pixel_size_um <- extract_im3_pixel_size_um(row$path)
      if (dim(cube)[3] != nrow(parameters$basis)) {
        stop("Spectral plane count differs from the representative image.")
      }
      components <- unmix_fixed_basis(cube, parameters$basis, parameters$nnls_iter)
    } else {
      components <- composite_components(read_composite_image(row$path))
    }

    if (!parameters$gfap_component %in% names(components)) {
      stop("GFAP component '", parameters$gfap_component, "' was not found.")
    }
    gfap <- components[[parameters$gfap_component]]
    dapi <- NULL
    if (nzchar(parameters$dapi_component) && parameters$dapi_component %in% names(components)) {
      dapi <- components[[parameters$dapi_component]]
    }
    af <- NULL
    if (nzchar(parameters$af_component) && parameters$af_component %in% names(components)) {
      af <- components[[parameters$af_component]]
    }
    corrected <- background_correct(gfap, parameters$background_sigma)
    mask <- make_gfap_mask(corrected, parameters$threshold, parameters$min_pixels)
    morphology <- analyze_gfap_skeleton(mask, pixel_size_um)
    metrics <- quantify_gfap(
      corrected, mask, dapi,
      pixel_size_um = pixel_size_um,
      morphology = morphology
    )
    if (isTRUE(export_options$add_scale_bar) && "composite" %in% export_options$items) {
      bar_spec <- scale_bar_spec(
        ncol(gfap), nrow(gfap), pixel_size_um, export_options$scale_bar_um
      )
      scale_bar_status <- bar_spec$status
    } else if (isTRUE(export_options$add_scale_bar)) {
      scale_bar_status <- "Not drawn: clean composite export was not selected"
    } else {
      scale_bar_status <- "Not requested"
    }
    write_batch_exports(
      row, components, gfap, dapi, af, corrected, mask,
      morphology, pixel_size_um, parameters, export_options, export_dir
    )
    cbind(
      data.frame(File = row$original_name, stringsAsFactors = FALSE),
      metrics,
      data.frame(
        Unmixing_mode = parameters$mode,
        GFAP_component = parameters$gfap_component,
        DAPI_component = parameters$dapi_component,
        Background_sigma = parameters$background_sigma,
        GFAP_threshold = parameters$threshold,
        Minimum_pixels = parameters$min_pixels,
        Pixel_size_um = pixel_size_um,
        Scale_bar_um = if (isTRUE(export_options$add_scale_bar)) export_options$scale_bar_um else NA_real_,
        Scale_bar_status = scale_bar_status,
        stringsAsFactors = FALSE
      )
    )
  }

  observeEvent(input$run_batch, {
    req(state$files)
    if (is.null(state$locked)) {
      showNotification("Lock the verified component basis and analysis settings before batch processing.",
                       type = "warning", duration = 8)
      return()
    }
    state$zip_path <- NULL
    parameters <- state$locked
    export_options <- list(
      items = input$export_items %||% character(),
      channel_format = input$channel_export_format %||% "tiff",
      jpeg_quality = input$jpeg_quality %||% 95,
      add_scale_bar = isTRUE(input$add_scale_bar),
      scale_bar_um = input$scale_bar_um %||% 100,
      scale_bar_color = input$scale_bar_color %||% "white",
      scale_bar_thickness = input$scale_bar_thickness %||% 6
    )
    export_dir <- tempfile(pattern = "GFAP_batch_exports_")
    if (!dir.create(export_dir, recursive = TRUE)) {
      showNotification("Could not create the temporary batch-export folder.", type = "error")
      return()
    }
    results <- vector("list", nrow(state$files))
    withProgress(message = "Batch processing", value = 0, {
      for (i in seq_len(nrow(state$files))) {
        incProgress(1 / nrow(state$files), detail = state$files$original_name[i])
        results[[i]] <- tryCatch(
          process_batch_file(
            state$files[i, , drop = FALSE],
            parameters,
            export_options,
            export_dir
          ),
          error = function(e) data.frame(
            File = state$files$original_name[i], Error = conditionMessage(e), stringsAsFactors = FALSE
          )
        )
      }
    })
    all_columns <- unique(unlist(lapply(results, names)))
    results <- lapply(results, function(x) {
      for (column in setdiff(all_columns, names(x))) x[[column]] <- NA
      x[all_columns]
    })
    state$results <- do.call(rbind, results)
    utils::write.csv(
      state$results,
      file.path(export_dir, "GFAP_quantification.csv"),
      row.names = FALSE
    )
    write_export_documentation(parameters, export_options, export_dir)
    state$export_dir <- export_dir
    state$exported_count <- length(list.files(export_dir, recursive = TRUE))
    state$zip_path <- tryCatch(
      create_batch_zip(export_dir),
      error = function(e) {
        showNotification(
          paste("Batch finished, but ZIP creation failed:", conditionMessage(e)),
          type = "error",
          duration = NULL
        )
        NULL
      }
    )
    if (!is.null(state$zip_path)) {
      showNotification("Batch finished and the ZIP download is ready.", type = "message")
    }
  })

  output$batch_status <- renderPrint({
    if (is.null(state$results)) cat("No batch analysis performed yet.")
    else {
      cat("Finished:", nrow(state$results), "image(s).\n")
      cat("The locked component basis and segmentation settings were used.\n")
      cat("Prepared", state$exported_count, "file(s).\n")
      if (is.null(state$zip_path)) cat("ZIP is not ready; see the error notification.")
      else cat("ZIP download is ready.")
    }
  })

  output$batch_table <- renderDT({
    req(state$results)
    datatable(state$results, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })

  output$download_results <- downloadHandler(
    filename = function() paste0("GFAP_quantification_", Sys.Date(), ".csv"),
    content = function(file) {
      req(state$results)
      utils::write.csv(state$results, file, row.names = FALSE)
    }
  )

  output$download_bundle_ui <- renderUI({
    if (is.null(state$zip_path) || !file.exists(state$zip_path)) {
      helpText("Run batch analysis to prepare the ZIP download.")
    } else {
      downloadButton("download_bundle", "Download results + images (ZIP)")
    }
  })

  output$download_bundle <- downloadHandler(
    filename = function() paste0("GFAP_batch_results_", Sys.Date(), ".zip"),
    contentType = "application/zip",
    content = function(file) {
      req(state$results, state$zip_path, file.exists(state$zip_path))
      if (!file.copy(state$zip_path, file, overwrite = TRUE)) {
        stop("The prepared ZIP could not be copied to the download location.")
      }
    }
  )
}

shinyApp(ui = ui, server = server)
