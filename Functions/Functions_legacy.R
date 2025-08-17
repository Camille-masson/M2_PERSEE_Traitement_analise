library(terra)
library(lubridate)

build_smod_fsca_stack <- function(alpage, YEAR,
                                  smod_hydro_tif,   # SMOD en DOY hydro (1 = 1 sept Y-1)
                                  out_dir,
                                  up_shape = NULL,  # ← NEW
                                  outfile  = NULL) {
  
  YEAR_chr <- as.character(YEAR)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(outfile)) {
    outfile <- file.path(out_dir,
                         paste0("SMOD_FSCA_stack_", alpage, "_", YEAR_chr, "_epsg2154.tif"))
  }
  
  ## 1) Lire SMOD & reprojeter en EPSG:2154
  smod_hydro <- rast(smod_hydro_tif)
  smod_hydro_2154 <- project(smod_hydro, "EPSG:2154", method = "near")
  
  ## (optionnel) lire l'UP et recadrer/masquer TOUT DE SUITE
  if (!is.null(up_shape) && file.exists(up_shape)) {
    up_vect_2154 <- project(vect(up_shape), "EPSG:2154")
    smod_hydro_2154 <- crop(smod_hydro_2154, up_vect_2154) |> mask(up_vect_2154)
  } else {
    up_vect_2154 <- NULL
  }
  
  ## 2) DOY hydro -> DOY calendrier
  start_hydro <- as.Date(paste0(as.integer(YEAR_chr) - 1, "-09-01"))  # 1er sept N-1
  jan1        <- as.Date(paste0(YEAR_chr, "-01-01"))
  dec31       <- as.Date(paste0(YEAR_chr, "-12-31"))
  
  hydro2cal <- function(x) {
    d   <- start_hydro + (as.integer(x) - 1L)
    doy <- as.integer(difftime(d, jan1, units = "days")) + 1L
    doy[d < jan1 | d > dec31] <- NA_integer_
    doy
  }
  
  smod_cal_2154 <- app(smod_hydro_2154, hydro2cal)
  names(smod_cal_2154) <- paste0("SMOD_DOYcal_", YEAR_chr)
  
  ## 3) FSCA par pixel (rang sur l'UP si fournie)
  vals <- values(smod_cal_2154)[, 1]
  N    <- length(vals)
  
  fsca_vals <- rep(NA_real_, N)
  valid_idx <- which(!is.na(vals))
  
  if (length(valid_idx) > 1) {
    valid_vals <- vals[valid_idx]
    rnk        <- rank(valid_vals, ties.method = "min")
    N_valid    <- length(valid_vals)
    fsca_vals[valid_idx] <- (N_valid - rnk) / (N_valid - 1)
  }
  
  fsca_r <- smod_cal_2154
  values(fsca_r) <- fsca_vals
  names(fsca_r)  <- "FSCA_frac"
  
  ## 4) Empilement & écriture
  stack_out <- c(smod_cal_2154, fsca_r)
  writeRaster(stack_out, outfile, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  
  cat("✓ Stack SMOD+FSCA (croppé UP si fourni) écrit :", outfile, "\n")
  invisible(outfile)
}

build_smod_fsca_stack <- function(
    alpage, YEAR,
    smod_hydro_tif,          # SMOD en DOY hydro (1 = 1 sept YEAR-1)
    out_dir,
    up_shape = NULL,         # shapefile UP_*_bis.shp (optionnel)
    outfile  = NULL,
    # ---- options de lissage ----
    smooth     = TRUE,       # activer/désactiver le lissage
    kernel     = c("gauss", "median", "mean", "circle"),
    radius_m   = 60,         # rayon visé (m) pour gauss/mean/circle
    round_DOY  = TRUE,       # arrondir DOY lissé à l'entier
    # ---- gestion des bords UP ----
    buffer_m   = NULL        # si NULL => = radius_m ; sinon taille du buffer (m)
) {
  suppressPackageStartupMessages(library(terra))
  kernel   <- match.arg(kernel)
  YEAR_chr <- as.character(YEAR)
  
  if (!file.exists(smod_hydro_tif))
    stop("Fichier SMOD introuvable : ", smod_hydro_tif)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(outfile)) {
    outfile <- file.path(out_dir,
                         paste0("SMOD_FSCA_stack_", alpage, "_", YEAR_chr, "_epsg2154.tif"))
  }
  if (is.null(buffer_m)) buffer_m <- radius_m
  
  ## 1) Lecture + reprojection EPSG:2154 -------------------------------------
  smod_hydro <- rast(smod_hydro_tif)
  smod_hydro_2154 <- project(smod_hydro, "EPSG:2154", method = "near")
  
  ## 2) Préparation UP + CLIP par buffer (évite effets de bord) --------------
  up_vect_2154 <- NULL
  smod_clip    <- smod_hydro_2154
  if (!is.null(up_shape) && file.exists(up_shape)) {
    up_vect_2154 <- project(vect(up_shape), "EPSG:2154")
    
    # on fabrique un buffer autour de l'UP pour garder des voisins pendant le lissage
    buf_geom <- buffer(up_vect_2154, width = buffer_m)
    
    # on RESTREINT l'emprise au buffer (crop), mais SANS masquer
    smod_clip <- crop(smod_hydro_2154, buf_geom)  # réduit l'emprise, garde les valeurs
  }
  
  ## 3) Conversion DOY hydro -> DOY calendrier --------------------------------
  jan1  <- as.Date(sprintf("%s-01-01", YEAR_chr))
  dec31 <- as.Date(sprintf("%s-12-31", YEAR_chr))
  start_hydro <- as.Date(sprintf("%d-09-01", as.integer(YEAR_chr) - 1))
  ndays <- as.integer(difftime(dec31, jan1, units = "days")) + 1L
  
  hydro2cal <- function(x) {
    x <- as.integer(x)
    d   <- start_hydro + (x - 1L)
    doy <- as.integer(difftime(d, jan1, units = "days")) + 1L
    doy[d < jan1 | d > dec31] <- NA_integer_
    doy
  }
  smod_cal <- app(smod_clip, hydro2cal)  # DOY calendrier (sur emprise buffer)
  
  ## 4) LISSAGE (moving window), AVANT masque UP ------------------------------
  smod_smooth <- smod_cal
  if (smooth) {
    px <- res(smod_cal)[1]                      # taille de pixel (m)
    
    if (kernel == "gauss") {
      # σ en mètres ; ≈ radius_m ≃ 3*σ
      sigma_m  <- max(radius_m / 3, px / 2)     # évite noyau dégénéré
      size_m   <- max(radius_m, 3 * sigma_m)    # force une fenêtre suffisante
      K <- focalMat(smod_cal, d = c(sigma_m, size_m), type = "Gauss")
      # secours si noyau trop petit
      if (!is.finite(sum(K)) || nrow(K) * ncol(K) == 1L) {
        sigma_m <- max(1.5 * px, sigma_m); size_m <- max(3 * px, size_m)
        K <- focalMat(smod_cal, d = c(sigma_m, size_m), type = "Gauss")
      }
      K <- K / sum(K, na.rm = TRUE)
      smod_smooth <- focal(smod_cal, w = K, fun = "sum",
                           na.policy = "omit", pad = TRUE)
      
    } else if (kernel == "median") {
      w <- 2 * ceiling(radius_m / px) + 1       # fenêtre impaire
      smod_smooth <- focal(smod_cal, w = w, fun = median,
                           na.policy = "omit", pad = TRUE)
      
    } else if (kernel == "mean") {
      w <- 2 * ceiling(radius_m / px) + 1
      W <- matrix(1, nrow = w, ncol = w); W <- W / sum(W)
      smod_smooth <- focal(smod_cal, w = W, fun = "sum",
                           na.policy = "omit", pad = TRUE)
      
    } else if (kernel == "circle") {
      W <- focalMat(smod_cal, d = radius_m, type = "circle")
      W <- W / sum(W, na.rm = TRUE)
      smod_smooth <- focal(smod_cal, w = W, fun = "sum",
                           na.policy = "omit", pad = TRUE)
    }
    
    if (round_DOY) smod_smooth <- round(smod_smooth)
    smod_smooth <- clamp(smod_smooth, lower = 1, upper = ndays, values = TRUE)
  }
  
  ## 5) Masque final à l’UP (valeurs uniquement DANS l’UP) --------------------
  if (!is.null(up_vect_2154)) {
    smod_final <- mask(smod_smooth, up_vect_2154)  # coupe net à l'UP
  } else {
    smod_final <- smod_smooth
  }
  names(smod_final) <- paste0("SMOD_", YEAR_chr)
  
  ## 6) FSCA (rang) CALCULÉ SUR L’UP UNIQUEMENT --------------------------------
  v  <- values(smod_final)[, 1]
  fs <- rep(NA_real_, length(v))
  idx <- which(!is.na(v))
  if (length(idx) > 1L) {
    vv  <- v[idx]
    rnk <- rank(vv, ties.method = "min")   # 1 = déneige le plus tôt
    N   <- length(vv)
    fs[idx] <- (N - rnk) / (N - 1)        # tôt → ~1 ; tard → ~0
  }
  fsca_r <- setValues(smod_final, fs)
  names(fsca_r) <- paste0("FSCA_", YEAR_chr)
  
  ## 7) Empilement & écriture --------------------------------------------------
  stack_out <- c(smod_final, fsca_r)       # ordre: FSCA_YYYY, SMOD_YYYY
  writeRaster(stack_out, outfile, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  
  cat("✓ Stack SMOD+FSCA (",
      if (smooth) paste0("smoothed ", kernel, ", radius=", radius_m, "m; ") else "no smoothing; ",
      if (!is.null(up_vect_2154)) paste0("buffer=", buffer_m, "m; masked to UP") else "no UP mask",
      ") écrit : ", outfile, "\n", sep = "")
  invisible(outfile)
}


build_pheno_dataset <- function(alpage, YEAR,
                                irg_tif, ndvi_tif,
                                habitat_tif, dah_tif,
                                smod_fsca_stack_tif,      # ← NOUVEAU : stack 2 bandes (SMOD + FSCA)
                                extra_tif,
                                charge_prev_rds,
                                first_use_rds,
                                lut_habitat_csv,
                                up_shape,
                                output_rds,
                                DOY_range = 60:365) {
  suppressPackageStartupMessages({
    library(terra);  library(dplyr);  library(tidyr); library(readr)
    library(stringr); library(purrr); library(conflicted)
  })
  conflict_prefer("select","dplyr", quiet=TRUE)
  conflict_prefer("filter","dplyr", quiet=TRUE)
  
  YEAR_chr <- as.character(YEAR)
  
  ## 0. Vérification des fichiers -------------------------------------------
  needed <- c(irg_tif, ndvi_tif, habitat_tif, dah_tif,
              smod_fsca_stack_tif,
              extra_tif, charge_prev_rds, first_use_rds,
              lut_habitat_csv, up_shape)
  miss <- needed[!file.exists(needed)]
  if (length(miss))
    stop("Fichiers manquants :\n", paste(" •", miss, collapse = "\n"))
  dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
  
  ## 1. Grille de référence --------------------------------------------------
  up_vect   <- vect(up_shape)                 # polygone UP
  meta_orig <- list(
    IRG   = rast(irg_tif)[[1]],
    NDVI  = rast(ndvi_tif)[[1]],
    HAB   = rast(habitat_tif),
    DAH   = rast(dah_tif),
    SMFS  = rast(smod_fsca_stack_tif),        # stack 2 bandes
    EXTRA = rast(extra_tif)[[1]]
  )
  # résolution la plus fine comme base (max nb de pixels -> min prod(res)? ici on garde comme avant)
  base_res <- res(meta_orig[[which.max(sapply(meta_orig, \(r) prod(res(r))))]])
  template <- rast(ext = ext(up_vect), resolution = base_res, crs = crs(up_vect))
  
  read_align <- function(file, method = "bilinear")
    rast(file) |>
    project(template, method = method) |>
    crop(up_vect) |>
    mask(up_vect)
  
  irg_stack   <- read_align(irg_tif,            "bilinear")
  ndvi_stack  <- read_align(ndvi_tif,           "bilinear")
  habitat_r   <- read_align(habitat_tif,        "near")
  dah_r       <- read_align(dah_tif,            "bilinear")
  smfs_stack  <- read_align(smod_fsca_stack_tif,"near")      # valeurs catégorielles / fraction
  extra_stack <- read_align(extra_tif,          "bilinear")
  
  # On suppose : bande 1 = SMOD, bande 2 = FSCA. Renomme dynamiquement :
  smod_r <- smfs_stack[[1]] ; names(smod_r) <- paste0("SMOD_", YEAR_chr)
  fsca_r <- smfs_stack[[2]] ; names(fsca_r) <- paste0("FSCA_", YEAR_chr)
  
  ## 2. Variables statiques --------------------------------------------------
  static_df <- as.data.frame(c(dah_r, fsca_r, smod_r, habitat_r),
                             cells = TRUE, na.rm = FALSE) |>
    setNames(c("cell",
               "dah",
               paste0("FSCA_", YEAR_chr),
               paste0("SMOD_", YEAR_chr),
               "habitat_code")) |>
    distinct(cell, .keep_all = TRUE)
  
  extra_df <- as.data.frame(extra_stack, cells = TRUE, na.rm = FALSE) |>
    distinct(cell, .keep_all = TRUE)
  
  static_df <- left_join(static_df, extra_df, by = "cell",
                         relationship = "one-to-one")
  
  habitat_lut <- read_csv2(lut_habitat_csv, col_types = "ic",
                           locale = locale(encoding = "latin1")) |>
    mutate(code = as.integer(code))
  
  ## 3‑A. Chargement total année n‑1 ----------------------------------------
  charge_df <- readRDS(charge_prev_rds) |> as_tibble()
  if (!"Charge" %in% names(charge_df))
    charge_df <- rename(charge_df, Charge = matches("charge", ignore.case = TRUE))
  
  if (!all(c("x", "y") %in% names(charge_df))) {
    xy0 <- xyFromCell(meta_orig$IRG, charge_df$cell)
    charge_df <- mutate(charge_df, x = xy0[, 1], y = xy0[, 2])
  }
  pts <- vect(cbind(charge_df$x, charge_df$y), crs = crs(meta_orig$IRG)) |>
    project(template)
  charge_df <- mutate(charge_df,
                      x = crds(pts)[, 1], y = crds(pts)[, 2],
                      cell = cellFromXY(template, cbind(x, y))) |>
    filter(!is.na(cell)) |>
    group_by(cell) |>
    summarise(Charge = mean(Charge), .groups = "drop")
  
  col_prev  <- paste0("Charge", as.integer(YEAR_chr) - 1L)
  charge_df <- mutate(charge_df, !!col_prev := Charge)
  
  ## 3‑B. Premier jour d’utilisation (année n) -------------------------------
  first_use_df <- readRDS(first_use_rds) |> as_tibble()
  pts_fu <- vect(cbind(first_use_df$x, first_use_df$y), crs = crs(meta_orig$IRG)) |>
    project(template)
  first_use_df <- mutate(first_use_df,
                         x = crds(pts_fu)[, 1], y = crds(pts_fu)[, 2],
                         cell = cellFromXY(template, cbind(x, y))) |>
    distinct(cell, .keep_all = TRUE) |>
    select(cell,
           first_day_use  = first_day,
           first_date_use = first_date)
  
  ## 4. Passage en format long IRG & NDVI ------------------------------------
  to_long <- function(stk, var)
    as.data.frame(stk, cells = TRUE, xy = TRUE, na.rm = FALSE) |>
    pivot_longer(-c(cell, x, y), names_to = "band", values_to = var) |>
    mutate(DOY = as.integer(stringr::str_extract(band, "\\d{3,}$"))) |>
    filter(DOY %in% DOY_range) |>
    select(-band)
  
  pheno_df <- inner_join(to_long(irg_stack,  "IRG"),
                         to_long(ndvi_stack, "NDVI"),
                         by = c("cell", "x", "y", "DOY"))
  
  ## 5. Fusion finale --------------------------------------------------------
  full_df <- pheno_df |>
    left_join(static_df,    by = "cell", relationship = "many-to-one") |>
    left_join(charge_df,    by = "cell", relationship = "many-to-one") |>
    left_join(first_use_df, by = "cell", relationship = "many-to-one") |>
    mutate(habitat_code = as.integer(habitat_code)) |>
    left_join(habitat_lut,
              by = c("habitat_code" = "code"),
              relationship = "many-to-one") |>
    rename(habitat_label = label) |>
    relocate(habitat_label, .after = habitat_code)
  
  ## 6. Export ---------------------------------------------------------------
  write_rds(full_df, output_rds, compress = "gz")
  cat("✓ Dataset écrit :", output_rds, "\n")
  
  invisible(full_df)
}




###############################################################################
# diag_dataset  – diagnostic PDF (IRG, NDVI, DAH, FSCA_YYYY, SMOD_YYYY,
#                                 AUCg_ONSET10_MAXD + PPI_extra + Charge)
###############################################################################
diag_dataset <- function(
    rds_path,
    up_shape   = NULL,
    sample_n   = 5e5,
    bins       = 60,
    output_dir = dirname(rds_path),
    palette    = c(IRG  = "#1b9e77", NDVI = "#d95f02",
                   fsca = "#7570b3", dah  = "#e7298a",
                   SMOD = "#66a61e",
                   AUCg_ONSET10_MAXD = "#8dd3c7",
                   Charge = "tomato",
                   MINV = "#4daf4a", MAXV = "#377eb8", AMPL = "#984ea3",
                   LENGTH = "#ff7f00", ONSET10 = "#e41a1c",
                   MaxSlope = "#999999",
                   GreenUpDur = "#a65628", GreenDownDur = "#f781bf",
                   AsymSlope = "#66c2a5",
                   LSLOPE = "#fc8d62", RSLOPE = "#8da0cb"),
    view  = FALSE,
    quiet = FALSE) {
  
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2)
    library(sf);    library(scales); library(grid)
    library(conflicted)
  })
  conflict_prefer("intersect", "base", quiet = TRUE)
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  stopifnot(file.exists(rds_path))
  names(palette) <- tolower(names(palette))
  
  # Couleur avec fallback + gestion préfixes FSCA_/SMOD_
  pcol <- function(v) {
    vn <- tolower(v)
    if (vn %in% names(palette)) return(palette[[vn]])
    if (startsWith(vn, "fsca"))  return(palette[["fsca"]])
    if (startsWith(vn, "smod"))  return(palette[["smod"]])
    NULL
  }
  
  ## 0. Chargement -----------------------------------------------------------
  df_full <- readRDS(rds_path)
  if (!is.data.frame(df_full))
    stop("Le fichier rds ne contient pas un data.frame.")
  cnames <- names(df_full); nrows <- nrow(df_full)
  
  ## Statistiques descriptives ----------------------------------------------
  if (!quiet) {
    message("Colonnes : ", paste(cnames, collapse = ", "))
    message("Nombre de pixels (lignes) : ", format(nrows, big.mark = " "))
    message("Nombre de variables (colonnes) : ", length(cnames))
    num_vars <- cnames[sapply(df_full, is.numeric)]
    for (v in num_vars) {
      vals <- df_full[[v]]
      message(sprintf(
        "  %s : mean=%.3f, median=%.3f, sd=%.3f, max=%.3f, NA=%d",
        v, mean(vals, na.rm = TRUE), median(vals, na.rm = TRUE),
        sd(vals, na.rm = TRUE), max(vals, na.rm = TRUE), sum(is.na(vals))))
    }
  }
  
  ## Colonne Charge éventuelle ----------------------------------------------
  charge_col <- cnames[tolower(cnames) == "charge"]
  
  ## 1. Étendue géographique -------------------------------------------------
  have_up <- !is.null(up_shape) && file.exists(up_shape)
  if (have_up) {
    up_sf <- st_read(up_shape, quiet = TRUE)
    bbox  <- st_bbox(up_sf) + c(-50,-50,50,50)
  } else {
    if (!all(c("x","y") %in% cnames))
      stop("Sans up_shape, x et y doivent exister.")
    bbox <- st_bbox(c(xmin = min(df_full$x), ymin = min(df_full$y),
                      xmax = max(df_full$x), ymax = max(df_full$y)),
                    crs = st_crs(2154)) + c(-50,-50,50,50)
  }
  poly_bbox <- st_as_sfc(bbox)
  
  ## 2. Échantillonnage ------------------------------------------------------
  # Colonnes dynamiques FSCA_YYYY / SMOD_YYYY
  fsca_cols <- grep("^FSCA_\\d{4}$", cnames, value = TRUE, ignore.case = TRUE)
  smod_cols <- grep("^SMOD_\\d{4}$", cnames, value = TRUE, ignore.case = TRUE)
  # Colonne AUCg_ONSET10_MAXD (si présente)
  aucg_col  <- intersect("AUCg_ONSET10_MAXD", cnames)
  
  vars_base <- c("IRG","NDVI","dah",
                 "MINV","MAXV","AMPL",
                 "LENGTH","ONSET10","MaxSlope",
                 "GreenUpDur","GreenDownDur","AsymSlope",
                 "LSLOPE","RSLOPE",
                 charge_col)
  vars <- unique(intersect(c(vars_base, fsca_cols, smod_cols, aucg_col), cnames))
  
  rows <- if (nrows <= sample_n) seq_len(nrows) else sample(nrows, sample_n)
  df   <- df_full[rows, vars, drop = FALSE]
  
  ## 3. Graphiques -----------------------------------------------------------
  plots <- list()
  # (a) carte d'emprise
  p_map <- ggplot() +
    { if (have_up) geom_sf(data = up_sf, fill = "grey85",
                           colour = "steelblue", linewidth = .5) } +
    geom_sf(data = poly_bbox, fill = NA, colour = "red3", linewidth = 1) +
    coord_sf(expand = FALSE) +
    labs(title = "Emprise du jeu de données",
         subtitle = basename(rds_path)) +
    theme_void()
  plots[[length(plots) + 1]] <- p_map
  
  # (b) histogrammes pour toutes les variables sauf Charge
  for (v in setdiff(vars, charge_col)) {
    if (!is.numeric(df[[v]])) next
    fill_col <- pcol(v) %||% "grey60"
    plots[[length(plots)+1]] <-
      ggplot(df, aes(.data[[v]])) +
      geom_histogram(bins = bins,
                     fill   = fill_col,
                     colour = "grey25") +
      labs(title     = paste("Distribution de", v),
           subtitle  = paste("Échantillon", format(nrow(df), big.mark = " ")),
           x = v, y  = "Fréquence") +
      theme_minimal(base_size = 9)
  }
  
  ## --- bloc Charge : barre + histo log10 (pixels Charge >= 1) --------------
  if (length(charge_col) == 1) {
    vC   <- charge_col
    
    df_bar <- df %>%
      mutate(etat = ifelse(!is.na(.data[[vC]]) & .data[[vC]] >= 1,
                           "Chargé (≥1)", "Non chargé (<1)"))
    
    plots[[length(plots)+1]] <-
      ggplot(df_bar, aes(etat)) +
      geom_bar(fill = pcol("charge") %||% "tomato") +
      labs(title="Pixels chargés / non chargés",
           x=NULL, y="Pixels") +
      theme_minimal(base_size = 9)
    
    df_pos <- dplyr::filter(df_bar, etat == "Chargé (≥1)")
    
    plots[[length(plots)+1]] <-
      ggplot(df_pos, aes(.data[[vC]])) +
      geom_histogram(bins = bins,
                     fill = pcol("charge") %||% "tomato", colour = "grey25") +
      scale_x_log10(labels = comma) +
      labs(title     = "Distribution de Charge (≥1 UA)",
           subtitle  = paste("Échantillon", format(nrow(df_pos), big.mark=" ")),
           x = "Charge (log10)", y = "Fréquence") +
      theme_minimal(base_size = 9)
  }
  
  ## 4. Export PDF -----------------------------------------------------------
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(output_dir,
                        paste0(tools::file_path_sans_ext(basename(rds_path)),
                               "_plots.pdf"))
  pdf(pdf_path, 7, 5)
  for (p in plots) print(p)    # ← ICI : on "print" les ggplot -> ils sortent bien dans le PDF
  dev.off()
  if (!quiet) message("✓ PDF écrit :", pdf_path)
  
  if (view) for (p in plots) { dev.new(7,5); print(p) }
  invisible(pdf_path)
}







merge_pheno_dataset <- function(year,
                               alpages,
                               root_dir   = getwd(),
                               output_dir = file.path(root_dir, "outputs"),
                               write_fst  = FALSE,
                               combo_name = "Alpe-Sud") {   # <‑‑ NEW
  ## 1. Liste unique des alpages ------------------------------------------------
  alpages <- sort(unique(alpages))
  
  ## 2. Lecture + filtre pour 1 alpage -----------------------------------------
  load_one <- function(ap) {
    in_file <- file.path(output_dir, "10. NDVI Effect",
                         paste0("dataset_", ap),
                         sprintf("dataset_pheno_LONG_%d_%s.rds", year, ap))
    if (!file.exists(in_file))
      stop("Fichier manquant : ", in_file, call. = FALSE)
    
    readRDS(in_file) |>
      filter(!is.na(IRG), !is.na(NDVI)) |>
      mutate(alpage = ap)
  }
  
  ## 3. Concaténation ----------------------------------------------------------
  dataset_all <- map_dfr(alpages, load_one)
  
  ## 4. Sauvegarde sous le nom “Alpe‑Sud” --------------------------------------
  out_dir <- file.path(output_dir, "10. NDVI Effect",
                       paste0("dataset_", combo_name))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_rds <- file.path(out_dir,
                       sprintf("dataset_pheno_LONG_%d_%s.rds", year, combo_name))
  saveRDS(dataset_all, out_rds)
  message("✓ Sauvegardé : ", out_rds)
  
  if (write_fst) {
    if (!requireNamespace("fst", quietly = TRUE))
      stop("Installez le package {fst} ou mettez write_fst = FALSE.")
    fst::write_fst(dataset_all, sub("\\.rds$", ".fst", out_rds))
  }
  
  invisible(dataset_all)
}







build_smod_fsca_stack <- function(
    alpage, YEAR,
    smod_hydro_tif, out_dir, up_shape = NULL, outfile = NULL,
    smooth = TRUE, kernel = c("gauss","median","mean","circle"),
    radius_m = 60, round_DOY = TRUE,
    buffer_m = NULL
){
  suppressPackageStartupMessages(library(terra))
  kernel   <- match.arg(kernel)
  YEAR_chr <- as.character(YEAR)
  
  if (!file.exists(smod_hydro_tif)) stop("Fichier SMOD introuvable : ", smod_hydro_tif)
  dir.create(out_dir, TRUE, FALSE)
  if (is.null(outfile))
    outfile <- file.path(out_dir, paste0("SMOD_FSCA_stack_", alpage, "_", YEAR_chr, "_epsg2154.tif"))
  if (is.null(buffer_m)) buffer_m <- radius_m
  
  ## 1) Lecture + reprojection
  smod_hydro      <- rast(smod_hydro_tif)
  smod_hydro_2154 <- project(smod_hydro, "EPSG:2154", method = "near")
  
  ## 2) Buffer autour de l'UP (pas de masque ici, donc pas de trous)
  up_vect_2154 <- NULL
  smod_clip    <- smod_hydro_2154
  if (!is.null(up_shape) && file.exists(up_shape)) {
    up_vect_2154 <- project(vect(up_shape), "EPSG:2154")
    smod_clip    <- crop(smod_hydro_2154, buffer(up_vect_2154, width = buffer_m))
  }
  
  ## 3) DOY hydro -> DOY calendrier (CIRCULAIRE => pas de NA)
  jan1        <- as.Date(sprintf("%s-01-01", YEAR_chr))
  ndays       <- as.integer(as.Date(sprintf("%s-12-31", YEAR_chr)) - jan1) + 1L
  start_hydro <- as.Date(sprintf("%d-09-01", as.integer(YEAR_chr) - 1))
  
  hydro2cal_wrap <- function(x){
    # x en 1..365/366 —> mappe dans 1..ndays SANS NA
    d <- start_hydro + (as.integer(x) - 1L)
    # décalage circulaire : ((d - jan1) %% ndays) + 1
    k <- as.integer(d - jan1) %% ndays
    k + 1L
  }
  smod_cal <- app(smod_clip, hydro2cal_wrap)
  
  ## 4) LISSAGE (zéro NA en entrée => zéro NA en sortie)
  smod_smooth <- smod_cal
  if (smooth){
    px <- res(smod_cal)[1]
    if (kernel == "gauss"){
      sigma_m <- max(radius_m/3, px/2); size_m <- max(radius_m, 3*sigma_m)
      K <- focalMat(smod_cal, d = c(sigma_m, size_m), type = "Gauss")
      K <- K / sum(K, na.rm = TRUE)
      smod_smooth <- focal(smod_cal, w = K, fun = "sum", pad = TRUE)
    } else if (kernel == "median"){
      w <- 2*ceiling(radius_m/px) + 1
      smod_smooth <- focal(smod_cal, w = w, fun = median, pad = TRUE)
    } else if (kernel == "mean"){
      w <- 2*ceiling(radius_m/px) + 1
      W <- matrix(1, w, w); W <- W / sum(W)
      smod_smooth <- focal(smod_cal, w = W, fun = "sum", pad = TRUE)
    } else { # circle
      W <- focalMat(smod_cal, d = radius_m, type = "circle")
      W <- W / sum(W, na.rm = TRUE)
      smod_smooth <- focal(smod_cal, w = W, fun = "sum", pad = TRUE)
    }
  }
  
  ## 5) Arrondi + bornage (TRONQUER, pas de NA)
  if (round_DOY) smod_smooth <- round(smod_smooth)
  smod_smooth <- clamp(smod_smooth, lower = 1, upper = ndays, values = FALSE)
  
  ## 6) Masque final à l’UP (optionnel) — l’intérieur de l’UP reste **plein**
  smod_final <- if (!is.null(up_vect_2154)) mask(smod_smooth, up_vect_2154) else smod_smooth
  names(smod_final) <- paste0("SMOD_", YEAR_chr)
  
  ## 7) FSCA (rang) pour info, calculé sur l’UP
  v  <- values(smod_final)[,1]
  fs <- rep(NA_real_, length(v))
  idx <- which(!is.na(v))
  if (length(idx) > 1L){
    rnk <- rank(v[idx], ties.method = "min"); N <- length(idx)
    fs[idx] <- (N - rnk) / (N - 1)
  }
  fsca_r <- setValues(smod_final, fs)
  names(fsca_r) <- paste0("FSCA_", YEAR_chr)
  
  ## 8) Écriture
  stack_out <- c(smod_final, fsca_r)
  writeRaster(stack_out, outfile, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  invisible(outfile)
}
