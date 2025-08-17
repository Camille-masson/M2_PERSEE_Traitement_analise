# ============================================================
# build_compact_dataset()
# Crée un dataset compact multi-années par pixel :
#  - Colonnes : cell, x, y, YEAR, DAH, habitat_code[, habitat_label],
#               chargement_median_<range>, SMOD, MAXV, GPROD
#  - Aligne tout sur la template (résolution la + grossière parmi sources)
#  - PPI_extra : ne garde que bandes 02 (MAXV) et 18 (GPROD)
#  - Chargement médian : NA -> 0
#  - Ecrit un .rds et retourne le data.frame (invisible)
# ============================================================
build_multi_year_dataset <- function(
    alpage,
    years,
    raster_dir,
    output_dir,
    median_range_tag = "2022-2024",
    out_tag = NULL,         # par défaut: "<min>-<max>"
    quiet = FALSE
) {
  stopifnot(length(alpage) == 1, length(years) >= 1)
  suppressPackageStartupMessages({ library(terra) })
  
  msg <- function(...) if (!quiet) message(...)
  
  # --- ENTREE : chemins ---
  UP_file     <- file.path(raster_dir, "UP",   paste0("UP_",  alpage, ".shp"))
  DAH_tif     <- file.path(raster_dir, "Alti", paste0("DAH_1_", alpage, ".tif"))
  habitat_tif <- file.path(raster_dir, "Classifications_fusion_ColorIndexed_sc1_landforms_mnh.tif")
  lut_habitat_csv <- file.path(raster_dir, "class_habitat.csv")
  
  median_tif <- file.path(output_dir, "5. Indicateurs_visualisation", "Chargement_representatif",
                          alpage, median_range_tag,
                          paste0("chargement_median_", alpage, "_", median_range_tag, ".tif"))
  
  smod_dir  <- file.path(output_dir, "8. Analysis_Climate", "SMOD")
  pheno_dir <- file.path(raster_dir, "Phenologie")
  smod_path  <- function(Y) file.path(smod_dir,  paste0("SMOD_FSCA_stack_", alpage, "_", Y, "_epsg2154.tif"))
  extra_path <- function(Y) file.path(pheno_dir, paste0("PPI_extra_",       alpage, "_", Y, ".tif"))
  
  # --- SORTIE ---
  if (is.null(out_tag)) out_tag <- paste0(min(years), "-", max(years))
  out_dir <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_", alpage))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_rds <- file.path(out_dir, paste0("dataset_legacy_multi-year_", alpage, "_", out_tag, ".rds"))
  
  # --- 1) Vérifs fichiers ---
  need_static <- c(UP_file, DAH_tif, habitat_tif, median_tif)
  miss_stat <- need_static[!file.exists(need_static)]
  if (length(miss_stat))
    stop("Fichiers manquants :\n", paste(" •", miss_stat, collapse = "\n"))
  
  need_years <- c(vapply(years, smod_path, ""), vapply(years, extra_path, ""))
  miss_years <- need_years[!file.exists(need_years)]
  if (length(miss_years))
    stop("Fichiers annuels manquants :\n", paste(" •", miss_years, collapse = "\n"))
  
  # --- 2) Template (coarsest) comme ta fonction ---
  up_vect <- vect(UP_file)
  meta_orig <- list(
    DAH    = rast(DAH_tif),
    MEDIAN = rast(median_tif),
    SMFS   = rast(smod_path(years[1])),        # stack 2 bandes
    EXTRA  = rast(extra_path(years[1]))[[1]]
  )
  base_res <- res(meta_orig[[ which.max(sapply(meta_orig, function(r) prod(res(r))) ) ]])
  template <- rast(ext = ext(up_vect), resolution = base_res, crs = crs(up_vect))
  msg(sprintf("[i] Template: %.2f x %.2f m", base_res[1], base_res[2]))
  
  # --- Helpers d'alignement ---
  read_align <- function(file, method = "bilinear") {
    r <- rast(file) |>
      project(template, method = method) |>
      crop(up_vect)
    if (!hasValues(r)) {
      if (nlyr(r) == 1) r <- setValues(r, rep(NA_real_, ncell(r)))
      else r <- setValues(r, matrix(NA_real_, nrow = ncell(r), ncol = nlyr(r)))
    }
    r <- mask(r, up_vect)
    r
  }
  
  # spéciale EXTRA : crop dans CRS source + fallback NAflag
  read_align_extra <- function(file, method = "bilinear") {
    r0 <- rast(file)
    up_src <- project(up_vect, crs(r0))
    r0c <- try(crop(r0, up_src), silent = TRUE)
    if (inherits(r0c, "try-error") || ncell(r0c) == 0) {
      na1 <- template * NA
      if (nlyr(r0) > 1) na_stack <- rast(replicate(nlyr(r0), na1, simplify = FALSE)) else na_stack <- na1
      names(na_stack) <- if (!is.null(names(r0))) names(r0) else paste0("extra", seq_len(nlyr(na_stack)))
      return(na_stack)
    }
    r1 <- project(r0c, template, method = method) |> crop(up_vect)
    if (!hasValues(r1)) {
      nf <- try(NAflag(r0), silent = TRUE)
      if (!inherits(nf, "try-error") && !is.na(nf)) {
        NAflag(r0) <- NA
        r0c2 <- crop(r0, up_src)
        r1 <- project(r0c2, template, method = method) |> crop(up_vect)
      }
    }
    if (!hasValues(r1)) {
      if (nlyr(r1) == 1) r1 <- setValues(r1, rep(NA_real_, ncell(r1))) else
        r1 <- setValues(r1, matrix(NA_real_, nrow = ncell(r1), ncol = nlyr(r1)))
    }
    r1 <- mask(r1, up_vect)
    r1
  }
  
  # --- 3) Statique : DAH + HABITAT + médian (NA->0) ---
  dah_r     <- read_align(DAH_tif,     "bilinear"); names(dah_r)     <- "DAH"
  habitat_r <- read_align(habitat_tif, "near");     names(habitat_r) <- "habitat_code"
  
  median_r  <- read_align(median_tif,  "bilinear")
  median_r  <- ifel(is.na(median_r), 0, median_r)  # NA -> 0
  med_col_name <- paste0("chargement_median_", median_range_tag)
  names(median_r) <- med_col_name
  
  base_df <- as.data.frame(c(dah_r, habitat_r, median_r), cells = TRUE, xy = TRUE, na.rm = FALSE)
  base_df$cell <- as.integer(base_df$cell)
  base_df <- base_df[!is.na(base_df$DAH), , drop = FALSE]
  
  # habitat label via LUT (longueur identique via match)
  base_df$habitat_code <- as.integer(round(base_df$habitat_code))
  if (file.exists(lut_habitat_csv)) {
    lut <- try(read.csv2(lut_habitat_csv, stringsAsFactors = FALSE), silent = TRUE)
    if (!inherits(lut, "try-error") && all(c("code","label") %in% names(lut))) {
      lut$code <- as.integer(lut$code)
      idx <- match(base_df$habitat_code, lut$code)
      base_df$habitat_label <- lut$label[idx]
    }
  }
  
  # --- 4) Boucle années : SMOD + EXTRA (MAXV & GPROD) ---
  out_list <- vector("list", length(years))
  
  for (i in seq_along(years)) {
    Y <- years[i]
    
    # SMOD
    smfs   <- read_align(smod_path(Y), "near")
    smod_r <- smfs[[1]]; names(smod_r) <- "SMOD"
    
    # EXTRA (ne garder que MAXV et GPROD)
    extra_al <- read_align_extra(extra_path(Y), "bilinear")
    nm_src   <- names(rast(extra_path(Y)))
    
    i_maxv_byname  <- if (length(nm_src)) which(grepl("^MAXV$", nm_src, ignore.case = TRUE)) else integer(0)
    i_gprod_byname <- if (length(nm_src)) which(grepl("AUCg[_-]?ONSET[_-]?10[_-]?MAXD", nm_src, ignore.case = TRUE) |
                                                  grepl("GPROD", nm_src, ignore.case = TRUE)) else integer(0)
    
    i_maxv  <- if (length(i_maxv_byname)) i_maxv_byname[1] else if (nlyr(extra_al) >= 2) 2 else NA_integer_
    i_gprod <- if (length(i_gprod_byname)) i_gprod_byname[1] else if (nlyr(extra_al) >= 18) 18 else NA_integer_
    
    maxv_r  <- if (is.na(i_maxv))  (template * NA) else extra_al[[i_maxv]]
    gprod_r <- if (is.na(i_gprod)) (template * NA) else extra_al[[i_gprod]]
    names(maxv_r)  <- "MAXV"
    names(gprod_r) <- "GPROD"
    
    ann_df <- as.data.frame(c(smod_r, maxv_r, gprod_r), cells = TRUE, na.rm = FALSE)
    ann_df$cell <- as.integer(ann_df$cell)
    
    dfY <- merge(base_df, ann_df, by = "cell", all.x = TRUE, sort = FALSE)
    dfY$YEAR <- Y
    
    ord <- c("cell","x","y","YEAR","DAH","habitat_code",
             if ("habitat_label" %in% names(dfY)) "habitat_label",
             med_col_name,"SMOD","MAXV","GPROD")
    ord <- ord[ord %in% names(dfY)]
    dfY <- dfY[, c(ord, setdiff(names(dfY), ord)), drop = FALSE]
    
    out_list[[i]] <- dfY
    msg(sprintf("  - %d : %d px (NA MAXV=%s%%, NA GPROD=%s%%)",
                Y, nrow(dfY),
                round(mean(is.na(dfY$MAXV))*100,1),
                round(mean(is.na(dfY$GPROD))*100,1)))
  }
  
  # --- 5) Empilement & export ---
  dataset <- do.call(rbind, out_list)
  saveRDS(dataset, out_rds)
  msg("✓ Dataset compact écrit : ", out_rds)
  
  invisible(dataset)
}

























# ============================================================
# diag_dataset_compact()
# Diagnostic visuel pour le dataset "compact" (build_compact_dataset)
# - Détecte automatiquement la colonne du chargement médian : ^chargement_median_
# - Graphs : emprise, histos numériques, barplot SMOD, focus charge médiane (>=1)
# - Ecrit un PDF dans output_dir, retourne le chemin du PDF (invisible)
# ============================================================
diag_dataset_legacy_multi_year<- function(
    rds_path,
    up_shape   = NULL,
    sample_n   = 5e5,
    bins       = 60,
    output_dir = dirname(rds_path),
    palette    = c(
      DAH   = "#e7298a",
      SMOD  = "#66a61e",
      MAXV  = "#377eb8",
      GPROD = "#8dd3c7",
      CHARGEMENT_MEDIAN = "tomato"   # préfixe pour la colonne chargement_median_*
    ),
    view  = FALSE,
    quiet = FALSE
) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2)
    library(sf);    library(scales); library(grid)
    library(conflicted)
  })
  conflict_prefer("intersect", "base", quiet = TRUE)
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  msg    <- function(...) if (!quiet) message(...)
  
  stopifnot(file.exists(rds_path))
  df_full <- readRDS(rds_path)
  if (!is.data.frame(df_full))
    stop("Le fichier rds ne contient pas un data.frame.")
  
  cnames <- names(df_full); nrows <- nrow(df_full)
  
  # Trouver la colonne du chargement médian (ex. "chargement_median_2022-2024")
  median_col <- grep("^chargement_median_", cnames, value = TRUE)
  if (length(median_col) == 0) {
    warning("Aucune colonne 'chargement_median_*' détectée.")
    median_col <- NA_character_
  } else if (length(median_col) > 1) {
    warning("Plusieurs colonnes 'chargement_median_*' détectées, on prend la première : ", median_col[1])
    median_col <- median_col[1]
  }
  
  # Petites stats console
  if (!quiet) {
    msg("Colonnes : ", paste(cnames, collapse = ", "))
    msg("Lignes : ", format(nrows, big.mark = " "))
    num_vars <- cnames[sapply(df_full, is.numeric)]
    for (v in num_vars) {
      vals <- df_full[[v]]
      msg(sprintf("  %s : mean=%.3f, median=%.3f, sd=%.3f, max=%.3f, NA=%d",
                  v, mean(vals, na.rm=TRUE), median(vals, na.rm=TRUE),
                  sd(vals, na.rm=TRUE), max(vals, na.rm=TRUE), sum(is.na(vals))))
    }
  }
  
  # Mappage couleurs (insensible à la casse + gestion préfixe chargement_median)
  names(palette) <- toupper(names(palette))
  pcol <- function(v) {
    vn <- toupper(v)
    if (vn %in% names(palette)) return(palette[[vn]])
    # préfixe CHARGEMENT_MEDIAN pour n'importe quelle colonne chargement_median_*
    if (grepl("^CHARGEMENT_MEDIAN", vn)) return(palette[["CHARGEMENT_MEDIAN"]] %||% "tomato")
    if (grepl("^SMOD", vn)) return(palette[["SMOD"]] %||% "#66a61e")
    NULL
  }
  
  ## 1) Étendue géographique
  have_up <- !is.null(up_shape) && file.exists(up_shape)
  if (have_up) {
    up_sf <- sf::st_read(up_shape, quiet = TRUE)
    bbox  <- sf::st_bbox(up_sf) + c(-50,-50,50,50)
  } else {
    if (!all(c("x","y") %in% cnames))
      stop("Sans up_shape, les colonnes x et y doivent exister.")
    bbox <- sf::st_bbox(c(xmin = min(df_full$x, na.rm=TRUE),
                          ymin = min(df_full$y, na.rm=TRUE),
                          xmax = max(df_full$x, na.rm=TRUE),
                          ymax = max(df_full$y, na.rm=TRUE)),
                        crs = sf::st_crs(2154)) + c(-50,-50,50,50)
  }
  poly_bbox <- sf::st_as_sfc(bbox)
  
  ## 2) Échantillonnage
  rows <- if (nrows <= sample_n) seq_len(nrows) else sample(nrows, sample_n)
  keep_cols <- c("x","y","YEAR","DAH","habitat_code","habitat_label","SMOD","MAXV","GPROD", median_col)
  keep_cols <- unique(keep_cols[keep_cols %in% cnames])
  df <- df_full[rows, keep_cols, drop = FALSE]
  
  ## 3) Graphiques
  plots <- list()
  
  # (a) carte d'emprise (optionnellement avec l'UP)
  p_map <- ggplot() +
    { if (have_up) geom_sf(data = up_sf, fill = "grey85",
                           colour = "steelblue", linewidth = .5) } +
    geom_sf(data = poly_bbox, fill = NA, colour = "red3", linewidth = 1) +
    coord_sf(expand = FALSE) +
    labs(title = "Emprise du jeu de données",
         subtitle = basename(rds_path)) +
    theme_void()
  plots[[length(plots) + 1]] <- p_map
  
  # (b) Histogrammes des numériques standard
  num_vars <- names(df)[sapply(df, is.numeric)]
  num_vars <- setdiff(num_vars, "SMOD")  # SMOD traité à part
  for (v in num_vars) {
    fill_col <- pcol(v) %||% "grey60"
    plots[[length(plots)+1]] <-
      ggplot(df, aes(.data[[v]])) +
      geom_histogram(bins = bins, fill = fill_col, colour = "grey25") +
      labs(title = paste("Distribution de", v),
           subtitle = paste("Échantillon", format(nrow(df), big.mark=" ")),
           x = v, y = "Fréquence") +
      theme_minimal(base_size = 9)
  }
  
  # (c) SMOD : barplot (si entier/catégoriel)
  if ("SMOD" %in% names(df)) {
    smod_df <- df %>% mutate(SMOD_class = factor(SMOD))
    plots[[length(plots)+1]] <-
      ggplot(smod_df, aes(SMOD_class)) +
      geom_bar(fill = pcol("SMOD") %||% "#66a61e") +
      labs(title = "Classes SMOD (toutes années confondues)",
           x = "SMOD", y = "Pixels") +
      theme_minimal(base_size = 9)
  }
  
  # (d) Focus chargement médian : bar Chargé / Non chargé + histo log10 (>0)
  if (!is.na(median_col) && median_col %in% names(df)) {
    vC <- median_col
    df_bar <- df %>%
      mutate(etat = ifelse(is.finite(.data[[vC]]) & .data[[vC]] >= 1,
                           "Chargé (≥1)", "Non chargé (<1)"))
    
    plots[[length(plots)+1]] <-
      ggplot(df_bar, aes(etat)) +
      geom_bar(fill = pcol("CHARGEMENT_MEDIAN") %||% "tomato") +
      labs(title="Pixels chargés / non chargés (chargement médian)",
           x=NULL, y="Pixels") +
      theme_minimal(base_size = 9)
    
    df_pos <- dplyr::filter(df_bar, .data[[vC]] > 0)
    if (nrow(df_pos) > 0) {
      plots[[length(plots)+1]] <-
        ggplot(df_pos, aes(.data[[vC]])) +
        geom_histogram(bins = bins,
                       fill = pcol("CHARGEMENT_MEDIAN") %||% "tomato",
                       colour = "grey25") +
        scale_x_log10(labels = scales::comma) +
        labs(title     = "Distribution du chargement médian (>0)",
             subtitle  = paste("Échantillon", format(nrow(df_pos), big.mark=" ")),
             x = paste0(vC, " (log10)"), y = "Fréquence") +
        theme_minimal(base_size = 9)
    }
  }
  
  ## 4) Export PDF
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(rds_path)), "_diag_compact.pdf")
  )
  grDevices::pdf(pdf_path, width = 7, height = 5)
  for (p in plots) print(p)
  grDevices::dev.off()
  msg("✓ PDF écrit :", pdf_path)
  
  if (view) for (p in plots) { dev.new(7,5); print(p) }
  invisible(pdf_path)
}













# ============================================================
# merge_compact_dataset_years()
# Fusionne les datasets compacts multi-années (déjà agrégés) de plusieurs alpages.
# - Cherche pour chaque alpage un "dataset_compact_<ap>_<min>-<max>.rds"
#   qui correspond exactement à range(years); sinon un fichier qui couvre ce range.
# - Filtre YEAR %in% years, ajoute 'alpage', aligne les colonnes et concatène.
# - Ecrit: outputs/10. NDVI Effect/dataset_<combo_name>/
#          dataset_compact_<combo_name>_<min>-<max>.rds
# ============================================================
merge_compact_dataset_years <- function(
    years,
    alpages,
    output_dir,
    combo_name = "Alpe-Sud",
    write_fst  = FALSE,
    quiet      = FALSE
) {
  stopifnot(length(years) >= 1, length(alpages) >= 1, dir.exists(output_dir))
  alpages <- sort(unique(alpages))
  yr_min  <- min(years); yr_max <- max(years)
  
  msg <- function(...) if (!quiet) message(...)
  
  pick_compact_file <- function(ap) {
    base_dir <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_", ap))
    if (!dir.exists(base_dir))
      stop("Dossier introuvable : ", base_dir, call. = FALSE)
    
    files <- list.files(base_dir,
                        pattern = sprintf("^dataset_legacy_multi-year_%s_.*\\.rds$", ap),
                        full.names = TRUE)
    if (length(files) == 0)
      stop("Aucun fichier compact trouvé pour '", ap, "' dans ", base_dir, call. = FALSE)
    
    # extraire l'intervalle yyyy-yyyy de chaque fichier
    get_range <- function(fn) {
      m <- regmatches(basename(fn),
                      regexpr("(\\d{4})-(\\d{4})", basename(fn)))
      if (length(m) == 0) return(c(NA_integer_, NA_integer_))
      as.integer(strsplit(m, "-")[[1]])
    }
    ranges <- t(vapply(files, get_range, integer(2)))
    colnames(ranges) <- c("min","max")
    
    # 1) priorité à l'intervalle EXACT
    exact <- which(ranges[, "min"] == yr_min & ranges[, "max"] == yr_max)
    if (length(exact) >= 1) return(files[exact[1]])
    
    # 2) sinon, un fichier qui couvre [yr_min, yr_max]
    covers <- which(ranges[, "min"] <= yr_min & ranges[, "max"] >= yr_max)
    if (length(covers) >= 1) {
      # prendre celui au "range" le plus serré
      widths <- ranges[covers, "max"] - ranges[covers, "min"]
      return(files[covers[which.min(widths)]])
    }
    
    stop("Aucun dataset compact pour '", ap, "' ne couvre [", yr_min, "-", yr_max,
         "]. Fichiers trouvés :\n  - ", paste(basename(files), collapse = "\n  - "), call. = FALSE)
  }
  
  load_one <- function(ap) {
    f <- pick_compact_file(ap)
    df <- readRDS(f)
    if (!is.data.frame(df))
      stop("Fichier corrompu (pas un data.frame) : ", f, call. = FALSE)
    if (!"YEAR" %in% names(df))
      stop("La colonne YEAR est absente dans : ", f, call. = FALSE)
    
    dfi <- df[df$YEAR %in% years, , drop = FALSE]
    if (nrow(dfi) == 0)
      msg(sprintf("[!] Aucun pixel YEAR ∈ [%d-%d] dans %s", yr_min, yr_max, basename(f)))
    
    dfi$alpage <- ap
    msg(sprintf("  - %s : %s | lignes gardées = %d",
                ap, basename(f), nrow(dfi)))
    dfi
  }
  
  lst <- lapply(alpages, load_one)
  
  # rbind avec remplissage des colonnes manquantes
  all_cols <- unique(unlist(lapply(lst, names)))
  lst2 <- lapply(lst, function(d) {
    miss <- setdiff(all_cols, names(d))
    for (m in miss) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  merged <- do.call(rbind, lst2)
  
  # sortie
  out_dir <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_", combo_name))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_rds <- file.path(out_dir, sprintf("dataset_legacy_multi-year_%s_%d-%d.rds", combo_name, yr_min, yr_max))
  saveRDS(merged, out_rds)
  msg("✓ Sauvegardé : ", out_rds)
  
  if (write_fst) {
    if (!requireNamespace("fst", quietly = TRUE))
      stop("Installez le package {fst} ou mettez write_fst = FALSE.")
    fst::write_fst(merged, sub("\\.rds$", ".fst", out_rds))
    msg("✓ Sauvegardé (fst) : ", sub("\\.rds$", ".fst", out_rds))
  }
  
  invisible(merged)
}
