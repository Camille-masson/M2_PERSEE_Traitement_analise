load_and_perc_by_vegetation_period <- function(pheno, charge) {
        charge_by_season <- data.frame(pheno[, c("x", "y")])
        charge_by_season$growing = 0
        charge_by_season$plateau = 0
        charge_by_season$senesc = 0
        for (d in unique(charge$day)) {
            charge_jour = charge[charge$day==d, ]
            growing = (d >= pheno$t0 & d <= pheno$t1)
            plateau = (d >= pheno$t1 & d <= pheno$t2)
            charge_by_season$growing = charge_by_season$growing + charge_jour$Charge * growing
            charge_by_season$plateau = charge_by_season$plateau + charge_jour$Charge * plateau
            charge_by_season$senesc = charge_by_season$senesc + charge_jour$Charge * (!growing & !plateau)
        }
        charge_by_season$total = charge_by_season$growing + charge_by_season$plateau + charge_by_season$senesc

        charge_by_season = charge_by_season[charge_by_season$total > 0, ]

        charge_by_season$perc_growing = 100 * charge_by_season$growing / charge_by_season$total
        charge_by_season$perc_plateau = 100 * charge_by_season$plateau / charge_by_season$total
        charge_by_season$perc_senesc <- 100 - charge_by_season$perc_growing - charge_by_season$perc_plateau

        return(charge_by_season)
}

get_charge_by_period <- function(charge, breakdays) {
    # Computes the flock load by periods
    # INPUTS
    #    charge : a data.frame with x, y, day and Charge columns
    #    breakdays : the julian days defining the begining and end of periods (the first period being anything before the first breakday,
    #                 the last on anithing after the last breakday)
    # OUPUTS
    #    a data.frame with x, y, day and period columns, periods going from 1 to length(breakdays)+1

    charge %>%
        mutate(period = findInterval(day, breakdays)+1) %>%
        group_by(x, y, period) %>%
        summarize(charge = sum(Charge)) %>%
        as.data.frame() %>%
    return()
}

get_chargement_by_period_and_vegetation_units <- function(charge_by_period, vegetation_units, period_labs) {
    # Computes the number of sheep.days per vegetation units and per period
    # INPUTS
    #    charge_by_period : a data.frame with x, y, day and period columns, periods ranging from 1 to length(period_labs)
    #    vegetation_units : a terra SpatVector containing the vegetation polygons, with a vegetation_type and attribute
    #    period_labs : the name of the different periods, ordered by time
    # OUPUTS
    #    a data.frame with x, y, period, vegetation_type, charge and area (of polygon) columns 
    pixel_surface = get_pixel_surface(charge_by_period)
    test = charge_by_period %>%
        terra::rast() %>%
        terra::extract(y= vegetation_units, fun = sum, na.rm = T, ID = F) %>%
        mutate(vegetation_type=vegetation_units$vegetation_type) %>%
        mutate(area=terra::expanse(vegetation_units, unit="ha")) %>%
        pivot_longer(cols=3:ncol(charge_by_period)-2,
                        names_to='period',
                        values_to='charge') %>%
        mutate(period = period_labs[as.numeric(period)]) %>%
        mutate(period = factor(period, levels = period_labs[length(period_labs):1])) %>%
        mutate(charge = charge*pixel_surface/10000) %>%
        as.data.frame() %>%
    return()
}



###  FIGURES ###
#**************#

plot_presence_perc_by_period_and_habitat <- function(data, period_title, period_labs, habitat_labs, title="") {
    data = data %>%
                count(vegetation_type, period, wt = charge, name = "charge")
    print(ggplot(data, aes(y = vegetation_type, x = charge, fill = period)) +
        geom_col(position = "fill") +
        scale_x_continuous(labels = scales::percent_format()) +
        scale_fill_viridis_d(option = "viridis", direction = -1, period_title,
                            breaks=levels(data$period)[nlevels(data$period):1], labels=period_labs) +
        xlab("Flock presence") +
        ylab("") +
        scale_y_discrete(breaks=habitat_labs$vegetation_type, labels=habitat_labs$lab) +
        ggtitle(title))
}

plot_charge_by_period_and_habitat <- function(data, period_title, period_labs, habitat_labs, title="") {
        print(ggplot(data, aes(y = vegetation_type, x = charge/area, fill = period)) +
        geom_boxplot(outlier.shape = NA, varwidth = F) +
        scale_fill_viridis_d(option = "viridis", direction = -1, period_title,
                            breaks=levels(data$period)[nlevels(data$period):1]) +
        xlab("Flock load (sheep.days/ha)") +
        ylab("") +
        scale_y_discrete(breaks=habitat_labs$vegetation_type, labels=habitat_labs$lab) +
        ggtitle(title) +
        coord_cartesian(xlim = c(0, quantile(data$charge/data$area, 0.98, na.rm=T)))) # contrairement à scale_x_continuous(xlim =...), coord_cartesian ne supprime pas les données hors graphe avant de calculer les statistiques
}





















# Fonction permettant de généré un rasteur du pourcentage de présence du troupeau
# en fonction des 4 stade phénologique identifié dans le script ressource : Plant_Phenology_Index_processing




rast_pourcentage_troupeau_phase_pheno <- function (daily_rds_file, pheno_tif_file, out_file){
  
  # —————————————————————————————————————————————————————————————
  # 1) Lecture et préparation des présences journalières
  # —————————————————————————————————————————————————————————————
  df <- readRDS(daily_rds_file)
  
  pres_df <- df %>%
    group_by(x, y, day) %>%
    summarise(totalCharge = sum(Charge, na.rm = TRUE), .groups = "drop") %>%
    mutate(present = as.integer(totalCharge > 1))
  
  # —————————————————————————————————————————————————————————————
  # 2) Chargement du raster de phénologie
  # —————————————————————————————————————————————————————————————
  pheno <- rast(pheno_tif_file)
  crs(pheno) <- "EPSG:2154"
  
  # extraire les DOY depuis "Phase_DOYXXX"
  doys <- as.integer(sub("Phase_DOY", "", names(pheno)))
  
  # —————————————————————————————————————————————————————————————
  # 3) Rasterisation ALIGNÉE : un layer/jour pour la présence
  # —————————————————————————————————————————————————————————————
  # On crée un cube 0/1 calé sur pheno
  base0 <- pheno[[1]] * 0
  pres_stack <- rep(base0, length(doys))
  names(pres_stack) <- paste0("day", doys)
  
  # Peupler avec les présences
  for(i in seq_along(doys)) {
    d <- doys[i]
    sub <- pres_df %>% filter(day == d)
    if(nrow(sub)==0) next
    vsub <- vect(sub, geom=c("x","y"), crs=crs(pheno))
    # rasterize sur la couche pheno[[i]] pour garder l'alignement
    r <- rasterize(vsub,
                   pheno[[i]],
                   field="present",
                   fun="max",
                   background=0)
    pres_stack[[i]] <- r
  }
  
  # —————————————————————————————————————————————————————————————
  # 4) Comptages et % par classe phénologique
  # —————————————————————————————————————————————————————————————
  classes <- c(pousse=1, plateau=2, deperiss=3, senescence=4)
  
  # (a) compter jours de présence PAR CLASSE
  counts <- lapply(classes, function(cl) {
    mask_cl <- pheno == cl         # binaire par jour, par pixel
    # multiplication couche à couche → 1 si (présent ET good class), sinon 0
    cst <- sum(mask_cl * pres_stack, na.rm=TRUE)
    names(cst) <- paste0("count_cl",cl)
    cst
  })
  
  # (b) total de jours présent (toutes classes)
  total_pres <- sum(pres_stack, na.rm=TRUE)
  
  # (c) % par classe = 100 * count_cl / total_pres
  percs <- mapply(function(r_count, nm){
    p <- (r_count / total_pres) * 100
    p[total_pres == 0] <- NA
    names(p) <- nm
    p
  }, counts, names(counts), SIMPLIFY=FALSE)
  
  # —————————————————————————————————————————————————————————————
  # 5) Écriture des TIFs finaux
  # —————————————————————————————————————————————————————————————
  for(nm in names(percs)) {
    writeRaster(percs[[nm]],
                filename  = out_files[[nm]],
                overwrite = TRUE)
  }
}


build_dataset_pheno_vec <- function(
    YEAR,
    alpage,
    pheno_tif_file,
    delta_max_tif_file,
    daily_rds_file,
    output_IRG_by_habitat,
    DOY_range   = 121:334,
    chunk_size  = 30,
    state_keep  = "Paturage"           # NULL → tous les états
){
  library(terra)
  library(data.table)
  library(pbapply)
  
  ## 1. Rasters ----------------------------------------------------------------
  ph  <- rast(pheno_tif_file)
  de  <- rast(delta_max_tif_file)
  
  stopifnot(ext(ph) == ext(de), res(ph) == res(de))
  
  avail_DOY <- as.integer(sub(".*DOY", "", names(ph)))
  layer_ids <- match(DOY_range, avail_DOY)
  if (anyNA(layer_ids))
    stop("DOY manquants : ",
         paste(DOY_range[is.na(layer_ids)], collapse = ", "))
  
  n_cells <- ncell(ph)
  
  ## 2. Coordonnées ------------------------------------------------------------
  coords <- as.data.table(crds(ph, df = TRUE))[, cell := .I]
  
  ## 3. Charges → vecteurs de lookup ------------------------------------------
  ch <- as.data.table(readRDS(daily_rds_file))
  setnames(ch,
           c("x","y","day","Charge","parc"),
           c("x","y","doy","charge","parc_nuit"))
  
  if (!is.null(state_keep) && "state" %in% names(ch))
    ch <- ch[state == state_keep]
  
  ch[, cell := terra::cellFromXY(ph, cbind(x, y))]
  ch <- ch[!is.na(cell)][ , doy := as.integer(doy)]
  
  ch <- ch[, .(
    charge    = sum(charge, na.rm = TRUE),
    parc_nuit = paste(unique(parc_nuit), collapse = ";")
  ), by = .(cell, doy)]
  
  ch[, id := paste(cell, doy, sep = ":")]
  charge_vec <- setNames(ch$charge,   ch$id)
  parc_vec   <- setNames(ch$parc_nuit,ch$id)
  
  ## 4. Boucle raster (sans melt) ---------------------------------------------
  res <- vector("list", ceiling(length(layer_ids)/chunk_size))
  pbapply::pboptions(type = "timer")
  
  for (chunk in seq_along(res)) {
    
    s <- (chunk-1)*chunk_size + 1
    e <- min(chunk*chunk_size, length(layer_ids))
    ly  <- layer_ids[s:e]
    dyy <- DOY_range[s:e]
    n_d <- length(dyy)
    
    delta_mat <- as.matrix(de[[ly]])
    pheno_mat <- as.matrix(ph[[ly]])
    
    ## Matrices (n_cells × n_d) → vecteurs ------------------------------------
    dt <- data.table(
      cell               = rep(seq_len(n_cells), times = n_d),
      doy                = rep(dyy, each = n_cells),
      delta_day_IRG_max  = as.vector(delta_mat),
      pheno_stage        = as.vector(pheno_mat)
    )
    
    ## Ajout charge / parc_nuit par lookup ------------------------------------
    dt[, id := paste(cell, doy, sep = ":")]
    dt[, charge    := charge_vec[id]]
    dt[, parc_nuit := parc_vec[id]]
    dt[, id := NULL]
    
    res[[chunk]] <- dt
  }
  
  ## 5. Assemblage & sauvegarde ------------------------------------------------
  dat <- rbindlist(res)
  dat <- merge(dat, coords, by = "cell", all.x = TRUE, sort = FALSE)
  
  setcolorder(dat,
              c("cell","x","y","doy",
                "delta_day_IRG_max","pheno_stage",
                "charge","parc_nuit"))
  
  saveRDS(dat, output_IRG_by_habitat, compress = "xz")
  message("✓ Dataset écrit : ", output_IRG_by_habitat,
          "\n   → ", format(object.size(dat), units = "auto"))
}















build_pheno_dataset <- function(alpage,
                                YEAR,
                                raster_dir,
                                output_dir,
                                DOY_range = 121:334,
                                lut_habitat_csv = file.path(raster_dir,
                                                            "class_habitat.csv"),
                                output_fst = NULL) {
  
  suppressPackageStartupMessages({
    library(terra)
    library(sf)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(fst)
    # (ligne conflicts_prefer supprimée)
  })
  
  # -------------------------------------------------------------------------
  # 0. Chemins d’accès -------------------------------------------------------
  phen_dir   <- file.path(raster_dir, "Phenologie")
  irg_file   <- file.path(phen_dir,
                          sprintf("IRG_season_%s_%d.tif",  alpage, YEAR))
  ndvi_file  <- file.path(phen_dir,
                          sprintf("NDVI_season_%s_%d.tif", alpage, YEAR))
  
  habitat_file <- file.path(raster_dir,
                            "Classifications_fusion_ColorIndexed_sc1_landforms_mnh.tif")
  dah_file     <- file.path(raster_dir, "Alti",
                            sprintf("DAH_1_%s.tif", alpage))
  fsca_file    <- file.path(output_dir, "8. Analysis_Climate", "SMOD",
                            sprintf("Fsca_%s.tif", alpage))
  
  charge_file  <- file.path(output_dir, "4. Chargements_Calcules",
                            sprintf("%d_%s", YEAR, alpage),
                            sprintf("total_%d_%s.rds", YEAR, alpage))
  
  # -------------------------------------------------------------------------
  # 1. Sortie ---------------------------------------------------------------
  if (is.null(output_fst)) {
    out_dir <- file.path(output_dir, "9. Analysis_Phenology", "data_joint")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    output_fst <- file.path(out_dir,
                            sprintf("dataset_pheno_LONG_%d_%s.fst", YEAR, alpage))
  }
  
  # -------------------------------------------------------------------------
  # 2. Lecture IRG / NDVI ----------------------------------------------------
  message("Lecture IRG  : ", irg_file)
  irg_stack  <- terra::rast(irg_file)
  
  message("Lecture NDVI : ", ndvi_file)
  ndvi_stack <- terra::rast(ndvi_file)
  
  if (!terra::compareGeom(irg_stack, ndvi_stack, stopOnError = FALSE))
    stop("IRG et NDVI ne partagent pas exactement la même grille.")
  
  template <- irg_stack[[1]]
  
  # -------------------------------------------------------------------------
  # 3. Couches statiques -----------------------------------------------------
  message("Lecture habitat, DAH, FSCA…")
  
  habitat_r <- terra::rast(habitat_file) |>
    terra::project(template, method = "near") |>
    terra::resample(template, "near")
  
  dah_r  <- terra::rast(dah_file)  |> terra::project(template) |> terra::resample(template, "bilinear")
  fsca_r <- terra::rast(fsca_file) |> terra::project(template) |> terra::resample(template, "bilinear")
  
  static_df <- terra::as.data.frame(
    c(dah_r, fsca_r, habitat_r),
    cells = TRUE, na.rm = FALSE
  ) |>
    rlang::set_names(c("cell", "dah", "fsca", "habitat_code"))
  
  habitat_lut <- readr::read_csv2(lut_habitat_csv,
                                  col_types = "ic",
                                  locale    = locale(encoding = "latin1")) |>
    dplyr::mutate(code = as.integer(code))
  
  if (!file.exists(charge_file))
    stop("Fichier chargement introuvable : ", charge_file)
  
  charge_df <- readRDS(charge_file) |> tibble::as_tibble()
  
  if (!"cell" %in% names(charge_df)) {
    message("↪ Ajout de 'cell' dans charge_df…")
    if (!all(c("x", "y") %in% names(charge_df)))
      stop("Le .rds doit contenir x et y pour calculer 'cell'.")
    charge_df <- charge_df |>
      dplyr::mutate(cell = terra::cellFromXY(template, cbind(x, y))) |>
      dplyr::filter(!is.na(cell))
  }
  
  charge_df <- charge_df |>
    dplyr::select(cell, Charge)   # plus de x / y
  
  # -------------------------------------------------------------------------
  # 4. Rasters journaliers → long -------------------------------------------
  message("Conversion IRG/NDVI en tables longues…")
  
  to_long <- function(r_stack, var) {
    terra::as.data.frame(r_stack, cells = TRUE, xy = TRUE, na.rm = FALSE) |>
      tidyr::pivot_longer(
        cols      = -c(cell, x, y),
        names_to  = "band",
        values_to = var
      ) |>
      dplyr::mutate(
        DOY = as.integer(stringr::str_extract(band, "\\d{3,}$"))
      ) |>
      dplyr::filter(DOY %in% DOY_range) |>
      dplyr::select(-band)
  }
  
  irg_df  <- to_long(irg_stack,  "IRG")
  ndvi_df <- to_long(ndvi_stack, "NDVI")
  
  pheno_df <- dplyr::inner_join(
    irg_df, ndvi_df,
    by = c("cell", "x", "y", "DOY")
  )
  
  rm(irg_df, ndvi_df) ; gc()
  
  # -------------------------------------------------------------------------
  # 5. Fusion finale ---------------------------------------------------------
  message("Fusion avec DAH, FSCA, habitat, chargement…")
  
  full_df <- pheno_df |>
    dplyr::left_join(static_df,  by = "cell") |>
    dplyr::left_join(charge_df,  by = "cell") |>
    dplyr::mutate(habitat_code = as.integer(habitat_code)) |>
    dplyr::left_join(habitat_lut, by = c("habitat_code" = "code")) |>
    dplyr::rename(habitat_label = label) |>
    dplyr::relocate(habitat_label, .after = habitat_code)
  
  # -------------------------------------------------------------------------
  # 6. Écriture .fst ---------------------------------------------------------
  message("Écriture du fichier FST : ", output_fst)
  fst::write_fst(full_df, output_fst, compress = 60)
  message("✓ Terminé : ", output_fst)
}














build_pheno_dataset <- function(alpage,
                                YEAR,
                                # chemins explicites ---------------------------------------------------
                                irg_tif,          # ".../IRG_season_<alp>_<yr>.tif"
                                ndvi_tif,         # ".../NDVI_season_<alp>_<yr>.tif"
                                habitat_tif,      # ".../Classifications_…tif"
                                dah_tif,          # ".../DAH_1_<alp>.tif"
                                fsca_tif,         # ".../Fsca_<alp>.tif"
                                charge_rds,       # ".../total_<yr>_<alp>.rds"
                                lut_habitat_csv,  # ".../class_habitat.csv"
                                # sortie ----------------------------------------------------------------
                                output_fst,
                                # options ---------------------------------------------------------------
                                DOY_range = 121:334) {
  
  suppressPackageStartupMessages({
    library(terra);  library(dplyr); library(tidyr)
    library(readr);  library(stringr); library(fst)
  })
  
  # --------------------------------------------------------------------- #
  # 0. Vérification des fichiers                                          #
  # --------------------------------------------------------------------- #
  needed <- c(irg_tif, ndvi_tif, habitat_tif,
              dah_tif, fsca_tif, charge_rds, lut_habitat_csv)
  miss <- needed[!file.exists(needed)]
  if (length(miss))
    stop("Fichiers manquants :\n", paste(" •", miss, collapse = "\n"), call. = FALSE)
  
  dir.create(dirname(output_fst), recursive = TRUE, showWarnings = FALSE)
  
  # --------------------------------------------------------------------- #
  # 1. Rasters dynamiques (IRG / NDVI)                                    #
  # --------------------------------------------------------------------- #
  message("Lecture IRG  : ", irg_tif)
  irg_stack  <- terra::rast(irg_tif)
  
  message("Lecture NDVI : ", ndvi_tif)
  ndvi_stack <- terra::rast(ndvi_tif)
  
  if (!terra::compareGeom(irg_stack, ndvi_stack, stopOnError = FALSE))
    stop("IRG et NDVI ne partagent pas la même grille.")
  
  template <- irg_stack[[1]]   # raster 10 m de référence
  
  # --------------------------------------------------------------------- #
  # 2. Rasters statiques (habitat, DAH, FSCA)                             #
  # --------------------------------------------------------------------- #
  message("Lecture habitat, DAH, FSCA…")
  
  habitat_r <- terra::rast(habitat_tif) |> 
    terra::project(template, "near") |> 
    terra::resample(template, "near")
  
  dah_r  <- terra::rast(dah_tif)  |> terra::project(template) |> terra::resample(template)
  fsca_r <- terra::rast(fsca_tif) |> terra::project(template) |> terra::resample(template)
  
  static_df <- terra::as.data.frame(
    c(dah_r, fsca_r, habitat_r),
    cells = TRUE, na.rm = FALSE
  ) |>
    rlang::set_names(c("cell", "dah", "fsca", "habitat_code"))
  
  # LUT habitat
  habitat_lut <- readr::read_csv2(lut_habitat_csv,
                                  col_types = "ic",
                                  locale    = locale(encoding = "latin1")) |>
    dplyr::mutate(code = as.integer(code))
  
  # --------------------------------------------------------------------- #
  # 3. Chargement pastoral                                                #
  # --------------------------------------------------------------------- #
  charge_df <- readRDS(charge_rds) |> tibble::as_tibble()
  
  if (!"cell" %in% names(charge_df)) {
    if (!all(c("x", "y") %in% names(charge_df)))
      stop("Le .rds doit contenir soit 'cell', soit 'x' et 'y'.")
    charge_df <- charge_df |>
      dplyr::mutate(cell = terra::cellFromXY(template, cbind(x, y))) |>
      dplyr::filter(!is.na(cell))
  }
  
  charge_df <- charge_df |>
    dplyr::select(cell, Charge)        # aucune collision possible
  
  # --------------------------------------------------------------------- #
  # 4. IRG / NDVI en format long                                          #
  # --------------------------------------------------------------------- #
  to_long <- function(stack, var) {
    terra::as.data.frame(stack, cells = TRUE, xy = TRUE, na.rm = FALSE) |>
      tidyr::pivot_longer(
        cols      = -c(cell, x, y),
        names_to  = "band",
        values_to = var
      ) |>
      dplyr::mutate(DOY = as.integer(stringr::str_extract(band, "\\d{3,}$"))) |>
      dplyr::filter(DOY %in% DOY_range) |>
      dplyr::select(-band)
  }
  
  irg_long  <- to_long(irg_stack,  "IRG")
  ndvi_long <- to_long(ndvi_stack, "NDVI")
  
  pheno_df  <- dplyr::inner_join(
    irg_long, ndvi_long,
    by = c("cell", "x", "y", "DOY")
  )
  
  rm(irg_long, ndvi_long) ; gc()
  
  # --------------------------------------------------------------------- #
  # 5. Fusion finale                                                      #
  # --------------------------------------------------------------------- #
  full_df <- pheno_df |>
    dplyr::left_join(static_df,  by = "cell") |>
    dplyr::left_join(charge_df,  by = "cell") |>
    dplyr::mutate(habitat_code = as.integer(habitat_code)) |>
    dplyr::left_join(habitat_lut, by = c("habitat_code" = "code")) |>
    dplyr::rename(habitat_label = label) |>
    dplyr::relocate(habitat_label, .after = habitat_code)
  
  n_pix <- dplyr::n_distinct(full_df$cell)
  message("Pixels uniques (≠ DOY) : ", format(n_pix, big.mark = " "))
  
  # --------------------------------------------------------------------- #
  # 6. Écriture .fst                                                      #
  # --------------------------------------------------------------------- #
  fst::write_fst(full_df, output_fst, compress = 60)
  message("✓ Fichier écrit : ", output_fst)
  
  invisible(full_df)
}



build_pheno_dataset <- function(alpage, YEAR,
                                irg_tif, ndvi_tif,
                                habitat_tif, dah_tif, fsca_tif,
                                charge_rds, lut_habitat_csv,
                                up_shape,            # UP EPSG:2154
                                output_fst,
                                DOY_range = 121:334) {
  
  suppressPackageStartupMessages({
    library(terra);  library(dplyr);  library(tidyr)
    library(readr);  library(stringr);library(fst)
  })
  
  `%||%` <- function(a,b) if(!is.null(a)) a else b
  needed <- c(irg_tif, ndvi_tif, habitat_tif, dah_tif, fsca_tif,
              charge_rds, lut_habitat_csv, up_shape)
  if (any(!file.exists(needed)))
    stop("Fichier manquant :\n", paste(needed[!file.exists(needed)], collapse="\n"))
  
  dir.create(dirname(output_fst), recursive=TRUE, showWarnings=FALSE)
  
  ## 1. UP et métadonnées raster ----------------------------------------
  up_vect_orig <- terra::vect(up_shape)              # EPSG 2154
  surf_up_km2  <- round(terra::expanse(up_vect_orig)/1e6, 2)
  message("Surface UP : ", surf_up_km2, " km²")
  
  # pixel surface pour chaque raster
  pix_surf <- \(r) prod(terra::res(r))
  
  meta <- list(
    IRG  = terra::rast(irg_tif)[[1]],
    NDVI = terra::rast(ndvi_tif)[[1]],
    HAB  = terra::rast(habitat_tif),
    DAH  = terra::rast(dah_tif),
    FSCA = terra::rast(fsca_tif)
  )
  sel_name <- names(which.max(sapply(meta, pix_surf)))
  template <- meta[[sel_name]]
  res_sel  <- terra::res(template)
  message("→ Résolution retenue : ", res_sel[1], " m × ", res_sel[2],
          " m (source ", sel_name, ")")
  
  ## 2. UP projeté dans le CRS template ---------------------------------
  up_vect <- terra::project(up_vect_orig, template)
  
  ## 3. Helper lecture/alignement ---------------------------------------
  read_align <- function(file, method) {
    terra::rast(file) |>
      terra::project(template, method = method) |>
      terra::crop(up_vect) |>
      terra::mask(up_vect)
  }
  
  irg_stack  <- read_align(irg_tif,  "bilinear")
  ndvi_stack <- read_align(ndvi_tif, "bilinear")
  habitat_r  <- read_align(habitat_tif, "near")
  dah_r      <- read_align(dah_tif,  "bilinear")
  fsca_r     <- read_align(fsca_tif, "bilinear")
  
  ## 4. Table statiques --------------------------------------------------
  static_df <- terra::as.data.frame(c(dah_r, fsca_r, habitat_r),
                                    cells=TRUE, na.rm=FALSE) |>
    rlang::set_names(c("cell","dah","fsca","habitat_code"))
  
  ## 5. LUT habitat ------------------------------------------------------
  habitat_lut <- read_csv2(lut_habitat_csv, col_types="ic",
                           locale=locale(encoding="latin1")) |>
    dplyr::mutate(code = as.integer(code))
  
  ## 6. Charge pastoral --------------------------------------------------
  charge_df <- readRDS(charge_rds) |> tibble::as_tibble()
  if (!"cell" %in% names(charge_df)) {
    charge_df <- charge_df |>
      dplyr::mutate(cell = terra::cellFromXY(template, cbind(x,y))) |>
      dplyr::filter(!is.na(cell))
  }
  charge_df <- charge_df |> dplyr::select(cell, Charge)
  
  ## 7. Long tables IRG/NDVI --------------------------------------------
  to_long <- function(stk, var){
    terra::as.data.frame(stk, cells=TRUE, xy=TRUE, na.rm=FALSE) |>
      tidyr::pivot_longer(-c(cell,x,y), names_to="band", values_to=var) |>
      dplyr::mutate(DOY = as.integer(stringr::str_extract(band,"\\d{3,}$"))) |>
      dplyr::filter(DOY %in% DOY_range) |>
      dplyr::select(-band)
  }
  irg_long  <- to_long(irg_stack,  "IRG")
  ndvi_long <- to_long(ndvi_stack, "NDVI")
  
  pheno_df <- dplyr::inner_join(irg_long, ndvi_long,
                                by=c("cell","x","y","DOY"))
  rm(irg_long, ndvi_long); gc()
  
  ## 8. Fusion finale ----------------------------------------------------
  full_df <- pheno_df |>
    dplyr::left_join(static_df,  by="cell") |>
    dplyr::left_join(charge_df,  by="cell") |>
    dplyr::mutate(habitat_code = as.integer(habitat_code)) |>
    dplyr::left_join(habitat_lut, by=c("habitat_code"="code")) |>
    dplyr::rename(habitat_label = label) |>
    dplyr::relocate(habitat_label, .after=habitat_code)
  
  ## 9. Diagnostic -------------------------------------------------------
  pix_area <- prod(res_sel)
  n_pix    <- dplyr::n_distinct(full_df$cell)
  theor_pix<- terra::expanse(up_vect)/pix_area
  message("Pixels uniques (≠ DOY) : ", format(n_pix, big.mark=" "))
  message("Pixels théoriques      : ", round(theor_pix))
  
  ## 10. Write FST -------------------------------------------------------
  fst::write_fst(full_df, output_fst, compress=60)
  message("✓ Fichier écrit : ", output_fst)
  
  invisible(full_df)
}




## NEW 

build_pheno_dataset <- function(alpage, YEAR,
                                irg_tif, ndvi_tif,
                                habitat_tif, dah_tif, fsca_tif,
                                charge_rds,             # .rds chargement
                                lut_habitat_csv,
                                up_shape,               # UP EPSG:2154
                                output_fst,
                                DOY_range = 121:334) {
  
  suppressPackageStartupMessages({
    library(terra);      library(dplyr);   library(tidyr)
    library(readr);      library(stringr); library(fst)
    library(conflicted)
  })
  # toujours les verbes dplyr plutôt que terra
  conflict_prefer("select",   "dplyr", quiet = TRUE)
  conflict_prefer("filter",   "dplyr", quiet = TRUE)
  conflict_prefer("mutate",   "dplyr", quiet = TRUE)
  conflict_prefer("summarise","dplyr", quiet = TRUE)
  
  ## 0. vérif fichiers ---------------------------------------------------
  needed <- c(irg_tif, ndvi_tif, habitat_tif, dah_tif, fsca_tif,
              charge_rds, lut_habitat_csv, up_shape)
  miss <- needed[!file.exists(needed)]
  if (length(miss))
    stop("Fichiers manquants :\n", paste(" •", miss, collapse = "\n"))
  dir.create(dirname(output_fst), recursive = TRUE, showWarnings = FALSE)
  
  ## 1. UP & choix template ---------------------------------------------
  up_vect_orig <- terra::vect(up_shape)                    # EPSG:2154
  cat("Surface UP :", round(terra::expanse(up_vect_orig)/1e6, 2), "km²\n")
  
  pixsurf <- \(r) prod(terra::res(r))
  meta <- list(
    IRG  = terra::rast(irg_tif)[[1]],
    NDVI = terra::rast(ndvi_tif)[[1]],
    HAB  = terra::rast(habitat_tif),
    DAH  = terra::rast(dah_tif),
    FSCA = terra::rast(fsca_tif)
  )
  sel_name <- names(which.max(sapply(meta, pixsurf)))      # maille la + grosse
  template <- meta[[sel_name]]
  cat("→ Template :", sel_name,
      "–", paste(terra::res(template), collapse = "×"), "m\n")
  
  up_vect <- terra::project(up_vect_orig, template)
  
  read_align <- function(f, m)
    terra::rast(f) |>
    terra::project(template, method = m) |>
    terra::crop(up_vect) |> terra::mask(up_vect)
  
  irg_stack  <- read_align(irg_tif,  "bilinear")
  ndvi_stack <- read_align(ndvi_tif, "bilinear")
  habitat_r  <- read_align(habitat_tif, "near")
  dah_r      <- read_align(dah_tif,   "bilinear")
  fsca_r     <- read_align(fsca_tif,  "bilinear")
  
  ## 2. Statique ---------------------------------------------------------
  static_df <- terra::as.data.frame(c(dah_r, fsca_r, habitat_r),
                                    cells = TRUE, na.rm = FALSE) |>
    rlang::set_names(c("cell","dah","fsca","habitat_code"))
  
  habitat_lut <- read_csv2(lut_habitat_csv, col_types = "ic",
                           locale = locale(encoding = "latin1")) |>
    mutate(code = as.integer(code))
  
  ## 3. Table chargement -------------------------------------------------
  charge_df <- readRDS(charge_rds) |> as_tibble()
  
  # colonne Charge (sensibilité casse)
  if (!"Charge" %in% names(charge_df)) {
    cand <- grep("charge", names(charge_df), ignore.case = TRUE, value = TRUE)
    if (length(cand))
      charge_df <- rename(charge_df, Charge = !!cand[1])
    else
      stop("Colonne 'Charge' introuvable dans le .rds.")
  }
  
  # x,y : les créer si on n’a qu’un cell (grille fine)
  if (!all(c("x","y") %in% names(charge_df))) {
    if (!"cell" %in% names(charge_df))
      stop("Le .rds doit contenir x/y ou cell (grille fine).")
    xy10 <- terra::xyFromCell(meta$IRG, charge_df$cell)
    charge_df <- mutate(charge_df, x = xy10[,1], y = xy10[,2])
  }
  
  # reprojection points -> CRS template
  pts <- terra::vect(cbind(charge_df$x, charge_df$y), crs = "EPSG:2154") |>
    terra::project(template)
  charge_df$x <- terra::crds(pts)[,1]
  charge_df$y <- terra::crds(pts)[,2]
  
  # nouveau cell sur la grille du template
  charge_df <- mutate(charge_df,
                      cell = terra::cellFromXY(template, cbind(x,y))) |>
    filter(!is.na(cell))
  
  # agrégation
  if ("DOY" %in% names(charge_df)) {
    charge_df <- group_by(charge_df, cell, DOY) |>
      summarise(Charge = mean(Charge, na.rm = TRUE), .groups="drop")
    by_cols <- c("cell","DOY")
  } else {
    charge_df <- group_by(charge_df, cell) |>
      summarise(Charge = mean(Charge, na.rm = TRUE), .groups="drop")
    by_cols <- "cell"
  }
  
  ## 4. IRG / NDVI -> long ----------------------------------------------
  to_long <- function(stk, var)
    terra::as.data.frame(stk, cells = TRUE, xy = TRUE, na.rm = FALSE) |>
    pivot_longer(-c(cell,x,y), names_to = "band", values_to = var) |>
    mutate(DOY = as.integer(str_extract(band, "\\d{3,}$"))) |>
    filter(DOY %in% DOY_range) |>
    select(-band)
  
  irg_long  <- to_long(irg_stack,  "IRG")
  ndvi_long <- to_long(ndvi_stack, "NDVI")
  pheno_df  <- inner_join(irg_long, ndvi_long,
                          by = c("cell","x","y","DOY"))
  rm(irg_long, ndvi_long); gc()
  
  ## 5. Fusion finale ----------------------------------------------------
  full_df <- pheno_df |>
    left_join(static_df, by = "cell") |>
    left_join(charge_df, by = by_cols) |>
    mutate(habitat_code = as.integer(habitat_code)) |>
    left_join(habitat_lut, by = c("habitat_code" = "code")) |>
    rename(habitat_label = label) |>
    relocate(habitat_label, .after = habitat_code)
  
  ## 6. Diagnostics & export --------------------------------------------
  n_pix <- n_distinct(full_df$cell)
  theor <- terra::expanse(up_vect) / pixsurf(template)
  cat("Pixels uniques :", format(n_pix, big.mark=" "),
      "(théorique ~", round(theor), ")\n")
  
  write_fst(full_df, output_fst, compress = 60)
  cat("✓ Écrit :", output_fst, "\n")
  
  invisible(full_df)
}










## NEW 2

#' Construire un dataset phénologique avec projection uniforme EPSG:2154 et export en .rds
#'
#' @param alpage   Identifiant de l'alpage (non utilisé directement ici)
#' @param YEAR     Année d'étude (non utilisé directement ici)
#' @param irg_tif  Chemin vers le fichier .tif IRG
#' @param ndvi_tif Chemin vers le fichier .tif NDVI
#' @param habitat_tif Chemin vers le .tif habitat
#' @param dah_tif  Chemin vers le .tif DAH
#' @param fsca_tif Chemin vers le .tif FSCA
#' @param charge_rds Chemin vers le .rds de charges (doit contenir x/y ou cell et éventuellement DOY)
#' @param lut_habitat_csv Chemin vers le .csv de correspondance habitat (code → label)
#' @param up_shape Chemin vers le shapefile de l'unité paysagère (EPSG:2154)
#' @param output_rds Chemin de sortie pour l'enregistrement du data.frame complet en .rds
#' @param DOY_range Vecteur de jours de l'année à considérer (par défaut 121:334)
build_pheno_dataset <- function(alpage, YEAR,
                                irg_tif, ndvi_tif,
                                habitat_tif, dah_tif, fsca_tif,
                                charge_rds,
                                lut_habitat_csv,
                                up_shape,
                                output_rds,
                                DOY_range = 60:365) {
  suppressPackageStartupMessages({
    library(terra)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(purrr)     #  <<— ajouté
    library(conflicted)
  })
  conflict_prefer("select",   "dplyr", quiet = TRUE)
  conflict_prefer("filter",   "dplyr", quiet = TRUE)
  conflict_prefer("mutate",   "dplyr", quiet = TRUE)
  conflict_prefer("summarise","dplyr", quiet = TRUE)
  
  ## 0. Vérification des fichiers -------------------------------------------------
  needed <- c(irg_tif, ndvi_tif, habitat_tif, dah_tif, fsca_tif,
              charge_rds, lut_habitat_csv, up_shape)
  miss <- needed[!file.exists(needed)]
  if (length(miss))
    stop("Fichiers manquants :\n", paste(" •", miss, collapse = "\n"))
  dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
  
  ## 1. Chargement UP et création du template EPSG:2154 -----------------------------
  up_vect <- terra::vect(up_shape)            # Doit être en EPSG:2154
  cat("Surface UP :", round(terra::expanse(up_vect)/1e6,2), "km²\n")
  
  # Chargement des rasters originaux
  meta_orig <- list(
    IRG  = terra::rast(irg_tif)[[1]],
    NDVI = terra::rast(ndvi_tif)[[1]],
    HAB  = terra::rast(habitat_tif),
    DAH  = terra::rast(dah_tif),
    FSCA = terra::rast(fsca_tif)
  )
  # Choix de la plus grande surface pixel (unités d'origine)
  pix_area <- function(r) prod(terra::res(r))
  sel_name <- names(which.max(sapply(meta_orig, pix_area)))
  base_res <- terra::res(meta_orig[[sel_name]])
  cat("→ Résolution choisie (", sel_name,") : ",
      paste(base_res, collapse="×"), " (unités du CRS)\n", sep="")
  
  # Template définissant EPSG:2154, étendue UP, résolution uniforme
  template <- terra::rast(ext = terra::ext(up_vect),
                          resolution = base_res,
                          crs = crs(up_vect))
  cat("Template CRS :", crs(template), "\n")
  
  # Fonction utilitaire pour lire, reprojeter, découper et masquer
  read_align <- function(file, method) {
    terra::rast(file) %>%
      terra::project(template, method = method) %>%
      terra::crop(up_vect) %>%
      terra::mask(up_vect)
  }
  
  # Reprojection de toutes les couches dans le même CRS et résolution
  irg_stack  <- read_align(irg_tif,  "bilinear")
  ndvi_stack <- read_align(ndvi_tif, "bilinear")
  habitat_r  <- read_align(habitat_tif, "near")
  dah_r      <- read_align(dah_tif,   "bilinear")
  fsca_r     <- read_align(fsca_tif,  "bilinear")
  
  ## 2. Extraction des variables statiques ----------------------------------------
  static_df <- terra::as.data.frame(c(dah_r, fsca_r, habitat_r),
                                    cells = TRUE, na.rm = FALSE) |>
    setNames(c("cell", "dah", "fsca", "habitat_code"))
  
  habitat_lut <- read_csv2(lut_habitat_csv, col_types = "ic",
                           locale = locale(encoding = "latin1")) %>%
    mutate(code = as.integer(code))
  
  ## 3. Chargement et préparation des données de charge ---------------------------
  charge_df <- readRDS(charge_rds) %>% as_tibble()
  if (!"Charge" %in% names(charge_df)) {
    cand <- grep("charge", names(charge_df), ignore.case = TRUE, value = TRUE)
    if (length(cand)) charge_df <- rename(charge_df, Charge = !!sym(cand[1]))
    else stop("Colonne 'Charge' introuvable dans le .rds")
  }
  if (!all(c("x","y") %in% names(charge_df))) {
    if (!"cell" %in% names(charge_df))
      stop("Le .rds doit contenir x/y ou cell (grille fine)")
    xy0 <- terra::xyFromCell(meta_orig$IRG, charge_df$cell)
    charge_df <- mutate(charge_df, x = xy0[,1], y = xy0[,2])
  }
  # Reprojection des points vers EPSG:2154 puis recalcule du cell
  pts <- terra::vect(cbind(charge_df$x, charge_df$y), crs = crs(meta_orig$IRG)) %>%
    terra::project(template)
  charge_df <- mutate(charge_df,
                      x    = crds(pts)[,1],
                      y    = crds(pts)[,2],
                      cell = terra::cellFromXY(template, cbind(x,y))) %>%
    filter(!is.na(cell))
  if ("DOY" %in% names(charge_df)) {
    charge_df <- group_by(charge_df, cell, DOY) %>%
      summarise(Charge = mean(Charge, na.rm = TRUE), .groups = "drop")
    by_cols <- c("cell", "DOY")
  } else {
    charge_df <- group_by(charge_df, cell) %>%
      summarise(Charge = mean(Charge, na.rm = TRUE), .groups = "drop")
    by_cols <- "cell"
  }
  
  ## 4. Format long pour IRG & NDVI ---------------------------------------------
  to_long <- function(stk, var) {
    terra::as.data.frame(stk, cells = TRUE, xy = TRUE, na.rm = FALSE) %>%
      pivot_longer(-c(cell, x, y), names_to = "band", values_to = var) %>%
      mutate(DOY = as.integer(str_extract(band, "\\d{3,}$"))) %>%
      filter(DOY %in% DOY_range) %>%
      select(-band)
  }
  irg_long  <- to_long(irg_stack,  "IRG")
  ndvi_long <- to_long(ndvi_stack, "NDVI")
  pheno_df  <- inner_join(irg_long, ndvi_long, by = c("cell","x","y","DOY"))
  rm(irg_long, ndvi_long); gc()
  
  ## 5. Fusion finale -----------------------------------------------------------
  full_df <- pheno_df %>%
    left_join(static_df, by = "cell") %>%
    left_join(charge_df, by = by_cols) %>%
    mutate(habitat_code = as.integer(habitat_code)) %>%
    left_join(habitat_lut, by = c("habitat_code" = "code")) %>%
    rename(habitat_label = label) %>%
    relocate(habitat_label, .after = habitat_code)
  
  ## 6. Export du résultat en .rds ------------------------------------------------
  n_pix <- n_distinct(full_df$cell)
  theor <- terra::expanse(up_vect) / pix_area(template)
  cat("Pixels uniques :", format(n_pix, big.mark = " "),
      "(théorique ~", round(theor), ")\n")
  write_rds(full_df, output_rds, compress = "gz")
  cat("✓ Écrit :", output_rds, "\n")
  
  invisible(full_df)
}














## NEW 3 : 


###############################################################################
# build_pheno_dataset  – prépare le dataset « long » pixel × DOY
#                       + 17 bandes PPI_extra (MINV, MAXV, …, AsymSlope)
###############################################################################
build_pheno_dataset <- function(alpage, YEAR,
                                irg_tif, ndvi_tif,
                                habitat_tif, dah_tif, fsca_tif,
                                extra_tif,                    # ← NEW
                                charge_rds,
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
  
  ## 0. Vérification des fichiers -------------------------------------------
  needed <- c(irg_tif, ndvi_tif, habitat_tif, dah_tif, fsca_tif,
              extra_tif,                          # ← NEW
              charge_rds, lut_habitat_csv, up_shape)
  miss <- needed[!file.exists(needed)]
  if (length(miss)) stop("Fichiers manquants :\n",
                         paste(" •", miss, collapse = "\n"))
  dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
  
  ## 1. Gabarit EPSG:2154 ----------------------------------------------------
  up_vect <- vect(up_shape)               # EPSG 2154
  meta_orig <- list(
    IRG  = rast(irg_tif)[[1]],
    NDVI = rast(ndvi_tif)[[1]],
    HAB  = rast(habitat_tif),
    DAH  = rast(dah_tif),
    FSCA = rast(fsca_tif),
    EXTRA= rast(extra_tif)[[1]]           # ← pour résolution de départ
  )
  base_res <- res(meta_orig[[which.max(sapply(meta_orig, \(r) prod(res(r))))]])
  template <- rast(ext=ext(up_vect), resolution=base_res, crs=crs(up_vect))
  
  read_align <- function(file, method)
    rast(file) |> project(template, method=method) |> crop(up_vect) |> mask(up_vect)
  
  irg_stack   <- read_align(irg_tif,   "bilinear")
  ndvi_stack  <- read_align(ndvi_tif,  "bilinear")
  habitat_r   <- read_align(habitat_tif,"near")
  dah_r       <- read_align(dah_tif,   "bilinear")
  fsca_r      <- read_align(fsca_tif,  "bilinear")
  extra_stack <- read_align(extra_tif, "bilinear")      # ← NEW
  
  ## 2. Variables statiques (dah, fsca, habitat + PPI_extra) -----------------
  static_df <- as.data.frame(c(dah_r, fsca_r, habitat_r),
                             cells=TRUE, na.rm=FALSE) |>
    setNames(c("cell","dah","fsca","habitat_code"))
  
  extra_df  <- as.data.frame(extra_stack, cells=TRUE, na.rm=FALSE)
  static_df <- left_join(static_df, extra_df, by="cell")   # ← fusion extra
  
  habitat_lut <- read_csv2(lut_habitat_csv, col_types="ic",
                           locale=locale(encoding="latin1")) |>
    mutate(code=as.integer(code))
  
  ## 3. Chargements ----------------------------------------------------------
  charge_df <- readRDS(charge_rds) |> as_tibble()
  if (!"Charge" %in% names(charge_df))
    charge_df <- rename(charge_df, Charge = matches("charge", ignore.case=TRUE))
  if (!all(c("x","y") %in% names(charge_df))) {
    xy0 <- xyFromCell(meta_orig$IRG, charge_df$cell)
    charge_df <- mutate(charge_df, x=xy0[,1], y=xy0[,2])
  }
  pts <- vect(cbind(charge_df$x, charge_df$y), crs=crs(meta_orig$IRG)) |>
    project(template)
  charge_df <- mutate(charge_df,
                      x=crds(pts)[,1], y=crds(pts)[,2],
                      cell=cellFromXY(template, cbind(x,y))) |>
    filter(!is.na(cell))
  if ("DOY" %in% names(charge_df)) {
    charge_df <- group_by(charge_df, cell, DOY) |>
      summarise(Charge=mean(Charge), .groups="drop")
    by_cols <- c("cell","DOY")
  } else {
    charge_df <- group_by(charge_df, cell) |>
      summarise(Charge=mean(Charge), .groups="drop")
    by_cols <- "cell"
  }
  
  ## 4. Format long IRG & NDVI ----------------------------------------------
  to_long <- function(stk, var)
    as.data.frame(stk, cells=TRUE, xy=TRUE, na.rm=FALSE) |>
    pivot_longer(-c(cell,x,y), names_to="band", values_to=var) |>
    mutate(DOY = as.integer(str_extract(band, "\\d{3,}$"))) |>
    filter(DOY %in% DOY_range) |>
    select(-band)
  
  pheno_df <- inner_join(to_long(irg_stack,"IRG"),
                         to_long(ndvi_stack,"NDVI"),
                         by=c("cell","x","y","DOY"))
  
  ## 5. Fusion finale --------------------------------------------------------
  full_df <- pheno_df |>
    left_join(static_df, by="cell") |>
    left_join(charge_df, by=by_cols) |>
    mutate(habitat_code = as.integer(habitat_code)) |>
    left_join(habitat_lut, by=c("habitat_code"="code")) |>
    rename(habitat_label = label) |>
    relocate(habitat_label, .after=habitat_code)
  
  ## 6. Export ---------------------------------------------------------------
  write_rds(full_df, output_rds, compress="gz")
  cat("✓ Dataset écrit :", output_rds, "\n")
  
  invisible(full_df)
}
























diag_dataset <- function(
    fst_path,
    up_shape   = NULL,          # shapefile UP EPSG:2154 ou NULL
    sample_n   = 5e5,
    bins       = 60,
    output_dir = dirname(fst_path),
    palette    = c(IRG  = "#1b9e77",
                   NDVI = "#d95f02",
                   fsca = "#7570b3",
                   dah  = "#e7298a",
                   Charge = "tomato"),
    view       = FALSE,
    quiet      = FALSE) {
  
  suppressPackageStartupMessages({
    library(fst);  library(dplyr);  library(tidyr)
    library(ggplot2);  library(sf);   library(scales);  library(grid)
  })
  
  stopifnot(file.exists(fst_path))
  names(palette) <- tolower(names(palette))
  pcol <- function(v) palette[[tolower(v)]]
  
  ## ── méta ------------------------------------------------------------
  meta   <- fst::metadata_fst(fst_path)
  cnames <- meta$columnNames
  nrows  <- meta$nrOfRows %||% meta$nrows
  charge_col <- cnames[ trimws(tolower(cnames)) == "charge" ]
  if (!length(charge_col)) charge_col <- ""
  
  ## ── emprise ---------------------------------------------------------
  have_up <- !is.null(up_shape) && file.exists(up_shape)
  if (have_up) {
    up_sf <- st_read(up_shape, quiet = TRUE)
    bbox  <- st_bbox(up_sf) + c(-50,-50,50,50)
  } else {
    xy   <- read_fst(fst_path, columns = c("x","y"))
    bbox <- st_bbox(c(xmin=min(xy$x), ymin=min(xy$y),
                      xmax=max(xy$x), ymax=max(xy$y)),
                    crs=st_crs(2154)) + c(-50,-50,50,50)
  }
  poly_bbox <- st_as_sfc(bbox)
  
  ## ── échantillon -----------------------------------------------------
  vars <- c("IRG","NDVI","fsca","dah", if (charge_col!="") charge_col)
  rows <- if (nrows<=sample_n) seq_len(nrows) else sort(sample(nrows,sample_n))
  df   <- read_fst(fst_path, columns=vars, from=min(rows), to=max(rows))
  df   <- df[rows - min(rows) + 1, ]
  
  ## ── liste de graphiques --------------------------------------------
  grobs <- list()
  
  # 1) carte
  g_map <- ggplot() +
    { if (have_up) geom_sf(data=up_sf, fill="grey85",
                           colour="steelblue", linewidth=.5) } +
    geom_sf(data = poly_bbox, fill = NA,
            colour = "red3", linewidth = 1) +
    coord_sf(expand = FALSE) +
    labs(title="Emprise du jeu de données",
         subtitle = basename(fst_path)) +
    theme_void()
  grobs[[1]] <- g_map
  
  # 2-5) histogrammes classiques
  for (v in c("IRG","NDVI","fsca","dah"))
    grobs[[length(grobs)+1]] <-
    ggplot(df, aes(.data[[v]])) +
    geom_histogram(bins=bins, fill=pcol(v), colour="grey25") +
    labs(title=paste("Distribution de",v),
         subtitle=paste0("Échantillon ",
                         format(nrow(df), big.mark=" ")," lignes"),
         x=v, y="Fréquence") +
    theme_minimal(10)
  
  # 6-7) Charge
  if (charge_col!="") {
    df_bar <- df |> mutate(etat = ifelse(!is.na(.data[[charge_col]]) &
                                           .data[[charge_col]]>=1,
                                         "Chargé (≥1)","Non chargé (<1)"))
    
    grobs[[length(grobs)+1]] <-
      ggplot(df_bar, aes(etat)) +
      geom_bar(fill=pcol("charge")) +
      scale_y_continuous(labels=comma) +
      labs(title="Pixels chargés / non chargés",
           x=NULL, y="Pixels") +
      theme_minimal(10)
    
    df_pos <- filter(df_bar, etat=="Chargé (≥1)")
    
    grobs[[length(grobs)+1]] <-
      ggplot(df_pos, aes(.data[[charge_col]])) +
      geom_histogram(bins=bins, fill=pcol("charge"),
                     colour="grey25") +
      scale_x_log10(labels=comma) +
      labs(title="Distribution de Charge (>0)",
           subtitle=paste("Échantillon",
                          format(nrow(df_pos),big.mark=" "),"lignes"),
           x="Charge (log10)", y="Fréquence") +
      theme_minimal(10)
  }
  
  ## ── export PDF ------------------------------------------------------
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  pdf_path <- file.path(output_dir,
                        paste0(tools::file_path_sans_ext(basename(fst_path)), "_plots.pdf"))
  
  pdf(pdf_path, width=7, height=5)
  grid.draw(grobs[[1]])                  # première page, pas de newpage
  if (length(grobs) > 1)
    for (g in grobs[-1]) { grid.newpage(); grid.draw(g) }
  dev.off()
  
  ## ── aperçu RStudio --------------------------------------------------
  if (view)
    for (g in grobs) { dev.new(width=7,height=5); grid.draw(g) }
  
  if (!quiet) message("✓ PDF écrit : ", pdf_path)
  invisible(pdf_path)
}





## NEW 1 : 

#' Diagnostic rapide d'un dataset phéno en .rds
#'
#' @param rds_path   Chemin vers le fichier .rds généré par build_pheno_dataset
#' @param up_shape   Chemin vers le shapefile UP (EPSG:2154) ou NULL pour autodétecter l'étendue
#' @param sample_n   Nombre de pixels à échantillonner (défaut 5e5)
#' @param bins       Nombre de classes pour les histogrammes (défaut 60)
#' @param output_dir Répertoire de sortie pour le PDF (défaut : même dossier que rds_path)
#' @param palette    Palette couleurs nommée pour les variables
#' @param view       Ouvrir les graphiques en extra-window (TRUE/FALSE)
#' @param quiet      Supprimer les messages d'avancement (TRUE/FALSE)
diag_dataset <- function(
    rds_path,
    up_shape   = NULL,
    sample_n   = 5e5,
    bins       = 60,
    output_dir = dirname(rds_path),
    palette    = c(IRG  = "#1b9e77",
                   NDVI = "#d95f02",
                   fsca = "#7570b3",
                   dah  = "#e7298a",
                   Charge = "tomato"),
    view       = FALSE,
    quiet      = FALSE) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr)
    library(ggplot2);  library(sf);   library(scales);  library(grid)
  })
  stopifnot(file.exists(rds_path))
  names(palette) <- tolower(names(palette))
  pcol <- function(v) palette[[tolower(v)]]
  
  ## 0. Chargement du data.frame depuis le .rds ------------------------------
  df_full <- readRDS(rds_path)
  if (!is.data.frame(df_full))
    stop("Le fichier rds ne contient pas un data.frame.")
  cnames <- names(df_full)
  nrows  <- nrow(df_full)
  # détecte colonne Charge si présente
  charge_col <- cnames[tolower(cnames)=="charge"] %||% character(0)
  
  ## 1. Définition de l'étendue géographique ----------------------------------
  have_up <- !is.null(up_shape) && file.exists(up_shape)
  if (have_up) {
    up_sf <- st_read(up_shape, quiet = TRUE)
    bbox  <- st_bbox(up_sf) + c(-50,-50,50,50)
  } else {
    if (!all(c("x","y") %in% cnames))
      stop("Sans up_shape, le data.frame doit contenir x et y.")
    coords <- df_full %>% select(x,y)
    bbox <- st_bbox(c(xmin=min(coords$x), ymin=min(coords$y),
                      xmax=max(coords$x), ymax=max(coords$y)),
                    crs=st_crs(2154)) + c(-50,-50,50,50)
  }
  poly_bbox <- st_as_sfc(bbox)
  
  ## 2. Échantillonnage ------------------------------------------------------
  vars <- c("IRG","NDVI","fsca","dah", charge_col)
  vars <- vars[vars %in% cnames]
  # échantillonnage aléatoire sur les indices
  rows <- if (nrows <= sample_n) seq_len(nrows)
  else sample(nrows, sample_n)
  df <- df_full[rows, vars]
  
  ## 3. Création des graphes -------------------------------------------------
  grobs <- list()
  # carte d'emprise
  g_map <- ggplot() +
    { if (have_up) geom_sf(data=up_sf, fill="grey85", colour="steelblue", linewidth=.5) } +
    geom_sf(data = poly_bbox, fill = NA, colour = "red3", linewidth = 1) +
    coord_sf(expand = FALSE) +
    labs(title="Emprise du jeu de données", subtitle = basename(rds_path)) +
    theme_void()
  grobs[[1]] <- g_map
  
  # histogrammes IRG, NDVI, fsca, dah
  for (v in c("IRG","NDVI","fsca","dah")[c("IRG","NDVI","fsca","dah") %in% names(df)]) {
    grobs[[length(grobs)+1]] <-
      ggplot(df, aes(.data[[v]])) +
      geom_histogram(bins = bins, fill = pcol(v), colour = "grey25") +
      labs(title = paste("Distribution de", v),
           subtitle = paste0("Échantillon ", format(nrow(df), big.mark=" "), " lignes"),
           x=v, y="Fréquence") +
      theme_minimal(10)
  }
  
  # diagnostics Charge si colonne présente
  if (length(charge_col)==1 && charge_col != "") {
    df_bar <- df %>% mutate(etat = ifelse(!is.na(.data[[charge_col]]) & .data[[charge_col]]>=1,
                                          "Chargé (≥1)", "Non chargé (<1)"))
    grobs[[length(grobs)+1]] <-
      ggplot(df_bar, aes(etat)) + geom_bar(fill=pcol("charge")) +
      scale_y_continuous(labels = comma) +
      labs(title="Pixels chargés / non chargés", x=NULL, y="Pixels") +
      theme_minimal(10)
    
    df_pos <- df_bar %>% filter(etat=="Chargé (≥1)")
    grobs[[length(grobs)+1]] <-
      ggplot(df_pos, aes(.data[[charge_col]])) +
      geom_histogram(bins=bins, fill=pcol("charge"), colour="grey25") +
      scale_x_log10(labels = comma) +
      labs(title="Distribution de Charge (>0)",
           subtitle = paste("Échantillon", format(nrow(df_pos), big.mark=" "), "lignes"),
           x="Charge (log10)", y="Fréquence") +
      theme_minimal(10)
  }
  
  ## 4. Export PDF -----------------------------------------------------------
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  pdf_path <- file.path(output_dir,
                        paste0(tools::file_path_sans_ext(basename(rds_path)), "_plots.pdf"))
  pdf(pdf_path, width = 7, height = 5)
  grid.draw(grobs[[1]])
  if (length(grobs) > 1)
    for (g in grobs[-1]) { grid.newpage(); grid.draw(g) }
  dev.off()
  
  if (!quiet) message("✓ PDF écrit : ", pdf_path)
  if (view) for (g in grobs) { dev.new(width=7, height=5); grid.draw(g) }
  invisible(pdf_path)
}










#' Diagnostic rapide d'un dataset phéno en .rds avec stats descriptives
#'
#' @param rds_path   Chemin vers le fichier .rds généré par build_pheno_dataset
#' @param up_shape   Chemin vers le shapefile UP (EPSG:2154) ou NULL pour autodétecter l'étendue
#' @param sample_n   Nombre de pixels à échantillonner (défaut 5e5)
#' @param bins       Nombre de classes pour les histogrammes (défaut 60)
#' @param output_dir Répertoire de sortie pour le PDF (défaut : même dossier que rds_path)
#' @param palette    Palette couleurs nommée pour les variables
#' @param view       Ouvrir les graphiques en extra-window (TRUE/FALSE)
#' @param quiet      Supprimer les messages d'avancement (TRUE/FALSE)
diag_dataset <- function(
    rds_path,
    up_shape   = NULL,
    sample_n   = 5e5,
    bins       = 60,
    output_dir = dirname(rds_path),
    palette    = c(IRG  = "#1b9e77",
                   NDVI = "#d95f02",
                   fsca = "#7570b3",
                   dah  = "#e7298a",
                   Charge = "tomato"),
    view       = FALSE,
    quiet      = FALSE) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr)
    library(ggplot2);  library(sf);   library(scales);  library(grid)
  })
  stopifnot(file.exists(rds_path))
  names(palette) <- tolower(names(palette))
  pcol <- function(v) palette[[tolower(v)]]
  
  ## 0. Chargement du data.frame depuis le .rds ------------------------------
  df_full <- readRDS(rds_path)
  if (!is.data.frame(df_full))
    stop("Le fichier rds ne contient pas un data.frame.")
  cnames <- names(df_full)
  nrows  <- nrow(df_full)
  
  ## 0.1 Statistiques descriptives -------------------------------------------
  if (!quiet) {
    message("Colonnes : ", paste(cnames, collapse=", "))
    message("Nombre de pixels (lignes) : ", nrows)
    message("Nombre de variables (colonnes) : ", length(cnames))
    # Stats pour variables numériques
    num_vars <- cnames[sapply(df_full, is.numeric)]
    for (v in num_vars) {
      vals <- df_full[[v]]
      message(sprintf("  %s : mean=%.2f, median=%.2f, sd=%.2f, NA=%d",
                      v,
                      mean(vals, na.rm=TRUE),
                      median(vals, na.rm=TRUE),
                      sd(vals, na.rm=TRUE),
                      sum(is.na(vals))))
    }
  }
  
  # détecte colonne Charge si présente
  charge_col <- cnames[tolower(cnames)=="charge"] %||% character(0)
  
  ## 1. Définition de l'étendue géographique ----------------------------------
  have_up <- !is.null(up_shape) && file.exists(up_shape)
  if (have_up) {
    up_sf <- st_read(up_shape, quiet = TRUE)
    bbox  <- st_bbox(up_sf) + c(-50,-50,50,50)
  } else {
    if (!all(c("x","y") %in% cnames))
      stop("Sans up_shape, le data.frame doit contenir x et y.")
    coords <- df_full %>% select(x,y)
    bbox <- st_bbox(c(xmin=min(coords$x), ymin=min(coords$y),
                      xmax=max(coords$x), ymax=max(coords$y)),
                    crs=st_crs(2154)) + c(-50,-50,50,50)
  }
  poly_bbox <- st_as_sfc(bbox)
  
  ## 2. Échantillonnage ------------------------------------------------------
  vars <- c("IRG","NDVI","fsca","dah", charge_col)
  vars <- vars[vars %in% cnames]
  rows <- if (nrows <= sample_n) seq_len(nrows) else sample(nrows, sample_n)
  df <- df_full[rows, vars]
  
  ## 3. Création des graphes -------------------------------------------------
  grobs <- list()
  # carte d'emprise
  g_map <- ggplot() +
    { if (have_up) geom_sf(data=up_sf, fill="grey85", colour="steelblue", linewidth=.5) } +
    geom_sf(data = poly_bbox, fill = NA, colour = "red3", linewidth = 1) +
    coord_sf(expand = FALSE) +
    labs(title="Emprise du jeu de données", subtitle = basename(rds_path)) +
    theme_void()
  grobs[[1]] <- g_map
  
  # histogrammes IRG, NDVI, fsca, dah
  for (v in c("IRG","NDVI","fsca","dah")[c("IRG","NDVI","fsca","dah") %in% names(df)]) {
    grobs[[length(grobs)+1]] <-
      ggplot(df, aes(.data[[v]])) +
      geom_histogram(bins = bins, fill = pcol(v), colour = "grey25") +
      labs(title = paste("Distribution de", v),
           subtitle = paste0("Échantillon ", format(nrow(df), big.mark=" "), " lignes"),
           x=v, y="Fréquence") +
      theme_minimal(10)
  }
  
  # diagnostics Charge si colonne présente
  if (length(charge_col)==1 && charge_col != "") {
    df_bar <- df %>% mutate(etat = ifelse(!is.na(.data[[charge_col]]) & .data[[charge_col]]>=1,
                                          "Chargé (≥1)", "Non chargé (<1)"))
    grobs[[length(grobs)+1]] <-
      ggplot(df_bar, aes(etat)) + geom_bar(fill=pcol("charge")) +
      scale_y_continuous(labels = comma) +
      labs(title="Pixels chargés / non chargés", x=NULL, y="Pixels") +
      theme_minimal(10)
    
    df_pos <- df_bar %>% filter(etat=="Chargé (≥1)")
    grobs[[length(grobs)+1]] <-
      ggplot(df_pos, aes(.data[[charge_col]])) +
      geom_histogram(bins=bins, fill=pcol("charge"), colour="grey25") +
      scale_x_log10(labels = comma) +
      labs(title="Distribution de Charge (>0)",
           subtitle = paste("Échantillon", format(nrow(df_pos), big.mark=" "), "lignes"),
           x="Charge (log10)", y="Fréquence") +
      theme_minimal(10)
  }
  
  ## 4. Export PDF -----------------------------------------------------------
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  pdf_path <- file.path(output_dir,
                        paste0(tools::file_path_sans_ext(basename(rds_path)), "_plots.pdf"))
  pdf(pdf_path, width = 7, height = 5)
  grid.draw(grobs[[1]])
  if (length(grobs) > 1)
    for (g in grobs[-1]) { grid.newpage(); grid.draw(g) }
  dev.off()
  
  if (!quiet) message("✓ PDF écrit : ", pdf_path)
  if (view) for (g in grobs) { dev.new(width=7, height=5); grid.draw(g) }
  invisible(pdf_path)
}













## NEW :

###############################################################################
# diag_dataset  – diagnostic PDF (IRG, NDVI, DAH, FSCA + PPI_extra + Charge)
###############################################################################
diag_dataset <- function(
    rds_path,
    up_shape   = NULL,
    sample_n   = 5e5,
    bins       = 60,
    output_dir = dirname(rds_path),
    palette    = c(IRG  = "#1b9e77", NDVI = "#d95f02",
                   fsca = "#7570b3", dah  = "#e7298a",
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
  
  stopifnot(file.exists(rds_path))
  names(palette) <- tolower(names(palette))
  pcol <- function(v) palette[[tolower(v)]]
  
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
  vars <- c("IRG","NDVI","fsca","dah",               # de base
            "MINV","MAXV","AMPL",                    # PPI_extra
            "LENGTH","ONSET10","MaxSlope",
            "GreenUpDur","GreenDownDur","AsymSlope",
            "LSLOPE","RSLOPE",
            charge_col)
  vars <- intersect(vars, cnames)
  rows <- if (nrows <= sample_n) seq_len(nrows) else sample(nrows, sample_n)
  df   <- df_full[rows, vars]
  
  ## 3. Graphiques -----------------------------------------------------------
  grobs <- list()
  # (a) carte d'emprise
  grobs[[1]] <- ggplot() +
    { if (have_up) geom_sf(data = up_sf, fill = "grey85",
                           colour = "steelblue", linewidth = .5) } +
    geom_sf(data = poly_bbox, fill = NA, colour = "red3", linewidth = 1) +
    coord_sf(expand = FALSE) +
    labs(title = "Emprise du jeu de données",
         subtitle = basename(rds_path)) +
    theme_void()
  
  # (b) histogrammes pour toutes les variables sauf Charge
  for (v in setdiff(vars, charge_col)) {
    if (!is.numeric(df[[v]])) next
    grobs[[length(grobs)+1]] <-
      ggplot(df, aes(.data[[v]])) +
      geom_histogram(bins = bins,
                     fill   = pcol(v) %||% "grey60",
                     colour = "grey25") +
      labs(title     = paste("Distribution de", v),
           subtitle  = paste("Échantillon", format(nrow(df), big.mark = " ")),
           x = v, y  = "Fréquence") +
      theme_minimal(9)
  }
  
  ## --- bloc Charge : barre + histo log10 (pixels Charge >= 1) ----------------
  if (length(charge_col) == 1) {
    vC   <- charge_col
    
    # barre Chargé / Non chargé (seuil 1 UA)
    df_bar <- df %>%
      mutate(etat = ifelse(!is.na(.data[[vC]]) & .data[[vC]] >= 1,
                           "Chargé (≥1)", "Non chargé (<1)"))
    
    grobs[[length(grobs)+1]] <-
      ggplot(df_bar, aes(etat)) +
      geom_bar(fill = pcol("charge")) +
      labs(title="Pixels chargés / non chargés",
           x=NULL, y="Pixels") +
      theme_minimal(9)
    
    # histogramme uniquement sur les pixels Charge >= 1
    df_pos <- filter(df_bar, etat == "Chargé (≥1)")
    
    grobs[[length(grobs)+1]] <-
      ggplot(df_pos, aes(.data[[vC]])) +
      geom_histogram(bins = bins,
                     fill = pcol("charge"), colour = "grey25") +
      scale_x_log10(labels = comma) +
      labs(title     = "Distribution de Charge (≥1 UA)",
           subtitle  = paste("Échantillon", format(nrow(df_pos), big.mark=" ")),
           x = "Charge (log10)", y = "Fréquence") +
      theme_minimal(9)
  }
  
  ## 4. Export PDF -----------------------------------------------------------
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(output_dir,
                        paste0(tools::file_path_sans_ext(basename(rds_path)),
                               "_plots.pdf"))
  pdf(pdf_path, 7, 5)
  for (i in seq_along(grobs)) {
    if (i > 1) grid.newpage()
    grid::grid.draw(grobs[[i]])
  }
  dev.off()
  if (!quiet) message("✓ PDF écrit :", pdf_path)
  
  if (view) for (g in grobs) { dev.new(7,5); grid::grid.draw(g) }
  invisible(pdf_path)
}








## create_zones_topo.R
## -----------------------------------------------------------
##  Version corrigée (avec reprojection automatique)
## -----------------------------------------------------------

#' Crée des zones topographiques classées "Faible / Moyen / Fort" selon la
#' charge, et exporte table, raster et shapefile.
#'
#' @param fst_in     Chemin du .fst contenant les colonnes x, y, Charge, DOY …
#' @param up_tif     Raster gabarit (chemin vers .tif ou tout format lu par terra).
#' @param table_out  Fichier de sortie (.fst ou .csv) récapitulant les pixels par classe.
#' @param ras_out    Fichier raster INT1U (1 = Faible, 2 = Moyen, 3 = Fort).
#' @param shp_out    (Optionnel) shapefile des pixels classés.
#' @param win_perc   Fenêtre (± %) autour de la médiane pour sélectionner DAH et FSCA.
#' @param method     "quantile" (par défaut) ou "kmeans" pour les seuils de Charge.
#' @param thr_low_manual, thr_high_manual  Seuils imposés (outrepasse *method*).
#' @param input_crs  Code EPSG ou wkt des coordonnées x-y si ≠ du raster gabarit ;
#'                   sinon NULL (on suppose déjà la même projection).
#' @return `data.table` récapitulative (invisible) + fichiers écrits sur disque.
#' @import data.table terra sf fst ggplot2 patchwork
#' @examples
#' create_zones_topo("pts.fst", "gabarit.tif", "tab.fst", "classes.tif",
#'                   shp_out = "classes.shp", input_crs = "EPSG:3857")
create_zones_topo <- function(
    fst_in,
    up_tif,
    table_out,
    ras_out,
    shp_out    = NULL,
    win_perc   = 0.10,
    method     = c("quantile", "kmeans"),
    thr_low_manual  = NULL,
    thr_high_manual = NULL,
    input_crs  = NULL
) {
  
  ## ---------------------------------------------------------
  ## 1. Préliminaires & paquets ------------------------------
  ## ---------------------------------------------------------
  method <- match.arg(method)
  
  suppressPackageStartupMessages({
    library(data.table); library(fst); library(dplyr)
    library(terra);      library(sf);  library(ggplot2)
    library(patchwork);  library(fs)
  })
  
  ## ---------------------------------------------------------
  ## 2. Lecture des données ----------------------------------
  ## ---------------------------------------------------------
  dt <- as.data.table(read_fst(fst_in))
  for (nm in list(c("x","y"), c("X","Y")))
    if (all(nm %in% names(dt))) { setnames(dt, nm, c("x","y")); break }
  
  if (!all(c("x","y","Charge","DOY","dah","fsca","habitat_label") %in% names(dt)))
    stop("❌ Colonnes obligatoires manquantes dans " , fst_in)
  
  ## ---------------------------------------------------------
  ## 3. Raster gabarit & reprojection ------------------------
  ## ---------------------------------------------------------
  gabarit <- rast(up_tif, subds = 1)
  
  # -- reprojection si input_crs précisé et différent du raster
  if (!is.null(input_crs) && !identical(crs(gabarit), input_crs)) {
    v <- terra::vect(dt[, .(x, y)], type = "points", crs = input_crs)
    v <- terra::project(v, gabarit)
    new_xy <- terra::geom(v)[, c("x", "y")]
    dt[, `:=`(x = new_xy[,1], y = new_xy[,2])]
  }
  
  ## ---------------------------------------------------------
  ## 4. Calcul de l'indice de cellule ------------------------
  ## ---------------------------------------------------------
  dt[, cell := terra::cellFromXY(gabarit, cbind(x, y))]
  dt <- dt[!is.na(cell)]
  
  if (!nrow(dt))
    stop("❌ Aucune donnée à l'intérieur du raster gabarit (après reprojection éventuelle).")
  
  ## point de référence (premier jour de l'année)
  dt0 <- dt[DOY == min(DOY)]
  habitats <- dt0[, .N, by = habitat_label][order(-N)]
  
  ## ---------------------------------------------------------
  ## 5. Objets de sortie initiaux ----------------------------
  ## ---------------------------------------------------------
  class_r <- gabarit; values(class_r) <- NA_integer_
  codes   <- c(Faible = 1L, Moyen = 2L, Fort = 3L)
  res_list <- list()
  all_pix  <- NULL
  
  ## ---------------------------------------------------------
  ## 6. Boucle habitats --------------------------------------
  ## ---------------------------------------------------------
  for (hab in habitats$habitat_label) {
    
    sub <- dt0[habitat_label == hab]
    if (nrow(sub) < 100) next
    
    # -- sélection selon DAH & FSCA --------------------------
    med_dah  <- median(sub$dah ,  na.rm = TRUE)
    med_fsca <- median(sub$fsca, na.rm = TRUE)
    
    sub <- sub[dah  %between% c(med_dah*(1-win_perc), med_dah*(1+win_perc)) &
                 fsca %between% c(med_fsca*(1-win_perc), med_fsca*(1+win_perc))]
    if (nrow(sub) < 100) next
    
    # -- seuils Charge --------------------------------------
    if (!is.null(thr_low_manual) && !is.null(thr_high_manual)) {
      thr_low  <- thr_low_manual
      thr_high <- thr_high_manual
    } else if (method == "quantile") {
      q <- quantile(sub$Charge, probs = c(1/3, 2/3), na.rm = TRUE)
      thr_low  <- q[1];  thr_high <- q[2]
    } else {
      km <- kmeans(sub$Charge, centers = 3)
      m  <- sort(tapply(sub$Charge, km$cluster, mean))
      thr_low  <- mean(m[1:2]);  thr_high <- mean(m[2:3])
    }
    
    # -- attribution des classes -----------------------------
    sub[, classe :=
          fifelse(Charge <= thr_low , "Faible",
                  fifelse(Charge >  thr_high, "Fort",   "Moyen"))]
    
    if (!any(sub$classe %in% c("Faible","Moyen","Fort"))) next
    
    # -- table de synthèse -----------------------------------
    tab <- sub[, .N, by = classe]
    tab <- merge(data.table(classe = names(codes)), tab,
                 by = "classe", all.x = TRUE)
    tab[is.na(N), N := 0L]
    
    res_list[[hab]] <- data.table(
      habitat     = hab,
      thr_low     = thr_low,
      thr_high    = thr_high,
      pix_faible  = tab[classe=="Faible", N],
      pix_moyen   = tab[classe=="Moyen",  N],
      pix_fort    = tab[classe=="Fort",   N]
    )
    
    # -- raster & stats globales -----------------------------
    sub[, code := codes[classe]]
    class_r[sub$cell] <- sub$code
    
    all_pix <- rbindlist(list(all_pix,
                              sub[, .(classe, dah, fsca)]),
                         use.names = TRUE)
  }
  
  ## ---------------------------------------------------------
  ## 7. Contrôle final --------------------------------------
  ## ---------------------------------------------------------
  if (is.null(all_pix) || !nrow(all_pix))
    stop("❌ Aucun pixel classé ; ajustez win_perc ou vérifiez 'Charge'.")
  
  ## ---------------------------------------------------------
  ## 8. Écriture de la table --------------------------------
  ## ---------------------------------------------------------
  res_tab <- rbindlist(res_list)
  dir_create(path_dir(table_out), recurse = TRUE)
  if (grepl("\\.fst$", table_out, ignore.case = TRUE)) {
    fst::write_fst(res_tab, table_out)
  } else {
    data.table::fwrite(res_tab, table_out)
  }
  message("✓ Tableau écrit : ", table_out)
  
  ## ---------------------------------------------------------
  ## 9. Écriture du raster ----------------------------------
  ## ---------------------------------------------------------
  terra::writeRaster(class_r, ras_out, datatype = "INT1U", overwrite = TRUE)
  message("✓ Raster écrit  : ", ras_out, " (1 Faible, 2 Moyen, 3 Fort)")
  
  ## ---------------------------------------------------------
  ## 10. Shapefile (optionnel) ------------------------------
  ## ---------------------------------------------------------
  if (!is.null(shp_out)) {
    dir_create(path_dir(shp_out), recurse = TRUE)
    pts <- which(!is.na(values(class_r)))
    if (length(pts)) {
      xy  <- terra::xyFromCell(class_r, pts)
      coords_ok <- complete.cases(xy)
      if (any(coords_ok)) {
        sf_pts <- sf::st_as_sf(
          data.frame(xy[coords_ok, ],
                     code = values(class_r)[pts][coords_ok]),
          coords = c("x", "y"), crs = terra::crs(class_r)
        )
        sf::st_write(sf_pts, shp_out, delete_dsn = TRUE, quiet = TRUE)
        message("✓ Shapefile écrit : ", shp_out)
      } else {
        warning("Tous les points sont hors emprise ; shapefile non créé.")
      }
    } else {
      warning("Aucun pixel classé ; shapefile non créé.")
    }
  }
  
  ## ---------------------------------------------------------
  ## 11. Boxplots -------------------------------------------
  ## ---------------------------------------------------------
  pal <- c(Faible = "khaki3", Moyen = "tan", Fort = "steelblue")
  p1 <- ggplot2::ggplot(all_pix, ggplot2::aes(classe, dah,  fill = classe)) +
    ggplot2::geom_boxplot(width = .65) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "DAH", x = "")
  p2 <- ggplot2::ggplot(all_pix, ggplot2::aes(classe, fsca, fill = classe)) +
    ggplot2::geom_boxplot(width = .65) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "FSCA", x = "")
  print(p1 / p2)
  
  invisible(res_tab)
}




















#' Détection de zones homogènes de traitement à partir d'un .rds phéno
#'
#' Cette fonction reprend la méthode DBSCAN appliquée à un dataset phénologique
#' exporté en .rds avec build_pheno_dataset, pour identifier des zones homogènes
#' en DAH/FSCA et caractériser la charge.
#'
#' @param rds_path    Chemin vers le .rds généré par build_pheno_dataset
#' @param up_shape    Chemin vers le shapefile UP (EPSG:2154) ou NULL pour autodétecter
#' @param thr_high    Seuil haut de charge (>thr_high = Fort)
#' @param thr_low     Seuil bas de charge (≤thr_low = Faible)
#' @param eps_m       Distance epsilon en mètres pour DBSCAN
#' @param minpts      Nombre minimum de points pour former un cluster
#' @param hab_exclus  Vecteur de labels d'habitats à exclure
#' @param ras_template Chemin vers un raster servant de template (pour reprojeter)
#' @param csv_out     Chemin du CSV de sortie pour zones
#' @param ras_out     Chemin du GeoTIFF de sortie des zones codées
#' @param shp_out     Chemin du Shapefile de sortie des points zonés
#'
detect_zones_pheno <- function(
    rds_path,
    up_shape        = NULL,
    thr_high        = 300,
    thr_low         =  50,
    eps_m           = 1000,
    minpts          = 20,
    hab_exclus      = c(NA, "Formations minérales"),
    ras_template    = NULL,
    csv_out,
    ras_out,
    shp_out
) {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(dbscan)
    library(terra)
    library(sf)
  })
  
  # 1. Chargement .rds
  stopifnot(file.exists(rds_path))
  df <- readRDS(rds_path)
  if (!is.data.frame(df)) stop("Le .rds ne contient pas un data.frame.")
  # fixer en data.table pour vitesse
  dt <- as.data.table(df)
  
  # 2. Prétraitement coordonnées
  # Renommer colonnes si besoin
  for (p in list(c("x","y"), c("x.x","y.x"), c("x.y","y.y"),
                 c("X","Y"), c("lon","lat"), c("longitude","latitude"))) {
    if (all(p %in% names(dt))) { setnames(dt, p, c("x","y")); break }
  }
  
  # Extraire première date (variables statiques)
  dt_first <- dt[DOY == min(DOY)]
  
  # 3. Sélection des 3 habitats majeurs
  top3 <- dt_first[!habitat_label %in% hab_exclus, .N, by=habitat_label]
  top3 <- top3[order(-N)][1:3, habitat_label]
  
  # 4. Fonction interne de création de zones via DBSCAN
  create_zones <- function(sub) {
    med_dah  <- median(sub$dah, na.rm=TRUE)
    med_fsca <- median(sub$fsca, na.rm=TRUE)
    win_dah  <- 0.10 * abs(med_dah)
    win_fsca <- 0.10 * abs(med_fsca)
    cand <- sub[dah  %between% list(med_dah - win_dah,  med_dah + win_dah) &
                  fsca %between% list(med_fsca - win_fsca, med_fsca + win_fsca)]
    if (nrow(cand) < minpts) return(NULL)
    coords <- as.matrix(cand[, .(x,y)])
    cl <- dbscan(coords, eps=eps_m, minPts=minpts)
    cand[, cluster := cl$cluster]
    return(cand[cluster > 0])
  }
  
  # 5. Traitement par habitat
  zones_list <- vector("list", length(top3))
  names(zones_list) <- top3
  for (hab in top3) {
    message("• Habitat : ", hab)
    sub <- dt_first[habitat_label == hab]
    z <- create_zones(sub)
    if (is.null(z)) { warning("  aucune zone pour ", hab); next }
    # classer faible/fort
    z[, traitement := fcase(
      Charge > thr_high, "Fort",
      Charge <= thr_low,  "Faible",
      default = NA_character_
    )]
    z <- z[!is.na(traitement)]
    # garder clusters mixtes
    tab <- z[, .N, by=.(cluster, traitement)]
    wide <- dcast(tab, cluster ~ traitement, value.var="N", fill=0)
    if (!all(c("Fort","Faible") %in% names(wide))) next
    keep <- wide[Fort >= minpts & Faible >= minpts, cluster]
    z <- z[cluster %in% keep]
    if (nrow(z)==0) next
    # renumérotation top3 clusters
    top_cl <- z[, .N, by=cluster][order(-N)][1:3]
    top_cl[, zone_id := seq_len(.N)]
    z <- merge(z, top_cl, by="cluster")
    z[, habitat_label := hab]
    zones_list[[hab]] <- z[, .(cell, x, y, zone_id, habitat_label, traitement, dah, fsca)]
  }
  
  # 6. Export CSV
  zones_all <- rbindlist(zones_list, use.names=TRUE)
  if (nrow(zones_all)==0) stop("Aucune zone valide détectée.")
  fwrite(zones_all, csv_out)
  message("✓ CSV exporté : ", csv_out)
  
  # 7. Export raster (& shape)
  if (!is.null(ras_template) && file.exists(ras_template)) {
    tmpl <- rast(ras_template)
  } else {
    stop("Un raster template valide est requis.")
  }
  zone_r <- rast(ext=ext(tmpl), res=res(tmpl), crs=crs(tmpl))
  values(zone_r) <- NA_integer_
  code_map <- c(Faible=1L, Fort=2L)
  valid <- !is.na(zones_all$cell)
  zones_all[, code := code_map[traitement]]
  zone_r[zones_all$cell[valid]] <- zones_all$code[valid]
  writeRaster(zone_r, ras_out, datatype="INT1U", overwrite=TRUE)
  message("✓ Raster exporté : ", ras_out, " (1=Faible, 2=Fort)")
  
  # shapefile points
  pts <- st_as_sf(zones_all[valid, .(x,y,code)], coords=c("x","y"), crs=crs(zone_r))
  st_write(pts, shp_out, delete_dsn=TRUE, quiet=TRUE)
  message("✓ Shapefile exporté : ", shp_out)
  
  # 8. Boxplots empilés
  library(ggplot2); library(patchwork)
  plot_dt <- zones_all[!is.na(dah) & !is.na(fsca)]
  pal <- c(Faible="khaki3", Fort="steelblue")
  theme_base <- theme_bw(base_size=11) +
    theme(axis.title.x=element_blank(), legend.position="none",
          panel.grid.major.y=element_line(colour="grey85", size=.3),
          panel.grid.minor.y=element_blank())
  p1 <- ggplot(plot_dt, aes(traitement, dah, fill=traitement)) +
    geom_boxplot(width=.65, outlier.shape=21, outlier.size=1) +
    scale_fill_manual(values=pal) +
    labs(title="Distribution DAH", y="DAH") + theme_base
  p2 <- ggplot(plot_dt, aes(traitement, fsca, fill=traitement)) +
    geom_boxplot(width=.65, outlier.shape=21, outlier.size=1) +
    scale_fill_manual(values=pal) +
    labs(title="Distribution FSCA", y="FSCA") + theme_base
  print(p1 / p2)
  
  invisible(list(csv=csv_out, raster=ras_out, shapefile=shp_out))
}





























#' Détection de zones homogènes de traitement à partir d'un .rds phéno
#'
#' Méthode DBSCAN pour identifier des zones homogènes en DAH/FSCA,
#' puis enrichissement avec IRG et NDVI pour chaque pixel et DOY.
#'
#' @param rds_path     Chemin vers le .rds généré par build_pheno_dataset
#' @param up_shape     Chemin vers le shapefile UP (EPSG:2154) ou NULL pour autodétecter
#' @param thr_high     Seuil haut de charge (>thr_high = Fort)
#' @param thr_low      Seuil bas de charge (≤thr_low = Faible)
#' @param eps_m        Distance epsilon en mètres pour DBSCAN
#' @param minpts       Nombre minimum de points pour former un cluster
#' @param hab_exclus   Vecteur de labels d'habitats à exclure
#' @param ras_template Chemin vers un raster servant de template (ou NULL pour dériver)
#' @param csv_out      Chemin du CSV de sortie pour zones (+ DOY, IRG, NDVI)
#' @param ras_out      Chemin du GeoTIFF de sortie des zones codées
#' @param shp_out      Chemin du Shapefile de sortie des points zonés
#'
detect_zones_pheno <- function(
    rds_path,
    up_shape     = NULL,
    thr_high     = 300,
    thr_low      =  50,
    eps_m        = 1000,
    minpts       = 20,
    hab_exclus   = c(NA, "Formations minérales"),
    ras_template = NULL,
    csv_out,
    ras_out,
    shp_out
) {
  suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(dbscan)
    library(terra); library(sf)
  })
  
  # 1. Chargement .rds et conversion en data.table
  stopifnot(file.exists(rds_path))
  df <- readRDS(rds_path)
  if (!is.data.frame(df)) stop("Le .rds ne contient pas un data.frame.")
  dt <- as.data.table(df)
  
  # 2. Harmonisation coordinates -> x,y
  for (p in list(c("x","y"), c("x.x","y.x"), c("x.y","y.y"),
                 c("X","Y"), c("lon","lat"), c("longitude","latitude"))) {
    if (all(p %in% names(dt))) { setnames(dt, p, c("x","y")); break }
  }
  # Namespace DOY, IRG, NDVI
  if (!all(c("DOY","IRG","NDVI") %in% names(dt)))
    stop("Le data.frame doit contenir DOY, IRG et NDVI.")
  
  # 3. Extraction de la première date pour filtration habitat/charge
  dt_first <- dt[DOY == min(DOY)]
  
  # 4. Sélection des 3 habitats majeurs
  top3 <- dt_first[!habitat_label %in% hab_exclus, .N, by=habitat_label][order(-N)][1:3, habitat_label]
  message("Habitats étudiés : ", paste(top3, collapse=", "))
  
  # 5. Statistiques avant clustering
  message("Habitats étudiés : ", paste(top3, collapse=", "))
  
  # 5. Statistiques avant clustering
  message("Statistiques avant clustering (total / >",thr_high," / ≤",thr_low,"):")
  for (hab in top3) {
    sub0 <- dt_first[habitat_label == hab]
    tot  <- nrow(sub0)
    hi   <- sum(sub0$Charge > thr_high, na.rm=TRUE)
    lo   <- sum(sub0$Charge <= thr_low, na.rm=TRUE)
    message(sprintf("  %s : total=%d, >%d=%d, ≤%d=%d",
                    hab, tot, thr_high, hi, thr_low, lo))
  }
  
  # 6. Fonction interne DBSCAN sur DAH & FSCA
  create_zones <- function(sub) {
    med_dah  <- median(sub$dah, na.rm=TRUE)
    med_fsca <- median(sub$fsca, na.rm=TRUE)
    win_dah  <- 0.20 * abs(med_dah)
    win_fsca <- 0.20 * abs(med_fsca)
    cand <- sub[dah  %between% list(med_dah - win_dah, med_dah + win_dah) &
                  fsca %between% list(med_fsca - win_fsca, med_fsca + win_fsca)]
    if (nrow(cand) < minpts) return(NULL)
    cl <- dbscan(as.matrix(cand[,.(x,y)]), eps=eps_m, minPts=minpts)
    cand[, cluster := cl$cluster]
    return(cand[cluster>0])
  }
  
  # 7. Boucle habitat -> zones
  zones_list <- list()
  for (hab in top3) {
    message("\n• Habitat : ", hab)
    sub <- dt_first[habitat_label==hab]
    z <- create_zones(sub)
    if (is.null(z)) {
      message("  → Aucune zone homogène."); next
    }
    message("  Candidats post-DAH/FSCA : ", nrow(z))
    # classification
    z[, traitement := fcase(
      Charge > thr_high, "Fort",
      Charge <= thr_low,  "Faible",
      default = NA_character_)
    ]
    tab0 <- z[, .N, by=traitement]
    message("  Répartition Fort/Faible : ",
            paste(tab0$traitement, tab0$N, sep="=", collapse=", "))
    z <- z[!is.na(traitement)]
    # clusters mixtes
    tab <- z[, .N, by=.(cluster, traitement)]
    wide <- dcast(tab, cluster ~ traitement, value.var="N", fill=0)
    if (!all(c("Fort","Faible") %in% names(wide))) {
      message("  → Pas de cluster mixte."); next
    }
    keep <- wide[Fort>=minpts & Faible>=minpts, cluster]
    z <- z[cluster %in% keep]
    if (nrow(z)==0) { message("  → Plus aucune zone après mixte."); next }
    message("  Zones retenues : ", length(unique(z$cluster)))
    # top3 clusters
    top_cl <- z[, .N, by=cluster][order(-N)][1:3]
    top_cl[, zone_id:=seq_len(.N)]
    z <- merge(z, top_cl, by="cluster")
    z[, habitat_label:=hab]
    zones_list[[hab]] <- z[,.(cell,x,y,zone_id,habitat_label,traitement)]
  }
  
  zones_all <- rbindlist(zones_list, use.names=TRUE)
  if (nrow(zones_all)==0) stop("Aucune zone valide détectée.")
  
  # 8. Enrichissement NDVI & IRG & DOY pour chaque pixel
  enrich <- dt[, .(cell, DOY, IRG, NDVI)]
  zones_full <- merge(zones_all, enrich, by="cell", allow.cartesian=TRUE)
  
  # 9. Export CSV final
  fwrite(zones_full, csv_out)
  message("✓ CSV exporté (avec DOY, IRG, NDVI) : ", csv_out)
  
  # 10. Construction du raster template si non fourni
  if (!is.null(ras_template) && file.exists(ras_template)) {
    tmpl <- rast(ras_template)
  } else {
    # dériver à partir des coordonnées
    dx <- median(diff(sort(unique(zones_full$x))))
    dy <- median(diff(sort(unique(zones_full$y))))
    ext0 <- ext(min(zones_full$x)-dx/2, max(zones_full$x)+dx/2,
                min(zones_full$y)-dy/2, max(zones_full$y)+dy/2)
    tmpl <- rast(ext=ext0, res=c(dx,dy), crs="EPSG:2154")
    message("→ Raster template dérivé avec res=",dx,"×",dy)
  }
  
  # 11. Génération de la colonne code si manquante
  if (!"code" %in% names(zones_full)) {
    code_map <- c(Faible=1L, Fort=2L)
    zones_full[, code := code_map[traitement]]
  }
  
  # 12. Export raster & shapefile
  # Calcul des indices de cellules sur le template à partir des coordonnées
  idx <- terra::cellFromXY(tmpl, cbind(zones_full$x, zones_full$y))
  # Vérifier les indices valides
  valid_idx <- which(!is.na(idx) & idx >= 1 & idx <= terra::ncell(tmpl))
  if (length(valid_idx) == 0) stop("Aucun pixel valide pour le raster." )
  idx <- idx[valid_idx]
  codes <- zones_full$code[valid_idx]
  zone_r <- tmpl
  values(zone_r) <- NA_integer_
  zone_r[idx] <- codes
  writeRaster(zone_r, ras_out, datatype = "INT1U", overwrite = TRUE)
  message("✓ Raster exporté : ", ras_out)
  # Shapefile : on recrée à partir des pixels valides
  pts <- st_as_sf(zones_full[valid_idx, .(x,y,code)], coords = c("x","y"), crs = crs(tmpl))
  st_write(pts, shp_out, delete_dsn = TRUE, quiet = TRUE)
  message("✓ Shapefile exporté : ", shp_out)
  
  # 12. Boxplots DAH/FSCA Boxplots DAH/FSCA
  library(ggplot2); library(patchwork)
  plot_dt <- merge(zones_full, dt[,.(cell,dah,fsca)], by="cell")
  theme_base <- theme_bw(base_size=11) +
    theme(axis.title.x=element_blank(), legend.position="none",
          panel.grid.major.y=element_line(colour="grey85", size=.3),
          panel.grid.minor.y=element_blank())
  pal <- c(Faible="khaki3", Fort="steelblue")
  p1 <- ggplot(plot_dt, aes(traitement, dah, fill=traitement)) +
    geom_boxplot(width=.65, outlier.shape=21, outlier.size=1)+ scale_fill_manual(values=pal)+
    labs(title="DAH", y="DAH") + theme_base
  p2 <- ggplot(plot_dt, aes(traitement, fsca, fill=traitement)) +
    geom_boxplot(width=.65, outlier.shape=21, outlier.size=1)+ scale_fill_manual(values=pal)+
    labs(title="FSCA", y="FSCA") + theme_base
  print(p1 / p2)
  
  invisible(list(csv=csv_out, raster=ras_out, shapefile=shp_out))
}





# detect_zones_pheno.R  ─────────────────────────────────────────────────────────
# Version révisée (juillet 2025, v2)
#   • detect_zones_pheno()     : sélection ±10 % DAH/FSCA, toutes classes
#   • run_detect_zones()       : wrapper 1↔N alpages, filtrage des ...
# -----------------------------------------------------------------------------

#' Détection de zones homogènes de traitement (tous habitats, 3 classes)
#'
#' Pour **chaque habitat** présent dans le .rds :
#'   1. Cherche la "fenêtre" DAH/FSCA (± `delta_pct`) qui maximise le nombre de pixels ;
#'   2. Conserve ces pixels ;
#'   3. Classe la charge : **Fort > 450 UA**, **Modéré 100–300 UA**, **Faible < 50 UA** ;
#'   4. Garde l'habitat si ≥ `minpts` pixels dans ≥ 2 classes ;
#'   5. Exporte CSV, raster, shapefile.
#'
#' @param rds_path      Chemin vers le .rds phéno (fusionné ou non)
#' @param thr_fort      Seuil Fort   (default = 450)
#' @param thr_faible    Seuil Faible (default = 50)
#' @param thr_mod_low   Borne basse Modéré (default = 100)
#' @param thr_mod_high  Borne haute Modéré (default = 300)
#' @param delta_pct     Tolérance ± en % sur DAH & FSCA (default = 0.10)
#' @param minpts        Minimum de pixels requis dans AU MOINS 2 classes
#' @param csv_out / ras_out / shp_out  Chemins de sortie
#' @import data.table terra sf
#' @export

detect_zones_pheno <- function(rds_path,
                               thr_fort      = 450,
                               thr_faible    = 50,
                               thr_mod_low   = 100,
                               thr_mod_high  = 300,
                               delta_pct     = 0.10,
                               minpts        = 1,
                               ras_template  = NULL,
                               csv_out,
                               ras_out,
                               shp_out) {
  suppressPackageStartupMessages({
    library(data.table); library(terra); library(sf)
  })
  
  stopifnot(file.exists(rds_path))
  dt <- as.data.table(readRDS(rds_path))
  if (!nrow(dt)) stop("Le .rds est vide.")
  
  # Harmoniser x / y ---------------------------------------------------------
  for (p in list(c("x","y"), c("X","Y"), c("lon","lat"), c("longitude","latitude"))) {
    if (all(p %in% names(dt))) { setnames(dt, p, c("x","y")); break }
  }
  needed <- c("x","y","DOY","IRG","NDVI","dah","fsca","Charge","habitat_label")
  stopifnot(all(needed %in% names(dt)))
  
  dt_first <- dt[DOY == min(DOY)]
  habitats <- dt_first[, unique(habitat_label)]
  message("Nb d'habitats trouvés : ", length(habitats))
  
  # ---- Chercher la meilleure fenêtre pour un habitat ----------------------
  best_window <- function(sub, delta) {
    q_dah  <- quantile(sub$dah,  probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
    q_fsca <- quantile(sub$fsca, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
    grid   <- data.table::CJ(c_dah = q_dah, c_fsca = q_fsca)
    # Nombre de pixels dans la fenêtre ±delta % autour de chaque couple centre
    grid[, n := mapply(function(cd, cf) {
      sum(abs(sub$dah  - cd) <= delta * abs(cd) &
            abs(sub$fsca - cf) <= delta * abs(cf), na.rm = TRUE)
    }, c_dah, c_fsca)]
    grid[which.max(n), .(c_dah, c_fsca)]
  }
  
  zones_list <- list()
  for (hab in habitats) {
    sub <- dt_first[habitat_label == hab]
    if (nrow(sub) < minpts) {
      message("Habitat ", hab, " ignoré (trop peu de pixels)"); next
    }
    
    win <- best_window(sub, delta_pct)
    cand <- sub[ abs(dah  - win$c_dah ) <= delta_pct * abs(win$c_dah) &
                   abs(fsca - win$c_fsca) <= delta_pct * abs(win$c_fsca) ]
    if (nrow(cand) < minpts) {
      message("Habitat ", hab, " : fenêtre vide"); next
    }
    
    # Classification --------------------------------------------------------
    cand[, traitement := fifelse(Charge >  thr_fort,                       "Fort",
                                 fifelse(Charge <  thr_faible,                   "Faible",
                                         fifelse(Charge >= thr_mod_low & Charge <= thr_mod_high,
                                                 "Modere", NA_character_)))]
    cand <- cand[!is.na(traitement)]
    
    # Garder si ≥ 2 classes avec minpts ------------------------------------
    if (cand[, .N, by = traitement][N >= minpts, .N] < 2) {
      message("Habitat ", hab, " : pas assez d'échantillons dans ≥2 classes");
      next
    }
    zones_list[[hab]] <- cand[, .(cell, x, y, habitat_label = hab, traitement,
                                  DOY, IRG, NDVI)]
  }
  
  zones_full <- rbindlist(zones_list, use.names = TRUE)
  if (!nrow(zones_full)) stop("Aucune zone retenue au final.")
  
  # Export CSV --------------------------------------------------------------
  fwrite(zones_full, csv_out)
  message("✓ CSV exporté : ", csv_out)
  
  # Raster template ---------------------------------------------------------
  if (!is.null(ras_template) && file.exists(ras_template)) {
    tmpl <- rast(ras_template)
  } else {
    dx <- median(diff(sort(unique(zones_full$x))))
    dy <- median(diff(sort(unique(zones_full$y))))
    ext0 <- ext(min(zones_full$x) - dx/2, max(zones_full$x) + dx/2,
                min(zones_full$y) - dy/2, max(zones_full$y) + dy/2)
    tmpl <- rast(ext = ext0, res = c(dx, dy), crs = "EPSG:2154")
  }
  
  # Codes 1/2/3 -------------------------------------------------------------
  zones_full[, code := c(Faible = 1L, Modere = 2L, Fort = 3L)[traitement]]
  idx   <- terra::cellFromXY(tmpl, as.matrix(zones_full[, .(x, y)]))
  valid <- which(!is.na(idx))
  zone_r <- tmpl; values(zone_r) <- NA_integer_
  zone_r[idx[valid]] <- zones_full$code[valid]
  writeRaster(zone_r, ras_out, datatype = "INT1U", overwrite = TRUE)
  message("✓ Raster exporté : ", ras_out)
  
  pts <- st_as_sf(zones_full[valid, .(x, y, code)], coords = c("x", "y"), crs = crs(tmpl))
  st_write(pts, shp_out, delete_dsn = TRUE, quiet = TRUE)
  message("✓ Shapefile exporté : ", shp_out)
  
  invisible(list(csv = csv_out, raster = ras_out, shapefile = shp_out))
}


# ---------------------------------------------------------------------------
# Wrapper 1↔N alpages
# ---------------------------------------------------------------------------

run_detect_zones <- function(alpages, YEAR, output_dir, ...) {
  if (!length(alpages)) stop("'alpages' vide")
  
  # -----------------------------------------------------------------
  # Construire (ou pointer) le .rds à analyser
  # -----------------------------------------------------------------
  if (length(alpages) == 1L) {
    tag       <- alpages[1]
    fused_rds <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_", tag),
                           sprintf("dataset_pheno_LONG_%d_%s.rds", YEAR, tag))
  } else {
    tag       <- paste(alpages, collapse = "_")
    fused_rds <- file.path(output_dir, "10. NDVI Effect",
                           sprintf("dataset_pheno_LONG_%d_fus_%s.rds", YEAR, tag))
    if (!file.exists(fused_rds)) {
      message("🔄 Fusion des datasets : ", tag)
      dt_all <- data.table::rbindlist(lapply(alpages, function(alp) {
        p <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_", alp),
                       sprintf("dataset_pheno_LONG_%d_%s.rds", YEAR, alp))
        d <- readRDS(p); d$alpage <- alp; d
      }), use.names = TRUE, fill = TRUE)
      saveRDS(dt_all, fused_rds)
    }
  }
  
  # -----------------------------------------------------------------
  # Construire les chemins de sortie
  # -----------------------------------------------------------------
  csv_out <- file.path(output_dir, "10. NDVI Effect",
                       sprintf("dataset_filtered_%s_%d.csv", tag, YEAR))
  ras_out <- file.path(output_dir, "10. NDVI Effect",
                       sprintf("zonesTraitement_%s_%d.tif", tag, YEAR))
  shp_out <- file.path(output_dir, "10. NDVI Effect",
                       sprintf("zonesTraitement_%s_%d.shp", tag, YEAR))
  
  # -----------------------------------------------------------------
  # Filtrer les arguments ... pour ne garder que ceux acceptés par
  # detect_zones_pheno(), afin d'éviter l'erreur « arguments inutilisés ».
  # -----------------------------------------------------------------
  dotl   <- list(...)
  accept <- names(formals(detect_zones_pheno))
  dotl   <- dotl[names(dotl) %in% accept]
  
  do.call(detect_zones_pheno, c(list(rds_path = fused_rds,
                                     csv_out  = csv_out,
                                     ras_out  = ras_out,
                                     shp_out  = shp_out),
                                dotl))
}









# detect_zones_pheno.R  ─────────────────────────────────────────────────────────
# Version v5 (juillet 2025)
#   • Corrige l’ambiguïté des identifiants de pixels quand plusieurs
#     alpages sont fusionnés (gid = alpage_cell).
#   • Extraction temporelle complète inchangée, mais sans mélange entre
#     alpages. Plus de valeurs négatives/vides liées à la confusion.
#   • API (arguments) inchangée par rapport à v4.
# -----------------------------------------------------------------------------

# ===================================================================== #
# 1.  FONCTION  detect_zones_pheno()                                    #
# ===================================================================== #

#' Détection de zones homogènes + extraction NDVI/IRG (toutes dates)
#'
#' Identique à la v4, à ceci près :
#'   * Ajout d’un identifiant global `gid` = paste(alpage, cell, sep="_")
#'     (ou simplement `cell` si pas d’alpage) pour empêcher tout mélange
#'     entre alpages lors de la fusion.
#'   * Toutes les sélections / fusions utilisent désormais `gid` et non
#'     plus `cell` seul.
#'
#' Le CSV de sortie contient :
#'   gid, cell, x, y, habitat_label, traitement, alpage, DOY, IRG, NDVI
#'
#' @inheritParams detect_zones_pheno (v4)
#' @import data.table terra sf
#' @export

detect_zones_pheno <- function(rds_path,
                               thr_fort      = 800,
                               thr_faible    = 50,
                               thr_mod_low   = 300,
                               thr_mod_high  = 600,
                               delta_pct     = 0.10,
                               minpts        = 20,
                               ras_template  = NULL,
                               csv_out,
                               ras_out,
                               shp_out) {
  suppressPackageStartupMessages({
    library(data.table); library(terra); library(sf)
  })
  
  # ------------------------------------------------------------------
  # Lecture et vérifications de base
  # ------------------------------------------------------------------
  if (!file.exists(rds_path)) stop("Fichier .rds introuvable : ", rds_path)
  dt <- as.data.table(readRDS(rds_path))
  if (!nrow(dt)) stop("Le .rds est vide.")
  
  # Harmoniser x/y ----------------------------------------------------
  for (p in list(c("x","y"), c("X","Y"), c("lon","lat"), c("longitude","latitude"))) {
    if (all(p %in% names(dt))) { setnames(dt, p, c("x","y")); break }
  }
  
  needed <- c("cell","x","y","DOY","IRG","NDVI","dah","fsca","Charge","habitat_label")
  if (!all(needed %in% names(dt)))
    stop("Colonnes manquantes : ", paste(setdiff(needed, names(dt)), collapse=","))
  
  has_alpage <- "alpage" %in% names(dt)
  if (!has_alpage) dt[, alpage := NA_character_]
  
  # Identifiant global (unique, même en fusion multi‑alpages) ----------
  dt[, gid := fifelse(is.na(alpage), as.character(cell), paste(alpage, cell, sep="_"))]
  
  # Snapshot première date --------------------------------------------
  dt_first <- dt[DOY == min(DOY)]
  habitats <- unique(dt_first$habitat_label)
  message("Nb d'habitats trouvés : ", length(habitats))
  
  # Helper : classification Fort / Modéré / Faible --------------------
  charge_class <- function(v) {
    fifelse(v >  thr_fort,      "Fort",
            fifelse(v <  thr_faible,  "Faible",
                    fifelse(v >= thr_mod_low & v <= thr_mod_high, "Modere", NA_character_)))
  }
  
  zones_list <- list()
  for (hab in habitats) {
    sub <- dt_first[habitat_label == hab]
    if (nrow(sub) < minpts) {
      message("Habitat ", hab, " ignoré (", nrow(sub), " pixels)"); next }
    
    # Grille de centres (déciles) -------------------------------------
    q_dah  <- quantile(sub$dah,  probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
    q_fsca <- quantile(sub$fsca, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
    grid   <- data.table::CJ(c_dah = q_dah, c_fsca = q_fsca)
    
    best_score  <- -Inf
    best_subset <- NULL
    
    for (k in seq_len(nrow(grid))) {
      cd <- grid$c_dah[k]; cf <- grid$c_fsca[k]
      cand <- sub[ abs(dah  - cd) <= delta_pct * abs(cd) &
                     abs(fsca - cf) <= delta_pct * abs(cf) ]
      if (nrow(cand) < 3*minpts) next
      cand[, traitement := charge_class(Charge)]
      cand <- cand[!is.na(traitement)]
      ok <- cand[, .N, by = traitement][traitement %in% c("Fort","Modere","Faible") & N >= minpts, .N] == 3
      if (!ok) next
      if (nrow(cand) > best_score) { best_score <- nrow(cand); best_subset <- cand }
    }
    
    if (is.null(best_subset)) {
      message("Habitat ", hab, " : aucune fenêtre valide"); next }
    
    # Pixels retenus
    gids_sel <- unique(best_subset$gid)
    
    # Métadonnées (une ligne par pixel)
    meta <- unique(best_subset[, .(gid, cell, x, y,
                                   habitat_label = hab,
                                   traitement,
                                   alpage)])
    
    # Time‑series complètes NDVI / IRG --------------------------------
    ts <- dt[gid %in% gids_sel, .(gid, DOY, IRG, NDVI, alpage)]
    
    # Fusion méta × ts
    zones_list[[hab]] <- merge(meta, ts, by = c("gid","alpage"), allow.cartesian = TRUE)
  }
  
  zones_full <- rbindlist(zones_list, use.names = TRUE)
  if (!nrow(zones_full)) stop("Aucune zone retenue au final.")
  
  # ------------------------------------------------------------------
  # Export CSV (toutes dates)
  # ------------------------------------------------------------------
  fwrite(zones_full, csv_out)
  message("✓ CSV exporté : ", csv_out)
  
  # ------------------------------------------------------------------
  # Raster + Shapefile (une seule ligne par pixel)
  # ------------------------------------------------------------------
  zones_meta <- unique(zones_full[, .(x, y, traitement)])
  zones_meta[, code := c(Faible=1L, Modere=2L, Fort=3L)[traitement]]
  
  if (!is.null(ras_template) && file.exists(ras_template)) {
    tmpl <- rast(ras_template)
  } else {
    dx <- median(diff(sort(unique(zones_meta$x))))
    dy <- median(diff(sort(unique(zones_meta$y))))
    ext0 <- ext(min(zones_meta$x)-dx/2, max(zones_meta$x)+dx/2,
                min(zones_meta$y)-dy/2, max(zones_meta$y)+dy/2)
    tmpl <- rast(ext = ext0, res = c(dx,dy), crs = "EPSG:2154")
  }
  
  idx <- terra::cellFromXY(tmpl, as.matrix(zones_meta[,.(x,y)]))
  zone_r <- tmpl; values(zone_r) <- NA_integer_
  zone_r[idx] <- zones_meta$code
  writeRaster(zone_r, ras_out, datatype="INT1U", overwrite=TRUE)
  message("✓ Raster exporté : ", ras_out)
  
  pts <- st_as_sf(zones_meta, coords=c("x","y"), crs = crs(tmpl))
  st_write(pts, shp_out, delete_dsn = TRUE, quiet = TRUE)
  message("✓ Shapefile exporté : ", shp_out)
  
  invisible(list(csv = csv_out, raster = ras_out, shapefile = shp_out))
}

# ===================================================================== #
# 2.  WRAPPER  run_detect_zones() (inchangé)                            #
# ===================================================================== #

run_detect_zones <- function(alpages, YEAR, output_dir, ...) {
  if (!length(alpages)) stop("Argument 'alpages' vide")
  
  if (length(alpages) == 1L) {
    tag       <- alpages[1]
    fused_rds <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_", tag),
                           sprintf("dataset_pheno_LONG_%d_%s.rds", YEAR, tag))
  } else {
    tag       <- paste(alpages, collapse = "_")
    fused_rds <- file.path(output_dir, "10. NDVI Effect",
                           sprintf("dataset_pheno_LONG_%d_fus_%s.rds", YEAR, tag))
    if (!file.exists(fused_rds)) {
      message("🔄 Fusion des datasets : ", tag)
      dt_all <- data.table::rbindlist(lapply(alpages, function(alp) {
        p <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_", alp),
                       sprintf("dataset_pheno_LONG_%d_%s.rds", YEAR, alp))
        if (!file.exists(p)) stop(".rds manquant : ", p)
        d <- readRDS(p); d$alpage <- alp; d
      }), use.names = TRUE, fill = TRUE)
      saveRDS(dt_all, fused_rds)
    }
  }
  
  csv_out <- file.path(output_dir, "10. NDVI Effect",
                       sprintf("dataset_filtered_%s_%d.csv", tag, YEAR))
  ras_out <- file.path(output_dir, "10. NDVI Effect",
                       sprintf("zonesTraitement_%s_%d.tif", tag, YEAR))
  shp_out <- file.path(output_dir, "10. NDVI Effect",
                       sprintf("zonesTraitement_%s_%d.shp", tag, YEAR))
  
  dotl   <- list(...)
  accept <- names(formals(detect_zones_pheno))
  dotl   <- dotl[names(dotl) %in% accept]
  
  do.call(detect_zones_pheno,
          c(list(rds_path = fused_rds,
                 csv_out  = csv_out,
                 ras_out  = ras_out,
                 shp_out  = shp_out),
            dotl))
}






























































optimize_loading_windows2 <- function(data,
                                      habitat_col  = "habitat_label",
                                      dah_col      = "dah",
                                      fsca_col     = "fsca",
                                      charge_col   = "Charge",
                                      min_pixels   = c(faible = 1L, moyen = 1L, fort = 0L),
                                      win_percents = c(0.10, 0.15, 0.20, 0.30),
                                      grid_probs   = seq(0.05, 0.95, 0.05)) {
  
  library(dplyr)
  
  ## 1. Pré‑traitement ----
  df <- data %>%
    filter(!is.na(.data[[habitat_col]]),
           !is.na(.data[[dah_col]]),
           !is.na(.data[[fsca_col]]),
           !is.na(.data[[charge_col]])) %>%
    mutate(
      load_class = case_when(
        .data[[charge_col]] <  50                       ~ "faible",
        .data[[charge_col]] >= 250 & .data[[charge_col]] <= 350 ~ "moyen",
        .data[[charge_col]] > 1000                      ~ "fort",
        TRUE                                            ~ NA_character_
      )) %>%
    filter(!is.na(load_class))
  
  habitats <- sort(unique(df[[habitat_col]]))
  full_mask <- rep(FALSE, nrow(df))
  recap <- vector("list", length(habitats))
  
  ## 2. Boucle habitat ----
  for (h in seq_along(habitats)) {
    
    sub <- df %>% filter(.data[[habitat_col]] == habitats[h])
    
    if (nrow(sub) == 0) next
    
    dah_centers  <- quantile(sub[[dah_col]],  probs = grid_probs, na.rm = TRUE)
    fsca_centers <- quantile(sub[[fsca_col]], probs = grid_probs, na.rm = TRUE)
    
    best_total  <- -1L
    best_counts <- c(faible = 0L, moyen = 0L, fort = 0L)
    best_mask_h <- rep(FALSE, nrow(sub))
    best_dah_c  <- NA_real_
    best_fsca_c <- NA_real_
    best_win    <- NA_real_
    
    ## 2a. On explore plusieurs tailles de fenêtre (± 10 %, 15 %, …) -----------
    for (w in win_percents) {
      for (dc in dah_centers) {
        dah_low <- (1 - w) * dc; dah_high <- (1 + w) * dc
        dah_ok  <- sub[[dah_col]] >= dah_low & sub[[dah_col]] <= dah_high
        if (!any(dah_ok)) next
        for (fc in fsca_centers) {
          fsca_low <- (1 - w) * fc; fsca_high <- (1 + w) * fc
          mask <- dah_ok & sub[[fsca_col]] >= fsca_low & sub[[fsca_col]] <= fsca_high
          if (!any(mask)) next
          
          counts <- c(faible = 0L, moyen = 0L, fort = 0L)
          tab <- table(sub$load_class[mask])
          counts[names(tab)] <- as.integer(tab)
          total <- sum(counts)
          
          ## 2b. Critère : respecter les minima, sinon maximiser la somme -------
          valid <- all(counts >= min_pixels)
          better <- if (valid) {
            # priorité : fenêtre valide avec total max
            (best_total < total) || (best_total < 0)
          } else if (all(best_counts < min_pixels)) {
            # aucune fenêtre valide trouvée encore : maximiser total
            total > best_total
          } else {
            FALSE
          }
          
          if (better) {
            best_total  <- total
            best_counts <- counts
            best_mask_h <- mask
            best_dah_c  <- dc
            best_fsca_c <- fc
            best_win    <- w
          }
        }
      }
      
      ## si on vient de trouver une fenêtre *valide*, on arrête d’élargir
      if (all(best_counts >= min_pixels)) break
    }
    
    full_mask[which(df[[habitat_col]] == habitats[h])] <- best_mask_h
    
    recap[[h]] <- data.frame(
      habitat      = habitats[h],
      dah_center   = best_dah_c,
      fsca_center  = best_fsca_c,
      win_pct      = best_win,
      pixels_total = best_total,
      n_faible     = best_counts["faible"],
      n_moyen      = best_counts["moyen"],
      n_fort       = best_counts["fort"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  
  recap <- bind_rows(recap)
  
  ## 3. Affichage ----
  cat("\n=== Résumé sélection par habitat ===\n")
  print(recap, row.names = FALSE)
  cat("------------------------------------\n",
      "Pixels retenus:", sum(full_mask), "sur", nrow(df), "\n\n")
  
  invisible(df[full_mask, ])
}











build_pheno_dataset <- function(alpage, YEAR,
                                irg_tif, ndvi_tif,
                                habitat_tif, dah_tif, fsca_tif,
                                extra_tif,                  # ← déjà présent
                                charge_prev_rds,            # ← année n‑1
                                first_use_rds,              # ← année n
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
  
  ## 0. Vérification des fichiers -------------------------------------------
  needed <- c(irg_tif, ndvi_tif, habitat_tif, dah_tif, fsca_tif,
              extra_tif, charge_prev_rds, first_use_rds,
              lut_habitat_csv, up_shape)
  miss <- needed[!file.exists(needed)]
  if (length(miss))
    stop("Fichiers manquants :\n", paste(" •", miss, collapse = "\n"))
  dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
  
  ## 1. Grille de référence EPSG:2154 ---------------------------------------
  up_vect   <- vect(up_shape)                 # polygone UP
  meta_orig <- list(
    IRG  = rast(irg_tif)[[1]],
    NDVI = rast(ndvi_tif)[[1]],
    HAB  = rast(habitat_tif),
    DAH  = rast(dah_tif),
    FSCA = rast(fsca_tif),
    EXTRA= rast(extra_tif)[[1]]
  )
  base_res <- res(meta_orig[[which.max(sapply(meta_orig, \(r) prod(res(r))))]])
  template <- rast(ext = ext(up_vect), resolution = base_res, crs = crs(up_vect))
  
  read_align <- function(file, method = "bilinear")
    rast(file) |>
    project(template, method = method) |>
    crop(up_vect) |>
    mask(up_vect)
  
  irg_stack   <- read_align(irg_tif,   "bilinear")
  ndvi_stack  <- read_align(ndvi_tif,  "bilinear")
  habitat_r   <- read_align(habitat_tif,"near")
  dah_r       <- read_align(dah_tif,   "bilinear")
  fsca_r      <- read_align(fsca_tif,  "bilinear")
  extra_stack <- read_align(extra_tif, "bilinear")
  
  ## 2. Variables statiques --------------------------------------------------
  static_df <- as.data.frame(c(dah_r, fsca_r, habitat_r),
                             cells = TRUE, na.rm = FALSE) |>
    setNames(c("cell", "dah", "fsca", "habitat_code")) |>
    distinct(cell, .keep_all = TRUE)               # ← unicité
  
  extra_df <- as.data.frame(extra_stack, cells = TRUE, na.rm = FALSE) |>
    distinct(cell, .keep_all = TRUE)               # ← unicité
  
  static_df <- left_join(static_df, extra_df, by = "cell",
                         relationship = "one-to-one")  # sûr qu’il reste unique
  
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
    summarise(Charge = mean(Charge), .groups = "drop") # ← 1 ligne/pixel
  
  col_prev  <- paste0("Charge", YEAR - 1)
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
    left_join(static_df,  by = "cell", relationship = "many-to-one") |>
    left_join(charge_df,  by = "cell", relationship = "many-to-one") |>
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
























build_dataset_pheno_vec <- function(
    YEAR,
    alpage,
    pheno_tif_file,
    delta_max_tif_file,
    daily_rds_file,
    output_IRG_by_habitat,
    DOY_range   = 121:334,
    chunk_size  = 30,
    state_keep  = "Paturage",   # NULL → tous les états
    align_report_csv = NULL      # chemin CSV optionnel pour un rapport d’alignement
){
  library(terra)
  library(data.table)
  library(pbapply)
  
  ## 1) Charger rasters --------------------------------------------------------
  ph <- rast(pheno_tif_file)
  de <- rast(delta_max_tif_file)
  
  if (!(ext(ph) == ext(de) && all(res(ph) == res(de))))
    stop("pheno_tif_file et delta_max_tif_file n'ont pas la même emprise/résolution.")
  
  ## DOY d'après les noms des couches (tout ce qu'il y a après 'DOY')
  get_doys <- function(nms) suppressWarnings(as.integer(sub(".*DOY", "", nms)))
  avail_DOY_ph <- get_doys(names(ph))
  avail_DOY_de <- get_doys(names(de))
  
  if (anyNA(avail_DOY_ph)) stop("Impossible d'extraire les DOY depuis 'pheno_tif_file'.")
  if (anyNA(avail_DOY_de)) stop("Impossible d'extraire les DOY depuis 'delta_max_tif_file'.")
  
  ## DOY à garder = ordre de DOY_range, mais présents dans ph & de
  DOY_keep <- DOY_range[DOY_range %in% avail_DOY_ph & DOY_range %in% avail_DOY_de]
  if (!length(DOY_keep))
    stop("Aucun jour commun entre DOY_range et les stacks fournis.")
  dropped <- setdiff(DOY_range, DOY_keep)
  if (length(dropped))
    warning("Jours ignorés (absents d'au moins un stack) : ", paste(dropped, collapse = ", "))
  
  ## Indices par raster
  idx_ph_all <- match(DOY_keep, avail_DOY_ph)  # positions dans ph
  idx_de_all <- match(DOY_keep, avail_DOY_de)  # positions dans de
  
  ## Rapport d’alignement (optionnel)
  if (!is.null(align_report_csv)) {
    rep_dt <- data.table(
      doy      = DOY_keep,
      layer_ph = names(ph)[idx_ph_all],
      layer_de = names(de)[idx_de_all]
    )
    fwrite(rep_dt, align_report_csv)
    message("• Rapport d’alignement écrit : ", align_report_csv)
  }
  
  n_cells  <- ncell(ph)
  
  ## 2) Coordonnées ------------------------------------------------------------
  coords <- as.data.table(crds(ph, df = TRUE))[, cell := .I]
  
  ## 3) Charges (lookup) -------------------------------------------------------
  ch <- as.data.table(readRDS(daily_rds_file))
  
  # Renommer si présents
  if ("day"    %in% names(ch)) setnames(ch, "day",    "doy")
  if ("Charge" %in% names(ch)) setnames(ch, "Charge", "charge")
  if ("parc"   %in% names(ch)) setnames(ch, "parc",   "parc_nuit")
  
  if (!is.null(state_keep) && "state" %in% names(ch))
    ch <- ch[state == state_keep]
  
  ch[, cell := terra::cellFromXY(ph, cbind(x, y))]
  ch <- ch[!is.na(cell)]
  ch[, doy := as.integer(doy)]
  
  ch <- ch[, .(
    charge    = sum(charge, na.rm = TRUE),
    parc_nuit = paste(unique(parc_nuit), collapse = ";")
  ), by = .(cell, doy)]
  
  ch[, id := paste(cell, doy, sep = ":")]
  charge_vec <- setNames(ch$charge,    ch$id)
  parc_vec   <- setNames(ch$parc_nuit, ch$id)
  
  ## 4) Boucle en chunks (alignée) --------------------------------------------
  n_days   <- length(DOY_keep)
  n_chunks <- ceiling(n_days / chunk_size)
  res      <- vector("list", n_chunks)
  
  pbapply::pboptions(type = "timer")
  for (chunk in seq_len(n_chunks)) {
    s <- (chunk - 1L) * chunk_size + 1L
    e <- min(chunk * chunk_size, n_days)
    
    dyy   <- DOY_keep[s:e]
    ly_ph <- idx_ph_all[s:e]
    ly_de <- idx_de_all[s:e]
    
    # Extraction matrices (ncell × n_d) ; les 2 doivent matcher
    pheno_mat <- as.matrix(ph[[ly_ph]])
    delta_mat <- as.matrix(de[[ly_de]])
    
    if (!identical(dim(pheno_mat), dim(delta_mat)))
      stop("Incohérence dimensions (phéno vs delta) sur le chunk ", chunk,
           " : ", paste(dim(pheno_mat), collapse="x"), " vs ",
           paste(dim(delta_mat), collapse="x"))
    
    if (nrow(pheno_mat) != n_cells)
      stop("nrow != n_cells : ", nrow(pheno_mat), " vs ", n_cells)
    
    n_d <- ncol(pheno_mat)
    stopifnot(length(as.vector(delta_mat)) == n_cells * n_d)
    
    ## IMPORTANT : pas de rep.int ici
    dt <- data.table(
      cell               = rep(seq_len(n_cells), times = n_d),  # 1..ncell répété pour chaque DOY
      doy                = rep(dyy, each = n_cells),            # chaque DOY répété pour toutes les cellules
      delta_day_IRG_max  = as.vector(delta_mat),
      pheno_stage        = as.vector(pheno_mat)
    )
    
    # Ajout charge / parc_nuit via lookups
    dt[, id := paste(cell, doy, sep = ":")]
    dt[, charge    := charge_vec[id]]
    dt[, parc_nuit := parc_vec[id]]
    dt[, id := NULL]
    
    res[[chunk]] <- dt
  }
  
  ## 5) Assemblage & sauvegarde ------------------------------------------------
  dat <- rbindlist(res, use.names = TRUE, fill = FALSE)
  dat <- merge(dat, coords, by = "cell", all.x = TRUE, sort = FALSE)
  
  setcolorder(dat, c("cell","x","y","doy",
                     "delta_day_IRG_max","pheno_stage",
                     "charge","parc_nuit"))
  
  saveRDS(dat, output_IRG_by_habitat, compress = "xz")
  message("✓ Dataset écrit : ", output_IRG_by_habitat,
          "\n   → ", format(object.size(dat), units = "auto"))
}
