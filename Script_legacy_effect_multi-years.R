#----------------------------------------------------------------------------##
##------------------------- EFFET D'HERITAGE ---------------------------------##
##----------------------------------------------------------------------------##

## Description : 
## Ce script et l'objet de mon sujet de M2, qui vise a étudier l'impacte du cha-
## rgement de l'année n-1, sur le fonctionnement de l'écosystème l'année n. Pour
## cela on utilisera le PPI pour qulifier le fonctionnement de l'écosystème en 
## ciblant certaine variable d'intéret comme MAXV, le ONSET10, GREEN-UP-DURATION
## IRG .... Nous regarderons donc sur un gradient de paturage quelle est la part
## expliqué par les paramètres topo-climatique de l'utilisation. 
##
## Pour s'affranchir des effets du paturage directe de l'années nous allons étu-
## dier que les paramètres avant et dont MAXV, la majorité du troupeau y paturant
## que bien plus tard. Et donc en connaissant la date de première utilisation du 
## pixel par le troupeau nous filterons tout ceux dont la première utilisation et
## antérieur a MAXV. 

## Dans ce premier cardre (Version 1 | 22/07/2025). Nous étudirons 3 alpages :
## Viso, Cayolles, Sanguinieres (appelé Alpe-Sud, au court de se script), avec
## comme année de référence n = 2023 et comme année n-1 pour le chargement 2022.


## Second Cadre (Version 2 | 12/08/2025)

## PARAMETRE : 

alpage = "Cayolle"
YEAR = "2022"
YEARS = 2022


#### 0. Transformation des données ####
#-------------------------------------#
if(TRUE){
  
  
  
  
  if ("package:raster" %in% search()) detach("package:raster", unload = TRUE)
  library(terra)
  terra::focal(terra::rast(matrix(1:9,3,3)), w = 3, fun = median, na.rm = TRUE)
  
  
  
  
  
  
  ## SNOW
  # Descriprion : préparation du smod de l'année, transformation de DOY en année
  # hydrique (1septembre) en année courante (1 janvier)
  if (TRUE){
  for (YEAR in YEARS){
  ## LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_legacy.R"))
  ## ENTREE
    
  snow_case <- file.path(raster_dir, "SNOW")
  
  SMOD_tif_file <- file.path(snow_case , paste0("SMOD_", alpage, "_", YEAR, ".tif"))
  
  case_UP_file <- file.path(raster_dir, "UP")
  UP_file      <- file.path(case_UP_file, paste0("UP_", alpage, "_bis.shp"))
  
  
  ## CODE
  stack_path <- build_smod_fsca_stack(
    alpage         = alpage,
    YEAR           = YEAR,
    smod_hydro_tif = SMOD_tif_file,
    out_dir        = SMOD_case,
    up_shape       = UP_file,   # masque final uniquement
    smooth         = TRUE,
    kernel         = "gauss",
    radius_m       = 60,
    round_DOY      = TRUE,
    buffer_m       = 60
  )
  
  }
  }
  
  
  ## Chargement médian
  if (FALSE) {
    suppressPackageStartupMessages({ library(raster) })
    
    ## ENTREE
    # Alpage et années à traiter
    alpage <- "Sanguiniere"       # <- change ici
    YEARS  <- 2022:2024        # <- ou c(2021,2023,2024)
    method_resample <- "bilinear"  # "bilinear" (lisse) ou "ngb" (entiers)
    
    stopifnot(exists("output_dir"))
    
    build_path <- function(Y) file.path(
      output_dir,
      "5. Indicateurs_visualisation", "Taux_chargement",
      paste0(Y, "_", alpage),
      paste0("total_", Y, "_", alpage, ".tif")
    )
    tifs <- vapply(YEARS, build_path, character(1))
    missing <- YEARS[!file.exists(tifs)]
    if (length(missing) > 0) stop("Fichiers manquants pour années: ", paste(missing, collapse=","))
    
    ## SORTIE
    range_tag <- paste0(min(YEARS), "-", max(YEARS))
    out_dir <- file.path(output_dir, "5. Indicateurs_visualisation", "Chargement_representatif", alpage, range_tag)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(out_dir, paste0("chargement_median_", alpage, "_", range_tag, ".tif"))
    
    ## CODE
    ras_list <- lapply(tifs, raster)
    ref <- ras_list[[1]]
    if (length(ras_list) > 1) {
      for (i in 2:length(ras_list)) {
        ok <- try(compareRaster(ref, ras_list[[i]], extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, stopiffalse=FALSE), silent=TRUE)
        if (!isTRUE(ok)) {
          message(sprintf("[i] Resample %d -> ref (%s)...", YEARS[i], method_resample))
          ras_list[[i]] <- resample(ras_list[[i]], ref, method = method_resample)
        }
      }
    }
    
    s <- stack(ras_list)
    
    # Masque robustesse : on enlève les cellules jamais utilisées (toutes années = 0 ou NA)
    mask_used_any <- calc(s, function(v) if (all(is.na(v)) || sum(v, na.rm=TRUE)==0) NA else 1)
    s_mask <- mask(s, mask_used_any)
    
    # Médiane multi-années
    median_r <- calc(s_mask, function(x) median(x, na.rm = TRUE))
    
    # Export
    writeRaster(median_r, out_file, datatype = "FLT4S", overwrite = TRUE)
    message("[ok] Raster médian exporté : ", out_file)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}

#### 1. Préparation du jeu de données   ####
#------------------------------------------#
if (TRUE){
  # Description : 
  # Préparation du dataset alpage par alpage regroupent pour chaque pixel : 
  # - Chargement total de l'année n-1 (script_analysis | PERSEE_Traitement_Catlog)
  # - Le premier jour d'utilisation du pixel de l'année n, seuil de l'utilisation 
  #   aux chargement 30 (Script_visualization | PERSEE_Traitement_Catlog)
  # - Le DAH
  # - Le SMOD
  # - L'habitat
  # - Valeur de PPI (Plant Phenology Idex), par jour (DOY)
  # - Valeur de IRG (Instantaneous Rate Green), par jour (DOY)
  # - Les 10 paramètres de la double logistique + des paramètres dérivé
  
  # Nécéssite : 
  # - raster carto habitat
  # - raster de DAH (raster/Alti/)
  # - raster enneigement (raster/snow/)
  # - stackraster du NDVI (raster/Phenologie/)
  # - stackraster de l'IRG (raster/Phenologie/)
  # - Taux de chargement total (rds) (1er test a changer avec by day)
  #
  # En sortie : 
  #
  # Un dataset au format long pixel (pixel répété par DOY) enregistrement sous fst
  # nécéssitte le package associé (fst)
  


    
    ## LIBRARY
    source(file.path(functions_dir, "Functions_legacy_multi_year.R"))
    
    
    ## ENTREE
    alpage <- "Sanguiniere"                 # <- change ici
    YEARS  <- 2017:2023              # <- années d'étude
    MEDIAN_RANGE_TAG <- "2022-2024"  # <- raster médian FIXE (utilisé pour toutes les années)
    
    stopifnot(exists("output_dir"), exists("raster_dir"))
    
    # Fichiers statiques
    UP_file     <- file.path(raster_dir, "UP",   paste0("UP_",  alpage, ".shp"))
    DAH_tif     <- file.path(raster_dir, "Alti", paste0("DAH_1_", alpage, ".tif"))
    habitat_tif <- file.path(raster_dir, "Classifications_fusion_ColorIndexed_sc1_landforms_mnh.tif")
    lut_habitat_csv <- file.path(raster_dir, "class_habitat.csv")
    
    # Chargement médian (FIXE)
    median_tif <- file.path(output_dir, "5. Indicateurs_visualisation", "Chargement_representatif",
                            alpage, MEDIAN_RANGE_TAG,
                            paste0("chargement_median_", alpage, "_", MEDIAN_RANGE_TAG, ".tif"))
    
    # Fichiers annuels (SMOD/EXTRA)
    smod_dir  <- file.path(output_dir, "8. Analysis_Climate", "SMOD")
    pheno_dir <- file.path(raster_dir, "Phenologie")
    smod_path  <- function(Y) file.path(smod_dir,  paste0("SMOD_FSCA_stack_", alpage, "_", Y, "_epsg2154.tif"))
    extra_path <- function(Y) file.path(pheno_dir, paste0("PPI_extra_",       alpage, "_", Y, ".tif"))
    
    ## SORTIE
    analysis_tag <- paste0(min(YEARS), "-", max(YEARS))
    out_dir <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_", alpage))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_rds <- file.path(out_dir, paste0("dataset_legacy_multi-year_", alpage, "_", analysis_tag, ".rds"))
    
   
    ## CODE 
    ## FONCTION 1 :
    
    build_multi_year_dataset(
      alpage,
      YEARS,
      raster_dir,
      output_dir,
      median_range_tag = MEDIAN_RANGE_TAG
    )
    
  
  ## FONCTION 2 :
  ## Controle du jeu de données (PDF dans le outpout du dataset)
  # Plot de la zone d'étude 
  # Histogramme de : DAH, FSCA, NDVI, IRG
  
    pdf_out <- diag_dataset_legacy_multi_year(
      rds_path   = out_rds,
      up_shape   = UP_file ,
      sample_n   = 2e5,
      bins       = 50,
      output_dir = out_dir
    )
  
  ## FONCTION 3 : 
  # Cocaténation des jeu de données pour regrouper les alpages, dans notre cas :
  # agrégation de Viso, Cayolle et Sanguiniere donnent le dataset Alpe-Sud
  alpages <- c("Viso", "Cayolle", "Sanguiniere")
    
  merge_compact_dataset_years(
      years      = 2017:2023,
      alpages    = alpages,
      output_dir = file.path(getwd(), "outputs"),
      combo_name = "Alpe-Sud",
      write_fst  = FALSE
    )
  
}

#### 2. Affinage du jeu de données      ####
#------------------------------------------#
if (TRUE){
  # Description : 
  # Filtre des pixel utilisé avant la date du MAXV, utilisation donc : 
  # - first_day_use
  # - MAXV
  
  ## PARAMETRE
  alpage = "Alpe-Sud"
  YEARs = "2017-2023"
  
  
  # LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_legacy.R"))
  library(dplyr)
  
  
  
  ## ENTREE
  out_dir <- file.path(output_dir, "10. NDVI Effect",
                       paste0("dataset_", alpage))
  
  dt_alpesud = file.path(out_dir, paste0("dataset_legacy_multi-year_", alpage, "_", YEARs, ".rds"))
  
  ## SORTIE
  dt_filtered = file.path(out_dir, paste0("df_filtered_legacy_multi-year_", YEARs, "_", alpage, ".rds"))

  
  
  
  ## CODE 
  
  # Filtre des pixel avec du minéral:
  if (TRUE){
    df_filtre <- readRDS(dt_alpesud) %>% 
      dplyr::filter(habitat_code != 18) %>% 
      dplyr::filter(!is.na(GPROD), !is.na(DAH), !is.na(SMOD), !is.na(`chargement_median_2022-2024`)) %>% 
      mutate(
        Charge_log = log1p(`chargement_median_2022-2024`),
        charge_bin = as.integer(`chargement_median_2022-2024` >= 1),)          # log(Charge + 1)
      
      
    saveRDS(df_filtre, dt_filtered) %>% 
      cat("Save dataset :", dt_filtered, "\n")
  }
  

}




#### 3. Analyse                         ####
#------------------------------------------#

if (TRUE ) {
  #### 3.1 Analyse préparatoire ###
  if (TRUE){
    ##Description :
    # distribution des données, analyse du chargement : quelle transformation
    # appliquer ? traitement du paturage en gradient ou classe ? 
    
    ## PARAMETRE
    alpage = "Alpe-Sud"
    YEARs = "2017-2023"
    
    # LIBRARY & FUNCTION
    source(file.path(functions_dir, "Functions_legacy.R"))
    library(ggplot2)
    library(dplyr)
    
    ## ENTREE
    out_dir <- file.path(output_dir, "10. NDVI Effect",
                         paste0("dataset_", alpage))
    
    dt_filtered = file.path(out_dir, paste0("df_filtered_legacy_multi-year_", YEARs, "_", alpage, ".rds"))
    
    dt_legacy = readRDS(dt_filtered)

    
    
    
    # Chargement : 
    ggplot(dt_legacy, aes(x = Charge_log)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~ YEAR, scales = "free_y") +
      labs(title = "Log du chargement par habitat")
    
    
    # DAH :
    ggplot(dt_legacy, aes(x = DAH)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~YEAR, scales = "free_y") +
      labs(title = "DAH par an")
    
        
    # SMOD :
    ggplot(dt_legacy, aes(x = SMOD)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~YEAR, scales = "free_y") +
      labs(title = "SMOD par an")
    
    
    # GPROD
    ggplot(dt_legacy, aes(x = GPROD)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~YEAR, scales = "free_y") +
      labs(title = "Production par an")

    
    
    
    
    
    ### PLOT RESULTAT 1.2 : Distribution SMOD et HABITAT
    
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(grid)
    
    # ---------- 0) Filtre + nettoyage ----------
    hab_focus      <- c(1, 9, 5, 6)
    hab_focus_chr  <- as.character(hab_focus)
    
    # Filtrer l'année (adapter si besoin)
    dt_plot <- dt_legacy %>%
      dplyr::filter(YEAR == 2022) %>%
      mutate(
        # harmoniser le type + l'ordre des habitats pour coller à la palette nommée
        habitat_code = factor(as.character(habitat_code),
                              levels = hab_focus_chr)
      )
    
    # Outliers (IQR) sur SMOD
    q   <- stats::quantile(dt_plot$SMOD, c(.25, .75), na.rm = TRUE)
    iqr <- diff(q); lo <- q[1] - 1.5*iqr; hi <- q[2] + 1.5*iqr
    dt_plot <- dt_plot %>% dplyr::filter(dplyr::between(SMOD, lo, hi))
    
    # ---------- 1) Bins communs + comptes ----------
    n_bins <- 30
    hb     <- hist(dt_plot$SMOD, breaks = n_bins, plot = FALSE)
    brks   <- hb$breaks
    binw   <- diff(brks)[1]
    
    # Étendue utile : de la 1re à la dernière barre non vide
    nz <- which(hb$counts > 0)
    xmin_use <- brks[min(nz)]
    xmax_use <- brks[max(nz) + 1]
    
    edges <- tibble::tibble(
      bin_id = seq_len(length(brks) - 1),
      xmin   = brks[-length(brks)],
      xmax   = brks[-1],
      xmid   = (xmin + xmax)/2
    )
    
    # Comptes globaux par bin
    bins_total <- dt_plot %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed = TRUE, all.inside = TRUE)) %>%
      count(bin_id, name = "count") %>%
      right_join(edges, by = "bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L))
    
    # Comptes par habitat d’intérêt (utiliser niveaux factorisés)
    bins_hab <- dt_plot %>%
      filter(habitat_code %in% hab_focus_chr) %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed = TRUE, all.inside = TRUE)) %>%
      count(habitat_code, bin_id, name = "count") %>%
      right_join(edges, by = "bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L)) %>%
      arrange(habitat_code, xmin)
    
    # Données pour contours "step" alignés sur bords des bins
    hab_step <- bins_hab %>%
      select(habitat_code, xmin, xmax, count) %>%
      pivot_longer(c(xmin, xmax), names_to = "edge", values_to = "x") %>%
      arrange(habitat_code, x) %>%
      transmute(habitat_code, x, y = count)
    
    # Seuils des 3 blocs (terciles) + snap sur brks
    bornes_brutes  <- stats::quantile(dt_plot$SMOD, c(0, 1/3, 2/3, 1), na.rm = TRUE)
    snap_to_break  <- function(v, brks) brks[which.min(abs(brks - v))]
    bornes         <- sort(unique(sapply(bornes_brutes, snap_to_break, brks = brks)))
    # bornes = c(xmin_classe1, coupure1, coupure2, xmax_classe3)
    
    # Graduations X identiques
    tick_by      <- 5
    x_breaks_all <- seq(floor(xmin_use/tick_by)*tick_by,
                        ceiling(xmax_use/tick_by)*tick_by, by = tick_by)
    
    # ---------- 2) Hauteurs/bandeau ----------
    bar_max     <- max(bins_total$count)
    strip_ratio <- 0.12
    gap_ratio   <- 0.05
    strip_h_abs <- bar_max * strip_ratio
    y_top       <- bar_max * (1 + gap_ratio) + strip_h_abs
    strip_ymin  <- y_top - strip_h_abs
    
    # Couleurs lisibles pour les 4 habitats (étiquettes = niveaux)
    cols_focus <- c(
      "1" = "#3B5BDB", # indigo
      "9" = "#E03131", # red
      "5" = "#F4A261", # orange
      "6" = "#1B9E77"  # teal
    )
    
    # Breaks sans répétition aux frontières
    panel_breaks <- function(xmin, xmax, side = c("left","middle","right")) {
      side <- match.arg(side)
      b <- x_breaks_all[x_breaks_all >= xmin & x_breaks_all <= xmax]
      if (side == "left")   b <- b[b <  xmax]
      if (side == "middle") b <- b[b >  xmin & b < xmax]
      if (side == "right")  b <- b[b >  xmin]
      b
    }
    
    # ---------- 3) Panneau avec bandeau collé tout en haut ----------
    make_panel <- function(xmin, xmax, titre, side, show_y = FALSE, show_x = FALSE) {
      ggplot() +
        # fond global
        geom_col(data = bins_total, aes(x = xmid, y = count),
                 width = binw, fill = "grey88", colour = "grey70", linewidth = 0.25) +
        # habitats d'intérêt (remplissage + contours)
        geom_col(data = bins_hab, aes(x = xmid, y = count, fill = habitat_code),
                 width = binw, position = "identity", alpha = 0.30, colour = NA) +
        geom_path(data = hab_step, aes(x = x, y = y, colour = habitat_code, group = habitat_code),
                  linewidth = 1.1, lineend = "butt") +
        
        # limites Y
        coord_cartesian(xlim = c(xmin, xmax), ylim = c(0, y_top), clip = "on") +
        
        # BANDEAU collé en haut
        annotate("rect", xmin = xmin, xmax = xmax, ymin = strip_ymin, ymax = Inf,
                 fill = "#E1E5EB", colour = NA) +
        annotate("segment", x = xmin, xend = xmax, y = strip_ymin, yend = strip_ymin,
                 colour = "black", linewidth = 1.1) +
        annotate("text", x = (xmin + xmax)/2, y = (strip_ymin + y_top)/2,
                 label = titre, fontface = "bold", size = 4.4, colour = "#222222") +
        
        # Axes & légende
        scale_x_continuous(breaks = panel_breaks(xmin, xmax, side = side),
                           minor_breaks = NULL, expand = expansion(add = 0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0))) +
        scale_fill_manual(values = cols_focus, name = "Habitats d’intérêt", drop = FALSE) +
        scale_colour_manual(values = cols_focus, guide = "none", drop = FALSE) +
        labs(x = if (show_x) NULL else NULL,
             y = if (show_y) "Effectif" else NULL) +
        theme_minimal(base_family = "Arial") +
        theme(
          panel.border       = element_rect(colour = "black", fill = NA, linewidth = 1.6),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.x       = element_text(margin = margin(t = 6)),
          axis.text.y        = if (show_y) element_text() else element_blank(),
          axis.ticks.y       = if (show_y) element_line() else element_blank(),
          plot.margin        = margin(0, -3, 0, -3, unit = "pt"),
          legend.position    = "right",
          legend.direction   = "vertical",
          legend.title       = element_text(face = "bold")
        )
    }
    
    # ---------- 4) Construire les 3 blocs ----------
    p1 <- make_panel(bornes[1], bornes[2], "Déneigement précoce", side = "left",  show_y = TRUE)  +
      theme(plot.margin = margin(0, -3, 0,  0, "pt"))
    p2 <- make_panel(bornes[2], bornes[3], "Déneigement moyen",   side = "middle") +
      theme(plot.margin = margin(0, -3, 0, -3, "pt"))
    p3 <- make_panel(bornes[3], bornes[4], "Déneigement tardif",  side = "right")  +
      theme(plot.margin = margin(0,  0, 0, -3, "pt"))
    
    # Largeurs proportionnelles à l’étendue de chaque classe
    w <- diff(bornes)
    
    # ---------- 5) Titre d’axe X global centré ----------
    year_label <- unique(dt_plot$YEAR)
    year_label <- if (length(year_label)) year_label[1] else "?"
    x_title <- ggplot() +
      theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = paste0("SMOD ", year_label),
               fontface = "bold", size = 4.8, colour = "#222222") +
      theme(plot.margin = margin(6, 0, 0, 0))
    
    # ---------- 6) Assemblage final ----------
    p_row <- (p1 | p2 | p3) + plot_layout(widths = w, guides = "collect")
    final_plot <- (p_row / x_title) +
      plot_layout(heights = c(1, 0.06)) +
      plot_annotation(
        title = paste0("SMOD ", year_label, " — histogramme global + habitats d’intérêt (3 blocs segmentés)"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.margin = margin(0, 0, 0, 0)
        )
      )
    
    # Affichage
    final_plot
    
    # Export PNG fond blanc
    ggsave(
      filename = file.path(inputs_case, paste0("SMOD_", year_label, "_HAB_segments.png")),
      plot     = final_plot,
      width    = 12, height = 5.5, dpi = 300, bg = "white"
    )
    
    
    
    
    
    ## V2 
    
    ### PLOT RESULTAT 1.2 : Distribution SMOD et HABITAT (légende par noms + rappel couleurs blocs)
    
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(grid)
    
    # ---------- 0) Filtre + nettoyage ----------
    hab_focus      <- c(1, 9, 5, 6)
    hab_focus_chr  <- as.character(hab_focus)
    
    # Mapping codes -> NOMS d'habitats (dans l'ordre voulu)
    hab_labels <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    
    # Palette pour les NOMS d'habitats
    cols_focus <- c(
      "P. nivales"                  = "#3B5BDB",  # indigo
      "P. thermiques écorchées"     = "#E03131",  # rouge
      "Queyrellins"                 = "#F4A261",  # orange
      "Nardaies denses du subalpin" = "#1B9E77"   # teal
    )
    
    # Couleurs des BLOCS "Déneigement"
    block_cols <- list(
      precoce = "#fed976",
      moyen   = "#41b6c4",
      tardif  = "#0c2c84"
    )
    
    # Filtrer l'année (adapter si besoin) et harmoniser habitat
    dt_plot <- dt_legacy %>%
      dplyr::filter(YEAR == 2023) %>%
      mutate(
        habitat_code   = as.character(habitat_code),
        habitat_label  = recode(habitat_code, !!!hab_labels),
        habitat_label  = factor(habitat_label, levels = names(cols_focus))
      )
    
    # Outliers (IQR) sur SMOD
    q   <- stats::quantile(dt_plot$SMOD, c(.25, .75), na.rm = TRUE)
    iqr <- diff(q); lo <- q[1] - 1.5*iqr; hi <- q[2] + 1.5*iqr
    dt_plot <- dt_plot %>% dplyr::filter(dplyr::between(SMOD, lo, hi))
    
    # ---------- 1) Bins communs + comptes ----------
    n_bins <- 25
    hb     <- hist(dt_plot$SMOD, breaks = n_bins, plot = FALSE)
    brks   <- hb$breaks
    binw   <- diff(brks)[1]
    
    # Étendue utile : de la 1re à la dernière barre non vide
    nz <- which(hb$counts > 0)
    xmin_use <- brks[min(nz)]
    xmax_use <- brks[max(nz) + 1]
    
    edges <- tibble::tibble(
      bin_id = seq_len(length(brks) - 1),
      xmin   = brks[-length(brks)],
      xmax   = brks[-1],
      xmid   = (xmin + xmax)/2
    )
    
    # Comptes globaux par bin
    bins_total <- dt_plot %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed = TRUE, all.inside = TRUE)) %>%
      count(bin_id, name = "count") %>%
      right_join(edges, by = "bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L))
    
    # Comptes par habitat (avec labels)
    bins_hab <- dt_plot %>%
      filter(habitat_label %in% names(cols_focus)) %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed = TRUE, all.inside = TRUE)) %>%
      count(habitat_label, bin_id, name = "count") %>%
      right_join(edges, by = "bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L)) %>%
      arrange(habitat_label, xmin)
    
    # Données pour contours "step" alignés sur bords des bins
    hab_step <- bins_hab %>%
      select(habitat_label, xmin, xmax, count) %>%
      pivot_longer(c(xmin, xmax), names_to = "edge", values_to = "x") %>%
      arrange(habitat_label, x) %>%
      transmute(habitat_label, x, y = count)
    
    # Seuils des 3 blocs (terciles) + snap sur brks
    bornes_brutes  <- stats::quantile(dt_plot$SMOD, c(0, 1/3, 2/3, 1), na.rm = TRUE)
    snap_to_break  <- function(v, brks) brks[which.min(abs(brks - v))]
    bornes         <- sort(unique(sapply(bornes_brutes, snap_to_break, brks = brks)))
    # bornes = c(xmin_classe1, coupure1, coupure2, xmax_classe3)
    
    # Graduations X identiques
    tick_by      <- 5
    x_breaks_all <- seq(floor(xmin_use/tick_by)*tick_by,
                        ceiling(xmax_use/tick_by)*tick_by, by = tick_by)
    
    # ---------- 2) Hauteurs/bandeau ----------
    bar_max     <- max(bins_total$count)
    strip_ratio <- 0.12
    gap_ratio   <- 0.05
    strip_h_abs <- bar_max * strip_ratio
    y_top       <- bar_max * (1 + gap_ratio) + strip_h_abs
    strip_ymin  <- y_top - strip_h_abs
    
    # Breaks sans répétition aux frontières
    panel_breaks <- function(xmin, xmax, side = c("left","middle","right")) {
      side <- match.arg(side)
      b <- x_breaks_all[x_breaks_all >= xmin & x_breaks_all <= xmax]
      if (side == "left")   b <- b[b <  xmax]
      if (side == "middle") b <- b[b >  xmin & b < xmax]
      if (side == "right")  b <- b[b >  xmin]
      b
    }
    
    # ---------- 3) Panneau avec bandeau coloré (rappel bloc) ----------
    make_panel <- function(xmin, xmax, titre, side, band_fill, show_y = FALSE, show_x = FALSE) {
      ggplot() +
        # fond global
        geom_col(data = bins_total, aes(x = xmid, y = count),
                 width = binw, fill = "grey88", colour = "grey70", linewidth = 0.25) +
        # habitats d'intérêt (remplissage + contours)
        geom_col(data = bins_hab, aes(x = xmid, y = count, fill = habitat_label),
                 width = binw, position = "identity", alpha = 0.30, colour = NA) +
        geom_path(data = hab_step, aes(x = x, y = y, colour = habitat_label, group = habitat_label),
                  linewidth = 1.1, lineend = "butt") +
        
        # limites Y
        coord_cartesian(xlim = c(xmin, xmax), ylim = c(0, y_top), clip = "on") +
        
        # BANDEAU collé en haut (couleur du bloc, discret)
        annotate("rect", xmin = xmin, xmax = xmax, ymin = strip_ymin, ymax = Inf,
                 fill = band_fill, alpha = 0.25, colour = NA) +
        annotate("segment", x = xmin, xend = xmax, y = strip_ymin, yend = strip_ymin,
                 colour = "black", linewidth = 1.1) +
        annotate("text", x = (xmin + xmax)/2, y = (strip_ymin + y_top)/2,
                 label = titre, fontface = "bold", size = 4.4, colour = "#222222") +
        
        # Axes & légende
        scale_x_continuous(breaks = panel_breaks(xmin, xmax, side = side),
                           minor_breaks = NULL, expand = expansion(add = 0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0))) +
        scale_fill_manual(values = cols_focus, name = "Habitat", drop = FALSE) +
        scale_colour_manual(values = cols_focus, guide = "none", drop = FALSE) +
        labs(x = if (show_x) NULL else NULL,
             y = if (show_y) "Effectif" else NULL) +
        theme_minimal(base_family = "Arial") +
        theme(
          panel.border       = element_rect(colour = "black", fill = NA, linewidth = 1.6),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.x       = element_text(margin = margin(t = 6)),
          axis.text.y        = if (show_y) element_text() else element_blank(),
          axis.ticks.y       = if (show_y) element_line() else element_blank(),
          plot.margin        = margin(0, -3, 0, -3, unit = "pt"),
          legend.position    = "right",
          legend.direction   = "vertical",
          legend.title       = element_text(face = "bold")
        )
    }
    
    # ---------- 4) Construire les 3 blocs ----------
    p1 <- make_panel(bornes[1], bornes[2], "Précoce", side = "left",
                     band_fill = block_cols$precoce, show_y = TRUE)  +
      theme(plot.margin = margin(0, -3, 0,  0, "pt"))
    p2 <- make_panel(bornes[2], bornes[3], "Moyen",   side = "middle",
                     band_fill = block_cols$moyen) +
      theme(plot.margin = margin(0, -3, 0, -3, "pt"))
    p3 <- make_panel(bornes[3], bornes[4], "Tardif",  side = "right",
                     band_fill = block_cols$tardif)  +
      theme(plot.margin = margin(0,  0, 0, -3, "pt"))
    
    # Largeurs proportionnelles à l’étendue de chaque classe
    w <- diff(bornes)
    
    # ---------- 5) Titre d’axe X global centré ----------
    year_label <- unique(dt_plot$YEAR)
    year_label <- if (length(year_label)) year_label[1] else "?"
    x_title <- ggplot() +
      theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = paste0("SMOD ", year_label),
               fontface = "bold", size = 4.8, colour = "#222222") +
      theme(plot.margin = margin(6, 0, 0, 0))
    
    # ---------- 6) Assemblage final ----------
    p_row <- (p1 | p2 | p3) + plot_layout(widths = w, guides = "collect")
    final_plot <- (p_row / x_title) +
      plot_layout(heights = c(1, 0.06)) +
      plot_annotation(
        title = paste0("SMOD ", year_label, " — histogramme global + habitats d’intérêt (3 blocs segmentés)"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.margin = margin(0, 0, 0, 0)
        )
      )
    
    # Affichage
    final_plot
    
    # Export PNG fond blanc
    ggsave(
      filename = file.path(out_dir, paste0("SMOD_", year_label, "_HAB_segments.png")),
      plot     = final_plot,
      width    = 10, height = 8, dpi = 300, bg = "white"
    )
    
    
    
    
    
    
    # === V3 (axes lisibles) + FIX Y=11000 constant & sans tick 10000 ===
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(grid)
    
    # -------- Paramètres globaux --------
    year_to_plot <- 2022
    x_min_force  <- 15
    x_max_force  <- 190
    n_bins       <- 25
    brks         <- seq(x_min_force, x_max_force, length.out = n_bins + 1)
    binw         <- diff(brks)[1]
    x_breaks_all <- seq(20, 190, by = 10)
    
    # Habitats (ordre & couleurs)
    hab_labels <- c("1"="P. nivales","9"="P. thermiques écorchées",
                    "5"="Queyrellins","6"="Nardaies denses du subalpin")
    cols_focus <- c("P. nivales"="#3B5BDB","P. thermiques écorchées"="#E03131",
                    "Queyrellins"="#F4A261","Nardaies denses du subalpin"="#1B9E77")
    
    # Bandeaux (légère transparence)
    block_cols <- list(precoce="#fed976", moyen="#41b6c4", tardif="#0c2c84")
    band_alpha <- 0.4
    
    # -------- Données --------
    dt_plot <- dt_legacy %>%
      filter(YEAR == year_to_plot) %>%
      mutate(habitat_code  = as.character(habitat_code),
             habitat_label = recode(habitat_code, !!!hab_labels, .default = NA_character_),
             habitat_label = factor(habitat_label, levels = names(cols_focus))) %>%
      filter(SMOD >= x_min_force, SMOD <= x_max_force)
    
    # (optionnel) filtre IQR
    q <- quantile(dt_plot$SMOD, c(.25,.75), na.rm=TRUE); iqr <- diff(q)
    dt_plot <- dt_plot %>% filter(between(SMOD, q[1]-1.5*iqr, q[2]+1.5*iqr))
    
    # -------- Bins & comptes --------
    edges <- tibble::tibble(
      bin_id = seq_len(length(brks)-1),
      xmin   = brks[-length(brks)],
      xmax   = brks[-1],
      xmid   = (brks[-length(brks)] + brks[-1]) / 2
    )
    bins_total <- dt_plot %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed=TRUE, all.inside=TRUE)) %>%
      count(bin_id, name="count") %>%
      right_join(edges, by="bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L))
    
    bins_hab <- dt_plot %>%
      filter(!is.na(habitat_label)) %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed=TRUE, all.inside=TRUE)) %>%
      count(habitat_label, bin_id, name="count") %>%
      right_join(edges, by="bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L)) %>%
      arrange(habitat_label, xmin)
    
    hab_step <- bins_hab %>%
      select(habitat_label, xmin, xmax, count) %>%
      pivot_longer(c(xmin, xmax), names_to="edge", values_to="x") %>%
      arrange(habitat_label, x) %>%
      transmute(habitat_label, x, y = count)
    
    # -------- Mise en page --------
    strip_ratio <- 0.12; gap_ratio <- 0.05
    
    # >>> Y FIXE = 11000 & pas de tick 10000 <<<
    y_top        <- 11000
    strip_h_abs  <- y_top * strip_ratio
    strip_ymin   <- y_top - strip_h_abs
    y_breaks     <- pretty(c(0, y_top), n = 6)
    y_breaks     <- y_breaks[y_breaks < y_top & y_breaks != 10000]  # ⬅ enlève 10000
    # <<< ------------------------------------ >>>
    
    panel_breaks <- function(xmin, xmax, side=c("left","middle","right")){
      side <- match.arg(side)
      b <- x_breaks_all[x_breaks_all >= xmin & x_breaks_all <= xmax]
      if (side=="left")   b <- b[b <  xmax]
      if (side=="middle") b <- b[b >  xmin & b < xmax]
      if (side=="right")  b <- b[b >  xmin]
      b
    }
    
    make_panel <- function(xmin, xmax, titre, side, band_fill, show_y=FALSE){
      ggplot() +
        geom_col(data=bins_total, aes(x=xmid, y=count),
                 width=binw, fill="grey88", colour="grey70", linewidth=.25) +
        geom_col(data=bins_hab, aes(x=xmid, y=count, fill=habitat_label),
                 width=binw, position="identity", alpha=.30, colour=NA) +
        geom_path(data=hab_step,
                  aes(x=x, y=y, colour=habitat_label, group=habitat_label),
                  linewidth=1.1, lineend="butt") +
        coord_cartesian(xlim=c(xmin,xmax), ylim=c(0, y_top), clip="on") +
        annotate("rect", xmin=xmin, xmax=xmax, ymin=strip_ymin, ymax=Inf,
                 fill=band_fill, alpha=band_alpha, colour=NA) +
        annotate("segment", x=xmin, xend=xmax, y=strip_ymin, yend=strip_ymin,
                 colour="black", linewidth=1.1) +
        annotate("text", x=(xmin+xmax)/2, y=(strip_ymin+y_top)/2,
                 label=titre, fontface="bold", size=5.3, colour="#222222") +
        scale_x_continuous(breaks=panel_breaks(xmin,xmax,side),
                           limits=c(x_min_force,x_max_force),
                           minor_breaks=NULL, expand=expansion(add=0)) +
        scale_y_continuous(breaks = y_breaks,
                           labels = scales::label_number(big.mark = " "),
                           limits = c(0, y_top),
                           expand = expansion(mult=c(0,0))) +
        scale_fill_manual(values=cols_focus, name="Habitat",
                          drop=FALSE, na.translate=FALSE) +
        scale_colour_manual(values=cols_focus, guide="none", drop=FALSE) +
        labs(x=NULL, y=if (show_y) "Effectif" else NULL) +
        theme_minimal(base_family="Arial", base_size=16) +
        theme(
          panel.border       = element_rect(colour="black", fill=NA, linewidth=1.6),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x        = element_text(size=14),
          axis.text.y        = if (show_y) element_text(size=14) else element_blank(),
          axis.title.y       = if (show_y) element_text(size=15, face="bold") else element_blank(),
          plot.margin        = margin(0,0,0,0,"pt"),
          legend.position    = "right",
          legend.title       = element_text(face="bold", size=14),
          legend.text        = element_text(size=12)
        )
    }
    
    # -------- Panneaux --------
    p1 <- make_panel(bornes[1], bornes[2], "Précoce", side="left",   band_fill=block_cols$precoce, show_y=TRUE)
    p2 <- make_panel(bornes[2], bornes[3], "Moyen",   side="middle", band_fill=block_cols$moyen)
    p3 <- make_panel(bornes[3], bornes[4], "Tardif",  side="right",  band_fill=block_cols$tardif)
    
    legend_col   <- guide_area(); legend_width <- 0.22
    p_row <- (p1 | p2 | p3 | legend_col) +
      plot_layout(widths = c(w / sum(w) * (1 - legend_width), legend_width),
                  guides = "collect") &
      theme(legend.position = "right",
            legend.key.height = unit(12,"pt"),
            legend.key.width  = unit(18,"pt"))
    
    x_title <- ggplot() + theme_void() +
      annotate("text", x=.5, y=.5, label=paste0("SMOD ", year_to_plot),
               fontface="bold", size=6, colour="#222") +
      theme(plot.margin = margin(8,0,0,0))
    
    final_plot <- (p_row / x_title) + plot_layout(heights=c(1,0.08)) +
      plot_annotation(
        title = paste0("Distribution du SMOD (", year_to_plot, ") — histogramme global + habitats d’intérêt"),
        theme = theme(plot.title = element_text(hjust=.5, face="bold", size=18))
      )
    
    final_plot
    
    
    
    
    ggsave(
      filename = file.path(out_dir, paste0("SMOD_", year_to_plot, "_HAB_segments.png")),
      plot     = final_plot,
      width    = 12, height = 8, dpi = 300, bg = "white"
    )
    
    
    
    
    
    
    
    
    
    
    
    
    # === V3 (readable axes) + FIX Y=11000 constant & remove 10000 tick ===
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(grid)
    
    # -------- Global parameters --------
    year_to_plot <- 2023
    x_min_force  <- 15
    x_max_force  <- 190
    n_bins       <- 25
    brks         <- seq(x_min_force, x_max_force, length.out = n_bins + 1)
    binw         <- diff(brks)[1]
    x_breaks_all <- seq(20, 190, by = 10)
    
    # Habitats (order & colors)
    hab_labels <- c("1"="P. nivales","9"="P. thermiques écorchées",
                    "5"="Queyrellins","6"="Nardaies denses du subalpin")
    cols_focus <- c("P. nivales"="#3B5BDB","P. thermiques écorchées"="#E03131",
                    "Queyrellins"="#F4A261","Nardaies denses du subalpin"="#1B9E77")
    
    # Bands (slight transparency)
    block_cols <- list(precoce="#fed976", moyen="#41b6c4", tardif="#0c2c84")
    band_alpha <- 0.4
    
    # -------- Data --------
    dt_plot <- dt_legacy %>%
      filter(YEAR == year_to_plot) %>%
      mutate(habitat_code  = as.character(habitat_code),
             habitat_label = recode(habitat_code, !!!hab_labels, .default = NA_character_),
             habitat_label = factor(habitat_label, levels = names(cols_focus))) %>%
      filter(SMOD >= x_min_force, SMOD <= x_max_force)
    
    # (optional) IQR filter
    q <- quantile(dt_plot$SMOD, c(.25,.75), na.rm=TRUE); iqr <- diff(q)
    dt_plot <- dt_plot %>% filter(between(SMOD, q[1]-1.5*iqr, q[2]+1.5*iqr))
    
    # -------- Bins & counts --------
    edges <- tibble::tibble(
      bin_id = seq_len(length(brks)-1),
      xmin   = brks[-length(brks)],
      xmax   = brks[-1],
      xmid   = (brks[-length(brks)] + brks[-1]) / 2
    )
    bins_total <- dt_plot %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed=TRUE, all.inside=TRUE)) %>%
      count(bin_id, name="count") %>%
      right_join(edges, by="bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L))
    
    bins_hab <- dt_plot %>%
      filter(!is.na(habitat_label)) %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed=TRUE, all.inside=TRUE)) %>%
      count(habitat_label, bin_id, name="count") %>%
      right_join(edges, by="bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L)) %>%
      arrange(habitat_label, xmin)
    
    hab_step <- bins_hab %>%
      dplyr::select(habitat_label, xmin, xmax, count) %>%
      pivot_longer(c(xmin, xmax), names_to="edge", values_to="x") %>%
      arrange(habitat_label, x) %>%
      transmute(habitat_label, x, y = count)
    
    # -------- Layout --------
    strip_ratio <- 0.12; gap_ratio <- 0.05
    
    # >>> FIXED Y = 11000 & remove tick 10000 <<<
    y_top        <- 11000
    strip_h_abs  <- y_top * strip_ratio
    strip_ymin   <- y_top - strip_h_abs
    y_breaks     <- pretty(c(0, y_top), n = 6)
    y_breaks     <- y_breaks[y_breaks < y_top & y_breaks != 10000]  # remove 10000
    # <<< ------------------------------------ >>>
    
    panel_breaks <- function(xmin, xmax, side=c("left","middle","right")){
      side <- match.arg(side)
      b <- x_breaks_all[x_breaks_all >= xmin & x_breaks_all <= xmax]
      if (side=="left")   b <- b[b <  xmax]
      if (side=="middle") b <- b[b >  xmin & b < xmax]
      if (side=="right")  b <- b[b >  xmin]
      b
    }
    
    make_panel <- function(xmin, xmax, title_txt, side, band_fill, show_y=FALSE){
      ggplot() +
        geom_col(data=bins_total, aes(x=xmid, y=count),
                 width=binw, fill="grey88", colour="grey70", linewidth=.25) +
        geom_col(data=bins_hab, aes(x=xmid, y=count, fill=habitat_label),
                 width=binw, position="identity", alpha=.30, colour=NA) +
        geom_path(data=hab_step,
                  aes(x=x, y=y, colour=habitat_label, group=habitat_label),
                  linewidth=1.1, lineend="butt") +
        coord_cartesian(xlim=c(xmin,xmax), ylim=c(0, y_top), clip="on") +
        annotate("rect", xmin=xmin, xmax=xmax, ymin=strip_ymin, ymax=Inf,
                 fill=band_fill, alpha=band_alpha, colour=NA) +
        annotate("segment", x=xmin, xend=xmax, y=strip_ymin, yend=strip_ymin,
                 colour="black", linewidth=1.1) +
        annotate("text", x=(xmin+xmax)/2, y=(strip_ymin+y_top)/2,
                 label=title_txt, fontface="bold", size=5.3, colour="#222222") +
        scale_x_continuous(breaks=panel_breaks(xmin,xmax,side),
                           limits=c(x_min_force,x_max_force),
                           minor_breaks=NULL, expand=expansion(add=0)) +
        scale_y_continuous(breaks = y_breaks,
                           labels = scales::label_number(big.mark = " "),
                           limits = c(0, y_top),
                           expand = expansion(mult=c(0,0))) +
        scale_fill_manual(values=cols_focus, name="Habitat",
                          drop=FALSE, na.translate=FALSE) +
        scale_colour_manual(values=cols_focus, guide="none", drop=FALSE) +
        labs(x=NULL, y=if (show_y) "Count" else NULL) +
        theme_minimal(base_family="Arial", base_size=16) +
        theme(
          panel.border       = element_rect(colour="black", fill=NA, linewidth=1.6),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x        = element_text(size=14),
          axis.text.y        = if (show_y) element_text(size=14) else element_blank(),
          axis.title.y       = if (show_y) element_text(size=15, face="plain") else element_blank(),  # not bold
          plot.margin        = margin(0,0,0,0,"pt"),
          legend.position    = "right",
          legend.title       = element_text(face="bold", size=14),
          legend.text        = element_text(size=12)
        )
    }
    
    # -------- Panels --------
    p1 <- make_panel(bornes[1], bornes[2], "Early", side="left",   band_fill=block_cols$precoce, show_y=TRUE)
    p2 <- make_panel(bornes[2], bornes[3], "Mid",   side="middle", band_fill=block_cols$moyen)
    p3 <- make_panel(bornes[3], bornes[4], "Late",  side="right",  band_fill=block_cols$tardif)
    
    legend_col   <- guide_area(); legend_width <- 0.22
    p_row <- (p1 | p2 | p3 | legend_col) +
      plot_layout(widths = c(w / sum(w) * (1 - legend_width), legend_width),
                  guides = "collect") &
      theme(legend.position = "right",
            legend.key.height = unit(12,"pt"),
            legend.key.width  = unit(18,"pt"))
    
    x_title <- ggplot() + theme_void() +
      annotate("text", x=.5, y=.5, label=paste0("SMOD ", year_to_plot),
               fontface="plain", size=6, colour="#222") +
      theme(plot.margin = margin(8,0,0,0))
    
    final_plot <- (p_row / x_title) + plot_layout(heights=c(1,0.08)) +
      plot_annotation(
        title = paste0("SMOD distribution (", year_to_plot, ") — global histogram + focus habitats"),
        theme = theme(plot.title = element_text(hjust=.5, face="bold", size=18))
      )
    
    final_plot
    
       # --- VERSION EXPORT SANS LÉGENDE (pour montage avec la carte) ----------------
    p_row_no_leg <- (p1 | p2 | p3) +                      # pas de colonne "guide_area"
      plot_layout(widths = w, guides = "keep") &          # on garde les guides internes
      theme(legend.position = "none")                     # ... mais on n'affiche rien
    
    final_plot_nolegend <- (p_row_no_leg / x_title) +
      plot_layout(heights = c(1, 0.06)) +
      plot_annotation(
        title = paste0("Distribution du SMOD (", year_to_plot,
                       ") — histogramme global + habitats d’intérêt"),
        theme = theme(
          plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
          plot.margin  = margin(4, 6, 4, 6)              # marges serrées pour gagner de la place
        )
      )
    
    # Variante encore plus compacte (SANS titre global) si tu veux maximiser la zone utile :
    final_plot_nolegend_notitle <- (p_row_no_leg / x_title) +
      plot_layout(heights = c(1, 0.06))
    
    
    # -----------------------------------------------------------------------------
    # si besoin : install.packages("ragg")
    p_out <- final_plot_nolegend + theme(plot.margin = margin(10,12,10,12))
    
    ggsave(
      filename = file.path(out_dir, sprintf("SMOD_%s_HAB_segments_noleg.png", year_to_plot)),
      plot     = p_out,
      width    = 230, height = 160, units = "mm",
      dpi      = 400,
      bg       = "white",
      device   = ragg::agg_png,
      limitsize = FALSE
    )
    # si besoin : install.packages("ragg")
    p_out2 <- final_plot_nolegend_notitle + theme(plot.margin = margin(8,10,8,10))
    
    ggsave(
      filename = file.path(out_dir, sprintf("SMOD_%s_HAB_segments_noleg_notitle.png", year_to_plot)),
      plot     = p_out2,
      width    = 215, height = 150, units = "mm",
      dpi      = 420,
      bg       = "white",
      device   = ragg::agg_png,
      limitsize = FALSE
    )
    
    # PDF vectoriel (si tu veux coller dans Inkscape/Illustrator/LaTeX)
    ggsave(
      filename = file.path(out_dir, sprintf("SMOD_%s_HAB_segments_noleg.pdf", year_to_plot)),
      plot     = final_plot_nolegend,
      width    = 230, height = 160, units = "mm",
      device   = cairo_pdf, bg = "white", useDingbats = FALSE
    )
    
    
    
    
    
    