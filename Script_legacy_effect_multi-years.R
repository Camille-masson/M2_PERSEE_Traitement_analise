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


source("config.R")

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
    
    dt_legacy_hab_V4 <- dt_legacy
    
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

    
    
  }
    
    
    #### 3.2 Figure SMOD + Habitat M2 ###
    if (TRUE){
  
    
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(grid)
    
    # -------- Paramètres globaux (à adapter si besoin) --------
    year_to_plot <- 2022
    
    # Domaine X + bins communs
    x_min_force  <- 15
    x_max_force  <- 190
    n_bins       <- 25
    brks         <- seq(x_min_force, x_max_force, length.out = n_bins + 1)
    binw         <- diff(brks)[1]
    x_breaks_all <- seq(20, 190, by = 10)  # graduations X pour les 3 panneaux
    
    # Habitats (ordre + couleurs)
    hab_labels <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    cols_focus <- c(
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    
    # Bandes (légère transparence)
    block_cols <- list(precoce = "#fed976", moyen = "#41b6c4", tardif = "#0c2c84")
    band_alpha <- 0.4
    
    # -------- Données --------
    dt_plot <- dt_legacy %>%
      filter(YEAR == year_to_plot) %>%
      mutate(
        habitat_code  = as.character(habitat_code),
        habitat_label = dplyr::recode(habitat_code, !!!hab_labels, .default = NA_character_),
        habitat_label = factor(habitat_label, levels = names(cols_focus))
      ) %>%
      filter(SMOD >= x_min_force, SMOD <= x_max_force)
    
    # (optionnel) filtre IQR sur SMOD
    q   <- stats::quantile(dt_plot$SMOD, c(.25, .75), na.rm = TRUE)
    iqr <- diff(q)
    dt_plot <- dt_plot %>%
      filter(dplyr::between(SMOD, q[1] - 1.5 * iqr, q[2] + 1.5 * iqr))
    
    # -------- Bins & comptes --------
    edges <- tibble::tibble(
      bin_id = seq_len(length(brks) - 1),
      xmin   = brks[-length(brks)],
      xmax   = brks[-1],
      xmid   = (brks[-length(brks)] + brks[-1]) / 2
    )
    
    bins_total <- dt_plot %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed = TRUE, all.inside = TRUE)) %>%
      count(bin_id, name = "count") %>%
      right_join(edges, by = "bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L))
    
    bins_hab <- dt_plot %>%
      filter(!is.na(habitat_label)) %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed = TRUE, all.inside = TRUE)) %>%
      count(habitat_label, bin_id, name = "count") %>%
      right_join(edges, by = "bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L)) %>%
      arrange(habitat_label, xmin)
    
    # Courbes "step" alignées aux bords de bin
    hab_step <- bins_hab %>%
      dplyr::select(habitat_label, xmin, xmax, count) %>%
      tidyr::pivot_longer(c(xmin, xmax), names_to = "edge", values_to = "x") %>%
      arrange(habitat_label, x) %>%
      transmute(habitat_label, x, y = count)
    
    # -------- Bornes des 3 blocs (terciles) "snappées" sur brks --------
    snap_to_break <- function(v, brks) brks[which.min(abs(brks - v))]
    bornes_brutes <- stats::quantile(dt_plot$SMOD, c(0, 1/3, 2/3, 1), na.rm = TRUE)
    bornes <- sort(unique(sapply(bornes_brutes, snap_to_break, brks = brks)))
    
    # Sécurisation si les snaps créent des doublons : fallback à 3 segments égaux
    if (length(bornes) < 4) {
      idx <- unique(pmax(1, pmin(length(brks),
                                 round(c(1, length(brks)/3, 2*length(brks)/3, length(brks))))))
      bornes <- sort(unique(brks[idx]))
    }
    # À ce stade, bornes = c(xmin, cut1, cut2, xmax)
    
    # -------- Mise en page / axes --------
    strip_ratio <- 0.12
    gap_ratio   <- 0.05
    
    # Y fixé à 11000 + suppression de la graduation 10000
    y_top       <- 11000
    strip_h_abs <- y_top * strip_ratio
    strip_ymin  <- y_top - strip_h_abs
    y_breaks    <- pretty(c(0, y_top), n = 6)
    y_breaks    <- y_breaks[y_breaks < y_top & y_breaks != 10000]
    
    panel_breaks <- function(xmin, xmax, side = c("left", "middle", "right")) {
      side <- match.arg(side)
      b <- x_breaks_all[x_breaks_all >= xmin & x_breaks_all <= xmax]
      if (side == "left")   b <- b[b <  xmax]
      if (side == "middle") b <- b[b >  xmin & b < xmax]
      if (side == "right")  b <- b[b >  xmin]
      b
    }
    
    make_panel <- function(xmin, xmax, title_txt, side, band_fill, show_y = FALSE) {
      ggplot() +
        geom_col(data = bins_total, aes(x = xmid, y = count),
                 width = binw, fill = "grey88", colour = "grey70", linewidth = .25) +
        geom_col(data = bins_hab, aes(x = xmid, y = count, fill = habitat_label),
                 width = binw, position = "identity", alpha = .30, colour = NA) +
        geom_path(data = hab_step,
                  aes(x = x, y = y, colour = habitat_label, group = habitat_label),
                  linewidth = 1.1, lineend = "butt") +
        coord_cartesian(xlim = c(xmin, xmax), ylim = c(0, y_top), clip = "on") +
        annotate("rect", xmin = xmin, xmax = xmax, ymin = strip_ymin, ymax = Inf,
                 fill = band_fill, alpha = band_alpha, colour = NA) +
        annotate("segment", x = xmin, xend = xmax, y = strip_ymin, yend = strip_ymin,
                 colour = "black", linewidth = 1.1) +
        annotate("text", x = (xmin + xmax) / 2, y = (strip_ymin + y_top) / 2,
                 label = title_txt, fontface = "bold", size = 5.3, colour = "#222222") +
        scale_x_continuous(breaks = panel_breaks(xmin, xmax, side),
                           limits = c(x_min_force, x_max_force),
                           minor_breaks = NULL, expand = expansion(add = 0)) +
        scale_y_continuous(breaks = y_breaks,
                           labels = scales::label_number(big.mark = " "),
                           limits = c(0, y_top),
                           expand = expansion(mult = c(0, 0))) +
        scale_fill_manual(values = cols_focus, name = "Habitat",
                          drop = FALSE, na.translate = FALSE) +
        scale_colour_manual(values = cols_focus, guide = "none", drop = FALSE) +
        labs(x = NULL, y = if (show_y) "Count" else NULL) +
        theme_minimal(base_family = "Arial", base_size = 16) +
        theme(
          panel.border       = element_rect(colour = "black", fill = NA, linewidth = 1.6),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x        = element_text(size = 14),
          axis.text.y        = if (show_y) element_text(size = 14) else element_blank(),
          axis.title.y       = if (show_y) element_text(size = 15, face = "plain") else element_blank(),
          plot.margin        = margin(0, 0, 0, 0, "pt"),
          legend.position    = "right",
          legend.title       = element_text(face = "bold", size = 14),
          legend.text        = element_text(size = 12)
        )
    }
    
    # -------- Panneaux --------
    p1 <- make_panel(bornes[1], bornes[2], "Early", side = "left",   band_fill = block_cols$precoce, show_y = TRUE)
    p2 <- make_panel(bornes[2], bornes[3], "Mid",   side = "middle", band_fill = block_cols$moyen)
    p3 <- make_panel(bornes[3], bornes[4], "Late",  side = "right",  band_fill = block_cols$tardif)
    
    # Largeurs relatives des panneaux proportionnelles à la largeur de chaque bloc
    w <- diff(bornes)                 # 3 valeurs (Early, Mid, Late)
    legend_width <- 0.22
    legend_col   <- guide_area()
    
    p_row <- (p1 | p2 | p3 | legend_col) +
      plot_layout(
        widths = c(w / sum(w) * (1 - legend_width), legend_width),
        guides = "collect"
      ) &
      theme(
        legend.position    = "right",
        legend.key.height  = unit(12, "pt"),
        legend.key.width   = unit(18, "pt")
      )
    
    # Titre de l'axe X
    x_title <- ggplot() + theme_void() +
      annotate("text", x = .5, y = .5, label = paste0("SMOD ", year_to_plot),
               fontface = "plain", size = 6, colour = "#222") +
      theme(plot.margin = margin(8, 0, 0, 0))
    
    # -------- Figure finale --------
    final_plot <- (p_row / x_title) + plot_layout(heights = c(1, 0.08)) +
      plot_annotation(
        title = paste0("SMOD distribution (", year_to_plot, ") — global histogram + focus habitats"),
        theme = theme(plot.title = element_text(hjust = .5, face = "bold", size = 18))
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
    }
    
  # si besoin :
  # install.packages("svglite")
  library(svglite)
  library(patchwork)
  library(ggplot2)
  
  # 1) Row sans légende
  p_row_no_leg <- (p1 | p2 | p3) +
    plot_layout(widths = w, guides = "keep") &     # structure inchangée
    theme(legend.position = "none")                 # légendes retirées partout
  
  # 2) Assemblage SANS titre global
  final_plot_nolegend_notitle <- (p_row_no_leg / x_title) +
    plot_layout(heights = c(1, 0.06)) &
    theme(plot.title = element_blank())             # au cas où un titre subsiste
  
  # 3) Export SVG (fond transparent pour montage carte ; mets "white" si besoin)
  p_out_svg <- final_plot_nolegend_notitle +
    theme(plot.margin = margin(8, 10, 8, 10))
  
  ggsave(
    filename = file.path(out_dir, sprintf("SMOD_%s_HAB_segments_noleg_notitle.svg", year_to_plot)),
    plot     = p_out_svg,
    width    = 215, height = 150, units = "mm",
    device   = svglite::svglite,
    bg       = "transparent",
    limitsize = FALSE
  )
  
  
  
  
  
  
  
  
  
  #### 3.2 Figure SMOD + Habitat M2 — Classes fixes ####
  if (TRUE){
    
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(grid)
    
    # -------- Paramètres globaux (à adapter si besoin) --------
    year_to_plot <- 2023
    
    # Domaine X + bins communs
    x_min_force  <- 15
    x_max_force  <- 190
    n_bins       <- 25
    brks         <- seq(x_min_force, x_max_force, length.out = n_bins + 1)
    binw         <- diff(brks)[1]
    x_breaks_all <- seq(20, 190, by = 10)  # graduations X pour les 3 panneaux
    
    # Habitats (ordre + couleurs)
    hab_labels <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    cols_focus <- c(
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    
    # Bandes (légère transparence)
    block_cols <- list(precoce = "#fed976", moyen = "#41b6c4", tardif = "#0c2c84")
    band_alpha <- 0.4
    
    # -------- Données --------
    dt_plot <- dt_legacy %>%
      filter(YEAR == year_to_plot) %>%
      mutate(
        habitat_code  = as.character(habitat_code),
        habitat_label = dplyr::recode(habitat_code, !!!hab_labels, .default = NA_character_),
        habitat_label = factor(habitat_label, levels = names(cols_focus))
      ) %>%
      filter(SMOD >= x_min_force, SMOD <= x_max_force)
    
    # (optionnel) filtre IQR sur SMOD
    q   <- stats::quantile(dt_plot$SMOD, c(.25, .75), na.rm = TRUE)
    iqr <- diff(q)
    dt_plot <- dt_plot %>%
      filter(dplyr::between(SMOD, q[1] - 1.5 * iqr, q[2] + 1.5 * iqr))
    
    # -------- Bins & comptes --------
    edges <- tibble::tibble(
      bin_id = seq_len(length(brks) - 1),
      xmin   = brks[-length(brks)],
      xmax   = brks[-1],
      xmid   = (brks[-length(brks)] + brks[-1]) / 2
    )
    
    bins_total <- dt_plot %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed = TRUE, all.inside = TRUE)) %>%
      count(bin_id, name = "count") %>%
      right_join(edges, by = "bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L))
    
    bins_hab <- dt_plot %>%
      filter(!is.na(habitat_label)) %>%
      mutate(bin_id = findInterval(SMOD, brks, rightmost.closed = TRUE, all.inside = TRUE)) %>%
      count(habitat_label, bin_id, name = "count") %>%
      right_join(edges, by = "bin_id") %>%
      mutate(count = tidyr::replace_na(count, 0L)) %>%
      arrange(habitat_label, xmin)
    
    # Courbes "step" alignées aux bords de bin
    hab_step <- bins_hab %>%
      dplyr::select(habitat_label, xmin, xmax, count) %>%
      tidyr::pivot_longer(c(xmin, xmax), names_to = "edge", values_to = "x") %>%
      arrange(habitat_label, x) %>%
      transmute(habitat_label, x, y = count)
    
    # -------- Bornes FIXES des 3 blocs (plus de terciles) --------
    # Bornes visuelles des panneaux aux seuils 120 et 150
    bornes <- c(x_min_force, 120, 150, x_max_force)
    
    # Variable de classe stricte (conforme à la demande) :
    # Early : SMOD < 120
    # Mid   : 121 ≤ SMOD ≤ 150
    # Late  : SMOD > 151
    # NB : 120 et 151 ne sont rattachés à aucune classe.
    dt_plot <- dt_plot %>%
      mutate(phase_fix = dplyr::case_when(
        SMOD < 120                ~ "Early",
        SMOD >= 121 & SMOD <= 150 ~ "Mid",
        SMOD > 151                ~ "Late",
        TRUE ~ NA_character_
      ))
    
    # -------- Mise en page / axes --------
    strip_ratio <- 0.12
    gap_ratio   <- 0.05
    
    # Y fixé à 11000 + suppression de la graduation 10000
    y_top       <- 11000
    strip_h_abs <- y_top * strip_ratio
    strip_ymin  <- y_top - strip_h_abs
    y_breaks    <- pretty(c(0, y_top), n = 6)
    y_breaks    <- y_breaks[y_breaks < y_top & y_breaks != 10000]
    
    panel_breaks <- function(xmin, xmax, side = c("left", "middle", "right")) {
      side <- match.arg(side)
      b <- x_breaks_all[x_breaks_all >= xmin & x_breaks_all <= xmax]
      if (side == "left")   b <- b[b <  xmax]
      if (side == "middle") b <- b[b >  xmin & b < xmax]
      if (side == "right")  b <- b[b >  xmin]
      b
    }
    
    make_panel <- function(xmin, xmax, title_txt, side, band_fill, show_y = FALSE) {
      ggplot() +
        geom_col(data = bins_total, aes(x = xmid, y = count),
                 width = binw, fill = "grey88", colour = "grey70", linewidth = .25) +
        geom_col(data = bins_hab, aes(x = xmid, y = count, fill = habitat_label),
                 width = binw, position = "identity", alpha = .30, colour = NA) +
        geom_path(data = hab_step,
                  aes(x = x, y = y, colour = habitat_label, group = habitat_label),
                  linewidth = 1.1, lineend = "butt") +
        coord_cartesian(xlim = c(xmin, xmax), ylim = c(0, y_top), clip = "on") +
        annotate("rect", xmin = xmin, xmax = xmax, ymin = strip_ymin, ymax = Inf,
                 fill = band_fill, alpha = band_alpha, colour = NA) +
        annotate("segment", x = xmin, xend = xmax, y = strip_ymin, yend = strip_ymin,
                 colour = "black", linewidth = 1.1) +
        annotate("text", x = (xmin + xmax) / 2, y = (strip_ymin + y_top) / 2,
                 label = title_txt, fontface = "bold", size = 5.3, colour = "#222222") +
        scale_x_continuous(breaks = panel_breaks(xmin, xmax, side),
                           limits = c(x_min_force, x_max_force),
                           minor_breaks = NULL, expand = expansion(add = 0)) +
        scale_y_continuous(breaks = y_breaks,
                           labels = scales::label_number(big.mark = " "),
                           limits = c(0, y_top),
                           expand = expansion(mult = c(0, 0))) +
        scale_fill_manual(values = cols_focus, name = "Habitat",
                          drop = FALSE, na.translate = FALSE) +
        scale_colour_manual(values = cols_focus, guide = "none", drop = FALSE) +
        labs(x = NULL, y = if (show_y) "Count" else NULL) +
        theme_minimal(base_family = "Arial", base_size = 16) +
        theme(
          panel.border       = element_rect(colour = "black", fill = NA, linewidth = 1.6),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x        = element_text(size = 14),
          axis.text.y        = if (show_y) element_text(size = 14) else element_blank(),
          axis.title.y       = if (show_y) element_text(size = 15, face = "plain") else element_blank(),
          plot.margin        = margin(0, 0, 0, 0, "pt"),
          legend.position    = "right",
          legend.title       = element_text(face = "bold", size = 14),
          legend.text        = element_text(size = 12)
        )
    }
    
    # -------- Panneaux --------
    p1 <- make_panel(bornes[1], bornes[2], "Early ",    side = "left",
                     band_fill = block_cols$precoce, show_y = TRUE)
    p2 <- make_panel(bornes[2], bornes[3], "Mid ",   side = "middle",
                     band_fill = block_cols$moyen)
    p3 <- make_panel(bornes[3], bornes[4], "Late ",     side = "right",
                     band_fill = block_cols$tardif)
    
    # Largeurs relatives des panneaux proportionnelles à la largeur de chaque bloc
    w <- diff(bornes)                 # 3 valeurs (Early, Mid, Late)
    legend_width <- 0.22
    legend_col   <- guide_area()
    
    p_row <- (p1 | p2 | p3 | legend_col) +
      plot_layout(
        widths = c(w / sum(w) * (1 - legend_width), legend_width),
        guides = "collect"
      ) &
      theme(
        legend.position    = "right",
        legend.key.height  = unit(12, "pt"),
        legend.key.width   = unit(18, "pt")
      )
    
    # Titre de l'axe X
    x_title <- ggplot() + theme_void() +
      annotate("text", x = .5, y = .5, label = paste0("SMOD ", year_to_plot),
               fontface = "plain", size = 6, colour = "#222") +
      theme(plot.margin = margin(8, 0, 0, 0))
    
    # -------- Figure finale --------
    final_plot <- (p_row / x_title) + plot_layout(heights = c(1, 0.08)) +
      plot_annotation(
        title = paste0("Distribution du SMOD (", year_to_plot, ") — histogramme global + habitats d’intérêt"),
        theme = theme(plot.title = element_text(hjust = .5, face = "bold", size = 18))
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
  }
  
  # si besoin :
  # install.packages("svglite")
  library(svglite)
  library(patchwork)
  library(ggplot2)
  
  # 1) Row sans légende
  p_row_no_leg <- (p1 | p2 | p3) +
    plot_layout(widths = w, guides = "keep") &     # structure inchangée
    theme(legend.position = "none")                 # légendes retirées partout
  
  # 2) Assemblage SANS titre global
  final_plot_nolegend_notitle <- (p_row_no_leg / x_title) +
    plot_layout(heights = c(1, 0.06)) &
    theme(plot.title = element_blank())             # au cas où un titre subsiste
  
  # 3) Export SVG (fond transparent pour montage carte ; mets "white" si besoin)
  p_out_svg <- final_plot_nolegend_notitle +
    theme(plot.margin = margin(8, 10, 8, 10))
  
  ggsave(
    filename = file.path(out_dir, sprintf("SMOD_%s_HAB_segments_noleg_notitle.svg", year_to_plot)),
    plot     = p_out_svg,
    width    = 215, height = 150, units = "mm",
    device   = svglite::svglite,
    bg       = "transparent",
    limitsize = FALSE
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #### 3.3 Analyse multi-year
  
  
  ## PARAMETRE
  alpage = "Alpe-Sud"
  ²YEARs = "2017-2023"
  
  # LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_legacy.R"))
  library(ggplot2)
  library(dplyr)
  
  ## ENTREE
  out_dir <- file.path(output_dir, "10. NDVI Effect",
                       paste0("dataset_", alpage))
  
  dt_filtered = file.path(out_dir, paste0("df_filtered_legacy_multi-year_", YEARs, "_", alpage, ".rds"))
  
  dt_legacy = readRDS(dt_filtered)
  
  
  summary(dt_legacy)
  str(dt_legacy) 
  
  
  
  
  
  library(mgcv)
  
  m00 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts"),
    data = dt_legacy, method = "fREML", discrete = TRUE, select = TRUE
  )
  
  
  gam.check(m00); summary(m00)$dev.expl
  
  summary(m00)
  
  
  m01 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")),
    data = dt_legacy, method = "fREML", discrete = TRUE, select = TRUE
  )
  
  summary(m01)
  AIC(m00 , m01)
  anova(m00, m01, test="Chisq")
  
  
  
  m02 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40),   # ou s(x_km, y_km, ...)
    data = dt_legacy, method="fREML", discrete=TRUE, select=TRUE
  )
  
  summary(m02)
  
  AIC(m00 , m01, m02)
  anova(m00, m01, m02 ,test="Chisq")
  gam.check(m02)
  
  
  
  dt_legacy$charge_bin <- factor(dt_legacy$charge_bin)
  dt_legacy$alpage <- factor(dt_legacy$alpage)
  
  m03 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      charge_bin +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40),
    data = dt_legacy, method="fREML", discrete=TRUE, select=FALSE
  )
  
  summary(m03)
  
  AIC(m00 , m01, m02, m03)
  
  
  # 3) Bon pattern "factor-by" (on met aussi le lisse de base)
  #    -> courbe de base + déviations par niveau :
  m04 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      charge_bin +                                  # intercepts par niveau
      s(Charge_log, bs="ts") +                      # lisse de base
      s(Charge_log, by=charge_bin, bs="ts") +       # déviations par niveau
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40),
    data = dt_legacy, method="fREML", discrete=TRUE, select=FALSE
  )
  
  summary(m04)
  
  AIC(m00 , m01, m02, m03, m04)
  

  
  
  m05 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(Charge_log, bs="ts") +                      # lisse de base
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40),
    data = dt_legacy, method="fREML", discrete=TRUE, select=FALSE
  )
  
  summary(m05)
  
  AIC(m00 , m01, m02, m03, m04, m05)
  
  
  
  dt_legacy$year_f <- factor(dt_legacy$YEAR)
  dt_legacy$alpage <- factor(dt_legacy$alpage)
  m06 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(Charge_log, bs="ts") +                      # lisse de base
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40)+
      year_f,
    data = dt_legacy, method="fREML", discrete=TRUE, select=FALSE
  )
  
  summary(m06)
  AIC(m00 , m01, m02, m03, m04, m05, m06)
  
  
  
  
  m07 <- bam(
    log(GPROD) ~ 
      s(DAH, k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(Charge_log, by=year_f, bs="ts") +   # effet Charge_log différent selon l’année
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      year_f,
    data = dt_legacy, method="fREML", discrete=TRUE
  )
  
  
  
  summary(m07)
  AIC(m00 , m01, m02, m03, m04, m05, m06, m07)
  
  
  
  
  
  m08 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(Charge_log, bs="ts") +                             # effet principal commun
      ti(Charge_log, year_f, bs=c("ts","re"), k=c(10, NA)) # déviation par année
    + s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      year_f,
    data = dt_legacy, method="fREML", discrete=TRUE
  )
  
  
  
  summary(m08)
  AIC(m00 , m01, m02, m03, m04, m05, m06, m07, m08)
  
  
  
  m07<- bam(
    log(GPROD) ~ 
      s(DAH, k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(Charge_log, by=year_f, bs="ts") +   # effet Charge_log différent selon l’année
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      year_f,
    data = dt_legacy, method="fREML", discrete=TRUE
  )
  
  
  
  library(dplyr)
  
  # 1) Récap des quantiles de Charge_log et covariables de ref par année
  ref_by_year <- dt_legacy %>%
    group_by(year_f) %>%
    summarise(
      c_lo = quantile(Charge_log, 0.10, na.rm = TRUE),
      c_hi = quantile(Charge_log, 0.60, na.rm = TRUE),
      DAH  = median(DAH,  na.rm = TRUE),
      SMOD = median(SMOD, na.rm = TRUE),
      x    = median(x,    na.rm = TRUE),
      y    = median(y,    na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2) Fixer l’alpage (un niveau existant). Tu peux mettre un alpage "référence".
  alp_ref <- levels(dt_legacy$alpage)[1]
  
  
  mk_newdat <- function(df, charge_value) {
    transform(df,
              Charge_log = charge_value,
              alpage = alp_ref)
  }
  
  new_lo <- mk_newdat(ref_by_year, ref_by_year$c_lo)
  new_hi <- mk_newdat(ref_by_year, ref_by_year$c_hi)
  
  
  library(mgcv)
  
  # Matrices de design
  X_lo <- predict(m07, newdata = new_lo, type = "lpmatrix")
  X_hi <- predict(m07, newdata = new_hi, type = "lpmatrix")
  
  # Contraste (une ligne par année)
  X_diff <- X_hi - X_lo
  
  # Coefficients et matrice de variance-covariance
  beta  <- coef(m07)
  Vb    <- vcov(m07)
  
  # Δ sur log(GPROD)
  delta <- as.vector(X_diff %*% beta)
  
  # SE(Δ)
  se_delta <- sqrt(rowSums((X_diff %*% Vb) * X_diff))
  
  # IC 95% sur log-échelle
  lo_delta <- delta - 1.96 * se_delta
  hi_delta <- delta + 1.96 * se_delta
  
  # Passage à l’échelle GPROD (ratio multiplicatif)
  ratio    <- exp(delta)
  lo_ratio <- exp(lo_delta)
  hi_ratio <- exp(hi_delta)
  
  effect_year <- ref_by_year %>%
    mutate(
      delta_log = delta,
      delta_se  = se_delta,
      delta_lo  = lo_delta,
      delta_hi  = hi_delta,
      ratio     = ratio,
      ratio_lo  = lo_ratio,
      ratio_hi  = hi_ratio
    ) %>%
    dplyr::select(year_f, c_lo, c_hi, delta_log, delta_se, delta_lo, delta_hi, ratio, ratio_lo, ratio_hi)
  
  
  library(ggplot2)
  
  ggplot(effect_year, aes(x = year_f, y = ratio)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ratio_lo, ymax = ratio_hi), width = 0.15) +
    geom_hline(yintercept = 1, linetype = 2) +
    labs(x = "Année", y = "Taille d'effet du pâturage (ratio Q90 vs Q10 de Charge_log)",
         title = "Effet du pâturage sur GPROD par année",
         subtitle = "Ratio = exp(Δ log(GPROD)) entre charge élevée (Q90) et faible (Q10)") +
    theme_minimal()
  
  
  
  
  
  
  
  # --- Pré-requis
  library(dplyr)
  library(mgcv)
  library(purrr)
  library(ggplot2)
  
  # 0) Quantiles globaux (comparables entre années)
  q <- quantile(dt_legacy$Charge_log, probs = c(0.10, 0.90), na.rm = TRUE)
  q_lo <- unname(q[1]); q_hi <- unname(q[2])
  
  # 1) Distribution de référence commune (intégration sur covariables)
  set.seed(1)
  base_design <- dt_legacy %>%
    select(DAH, SMOD, x, y) %>%
    sample_n(size = min(5000, n()), replace = FALSE)
  
  lev_year <- levels(dt_legacy$year_f)
  alp_ref  <- levels(dt_legacy$alpage)[1]
  
  # 2) Extraire beta & Vb, en tenant compte de l'incertitude de lissage
  beta <- coef(m07)
  Vb   <- vcov(m07, unconditional = TRUE)
  
  # 3) On exclut les effets alpage & spatial pour un effet population-level
  term_excl <- c("s(alpage)", "s(x,y)")
  
  compute_year_contrast <- function(year_level, q_lo, q_hi, use_within_year = FALSE) {
    base <- if (use_within_year) {
      dt_legacy %>% filter(year_f == year_level) %>% select(DAH, SMOD, x, y)
    } else {
      base_design
    }
    if (nrow(base) == 0) return(NULL)
    
    make_new <- function(cval) {
      mutate(base,
             Charge_log = cval,
             year_f     = factor(year_level, levels = lev_year),
             alpage     = alp_ref
      )
    }
    new_lo <- make_new(q_lo)
    new_hi <- make_new(q_hi)
    
    X_lo <- predict(m07, newdata = new_lo, type = "lpmatrix", exclude = term_excl)
    X_hi <- predict(m07, newdata = new_hi, type = "lpmatrix", exclude = term_excl)
    
    Xbar <- colMeans(X_hi - X_lo)                     # moyenne du contraste
    delta <- as.numeric(Xbar %*% beta)                # Δ log(GPROD)
    se    <- sqrt(sum((Xbar %*% Vb) * Xbar))          # SE(Δ)
    
    tibble(
      year_f   = year_level,
      delta_log = delta,
      delta_se  = se,
      delta_lo  = delta - 1.96*se,
      delta_hi  = delta + 1.96*se,
      ratio     = exp(delta),
      ratio_lo  = exp(delta - 1.96*se),
      ratio_hi  = exp(delta + 1.96*se)
    )
  }
  
  # 4) Effet marginal standardisé (comparables entre années)
  effect_year_std <- map_dfr(lev_year, compute_year_contrast, q_lo, q_hi, use_within_year = FALSE)
  
  # (Option) Effet marginal "dans l'année" (intègre sur la distribution propre à chaque année)
  effect_year_inyear <- map_dfr(lev_year, compute_year_contrast, q_lo, q_hi, use_within_year = TRUE)
  
  # 5) Graphique recommandé : contraste standardisé (population-level)
  ggplot(effect_year_std, aes(x = year_f, y = ratio)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ratio_lo, ymax = ratio_hi), width = 0.15) +
    geom_hline(yintercept = 1, linetype = 2) +
    labs(
      x = "Année",
      y = "Taille d'effet de la charge (Q90 vs Q10, ratio sur GPROD)",
      title = "Effet marginal standardisé de la charge pastorale (Charge_log) par année",
      subtitle = paste0("Contraste global ", round(q_lo, 2), " → ", round(q_hi, 2),
                        " ; covariables intégrées ; effets alpage & spatial exclus")
    ) +
    theme_minimal()
  
  
  
  
  
  
  
  ### APPROCHE HABITATS : 
  
  
  
  
  library(dplyr)
  library(mgcv)
  
  ## — mapping des 4 habitats ciblés
  hab_map <- c("1"="P. nivales",
               "9"="P. thermiques écorchées",
               "5"="Queyrellins",
               "6"="Nardaies denses du subalpin")
  hab_codes <- names(hab_map)
  
  ## — restreindre le jeu de données à ces 4 habitats (on peut élargir plus tard si besoin)
  dt_hab <- dt_legacy %>%
    filter(habitat_code %in% hab_codes) %>%
    mutate(
      habitat4 = factor(habitat_code, levels = hab_codes, labels = unname(hab_map)),
      year_f   = factor(year_f)  # au cas où
    )
  
  ## — facteur interaction Année×Habitat (utile si tu veux un terme d'interaction)
  dt_hab <- dt_hab %>%
    mutate(year_hab = interaction(year_f, habitat4, drop = TRUE))
  
  ## — modèle HGAM avec :
  ##    - lissage global de Charge_log
  ##    - déviations par année (fs)
  ##    - déviations par habitat (fs)
  ##    - (optionnel) déviations spécifiques Année×Habitat (fs) si assez de données
  m08 <- bam(
    log(GPROD) ~
      s(DAH,  k = 10, bs = "ts") +
      s(SMOD, k = 10, bs = "ts") +
      ti(DAH, SMOD, k = c(20,20), bs = c("ts","ts")) +
      s(Charge_log, k = 10, bs = "tp") +                 # global
      s(Charge_log, year_f,   bs = "fs", k = 6, m = 1) + # déviation par année
      s(Charge_log, habitat4, bs = "fs", k = 6, m = 1) + # déviation par habitat
      # --- décommente la ligne suivante seulement si tu as assez d'obs dans
      # --- plusieurs combos Année×Habitat ; sinon laisse le pénaliser à zéro :
      # s(Charge_log, year_hab, bs = "fs", k = 5, m = 1) +
      s(alpage, bs = "re") +
      s(x, y, bs = "tp", k = 35) +
      year_f + habitat4,
    data = dt_hab, method = "fREML", discrete = TRUE
  )
  
  
  
  
  
  
  library(purrr)
  library(ggplot2)
  
  ## --- 0) Quantiles globaux (sur les données des 4 habitats)
  q <- quantile(dt_hab$Charge_log, probs = c(0.10, 0.90), na.rm = TRUE)
  q_lo <- unname(q[1]); q_hi <- unname(q[2])
  
  ## --- 1) Distribution de référence commune pour l’intégration
  set.seed(1)
  base_design <- dt_hab %>%
    select(DAH, SMOD, x, y) %>%
    sample_n(size = min(5000, n()), replace = FALSE)
  
  lev_year <- levels(dt_hab$year_f)
  lev_hab  <- levels(dt_hab$habitat4)
  alp_ref  <- levels(dt_hab$alpage)[1]
  
  ## --- 2) Coefs & Vb (incertitude de lissage incluse)
  beta <- coef(m08)
  Vb   <- vcov(m08, unconditional = TRUE)
  
  ## --- 3) On veut un effet “population-level” → exclure alpage & spatial
  term_excl <- c("s(alpage)", "s(x,y)")
  
  ## --- 4) Helpers
  mk_new <- function(df, cval, year_level, hab_level) {
    mutate(df,
           Charge_log = cval,
           year_f     = factor(year_level, levels = lev_year),
           habitat4   = factor(hab_level, levels = lev_hab),
           year_hab   = interaction(year_f, habitat4, drop = TRUE),
           alpage     = alp_ref
    )
  }
  
  ## — contraste moyen pour (année donnée, moyenné sur 4 habitats à parts égales)
  contrast_year_total <- function(year_level) {
    # empile la base pour les 4 habitats, poids égaux
    base_all_hab <- map_dfr(lev_hab, ~ base_design %>% mutate(habitat4 = .x))
    new_lo <- mk_new(base_all_hab, q_lo, year_level, base_all_hab$habitat4)
    new_hi <- mk_new(base_all_hab, q_hi, year_level, base_all_hab$habitat4)
    
    X_lo <- predict(m08, newdata = new_lo, type = "lpmatrix", exclude = term_excl)
    X_hi <- predict(m08, newdata = new_hi, type = "lpmatrix", exclude = term_excl)
    
    # moyenne simple sur toutes les lignes (égale pour chaque habitat)
    Xbar <- colMeans(X_hi - X_lo)
    delta <- as.numeric(Xbar %*% beta)
    se    <- sqrt(sum((Xbar %*% Vb) * Xbar))
    
    tibble::tibble(
      year_f    = year_level,
      delta_log = delta,
      delta_se  = se,
      delta_lo  = delta - 1.96*se,
      delta_hi  = delta + 1.96*se,
      ratio     = exp(delta),
      ratio_lo  = exp(delta - 1.96*se),
      ratio_hi  = exp(delta + 1.96*se)
    )
  }
  
  ## — contraste pour (année, habitat) spécifique
  #    (on calcule seulement sur les combos observés pour éviter des extrapolations)
  combos_obs <- dt_hab %>% count(year_f, habitat4) %>% filter(n > 0) %>%
    transmute(year_f = as.character(year_f), habitat4 = as.character(habitat4))
  
  contrast_year_hab <- function(year_level, hab_level) {
    base_h <- base_design %>% mutate(habitat4 = hab_level)
    new_lo <- mk_new(base_h, q_lo, year_level, hab_level)
    new_hi <- mk_new(base_h, q_hi, year_level, hab_level)
    
    X_lo <- predict(m08, newdata = new_lo, type = "lpmatrix", exclude = term_excl)
    X_hi <- predict(m08, newdata = new_hi, type = "lpmatrix", exclude = term_excl)
    
    Xbar <- colMeans(X_hi - X_lo)
    delta <- as.numeric(Xbar %*% beta)
    se    <- sqrt(sum((Xbar %*% Vb) * Xbar))
    
    tibble::tibble(
      year_f    = year_level,
      habitat4  = hab_level,
      delta_log = delta,
      delta_se  = se,
      delta_lo  = delta - 1.96*se,
      delta_hi  = delta + 1.96*se,
      ratio     = exp(delta),
      ratio_lo  = exp(delta - 1.96*se),
      ratio_hi  = exp(delta + 1.96*se)
    )
  }
  
  ## --- 5) Calculs
  effect_year_total <- purrr::map_dfr(lev_year, contrast_year_total)
  
  effect_year_hab <- purrr::pmap_dfr(
    list(combos_obs$year_f, combos_obs$habitat4),
    contrast_year_hab
  ) %>%
    mutate(
      year_f   = factor(year_f, levels = lev_year),
      habitat4 = factor(habitat4, levels = lev_hab)
    )
  
  
  
  p_total <- ggplot(effect_year_total, aes(x = year_f, y = ratio)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ratio_lo, ymax = ratio_hi), width = 0.15) +
    geom_hline(yintercept = 1, linetype = 2) +
    labs(
      x = "Année",
      y = "Taille d'effet de la charge (Q90 vs Q10, ratio sur GPROD)",
      title = "Effet marginal standardisé de la charge (tous habitats, poids égaux)",
      subtitle = paste0("Contraste global ", round(q_lo, 2), " → ", round(q_hi, 2),
                        " ; covariables intégrées ; effets alpage & spatial exclus")
    ) +
    theme_minimal()
  p_total
  
  
  
  
  
  
  p_by_hab <- ggplot(effect_year_hab, aes(x = year_f, y = ratio)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ratio_lo, ymax = ratio_hi), width = 0.15) +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_wrap(~ habitat4) +
    labs(
      x = "Année",
      y = "Taille d'effet de la charge (Q90 vs Q10, ratio sur GPROD)",
      title = "Effet marginal standardisé par habitat",
      subtitle = paste0("Contraste global ", round(q_lo, 2), " → ", round(q_hi, 2),
                        " ; covariables intégrées ; effets alpage & spatial exclus")
    ) +
    theme_minimal()
  p_by_hab
  
  
  
  
  
  
  
  
  library(dplyr)
  library(mgcv)
  library(purrr)
  library(ggplot2)
  
  # --- mapping des 4 habitats ciblés
  hab_map <- c("1"="P. nivales",
               "9"="P. thermiques écorchées",
               "5"="Queyrellins",
               "6"="Nardaies denses du subalpin")
  hab_codes <- names(hab_map)
  
  # --- sous-données : uniquement ces 4 habitats (et colonnes utiles)
  dt_hab <- dt_legacy %>%
    filter(habitat_code %in% hab_codes) %>%
    mutate(
      habitat4 = factor(habitat_code, levels = hab_codes, labels = unname(hab_map)),
      year_f   = factor(year_f),
      alpage   = factor(alpage)
    ) %>%
    filter(is.finite(GPROD), is.finite(DAH), is.finite(SMOD),
           is.finite(Charge_log), is.finite(x), is.finite(y)) %>%
    droplevels()
  
  # --- k "sûrs" et simples (petits) pour éviter l'erreur k>unicité
  min_u_year <- dt_hab %>% group_by(year_f)   %>% summarise(nu=n_distinct(Charge_log), .groups="drop") %>% pull(nu) %>% min(na.rm=TRUE)
  min_u_hab  <- dt_hab %>% group_by(habitat4) %>% summarise(nu=n_distinct(Charge_log), .groups="drop") %>% pull(nu) %>% min(na.rm=TRUE)
  
  k_year <- max(4, min(8, min_u_year - 1))
  k_hab  <- max(4, min(8, min_u_hab  - 1))
  
  # --- modèle: lissage par année + lissage par habitat (additifs)
  m07_hab <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(Charge_log, by = year_f,   bs="ts", k = k_year) +   # lissage par année
      s(Charge_log, by = habitat4, bs="ts", k = k_hab ) +   # lissage par habitat
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=30) +
      year_f + habitat4,
    data     = dt_hab,
    method   = "fREML",
    select   = TRUE,
    discrete = TRUE
  )
  
  
  
  
  
  
  library(dplyr); library(mgcv)
  
  # facteur Année×Habitat
  dt_hab <- dt_hab %>%
    mutate(year_hab = interaction(year_f, habitat4, drop = TRUE))
  
  # petit k sûr pour l’interaction (évite l’erreur k>unicité)
  min_u_yh <- dt_hab %>% group_by(year_hab) %>% summarise(nu=n_distinct(Charge_log), .groups="drop") %>% pull(nu)
  k_fs_yh  <- max(3, min(5, max(3, min(min_u_yh, na.rm=TRUE) - 1)))
  
  # modèle : baseline + déviation Année×Habitat (on garde le reste inchangé)
  m07_hab2 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(Charge_log, k = 8, bs="tp") +                     # tendance globale
      s(Charge_log, year_hab, bs="fs", k = k_fs_yh, m=1) +# déviation Année×Habitat  ← NEW
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=30) +
      year_f + habitat4,
    data = dt_hab, method="fREML", select=TRUE, discrete=TRUE
  )
  # Si l’EDFi de s(Charge_log, year_hab, ...) ~ 0, c’est que les données n’appuient pas des différences fortes.
  # summary(m07_hab2)$s.tabl
  
  
  m07_hab <- m07hab
  
  
  
  
  # --- quantiles globaux (comparables)
  q <- quantile(dt_hab$Charge_log, probs = c(0, 0.5), na.rm = TRUE)
  q_lo <- unname(q[1]); q_hi <- unname(q[2])
  
  # --- base d'intégration commune pour lisser DAH/SMOD/x/y
  set.seed(1)
  base_design <- dt_hab %>%
    select(DAH, SMOD, x, y) %>%
    sample_n(size = min(5000, n()), replace = FALSE)
  
  lev_year <- levels(dt_hab$year_f)
  lev_hab  <- levels(dt_hab$habitat4)
  alp_ref  <- levels(dt_hab$alpage)[1]
  
  # --- on exclut les effets alpage & spatial dans la prédiction (effet "population-level")
  term_excl <- c("s(alpage)", "s(x,y)")
  
  mk_new <- function(df, cval, yy, hab) {
    mutate(df,
           Charge_log = cval,
           year_f     = factor(yy,  levels = lev_year),
           habitat4   = factor(hab, levels = lev_hab),
           alpage     = alp_ref
    )
  }
  
  # --- contraste pour une année et un groupe (habitat ou "Total")
  contrast_one <- function(yy, group = c("Total", as.character(lev_hab))) {
    group <- match.arg(group, c("Total", as.character(lev_hab)))
    
    if (group == "Total") {
      # poids égaux sur les 4 habitats
      base_all <- map_dfr(lev_hab, ~ base_design %>% mutate(habitat4 = .x))
    } else {
      base_all <- base_design %>% mutate(habitat4 = factor(group, levels = lev_hab))
    }
    
    new_lo <- mk_new(base_all, q_lo, yy, base_all$habitat4)
    new_hi <- mk_new(base_all, q_hi, yy, base_all$habitat4)
    
    X_lo <- predict(m07_hab, newdata = new_lo, type = "lpmatrix", exclude = term_excl)
    X_hi <- predict(m07_hab, newdata = new_hi, type = "lpmatrix", exclude = term_excl)
    
    Xbar  <- colMeans(X_hi - X_lo)
    b     <- coef(m07_hab)
    V     <- vcov(m07_hab, unconditional = TRUE)
    
    delta <- as.numeric(Xbar %*% b)
    se    <- sqrt(sum((Xbar %*% V) * Xbar))
    
    tibble::tibble(
      year_f   = yy,
      group    = group,
      delta_log = delta,
      delta_lo  = delta - 1.96*se,
      delta_hi  = delta + 1.96*se,
      ratio     = exp(delta),
      ratio_lo  = exp(delta - 1.96*se),
      ratio_hi  = exp(delta + 1.96*se)
    )
  }
  
  # --- calcul complet: pour chaque année -> Total + 4 habitats
  eff <- map_dfr(lev_year, function(yy) {
    bind_rows(
      contrast_one(yy, "Total"),
      map_dfr(as.list(lev_hab), ~ contrast_one(yy, .x))
    )
  }) %>%
    mutate(group = factor(group, levels = c("Total", as.character(lev_hab))))
  
  # --- plot: 5 points par année (Total + 4 habitats)
  ggplot(eff, aes(x = group, y = ratio)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ratio_lo, ymax = ratio_hi), width = 0.15) +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_wrap(~ year_f, scales = "free_y") +
    labs(
      x = NULL,
      y = "Taille d'effet de la charge (Q90 vs Q10, ratio sur GPROD)",
      title = "Effet du pâturage par année — Total et 4 habitats",
      subtitle = paste0("Contraste global ", round(q_lo,2), " → ", round(q_hi,2),
                        " ; covariables intégrées ; effets alpage & spatial exclus\n",
                        "Total = moyenne à parts égales des 4 habitats")
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 20, hjust = 1))
  
  
  
  
  
  
  
  
  
  # --- PRÉ-REQUIS : m07_hab déjà fit comme dans le message précédent
  # m07_hab <- bam(... s(Charge_log, by=year_f) + s(Charge_log, by=habitat4) + ...)
  
  library(dplyr)
  library(purrr)
  library(mgcv)
  library(ggplot2)
  
  # 1) Petites aides
  lev_year <- levels(dt_hab$year_f)
  lev_hab  <- levels(dt_hab$habitat4)
  alp_ref  <- levels(dt_hab$alpage)[1]
  
  # base d'intégration commune (DAH/SMOD/x/y)
  set.seed(1)
  base_design <- dt_hab %>%
    select(DAH, SMOD, x, y) %>%
    sample_n(size = min(4000, n()), replace = FALSE)
  
  # on exclut ces termes pour des courbes "population-level"
  term_excl <- c("s(alpage)", "s(x,y)")
  
  # grille de Charge_log (évite l'extrapolation)
  xr <- quantile(dt_hab$Charge_log, c(0.05, 0.95), na.rm = TRUE)
  grid_c <- as.numeric(seq(xr[1], xr[2], length.out = 120))
  
  # 2) Fonction qui calcule la courbe moyenne par habitat
  curve_for_hab <- function(hab, year_mode = c("average","fix"), year_fix = NULL) {
    year_mode <- match.arg(year_mode)
    # base: répète toutes les années (moyenne égale) ou fixe une année
    base0 <- base_design %>% mutate(habitat4 = factor(hab, levels = lev_hab))
    baseY <- if (year_mode == "average") {
      map_dfr(lev_year, ~ base0 %>% mutate(year_f = factor(.x, levels = lev_year)))
    } else {
      base0 %>% mutate(year_f = factor(year_fix, levels = lev_year))
    }
    
    beta <- coef(m07_hab)
    Vb   <- vcov(m07_hab, unconditional = TRUE)
    
    map_dfr(grid_c, function(cval) {
      newd <- baseY %>% mutate(
        Charge_log = cval,
        alpage     = alp_ref
      )
      # lpmatrix puis moyenne (delta-method)
      X <- predict(m07_hab, newdata = newd, type = "lpmatrix", exclude = term_excl)
      Xbar   <- colMeans(X)
      eta    <- as.numeric(Xbar %*% beta)
      se_eta <- sqrt(sum((Xbar %*% Vb) * Xbar))
      
      tibble::tibble(
        habitat4 = hab,
        Charge_log = cval,
        eta  = eta,                 # log-échelle
        lo   = eta - 1.96*se_eta,
        hi   = eta + 1.96*se_eta,
        mu   = exp(eta),            # échelle GPROD
        mu_lo= exp(lo),
        mu_hi= exp(hi)
      )
    })
  }
  
  # 3) Calcule pour les 4 habitats (moyenne sur les années)
  curves_hab <- map_dfr(lev_hab, curve_for_hab, year_mode = "average")
  
  # 4) Plot — courbe GPROD ~ Charge_log par habitat
  ggplot(curves_hab, aes(x = Charge_log, y = mu)) +
    geom_ribbon(aes(ymin = mu_lo, ymax = mu_hi), alpha = 0.15) +
    geom_line(size = 1) +
    facet_wrap(~ habitat4, scales = "free_y") +
    labs(
      x = "Charge (log)",
      y = "GPROD (prédit, échelle originale)",
      title = "Relation GPROD ~ Charge (log) par habitat",
      subtitle = "Courbes marginales (moyenne sur les années) — covariables intégrées, effets alpage & spatial exclus"
    ) +
    theme_minimal()
  
  
  
  
  
  
  
  
  
  
  m05 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      charge_bin +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40),   # ou s(x_km, y_km, ...)
    data = dt_legacy, method="fREML", discrete=TRUE, select=TRUE
  )
  
  summary(m04)
  
  AIC(m00 , m01, m02, m03, m04)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  library(mgcv)
  dt_legacy$year_f <- factor(dt_legacy$YEAR)
  dt_legacy$alpage <- factor(dt_legacy$alpage)
  
  
  
  
  
  
  m0 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40),
    data = dt_legacy, method="fREML", discrete=TRUE, select=TRUE
  )
  
  
  m1 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +charge_bin +                                      # décalage d'intercept (chargé vs non)
      s(Charge_log, by=charge_bin, bs="ts", k=10) +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40),      # pente lissée active si charge_bin==1
    data = dt_legacy, method="fREML", discrete=TRUE, select=TRUE
  )
  
  
  
  
  anova(m0, m1, test="Chisq"); AIC(m0, m1)
  dev0 <- deviance(m0); dev1 <- deviance(m1)
  partial_R2_charge <- (dev0 - dev1) / dev0; partial_R2_charge
  
  
  
  
  
  
  m2 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      charge_bin +
      s(Charge_log, by=charge_bin, bs="ts", k=10) +
      year_f,                                           # effets fixes par année
    data = dt_legacy, method="fREML", discrete=TRUE, select=TRUE
  )
  
  
  
  
  
  
  library(mgcv)
  
  # Prépa simple
  dt_legacy$charge_bin <- as.numeric(dt_legacy$charge_bin)  # 0/1
  dt_legacy$year_f <- relevel(factor(dt_legacy$year_f), ref = "2017")  # choisis 2017 comme référence
  
  # Modèle : topo-climat (GAM) + années fixes + intercept pâturé vs non
  #          + pente linéaire de Charge chez les pâturés
  #          + déviation de pente par année (interactions fixes)
  m2_est <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      year_f +                       # effets fixes moyens par année
      charge_bin +                   # intercept pâturé vs non
      Charge_log:charge_bin +        # pente globale chez les pâturés
      year_f:Charge_log:charge_bin,  # déviation de pente par année (fixe)
    data = dt_legacy, method="fREML", discrete=TRUE, select=TRUE
  )
  summary(m2_est)
  
  AIC(m2, m2_est)
  
  
  library(mgcv)
  library(ggplot2)
  
  # ---------- Prépa ----------
  dt_legacy$year_f     <- factor(dt_legacy$year_f)
  dt_legacy$charge_bin <- as.numeric(dt_legacy$charge_bin)   # 0/1
  dt_legacy$Charge_g   <- dt_legacy$Charge_log * dt_legacy$charge_bin
  years_all <- levels(dt_legacy$year_f)
  
  hab_map   <- c("1"="P. nivales",
                 "9"="P. thermiques écorchées",
                 "5"="Queyrellins",
                 "6"="Nardaies denses du subalpin")
  hab_codes <- names(hab_map)
  
  # ---------- Helpers pour dimensionner k en sécurité ----------
  n_uniq <- function(x) length(unique(x))
  n_uniq_xy <- function(d) length(unique(interaction(d$x, d$y, drop = TRUE)))
  
  safe_k_1d <- function(x, target=10, kmin=3) {
    ku <- n_uniq(x) - 1L
    max(kmin, min(target, ku))
  }
  safe_k_2d <- function(d, target=40, kmin=5) {
    ku <- n_uniq_xy(d) - 1L
    max(kmin, min(target, ku))
  }
  
  # ---------- Modèle : pente linéaire de Charge_g, pente qui varie par année ----------
  fit_lin_year <- function(dat){
    dat <- droplevels(dat)
    
    # Garde-fous
    has_grazing <- any(dat$charge_bin == 1, na.rm=TRUE)
    has_var     <- sd(dat$Charge_g[dat$charge_bin==1], na.rm=TRUE) > 0
    if(!has_grazing || !has_var) return(NULL)
    
    # k adaptés au sous-échantillon
    k_DAH  <- safe_k_1d(dat$DAH,  target=10, kmin=3)
    k_SMOD <- safe_k_1d(dat$SMOD, target=10, kmin=3)
    k_xy   <- safe_k_2d(dat,      target=40, kmin=5)
    k_ti1  <- safe_k_1d(dat$DAH,  target= min(20, k_DAH),  kmin=3)
    k_ti2  <- safe_k_1d(dat$SMOD, target= min(20, k_SMOD), kmin=3)
    
    bam(
      log(GPROD) ~ 
        s(DAH,  k = k_DAH,  bs = "ts") +
        s(SMOD, k = k_SMOD, bs = "ts") +
        ti(DAH, SMOD, k = c(k_ti1, k_ti2), bs = c("ts","ts")) +
        s(alpage, bs = "re") +
        s(x, y, bs = "tp", k = k_xy) +
        year_f +                 # effet moyen par année (capte les années "tôt/tard")
        charge_bin +             # intercept pâturé vs non
        Charge_g +               # pente moyenne (année de référence)
        year_f:Charge_g,         # déviation de pente par année (effets fixes)
      data = dat, method = "fREML",
      discrete = FALSE,          # plus sûr sur petits n
      select   = TRUE
    )
  }
  
  # ---------- Extraction des pentes par année -> % (+ IC95%) ----------
  extract_year_slopes <- function(mod, dat, label){
    if (is.null(mod)) {
      return(data.frame(year_f=years_all, group=label, effect_pct=NA, lwr=NA, upr=NA))
    }
    b  <- coef(mod); V <- vcov(mod); pn <- names(b)
    yrs_present <- levels(droplevels(dat$year_f))
    
    ix_base <- which(pn == "Charge_g")
    if(length(ix_base)!=1) stop("Coef 'Charge_g' introuvable (vérifie la formule).")
    
    res <- lapply(years_all, function(y){
      if(!(y %in% yrs_present)) {
        return(data.frame(year_f=y, group=label, effect_pct=NA, lwr=NA, upr=NA))
      }
      idx <- which(pn == paste0("year_f", y, ":Charge_g"))
      L <- rep(0, length(b)); L[ix_base] <- 1
      if(length(idx)==1) L[idx] <- 1
      
      est <- as.numeric(crossprod(L, b))
      se  <- sqrt(as.numeric(t(L) %*% V %*% L))
      
      data.frame(
        year_f = y, group = label,
        effect_pct = (exp(est) - 1) * 100,
        lwr        = (exp(est - 1.96*se) - 1) * 100,
        upr        = (exp(est + 1.96*se) - 1) * 100
      )
    })
    do.call(rbind, res)
  }
  
  # ---------- 1) Global ----------
  m_global   <- fit_lin_year(dt_legacy)
  res_global <- extract_year_slopes(m_global, dt_legacy, "Global")
  
  # ---------- 2) Par habitat ----------
  res_habs <- lapply(hab_codes, function(code){
    dat_h <- subset(dt_legacy, as.character(habitat_code) == code)
    dat_h$year_f <- factor(dat_h$year_f, levels = years_all)
    m_h   <- fit_lin_year(dat_h)
    extract_year_slopes(m_h, dat_h, hab_map[[code]])
  })
  res_habs <- do.call(rbind, res_habs)
  
  # ---------- 3) Combine & plot : 5 barres par année ----------
  res_all <- rbind(res_global, res_habs)
  res_all$group <- factor(res_all$group, levels=c("Global", unname(hab_map)))
  res_plot <- na.omit(res_all)
  
  gg <- ggplot(res_plot, aes(x=year_f, y=effect_pct, fill=group)) +
    geom_col(position=position_dodge(width=0.8), width=0.72) +
    geom_errorbar(aes(ymin=lwr, ymax=upr),
                  position=position_dodge(width=0.8), width=0.2) +
    labs(x="Année",
         y="Effet du chargement sur GPROD (%) par +1 de Charge_log",
         fill="Groupe",
         title="Taille d'effet (estimate) du chargement – Global et par habitat",
         subtitle="Effets partiels ajustés pour DAH, SMOD, DAH×SMOD, spatial (s(x,y)) et année (additive); k adaptés par habitat") +
    theme_minimal(base_size = 12)
  print(gg)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ## --- extraire pentes par année (+ IC) à partir de m2_est ---
  b <- coef(m2_est); V <- vcov(m2_est); pname <- names(b)
  years <- levels(dt_legacy$year_f)  # 2017 est la ref
  
  ix_base <- grep("^Charge_log:charge_bin$", pname)
  
  ix_year <- sapply(years[-1], function(y) {
    i <- grep(paste0("(^|:)", y, "($|:)"), pname, perl=TRUE)     # coef contenant l'année
    i <- i[grepl("Charge_log", pname[i]) & grepl("charge_bin", pname[i])]
    if (length(i) == 0) NA_integer_ else i
  }, simplify = FALSE)
  
  make_L <- function(y) {
    L <- rep(0, length(b)); L[ix_base] <- 1
    if (y != years[1] && !is.na(ix_year[[y]])) L[ix_year[[y]]] <- 1
    L
  }
  
  res <- do.call(rbind, lapply(years, function(y){
    L <- make_L(y)
    est <- as.numeric(crossprod(L, b))
    se  <- sqrt(as.numeric(t(L) %*% V %*% L))
    data.frame(
      year_f = y,
      effect_pct = (exp(est)-1)*100,
      lwr = (exp(est-1.96*se)-1)*100,
      upr = (exp(est+1.96*se)-1)*100
    )
  }))
  
  # --- barplot simple + IC95% ---
  op <- par(mar=c(4,5,2,1))
  bp <- barplot(res$effect_pct, names.arg=res$year_f, las=1,
                ylab="Effet du chargement sur GPROD (%), par +1 de Charge_log",
                ylim=range(c(0, res$lwr, res$upr)), border=NA)
  arrows(x0=bp, y0=res$lwr, x1=bp, y1=res$upr, angle=90, code=3, length=0.04)
  box(); par(op)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Assure-toi que charge_bin est bien numérique 0/1 (pas facteur)
  dt_legacy$charge_bin <- as.numeric(dt_legacy$charge_bin)
  
  m3_fix <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      year_f +                      # baseline interannuel
      charge_bin +                  # décalage d'intercept (grazed vs non)
      s(Charge_log, by=charge_bin, bs="ts", k=8) +               # effet moyen de dose, seulement si pâturé
      s(Charge_log, year_f, by=charge_bin, bs="fs", k=6, m=2),   # déviation par année, seulement si pâturé
    data = dt_legacy, method="fREML",
    # <<< astuce convergence >>>
    discrete = FALSE,                          # d'abord sans 'discrete' pour stabiliser
    select   = TRUE,
    gamma    = 1.2,                            # pénalisation un peu plus forte = plus stable
    optimizer = c("efs","newton")
  )
  
  
  summary(m3_fix)
  AIC(m0, m1, m2, m3_fix)
  
  
  # 0) Prépa : activer la charge seulement sur les pixels pâturés
  dt_legacy$charge_bin <- as.numeric(dt_legacy$charge_bin)  # 0/1
  dt_legacy$Charge_log_grazed <- dt_legacy$Charge_log * dt_legacy$charge_bin
  
  # 1) Modèle base + effet moyen non-linéaire de la charge + pente (linéaire) qui varie par année
  library(mgcv)
  
  m3_vc <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      year_f +                              # différences moyennes entre années
      charge_bin +                          # intercept grazed vs non
      s(Charge_log, by=charge_bin, bs="ts", k=8) +  # courbe dose-réponse moyenne (seulement si pâturé)
      s(year_f, by=Charge_log_grazed, bs="re"),     # *random slope* de Charge par année (seulement si pâturé)
    data = dt_legacy, method="fREML",
    discrete = FALSE,            # d’abord sans discretisation pour stabiliser
    select   = TRUE,
    gamma    = 1.2,
    optimizer = c("efs","newton")
  )
  
  
  summary(m3_vc)
  AIC(m0, m1, m2, m3_fix, m3_vc)
  
  
  
  
  
  # --- QUANTILES DE CHARGE CHEZ LES PATURES (base R) ---
  qs <- quantile(dt_legacy$Charge_log[dt_legacy$charge_bin == 1],
                 probs = c(0.25, 0.75), na.rm = TRUE)
  Q25 <- qs[[1]]; Q75 <- qs[[2]]
  
  # --- TERME A GARDER: random slope annuel ---
  lab_smooth <- vapply(m3_vc$smooth, function(s) s$label, character(1))
  keep_label <- lab_smooth[grepl("s\\(year_f\\).*Charge_log_grazed", lab_smooth)]
  
  # Exclure tous les autres lissages (on ne garde QUE le random slope annuel)
  exclude_labels <- setdiff(lab_smooth, keep_label)
  
  # --- NEWDATA: 2 lignes par année (Q25 et Q75), le reste n'a pas d'importance car exclu ---
  years <- levels(dt_legacy$year_f)
  alp0  <- if (is.factor(dt_legacy$alpage)) levels(dt_legacy$alpage)[1] else unique(dt_legacy$alpage)[1]
  
  nd <- do.call(rbind, lapply(years, function(y) {
    data.frame(
      year_f = factor(c(y, y), levels = years),
      # ces colonnes sont requises par la formule, mais leurs effets seront exclus :
      DAH = 0, SMOD = 0, x = 0, y = 0, alpage = alp0,
      charge_bin = 1,
      Charge_log = c(Q25, Q75),
      Charge_log_grazed = c(Q25, Q75)
    )
  }))
  nd$pair <- rep(seq_along(years), each = 2)
  
  # --- CONTRASTE (Q75 - Q25) EN NE GARDANT QUE LE TERME random slope annuel ---
  X <- predict(m3_vc, newdata = nd, type = "lpmatrix", exclude = exclude_labels)
  b <- coef(m3_vc); V <- vcov(m3_vc)
  
  split_rows <- split(seq_len(nrow(nd)), nd$pair)
  eff <- do.call(rbind, lapply(seq_along(split_rows), function(i) {
    id <- split_rows[[i]]
    cvec <- X[id[2], , drop = FALSE] - X[id[1], , drop = FALSE]  # Q75 - Q25
    est  <- as.numeric(cvec %*% b)
    se   <- sqrt(as.numeric(cvec %*% V %*% t(cvec)))
    data.frame(
      year_f = years[i],
      effect_pct = (exp(est) - 1) * 100,
      lwr = (exp(est - 1.96 * se) - 1) * 100,
      upr = (exp(est + 1.96 * se) - 1) * 100
    )
  }))
  
  # --- BARPLOT ---
  op <- par(mar = c(4, 5, 2, 1))
  bp <- barplot(eff$effect_pct, names.arg = eff$year_f, las = 1,
                ylab = "Effet interannuel du chargement sur GPROD (%)\n[Q75 vs Q25, terme 'random slope' seul]",
                ylim = range(c(0, eff$upr)), border = NA)
  arrows(x0 = bp, y0 = eff$lwr, x1 = bp, y1 = eff$upr, angle = 90, code = 3, length = 0.04)
  box()
  par(op)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  anova(m1, m2, test="Chisq"); AIC(m1, m2)
  summary(m2)$p.table[grep("^year_f", rownames(summary(m2)$p.table)), ]
  
  
  
  
  m3A <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD, k=10, bs="ts") +
      ti(DAH, SMOD, k=c(20,20), bs=c("ts","ts")) +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      year_f + 
      charge_bin +
      s(Charge_log, by=charge_bin, bs="ts", k=10) +  # forme moyenne quand chargé
      charge_bin:year_f +                            # intercept "chargé" qui change par année
      Charge_log:year_f:charge_bin,                  # pente (linéaire) qui change par année (uniquement chargé)
    data = dt_legacy, method="fREML", discrete=TRUE, select=TRUE
  )
  
  # Test de variation par année de l’effet charge :
  anova(m2, m3A, test="Chisq"); AIC(m2, m3A)
  
  
  pt <- summary(m3A)$p.table
  rn <- rownames(pt)
  rn
  
  idx_intercepts <- grepl("year_f", rn) & grepl("charge_bin", rn) & !grepl("Charge_log", rn)
  intercepts_yr <- pt[idx_intercepts, , drop=FALSE]
  intercepts_yr
  
  idx_slopes <- grepl("year_f", rn) & grepl("charge_bin", rn) & grepl("Charge_log", rn)
  slopes_yr <- pt[idx_slopes, , drop=FALSE]
  slopes_yr
  
  
  # Décomposition SMOD
  SMOD_mean_by_cell <- tapply(dt_legacy$SMOD, dt_legacy$cell, mean, na.rm=TRUE)
  dt_legacy$SMOD_mean_cell <- SMOD_mean_by_cell[as.character(dt_legacy$cell)]
  dt_legacy$SMOD_anom      <- dt_legacy$SMOD - dt_legacy$SMOD_mean_cell  # <0 = déneigement précoce / faible neige
  
  # Modèle "stress neige × charge"
  m_stress <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD_mean_cell, bs="ts") + 
      s(SMOD_anom,      bs="ts") +
      ti(DAH, SMOD_anom, k=c(20,20), bs=c("ts","ts")) +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      year_f +
      charge_bin +
      s(Charge_log, by=charge_bin, bs="ts", k=10) +
      ti(Charge_log, SMOD_anom, bs=c("ts","ts"), k=c(10,10)),  # <-- test clé
    data = dt_legacy, method="fREML", discrete=TRUE, select=TRUE
  )
  
  AIC(m2, m3A, m_stress)
  
  summary(m_stress)
  
  
  
  
  
  
  library(mgcv)
  
  # Facteurs
  dt_legacy$year_f <- factor(dt_legacy$YEAR)
  dt_legacy$alpage <- factor(dt_legacy$alpage)
  
  # Décomposition SMOD (base R, robuste)
  SMOD_mean_by_cell <- tapply(dt_legacy$SMOD, dt_legacy$cell, mean, na.rm=TRUE)
  dt_legacy$SMOD_mean_cell <- SMOD_mean_by_cell[as.character(dt_legacy$cell)]
  dt_legacy$SMOD_anom      <- dt_legacy$SMOD - dt_legacy$SMOD_mean_cell
  
  
  
  mS0 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD_mean_cell, bs="ts") +      # structure spatiale neige
      s(SMOD_anom,      bs="ts") +      # stress inter-annuel
      ti(DAH, SMOD_anom, k=c(20,20), bs=c("ts","ts")) +  # optionnel mais pertinent
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      charge_bin +
      s(Charge_log, by=charge_bin, bs="ts", k=10) +
      year_f,
    data = dt_legacy, method="fREML", discrete=TRUE, select=TRUE
  )
  
  
  mS1 <- bam(
    log(GPROD) ~ 
      s(DAH,  k=10, bs="ts") +
      s(SMOD_mean_cell, bs="ts") +
      s(SMOD_anom,      bs="ts") +
      ti(DAH, SMOD_anom, k=c(20,20), bs=c("ts","ts")) +
      s(alpage, bs="re") +
      s(x, y, bs="tp", k=40) +
      year_f +                                  # garde les écarts annuels résiduels
      charge_bin +
      s(Charge_log, by=charge_bin, bs="ts", k=10) +             # effet moyen (chargé)
      ti(Charge_log, SMOD_anom, bs=c("ts","ts"), k=c(10,10),    # <-- interaction clé
         by = charge_bin),
    data = dt_legacy, method="fREML", discrete=TRUE, select=TRUE
  )
  
  # Comparaison
  anova(mS0, mS1, test="Chisq"); AIC(mS0, mS1)
  
  
  
  summary(mS1)$s.table   # regarde la ligne 'ti(Charge_log,SMOD_anom):charge_bin...'
  gam.check(mS1)         # k-index; n'augmente k que si k-index<1 ET edf≈k'
  concurvity(mS1)    
  
  
  
  
  # Niveaux d'anomalie : précoce / médian / tardif
  q_anom <- quantile(dt_legacy$SMOD_anom, c(0.1, 0.5, 0.9), na.rm=TRUE)
  names(q_anom) <- c("precoce","moyenne","tardive")
  
  # Grille de charge (régime chargé)
  g <- expand.grid(
    Charge_log = seq(
      quantile(dt_legacy$Charge_log[dt_legacy$charge_bin==1], 0.02, na.rm=TRUE),
      quantile(dt_legacy$Charge_log[dt_legacy$charge_bin==1], 0.98, na.rm=TRUE),
      length.out=100
    ),
    SMOD_anom = q_anom
  )
  
  # Valeurs de référence (on neutralise tout le reste)
  g$charge_bin <- 1
  g$DAH  <- median(dt_legacy$DAH,  na.rm=TRUE)
  g$SMOD_mean_cell <- median(dt_legacy$SMOD_mean_cell, na.rm=TRUE)
  g$alpage <- levels(dt_legacy$alpage)[1]
  g$x <- median(dt_legacy$x); g$y <- median(dt_legacy$y)
  g$year_f <- levels(dt_legacy$year_f)[1]
  
  # Contraste "chargé" - "non chargé"
  g_on  <- g
  g_off <- g; g_off$charge_bin <- 0; g_off$Charge_log <- 0
  
  p_on  <- predict(mS1, newdata=g_on,  type="link")
  p_off <- predict(mS1, newdata=g_off, type="link")
  
  eff <- transform(g,
                   eff_log = p_on - p_off,
                   mult    = exp(p_on - p_off)  # facteur multiplicatif sur GPROD
  )
  # → trace 'mult' ~ Charge_log pour chaque niveau d'anomalie (précoce/moyenne/tardive).
  
  
  
  
  
  # Anomalie médiane de neige par année (observée)
  anom_year <- tapply(dt_legacy$SMOD_anom, dt_legacy$year_f, median, na.rm=TRUE)
  
  # Une grille de Charge_log et contraste par année
  gY <- data.frame(
    Charge_log = seq(
      quantile(dt_legacy$Charge_log[dt_legacy$charge_bin==1], 0.25, na.rm=TRUE),
      quantile(dt_legacy$Charge_log[dt_legacy$charge_bin==1], 0.75, na.rm=TRUE),
      length.out=5
    ),
    charge_bin = 1,
    DAH  = median(dt_legacy$DAH,  na.rm=TRUE),
    SMOD_mean_cell = median(dt_legacy$SMOD_mean_cell, na.rm=TRUE),
    alpage = levels(dt_legacy$alpage)[1],
    x = median(dt_legacy$x), y = median(dt_legacy$y)
  )
  
  yrs <- levels(dt_legacy$year_f)
  eff_year <- do.call(rbind, lapply(yrs, function(y){
    g_on  <- gY; g_off <- gY
    g_on$year_f  <- g_off$year_f <- factor(y, levels=yrs)
    g_on$SMOD_anom <- g_off$SMOD_anom <- as.numeric(anom_year[y])
    g_off$charge_bin <- 0; g_off$Charge_log <- 0
    p_on  <- predict(mS1, newdata=g_on,  type="link")
    p_off <- predict(mS1, newdata=g_off, type="link")
    data.frame(year=y, Charge_log=gY$Charge_log,
               effet_log=p_on-p_off, facteur=exp(p_on-p_off))
  }))
  eff_year
  
  
  
  
  
  
  # Choisir 3 niveaux d'anomalie : précoce (q10), moyenne (q50), tardive (q90)
  q_anom <- quantile(dt_legacy$SMOD_anom, c(0.1, 0.5, 0.9), na.rm=TRUE)
  names(q_anom) <- c("precoce","moyenne","tardive")
  
  # Grille de charge (régime chargé), avec covariables fixées
  g <- expand.grid(
    Charge_log = seq(
      quantile(dt_legacy$Charge_log[dt_legacy$charge_bin==1], 0.25, na.rm=TRUE),
      quantile(dt_legacy$Charge_log[dt_legacy$charge_bin==1], 0.75, na.rm=TRUE),
      length.out=5
    ),
    SMOD_anom = q_anom
  )
  g$charge_bin <- 1
  g$DAH  <- median(dt_legacy$DAH,  na.rm=TRUE)
  g$SMOD_mean_cell <- median(dt_legacy$SMOD_mean_cell, na.rm=TRUE)
  g$alpage <- levels(dt_legacy$alpage)[1]
  g$x <- median(dt_legacy$x); g$y <- median(dt_legacy$y)
  g$year_f <- levels(dt_legacy$year_f)[1]
  
  # Contraste "chargé" - "non chargé"
  g_on  <- g
  g_off <- g; g_off$charge_bin <- 0; g_off$Charge_log <- 0
  
  p_on  <- predict(mS1, newdata=g_on,  type="link")
  p_off <- predict(mS1, newdata=g_off, type="link")
  
  res <- transform(g,
                   effet_log = p_on - p_off,
                   facteur   = exp(p_on - p_off)   # multiplicateur sur GPROD
  )
  res[order(res$SMOD_anom, res$Charge_log), ]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # =========================================================
  # Packages
  # =========================================================
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  
  # =========================================================
  # 0) Standardisation des noms / colonnes
  #    (on ne touche pas Charge_log / charge_bin déjà présents)
  # =========================================================
  stopifnot(all(c("YEAR","alpage","GPROD","DAH","SMOD",
                  "Charge_log","charge_bin","habitat_code","habitat_label") %in% names(dt_legacy)))
  
  data_all <- dt_legacy %>%
    transmute(
      YEAR,
      alpage           = factor(alpage),
      AUCg             = GPROD,             # même indicateur que dans data_mod
      dah              = DAH,
      SMOD_2023        = SMOD,
      Charge_log       = Charge_log,
      charge_bin       = as.integer(ifelse(is.na(charge_bin), Charge_log >= log(2), charge_bin)),
      habitat_code     = as.integer(habitat_code),
      habitat_label    = habitat_label
    )
  
  # 4 habitats focus
  hab_labels_map <- c(
    "1" = "P. nivales",
    "9" = "P. thermiques écorchées",
    "5" = "Queyrellins",
    "6" = "Nardaies denses du subalpin"
  )
  focus_codes  <- c(1,9,5,6)
  focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
  
  # =========================================================
  # 1) Formule modèle (identique à avant)
  # =========================================================
  form_by_alp <- formula(
    log(AUCg) ~
      s(dah,       k = 10, bs = "ts") +
      s(SMOD_2023, k = 10, bs = "ts") +
      ti(dah, SMOD_2023, k = c(20, 20)) +
      s(alpage, bs = "re") +
      charge_bin +
      s(Charge_log, alpage, bs = "fs", k = 10, by = charge_bin) +
      ti(Charge_log, SMOD_2023, k = c(15, 15), by = charge_bin)
  )
  
  # =========================================================
  # 2) Fit par année (global + 4 habitats) avec bam
  # =========================================================
  n_th <- tryCatch(max(1L, parallel::detectCores(TRUE) - 1L), error = function(e) 1L)
  
  fit_fun <- function(dat) {
    dat <- dat %>%
      mutate(
        alpage     = factor(alpage),                    # refactorise au subset (crucial)
        charge_bin = as.integer(Charge_log >= log(2))   # gating cohérent
      )
    bam(form_by_alp, data = dat, method = "fREML",
        select = TRUE, discrete = TRUE, nthreads = n_th)
  }
  
  # listes par année
  years <- sort(unique(data_all$YEAR))
  d_by_year <- split(data_all, data_all$YEAR)
  
  # modèles globaux par année
  mglob <- map(d_by_year, fit_fun)
  
  # modèles par habitat (4 focus) par année
  mhab  <- map(d_by_year, function(dy) {
    map(setNames(as.list(focus_codes), as.character(focus_codes)), function(hc) {
      di <- dy %>% filter(habitat_code == hc)
      if (nrow(di) < 50) return(NULL)
      fit_fun(di)
    })
  })
  
  # =========================================================
  # 3) Effet Δ% entre Q25 -> Q75 de Charge_log (par an + groupes)
  # =========================================================
  delta_q25_q75 <- function(model, data_in, group_label = "Total") {
    d <- data_in
    pos <- d$Charge_log[d$Charge_log > 0]
    if (length(pos) < 10) {
      return(tibble::tibble(
        group = group_label,
        effect_log = NA_real_, se_log = NA_real_,
        lwr_log = NA_real_, upr_log = NA_real_,
        effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
        sig = FALSE, q_low = NA_real_, q_high = NA_real_, n = nrow(d)
      ))
    }
    ql <- as.numeric(quantile(pos, 0.25, na.rm = TRUE))
    qh <- as.numeric(quantile(pos, 0.75, na.rm = TRUE))
    if (!is.finite(ql) || !is.finite(qh) || ql >= qh) {
      ql <- as.numeric(quantile(pos, 0.10, na.rm = TRUE))
      qh <- as.numeric(quantile(pos, 0.90, na.rm = TRUE))
    }
    
    nd_lo <- d; nd_lo$Charge_log <- ql; nd_lo$charge_bin <- as.integer(ql >= log(2))
    nd_hi <- d; nd_hi$Charge_log <- qh; nd_hi$charge_bin <- as.integer(qh >= log(2))
    
    Xhi  <- predict(model, newdata = nd_hi, type = "lpmatrix")
    Xlo  <- predict(model, newdata = nd_lo, type = "lpmatrix")
    dbar <- colMeans(Xhi - Xlo)
    
    b <- coef(model); V <- vcov(model)
    eff <- as.numeric(drop(dbar %*% b))
    se  <- sqrt(as.numeric(drop(dbar %*% V %*% dbar)))
    lwr <- eff - 1.96*se; upr <- eff + 1.96*se
    
    tibble::tibble(
      group = group_label,
      effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
      effect_perc = 100*(exp(eff)-1),
      lwr_perc    = 100*(exp(lwr)-1),
      upr_perc    = 100*(exp(upr)-1),
      sig         = (lwr > 0) | (upr < 0),
      q_low = ql, q_high = qh, n = nrow(d)
    )
  }
  
  # collecte des effets pour toutes les années
  eff_res <- purrr::imap_dfr(d_by_year, function(dy, yy) {
    # total
    eg <- delta_q25_q75(mglob[[yy]], dy, group_label = "Total") %>%
      mutate(YEAR = as.integer(yy))
    # 4 habitats
    eh <- purrr::imap_dfr(mhab[[yy]], function(mod, code_chr) {
      lab <- hab_labels_map[[code_chr]]
      di  <- dy %>% filter(habitat_code == as.integer(code_chr))
      if (is.null(mod)) {
        tibble::tibble(
          group = lab, effect_log = NA_real_, se_log = NA_real_,
          lwr_log = NA_real_, upr_log = NA_real_,
          effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
          sig = FALSE, q_low = NA_real_, q_high = NA_real_, n = nrow(di),
          YEAR = as.integer(yy)
        )
      } else {
        delta_q25_q75(mod, di, group_label = lab) %>% mutate(YEAR = as.integer(yy))
      }
    })
    bind_rows(eg, eh)
  })
  
  # =========================================================
  # 4) Barplot (5 barres par année)
  # =========================================================
  eff_res <- eff_res %>%
    mutate(
      YEAR  = factor(YEAR, levels = years),
      group = factor(group, levels = c("Total", focus_labels))
    )
  
  cols_focus <- c(
    "Total"                        = "grey60",
    "P. nivales"                   = "#3B5BDB",
    "P. thermiques écorchées"      = "#E03131",
    "Queyrellins"                  = "#F4A261",
    "Nardaies denses du subalpin"  = "#1B9E77"
  )
  pal <- cols_focus[levels(eff_res$group)]
  
  p_mag_all <- ggplot(eff_res, aes(x = YEAR, y = effect_perc, fill = group)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"),
             width = 0.85, aes(alpha = ifelse(sig, 1, 0.4))) +
    geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc),
                  position = position_dodge2(width = 0.9, preserve = "single"),
                  width = 0.22) +
    scale_fill_manual(values = pal, name = "Groupe", drop = FALSE) +
    scale_alpha_identity() +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = "Effet de la charge sur la production (Δ% entre Q25→Q75 de Charge_log)",
      subtitle = "Total = modèle global annuel ; Habitats = modèles annuels séparés (codes 1, 9, 5, 6). Barres = IC95%.",
      x = "Année", y = "Variation (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")
  
  print(p_mag_all)
  # ggsave("barplot_magnitude_charge_2017_2023.png", p_mag_all, width = 12, height = 6.5, dpi = 300)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # =========================================================
  # EXTRACTION des effets (robuste) 2017–2023
  # =========================================================
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  
  eff_res <- purrr::imap_dfr(d_by_year, function(dy, yy) {
    # ----- Total (global annuel)
    eg <- delta_q25_q75_robust(mglob[[yy]], dy, group_label = "Total") %>%
      mutate(YEAR = as.integer(yy))
    
    # ----- 4 habitats focus (modèles séparés)
    eh <- purrr::imap_dfr(mhab[[yy]], function(mod, code_chr) {
      lab <- hab_labels_map[[code_chr]]
      di  <- dy %>% filter(habitat_code == as.integer(code_chr))
      if (is.null(mod)) {
        tibble::tibble(
          group = lab,
          effect_log = NA_real_, se_log = NA_real_,
          lwr_log = NA_real_, upr_log = NA_real_,
          effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
          sig = FALSE, q_low = NA_real_, q_high = NA_real_,
          n = nrow(di), kept_share = NA_real_, kept_n = nrow(di),
          adjusted = TRUE, YEAR = as.integer(yy)
        )
      } else {
        delta_q25_q75_robust(mod, di, group_label = lab) %>% mutate(YEAR = as.integer(yy))
      }
    })
    
    bind_rows(eg, eh)
  })
  
  # Mise en forme
  eff_res <- eff_res %>%
    mutate(
      YEAR  = factor(YEAR, levels = years),
      group = factor(group, levels = c("Total", focus_labels)),
      # Transparence : 1 si significatif, 0.4 sinon
      alpha_plot = ifelse(sig, 1, 0.4),
      # Type d'errorbar : "dashed" si quantiles réajustés (par manque de support)
      adj_flag   = factor(ifelse(isTRUE(adjusted), "Ajusté", "Nominal"),
                          levels = c("Nominal","Ajusté"))
    )
  
  # Palette
  cols_focus <- c(
    "Total"                        = "grey60",
    "P. nivales"                   = "#3B5BDB",
    "P. thermiques écorchées"      = "#E03131",
    "Queyrellins"                  = "#F4A261",
    "Nardaies denses du subalpin"  = "#1B9E77"
  )
  pal <- cols_focus[levels(eff_res$group)]
  
  # =========================================================
  # BARPLOT multi-années (Total + 4 habitats)
  # =========================================================
  p_mag_all <- ggplot(eff_res,
                      aes(x = YEAR, y = effect_perc, fill = group)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"),
             width = 0.85, aes(alpha = alpha_plot)) +
    geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc, linetype = adj_flag),
                  position = position_dodge2(width = 0.9, preserve = "single"),
                  width = 0.22) +
    scale_fill_manual(values = pal, name = "Groupe", drop = FALSE) +
    scale_alpha_identity() +
    scale_linetype_manual(values = c("solid","dashed"), name = "Bornes utilisées") +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = "Effet de la charge sur la production — 2017–2023",
      subtitle = "Δ% entre Q25 et Q75 de Charge_log (méthode robuste au support).\nBarres = IC95%. Lignes pointillées = quantiles réajustés (support insuffisant).",
      x = "Année", y = "Variation (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")
  
  print(p_mag_all)
  # ggsave("barplot_magnitude_charge_2017_2023_robuste.png", p_mag_all, width = 12, height = 6.5, dpi = 300)
  
  # =========================================================
  # BARPLOT 2023 uniquement (zoom)
  # =========================================================
  eff_2023 <- eff_res %>% filter(YEAR == max(YEAR))
  
  p_mag_2023 <- ggplot(eff_2023, aes(x = group, y = effect_perc, fill = group)) +
    geom_col(width = 0.80, aes(alpha = alpha_plot)) +
    geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc, linetype = adj_flag), width = 0.20) +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    scale_alpha_identity() +
    scale_linetype_manual(values = c("solid","dashed"), name = NULL) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = "Effet de la charge sur la production — 2023",
      subtitle = expression(Delta~"% entre "~Q[25]~" et "~Q[75]~" de "~Charge[log]~"(robuste)"),
      x = NULL, y = "Variation (%) vs Q25"
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(size = 12))
  
  print(p_mag_2023)
  # ggsave("barplot_magnitude_charge_2023_robuste.png", p_mag_2023, width = 9, height = 5, dpi = 300)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # =========================================================
  # Packages
  # =========================================================
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  
  # =========================================================
  # Données : dt_legacy  (colonnes vues dans ton summary)
  # =========================================================
  stopifnot(all(c("YEAR","alpage","GPROD","DAH","SMOD",
                  "Charge_log","charge_bin","habitat_code","habitat_label") %in% names(dt_legacy)))
  
  data_all <- dt_legacy %>%
    transmute(
      YEAR,
      alpage        = factor(alpage),
      AUCg          = GPROD,
      dah           = DAH,
      SMOD_2023     = SMOD,
      Charge_log    = Charge_log,
      charge_bin    = as.integer(ifelse(is.na(charge_bin), Charge_log >= log(2), charge_bin)),
      habitat_code  = as.integer(habitat_code),
      habitat_label = habitat_label
    )
  
  # =========================================================
  # Habitats focus
  # =========================================================
  hab_labels_map <- c(
    "1" = "P. nivales",
    "9" = "P. thermiques écorchées",
    "5" = "Queyrellins",
    "6" = "Nardaies denses du subalpin"
  )
  focus_codes  <- c(1,9,5,6)
  focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
  
  # =========================================================
  # Formule modèle (identique à avant, courbes par alpage)
  # =========================================================
  form_by_alp <- formula(
    log(AUCg) ~
      s(dah,       k = 10, bs = "ts") +
      s(SMOD_2023, k = 10, bs = "ts") +
      ti(dah, SMOD_2023, k = c(20,20)) +
      s(alpage, bs = "re") +
      charge_bin +
      s(Charge_log, alpage, bs = "fs", k = 10, by = charge_bin) +
      ti(Charge_log, SMOD_2023, k = c(15,15), by = charge_bin)
  )
  
  # =========================================================
  # Fit par année (global + 4 habitats)
  # =========================================================
  n_th <- tryCatch(max(1L, parallel::detectCores(TRUE) - 1L), error = function(e) 1L)
  fit_fun <- function(dat){
    dat <- dat %>% mutate(alpage = factor(alpage),
                          charge_bin = as.integer(Charge_log >= log(2)))
    bam(form_by_alp, data = dat, method = "fREML", select = TRUE,
        discrete = TRUE, nthreads = n_th)
  }
  
  years <- sort(unique(data_all$YEAR))
  d_by_year <- split(data_all, data_all$YEAR)
  mglob <- map(d_by_year, fit_fun)
  mhab  <- map(d_by_year, function(dy){
    map(setNames(as.list(focus_codes), as.character(focus_codes)), function(hc){
      di <- dy %>% filter(habitat_code == hc)
      if (nrow(di) < 50) return(NULL)
      fit_fun(di)
    })
  })
  
  # =========================================================
  # Aides : +10% et robust IQR
  # =========================================================
  log_plus_p <- function(x_log, p = 0.10) log1p((1+p) * pmax(expm1(x_log), 0))
  
  effect_plus10 <- function(model, data_in, group_label, p = 0.10){
    d <- data_in %>% filter(is.finite(Charge_log), Charge_log > 0)
    if (nrow(d) < 10) return(tibble(
      group = group_label, method = "Plus10 (fallback)",
      effect_log=NA_real_, se_log=NA_real_, lwr_log=NA_real_, upr_log=NA_real_,
      effect_perc=NA_real_, lwr_perc=NA_real_, upr_perc=NA_real_,
      sig=FALSE, q_low=NA_real_, q_high=NA_real_, p_lo_used=NA_real_, p_hi_used=NA_real_,
      kept_share=NA_real_, kept_n=nrow(d)
    ))
    nd0 <- d
    nd1 <- d; nd1$Charge_log <- log_plus_p(d$Charge_log, p)
    nd0$charge_bin <- as.integer(nd0$Charge_log >= log(2))
    nd1$charge_bin <- as.integer(nd1$Charge_log >= log(2))
    
    X1 <- predict(model, newdata = nd1, type="lpmatrix")
    X0 <- predict(model, newdata = nd0, type="lpmatrix")
    dbar <- colMeans(X1 - X0)
    
    b <- coef(model); V <- vcov(model)
    eff <- as.numeric(drop(dbar %*% b))
    se  <- sqrt(as.numeric(drop(dbar %*% V %*% dbar)))
    lwr <- eff - 1.96*se; upr <- eff + 1.96*se
    
    tibble(
      group = group_label, method = "Plus10 (fallback)",
      effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
      effect_perc = 100*(exp(eff)-1),
      lwr_perc    = 100*(exp(lwr)-1),
      upr_perc    = 100*(exp(upr)-1),
      sig = (lwr > 0) | (upr < 0),
      q_low = NA_real_, q_high = NA_real_,
      p_lo_used = NA_real_, p_hi_used = NA_real_,
      kept_share = nrow(d)/nrow(data_in), kept_n = nrow(d)
    )
  }
  
  # ---- IQR robuste : panel-safe + densité locale + resserrement si besoin ----
  delta_q25_q75_robust <- function(model, data_in, group_label,
                                   qlo_init = 0.25, qhi_init = 0.75,
                                   panel_trim = 0.05, width_x = 0.25,
                                   n_min_loc = 30, prop_min_loc = 0.01,
                                   min_share_rows = 0.30) {
    d0 <- data_in %>% filter(is.finite(Charge_log), Charge_log > 0)
    n_tot <- nrow(d0)
    if (n_tot < 10) return(effect_plus10(model, data_in, group_label))  # bascule direct
    
    # bornes par alpage (pour éviter les queues vides par panel)
    rng_by_alp <- d0 %>%
      group_by(alpage) %>%
      summarise(lo = quantile(Charge_log, panel_trim, na.rm = TRUE),
                hi = quantile(Charge_log, 1 - panel_trim, na.rm = TRUE),
                .groups = "drop")
    
    grid_p <- seq(qlo_init, 0.50, by = 0.05)     # 25→50
    pick <- NULL
    for (p_lo in grid_p) {
      p_hi <- 1 - p_lo
      ql <- as.numeric(quantile(d0$Charge_log, p_lo, na.rm = TRUE))
      qh <- as.numeric(quantile(d0$Charge_log, p_hi, na.rm = TRUE))
      if (!is.finite(ql) || !is.finite(qh) || ql >= qh) next
      
      ok_panels <- rng_by_alp %>% filter(ql >= lo, ql <= hi, qh >= lo, qh <= hi) %>% pull(alpage)
      d <- d0 %>% filter(alpage %in% ok_panels)
      if (nrow(d) == 0) next
      
      kept_share <- nrow(d) / n_tot
      n_loc_l <- sum(abs(d$Charge_log - ql) <= width_x)
      n_loc_h <- sum(abs(d$Charge_log - qh) <= width_x)
      ok_loc  <- (n_loc_l >= n_min_loc && n_loc_h >= n_min_loc &&
                    n_loc_l / n_tot >= prop_min_loc && n_loc_h / n_tot >= prop_min_loc)
      
      if (kept_share >= min_share_rows && ok_loc) {
        pick <- list(d=d, ql=ql, qh=qh, p_lo=p_lo, p_hi=p_hi, kept_share=kept_share); break
      }
    }
    
    if (is.null(pick)) {
      # si on n'a rien, on bascule proprement sur +10%
      return(effect_plus10(model, data_in, group_label))
    }
    
    d  <- pick$d; ql <- pick$ql; qh <- pick$qh
    
    nd_lo <- d; nd_lo$Charge_log <- ql; nd_lo$charge_bin <- as.integer(ql >= log(2))
    nd_hi <- d; nd_hi$Charge_log <- qh; nd_hi$charge_bin <- as.integer(qh >= log(2))
    
    Xhi <- predict(model, newdata = nd_hi, type = "lpmatrix")
    Xlo <- predict(model, newdata = nd_lo, type = "lpmatrix")
    dbar <- colMeans(Xhi - Xlo)
    
    b <- coef(model); V <- vcov(model)
    eff <- as.numeric(drop(dbar %*% b))
    se  <- sqrt(as.numeric(drop(dbar %*% V %*% dbar)))
    lwr <- eff - 1.96*se; upr <- eff + 1.96*se
    
    tibble(
      group = group_label,
      method = if (abs(pick$p_lo - qlo_init) < 1e-9) "IQR (Q25–Q75)" else sprintf("IQR ajusté (Q%.0f–Q%.0f)", 100*pick$p_lo, 100*pick$p_hi),
      effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
      effect_perc = 100*(exp(eff)-1),
      lwr_perc    = 100*(exp(lwr)-1),
      upr_perc    = 100*(exp(upr)-1),
      sig = (lwr > 0) | (upr < 0),
      q_low = ql, q_high = qh,
      p_lo_used = pick$p_lo, p_hi_used = pick$p_hi,
      kept_share = pick$kept_share, kept_n = nrow(d)
    )
  }
  
  # =========================================================
  # Extraction des effets pour 2017–2023
  # =========================================================
  safe_effect <- function(mod, dat, label){
    delta_q25_q75_robust(mod, dat, group_label = label)
  }
  
  eff_res <- imap_dfr(d_by_year, function(dy, yy){
    eg <- safe_effect(mglob[[yy]], dy, "Total") %>% mutate(YEAR = as.integer(yy))
    eh <- imap_dfr(mhab[[yy]], function(mod, code_chr){
      lab <- hab_labels_map[[code_chr]]
      di  <- dy %>% filter(habitat_code == as.integer(code_chr))
      if (is.null(mod)) {
        tibble(group = lab, method = "insuffisant",
               effect_log=NA_real_, se_log=NA_real_, lwr_log=NA_real_, upr_log=NA_real_,
               effect_perc=NA_real_, lwr_perc=NA_real_, upr_perc=NA_real_,
               sig=FALSE, q_low=NA_real_, q_high=NA_real_,
               p_lo_used=NA_real_, p_hi_used=NA_real_,
               kept_share=NA_real_, kept_n=nrow(di), YEAR = as.integer(yy))
      } else {
        safe_effect(mod, di, lab) %>% mutate(YEAR = as.integer(yy))
      }
    })
    bind_rows(eg, eh)
  })
  
  # =========================================================
  # Graphs
  # =========================================================
  eff_res <- eff_res %>%
    mutate(
      YEAR  = factor(YEAR, levels = years),
      group = factor(group, levels = c("Total", focus_labels)),
      alpha_plot = ifelse(sig, 1, 0.4),
      lty_plot   = ifelse(grepl("ajust", method), "dashed",
                          ifelse(grepl("Plus10", method), "twodash", "solid"))
    )
  
  cols_focus <- c(
    "Total"                        = "grey60",
    "P. nivales"                   = "#3B5BDB",
    "P. thermiques écorchées"      = "#E03131",
    "Queyrellins"                  = "#F4A261",
    "Nardaies denses du subalpin"  = "#1B9E77"
  )
  pal <- cols_focus[levels(eff_res$group)]
  
  # Multi-années
  p_mag_all <- ggplot(eff_res, aes(x = YEAR, y = effect_perc, fill = group)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"),
             width = 0.85, aes(alpha = alpha_plot)) +
    geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc, linetype = lty_plot),
                  position = position_dodge2(width = 0.9, preserve = "single"),
                  width = 0.22) +
    scale_fill_manual(values = pal, name = "Groupe", drop = FALSE) +
    scale_alpha_identity() +
    scale_linetype_identity(name = "Méthode",
                            guide = guide_legend(override.aes = list(color="black"))) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = "Effet de la charge sur la production (2017–2023)",
      subtitle = "Δ% entre Q25 et Q75 de Charge_log si support suffisant ; sinon IQR ajusté, sinon +10% (fallback).\nBarres = IC95%.",
      x = "Année", y = "Variation (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")
  
  print(p_mag_all)
  # ggsave("barplot_magnitude_charge_2017_2023_robuste.png", p_mag_all, width = 12, height = 6.5, dpi = 300)
  
  # Zoom 2023
  eff_2023 <- eff_res %>% filter(YEAR == max(YEAR))
  p_mag_2023 <- ggplot(eff_2023, aes(x = group, y = effect_perc, fill = group)) +
    geom_col(width = 0.80, aes(alpha = alpha_plot)) +
    geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc, linetype = lty_plot), width = 0.20) +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    scale_alpha_identity() +
    scale_linetype_identity() +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = "Effet de la charge sur la production — 2023",
      subtitle = "Méthode robuste au support (IQR / IQR ajusté / +10% fallback)",
      x = NULL, y = "Variation (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(size = 12))
  
  print(p_mag_2023)
  # ggsave("barplot_magnitude_charge_2023_robuste.png", p_mag_2023, width = 9, height = 5, dpi = 300)
  
  # (Option) Inspecter les cas pauvres en support
  # eff_res %>% arrange(YEAR, group) %>% select(YEAR, group, method, kept_n, kept_share, p_lo_used, p_hi_used, q_low, q_high) %>% print(n=Inf)
  
  
  
  
  
  
  
  
  
  
  
  # =========================================================
  # Packages
  # =========================================================
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(mgcv)
  
  # =========================================================
  # Helpers
  # =========================================================
  thr <- log(2)                     # seuil de charge (Charge = 1)
  log_plus_p <- function(x_log, p = 0.10) log1p((1+p) * pmax(expm1(x_log), 0))
  
  # Fallback: +10% mais UNIQUEMENT sur le POST-SEUIL
  effect_plus10_post <- function(model, data_in, group_label, p = 0.10, delta_thr = 0.02){
    d <- data_in %>% filter(is.finite(Charge_log), Charge_log >= thr + delta_thr)
    if (nrow(d) < 10) return(tibble::tibble(
      group = group_label, method = "Post +10% (fallback, insuffisant)",
      effect_log=NA_real_, se_log=NA_real_, lwr_log=NA_real_, upr_log=NA_real_,
      effect_perc=NA_real_, lwr_perc=NA_real_, upr_perc=NA_real_,
      sig=FALSE, from=NA_real_, to=NA_real_, kept_share=NA_real_, kept_n=nrow(d)
    ))
    
    nd0 <- d
    nd1 <- d; nd1$Charge_log <- log_plus_p(d$Charge_log, p)
    nd0$charge_bin <- as.integer(nd0$Charge_log >= thr)
    nd1$charge_bin <- as.integer(nd1$Charge_log >= thr)
    
    X1 <- predict(model, newdata=nd1, type="lpmatrix")
    X0 <- predict(model, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1 - X0)
    
    b <- coef(model); V <- vcov(model)
    eff <- as.numeric(drop(dbar %*% b))
    se  <- sqrt(as.numeric(drop(dbar %*% V %*% dbar)))
    lwr <- eff - 1.96*se; upr <- eff + 1.96*se
    
    tibble::tibble(
      group = group_label, method = "Post +10% (fallback)",
      effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
      effect_perc = 100*(exp(eff)-1),
      lwr_perc    = 100*(exp(lwr)-1),
      upr_perc    = 100*(exp(upr)-1),
      sig = (lwr > 0) | (upr < 0),
      from = NA_real_, to = NA_real_,
      kept_share = nrow(d)/nrow(data_in), kept_n = nrow(d)
    )
  }
  
  # Choix des deux points: "juste après seuil" -> p_to du POST-SEUIL
  pick_post_points <- function(d, p_to = 0.60, delta_thr = 0.02, cap_hi = 0.95){
    post <- d$Charge_log[d$Charge_log >= thr + delta_thr]
    if (!length(post)) return(c(from = NA_real_, to = NA_real_))
    x_from <- thr + delta_thr
    # on borne le haut pour éviter les queues : min(Qp_to, Q95 du post)
    x_to <- min(as.numeric(quantile(post, p_to, na.rm=TRUE)),
                as.numeric(quantile(post, cap_hi, na.rm=TRUE)))
    c(from = x_from, to = x_to)
  }
  
  # Contraste post-seuil robuste (panel-safe + densité locale), sinon fallback +10%
  effect_post_range_robust <- function(model, data_in, group_label,
                                       p_to = 0.60, delta_thr = 0.02, cap_hi = 0.95,
                                       panel_trim = 0.05, width_x = 0.25,
                                       n_min_loc = 30, prop_min_loc = 0.01,
                                       min_share_rows = 0.30){
    d0 <- data_in %>% filter(is.finite(Charge_log))
    if (nrow(d0) < 10) return(effect_plus10_post(model, data_in, group_label, delta_thr = delta_thr))
    
    # points à comparer
    pts <- pick_post_points(d0, p_to = p_to, delta_thr = delta_thr, cap_hi = cap_hi)
    x_from <- pts["from"]; x_to <- pts["to"]
    if (!is.finite(x_from) || !is.finite(x_to) || x_to <= x_from) {
      return(effect_plus10_post(model, data_in, group_label, delta_thr = delta_thr))
    }
    
    # bornes par alpage (éviter les queues propres à un panel)
    rng_by_alp <- d0 %>%
      group_by(alpage) %>%
      summarise(lo = quantile(Charge_log, panel_trim, na.rm = TRUE),
                hi = quantile(Charge_log, 1 - panel_trim, na.rm = TRUE),
                .groups = "drop")
    
    ok_panels <- rng_by_alp %>%
      filter(x_from >= lo, x_from <= hi, x_to >= lo, x_to <= hi) %>%
      pull(alpage)
    
    d <- d0 %>% filter(alpage %in% ok_panels)
    kept_share <- nrow(d)/max(nrow(d0),1)
    
    # densité locale près des deux points (fenêtre |x-x0|<=width_x)
    n_loc_from <- sum(abs(d$Charge_log - x_from) <= width_x)
    n_loc_to   <- sum(abs(d$Charge_log - x_to)   <= width_x)
    ok_loc <- (n_loc_from >= n_min_loc && n_loc_to >= n_min_loc &&
                 n_loc_from/nrow(d0) >= prop_min_loc && n_loc_to/nrow(d0) >= prop_min_loc)
    
    if (kept_share < min_share_rows || !ok_loc) {
      return(effect_plus10_post(model, data_in, group_label, delta_thr = delta_thr))
    }
    
    # lpmatrix diff + gating cohérent
    nd_a <- d; nd_a$Charge_log <- x_from; nd_a$charge_bin <- as.integer(x_from >= thr)
    nd_b <- d; nd_b$Charge_log <- x_to;   nd_b$charge_bin <- as.integer(x_to   >= thr)
    
    Xb <- predict(model, newdata = nd_b, type = "lpmatrix")
    Xa <- predict(model, newdata = nd_a, type = "lpmatrix")
    dbar <- colMeans(Xb - Xa)
    
    b <- coef(model); V <- vcov(model)
    eff <- as.numeric(drop(dbar %*% b))
    se  <- sqrt(as.numeric(drop(dbar %*% V %*% dbar)))
    lwr <- eff - 1.96*se; upr <- eff + 1.96*se
    
    tibble::tibble(
      group = group_label, method = sprintf("Post (%.0f%% → %.0f%% post-seuil)", 0, 100*p_to),
      effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
      effect_perc = 100*(exp(eff)-1),
      lwr_perc    = 100*(exp(lwr)-1),
      upr_perc    = 100*(exp(upr)-1),
      sig = (lwr > 0) | (upr < 0),
      from = x_from, to = x_to,
      kept_share = kept_share, kept_n = nrow(d)
    )
  }
  
  # =========================================================
  # Extraction multi-années : Total + 4 habitats, post-seuil
  #    Choisir p_to = 0.60 (ou 0.50) et delta_thr (ex. 0.02)
  # =========================================================
  compute_post_all <- function(d_by_year, mglob, mhab,
                               p_to = 0.60, delta_thr = 0.02){
    imap_dfr(d_by_year, function(dy, yy){
      # Global annuel
      eg <- effect_post_range_robust(mglob[[yy]], dy, "Total",
                                     p_to = p_to, delta_thr = delta_thr) %>%
        mutate(YEAR = as.integer(yy))
      
      # 4 habitats focus
      eh <- imap_dfr(mhab[[yy]], function(mod, code_chr){
        lab <- hab_labels_map[[code_chr]]
        di  <- dy %>% filter(habitat_code == as.integer(code_chr))
        if (is.null(mod)) {
          tibble::tibble(
            group=lab, method="Insuffisant",
            effect_log=NA_real_, se_log=NA_real_, lwr_log=NA_real_, upr_log=NA_real_,
            effect_perc=NA_real_, lwr_perc=NA_real_, upr_perc=NA_real_,
            sig=FALSE, from=NA_real_, to=NA_real_, kept_share=NA_real_, kept_n=nrow(di),
            YEAR = as.integer(yy)
          )
        } else {
          effect_post_range_robust(mod, di, lab, p_to = p_to, delta_thr = delta_thr) %>%
            mutate(YEAR = as.integer(yy))
        }
      })
      bind_rows(eg, eh)
    })
  }
  
  # ---- Lance le calcul (choisis p_to) ----
  # p_to = 0.60 (ou 0.50)
  eff_post <- compute_post_all(d_by_year, mglob, mhab, p_to = 0.60, delta_thr = 0.02)
  
  # =========================================================
  # Graph — 5 barres par année (Post-seuil)
  # =========================================================
  eff_post <- eff_post %>%
    mutate(
      YEAR  = factor(YEAR, levels = sort(unique(as.integer(YEAR)))),
      group = factor(group, levels = c("Total",
                                       "P. nivales","P. thermiques écorchées",
                                       "Queyrellins","Nardaies denses du subalpin")),
      alpha_plot = ifelse(sig, 1, 0.4),
      lty_plot   = ifelse(grepl("fallback", method, ignore.case = TRUE), "twodash", "solid")
    )
  
  cols_focus <- c(
    "Total"                        = "grey60",
    "P. nivales"                   = "#3B5BDB",
    "P. thermiques écorchées"      = "#E03131",
    "Queyrellins"                  = "#F4A261",
    "Nardaies denses du subalpin"  = "#1B9E77"
  )
  pal <- cols_focus[levels(eff_post$group)]
  
  p_post <- ggplot(eff_post, aes(x = YEAR, y = effect_perc, fill = group)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"),
             width = 0.85, aes(alpha = alpha_plot)) +
    geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc, linetype = lty_plot),
                  position = position_dodge2(width = 0.9, preserve = "single"),
                  width = 0.22) +
    scale_fill_manual(values = pal, name = "Groupe", drop = FALSE) +
    scale_alpha_identity() +
    scale_linetype_identity(name = "Méthode",
                            guide = guide_legend(override.aes = list(color="black"))) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = "Effet de la charge (POST-SEUIL) — 2017–2023",
      subtitle = "Contraste « juste après seuil → p_to du post-seuil ». Fallback = +10% post si support insuffisant. Barres = IC95%.",
      x = "Année", y = "Variation (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")
  
  print(p_post)
  # ggsave("barplot_post_only_2017_2023.png", p_post, width = 12, height = 6.5, dpi = 300)
  
  # ---- Zoom 2023 (facultatif) ----
  eff_post_2023 <- eff_post %>% filter(YEAR == max(YEAR))
  p_post_2023 <- ggplot(eff_post_2023, aes(x = group, y = effect_perc, fill = group)) +
    geom_col(width = 0.80, aes(alpha = alpha_plot)) +
    geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc, linetype = lty_plot), width = 0.20) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_alpha_identity() + scale_linetype_identity() +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(title = "Effet de la charge (POST-SEUIL) — 2023",
         subtitle = "Contraste « seuil+δ → p_to post »",
         x = NULL, y = "Variation (%)") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(size = 12))
  print(p_post_2023)
  # ggsave("barplot_post_only_2023.png", p_post_2023, width = 9, height = 5, dpi = 300)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # =========================================================
  # Packages
  # =========================================================
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  
  # =========================================================
  # Données : dt_legacy (2017–2023)
  # =========================================================
  stopifnot(all(c("YEAR","alpage","GPROD","DAH","SMOD",
                  "Charge_log","charge_bin","habitat_code","habitat_label") %in% names(dt_legacy)))
  
  data_all <- dt_legacy %>%
    transmute(
      YEAR,
      alpage        = factor(alpage),
      AUCg          = GPROD,
      dah           = DAH,
      SMOD_2023     = SMOD,
      Charge_log    = Charge_log,
      # sécurisation charge_bin (0/1)
      charge_bin    = as.integer(ifelse(is.na(charge_bin),
                                        (Charge_log >= log(2)), charge_bin)),
      habitat_code  = as.integer(habitat_code),
      habitat_label = habitat_label
    )
  
  # 4 habitats focus
  hab_labels_map <- c(
    "1" = "P. nivales",
    "9" = "P. thermiques écorchées",
    "5" = "Queyrellins",
    "6" = "Nardaies denses du subalpin"
  )
  focus_codes  <- c(1,9,5,6)
  focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
  
  # =========================================================
  # Modèle "binaire" (charge_bin seulement)
  # - base topo + neige
  # - aléatoire d'intercept par alpage
  # - effet paramétrique charge_bin (0→1)
  # - aléatoire "slope" charge_bin par alpage (random slope)
  # =========================================================
  form_bin <- formula(
    log(AUCg) ~
      s(dah,       k = 10, bs = "ts") +
      s(SMOD_2023, k = 10, bs = "ts") +
      ti(dah, SMOD_2023, k = c(20,20)) +
      s(alpage, bs = "re") +                       # intercept alpage
      charge_bin +                                  # effet moyen 0→1
      s(alpage, bs = "re", by = charge_bin)         # random slope 0→1 par alpage
  )
  
  # =========================================================
  # Fit par année (global + 4 habitats) avec bam
  # =========================================================
  n_th <- tryCatch(max(1L, parallel::detectCores(TRUE) - 1L), error = function(e) 1L)
  
  fit_fun <- function(dat){
    # refactoriser alpage dans le subset + regénérer charge_bin proprement
    dat <- dat %>% mutate(
      alpage     = factor(alpage),
      charge_bin = as.integer(Charge_log >= log(2))
    )
    # garde-fou : faut des 0 ET des 1 sinon le modèle binaire n'est pas identifiable
    if (length(unique(dat$charge_bin)) < 2L) return(NULL)
    
    bam(form_bin, data = dat, method = "fREML",
        select = TRUE, discrete = TRUE, nthreads = n_th)
  }
  
  years <- sort(unique(data_all$YEAR))
  d_by_year <- split(data_all, data_all$YEAR)
  
  # modèles globaux
  mglob <- map(d_by_year, fit_fun)
  
  # modèles par habitat
  mhab  <- map(d_by_year, function(dy){
    map(setNames(as.list(focus_codes), as.character(focus_codes)), function(hc){
      di <- dy %>% filter(habitat_code == hc)
      if (nrow(di) < 50) return(NULL)
      fit_fun(di)
    })
  })
  
  # =========================================================
  # Effet "0 → 1" (chargé vs non chargé) en % + IC95%
  # - On clone le data frame 2x : charge_bin=0 et charge_bin=1
  # - Différence du prédicteur moyen (lpmatrix), puis exp(Δ)-1 en %
  # - Garde-fou : au moins Nmin par état
  # =========================================================
  effect_bin01 <- function(model, data_in, group_label = "Total",
                           nmin_per_state = 30, prop_min_state = 0.05){
    d <- data_in %>% mutate(charge_bin = as.integer(Charge_log >= log(2)))
    n0 <- sum(d$charge_bin == 0, na.rm = TRUE)
    n1 <- sum(d$charge_bin == 1, na.rm = TRUE)
    if (is.null(model) || n0 < nmin_per_state || n1 < nmin_per_state ||
        (n0/nrow(d) < prop_min_state) || (n1/nrow(d) < prop_min_state)) {
      return(tibble(
        group = group_label, method = "0→1 (insuffisant)",
        effect_log = NA_real_, se_log = NA_real_, lwr_log = NA_real_, upr_log = NA_real_,
        effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
        sig = FALSE, n = nrow(d), n0 = n0, n1 = n1
      ))
    }
    
    nd0 <- d; nd0$charge_bin <- 0L
    nd1 <- d; nd1$charge_bin <- 1L
    
    X1 <- predict(model, newdata = nd1, type = "lpmatrix")
    X0 <- predict(model, newdata = nd0, type = "lpmatrix")
    dbar <- colMeans(X1 - X0)
    
    b <- coef(model); V <- vcov(model)
    eff <- as.numeric(drop(dbar %*% b))
    se  <- sqrt(as.numeric(drop(dbar %*% V %*% dbar)))
    lwr <- eff - 1.96*se; upr <- eff + 1.96*se
    
    tibble(
      group = group_label, method = "0→1",
      effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
      effect_perc = 100*(exp(eff)-1),
      lwr_perc    = 100*(exp(lwr)-1),
      upr_perc    = 100*(exp(upr)-1),
      sig = (lwr > 0) | (upr < 0),
      n = nrow(d), n0 = n0, n1 = n1
    )
  }
  
  # =========================================================
  # Extraction multi-années (Total + 4 habitats)
  # =========================================================
  eff_bin <- purrr::imap_dfr(d_by_year, function(dy, yy){
    # global
    eg <- effect_bin01(mglob[[yy]], dy, "Total") %>% mutate(YEAR = as.integer(yy))
    # habitats
    eh <- purrr::imap_dfr(mhab[[yy]], function(mod, code_chr){
      lab <- hab_labels_map[[code_chr]]
      di  <- dy %>% filter(habitat_code == as.integer(code_chr))
      effect_bin01(mod, di, lab) %>% mutate(YEAR = as.integer(yy))
    })
    bind_rows(eg, eh)
  })
  
  # =========================================================
  # Graphs
  # =========================================================
  eff_bin <- eff_bin %>%
    mutate(
      YEAR  = factor(YEAR, levels = years),
      group = factor(group, levels = c("Total", focus_labels)),
      alpha_plot = ifelse(sig, 1, 0.4),
      lty_plot   = ifelse(grepl("insuffisant", method, ignore.case = TRUE), "twodash", "solid")
    )
  
  cols_focus <- c(
    "Total"                        = "grey60",
    "P. nivales"                   = "#3B5BDB",
    "P. thermiques écorchées"      = "#E03131",
    "Queyrellins"                  = "#F4A261",
    "Nardaies denses du subalpin"  = "#1B9E77"
  )
  pal <- cols_focus[levels(eff_bin$group)]
  
  # ----- Multi-années -----
  p_bin_all <- ggplot(eff_bin, aes(x = YEAR, y = effect_perc, fill = group)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"),
             width = 0.85, aes(alpha = alpha_plot)) +
    geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc, linetype = lty_plot),
                  position = position_dodge2(width = 0.9, preserve = "single"),
                  width = 0.22) +
    scale_fill_manual(values = pal, name = "Groupe", drop = FALSE) +
    scale_alpha_identity() +
    scale_linetype_identity(name = "Méthode",
                            guide = guide_legend(override.aes = list(color="black"))) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = "Effet binaire de la charge (0 → 1) — 2017–2023",
      subtitle = "Modèle GAM : charge traitée en binaire. Effet = Δ% (exp(Δlog)−1) entre états non chargé et chargé.\nBarres = IC95%. 'twodash' = données insuffisantes (pas d'estimation).",
      x = "Année", y = "Variation (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")
  
  print(p_bin_all)
  # ggsave("barplot_charge_binaire_2017_2023.png", p_bin_all, width = 12, height = 6.5, dpi = 300)
  
  # ----- Zoom 2023 -----
  eff_bin_2023 <- eff_bin %>% filter(YEAR == max(YEAR))
  p_bin_2023 <- ggplot(eff_bin_2023, aes(x = group, y = effect_perc, fill = group)) +
    geom_col(width = 0.80, aes(alpha = alpha_plot)) +
    geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc, linetype = lty_plot), width = 0.20) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_alpha_identity() + scale_linetype_identity() +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(title = "Effet binaire de la charge (0 → 1) — 2023",
         x = NULL, y = "Variation (%)") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(size = 12))
  print(p_bin_2023)
  # ggsave("barplot_charge_binaire_2023.png", p_bin_2023, width = 9, height = 5, dpi = 300)
  
  
  
  
  
  
  
  # ===== Modèles annuels combinés (global + déviations par habitat) =====
  library(dplyr)
  library(mgcv)
  library(parallel)
  
  nthreads <- max(1, parallel::detectCores() - 1)
  
  # Sécurité réponse : si AUCg absent => utiliser GPROD
  if (!"AUCg" %in% names(dt_legacy) && "GPROD" %in% names(dt_legacy)) {
    dt_legacy$AUCg <- dt_legacy$GPROD
  }
  dt_legacy$AUCg <- pmax(dt_legacy$AUCg, .Machine$double.eps)
  
  # ========= 2017 =========
  d2017 <- dt_legacy %>%
    filter(YEAR == 2017) %>%
    mutate(
      DAH_std        = as.numeric(scale(DAH)),
      SMOD_std       = as.numeric(scale(SMOD)),
      SMOD_2023_std  = SMOD_std,
      alpage         = factor(alpage),
      hab_all        = factor(as.character(habitat_code)),
      h1 = as.numeric(habitat_code == 1),
      h9 = as.numeric(habitat_code == 9),
      h5 = as.numeric(habitat_code == 5),
      h6 = as.numeric(habitat_code == 6)
    ) %>% filter(is.finite(Charge_log))
  
  mod_2017 <- bam(
    log(AUCg) ~
      s(DAH_std, k=10, bs="ts") +
      s(SMOD_2023_std, k=10, bs="ts") +
      ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
      hab_all + s(alpage, bs="re") +
      # ----- Gradient GLOBAL + déviations pour 1,9,5,6 -----
    s(Charge_log, k=10, bs="ts") +
      s(Charge_log, k=10, bs="ts", by=h1) +
      s(Charge_log, k=10, bs="ts", by=h9) +
      s(Charge_log, k=10, bs="ts", by=h5) +
      s(Charge_log, k=10, bs="ts", by=h6) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15)) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h1) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h9) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h5) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h6),
    data=d2017, method="fREML", select=TRUE,
    discrete=TRUE, nthreads=nthreads, na.action=na.exclude
  )
  
  # ========= 2018 =========
  d2018 <- dt_legacy %>%
    filter(YEAR == 2018) %>%
    mutate(
      DAH_std = as.numeric(scale(DAH)),
      SMOD_std = as.numeric(scale(SMOD)),
      SMOD_2023_std = SMOD_std,
      alpage = factor(alpage),
      hab_all = factor(as.character(habitat_code)),
      h1 = as.numeric(habitat_code == 1),
      h9 = as.numeric(habitat_code == 9),
      h5 = as.numeric(habitat_code == 5),
      h6 = as.numeric(habitat_code == 6)
    ) %>% filter(is.finite(Charge_log))
  
  mod_2018 <- bam(
    log(AUCg) ~
      s(DAH_std, k=10, bs="ts") +
      s(SMOD_2023_std, k=10, bs="ts") +
      ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
      hab_all + s(alpage, bs="re") +
      s(Charge_log, k=10, bs="ts") +
      s(Charge_log, k=10, bs="ts", by=h1) +
      s(Charge_log, k=10, bs="ts", by=h9) +
      s(Charge_log, k=10, bs="ts", by=h5) +
      s(Charge_log, k=10, bs="ts", by=h6) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15)) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h1) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h9) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h5) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h6),
    data=d2018, method="fREML", select=TRUE,
    discrete=TRUE, nthreads=nthreads, na.action=na.exclude
  )
  
  # ========= 2019 =========
  d2019 <- dt_legacy %>%
    filter(YEAR == 2019) %>%
    mutate(
      DAH_std = as.numeric(scale(DAH)),
      SMOD_std = as.numeric(scale(SMOD)),
      SMOD_2023_std = SMOD_std,
      alpage = factor(alpage),
      hab_all = factor(as.character(habitat_code)),
      h1 = as.numeric(habitat_code == 1),
      h9 = as.numeric(habitat_code == 9),
      h5 = as.numeric(habitat_code == 5),
      h6 = as.numeric(habitat_code == 6)
    ) %>% filter(is.finite(Charge_log))
  
  mod_2019 <- bam(
    log(AUCg) ~
      s(DAH_std, k=10, bs="ts") +
      s(SMOD_2023_std, k=10, bs="ts") +
      ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
      hab_all + s(alpage, bs="re") +
      s(Charge_log, k=10, bs="ts") +
      s(Charge_log, k=10, bs="ts", by=h1) +
      s(Charge_log, k=10, bs="ts", by=h9) +
      s(Charge_log, k=10, bs="ts", by=h5) +
      s(Charge_log, k=10, bs="ts", by=h6) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15)) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h1) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h9) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h5) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h6),
    data=d2019, method="fREML", select=TRUE,
    discrete=TRUE, nthreads=nthreads, na.action=na.exclude
  )
  
  # ========= 2020 =========
  d2020 <- dt_legacy %>%
    filter(YEAR == 2020) %>%
    mutate(
      DAH_std = as.numeric(scale(DAH)),
      SMOD_std = as.numeric(scale(SMOD)),
      SMOD_2023_std = SMOD_std,
      alpage = factor(alpage),
      hab_all = factor(as.character(habitat_code)),
      h1 = as.numeric(habitat_code == 1),
      h9 = as.numeric(habitat_code == 9),
      h5 = as.numeric(habitat_code == 5),
      h6 = as.numeric(habitat_code == 6)
    ) %>% filter(is.finite(Charge_log))
  
  mod_2020 <- bam(
    log(AUCg) ~
      s(DAH_std, k=10, bs="ts") +
      s(SMOD_2023_std, k=10, bs="ts") +
      ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
      hab_all + s(alpage, bs="re") +
      s(Charge_log, k=10, bs="ts") +
      s(Charge_log, k=10, bs="ts", by=h1) +
      s(Charge_log, k=10, bs="ts", by=h9) +
      s(Charge_log, k=10, bs="ts", by=h5) +
      s(Charge_log, k=10, bs="ts", by=h6) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15)) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h1) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h9) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h5) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h6),
    data=d2020, method="fREML", select=TRUE,
    discrete=TRUE, nthreads=nthreads, na.action=na.exclude
  )
  
  # ========= 2021 =========
  d2021 <- dt_legacy %>%
    filter(YEAR == 2021) %>%
    mutate(
      DAH_std = as.numeric(scale(DAH)),
      SMOD_std = as.numeric(scale(SMOD)),
      SMOD_2023_std = SMOD_std,
      alpage = factor(alpage),
      hab_all = factor(as.character(habitat_code)),
      h1 = as.numeric(habitat_code == 1),
      h9 = as.numeric(habitat_code == 9),
      h5 = as.numeric(habitat_code == 5),
      h6 = as.numeric(habitat_code == 6)
    ) %>% filter(is.finite(Charge_log))
  
  mod_2021 <- bam(
    log(AUCg) ~
      s(DAH_std, k=10, bs="ts") +
      s(SMOD_2023_std, k=10, bs="ts") +
      ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
      hab_all + s(alpage, bs="re") +
      s(Charge_log, k=10, bs="ts") +
      s(Charge_log, k=10, bs="ts", by=h1) +
      s(Charge_log, k=10, bs="ts", by=h9) +
      s(Charge_log, k=10, bs="ts", by=h5) +
      s(Charge_log, k=10, bs="ts", by=h6) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15)) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h1) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h9) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h5) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h6),
    data=d2021, method="fREML", select=TRUE,
    discrete=TRUE, nthreads=nthreads, na.action=na.exclude
  )
  
  # ========= 2022 =========
  d2022 <- dt_legacy %>%
    filter(YEAR == 2022) %>%
    mutate(
      DAH_std = as.numeric(scale(DAH)),
      SMOD_std = as.numeric(scale(SMOD)),
      SMOD_2023_std = SMOD_std,
      alpage = factor(alpage),
      hab_all = factor(as.character(habitat_code)),
      h1 = as.numeric(habitat_code == 1),
      h9 = as.numeric(habitat_code == 9),
      h5 = as.numeric(habitat_code == 5),
      h6 = as.numeric(habitat_code == 6)
    ) %>% filter(is.finite(Charge_log))
  
  mod_2022 <- bam(
    log(AUCg) ~
      s(DAH_std, k=10, bs="ts") +
      s(SMOD_2023_std, k=10, bs="ts") +
      ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
      hab_all + s(alpage, bs="re") +
      s(Charge_log, k=10, bs="ts") +
      s(Charge_log, k=10, bs="ts", by=h1) +
      s(Charge_log, k=10, bs="ts", by=h9) +
      s(Charge_log, k=10, bs="ts", by=h5) +
      s(Charge_log, k=10, bs="ts", by=h6) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15)) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h1) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h9) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h5) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h6),
    data=d2022, method="fREML", select=TRUE,
    discrete=TRUE, nthreads=nthreads, na.action=na.exclude
  )
  
  # ========= 2023 =========
  d2023 <- dt_legacy %>%
    filter(YEAR == 2023) %>%
    mutate(
      DAH_std = as.numeric(scale(DAH)),
      SMOD_std = as.numeric(scale(SMOD)),
      SMOD_2023_std = SMOD_std,
      alpage = factor(alpage),
      hab_all = factor(as.character(habitat_code)),
      h1 = as.numeric(habitat_code == 1),
      h9 = as.numeric(habitat_code == 9),
      h5 = as.numeric(habitat_code == 5),
      h6 = as.numeric(habitat_code == 6)
    ) %>% filter(is.finite(Charge_log))
  
  mod_2023 <- bam(
    log(AUCg) ~
      s(DAH_std, k=10, bs="ts") +
      s(SMOD_2023_std, k=10, bs="ts") +
      ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
      hab_all + s(alpage, bs="re") +
      s(Charge_log, k=10, bs="ts") +
      s(Charge_log, k=10, bs="ts", by=h1) +
      s(Charge_log, k=10, bs="ts", by=h9) +
      s(Charge_log, k=10, bs="ts", by=h5) +
      s(Charge_log, k=10, bs="ts", by=h6) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15)) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h1) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h9) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h5) +
      ti(Charge_log, SMOD_2023_std, k=c(15,15), by=h6),
    data=d2023, method="fREML", select=TRUE,
    discrete=TRUE, nthreads=nthreads, na.action=na.exclude
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # === % de variation de GPROD pour +10% de charge (exact, log1p) ===
  library(dplyr); library(ggplot2); library(mgcv); library(tidyr)
  
  hab_labels_map <- c("1"="P. nivales","9"="P. thermiques écorchées",
                      "5"="Queyrellins","6"="Nardaies denses du subalpin")
  hab_levels <- unname(hab_labels_map[c("1","9","5","6")])
  cols_focus <- c("Total"="grey60","P. nivales"="#3B5BDB","P. thermiques écorchées"="#E03131",
                  "Queyrellins"="#F4A261","Nardaies denses du subalpin"="#1B9E77")
  
  # Modèles & data déjà en mémoire (ceux que tu as fit) :
  mglob <- list(`2017`=mod_2017_global, `2018`=mod_2018_global, `2019`=mod_2019_global,
                `2020`=mod_2020_global, `2021`=mod_2021_global, `2022`=mod_2022_global, `2023`=mod_2023_global)
  dglob <- list(`2017`=d2017_all,      `2018`=d2018_all,      `2019`=d2019_all,
                `2020`=d2020_all,      `2021`=d2021_all,      `2022`=d2022_all,      `2023`=d2023_all)
  mhab  <- list(`2017`=mod_2017, `2018`=mod_2018, `2019`=mod_2019,
                `2020`=mod_2020, `2021`=mod_2021, `2022`=mod_2022, `2023`=mod_2023)
  dhab  <- list(`2017`=d2017,    `2018`=d2018,    `2019`=d2019,
                `2020`=d2020,    `2021`=d2021,    `2022`=d2022,    `2023`=d2023)
  
  raw_col <- "chargement_median_2022-2024"  # ta charge brute
  eff10 <- tibble::tibble(YEAR=integer(), group=character(),
                          pct10=numeric(), lwr=numeric(), upr=numeric(), n=integer())
  
  for (yy in names(mglob)) {
    # ----- Global -----
    mg <- mglob[[yy]]; dg <- dglob[[yy]]
    d0 <- dg; d1 <- dg
    d0$Charge_log <- log1p(pmax(dg[[raw_col]], 0))
    d1$Charge_log <- log1p(pmax(dg[[raw_col]]*1.10, 0))
    X1 <- predict(mg, newdata=d1, type="lpmatrix")
    X0 <- predict(mg, newdata=d0, type="lpmatrix")
    dbar <- colMeans(X1 - X0); b <- coef(mg); V <- vcov(mg)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    pct <- 100*(exp(est)-1); lwr <- 100*(exp(est-1.96*se)-1); upr <- 100*(exp(est+1.96*se)-1)
    eff10 <- bind_rows(eff10, tibble::tibble(YEAR=as.integer(yy), group="Total",
                                             pct10=pct, lwr=lwr, upr=upr, n=nrow(dg)))
    # ----- Habitats -----
    mh <- mhab[[yy]]; dh <- dhab[[yy]]
    for (h in hab_levels) {
      dH <- dh %>% filter(hab4==h)
      if (nrow(dH)<30) next
      d0 <- dH; d1 <- dH
      d0$Charge_log <- log1p(pmax(dH[[raw_col]], 0))
      d1$Charge_log <- log1p(pmax(dH[[raw_col]]*1.10, 0))
      X1 <- predict(mh, newdata=d1, type="lpmatrix")
      X0 <- predict(mh, newdata=d0, type="lpmatrix")
      dbar <- colMeans(X1 - X0); b <- coef(mh); V <- vcov(mh)
      est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
      pct <- 100*(exp(est)-1); lwr <- 100*(exp(est-1.96*se)-1); upr <- 100*(exp(est+1.96*se)-1)
      eff10 <- bind_rows(eff10, tibble::tibble(YEAR=as.integer(yy), group=h,
                                               pct10=pct, lwr=lwr, upr=upr, n=nrow(dH)))
    }
  }
  
  eff10 <- eff10 %>% mutate(YEAR=factor(YEAR, levels=2017:2023),
                            group=factor(group, levels=c("Total", hab_levels)))
  
  pal <- cols_focus[levels(eff10$group)]
  p_elast <- ggplot(eff10, aes(YEAR, pct10, fill=group)) +
    geom_col(position=position_dodge2(width=0.9, preserve="single"), width=0.85) +
    geom_errorbar(aes(ymin=lwr, ymax=upr),
                  position=position_dodge2(width=0.9, preserve="single"), width=0.22) +
    scale_fill_manual(values=pal, name="Groupe", drop=FALSE) +
    geom_hline(yintercept=0, linewidth=0.5) +
    labs(title="Élasticité moyenne : effet d’un +10% de charge",
         subtitle="% de variation de la production attendue (GPROD) — IC95%",
         x="Année", y="% pour +10% de charge") +
    theme_minimal(base_size=14) +
    theme(panel.grid.minor=element_blank(), legend.position="right")
  
  p_elast
  
  
  
  
  
  
  
  
  
  
  
  
  # === % de variation de GPROD pour un doublement (médiane -> 2× médiane) ===
  eff2x <- tibble::tibble(YEAR=integer(), group=character(),
                          pct2x=numeric(), lwr=numeric(), upr=numeric(), n=integer())
  
  for (yy in names(mglob)) {
    # quantiles sur la charge brute (pour clamp)
    dg <- dglob[[yy]]
    lo <- quantile(pmax(dg[[raw_col]],0), 0.01, na.rm=TRUE)
    hi <- quantile(pmax(dg[[raw_col]],0), 0.99, na.rm=TRUE)
    
    # ----- Global -----
    mg <- mglob[[yy]]
    med <- as.numeric(median(pmax(dg[[raw_col]],0), na.rm=TRUE))
    c0 <- max(lo, min(med, hi))
    c1 <- max(lo, min(2*med, hi))
    nd0 <- dg; nd1 <- dg
    nd0$Charge_log <- log1p(c0); nd1$Charge_log <- log1p(c1)
    X1 <- predict(mg, newdata=nd1, type="lpmatrix")
    X0 <- predict(mg, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1 - X0); b <- coef(mg); V <- vcov(mg)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    pct <- 100*(exp(est)-1); lwr <- 100*(exp(est-1.96*se)-1); upr <- 100*(exp(est+1.96*se)-1)
    eff2x <- bind_rows(eff2x, tibble::tibble(YEAR=as.integer(yy), group="Total",
                                             pct2x=pct, lwr=lwr, upr=upr, n=nrow(dg)))
    
    # ----- Habitats -----
    mh <- mhab[[yy]]; dh <- dhab[[yy]]
    for (h in hab_levels) {
      dH <- dh %>% filter(hab4==h); if (nrow(dH)<30) next
      medH <- as.numeric(median(pmax(dH[[raw_col]],0), na.rm=TRUE))
      c0 <- max(lo, min(medH, hi)); c1 <- max(lo, min(2*medH, hi))
      nd0 <- dH; nd1 <- dH
      nd0$Charge_log <- log1p(c0); nd1$Charge_log <- log1p(c1)
      X1 <- predict(mh, newdata=nd1, type="lpmatrix")
      X0 <- predict(mh, newdata=nd0, type="lpmatrix")
      dbar <- colMeans(X1 - X0); b <- coef(mh); V <- vcov(mh)
      est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
      pct <- 100*(exp(est)-1); lwr <- 100*(exp(est-1.96*se)-1); upr <- 100*(exp(est+1.96*se)-1)
      eff2x <- bind_rows(eff2x, tibble::tibble(YEAR=as.integer(yy), group=h,
                                               pct2x=pct, lwr=lwr, upr=upr, n=nrow(dH)))
    }
  }
  
  eff2x <- eff2x %>% mutate(YEAR=factor(YEAR, levels=2017:2023),
                            group=factor(group, levels=c("Total", hab_levels)))
  pal <- cols_focus[levels(eff2x$group)]
  p_2x <- ggplot(eff2x, aes(YEAR, pct2x, fill=group)) +
    geom_col(position=position_dodge2(width=0.9, preserve="single"), width=0.85) +
    geom_errorbar(aes(ymin=lwr, ymax=upr),
                  position=position_dodge2(width=0.9, preserve="single"), width=0.22) +
    scale_fill_manual(values=pal, name="Groupe", drop=FALSE) +
    geom_hline(yintercept=0, linewidth=0.5) +
    labs(title="Effet d’un doublement de la charge (médiane → 2× médiane)",
         subtitle="% de variation de la production attendue (GPROD) — IC95%",
         x="Année", y="% pour un doublement de la charge") +
    theme_minimal(base_size=14) +
    theme(panel.grid.minor=element_blank(), legend.position="right")
  
  p_2x
  
  
  
  
  
  # ===== Effet d'un doublement de charge (médiane -> 2× médiane)
  # ===== avec SMOD fixé (global=mediane annuelle / habitat=mediane habitat×année) =====
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  
  # Libellés / couleurs
  hab_labels_map <- c("1"="P. nivales","9"="P. thermiques écorchées",
                      "5"="Queyrellins","6"="Nardaies denses du subalpin")
  hab_levels <- unname(hab_labels_map[c("1","9","5","6")])
  cols_focus <- c("Total"="grey60","P. nivales"="#3B5BDB","P. thermiques écorchées"="#E03131",
                  "Queyrellins"="#F4A261","Nardaies denses du subalpin"="#1B9E77")
  
  # Modèles & data déjà fités
  mglob <- list(`2017`=mod_2017_global, `2018`=mod_2018_global, `2019`=mod_2019_global,
                `2020`=mod_2020_global, `2021`=mod_2021_global, `2022`=mod_2022_global, `2023`=mod_2023_global)
  dglob <- list(`2017`=d2017_all,      `2018`=d2018_all,      `2019`=d2019_all,
                `2020`=d2020_all,      `2021`=d2021_all,      `2022`=d2022_all,      `2023`=d2023_all)
  mhab  <- list(`2017`=mod_2017, `2018`=mod_2018, `2019`=mod_2019,
                `2020`=mod_2020, `2021`=mod_2021, `2022`=mod_2022, `2023`=mod_2023)
  dhab  <- list(`2017`=d2017,    `2018`=d2018,    `2019`=d2019,
                `2020`=d2020,    `2021`=d2021,    `2022`=d2022,    `2023`=d2023)
  
  raw_col <- "chargement_median_2022-2024"   # charge brute (nom avec tirets OK via [[raw_col]])
  
  res2x <- tibble::tibble(
    YEAR=integer(), group=character(),
    pct2x=numeric(), lwr=numeric(), upr=numeric(), n=integer(),
    smod_ref=numeric(), c0=numeric(), c1=numeric()
  )
  
  for (yy in names(mglob)) {
    
    # ---------- GLOBAL (SMOD fixé à la médiane annuelle) ----------
    mg <- mglob[[yy]]; dg <- dglob[[yy]]
    # borne de sécurité sur charge (pour éviter extrapolation délirante)
    lo_all <- quantile(pmax(dg[[raw_col]],0), 0.01, na.rm=TRUE)
    hi_all <- quantile(pmax(dg[[raw_col]],0), 0.99, na.rm=TRUE)
    
    med_all <- as.numeric(median(pmax(dg[[raw_col]],0), na.rm=TRUE))
    c0 <- max(lo_all, min(med_all, hi_all))
    c1 <- max(lo_all, min(2*med_all, hi_all))
    
    # SMOD : référence annuelle puis standardisation comme au fit global
    smod_ref_year <- median(dg$SMOD, na.rm=TRUE)
    muS_g <- mean(dg$SMOD, na.rm=TRUE); sdS_g <- sd(dg$SMOD, na.rm=TRUE)
    smod_std_ref_g <- if (is.finite(sdS_g) && sdS_g>0) (smod_ref_year - muS_g)/sdS_g else 0
    
    nd0 <- dg; nd1 <- dg
    nd0$Charge_log <- log1p(c0); nd1$Charge_log <- log1p(c1)
    nd0$SMOD_2023_std <- smod_std_ref_g; nd1$SMOD_2023_std <- smod_std_ref_g
    
    X1 <- predict(mg, newdata=nd1, type="lpmatrix")
    X0 <- predict(mg, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1 - X0); b <- coef(mg); V <- vcov(mg)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    pct <- 100*(exp(est)-1); lwr <- 100*(exp(est-1.96*se)-1); upr <- 100*(exp(est+1.96*se)-1)
    
    res2x <- bind_rows(res2x, tibble::tibble(
      YEAR=as.integer(yy), group="Total",
      pct2x=pct, lwr=lwr, upr=upr, n=nrow(dg),
      smod_ref=smod_ref_year, c0=c0, c1=c1
    ))
    
    # ---------- HABITATS (SMOD fixé à la médiane habitat×année) ----------
    mh <- mhab[[yy]]; dh <- dhab[[yy]]
    
    # Standardisation utilisée dans le modèle habitat (faite sur dh complet des 4 habitats)
    muS_h <- mean(dh$SMOD, na.rm=TRUE); sdS_h <- sd(dh$SMOD, na.rm=TRUE)
    
    for (h in hab_levels) {
      dH <- dh %>% filter(hab4==h)
      if (nrow(dH) < 30) next
      
      # bornes & médiane de charge pour cet habitat×année
      lo_h <- quantile(pmax(dH[[raw_col]],0), 0.01, na.rm=TRUE)
      hi_h <- quantile(pmax(dH[[raw_col]],0), 0.99, na.rm=TRUE)
      med_h <- as.numeric(median(pmax(dH[[raw_col]],0), na.rm=TRUE))
      c0 <- max(lo_h, min(med_h, hi_h))
      c1 <- max(lo_h, min(2*med_h, hi_h))
      
      # SMOD de référence = médiane habitat×année,
      # transformée avec la même standardisation que celle du fit habitat (mu/sd de dh)
      smod_ref_h <- median(dH$SMOD, na.rm=TRUE)
      smod_std_ref_h <- if (is.finite(sdS_h) && sdS_h>0) (smod_ref_h - muS_h)/sdS_h else 0
      
      nd0 <- dH; nd1 <- dH
      nd0$Charge_log <- log1p(c0); nd1$Charge_log <- log1p(c1)
      nd0$SMOD_2023_std <- smod_std_ref_h; nd1$SMOD_2023_std <- smod_std_ref_h
      
      X1 <- predict(mh, newdata=nd1, type="lpmatrix")
      X0 <- predict(mh, newdata=nd0, type="lpmatrix")
      dbar <- colMeans(X1 - X0); b <- coef(mh); V <- vcov(mh)
      est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
      pct <- 100*(exp(est)-1); lwr <- 100*(exp(est-1.96*se)-1); upr <- 100*(exp(est+1.96*se)-1)
      
      res2x <- bind_rows(res2x, tibble::tibble(
        YEAR=as.integer(yy), group=h,
        pct2x=pct, lwr=lwr, upr=upr, n=nrow(dH),
        smod_ref=smod_ref_h, c0=c0, c1=c1
      ))
    }
  }
  
  # ---- Barplot : 5 barres par année ----
  res2x <- res2x %>%
    mutate(YEAR=factor(YEAR, levels=2017:2023),
           group=factor(group, levels=c("Total", hab_levels)),
           sig = (lwr>0)|(upr<0))
  
  pal <- cols_focus[levels(res2x$group)]
  
  p_2x_fixsmod <- ggplot(res2x, aes(YEAR, pct2x, fill=group)) +
    geom_col(position=position_dodge2(width=0.9, preserve="single"),
             width=0.85, aes(alpha=ifelse(sig, 1, 0.4))) +
    geom_errorbar(aes(ymin=lwr, ymax=upr),
                  position=position_dodge2(width=0.9, preserve="single"),
                  width=0.22) +
    scale_fill_manual(values=pal, name="Groupe", drop=FALSE) +
    scale_alpha_identity() +
    geom_hline(yintercept=0, linewidth=0.5) +
    labs(
      title = "Effet d’un doublement de la charge (médiane → 2× médiane)",
      subtitle = "SMOD fixé : médiane annuelle (Total) / médiane habitat×année (barres habitat). Barres = IC95%.",
      x = "Année", y = "% de variation de la production (GPROD)"
    ) +
    theme_minimal(base_size=14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")
  
  p_2x_fixsmod
  
  # (optionnel) voir les références utilisées :
  # res2x %>% arrange(YEAR, group) %>% select(YEAR, group, smod_ref, c0, c1, n) %>% print(n=Inf)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ========= Préambule commun =========
  library(dplyr)
  library(mgcv)
  library(parallel)
  
  # Threads
  n_th <- max(1, parallel::detectCores() - 1)
  
  # AUCg de secours (si manquant)
  if (!"AUCg" %in% names(dt_legacy) && "GPROD" %in% names(dt_legacy)) {
    dt_legacy$AUCg <- dt_legacy$GPROD
  }
  dt_legacy$AUCg <- pmax(dt_legacy$AUCg, .Machine$double.eps)
  
  # ========= 2017 =========
  data_2017 <- dt_legacy %>%
    filter(YEAR == 2017) %>%
    transmute(
      AUCg = AUCg,
      dah = as.numeric(DAH),
      SMOD_2023 = as.numeric(SMOD),
      alpage = factor(alpage),
      Charge_log = as.numeric(Charge_log),      # log1p(charge) déjà calculé en amont
      charge_bin = factor(charge_bin)           # pour by=charge_bin (séparateur 0/1)
    ) %>% filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023),
                 is.finite(Charge_log))
  
  mod_2017_plus <- bam(
    log(AUCg) ~
      s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
      ti(dah, SMOD_2023, k=c(20,20)) +
      s(alpage, bs="re") +
      charge_bin +
      s(Charge_log, k=10, bs="ts", by=charge_bin) +
      ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin) +
      s(alpage, bs="re", by=charge_bin) +
      ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
    data = data_2017, method = "fREML", select = TRUE,
    discrete = TRUE, nthreads = n_th, na.action = na.exclude
  )
  
  # ========= 2018 =========
  data_2018 <- dt_legacy %>%
    filter(YEAR == 2018) %>%
    transmute(
      AUCg = AUCg,
      dah = as.numeric(DAH),
      SMOD_2023 = as.numeric(SMOD),
      alpage = factor(alpage),
      Charge_log = as.numeric(Charge_log),
      charge_bin = factor(charge_bin)
    ) %>% filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023),
                 is.finite(Charge_log))
  
  mod_2018_plus <- bam(
    log(AUCg) ~
      s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
      ti(dah, SMOD_2023, k=c(20,20)) +
      s(alpage, bs="re") +
      charge_bin +
      s(Charge_log, k=10, bs="ts", by=charge_bin) +
      ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin) +
      s(alpage, bs="re", by=charge_bin) +
      ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
    data = data_2018, method = "fREML", select = TRUE,
    discrete = TRUE, nthreads = n_th, na.action = na.exclude
  )
  
  # ========= 2019 =========
  data_2019 <- dt_legacy %>%
    filter(YEAR == 2019) %>%
    transmute(
      AUCg = AUCg,
      dah = as.numeric(DAH),
      SMOD_2023 = as.numeric(SMOD),
      alpage = factor(alpage),
      Charge_log = as.numeric(Charge_log),
      charge_bin = factor(charge_bin)
    ) %>% filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023),
                 is.finite(Charge_log))
  
  mod_2019_plus <- bam(
    log(AUCg) ~
      s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
      ti(dah, SMOD_2023, k=c(20,20)) +
      s(alpage, bs="re") +
      charge_bin +
      s(Charge_log, k=10, bs="ts", by=charge_bin) +
      ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin) +
      s(alpage, bs="re", by=charge_bin) +
      ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
    data = data_2019, method = "fREML", select = TRUE,
    discrete = TRUE, nthreads = n_th, na.action = na.exclude
  )
  
  # ========= 2020 =========
  data_2020 <- dt_legacy %>%
    filter(YEAR == 2020) %>%
    transmute(
      AUCg = AUCg,
      dah = as.numeric(DAH),
      SMOD_2023 = as.numeric(SMOD),
      alpage = factor(alpage),
      Charge_log = as.numeric(Charge_log),
      charge_bin = factor(charge_bin)
    ) %>% filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023),
                 is.finite(Charge_log))
  
  mod_2020_plus <- bam(
    log(AUCg) ~
      s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
      ti(dah, SMOD_2023, k=c(20,20)) +
      s(alpage, bs="re") +
      charge_bin +
      s(Charge_log, k=10, bs="ts", by=charge_bin) +
      ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin) +
      s(alpage, bs="re", by=charge_bin) +
      ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
    data = data_2020, method = "fREML", select = TRUE,
    discrete = TRUE, nthreads = n_th, na.action = na.exclude
  )
  
  # ========= 2021 =========
  data_2021 <- dt_legacy %>%
    filter(YEAR == 2021) %>%
    transmute(
      AUCg = AUCg,
      dah = as.numeric(DAH),
      SMOD_2023 = as.numeric(SMOD),
      alpage = factor(alpage),
      Charge_log = as.numeric(Charge_log),
      charge_bin = factor(charge_bin)
    ) %>% filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023),
                 is.finite(Charge_log))
  
  mod_2021_plus <- bam(
    log(AUCg) ~
      s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
      ti(dah, SMOD_2023, k=c(20,20)) +
      s(alpage, bs="re") +
      charge_bin +
      s(Charge_log, k=10, bs="ts", by=charge_bin) +
      ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin) +
      s(alpage, bs="re", by=charge_bin) +
      ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
    data = data_2021, method = "fREML", select = TRUE,
    discrete = TRUE, nthreads = n_th, na.action = na.exclude
  )
  
  # ========= 2022 =========
  data_2022 <- dt_legacy %>%
    filter(YEAR == 2022) %>%
    transmute(
      AUCg = AUCg,
      dah = as.numeric(DAH),
      SMOD_2023 = as.numeric(SMOD),
      alpage = factor(alpage),
      Charge_log = as.numeric(Charge_log),
      charge_bin = factor(charge_bin)
    ) %>% filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023),
                 is.finite(Charge_log))
  
  mod_2022_plus <- bam(
    log(AUCg) ~
      s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
      ti(dah, SMOD_2023, k=c(20,20)) +
      s(alpage, bs="re") +
      charge_bin +
      s(Charge_log, k=10, bs="ts", by=charge_bin) +
      ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin) +
      s(alpage, bs="re", by=charge_bin) +
      ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
    data = data_2022, method = "fREML", select = TRUE,
    discrete = TRUE, nthreads = n_th, na.action = na.exclude
  )
  
  # ========= 2023 =========
  data_2023 <- dt_legacy %>%
    filter(YEAR == 2023) %>%
    transmute(
      AUCg = AUCg,
      dah = as.numeric(DAH),
      SMOD_2023 = as.numeric(SMOD),
      alpage = factor(alpage),
      Charge_log = as.numeric(Charge_log),
      charge_bin = factor(charge_bin)
    ) %>% filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023),
                 is.finite(Charge_log))
  
  mod_2023_plus <- bam(
    log(AUCg) ~
      s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
      ti(dah, SMOD_2023, k=c(20,20)) +
      s(alpage, bs="re") +
      charge_bin +
      s(Charge_log, k=10, bs="ts", by=charge_bin) +
      ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin) +
      s(alpage, bs="re", by=charge_bin) +
      ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
    data = data_2023, method = "fREML", select = TRUE,
    discrete = TRUE, nthreads = n_th, na.action = na.exclude
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ==== Effet "×2 robuste" (Q40 -> min(2×, Q90)) par année × alpage ====
  # Hypothèses : les modèles mod_2017_plus ... mod_2023_plus existent
  # et data_2017 ... data_2023 aussi (comme plus haut).
  
  library(dplyr); library(tidyr); library(ggplot2); library(mgcv)
  
  mplus <- list(`2017`=mod_2017_plus, `2018`=mod_2018_plus, `2019`=mod_2019_plus,
                `2020`=mod_2020_plus, `2021`=mod_2021_plus, `2022`=mod_2022_plus, `2023`=mod_2023_plus)
  
  dlist <- list(`2017`=data_2017, `2018`=data_2018, `2019`=data_2019,
                `2020`=data_2020, `2021`=data_2021, `2022`=data_2022, `2023`=data_2023)
  
  raw_col <- "chargement_median_2022-2024"
  
  # Choisis (ou fixe) tes 3 alpages
  if (!exists("alpages_focus")) {
    alpages_focus <- dt_legacy %>%
      dplyr::filter(YEAR %in% 2017:2023) %>%
      dplyr::count(alpage, YEAR) %>%
      dplyr::group_by(alpage) %>% dplyr::summarise(n_years=dplyr::n_distinct(YEAR), n=sum(n), .groups="drop") %>%
      dplyr::arrange(dplyr::desc(n)) %>% dplyr::slice_head(n=3) %>% dplyr::pull(alpage)
  }
  
  # Réglages "robustes"
  p_anchor <- 0.40   # ancre au centre (Q40)
  p_cap    <- 0.90   # plafond (Q90)
  minN     <- 30     # minimum d'observations chargées par barre
  flagN    <- 150    # seuil d'appui pour griser si faible
  
  res_alp <- tibble::tibble(
    YEAR=integer(), alpage=character(),
    pct2x=numeric(), lwr=numeric(), upr=numeric(), n=integer(),
    c0=numeric(), c1=numeric(), mult_eff=numeric(), n_high=numeric()
  )
  
  for (yy in names(mplus)) {
    modY <- mplus[[yy]]; datY <- dlist[[yy]]
    if (is.null(modY) || is.null(datY)) next
    alp_levels <- levels(datY$alpage)
    
    rawY <- dt_legacy %>%
      dplyr::filter(YEAR == as.integer(yy)) %>%
      dplyr::select(alpage, tidyselect::all_of(raw_col), Charge_log, charge_bin, DAH, SMOD) %>%
      dplyr::mutate(alpage = factor(alpage, levels = alp_levels))
    
    for (alp in alpages_focus) {
      rH <- rawY %>% dplyr::filter(alpage == alp, charge_bin == 1)
      dH <- datY  %>% dplyr::filter(alpage == alp, charge_bin == 1)
      if (nrow(dH) < minN || sum(is.finite(rH[[raw_col]])) < minN) next
      
      v  <- pmax(rH[[raw_col]], 0)
      lo <- stats::quantile(v, 0.01, na.rm=TRUE)
      hi <- stats::quantile(v, 0.99, na.rm=TRUE)
      c0 <- as.numeric(stats::quantile(v, p_anchor, na.rm=TRUE))
      cap<- as.numeric(stats::quantile(v, p_cap,    na.rm=TRUE))
      c0 <- max(lo, min(c0, hi))
      c1 <- min(max(lo, 2*c0), cap, hi)
      mult_eff <- c1 / c0
      n_high   <- sum(v >= c1, na.rm=TRUE)
      
      nd0 <- dH; nd1 <- dH
      nd0$Charge_log <- log1p(c0); nd1$Charge_log <- log1p(c1)
      
      X1 <- predict(modY, newdata=nd1, type="lpmatrix")
      X0 <- predict(modY, newdata=nd0, type="lpmatrix")
      dbar <- colMeans(X1 - X0); b <- coef(modY); V <- vcov(modY)
      est <- as.numeric(dbar %*% b)
      se  <- sqrt(drop(dbar %*% V %*% dbar))
      
      pct <- 100*(exp(est)-1)
      lwr <- 100*(exp(est-1.96*se)-1)
      upr <- 100*(exp(est+1.96*se)-1)
      
      res_alp <- dplyr::bind_rows(res_alp, tibble::tibble(
        YEAR=as.integer(yy), alpage=as.character(alp),
        pct2x=pct, lwr=lwr, upr=upr, n=nrow(dH),
        c0=c0, c1=c1, mult_eff=mult_eff, n_high=n_high
      ))
    }
  }
  
  # Plot (alpha réduit si appui faible au-dessus de c1)
  res_alp <- res_alp %>%
    dplyr::mutate(
      YEAR=factor(YEAR, levels=2017:2023),
      alpage=factor(alpage, levels=alpages_focus),
      sig=(lwr>0)|(upr<0),
      alpha_bar = ifelse(n_high < flagN, 0.45, 1)  # peu de pixels dans la zone haute -> griser
    )
  
  pal_alp <- setNames(RColorBrewer::brewer.pal(length(alpages_focus), "Set2"), alpages_focus)
  
  p_2x_robuste <- ggplot(res_alp, aes(YEAR, pct2x, fill=alpage)) +
    geom_col(position=position_dodge2(width=0.9, preserve="single"),
             width=0.85, aes(alpha=alpha_bar)) +
    geom_errorbar(aes(ymin=lwr, ymax=upr),
                  position=position_dodge2(width=0.9, preserve="single"), width=0.22) +
    scale_fill_manual(values=pal_alp, name="Alpage", drop=FALSE) +
    scale_alpha_identity() +
    geom_hline(yintercept=0, linewidth=0.5) +
    labs(
      title="Effet '×2' ROBUSTE de la charge : Q40 → min(2×, Q90)",
      subtitle="Sans fixer SMOD/DAH. Barres = IC95%. Transparence si peu de pixels ≥ c1.",
      x="Année", y="% de variation de la production (GPROD)"
    ) +
    theme_minimal(base_size=14) +
    theme(panel.grid.minor=element_blank(), legend.position="right")
  
  p_2x_robuste
  
  # (facultatif) tableau des multiplicateurs effectifs et appuis :
  # res_alp %>% dplyr::arrange(YEAR, alpage) %>% dplyr::select(YEAR, alpage, mult_eff, n_high, n, c0, c1) %>% print(n=Inf)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ===== Effet +10% de charge, SMOD restreint à la classe Mid [121,150], par Année × Alpage =====
  library(dplyr); library(ggplot2); library(mgcv)
  
  # Modèles & jeux annuels (issus des fits précédents)
  mplus <- list(`2017`=mod_2017_plus, `2018`=mod_2018_plus, `2019`=mod_2019_plus,
                `2020`=mod_2020_plus, `2021`=mod_2021_plus, `2022`=mod_2022_plus, `2023`=mod_2023_plus)
  
  dlist <- list(`2017`=data_2017, `2018`=data_2018, `2019`=data_2019,
                `2020`=data_2020, `2021`=data_2021, `2022`=data_2022, `2023`=data_2023)
  
  # Alpages à afficher (si tu veux forcer : alpages_focus <- c("Cayolle","Sanguiniere","Viso"))
  if (!exists("alpages_focus")) {
    alpages_focus <- dt_legacy %>%
      dplyr::filter(YEAR %in% 2017:2023, SMOD >= 121, SMOD <= 150) %>%
      dplyr::count(alpage) %>% dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::slice_head(n=3) %>% dplyr::pull(alpage)
  }
  
  # Seuil de taille d'échantillon par barre (éviter des IC absurdes)
  minN <- 30
  
  res_mid10 <- tibble::tibble(
    YEAR=integer(), alpage=character(),
    pct10=numeric(), lwr=numeric(), upr=numeric(), n=integer()
  )
  
  for (yy in names(mplus)) {
    modY <- mplus[[yy]]; datY <- dlist[[yy]]
    if (is.null(modY) || is.null(datY)) next
    
    for (alp in alpages_focus) {
      # Pixels CHARGÉS et SMOD dans la classe MID (121–150)
      dH <- datY %>%
        dplyr::filter(alpage == alp, charge_bin == 1, SMOD_2023 >= 121, SMOD_2023 <= 150)
      
      if (nrow(dH) < minN) next
      
      # +10% SUR LA CHARGE BRUTE (on repart de log1p -> inverse = exp(.)-1)
      ch0 <- pmax(exp(dH$Charge_log) - 1, 0)
      ch1 <- ch0 * 1.10
      
      nd0 <- dH; nd1 <- dH
      nd0$Charge_log <- log1p(ch0)
      nd1$Charge_log <- log1p(ch1)
      # On laisse SMOD_2023 AUX VALEURS OBSERVÉES DANS [121,150]; DAH inchangé.
      
      # Δη moyen + IC95% (delta method avec V(β))
      X1 <- predict(modY, newdata = nd1, type = "lpmatrix")
      X0 <- predict(modY, newdata = nd0, type = "lpmatrix")
      dbar <- colMeans(X1 - X0); b <- coef(modY); V <- vcov(modY)
      est <- as.numeric(dbar %*% b)
      se  <- sqrt(drop(dbar %*% V %*% dbar))
      
      # Retour en % sur GPROD (réponse = log(AUCg))
      pct <- 100*(exp(est) - 1)
      lwr <- 100*(exp(est - 1.96*se) - 1)
      upr <- 100*(exp(est + 1.96*se) - 1)
      
      res_mid10 <- dplyr::bind_rows(res_mid10, tibble::tibble(
        YEAR = as.integer(yy), alpage = as.character(alp),
        pct10 = pct, lwr = lwr, upr = upr, n = nrow(dH)
      ))
    }
  }
  
  # ----- Plot : 3 barres par année (effet +10% ; SMOD ∈ [121,150]) -----
  res_mid10 <- res_mid10 %>%
    dplyr::mutate(
      YEAR   = factor(YEAR, levels = 2017:2023),
      alpage = factor(alpage, levels = alpages_focus),
      sig    = (lwr > 0) | (upr < 0)
    )
  
  # Palette simple pour 3 alpages
  pal_alp <- setNames(c("#66C2A5","#FC8D62","#8DA0CB")[seq_along(alpages_focus)], alpages_focus)
  
  p_mid10_class <- ggplot(res_mid10, aes(YEAR, pct10, fill=alpage)) +
    geom_col(position=position_dodge2(width=0.9, preserve="single"),
             width=0.85, aes(alpha=ifelse(sig, 1, 0.45))) +
    geom_errorbar(aes(ymin=lwr, ymax=upr),
                  position=position_dodge2(width=0.9, preserve="single"),
                  width=0.22) +
    scale_fill_manual(values=pal_alp, name="Alpage", drop=FALSE) +
    scale_alpha_identity() +
    geom_hline(yintercept=0, linewidth=0.5) +
    labs(
      title   = "Élasticité +10% de charge (SMOD ∈ [121,150])",
      subtitle= "Gradient uniquement (charge_bin=1). Barres = IC95%.",
      x="Année", y="% de variation de la production (GPROD) pour +10% de charge"
    ) +
    theme_minimal(base_size=14) +
    theme(panel.grid.minor=element_blank(),
          legend.position="right")
  
  p_mid10_class
  
  # (optionnel) voir les effectifs par barre :
  # res_mid10 %>% dplyr::arrange(YEAR, alpage) %>% dplyr::select(YEAR, alpage, n, pct10, lwr, upr) %>% print(n=Inf)
    
}

































# =========================================================
# 0) Préambule (tu as déjà ces 7 modèles "mod_YYYY_plus")
# =========================================================
library(dplyr); library(mgcv); library(ggplot2); library(parallel)

n_th <- max(1, parallel::detectCores() - 1)

if (!"AUCg" %in% names(dt_legacy) && "GPROD" %in% names(dt_legacy)) {
  dt_legacy$AUCg <- dt_legacy$GPROD
}
dt_legacy$AUCg <- pmax(dt_legacy$AUCg, .Machine$double.eps)

hab_labels_map <- c(
  "1" = "P. nivales",
  "9" = "P. thermiques écorchées",
  "5" = "Queyrellins",
  "6" = "Nardaies denses du subalpin"
)

# =========================================================
# 1) Modèles PAR HABITAT (4 classes) — un par année
#    (les effets liés à Charge_log varient selon hab4 et charge_bin)
# =========================================================

# ---------- 2017 ----------
d2017_hab <- dt_legacy %>%
  dplyr::filter(YEAR == 2017, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(
    AUCg = AUCg,
    dah  = as.numeric(DAH),
    SMOD_2023 = as.numeric(SMOD),
    alpage = factor(alpage),
    Charge_log = as.numeric(Charge_log),
    charge_bin = factor(charge_bin),
    hab4 = factor(as.character(habitat_code),
                  levels = c("1","9","5","6"),
                  labels = unname(hab_labels_map[c("1","9","5","6")]))
  ) %>% dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2017_hab_plus <- bam(
  log(AUCg) ~
    s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
    ti(dah, SMOD_2023, k=c(20,20)) +
    hab4 + s(alpage, bs="re") +
    charge_bin +
    s(Charge_log, k=10, bs="ts", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, SMOD_2023, k=c(15,15), by=interaction(hab4, charge_bin)) +
    s(alpage, bs="re", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, alpage, bs=c("tp","re"), by=interaction(hab4, charge_bin), k=c(6, NA)),
  data = d2017_hab, method="fREML", select=TRUE,
  discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# ---------- 2018 ----------
d2018_hab <- dt_legacy %>%
  dplyr::filter(YEAR == 2018, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(
    AUCg = AUCg, dah = as.numeric(DAH), SMOD_2023 = as.numeric(SMOD),
    alpage = factor(alpage), Charge_log = as.numeric(Charge_log),
    charge_bin = factor(charge_bin),
    hab4 = factor(as.character(habitat_code),
                  levels=c("1","9","5","6"),
                  labels=unname(hab_labels_map[c("1","9","5","6")]))
  ) %>% dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2018_hab_plus <- bam(
  log(AUCg) ~
    s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
    ti(dah, SMOD_2023, k=c(20,20)) +
    hab4 + s(alpage, bs="re") +
    charge_bin +
    s(Charge_log, k=10, bs="ts", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, SMOD_2023, k=c(15,15), by=interaction(hab4, charge_bin)) +
    s(alpage, bs="re", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, alpage, bs=c("tp","re"), by=interaction(hab4, charge_bin), k=c(6, NA)),
  data = d2018_hab, method="fREML", select=TRUE,
  discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# ---------- 2019 ----------
d2019_hab <- dt_legacy %>%
  dplyr::filter(YEAR == 2019, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(
    AUCg = AUCg, dah = as.numeric(DAH), SMOD_2023 = as.numeric(SMOD),
    alpage = factor(alpage), Charge_log = as.numeric(Charge_log),
    charge_bin = factor(charge_bin),
    hab4 = factor(as.character(habitat_code),
                  levels=c("1","9","5","6"),
                  labels=unname(hab_labels_map[c("1","9","5","6")]))
  ) %>% dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2019_hab_plus <- bam(
  log(AUCg) ~
    s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
    ti(dah, SMOD_2023, k=c(20,20)) +
    hab4 + s(alpage, bs="re") +
    charge_bin +
    s(Charge_log, k=10, bs="ts", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, SMOD_2023, k=c(15,15), by=interaction(hab4, charge_bin)) +
    s(alpage, bs="re", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, alpage, bs=c("tp","re"), by=interaction(hab4, charge_bin), k=c(6, NA)),
  data = d2019_hab, method="fREML", select=TRUE,
  discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# ---------- 2020 ----------
d2020_hab <- dt_legacy %>%
  dplyr::filter(YEAR == 2020, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(
    AUCg = AUCg, dah = as.numeric(DAH), SMOD_2023 = as.numeric(SMOD),
    alpage = factor(alpage), Charge_log = as.numeric(Charge_log),
    charge_bin = factor(charge_bin),
    hab4 = factor(as.character(habitat_code),
                  levels=c("1","9","5","6"),
                  labels=unname(hab_labels_map[c("1","9","5","6")]))
  ) %>% dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2020_hab_plus <- bam(
  log(AUCg) ~
    s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
    ti(dah, SMOD_2023, k=c(20,20)) +
    hab4 + s(alpage, bs="re") +
    charge_bin +
    s(Charge_log, k=10, bs="ts", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, SMOD_2023, k=c(15,15), by=interaction(hab4, charge_bin)) +
    s(alpage, bs="re", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, alpage, bs=c("tp","re"), by=interaction(hab4, charge_bin), k=c(6, NA)),
  data = d2020_hab, method="fREML", select=TRUE,
  discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# ---------- 2021 ----------
d2021_hab <- dt_legacy %>%
  dplyr::filter(YEAR == 2021, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(
    AUCg = AUCg, dah = as.numeric(DAH), SMOD_2023 = as.numeric(SMOD),
    alpage = factor(alpage), Charge_log = as.numeric(Charge_log),
    charge_bin = factor(charge_bin),
    hab4 = factor(as.character(habitat_code),
                  levels=c("1","9","5","6"),
                  labels=unname(hab_labels_map[c("1","9","5","6")]))
  ) %>% dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2021_hab_plus <- bam(
  log(AUCg) ~
    s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
    ti(dah, SMOD_2023, k=c(20,20)) +
    hab4 + s(alpage, bs="re") +
    charge_bin +
    s(Charge_log, k=10, bs="ts", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, SMOD_2023, k=c(15,15), by=interaction(hab4, charge_bin)) +
    s(alpage, bs="re", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, alpage, bs=c("tp","re"), by=interaction(hab4, charge_bin), k=c(6, NA)),
  data = d2021_hab, method="fREML", select=TRUE,
  discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# ---------- 2022 ----------
d2022_hab <- dt_legacy %>%
  dplyr::filter(YEAR == 2022, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(
    AUCg = AUCg, dah = as.numeric(DAH), SMOD_2023 = as.numeric(SMOD),
    alpage = factor(alpage), Charge_log = as.numeric(Charge_log),
    charge_bin = factor(charge_bin),
    hab4 = factor(as.character(habitat_code),
                  levels=c("1","9","5","6"),
                  labels=unname(hab_labels_map[c("1","9","5","6")]))
  ) %>% dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2022_hab_plus <- bam(
  log(AUCg) ~
    s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
    ti(dah, SMOD_2023, k=c(20,20)) +
    hab4 + s(alpage, bs="re") +
    charge_bin +
    s(Charge_log, k=10, bs="ts", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, SMOD_2023, k=c(15,15), by=interaction(hab4, charge_bin)) +
    s(alpage, bs="re", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, alpage, bs=c("tp","re"), by=interaction(hab4, charge_bin), k=c(6, NA)),
  data = d2022_hab, method="fREML", select=TRUE,
  discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# ---------- 2023 ----------
d2023_hab <- dt_legacy %>%
  dplyr::filter(YEAR == 2023, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(
    AUCg = AUCg, dah = as.numeric(DAH), SMOD_2023 = as.numeric(SMOD),
    alpage = factor(alpage), Charge_log = as.numeric(Charge_log),
    charge_bin = factor(charge_bin),
    hab4 = factor(as.character(habitat_code),
                  levels=c("1","9","5","6"),
                  labels=unname(hab_labels_map[c("1","9","5","6")]))
  ) %>% dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2023_hab_plus <- bam(
  log(AUCg) ~
    s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
    ti(dah, SMOD_2023, k=c(20,20)) +
    hab4 + s(alpage, bs="re") +
    charge_bin +
    s(Charge_log, k=10, bs="ts", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, SMOD_2023, k=c(15,15), by=interaction(hab4, charge_bin)) +
    s(alpage, bs="re", by=interaction(hab4, charge_bin)) +
    ti(Charge_log, alpage, bs=c("tp","re"), by=interaction(hab4, charge_bin), k=c(6, NA)),
  data = d2023_hab, method="fREML", select=TRUE,
  discrete=TRUE, nthreads=n_th, na.action=na.exclude
)




# Registres pour accès facile
mplus_global <- list(`2017`=mod_2017_plus, `2018`=mod_2018_plus, `2019`=mod_2019_plus,
                     `2020`=mod_2020_plus, `2021`=mod_2021_plus, `2022`=mod_2022_plus, `2023`=mod_2023_plus)

d_global <- list(`2017`=data_2017, `2018`=data_2018, `2019`=data_2019,
                 `2020`=data_2020, `2021`=data_2021, `2022`=data_2022, `2023`=data_2023)

mplus_hab <- list(`2017`=mod_2017_hab_plus, `2018`=mod_2018_hab_plus, `2019`=mod_2019_hab_plus,
                  `2020`=mod_2020_hab_plus, `2021`=mod_2021_hab_plus, `2022`=mod_2022_hab_plus, `2023`=mod_2023_hab_plus)

d_hab <- list(`2017`=d2017_hab, `2018`=d2018_hab, `2019`=d2019_hab,
              `2020`=d2020_hab, `2021`=d2021_hab, `2022`=d2022_hab, `2023`=d2023_hab)

# ===== Effet +10% de charge, SANS fixer SMOD (marginalisé), par Année × Alpage × Groupe =====
library(dplyr); library(ggplot2); library(mgcv)

# Registres (venant du bloc modèles que tu viens d’exécuter)
# mplus_global : mod_YYYY_plus ; d_global : data_YYYY
# mplus_hab    : mod_YYYY_hab_plus ; d_hab : dYYYY_hab
stopifnot(exists("mplus_global"), exists("d_global"), exists("mplus_hab"), exists("d_hab"))

# Groupes d'habitats (labels)
hab_labels_map <- c(
  "1" = "P. nivales",
  "9" = "P. thermiques écorchées",
  "5" = "Queyrellins",
  "6" = "Nardaies denses du subalpin"
)
hab_levels <- unname(hab_labels_map[c("1","9","5","6")])

# Alpages à afficher (ordre souhaité) — ajuste si besoin
alpages_focus <- c("Viso","Cayolle","Sanguiniere")

minN <- 20  # seuil minimum d'observations par barre

res <- tibble::tibble(
  YEAR=integer(), alpage=character(), groupe=character(),
  pct10=numeric(), lwr=numeric(), upr=numeric(), n=integer()
)

for (yy in names(mplus_global)) {
  # ----- TOTAL (tous habitats; modèle global) -----
  modG <- mplus_global[[yy]]; datG <- d_global[[yy]]
  for (alp in alpages_focus) {
    dG <- datG %>%
      dplyr::filter(alpage == alp, charge_bin == 1)   # SMOD/DAH laissés tels qu’observés
    if (nrow(dG) < minN) next
    
    ch0 <- pmax(exp(dG$Charge_log) - 1, 0)  # inverse de log1p
    ch1 <- ch0 * 1.10
    
    nd0 <- dG; nd1 <- dG
    nd0$Charge_log <- log1p(ch0)
    nd1$Charge_log <- log1p(ch1)
    
    X1 <- predict(modG, newdata=nd1, type="lpmatrix")
    X0 <- predict(modG, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1 - X0); b <- coef(modG); V <- vcov(modG)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    
    res <- dplyr::bind_rows(res, tibble::tibble(
      YEAR=as.integer(yy), alpage=alp, groupe="Total",
      pct10=100*(exp(est)-1),
      lwr =100*(exp(est-1.96*se)-1),
      upr =100*(exp(est+1.96*se)-1),
      n=nrow(dG)
    ))
  }
  
  # ----- PAR HABITAT (modèle hab4) -----
  modH <- mplus_hab[[yy]]; datH <- d_hab[[yy]]
  for (alp in alpages_focus) {
    for (lab in hab_levels) {
      dH <- datH %>%
        dplyr::filter(alpage == alp, hab4 == lab, charge_bin == 1)
      if (nrow(dH) < minN) next
      
      ch0 <- pmax(exp(dH$Charge_log) - 1, 0)
      ch1 <- ch0 * 1.10
      nd0 <- dH; nd1 <- dH
      nd0$Charge_log <- log1p(ch0)
      nd1$Charge_log <- log1p(ch1)
      
      X1 <- predict(modH, newdata=nd1, type="lpmatrix")
      X0 <- predict(modH, newdata=nd0, type="lpmatrix")
      dbar <- colMeans(X1 - X0); b <- coef(modH); V <- vcov(modH)
      est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
      
      res <- dplyr::bind_rows(res, tibble::tibble(
        YEAR=as.integer(yy), alpage=alp, groupe=lab,
        pct10=100*(exp(est)-1),
        lwr =100*(exp(est-1.96*se)-1),
        upr =100*(exp(est+1.96*se)-1),
        n=nrow(dH)
      ))
    }
  }
}

# ----- Plot facetté (par alpage, en lignes) -----
res_plot <- res %>%
  dplyr::mutate(
    YEAR  = factor(YEAR, levels=2017:2023),
    alpage= factor(alpage, levels=alpages_focus),
    groupe= factor(groupe, levels=c("Total", hab_levels)),
    sig   = (lwr > 0) | (upr < 0)
  )

cols <- c(
  "Total"                        = "grey60",
  "P. nivales"                   = "#3B5BDB",
  "P. thermiques écorchées"      = "#E03131",
  "Queyrellins"                  = "#F4A261",
  "Nardaies denses du subalpin"  = "#1B9E77"
)

p_final_noFix <- ggplot(res_plot, aes(x=YEAR, y=pct10, fill=groupe)) +
  geom_col(position=position_dodge2(width=0.9, preserve="single"),
           width=0.86, aes(alpha=ifelse(sig,1,0.45))) +
  geom_errorbar(aes(ymin=lwr, ymax=upr),
                position=position_dodge2(width=0.9, preserve="single"),
                width=0.22) +
  scale_fill_manual(values=cols, name=NULL, drop=FALSE) +
  scale_alpha_identity() +
  geom_hline(yintercept=0, linewidth=0.5) +
  labs(
    title   = "Élasticité +10% de charge — SMOD non fixé (marginalisé)",
    subtitle= "Par année × alpage. Groupes = Total + 4 habitats. Barres = IC95%; transparence = non significatif.",
    x="Année", y="% variation de GPROD pour +10% de charge"
  ) +
  facet_grid(rows = vars(alpage)) +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right")

p_final_noFix

# (optionnel) inspecter les effectifs par barre
# res_plot %>% arrange(alpage, YEAR, groupe) %>% select(alpage, YEAR, groupe, n) %>% print(n=Inf)










# ================================
# Commun
# ================================
library(dplyr)
library(ggplot2)
library(mgcv)
library(tibble)

stopifnot(exists("mplus_global"), exists("d_global"),
          exists("mplus_hab"),    exists("d_hab"))

hab_labels_map <- c(
  "1" = "P. nivales",
  "9" = "P. thermiques écorchées",
  "5" = "Queyrellins",
  "6" = "Nardaies denses du subalpin"
)
hab_levels <- unname(hab_labels_map[c("1","9","5","6")])

cols <- c(
  "Total"                        = "grey60",
  "P. nivales"                   = "#3B5BDB",
  "P. thermiques écorchées"      = "#E03131",
  "Queyrellins"                  = "#F4A261",
  "Nardaies denses du subalpin"  = "#1B9E77"
)

minN <- 20  # minimum d'observations par barre
years <- names(mplus_global)  # "2017"..."2023"

# ================================
# GRAPH 1 : 3 alpages ensemble
# ================================
alpages_all3 <- c("Viso","Cayolle","Sanguiniere")

res_all3 <- tibble(YEAR=integer(), groupe=character(),
                   pct10=numeric(), lwr=numeric(), upr=numeric(), n=integer())

for (yy in years) {
  ## Total (modèle global)
  modG <- mplus_global[[yy]]; datG <- d_global[[yy]]
  dG <- datG %>% dplyr::filter(alpage %in% alpages_all3, charge_bin == 1)
  if (nrow(dG) >= minN) {
    ch0 <- pmax(exp(dG$Charge_log) - 1, 0); ch1 <- ch0 * 1.10
    nd0 <- dG; nd1 <- dG; nd0$Charge_log <- log1p(ch0); nd1$Charge_log <- log1p(ch1)
    X1 <- predict(modG, newdata=nd1, type="lpmatrix")
    X0 <- predict(modG, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1 - X0); b <- coef(modG); V <- vcov(modG)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    res_all3 <- dplyr::bind_rows(res_all3, tibble(
      YEAR=as.integer(yy), groupe="Total",
      pct10=100*(exp(est)-1), lwr=100*(exp(est-1.96*se)-1),
      upr=100*(exp(est+1.96*se)-1), n=nrow(dG)
    ))
  }
  
  ## 4 habitats (modèle par habitat)
  modH <- mplus_hab[[yy]]; datH <- d_hab[[yy]]
  for (lab in hab_levels) {
    dH <- datH %>% dplyr::filter(alpage %in% alpages_all3, hab4 == lab, charge_bin == 1)
    if (nrow(dH) < minN) next
    ch0 <- pmax(exp(dH$Charge_log) - 1, 0); ch1 <- ch0 * 1.10
    nd0 <- dH; nd1 <- dH; nd0$Charge_log <- log1p(ch0); nd1$Charge_log <- log1p(ch1)
    X1 <- predict(modH, newdata=nd1, type="lpmatrix")
    X0 <- predict(modH, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1 - X0); b <- coef(modH); V <- vcov(modH)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    res_all3 <- dplyr::bind_rows(res_all3, tibble(
      YEAR=as.integer(yy), groupe=lab,
      pct10=100*(exp(est)-1), lwr=100*(exp(est-1.96*se)-1),
      upr=100*(exp(est+1.96*se)-1), n=nrow(dH)
    ))
  }
}

res_all3 <- res_all3 %>%
  dplyr::mutate(
    YEAR  = factor(YEAR, levels = 2017:2023),
    groupe= factor(groupe, levels = c("Total", hab_levels)),
    sig   = (lwr > 0) | (upr < 0)
  )

p_all3 <- ggplot(res_all3, aes(YEAR, pct10, fill=groupe)) +
  geom_col(position=position_dodge2(width=0.9, preserve="single"),
           width=0.86, aes(alpha=ifelse(sig,1,0.45))) +
  geom_errorbar(aes(ymin=lwr, ymax=upr),
                position=position_dodge2(width=0.9, preserve="single"),
                width=0.22) +
  scale_fill_manual(values=cols, name=NULL, drop=FALSE) +
  scale_alpha_identity() +
  geom_hline(yintercept=0, linewidth=0.5) +
  labs(
    title   = "Élasticité +10% de charge — 3 alpages ensemble (SMOD non fixé)",
    subtitle= "Groupes = Total + 4 habitats. Barres = IC95%; transparence = non significatif.",
    x="Année", y="% variation de GPROD pour +10% de charge"
  ) +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(), legend.position="right")

p_all3
# ggsave("bar_all3.png", p_all3, width=10, height=6, dpi=300)

# ================================
# GRAPH 2 : Cayolle + Sanguiniere (sans Viso)
# ================================
alpages_2 <- c("Cayolle","Sanguiniere")

res_2 <- tibble(YEAR=integer(), groupe=character(),
                pct10=numeric(), lwr=numeric(), upr=numeric(), n=integer())

for (yy in years) {
  ## Total
  modG <- mplus_global[[yy]]; datG <- d_global[[yy]]
  dG <- datG %>% dplyr::filter(alpage %in% alpages_2, charge_bin == 1)
  if (nrow(dG) >= minN) {
    ch0 <- pmax(exp(dG$Charge_log) - 1, 0); ch1 <- ch0 * 1.10
    nd0 <- dG; nd1 <- dG; nd0$Charge_log <- log1p(ch0); nd1$Charge_log <- log1p(ch1)
    X1 <- predict(modG, newdata=nd1, type="lpmatrix")
    X0 <- predict(modG, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1 - X0); b <- coef(modG); V <- vcov(modG)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    res_2 <- dplyr::bind_rows(res_2, tibble(
      YEAR=as.integer(yy), groupe="Total",
      pct10=100*(exp(est)-1), lwr=100*(exp(est-1.96*se)-1),
      upr=100*(exp(est+1.96*se)-1), n=nrow(dG)
    ))
  }
  
  ## 4 habitats
  modH <- mplus_hab[[yy]]; datH <- d_hab[[yy]]
  for (lab in hab_levels) {
    dH <- datH %>% dplyr::filter(alpage %in% alpages_2, hab4 == lab, charge_bin == 1)
    if (nrow(dH) < minN) next
    ch0 <- pmax(exp(dH$Charge_log) - 1, 0); ch1 <- ch0 * 1.10
    nd0 <- dH; nd1 <- dH; nd0$Charge_log <- log1p(ch0); nd1$Charge_log <- log1p(ch1)
    X1 <- predict(modH, newdata=nd1, type="lpmatrix")
    X0 <- predict(modH, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1 - X0); b <- coef(modH); V <- vcov(modH)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    res_2 <- dplyr::bind_rows(res_2, tibble(
      YEAR=as.integer(yy), groupe=lab,
      pct10=100*(exp(est)-1), lwr=100*(exp(est-1.96*se)-1),
      upr=100*(exp(est+1.96*se)-1), n=nrow(dH)
    ))
  }
}

res_2 <- res_2 %>%
  dplyr::mutate(
    YEAR  = factor(YEAR, levels = 2017:2023),
    groupe= factor(groupe, levels = c("Total", hab_levels)),
    sig   = (lwr > 0) | (upr < 0)
  )

p_two <- ggplot(res_2, aes(YEAR, pct10, fill=groupe)) +
  geom_col(position=position_dodge2(width=0.9, preserve="single"),
           width=0.86, aes(alpha=ifelse(sig,1,0.45))) +
  geom_errorbar(aes(ymin=lwr, ymax=upr),
                position=position_dodge2(width=0.9, preserve="single"),
                width=0.22) +
  scale_fill_manual(values=cols, name=NULL, drop=FALSE) +
  scale_alpha_identity() +
  geom_hline(yintercept=0, linewidth=0.5) +
  labs(
    title   = "Élasticité +10% de charge — Cayolle + Sanguiniere (SMOD non fixé)",
    subtitle= "Groupes = Total + 4 habitats. Barres = IC95%; transparence = non significatif.",
    x="Année", y="% variation de GPROD pour +10% de charge"
  ) +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(), legend.position="right")

p_two
# ggsave("bar_cay_sang.png", p_two, width=10, height=6, dpi=300)

    














































































# =========================
# 0) Préambule & réglages
# =========================
library(dplyr); library(ggplot2); library(mgcv); library(parallel); library(tibble)

n_th <- max(1, parallel::detectCores() - 1)

# AUCg de secours
if (!"AUCg" %in% names(dt_legacy) && "GPROD" %in% names(dt_legacy)) {
  dt_legacy$AUCg <- dt_legacy$GPROD
}
dt_legacy$AUCg <- pmax(dt_legacy$AUCg, .Machine$double.eps)

# Alpages ciblés
alpages_all3  <- c("Viso","Cayolle","Sanguiniere")
alpages_noViso<- c("Cayolle","Sanguiniere")

# 4 habitats
hab_labels_map <- c("1"="P. nivales","9"="P. thermiques écorchées","5"="Queyrellins","6"="Nardaies denses du subalpin")
hab_levels <- unname(hab_labels_map[c("1","9","5","6")])

# =========================
# 1) Jeux & MODÈLES / ANNÉE
#    (4 modèles par année)
# =========================

# ---------- 2017 ----------
data_2017_all3 <- dt_legacy %>%
  dplyr::filter(YEAR==2017, alpage %in% alpages_all3) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log),
                   charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2017_ref_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") +
    charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2017_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

d2017_hab_all3 <- dt_legacy %>%
  dplyr::filter(YEAR==2017, alpage %in% alpages_all3, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log),
                   charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code),
                               levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2017_hab_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") +
    charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2017_hab_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

data_2017_noViso <- dt_legacy %>%
  dplyr::filter(YEAR==2017, alpage %in% alpages_noViso) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log),
                   charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2017_ref_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") +
    charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2017_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

d2017_hab_noViso <- dt_legacy %>%
  dplyr::filter(YEAR==2017, alpage %in% alpages_noViso, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log),
                   charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code),
                               levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2017_hab_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") +
    charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2017_hab_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# -------- 2018 --------
data_2018_all3 <- dt_legacy %>% dplyr::filter(YEAR==2018, alpage %in% alpages_all3) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log),
                   charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2018_ref_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") +
    charge_bin + s(Charge_log,k=10,bs="ts",by=charge_bin) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2018_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

d2018_hab_all3 <- dt_legacy %>% dplyr::filter(YEAR==2018, alpage %in% alpages_all3, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log),
                   charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2018_hab_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") +
    charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2018_hab_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

data_2018_noViso <- dt_legacy %>% dplyr::filter(YEAR==2018, alpage %in% alpages_noViso) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log),
                   charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2018_ref_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") +
    charge_bin + s(Charge_log,k=10,bs="ts",by=charge_bin) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2018_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

d2018_hab_noViso <- dt_legacy %>% dplyr::filter(YEAR==2018, alpage %in% alpages_noViso, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log),
                   charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))

mod_2018_hab_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") +
    charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2018_hab_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# -------- 2019 --------
data_2019_all3 <- dt_legacy %>% dplyr::filter(YEAR==2019, alpage %in% alpages_all3) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log),
                   charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2019_ref_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") +
    charge_bin + s(Charge_log,k=10,bs="ts",by=charge_bin) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2019_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
d2019_hab_all3 <- dt_legacy %>% dplyr::filter(YEAR==2019, alpage %in% alpages_all3, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2019_hab_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") +
    charge_bin + s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2019_hab_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
data_2019_noViso <- dt_legacy %>% dplyr::filter(YEAR==2019, alpage %in% alpages_noViso) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2019_ref_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) + ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2019_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
d2019_hab_noViso <- dt_legacy %>% dplyr::filter(YEAR==2019, alpage %in% alpages_noViso, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2019_hab_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2019_hab_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# -------- 2020 --------
data_2020_all3 <- dt_legacy %>% dplyr::filter(YEAR==2020, alpage %in% alpages_all3) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2020_ref_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) + ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2020_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
d2020_hab_all3 <- dt_legacy %>% dplyr::filter(YEAR==2020, alpage %in% alpages_all3, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2020_hab_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2020_hab_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
data_2020_noViso <- dt_legacy %>% dplyr::filter(YEAR==2020, alpage %in% alpages_noViso) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2020_ref_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) + ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2020_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
d2020_hab_noViso <- dt_legacy %>% dplyr::filter(YEAR==2020, alpage %in% alpages_noViso, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2020_hab_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2020_hab_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# -------- 2021 --------
data_2021_all3 <- dt_legacy %>% dplyr::filter(YEAR==2021, alpage %in% alpages_all3) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD),
                   alpage=factor(alpage), Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2021_ref_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) + ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2021_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
d2021_hab_all3 <- dt_legacy %>% dplyr::filter(YEAR==2021, alpage %in% alpages_all3, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2021_hab_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2021_hab_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
data_2021_noViso <- dt_legacy %>% dplyr::filter(YEAR==2021, alpage %in% alpages_noViso) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2021_ref_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) + ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2021_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
d2021_hab_noViso <- dt_legacy %>% dplyr::filter(YEAR==2021, alpage %in% alpages_noViso, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2021_hab_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2021_hab_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# -------- 2022 --------
data_2022_all3 <- dt_legacy %>% dplyr::filter(YEAR==2022, alpage %in% alpages_all3) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2022_ref_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) + ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2022_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
d2022_hab_all3 <- dt_legacy %>% dplyr::filter(YEAR==2022, alpage %in% alpages_all3, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2022_hab_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2022_hab_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
data_2022_noViso <- dt_legacy %>% dplyr::filter(YEAR==2022, alpage %in% alpages_noViso) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2022_ref_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) + ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2022_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
d2022_hab_noViso <- dt_legacy %>% dplyr::filter(YEAR==2022, alpage %in% alpages_noViso, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2022_hab_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2022_hab_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# -------- 2023 --------
data_2023_all3 <- dt_legacy %>% dplyr::filter(YEAR==2023, alpage %in% alpages_all3) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2023_ref_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) + ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2023_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
d2023_hab_all3 <- dt_legacy %>% dplyr::filter(YEAR==2023, alpage %in% alpages_all3, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2023_hab_all3 <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2023_hab_all3, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
data_2023_noViso <- dt_legacy %>% dplyr::filter(YEAR==2023, alpage %in% alpages_noViso) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2023_ref_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=charge_bin) + ti(Charge_log,SMOD_2023,k=c(15,15),by=charge_bin),
  data=data_2023_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)
d2023_hab_noViso <- dt_legacy %>% dplyr::filter(YEAR==2023, alpage %in% alpages_noViso, habitat_code %in% c(1,9,5,6)) %>%
  dplyr::transmute(AUCg, dah=as.numeric(DAH), SMOD_2023=as.numeric(SMOD), alpage=factor(alpage),
                   Charge_log=as.numeric(Charge_log), charge_bin=factor(charge_bin),
                   hab4=factor(as.character(habitat_code), levels=c("1","9","5","6"), labels=hab_levels)) %>%
  dplyr::filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
mod_2023_hab_noViso <- bam(
  log(AUCg) ~ s(dah,k=10,bs="ts") + s(SMOD_2023,k=10,bs="ts") +
    ti(dah,SMOD_2023,k=c(20,20)) + hab4 + s(alpage,bs="re") + charge_bin +
    s(Charge_log,k=10,bs="ts",by=interaction(hab4,charge_bin)) +
    ti(Charge_log,SMOD_2023,k=c(15,15),by=interaction(hab4,charge_bin)),
  data=d2023_hab_noViso, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, na.action=na.exclude
)

# Registres pour itérer proprement
m_all3_ref <- list(`2017`=mod_2017_ref_all3, `2018`=mod_2018_ref_all3, `2019`=mod_2019_ref_all3,
                   `2020`=mod_2020_ref_all3, `2021`=mod_2021_ref_all3, `2022`=mod_2022_ref_all3, `2023`=mod_2023_ref_all3)
d_all3_ref <- list(`2017`=data_2017_all3, `2018`=data_2018_all3, `2019`=data_2019_all3,
                   `2020`=data_2020_all3, `2021`=data_2021_all3, `2022`=data_2022_all3, `2023`=data_2023_all3)

m_all3_hab <- list(`2017`=mod_2017_hab_all3, `2018`=mod_2018_hab_all3, `2019`=mod_2019_hab_all3,
                   `2020`=mod_2020_hab_all3, `2021`=mod_2021_hab_all3, `2022`=mod_2022_hab_all3, `2023`=mod_2023_hab_all3)
d_all3_hab <- list(`2017`=d2017_hab_all3, `2018`=d2018_hab_all3, `2019`=d2019_hab_all3,
                   `2020`=d2020_hab_all3, `2021`=d2021_hab_all3, `2022`=d2022_hab_all3, `2023`=d2023_hab_all3)

m_noV_ref   <- list(`2017`=mod_2017_ref_noViso, `2018`=mod_2018_ref_noViso, `2019`=mod_2019_ref_noViso,
                    `2020`=mod_2020_ref_noViso, `2021`=mod_2021_ref_noViso, `2022`=mod_2022_ref_noViso, `2023`=mod_2023_ref_noViso)
d_noV_ref   <- list(`2017`=data_2017_noViso, `2018`=data_2018_noViso, `2019`=data_2019_noViso,
                    `2020`=data_2020_noViso, `2021`=data_2021_noViso, `2022`=data_2022_noViso, `2023`=data_2023_noViso)

m_noV_hab   <- list(`2017`=mod_2017_hab_noViso, `2018`=mod_2018_hab_noViso, `2019`=mod_2019_hab_noViso,
                    `2020`=mod_2020_hab_noViso, `2021`=mod_2021_hab_noViso, `2022`=mod_2022_hab_noViso, `2023`=mod_2023_hab_noViso)
d_noV_hab   <- list(`2017`=d2017_hab_noViso, `2018`=d2018_hab_noViso, `2019`=d2019_hab_noViso,
                    `2020`=d2020_hab_noViso, `2021`=d2021_hab_noViso, `2022`=d2022_hab_noViso, `2023`=d2023_hab_noViso)

# =========================
# 2) Effet +10% (marginalisé) & GRAPHE 1 (3 alpages)
# =========================
minN <- 20
years <- names(m_all3_ref)

res_all3 <- tibble(YEAR=integer(), groupe=character(),
                   pct10=numeric(), lwr=numeric(), upr=numeric(), n=integer())

for (yy in years) {
  # TOTAL
  modG <- m_all3_ref[[yy]]; datG <- d_all3_ref[[yy]]
  dG <- datG %>% dplyr::filter(charge_bin==1)
  if (nrow(dG) >= minN) {
    ch0 <- pmax(exp(dG$Charge_log)-1, 0); ch1 <- ch0*1.10
    nd0 <- dG; nd1 <- dG; nd0$Charge_log <- log1p(ch0); nd1$Charge_log <- log1p(ch1)
    X1 <- predict(modG, newdata=nd1, type="lpmatrix"); X0 <- predict(modG, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1-X0); b <- coef(modG); V <- vcov(modG)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    res_all3 <- dplyr::bind_rows(res_all3, tibble(
      YEAR=as.integer(yy), groupe="Total",
      pct10=100*(exp(est)-1), lwr=100*(exp(est-1.96*se)-1), upr=100*(exp(est+1.96*se)-1),
      n=nrow(dG)
    ))
  }
  # 4 HABITATS
  modH <- m_all3_hab[[yy]]; datH <- d_all3_hab[[yy]]
  for (lab in hab_levels) {
    dH <- datH %>% dplyr::filter(hab4==lab, charge_bin==1)
    if (nrow(dH) < minN) next
    ch0 <- pmax(exp(dH$Charge_log)-1, 0); ch1 <- ch0*1.10
    nd0 <- dH; nd1 <- dH; nd0$Charge_log <- log1p(ch0); nd1$Charge_log <- log1p(ch1)
    X1 <- predict(modH, newdata=nd1, type="lpmatrix"); X0 <- predict(modH, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1-X0); b <- coef(modH); V <- vcov(modH)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    res_all3 <- dplyr::bind_rows(res_all3, tibble(
      YEAR=as.integer(yy), groupe=lab,
      pct10=100*(exp(est)-1), lwr=100*(exp(est-1.96*se)-1), upr=100*(exp(est+1.96*se)-1),
      n=nrow(dH)
    ))
  }
}

res_all3 <- res_all3 %>%
  dplyr::mutate(YEAR=factor(YEAR, levels=2017:2023),
                groupe=factor(groupe, levels=c("Total", hab_levels)),
                sig=(lwr>0)|(upr<0))

cols <- c("Total"="grey60",
          "P. nivales"="#3B5BDB","P. thermiques écorchées"="#E03131",
          "Queyrellins"="#F4A261","Nardaies denses du subalpin"="#1B9E77")

p_all3 <- ggplot(res_all3, aes(YEAR, pct10, fill=groupe)) +
  geom_col(position=position_dodge2(width=0.9, preserve="single"),
           width=0.86, aes(alpha=ifelse(sig,1,0.45))) +
  geom_errorbar(aes(ymin=lwr, ymax=upr),
                position=position_dodge2(width=0.9, preserve="single"), width=0.22) +
  scale_fill_manual(values=cols, name=NULL, drop=FALSE) +
  scale_alpha_identity() +
  geom_hline(yintercept=0, linewidth=0.5) +
  labs(title="Élasticité +10% de charge — 3 alpages (SMOD non fixé)",
       subtitle="Groupes = Total + 4 habitats. Barres = IC95%; transparence = non significatif.",
       x="Année", y="% variation de GPROD pour +10% de charge") +
  theme_minimal(base_size=14) + theme(panel.grid.minor=element_blank(), legend.position="right")

p_all3
# ggsave("bar_all3_alpages.png", p_all3, width=10, height=6, dpi=300)

# =========================
# 3) Effet +10% & GRAPHE 2 (sans Viso)
# =========================
res_noV <- tibble(YEAR=integer(), groupe=character(),
                  pct10=numeric(), lwr=numeric(), upr=numeric(), n=integer())

for (yy in names(m_noV_ref)) {
  # TOTAL
  modG <- m_noV_ref[[yy]]; datG <- d_noV_ref[[yy]]
  dG <- datG %>% dplyr::filter(charge_bin==1)
  if (nrow(dG) >= minN) {
    ch0 <- pmax(exp(dG$Charge_log)-1, 0); ch1 <- ch0*1.10
    nd0 <- dG; nd1 <- dG; nd0$Charge_log <- log1p(ch0); nd1$Charge_log <- log1p(ch1)
    X1 <- predict(modG, newdata=nd1, type="lpmatrix"); X0 <- predict(modG, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1-X0); b <- coef(modG); V <- vcov(modG)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    res_noV <- dplyr::bind_rows(res_noV, tibble(
      YEAR=as.integer(yy), groupe="Total",
      pct10=100*(exp(est)-1), lwr=100*(exp(est-1.96*se)-1), upr=100*(exp(est+1.96*se)-1),
      n=nrow(dG)
    ))
  }
  # 4 HABITATS
  modH <- m_noV_hab[[yy]]; datH <- d_noV_hab[[yy]]
  for (lab in hab_levels) {
    dH <- datH %>% dplyr::filter(hab4==lab, charge_bin==1)
    if (nrow(dH) < minN) next
    ch0 <- pmax(exp(dH$Charge_log)-1, 0); ch1 <- ch0*1.10
    nd0 <- dH; nd1 <- dH; nd0$Charge_log <- log1p(ch0); nd1$Charge_log <- log1p(ch1)
    X1 <- predict(modH, newdata=nd1, type="lpmatrix"); X0 <- predict(modH, newdata=nd0, type="lpmatrix")
    dbar <- colMeans(X1-X0); b <- coef(modH); V <- vcov(modH)
    est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
    res_noV <- dplyr::bind_rows(res_noV, tibble(
      YEAR=as.integer(yy), groupe=lab,
      pct10=100*(exp(est)-1), lwr=100*(exp(est-1.96*se)-1), upr=100*(exp(est+1.96*se)-1),
      n=nrow(dH)
    ))
  }
}

res_noV <- res_noV %>%
  dplyr::mutate(YEAR=factor(YEAR, levels=2017:2023),
                groupe=factor(groupe, levels=c("Total", hab_levels)),
                sig=(lwr>0)|(upr<0))

p_noViso <- ggplot(res_noV, aes(YEAR, pct10, fill=groupe)) +
  geom_col(position=position_dodge2(width=0.9, preserve="single"),
           width=0.86, aes(alpha=ifelse(sig,1,0.45))) +
  geom_errorbar(aes(ymin=lwr, ymax=upr),
                position=position_dodge2(width=0.9, preserve="single"), width=0.22) +
  scale_fill_manual(values=cols, name=NULL, drop=FALSE) +
  scale_alpha_identity() +
  geom_hline(yintercept=0, linewidth=0.5) +
  labs(title="Élasticité +10% de charge — sans Viso (SMOD non fixé)",
       subtitle="Groupes = Total + 4 habitats. Barres = IC95%; transparence = non significatif.",
       x="Année", y="% variation de GPROD pour +10% de charge") +
  theme_minimal(base_size=14) + theme(panel.grid.minor=element_blank(), legend.position="right")

p_noViso
# ggsave("bar_noViso.png", p_noViso, width=10, height=6, dpi=300)

    







# ==== Packages ====
library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(parallel)

# ==== Préparation ====
# AUCg de secours
if (!"AUCg" %in% names(dt_legacy) && "GPROD" %in% names(dt_legacy)) {
  dt_legacy$AUCg <- dt_legacy$GPROD
}
dt_legacy$AUCg <- pmax(dt_legacy$AUCg, .Machine$double.eps)

years  <- 2017:2023
n_th   <- max(1, parallel::detectCores() - 1)
x_thr  <- log(2)  # seuil charge=1

# Classes SMOD
class_bounds <- tibble::tibble(
  groupe = factor(c("Total",
                    "Déneigement : Précoce (≤120)",
                    "Déneigement : Moyen (121–150)",
                    "Déneigement : Tardif (≥151)"),
                  levels = c("Total",
                             "Déneigement : Précoce (≤120)",
                             "Déneigement : Moyen (121–150)",
                             "Déneigement : Tardif (≥151)")),
  z_min = c(-Inf, -Inf, 121, 151),
  z_max = c( Inf, 120 , 150,  Inf)
)

pal <- c("Total" = "grey50",
         "Déneigement : Précoce (≤120)" = "#fed976",
         "Déneigement : Moyen (121–150)"= "#41b6c4",
         "Déneigement : Tardif (≥151)"  = "#0c2c84")

# ==== Fonction utilitaire : effet +10% pour un sous-ensemble ====
effet_plus10pct <- function(mod, df_sub) {
  if (nrow(df_sub) < 50) return(c(NA_real_, NA_real_, NA_real_, n = nrow(df_sub)))  # trop peu
  
  # 1) données "avant" (état observé)
  nd1 <- df_sub %>%
    mutate(
      charge_bin = factor(as.integer(Charge_log >= x_thr),
                          levels = levels(model.frame(mod)$charge_bin))
    )
  
  # 2) données "après" (+10% de charge brute -> Charge_log2)
  C      <- pmax(exp(nd1$Charge_log) - 1, 0)
  C2     <- 1.1 * C
  nd2    <- nd1
  nd2$Charge_log <- log1p(C2)
  nd2$charge_bin <- factor(as.integer(nd2$Charge_log >= x_thr),
                           levels = levels(model.frame(mod)$charge_bin))
  
  # Matrices de design sur l’échelle lien
  X1 <- predict(mod, newdata = nd1, type = "lpmatrix")
  X2 <- predict(mod, newdata = nd2, type = "lpmatrix")
  Xd_bar <- colMeans(X2 - X1)                  # moyenne des deltas de lignes
  beta   <- coef(mod)
  Vp     <- vcov(mod)
  
  dEta_bar <- as.numeric(Xd_bar %*% beta)
  se_dEta  <- sqrt(as.numeric(Xd_bar %*% Vp %*% Xd_bar))
  
  # Passage en pourcentage
  eff  <- 100 * (exp(dEta_bar) - 1)
  se_p <- 100 * exp(dEta_bar) * se_dEta
  lo   <- eff - 1.96 * se_p
  hi   <- eff + 1.96 * se_p
  
  c(eff = eff, lo = lo, hi = hi, n = nrow(df_sub))
}

# ==== Boucle années : fit + extraction effets ====
res <- list()

for (yy in years) {
  d_yr <- dt_legacy %>%
    filter(YEAR == yy) %>%
    transmute(
      AUCg,
      dah        = as.numeric(DAH),
      SMOD_2023  = as.numeric(SMOD),
      alpage     = factor(alpage),
      Charge_log = as.numeric(Charge_log),
      charge_bin = factor(charge_bin)
    ) %>%
    filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
  
  # Fit (un modèle global par année, alpage en aléatoire)
  mod <- bam(
    log(AUCg) ~
      s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
      ti(dah, SMOD_2023, k=c(20,20)) +
      s(alpage, bs="re") +
      charge_bin +
      s(Charge_log, k=10, bs="ts", by=charge_bin) +
      ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin),
    data = d_yr, method = "fREML", select = TRUE,
    discrete = TRUE, nthreads = n_th, na.action = na.exclude
  )
  
  # Sous-ensembles : Total + classes
  out_year <- purrr::map_dfr(seq_len(nrow(class_bounds)), function(i){
    zmin <- class_bounds$z_min[i]; zmax <- class_bounds$z_max[i]
    grp  <- class_bounds$groupe[i]
    
    df_sub <- if (grp == "Total") d_yr else d_yr %>% dplyr::filter(SMOD_2023 >= zmin, SMOD_2023 <= zmax)
    
    eff <- effet_plus10pct(mod, df_sub)
    
    tibble::tibble(
      YEAR   = yy,
      groupe = grp,
      eff    = eff["eff"],
      lo     = eff["lo"],
      hi     = eff["hi"],
      n      = as.integer(eff["n"])
    )
  })
  
  res[[as.character(yy)]] <- out_year
}

eff_df <- dplyr::bind_rows(res) %>%
  mutate(
    groupe = factor(groupe, levels = levels(class_bounds$groupe)),
    sig    = !is.na(eff) & lo * hi > 0,        # IC ne recouvre pas 0
    alpha  = ifelse(is.na(eff) | n < 150, 0.35, 1)  # transparence si peu de pixels
  )

# ==== PLOT ====
gg <- ggplot(eff_df, aes(x = factor(YEAR), y = eff, fill = groupe, alpha = alpha)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "grey20") +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.8), width = 0.35, linewidth = 0.5) +
  scale_fill_manual(values = pal, name = "Groupe") +
  scale_alpha_identity() +
  labs(
    title = "Élasticité +10% de charge — SMOD non fixé (par classes)",
    subtitle = "Barres = IC95%. Transparence si n < 150. Modèle : 1 par année, alpage en effet aléatoire.",
    x = "Année",
    y = "% variation de GPROD pour +10% de charge"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(gg)

# Export (optionnel)
# ggsave("Barplot_effet_plus10pct_par_classe.png", gg, width = 12, height = 6, dpi = 300)













# =========================
# Barplot +10% par classes, aligné aux courbes (modèles PAR ALPAGE)
# =========================
library(dplyr); library(tidyr); library(ggplot2); library(mgcv); library(purrr); library(parallel)

# Secours AUCg
if (!"AUCg" %in% names(dt_legacy) && "GPROD" %in% names(dt_legacy)) dt_legacy$AUCg <- dt_legacy$GPROD
dt_legacy$AUCg <- pmax(dt_legacy$AUCg, .Machine$double.eps)

years <- 2017:2023
alpages_focus <- c("Viso","Cayolle","Sanguiniere")  # adapte si besoin
n_th <- max(1, parallel::detectCores()-1)
x_thr <- log(2)

# Classes SMOD
class_bounds <- tibble::tibble(
  groupe = factor(c("Total", "Déneigement : Précoce (≤120)", "Déneigement : Moyen (121–150)", "Déneigement : Tardif (≥151)"),
                  levels = c("Total","Déneigement : Précoce (≤120)","Déneigement : Moyen (121–150)","Déneigement : Tardif (≥151)")),
  z_min = c(-Inf, -Inf, 121, 151),
  z_max = c( Inf, 120 , 150,  Inf)
)

pal <- c("Total"="grey50",
         "Déneigement : Précoce (≤120)"="#fed976",
         "Déneigement : Moyen (121–150)"="#41b6c4",
         "Déneigement : Tardif (≥151)"  ="#0c2c84")

# Effet +10% (retourne dEta moyen et son SE, pas encore en %)
deta_plus10pct <- function(mod, df_sub){
  if (nrow(df_sub) < 50) return(c(deta=NA_real_, se=NA_real_, n=nrow(df_sub)))
  
  # On restreint aux pixels chargés (après seuil) et dans la plage observée de charge
  qs  <- quantile(df_sub$Charge_log, c(.01,.99), na.rm=TRUE)
  df1 <- df_sub %>% filter(charge_bin==1, between(Charge_log, qs[1], qs[2]))
  if (nrow(df1) < 30) return(c(deta=NA_real_, se=NA_real_, n=nrow(df1)))
  
  # Avant/après +10% de charge brute
  C   <- pmax(exp(df1$Charge_log)-1, 0)
  C2  <- 1.1 * C
  nd1 <- df1
  nd2 <- df1; nd2$Charge_log <- log1p(C2)
  nd2$charge_bin <- factor(as.integer(nd2$Charge_log >= x_thr),
                           levels = levels(model.frame(mod)$charge_bin))
  
  X1 <- predict(mod, newdata=nd1, type="lpmatrix")
  X2 <- predict(mod, newdata=nd2, type="lpmatrix")
  Xd_bar <- colMeans(X2 - X1)
  
  beta <- coef(mod); Vp <- vcov(mod)
  deta <- as.numeric(Xd_bar %*% beta)
  se   <- sqrt(as.numeric(Xd_bar %*% Vp %*% Xd_bar))
  c(deta=deta, se=se, n=nrow(df1))
}

# Fit par ANNÉE × ALPAGE, calcule effets par CLASSE, puis AGRÈGE entre alpages
res <- list()

for (yy in years) {
  d_yr <- dt_legacy %>%
    filter(YEAR==yy, alpage %in% alpages_focus) %>%
    transmute(
      AUCg, dah = as.numeric(DAH),
      SMOD_2023 = as.numeric(SMOD),
      alpage = factor(alpage, levels=alpages_focus),
      Charge_log = as.numeric(Charge_log),
      charge_bin = factor(charge_bin)
    ) %>%
    filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
  
  # Modèles par alpage (même structure que tes courbes)
  mods <- lapply(split(d_yr, d_yr$alpage), function(dd){
    if (nrow(dd) < 200) return(NULL)
    bam(
      log(AUCg) ~
        s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
        ti(dah, SMOD_2023, k=c(20,20)) +
        charge_bin +
        s(Charge_log, k=10, bs="ts", by=charge_bin) +
        ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin),
      data=dd, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=n_th, na.action=na.exclude
    )
  })
  
  # Effets par alpage × classe (on pèse par n)
  by_alp <- purrr::imap_dfr(mods, function(m, alp){
    if (is.null(m)) return(tibble())
    dd <- filter(d_yr, alpage==alp)
    
    purrr::map_dfr(seq_len(nrow(class_bounds)), function(i){
      g   <- class_bounds$groupe[i]
      zlo <- class_bounds$z_min[i]; zhi <- class_bounds$z_max[i]
      sub <- if (g=="Total") dd else dd %>% filter(SMOD_2023>=zlo, SMOD_2023<=zhi)
      
      out <- deta_plus10pct(m, sub)
      tibble(alpage=alp, YEAR=yy, groupe=g,
             deta=out["deta"], se_deta=out["se"], n=as.integer(out["n"]))
    })
  })
  
  # Agrégation entre alpages (poids = n)
  agg <- by_alp %>%
    group_by(YEAR, groupe) %>%
    summarise(
      wsum = sum(n, na.rm=TRUE),
      deta = sum(n * deta, na.rm=TRUE) / pmax(wsum, 1),
      var_ = sum((n^2) * (se_deta^2), na.rm=TRUE) / pmax(wsum, 1)^2,
      n    = wsum,
      .groups="drop"
    ) %>%
    mutate(
      eff = 100 * (exp(deta) - 1),
      se_eff = 100 * exp(deta) * sqrt(var_),
      lo  = eff - 1.96*se_eff,
      hi  = eff + 1.96*se_eff
    )
  
  res[[as.character(yy)]] <- agg
}

eff_df <- bind_rows(res) %>%
  mutate(
    groupe = factor(groupe, levels=levels(class_bounds$groupe)),
    alpha  = ifelse(is.na(eff) | n < 150, 0.35, 1)
  )

# ---- PLOT
gg <- ggplot(eff_df, aes(x=factor(YEAR), y=eff, fill=groupe, alpha=alpha)) +
  geom_hline(yintercept=0, colour="grey50") +
  geom_col(position=position_dodge(0.8), width=0.72, colour="grey25") +
  geom_errorbar(aes(ymin=lo, ymax=hi), position=position_dodge(0.8), width=0.34, linewidth=0.5) +
  scale_fill_manual(values=pal, name="Groupe") +
  scale_alpha_identity() +
  labs(
    title="Élasticité +10% de charge — agrégé d’après modèles par alpage (SMOD non fixé)",
    subtitle="Barres = IC95%. Moyenne après-seuil (charge_bin=1), pondérée par n. Transparence si n<150.",
    x="Année", y="% variation de GPROD pour +10% de charge"
  ) +
  theme_minimal(base_size=13) +
  theme(legend.position="right", panel.grid.minor=element_blank())
print(gg)







# =========================
# Barplot +10% par classes — versions "alignées courbes"
# =========================
library(dplyr); library(tidyr); library(ggplot2); library(mgcv); library(purrr); library(parallel)

# Secours AUCg
if (!"AUCg" %in% names(dt_legacy) && "GPROD" %in% names(dt_legacy)) dt_legacy$AUCg <- dt_legacy$GPROD
dt_legacy$AUCg <- pmax(dt_legacy$AUCg, .Machine$double.eps)

years <- 2017:2023
alpages_focus <- c("Viso","Cayolle","Sanguiniere")   # adapte si besoin
n_th <- max(1, parallel::detectCores()-1)
x_thr <- log(2)

# Classes SMOD
class_bounds <- tibble::tibble(
  groupe = factor(c("Total",
                    "Déneigement : Précoce (≤120)",
                    "Déneigement : Moyen (121–150)",
                    "Déneigement : Tardif (≥151)"),
                  levels = c("Total",
                             "Déneigement : Précoce (≤120)",
                             "Déneigement : Moyen (121–150)",
                             "Déneigement : Tardif (≥151)")),
  z_min = c(-Inf, -Inf, 121, 151),
  z_max = c( Inf, 120 , 150,  Inf)
)

pal <- c("Total"="grey50",
         "Déneigement : Précoce (≤120)"="#fed976",
         "Déneigement : Moyen (121–150)"="#41b6c4",
         "Déneigement : Tardif (≥151)"  ="#0c2c84")

# -------- Effet +10% au NIVEAU PIXEL, sur échelle réponse, moyenne + delta-method
eff_plus10_pixel_mean <- function(mod, df_sub, trim=c(.01,.99)){
  if (nrow(df_sub) < 50) return(c(eff=NA_real_, lo=NA_real_, hi=NA_real_, n=nrow(df_sub)))
  
  # pixels après seuil et dans la plage observée (trim 1-99%)
  qs  <- quantile(df_sub$Charge_log, trim, na.rm=TRUE)
  D   <- df_sub %>% filter(charge_bin==1, between(Charge_log, qs[1], qs[2]))
  if (nrow(D) < 30) return(c(eff=NA_real_, lo=NA_real_, hi=NA_real_, n=nrow(D)))
  
  # Avant / Après +10% (charge brute)
  C    <- pmax(exp(D$Charge_log)-1, 0)
  C2   <- 1.1 * C
  nd1  <- D
  nd2  <- D; nd2$Charge_log <- log1p(C2)
  nd2$charge_bin <- factor(as.integer(nd2$Charge_log >= x_thr),
                           levels = levels(model.frame(mod)$charge_bin))
  
  X1 <- predict(mod, nd1, type="lpmatrix")
  X2 <- predict(mod, nd2, type="lpmatrix")
  Xd <- X2 - X1
  
  beta <- coef(mod); Vp <- vcov(mod)
  dEta <- as.numeric(Xd %*% beta)     # Δη_i
  r_i  <- exp(dEta) - 1               # variation relative (réponse) par pixel
  
  # moyenne en %
  eff <- 100 * mean(r_i)
  
  # delta-method pour la moyenne de exp(Δη) :
  # ∂/∂β mean(exp(Δη_i)) = mean(exp(Δη_i) * Xd_i)
  g   <- colMeans(exp(dEta) * Xd) * 100
  se  <- sqrt(as.numeric(t(g) %*% Vp %*% g))
  
  lo <- eff - 1.96*se
  hi <- eff + 1.96*se
  
  c(eff=eff, lo=lo, hi=hi, n=nrow(D))
}

# -------- Fit par ANNÉE × ALPAGE, calcul par CLASSE, puis agrégation entre alpages
make_barplot_df <- function(weight_scheme=c("n","equal")) {
  weight_scheme <- match.arg(weight_scheme)
  res <- list()
  
  for (yy in years) {
    d_yr <- dt_legacy %>%
      filter(YEAR==yy, alpage %in% alpages_focus) %>%
      transmute(
        AUCg, dah = as.numeric(DAH),
        SMOD_2023 = as.numeric(SMOD),
        alpage = factor(alpage, levels=alpages_focus),
        Charge_log = as.numeric(Charge_log),
        charge_bin = factor(charge_bin)
      ) %>%
      filter(is.finite(AUCg), is.finite(dah), is.finite(SMOD_2023), is.finite(Charge_log))
    
    # Modèles par alpage (comme les courbes)
    mods <- lapply(split(d_yr, d_yr$alpage), function(dd){
      if (nrow(dd) < 200) return(NULL)
      bam(
        log(AUCg) ~
          s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
          ti(dah, SMOD_2023, k=c(20,20)) +
          charge_bin +
          s(Charge_log, k=10, bs="ts", by=charge_bin) +
          ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin),
        data=dd, method="fREML", select=TRUE,
        discrete=TRUE, nthreads=n_th, na.action=na.exclude
      )
    })
    
    # Effets alpage × classe
    by_alp <- purrr::imap_dfr(mods, function(m, alp){
      if (is.null(m)) return(tibble())
      dd <- filter(d_yr, alpage==alp)
      
      purrr::map_dfr(seq_len(nrow(class_bounds)), function(i){
        g   <- class_bounds$groupe[i]
        zlo <- class_bounds$z_min[i]; zhi <- class_bounds$z_max[i]
        sub <- if (g=="Total") dd else dd %>% filter(SMOD_2023>=zlo, SMOD_2023<=zhi)
        
        out <- eff_plus10_pixel_mean(m, sub)
        tibble(alpage=alp, YEAR=yy, groupe=g,
               eff=out["eff"], lo=out["lo"], hi=out["hi"], n=as.integer(out["n"]))
      })
    })
    
    # Agrégation entre alpages
    agg <- by_alp %>%
      group_by(YEAR, groupe) %>%
      { if (weight_scheme=="n") {
        summarise(., eff = weighted.mean(eff, w = n, na.rm=TRUE),
                  # delta pour la moyenne pondérée : approx via variances individuelles
                  se  = sqrt(weighted.mean(((hi-eff)/1.96)^2, w=n, na.rm=TRUE)),
                  n   = sum(n, na.rm=TRUE), .groups="drop")
      } else {
        summarise(., eff = mean(eff, na.rm=TRUE),
                  se  = sqrt(mean(((hi-eff)/1.96)^2, na.rm=TRUE)),
                  n   = sum(n, na.rm=TRUE), .groups="drop")
      }} %>%
      mutate(lo = eff - 1.96*se, hi = eff + 1.96*se)
    
    res[[as.character(yy)]] <- agg
  }
  
  bind_rows(res) %>%
    mutate(groupe = factor(groupe, levels=levels(class_bounds$groupe)),
           alpha  = ifelse(is.na(eff) | n<150, 0.35, 1))
}

# --- Version globale (poids = n, i.e. « tout alpages confondus »)
eff_nweight  <- make_barplot_df(weight_scheme = "n")
# --- Version « alpages égaux » (pour éviter que Viso domine)
eff_eqweight <- make_barplot_df(weight_scheme = "equal")

plot_bar <- function(df, titre){
  ggplot(df, aes(x=factor(YEAR), y=eff, fill=groupe, alpha=alpha)) +
    geom_hline(yintercept=0, colour="grey55") +
    geom_col(position=position_dodge(0.8), width=0.72, colour="grey25") +
    geom_errorbar(aes(ymin=lo, ymax=hi), position=position_dodge(0.8), width=0.34, linewidth=0.5) +
    scale_fill_manual(values=pal, name="Groupe") +
    scale_alpha_identity() +
    labs(title=titre,
         subtitle="Moyenne sur échelle réponse, pixels chargés (charge_bin=1), plage 1–99% de charge.\nIC95% via delta-method; transparence si n<150.",
         x="Année", y="% variation de GPROD pour +10% de charge") +
    theme_minimal(base_size=13) +
    theme(legend.position="right", panel.grid.minor=element_blank())
}

p1 <- plot_bar(eff_nweight,  "Élasticité +10% — agrégation pondérée par n (tous alpages)")
p2 <- plot_bar(eff_eqweight, "Élasticité +10% — agrégation alpages égaux")

print(p1); print(p2)











}
