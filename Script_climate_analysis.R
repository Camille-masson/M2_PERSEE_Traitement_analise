####------------- SCRIPT : Variations climatiques inter-annuel -------------####
## Script rédiger dans le cadre du rapport de stage de M2 visant a annalyser 
## l'influence des variations climatiques inter-annuel sur le schéma de paturage
## de l'alapge. 

## Observer d'eventuel report de charge entre des années contrastées (printemps
## précoce ou tardif)


#### 0. LIBRARIES AND CONSTANTS ####
#----------------------------------#

gc()
# Chargement de la configuration
source("config.R")
source(file.path(functions_dir, "Functions_filtering.R")) 


# Définition de l'année d'analyse
YEAR <- 2023
TYPE <- "catlog" #Type de données d'entrée (CATLOG, OFB )
alpage <- "Cayolle"
alpages <- c("Cayolle","Viso","Sanguiniere")
# Liste complète des alpages 2023 : "Cayolle", "Crouzet", "Grande-Cabane", "Lanchatra", "Rouanette", "Sanguiniere", "Vacherie-de-Roubion", "Viso"
# Liste complète des alpages 2022 : "Cayolle", "Combe-Madame", "Grande-Fesse", "Jas-des-Lievres", "Lanchatra", "Pelvas", "Sanguiniere", "Viso"

ALPAGES_TOTAL <- list(
  "9999" = c("Alpage_demo"),
  "2013" = c("Combe-Madame"),
  "2014" = c("Combe-Madame"),
  "2015" = c("Combe-Madame"),
  "2016" = c("Combe-Madame"),
  "2017" = c("Combe-Madame"),
  "2018" = c("Ane-et-Buyant", "Bedina", "Pesee", "Sept-Laux"),
  "2019" = c("Ane-et-Buyant", "Bedina", "Pesee", "Sept-Laux"),
  "2020" = c("Ane-et-Buyant", "Bedina", "Pesee","Rieuxclaret", "Sept-Laux"),
  "2021" = c("Ane-et-Buyant", "Bedina", "Pesee","Combe-Madame", "Sept-Laux"),
  "2022" = c("Ane-et-Buyant", "Bedina", "Cayolle", "Combe-Madame", "Grande-Fesse", "Jas-des-Lievres", "Lanchatra", "Pelvas","Pesee", "Sanguiniere","Sept-Laux", "Viso"),
  "2023" = c("Ane-et-Buyant", "Bedina", "Cayolle", "Crouzet", "Combe", "Combe-Madame", "Grande-Cabane", "Lanchatra", "Pesee", "Rouanette", "Sanguiniere", "Sept-Laux", "Vacherie-de-Roubion", "Viso"),
  "2024" = c("Viso", "Cayolle", "Sanguiniere")
)
ALPAGES <- ALPAGES_TOTAL[[as.character(YEAR)]]



#### 1. Préparation des données de SMOD ####
#------------------------------------------#
if (TRUE) {
  # Crop du smod par alpage 
  # Calcul du SMOD moyen sur 10 ans (2023 - 2013)
  # Calcul du fsca (fraction of snow cover)
  
  year_1 = "1999-09-01"
  year_2 = "2023-08-31"
  
  # LIBRARY & FUNCTION
  library(raster)
  library(dplyr)  
  library(stringr) 
  library(sf)
  library(tools)
  library(ggplot2)
  library(terra)    # pour lire les rasters SMOD
  library(dplyr)
  library(lubridate)
  source(file.path(functions_dir, "Functions_traitement_smod.R"))

  
  # ENTREE
  
  # Un dossier contenant les SMOD
  case_SMOD_file = file.path(raster_dir, "Snow")
  
  # Un .CSV : correspondance entre les dalles du SMOD et les alpages
  table_corresp_file = file.path(case_SMOD_file, "SMOD_ID_by_alpage.csv")
  
  # Un dossier contenant les ratsers des Unités Pastorales (UP)
  case_UP_file = file.path(raster_dir, "UP")
  
  # Un .SHP avec les Unités pastorales UP
  UP_file = file.path(case_UP_file, "v1_bd_shape_up_inra_2012_2014_2154_all_emprise.shp")
  
  # Un dossier contenant les Infos sur les alpages
  raw_data_dir = file.path(data_dir,paste0("Colliers_",YEAR,"_brutes"))
  # Un data.frame contenant les dates de pose et de retrait des colliers
  alpage_info_file <- file.path(raw_data_dir, paste0(YEAR,"_infos_alpages.csv"))
  
  # SORTIE
  
  #Création du dossier de sortie des indicateur pour la visualistaion
  output_clim_case <- file.path(output_dir, "8. Analysis_Climate")
  if (!dir.exists(output_clim_case)) {
    dir.create(output_clim_case, recursive = TRUE)
  }
  #Création du sous-dossier Indicateur traitée : Chargement
  output_SMOD_case <- file.path(output_clim_case, "SMOD")
  if (!dir.exists(output_SMOD_case)) {
    dir.create(output_SMOD_case, recursive = TRUE)
  }
  
  
  
  
  # CODE
  
  # Pipeline du FSCA : 
  # 1. Préparation des SMOD
  process_smod_for_alpage(alpage, year_1, year_2,case_SMOD_file,table_corresp_file,alpage_info_file,UP_file, output_SMOD_case)
  
  # 2. Calcul du SMOD médian
  median_file <- median_smod_for_alpage(alpage, output_SMOD_case)
  
  # 3. Calcul du FSCA
  fsca_file <- compute_fsca_for_alpage(alpage, output_SMOD_case)
  
  
  # Pipeline du Enneigement median (par : Alpage, Année)
  # Utilise la partie 1 du pipeline precédent : 1. Préparation des SMOD
  
  dataset_smod_median_alpage_year(output_clim_case, output_SMOD_case)
  
  
  
}
  
#### 2. Création du jeu de données : Alti, FSCA, Présence ####
#------------------------------------------------------------#
if (TRUE) {
  # Pourchaque pixel a 20 mètres attribution  : 
  # - FSCA
  # - Altitude
  # - Chargement / Densité de kernel
  # Donc soit l'utilisation par le SHP des polygones de kernel
  # Soit le taux de chargement par quizaine (méthode préféré pour les plots)
  
 
  # LIBRARY & FUNCTION
  library(raster)
  library(dplyr)  
  library(stringr) 
  library(sf)
  library(tools)
  library(ggplot2)
  source(file.path(functions_dir, "Functions_traitement_smod.R"))
  
  # ENTREE
  # Un .TIF du SMOD Médian sur 10 ans 
  clim_case <- file.path(output_dir, "8. Analysis_Climate")
  SMOD_case <- file.path(clim_case, "SMOD")
  SMOD_tif_file <- file.path(SMOD_case, paste0("Fsca_",alpage,".tif"))
  
  
  # Un .TIF du DEM de l'alpage
  Alti_case <- file.path(raster_dir, "Alti")
  MNT_tif_file <- file.path(Alti_case, paste0(alpage,"_MNT.tif"))
  
  
  # Un .SHP de l'utilisation par quinzaine
  visu_case <- file.path(output_dir, "5. Indicateurs_visualisation")
  polygon_case <- file.path(visu_case, "Utilisation_par_quinzaine")
  polygon_use_shp <- file.path(polygon_case, paste0("Use_polygon_", YEAR, "_", alpage, ".shp"))
  
  
  
  
  
  # Des .TIF du chargement par quizaine
  #Création du dossier de sortie des indicateur pour la visualistaion
  visu_case <- file.path(output_dir, "5. Indicateurs_visualisation")
  chargement_case <- file.path(visu_case, "Taux_chargement")
  
  
  case_alpage <- file.path(chargement_case, paste0(YEAR,"_",alpage))
  load_16_30_jun_tif_file <- file.path(case_alpage, paste0("by_quinzaine_",YEAR,"_",alpage,"_16_30_jun.tif"))
  load_1_15_jul_tif_file <- file.path(case_alpage, paste0("by_quinzaine_",YEAR,"_",alpage,"_1_15_jul.tif"))
  load_16_30_jul_tif_file <- file.path(case_alpage, paste0("by_quinzaine_",YEAR,"_",alpage,"_16_30_jul.tif"))
  load_1_15_aou_tif_file <- file.path(case_alpage, paste0("by_quinzaine_",YEAR,"_",alpage,"_1_15_aou.tif"))
  load_16_30_aou_tif_file <- file.path(case_alpage, paste0("by_quinzaine_",YEAR,"_",alpage,"_16_30_aou.tif"))
  load_1_15_sep_tif_file <- file.path(case_alpage, paste0("by_quinzaine_",YEAR,"_",alpage,"_1_15_sep.tif"))
  load_apres16_sep_tif_file <- file.path(case_alpage, paste0("by_quinzaine_",YEAR,"_",alpage,"_apres16_sep.tif"))
  
  
  # SORTIE
  #Création du dossier de sortie des indicateur pour la visualistaion
  output_clim_case <- file.path(output_dir, "8. Analysis_Climate")
  if (!dir.exists(output_clim_case)) {
    dir.create(output_clim_case, recursive = TRUE)
  }
  #Création du sous-dossier Indicateur traitée : Chargement
  output_data_case <- file.path(output_clim_case, "Data_Use_Fsca_Alti")
  if (!dir.exists(output_data_case)) {
    dir.create(output_data_case, recursive = TRUE)
  }
  
  # Un .RDS des données final
  output_clim_data_rds_file_shp <- file.path(output_data_case, paste0("Use_Fsca_Alti_",alpage,"_by_shp.rds"))
  # Un .RDS des données final
  output_clim_data_rds_file_raster <- file.path(output_data_case, paste0("Use_Fsca_Alti_",alpage,"_by_raster.rds"))
  # Un .RDS des données final
  output_clim_data_rds_file_raster_parc <- file.path(output_data_case, paste0("Use_Fsca_Alti_",alpage,"_by_raster_and_parc.rds"))
  output_rds_file_parc <- file.path(output_data_case, paste0("Use_Fsca_Alti_", alpage, "_by_raster_and_parc.rds"))
  
  # CODE 
  
  # Création du data.frame (sortie .RDS), basés sur les polygone de densité
  df_final <- generate_presence_data_by_quinzaine(SMOD_tif_file, MNT_tif_file, polygon_case,
                                                  years = 2022:2024, alpage , 
                                                  output_clim_data_rds_file = output_clim_data_rds_file_shp)
  
  
  # Création du data.frame (sortie .RDS), basés sur le chargement (raster)
  df_final <- generate_loading_data_by_quinzaine(SMOD_tif_file, MNT_tif_file, chargement_case, 
                                                 years = 2022:2024, alpage, 
                                                 output_rds_file = output_clim_data_rds_file_raster,
                                                 threshold = 10,res_raster = 10)
  
  
  
  df_final <- generate_loading_data_by_parc(SMOD_tif_file, MNT_tif_file, chargement_case,
                                            years= c(2022,2023,2024),alpage ,
                                            output_rds_file = output_clim_data_rds_file_raster_parc,
                                            threshold = 10, res_raster = 10
  )
  
  
  
  
}

#### 3. Création du plot de l'utilisation en fonction du climat ####
#------------------------------------------------------------------#
if (FALSE) {
  # Réglage :
  # - Réglages des points, seuil : 10 à 1000 ; Gradient de couleur log (cap à 500) par période
  # - Gradient basés sur l'intensité du chargement
  # - Polygon d'utilisation Calculé sur seuil chargement 100 à 1000 (polygone a 70%)
  # 
  # Les différents plots :
  # - Ellipse
  # - Point 
  # - Violin en fonction fsca
  # - Violin en focntion alti
  
  # LIBRARY & FUNCTION
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(ellipse)
  library(ggnewscale)
  library(scales)
  source(file.path(functions_dir, "Functions_plot_snow_ndvi.R"))
  
  #Paramètres d'alpages
  alpages <- "Sanguiniere"
  
  # ENTREE
  # Dossier général pour l'analyse climatique
  clim_case <- file.path(output_dir, "8. Analysis_Climate")
  # Création des dossiers si nécessaire
  if (!dir.exists(clim_case)) {
    dir.create(clim_case, recursive = TRUE)
  }
  # Dossier d'entrée contenant les fichiers RDS
  data_case <- file.path(clim_case, "Data_Use_Fsca_Alti")
  if (!dir.exists(data_case)) {
    dir.create(data_case, recursive = TRUE)
  }
  
  
  # SORTIE
  # Dossier de sortie pour les graphiques
  output_plot_case <- file.path(clim_case, "Graphique")
  if (!dir.exists(output_plot_case)) {
    dir.create(output_plot_case, recursive = TRUE)
  }
  
  
  
  # CODE
  
  # Plot avec les ellipses
  if(FALSE){
  plot_fsca_alti_elypse(
    alpages       = alpages,
    data_dir      = data_case,
    output_dir    = output_plot_case,
    years_to_use  = c(2022, 2023,2024),
    ellipse_level = 0.6)
    }
  
  # Plot avec les points
  if(TRUE){
  plot_fsca_alti_points(
    alpages, 
    data_dir = data_case, 
    output_dir = output_plot_case, 
    years_to_use = c(2022, 2023, 2024))
  }
  
  
  
  
  source(file.path(functions_dir, "Functions_plot_snow_ndvi.R"))
  
  
  
  
  plot_fsca_alti_points_parc_verif(
    alpages       = alpages,
    data_dir      = data_case,
    output_dir    = output_plot_case,
    years_to_use  = c(2022,2023, 2024)
  )
  
  
 
  plot_fsca_alti_points_parc(
    alpages       = alpages,
    data_dir      = data_case,
    output_dir    = output_plot_case,
    years_to_use  = c( 2022, 2023, 2024)
  )
  
  
  
  
  
  # Plot violin avec le fsca
  if(TRUE){
   plot_violin_fsca(
     alpages, 
     data_dir = data_case, 
     output_dir = output_plot_case, 
     years_to_use = c(2022, 2023, 2024))
   }
  
  
  
  
  # Plot violin avec l'altitude
  if(TRUE){
   plot_violin_alti(
     alpages,
     data_dir   = data_case,
     output_dir = output_plot_case,
     years_to_use = c(2022, 2023, 2024))
   }
  
  
  
  
  
  
}
  
#### 4. Analyse de la Phénologie ####
#-----------------------------------#
if (TRUE){
  ### 4.1 IRG par type habitat ###
  #------------------------------#
  
  # ---------------------------------------------
  # 4.1 IRG par type habitat
  # ---------------------------------------------
  # Description :
  # Préparation d'un jeu de données avec un cycle IRG par habitat et par an.
  # Chaque pixel aura une valeur IRG par jour (DOY) et son type d'habitat.
  # La table finale est sauvegardée au format .rds pour les analyses ultérieures
  # (courbes IRG par habitat).
  #
  # Entrées nécessaires :
  #  - Un stack multi-bandes IRG (.tif) : chaque bande = un jour de l'année (DOY).
  #    Exemple de nom de bande : "IRG_DOY121".
  #  - Raster de type d'habitat (pixel 10 m, valeurs entières catégorielles).
  #  - (Optionnel) Shapefile des Unités Pastorales (UP) pour réaliser un crop/mask.
  #
  # Sortie :
  #  - Un fichier .rds (output_IRG_by_habitat) contenant un data.frame long
  #    (colonnes : cell, x, y, habitat, DOY, IRG).
  
  
  # LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_traitement_smod.R"))
  
  # ENTREE
  #Dossier de phénologie 
  case_phenologie = file.path(raster_dir, "Phenologie")
  if (!dir.exists(case_phenologie)) {
    dir.create(case_phenologie, recursive = TRUE)
  }
  
  # Un .TIF de l'IRG catégorisé par jour et pixel de 10 mètre
  pheno_tif_file = file.path(case_phenologie, paste0("IRG_season",alpage,"_",YEAR,".tif"))
  
  
  # Un dossier contenant carte de végétation
  carto_file = file.path(raster_dir, "Classifications_fusion_ColorIndexed_sc1_landforms_mnh.tif")
  
  # Un dossier contenant les ratsers des Unités Pastorales (UP)
  case_UP_file = file.path(raster_dir, "UP")
  # Un .SHP avec les Unités pastorales UP
  UP_file = file.path(case_UP_file, "v1_bd_shape_up_inra_2012_2014_2154_all_emprise.shp")
  
  
  # SORTIE
  # SORTIE
  # Dossiers de sortie
  out_dir <- file.path(output_dir, "9. Analysis_Phenology", "data_IRG")
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  # Un .csv avec pour chaque pixel un valeur IRG par habitat
  output_IRG_by_habitat <- file.path(out_dir, paste0("IRG_by_habitat_",YEAR,"_",alpage,".csv"))
  
  
    
    
  # CODE  
  
  # Import de la carte de veget
  # Création du polygone de découpe à partir de UP
  # Récupérer l'UP 
  UP <- get_UP_shp(alpage, alpage_info_file, UP_file)
  # Conversion de UP (sf) en SpatVector (terra) et buffer de 300 m
  UP_vect <- terra::vect(UP)
  cropping_polygon <- terra::buffer(UP_vect, 300)
  # Typo de veget
  vegetation_typology_name <- get_alpage_info(alpage, alpage_info_file, "typologie_vegetation")
  
  # On aligne la rasterisation sur la grille contenue dans df_filtered (colonnes x et y)
  vegetation_df <- get_vegetation_rasterized(
    vegetation_file = carto_file,
    vegetation_typology_name = vegetation_typology_name,
    grid = df_filtered[c("x", "y")],
    cropping_polygon = cropping_polygon
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ---------------------------------------------
  # 4.1 IRG par type habitat (v2 – alignement robuste)
  # ---------------------------------------------
  # Description :
  # Préparation d'un jeu de données avec un cycle IRG par habitat et par an.
  # Chaque pixel aura une valeur IRG par jour (DOY) et son type d'habitat.
  # ---------------------------------------------
  # Entrées :
  #  - Stack IRG multi‑bandes (.tif) : noms « IRG_DOY121 », etc.
  #  - Raster habitat (catégoriel)
  #  - Shapefile UP (pour crop/mask)
  # ---------------------------------------------
  # Sortie : .rds long (cell, x, y, habitat, DOY, IRG)
  # ---------------------------------------------
  # Paramètres utilisateur ------------------------
  YEAR    <- 2021          # année analysée
  alpage  <- "Cayolle"       # alpage
  TYPE    <- "catlog"     # source données
  
  # Définissez functions_dir, raster_dir, output_dir avant exécution.
  # ---------------------------------------------------------------------------
  # Librairies ----------------------------------------------------------------
  library(terra)      # >= 1.7.0
  library(sf)
  library(stringr)
  library(tidyverse)
  
  # ---------------------------------------------------------------------------
  # Chemins d'accès -----------------------------------------------------------
  case_phenologie <- file.path(raster_dir, "Phenologie")
  pheno_tif_file  <- file.path(case_phenologie,
                               sprintf("IRG_season_%s_%d.tif", alpage, YEAR))
  
  carto_file      <- file.path(raster_dir,
                               "Classifications_fusion_ColorIndexed_sc1_landforms_mnh.tif")
  
  # LUT code ↔ label pour les habitats (CSV séparateur ;) ---------------------
  class_habitat_file <- file.path(raster_dir, "class_habitat.csv")
  
  case_UP_file    <- file.path(raster_dir, "UP")
  UP_file         <- file.path(case_UP_file,
                               "v1_bd_shape_up_inra_2012_2014_2154_all_emprise.shp")
  
  out_dir <- file.path(output_dir, "9. Analysis_Phenology", "data_IRG")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  output_IRG_by_habitat <- file.path(out_dir,
                                     sprintf("IRG_by_habitat_%d_%s.rds", YEAR, alpage))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  output_IRG_by_habitat <- file.path(out_dir,
                                     sprintf("IRG_by_habitat_%d_%s.rds", YEAR, alpage))
  
  # ---------------------------------------------------------------------------
  # 1. Lecture des rasters ----------------------------------------------------
  message("Lecture du stack IRG : ", pheno_tif_file)
  irg_stack <- terra::rast(pheno_tif_file)
  
  message("Lecture de la carte d'habitat : ", carto_file)
  habitat_r <- terra::rast(carto_file)
  
  # Harmonisation CRS ---------------------------------------------------------
  if (!is.na(terra::crs(habitat_r)) &&
      !terra::crs(irg_stack) == terra::crs(habitat_r)) {
    habitat_r <- terra::project(habitat_r, terra::crs(irg_stack), method = "near")
  }
  
  # ---------------------------------------------------------------------------
  # 2. Crop / mask par Unité Pastorale (optionnel) ----------------------------
  if (file.exists(UP_file)) {
    message("Application du masque UP (", alpage, ")")
    up_vect <- terra::vect(UP_file)
    
    # nom de champ contenant l'alpage
    alpage_field <- c("ALPAGE", "NOM_ALPA", "ALPAGE_NM")[
      c("ALPAGE", "NOM_ALPA", "ALPAGE_NM") %in% names(up_vect)][1]
    
    up_sel <- if (!is.na(alpage_field)) {
      up_vect[ up_vect[[alpage_field]] == alpage, ]
    } else up_vect
    
    up_sel <- terra::project(up_sel, terra::crs(irg_stack))
    
    # IRG : on conserve le grid d'origine comme référence
    irg_stack <- irg_stack |> terra::crop(up_sel) |> terra::mask(up_sel)
    # Habitat : découpe brute (grille propre), on resamplera ensuite
    habitat_r <- habitat_r |> terra::crop(up_sel) |> terra::mask(up_sel)
  }
  
  # ---------------------------------------------------------------------------
  # 3. Alignement habitat -> grille IRG ---------------------------------------
  #   On impose STRICTEMENT la grille (extent, nrow, ncol, res) de irg_stack.
  message("Alignement de la carte d'habitat sur la grille IRG…")
  habitat_r <- terra::resample(habitat_r, irg_stack, method = "near")
  
  # Vérif finale --------------------------------------------------------------
  if (!terra::compareGeom(irg_stack, habitat_r, stopOnError = FALSE)) {
    stop("Les géométries de IRG et habitat ne correspondent toujours pas – ",
         "vérifiez resolution/extent dans QGIS !")
  }
  
  # ---------------------------------------------------------------------------
  # 4. Renommage des bandes IRG avec le DOY -----------------------------------
  orig_names <- names(irg_stack)
  doy <- stringr::str_extract(orig_names, "\\d{3,}")
  if (any(is.na(doy)))
    stop("Les noms de bandes IRG doivent contenir le DOY (ex. 'IRG_DOY121').")
  
  names(irg_stack) <- sprintf("DOY%03d", as.integer(doy))
  
  # ---------------------------------------------------------------------------
  # 5. Fusion raster & passage au format long ---------------------------------
  
  # Lecture LUT habitat -------------------------------------------------------
  if (!file.exists(class_habitat_file))
    stop("Fichier LUT habitat introuvable : ", class_habitat_file)
  
  habitat_lut <- read.csv(class_habitat_file, sep = ";", stringsAsFactors = FALSE,
                          colClasses = c("code" = "integer", "label" = "character")) %>%
    dplyr::mutate(code = as.integer(code))
  
  names(habitat_r) <- "habitat_code"
  combined <- c(irg_stack, habitat_r)
  
  message("Conversion en table longue – patience…")
  combined_df <- as.data.frame(combined, xy = TRUE, cells = TRUE, na.rm = FALSE)
  
  # Assure la même classe pour la jointure
  combined_df$habitat_code <- as.integer(combined_df$habitat_code)
  
  long_df <- combined_df %>%
    tidyr::pivot_longer(
      cols      = starts_with("DOY"),
      names_to  = "DOY",
      values_to = "IRG"
    ) %>%
    dplyr::mutate(
      DOY = as.integer(stringr::str_extract(DOY, "[0-9]{3}"))
    ) %>%
    dplyr::left_join(habitat_lut, by = c("habitat_code" = "code")) %>%
    dplyr::mutate(
      habitat_label = factor(label, levels = habitat_lut$label)
    ) %>%
    dplyr::select(cell, x, y, habitat_code, habitat_label, DOY, IRG) %>%
    dplyr::arrange(cell, DOY)
  
  # ---------------------------------------------------------------------------
  # 6. Sauvegarde --------------------------------------------------------------
  message("Sauvegarde du tableau : ", output_IRG_by_habitat)
  # Utilise "gz" (plus universel) ; changez en "xz" si votre R le gère
  # Utilise "gzip" (ou TRUE) – compatible sur toutes les installations R
  saveRDS(long_df, output_IRG_by_habitat, compress = "gzip")
  
  message("Terminé ! Fichier créé : ", basename(output_IRG_by_habitat))
  
  
  
  
  data <- readRDS(output_IRG_by_habitat)
  
  

lut <- read_csv2(file.path(raster_dir, "class_habitat.csv"),
                 col_types = "ic",               # int ; char
                 locale    = locale(encoding = "latin1"))  # ← change en "UTF-8" si besoin




# ───────────────────────────────────────────────
# 1. Nettoyage  & libellés d’habitat
# ───────────────────────────────────────────────
plot_df <- data %>% 
  filter(!is.na(IRG),                      # enlève les valeurs IRG manquantes
         !is.na(habitat_code),             # enlève lignes sans code
         habitat_code != 0) %>%            # code 0 = fond / nodata
  left_join(lut, by = c("habitat_code" = "code")) %>% 
  mutate(
    habitat_label = coalesce(label, paste0("code_", habitat_code)),
    habitat_label = factor(habitat_label)          # garde l’ordre tel quel
  ) %>% 
  select(habitat_label, DOY, IRG)

# ───────────────────────────────────────────────
# 2. Résumé jour × habitat  (médiane + IC 95 %)
# ───────────────────────────────────────────────
summary_IRG <- plot_df %>% 
  group_by(habitat_label, DOY) %>% 
  summarise(
    median_IRG = median(IRG),
    low95      = quantile(IRG, .025),
    high95     = quantile(IRG, .975),
    .groups    = "drop"
  )

# Limites Y communes
y_lim <- range(summary_IRG$low95, summary_IRG$high95, na.rm = TRUE)

# ───────────────────────────────────────────────
# 3. Graphique
# ───────────────────────────────────────────────
p <- ggplot(summary_IRG, aes(DOY)) +
  geom_ribbon(aes(ymin = low95, ymax = high95),
              fill = "grey70", alpha = .35) +
  geom_line(aes(y = median_IRG),
            colour = "steelblue", linewidth = .6) +
  facet_wrap(~ habitat_label, ncol = 4) +           # ← ajuste ncol
  scale_x_continuous(breaks = seq(0, 360, 30),
                     minor_breaks = seq(0, 360, 10)) +
  scale_y_continuous(limits = y_lim) +
  labs(x = "Jour de l'année (DOY)",
       y = "IRG (médiane)",
       title   = "",
       caption = "") +
  theme_bw(base_size = 9) +
  theme(
    strip.text       = element_text(size = 9, face = "bold"),
    plot.title       = element_text(size = 12, face = "bold", hjust = .5),
    axis.title       = element_text(size = 10, face = "bold"),
    axis.text.x      = element_text(size = 7),
    axis.text.y      = element_text(size = 7),
    panel.grid.minor = element_line(linewidth = .2, linetype = "dotted")
  )

print(p)



}

if (TRUE){
  
  ### 4.2 Préparation du jeu de données phéno ###
  #---------------------------------------------#
  
  # Description : 
  # Préparation du jeu pour la phénologie; le jeu de données est issue du 
  # Script ressource : Plant_Phenology_Index_processing, et identifie les 4 phases
  # Phéno a partir du PPI puis de l'IRG
  
  # Fonction 1, correspond au jeu de données pour obtenir le pouventage du troupeau 
  # présent par état de végétation (stade phéno : pousse, maturé, dépérissant, senescent)
  
  # Nécéssite : 
  # - Le rds du chargement par jour (part 4, script 1)
  # - La phénologie : IRG (catégorisé selon 4 classe : 1 = pousse, 2 = plateau,
  #   3 = dépérissement, 4 = sénéssence) par jour et pixel de 10 mètres
  
  
  # LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_traitement_smod.R"))
  
  # ENTREE
  #Dossier contenant les sous dossier des chargement
  case_flock_file = file.path(output_dir, "4. Chargements_Calcules")
  #Dossier contenant les fichiers du tot de chargement
  case_flock_alpage_file = file.path(case_flock_file,paste0(YEAR,"_",alpage))
  
  # Un .RDS par alpage contenant les charges journalières
  daily_rds_file = file.path(case_flock_alpage_file, paste0("by_day_and_state_",YEAR,"_",alpage,".rds"))
  
  
  #Dossier de phénologie 
  case_phenologie = file.path(raster_dir, "Phenologie")
  if (!dir.exists(case_phenologie)) {
    dir.create(case_phenologie, recursive = TRUE)
  }
  
  # Un .TIF de l'IRG catégorisé par jour et pixel de 10 mètre
  pheno_tif_file = file.path(case_phenologie, paste0("Phenology_Phase_",alpage,"_",YEAR,".tif"))
  
  
  # SORTIE
  # Dossiers de sortie
  out_dir <- file.path(output_dir, "9. Analysis_Phenology", "data_phenology")
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  
  out_files <- list(
    pousse       = file.path(out_dir, paste0("Pourcentage_presence_stade_croissance_",   alpage,"_",YEAR,".tif")),
    plateau      = file.path(out_dir, paste0("Pourcentage_presence_stade_plateau_",      alpage,"_",YEAR,".tif")),
    deperiss    = file.path(out_dir, paste0("Pourcentage_presence_stade_deperiss_",    alpage,"_",YEAR,".tif")),
    senescence   = file.path(out_dir, paste0("Pourcentage_presence_stade_senescence_",  alpage,"_",YEAR,".tif"))
  )
  
  # CODE
  rast_pourcentage_troupeau_phase_pheno(daily_rds_file, pheno_tif_file, out_file)

}

if (TRUE) {
  
  ### 4.2 Préparation du jeu de données delta évolution en fonction IRGmax ###
  #--------------------------------------------------------------------------#
  
  # Description : 
  # Préparation du jeu pour le delta d'utilisation en fonction du max IRG;
  
  
  
  # LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_phenology.R"))
  
  
  
  # ENTREE
  #Dossier contenant les sous dossier des chargement
  case_flock_file = file.path(output_dir, "4. Chargements_Calcules")
  #Dossier contenant les fichiers du tot de chargement
  case_flock_alpage_file = file.path(case_flock_file,paste0(YEAR,"_",alpage))
  
  # Un .RDS par alpage contenant les charges journalières et par parc de nuit
  daily_rds_file = file.path(case_flock_alpage_file, paste0("by_park_day_and_state_transition_filtered_",YEAR,"_",alpage,".rds"))
  
  
  #Dossier de phénologie 
  case_phenologie = file.path(raster_dir, "Phenologie")
  if (!dir.exists(case_phenologie)) {
    dir.create(case_phenologie, recursive = TRUE)
  }
  
  # Un .TIF de l'IRG catégorisé par Stade Phénologique par jour et pixel de 10 mètres
  pheno_tif_file = file.path(case_phenologie, paste0("Phenology_Phase_",alpage,"_",YEAR,".tif"))
  
  # Un .TIF du delta de l'IRGmax en jour par jour et par pixel de 10 mètres
  delta_max_tif_file = file.path(case_phenologie, paste0("DeltaDay_from_MAXD_",alpage,"_",YEAR,".tif"))
  
  # SORTIE
  out_dir <- file.path(output_dir, "9. Analysis_Phenology", "data_IRG")
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  # Un .csv avec pour chaque pixel un valeur IRG par habitat
  output_IRG_by_habitat <- file.path(out_dir, paste0("dataset_delta_IRGmax",YEAR,"_",alpage,".rds"))
  
  
  # CODE
  
  ## ──────────────────────────────────────────────────────────────────────────────
  ##  4.2  Construction du data-set « Delta–IRG & Phéno & Charge »
  ## ──────────────────────────────────────────────────────────────────────────────
  library(terra)      # rasters
  library(data.table) # manip ultra-rapide
  library(glue)
  library(pbapply)    # barre de progression
  
  build_dataset_pheno_vec(
    YEAR                  = YEAR,
    alpage                = alpage,
    pheno_tif_file        = pheno_tif_file,
    delta_max_tif_file    = delta_max_tif_file,
    daily_rds_file        = daily_rds_file,
    output_IRG_by_habitat = output_IRG_by_habitat,
    DOY_range             = 121:334,
    chunk_size            = 30
  )
  
  
  
  

  
  
  
  
  ###############################################################################
  # VERSION 1 Δ-IRGmax quotidien : ruban = quantiles 2.5–97.5 %                          #
  ###############################################################################
  
  library(data.table)
  library(ggplot2)
  library(scales)
  library(glue)
  
  dt_2022 <- readRDS(file.path(out_dir, paste0("dataset_delta_IRGmax2022_",alpage,".rds")))
  dt_2023 <- readRDS(file.path(out_dir, paste0("dataset_delta_IRGmax2023_",alpage,".rds")))
  
  
  ## 1. Pixels utilisés (charge ≥ 1) -------------------------------------------
  dt_used <- dt[!is.na(charge) & charge >= 30]
  
  ## 2. Statistiques journalières ----------------------------------------------
  daily <- dt_used[
    , .(
      mean_delta = mean(delta_day_IRG_max, na.rm = TRUE),
      q2_5       = quantile(delta_day_IRG_max, 0.05, na.rm = TRUE),
      q97_5      = quantile(delta_day_IRG_max, 0.95, na.rm = TRUE),
      mean_pheno = mean(pheno_stage,           na.rm = TRUE)
    ),
    by = doy][order(doy)]
  
  ## 3. Palette gradient --------------------------------------------------------
  pheno_cols <- c("#008000", "#55C900", "#E6E600", "#FF7F00", "#4B2200")
  
  ## 4. Graphique ---------------------------------------------------------------
  ggplot(daily, aes(doy, mean_delta)) +
    geom_ribbon(aes(ymin = q2_5, ymax = q97_5),
                fill = "grey70", alpha = 0.3) +
    geom_line(aes(colour = mean_pheno),
              linewidth = 1.8, lineend = "round") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
    scale_colour_gradientn(
      "Stade phéno moyen",
      colours = pheno_cols,
      values  = rescale(c(1, 1.7, 2.5, 3.3, 4)),
      limits  = c(1, 4),
      breaks  = c(1, 2, 3, 4),
      labels  = c("Pousse", "Plateau", "Déclin", "Sénescence")
    ) +
    scale_x_continuous("Jour de l’année (DOY)",
                       breaks = seq(min(daily$doy),
                                    max(daily$doy), by = 10)) +
    labs(
      y = "Δ IRGmax (Day)",
      title = glue(""),
      subtitle = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  
  
  
  ### VERSION 2 ###
  
  library(data.table)
  library(ggplot2)
  library(scales)
  library(glue)
  
  
  
  
  
  
  
  
  # ---------------- Params “look & feel” ----------------
  smooth_k   <- 2      # fenêtre de lissage (jours) : 7–11 marche bien
  charge_min <- 20      # seuil de pixels utilisés
  pal_year   <- c(`2022`="#c0ccdf", `2023`="#c77f64")  # bleu & orange propres
  
  

  
  
  # --- Chargement ---------------------------------------------------------------
  dt_2022 <- as.data.table(readRDS(file.path(out_dir, paste0("dataset_delta_IRGmax2022_", alpage, ".rds"))))
  dt_2023 <- as.data.table(readRDS(file.path(out_dir, paste0("dataset_delta_IRGmax2023_", alpage, ".rds"))))
  
  dt_2022[, year := 2022L]
  dt_2023[, year := 2023L]
  
  dt <- rbindlist(list(dt_2022, dt_2023), use.names = TRUE, fill = TRUE)
  
  # --- Filtre pixels utilisés ---------------------------------------------------
  # (tu peux ajuster le seuil si besoin)
  dt_used <- dt[!is.na(charge) & charge >= 30]
  
  # --- Statistiques journalières par année -------------------------------------
  # Stat journalières par année (Δ moyen + IQR)
  daily <- dt_used[
    , .(
      mean_delta = mean(delta_day_IRG_max, na.rm = TRUE),
      q25        = quantile(delta_day_IRG_max, 0.1, na.rm = TRUE),
      q75        = quantile(delta_day_IRG_max, 0.9, na.rm = TRUE),
      n          = .N
    ),
    by = .(year, doy)
  ][order(year, doy)]
  
  # Lissage (moyenne mobile centrée)
  daily[, mean_smooth := frollmean(mean_delta, n = smooth_k, align = "center", na.rm = TRUE), by = year]
  daily[, q25_smooth  := frollmean(q25,        n = smooth_k, align = "center", na.rm = TRUE), by = year]
  daily[, q75_smooth  := frollmean(q75,        n = smooth_k, align = "center", na.rm = TRUE), by = year]
  
 
  
  
  
  # Breaks/labels Y robustes (assure des valeurs > 0 en plus de 0)
  y_breaks2 <- sort(unique(pretty(c(y_lower, y_upper), n = 6)))
  if (!0 %in% y_breaks2) y_breaks2 <- sort(c(y_breaks2, 0))
  
  lab_y2 <- function(b) {
    out <- scales::number(b, accuracy = 1)
    out[b == 0] <- "MAXV"
    out
  }
  
  # --- Plot --------------------------------------------------------------------
  p <- ggplot(daily, aes(x = doy, group = factor(year))) +
    geom_ribbon(aes(ymin = q25_smooth, ymax = q75_smooth, fill = factor(year)),
                alpha = 0.12, colour = NA, na.rm = TRUE) +
    geom_line(aes(y = mean_smooth, colour = factor(year)),
              linewidth = 1.6, lineend = "round", na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1, colour = "grey20") +
    scale_color_manual(
      name   = "Grazing resource used by the herd :",
      values = pal_year,
      labels = c("2022", "2023")
    ) +
    scale_fill_manual(values = pal_year, guide = "none") +
    scale_x_continuous("Day of year (DOY)",
                       breaks = pretty(daily$doy, n = 8),
                       expand = expansion(mult = c(0.005, 0.02))) +
    scale_y_continuous("Delta days from vegetation peak (MAXV)",
                       limits = c(y_lower, y_upper),
                       breaks = y_breaks2,          # << force des graduations Y
                       labels = lab_y2,             # << 0 -> MAXV, le reste numérique
                       expand = expansion(mult = c(0.02, 0.08))) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position      = c(0.985, 0.98),
      legend.justification = c(1, 1),
      legend.direction     = "vertical",
      legend.background    = element_rect(fill = scales::alpha("white", 0.85), colour = "grey80"),
      legend.key.height    = unit(10, "pt"),
      legend.title         = element_text(size = 10, face = "bold"),
      legend.text          = element_text(size = 10),
      
      # Axes + ticks bien visibles
      axis.text.x          = element_text(size = 10),
      axis.text.y          = element_text(size = 10),
      axis.ticks.x         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.y         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length    = unit(3, "pt"),
      
      panel.grid.minor     = element_blank(),
      panel.grid.major.x   = element_line(linewidth = 0.25),
      panel.grid.major.y   = element_line(linewidth = 0.25),
      
      # Cadre panneau + cadre global
      panel.border         = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      
      plot.title           = element_blank(),
      plot.margin          = margin(4, 8, 2, 8)
    )
  
  print(p)
  
  
  
  # Export as a low-height banner (nearly A4 width)
  ggsave(file.path(out_dir, glue("delta_IRGmax_{alpage}_2022_2023_banner_en.svg")),
         p, width = 11.2, height = 5.5, dpi = 300)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ### VERSION 2 — corrigée & robuste ###
  
  library(data.table)
  library(ggplot2)
  library(scales)
  library(glue)
  
  # ---------------- Params “look & feel” ----------------
  smooth_k   <- 2            # impair pour un vrai centrage; 7–11 marche bien
  charge_min <- 20           # seuil de pixels utilisés (cohérent partout)
  pal_year   <- c(`2022`="#c0ccdf", `2023`="#c77f64")  # orange & bleu doux
  
  # Forcer smooth_k impair si besoin
  if (smooth_k %% 2 == 0) smooth_k <- smooth_k + 1L
  
  # --- Chargement ---------------------------------------------------------------
  dt_2022 <- as.data.table(readRDS(file.path(out_dir, paste0("dataset_delta_IRGmax2022_", alpage, ".rds"))))
  dt_2023 <- as.data.table(readRDS(file.path(out_dir, paste0("dataset_delta_IRGmax2023_", alpage, ".rds"))))
  
  dt_2022[, year := 2022L]
  dt_2023[, year := 2023L]
  
  dt <- rbindlist(list(dt_2022, dt_2023), use.names = TRUE, fill = TRUE)
  
  # --- Filtre pixels utilisés ---------------------------------------------------
  dt_used <- dt[!is.na(charge) & charge >= charge_min]
  
  # --- Statistiques journalières par année -------------------------------------
  # Δ moyen + quantiles 10/90, avec garde-fous si peu d'observations
  n_min_grp <- charge_min
  daily <- dt_used[
    , .(
      mean_delta = mean(delta_day_IRG_max, na.rm = TRUE),
      q10        = if (.N >= n_min_grp) quantile(delta_day_IRG_max, 0.10, na.rm = TRUE) else NA_real_,
      q90        = if (.N >= n_min_grp) quantile(delta_day_IRG_max, 0.90, na.rm = TRUE) else NA_real_,
      n          = .N
    ),
    by = .(year, doy)
  ][order(year, doy)]
  
  # --- Lissage (moyenne mobile centrée) ----------------------------------------
  daily[, mean_smooth := frollmean(mean_delta, n = smooth_k, align = "center", na.rm = TRUE), by = year]
  daily[, q10_smooth  := frollmean(q10,        n = smooth_k, align = "center", na.rm = TRUE), by = year]
  daily[, q90_smooth  := frollmean(q90,        n = smooth_k, align = "center", na.rm = TRUE), by = year]
  
  # --- Plage Y & ticks robustes -------------------------------------------------
  yrng    <- range(c(daily$mean_smooth, daily$q10_smooth, daily$q90_smooth), na.rm = TRUE)
  if (!all(is.finite(yrng))) yrng <- c(-1, 1)  # fallback si tout NA
  pad     <- diff(yrng) * 0.06
  y_lower <- floor(yrng[1] - pad)
  y_upper <- ceiling(yrng[2] + pad)
  
  y_breaks2 <- pretty(c(y_lower, y_upper), n = 6)
  if (!0 %in% y_breaks2) y_breaks2 <- sort(c(y_breaks2, 0))
  
  lab_y2 <- function(b) {
    out <- number(b, accuracy = 1)
    out[b == 0] <- "MAXV"
    out
  }
  
  # --- Plot --------------------------------------------------------------------
  # --- Plot (avec -20 sur Y + label "Delta") -----------------------------------
  # on force -20 dans les limites et dans les graduations
  y_lower2  <- min(y_lower, -20)
  y_breaks2 <- sort(unique(c(y_breaks2, -20, 0)))
  
  p <- ggplot(daily[n >= n_min_grp], aes(x = doy, group = factor(year))) +
    geom_ribbon(aes(ymin = q10_smooth, ymax = q90_smooth, fill = factor(year)),
                alpha = 0.18, colour = NA, na.rm = TRUE) +
    geom_line(aes(y = mean_smooth, colour = factor(year)),
              linewidth = 1.4, lineend = "round", na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "longdash", linewidth = 0.9, colour = "grey20") +
    scale_color_manual(
      name   = "Grazing resource use by the herd :",
      values = pal_year,
      labels = c("2022", "2023")
    ) +
    scale_fill_manual(values = pal_year, guide = "none") +
    scale_x_continuous("Day of year (DOY)",
                       breaks = pretty(unique(daily$doy), n = 8),
                       expand = expansion(mult = c(0.005, 0.02))) +
    scale_y_continuous("Delta days from vegetation peak (MAXV)",
                       limits = c(y_lower2, y_upper),
                       breaks = y_breaks2,
                       labels = lab_y2,         # 0 reste étiqueté "MAXV"
                       expand = expansion(mult = c(0.02, 0.08))) +
    theme_minimal(base_size = 15) +
    theme(
      legend.position      = c(0.985, 0.98),
      legend.justification = c(1, 1),
      legend.direction     = "vertical",
      legend.background    = element_rect(fill = scales::alpha("white", 0.85), colour = "grey80"),
      legend.key.height    = unit(10, "pt"),
      legend.title         = element_text(size = 12, face = "bold"),
      legend.text          = element_text(size = 12),
      axis.text.x          = element_text(size = 12),
      axis.text.y          = element_text(size = 12),
      axis.ticks.x         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.y         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length    = unit(3, "pt"),
      panel.grid.minor     = element_blank(),
      panel.grid.major.x   = element_line(linewidth = 0.25),
      panel.grid.major.y   = element_line(linewidth = 0.25),
      panel.border         = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      plot.title           = element_blank(),
      plot.margin          = margin(4, 8, 2, 8)
    )
  
  print(p)
  
  
  # --- Export -------------------------------------------------------------------
  outfile_png <- file.path(out_dir, glue("delta_IRGmax_{alpage}_2022_2023_banner_en.svg"))
  ggsave(outfile_png, p, width = 6.2, height = 5.5, dpi = 300)
  cat("✅ Figure écrite :", outfile_png, "\n")
  
  
  
  
  
  
  
  
  
  
  
  ###############################################################################
  #  Δ-IRGmax quotidien – ruban (2.5–97.5 %) lissé par moyenne mobile 7 jours   #
  ###############################################################################
  
  library(data.table)
  library(zoo)        # rollapply
  library(ggplot2)
  library(scales)
  library(glue)
  
  dt <- readRDS(output_IRG_by_habitat)
  
  ## 1. Pixels “utilisés” (charge ≥ 1) ------------------------------------------
  dt_used <- dt[!is.na(charge) & charge >= 10]
  
  ## 2. Statistiques journalières (brutes) --------------------------------------
  daily <- dt_used[
    , .(
      mean_delta = mean(delta_day_IRG_max, na.rm = TRUE),
      q2_5       = quantile(delta_day_IRG_max, 0.025, na.rm = TRUE),
      q97_5      = quantile(delta_day_IRG_max, 0.975, na.rm = TRUE),
      mean_pheno = mean(pheno_stage,           na.rm = TRUE)
    ),
    by = doy][order(doy)]
  
  ## 3. Lissage *uniquement* de l’enveloppe 2.5–97.5 % --------------------------
  k <- 7        # largeur fenêtre glissante (jours) — ajuste si besoin
  
  daily[, q2_5_smooth  := rollapply(q2_5,  k, mean, align = "center", fill = NA)]
  daily[, q97_5_smooth := rollapply(q97_5, k, mean, align = "center", fill = NA)]
  
  ## 4. Palette gradient continue ----------------------------------------------
  pheno_cols <- c(
    "#008000", "#55C900", "#E6E600", "#FF7F00", "#4B2200"
  )
  
  ## 5. Graphique ---------------------------------------------------------------
  ggplot(daily, aes(doy)) +
    # ruban lissé (enveloppe)
    geom_ribbon(aes(ymin = q2_5_smooth, ymax = q97_5_smooth),
                fill = "grey70", alpha = 0.35, na.rm = TRUE) +
    # courbe journalière brute (pas lissée)
    geom_line(aes(y = mean_delta, colour = mean_pheno),
              linewidth = 1.6, lineend = "round", na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
    scale_colour_gradientn(
      "Stade phéno moyen",
      colours = pheno_cols,
      values  = rescale(c(1, 1.7, 2.5, 3.3, 4)),
      limits  = c(1, 4),
      breaks  = 1:4,
      labels  = c("Pousse", "Plateau", "Déclin", "Sénescence")
    ) +
    scale_x_continuous(
      "Jour de l’année (DOY)",
      breaks = seq(min(daily$doy), max(daily$doy), by = 10)
    ) +
    labs(
      y = "Δ jours (moyenne ± quantiles lissés 2.5–97.5 %)",
      title = glue("Δ-IRG_max quotidien – {alpage} {YEAR}"),
      subtitle = "Pixels utilisés : charge ≥ 1  |  Couleur = stade phéno moyen (continu)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  library(data.table)
  
  ## 1. Charge ton .rds  (juste la partie "charges")
  charges <- as.data.table(readRDS(daily_rds_file))
  
  setnames(charges,
           c("x","y","day","Charge","parc"),
           c("x","y","doy","charge","parc_nuit"))
  
  # éventuel filtrage :
  charges <- charges[state == "Paturage"]
  
  # projection sur la grille (reproduis la même ligne que dans la fonction)
  library(terra)
  pheno_stack <- rast(pheno_tif_file)
  charges[, cell := terra::cellFromXY(pheno_stack, cbind(x, y))]
  charges <- charges[!is.na(cell)]
  charges[, doy := as.integer(doy)]
  
  # agrégation stricte
  charges <- charges[, .(
    charge    = sum(charge, na.rm = TRUE),
    parc_nuit = paste(unique(parc_nuit), collapse = ";")
  ), by = .(cell, doy)]
  
  ## 2. Vérifie l’unicité
  dup_n <- anyDuplicated(charges, by = c("cell", "doy"))
  dup_n
  
  
  
  
  
  
  
  
  
}

if (TRUE) {
  
  ### 4.4 Préparation du jeu de données impacte pasto sur le NDVI et IRGmax ###
  #--------------------------------------------------------------------------#
  
  # Description : 
  # Préparation du dataset pour evalué si pour un meme habitats dans les meme 
  # contrainte stationnelle (Enneigement, DAH) ; si pour 3 classe de chargement :
  # 
  # - Fort (chargement > 500)
  # - Moyen (chargement = 200)
  # - Faible (chargement < 10)
  #
  # On oberve un impacte sur la strcture et le fonctionnement des communautés végétal a travers
  # le NDVI et l'IRG
  #
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
  
  
  
  
  # LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_phenology.R"))
  
  # ENTREE
  #Dossier de phénologie 
  case_phenologie = file.path(raster_dir, "Phenologie")
  if (!dir.exists(case_phenologie)) {
    dir.create(case_phenologie, recursive = TRUE)
  }
  
  # Un .TIF de l'IRG catégorisé par jour et pixel de 20 mètre
  IRG_tif_file = file.path(case_phenologie, paste0("IRG_season_",alpage,"_",YEAR,".tif"))
 
  # Un .TIF de l'IRG catégorisé par jour et pixel de 20 mètre
  NDVI_tif_file = file.path(case_phenologie, paste0("NDVI_season_",alpage,"_",YEAR,".tif"))
  
  
  # Un dossier contenant carte de végétation
  carto_file = file.path(raster_dir, "Classifications_fusion_ColorIndexed_sc1_landforms_mnh.tif")
 
  
  # Un dossier contenant les ratsers des Unités Pastorales (UP)
  case_UP_file = file.path(raster_dir, "UP")
  # Un .SHP avec les Unités pastorales UP
  UP_file = file.path(case_UP_file,paste0 ("UP_",alpage,".shp"))
  
  
  #Dossier contenant les sous dossier des chargement
  case_flock_file = file.path(output_dir, "4. Chargements_Calcules")
  #Dossier contenant les fichiers du tot de chargement
  case_flock_alpage_file = file.path(case_flock_file,paste0(YEAR,"_",alpage))
  
  # Un .RDS par alpage contenant les charges journalières et par parc de nuit
  total_rds_file = file.path(case_flock_alpage_file, paste0("total_",YEAR,"_",alpage,".rds"))
  
  
  # Un .TIF du SMOD Médian sur 10 ans 
  clim_case <- file.path(output_dir, "8. Analysis_Climate")
  SMOD_case <- file.path(clim_case, "SMOD")
  SMOD_tif_file <- file.path(SMOD_case, paste0("Fsca_",alpage,".tif"))
  
  # Un .TIF du DAH de l'alpage
  Alti_case <- file.path(raster_dir, "Alti")
  DAH_tif_file <- file.path(Alti_case, paste0("DAH_1_",alpage,".tif"))
  
  lut_csv        <- file.path(raster_dir, "class_habitat.csv")
  
  
  # Un .TIF avec les extra du PPI
  EXTRA_tif = file.path(case_phenologie, paste0("PPI_extra_",alpage,"_",YEAR,".tif"))
 
  
  ## SORTIE
  # Dossier de sortie
  out_dataset       <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_",alpage))
  dir.create(out_dataset, recursive=TRUE, showWarnings=FALSE)
  # Un .fst contetant l'ensemble du dataset 
  out_rds <- file.path(out_dataset, paste0("dataset_pheno_LONG_", YEAR, "_", alpage, ".rds"))
  
  ## CODE
  
  

  
  
  ## FONCTION 1 :
  ## Création du dataset (X, Y, DAH, FSCA, NDVI, IRG, HABITAT, CHARGEMENT)
  build_pheno_dataset(
    alpage   = "Viso",
    YEAR     = 2022,
    irg_tif  = IRG_tif_file,
    ndvi_tif = NDVI_tif_file,
    habitat_tif = carto_file,
    dah_tif     = DAH_tif_file,
    fsca_tif    = SMOD_tif_file,
    extra_tif   = EXTRA_tif,            # ← NEW
    charge_rds  = total_rds_file,
    lut_habitat_csv = lut_csv,
    up_shape  = UP_file,
    output_rds = out_rds
  )
  
  
  ## FONCTION 2 :
  ## Controle du jeu de données (PDF dans le outpout du dataset)
  # Plot de la zone d'étude 
  # Histogramme de : DAH, FSCA, NDVI, IRG
  diag_dataset(
    rds_path = out_rds,
    up_shape = UP_file,   # ou NULL pour bbox seule
    view     = FALSE     # pour voir les graphes dans RStudio
  )
  
  
  
  
  
  ## FONCTION 3 : 
  # Séléction des zones ayant les meme conditions de DAH, FSCA, et Habitat. 
  # Et identification dans ces zones d'un gradient de paturage : 
  # - Faible (inf. 50)
  # - Moyen (a faire)
  # - Fort (sup. 300)
  
  ## ENTREE
  template_tif <- file.path(case_phenologie, paste0("NDVI_season_",alpage,"_",YEAR,".tif"))
  
  ## SORTIE
  # Un . TIF de la zone d'étude séléctioné par la fonction
  ras_out <- file.path(out_dataset, "zonesTraitement.tif")
  # Un . SHP de la zone d'étude séléctionné par la fonction
  shp_out <- file.path(out_dataset, "zonesTraitement.shp")
  
  # Le DataSet de sortie apres séléction des meme 
  
  csv_out <- file.path(out_dataset, paste0("dataset_filtered_",alpage,"_",YEAR,".csv"))
  
  
  
  
  
  
  
  # Définition de l'année d'analyse
  YEAR <- 2022
  TYPE <- "catlog" #Type de données d'entrée (CATLOG, OFB )
  alpage <- "Viso"
  
  
  
  alpages <- c("Cayolle","Viso","Sanguiniere")
  
  
  
  
  
  
  
  
  
  detect_zones_pheno(
    rds_path = merge_rds_path,
    up_shape        = NULL,
    thr_high        = 200,
    thr_low         =  50,
    eps_m           = 10000,
    minpts          = 10,
    hab_exclus      = c(NA, "Formations minérales"),
    ras_template    = NULL,
    csv_out,
    ras_out,
    shp_out
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  source(file.path(functions_dir, "Functions_phenology.R"))
  
  run_detect_zones(
    alpages    = "Cayolle",  # ou un seul
    YEAR       = 2022,
    output_dir = output_dir,
    delta_pct  = 0.10,
    minpts     = 1
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    #' Génération de graphiques NDVI par traitement pour chaque habitat
  #'
  #' Lecture directe du CSV zonesTraitement.csv.gz (qui contient DOY, IRG, NDVI)
  #' et production d'un graphique par habitat montrant la médiane et les intervalles
  #' (20%-80%) de NDVI au cours du temps pour chaque traitement.
  #'
  # ────────────────────────────────────────────────────────────────────────────────
  # 0. Packages --------------------------------------------------------------------
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
  })
  
  # ────────────────────────────────────────────────────────────────────────────────
  # 1. Lecture des données ---------------------------------------------------------
  csv_zones <- "C:/Users/massocam/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/outputs/10. NDVI Effect/dataset_filtered_Cayolle_Viso_Sanguiniere_2022.csv"
  dt_zones  <- read.csv(csv_zones)
  
  req <- c("habitat_label", "traitement", "DOY", "NDVI", "IRG")
  stopifnot(all(req %in% names(dt_zones)))
  
  # ────────────────────────────────────────────────────────────────────────────────
  # 2. Paramètres généraux ---------------------------------------------------------
  habitats <- unique(dt_zones$habitat_label)
  cols     <- c(Fort = "steelblue", Faible = "khaki3")
  labs_trt <- c(Fort = "high grazing intensity",
                Faible = "low grazing intensity")
  
  plots <- vector("list", length(habitats))
  names(plots) <- habitats
  
  # ────────────────────────────────────────────────────────────────────────────────
  # 3. Boucle par habitat ----------------------------------------------------------
  for (hab in habitats) {
    
    df_h <- dt_zones %>% filter(habitat_label == hab)
    
    # ── 3.1 NDVI : médiane & bande 20-80 %
    ndvi_stats <- df_h %>%
      group_by(DOY, traitement) %>%
      summarise(
        med_ndvi   = median(NDVI, na.rm = TRUE),
        low80_ndvi = quantile(NDVI, 0.20, na.rm = TRUE),
        hi80_ndvi  = quantile(NDVI, 0.80, na.rm = TRUE),
        .groups = "drop"
      )
    
    # ── 3.2 IRG : jour & valeur du pic
    irg_peaks <- df_h %>%
      group_by(DOY, traitement) %>%
      summarise(med_irg = median(IRG, na.rm = TRUE), .groups = "drop") %>%
      group_by(traitement) %>%
      slice_max(order_by = med_irg, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      rename(DOY_peak = DOY, irg_max = med_irg)
    
    # ── 3.3 Bornes NDVI et padding
    y_min <- min(ndvi_stats$low80_ndvi, na.rm = TRUE)
    y_max <- max(ndvi_stats$hi80_ndvi , na.rm = TRUE)
    pad   <- 0.10 * (y_max - y_min)          # marge pour respirer
    
    
    
    # ── 3.4 Mise à l’échelle IRG → NDVI (harmonisée)
    target_frac <- 0.90                       # 90 % de la hauteur NDVI
    baseline_y  <- y_min                      # pied des barres
    peak_y      <- baseline_y + (y_max - y_min) * target_frac
    sf          <- (peak_y - baseline_y) / max(irg_peaks$irg_max)
    
    # Limite supérieure : on laisse une marge confortable (30 % du pad)
    y_max_all   <- max(y_max, peak_y) + pad * 0.30   # ← nouvelle ligne
    
    # ── 3.5 Construction du graphique
    p <- ggplot() +
      # NDVI (ruban + médiane)
      geom_ribbon(
        data = ndvi_stats,
        aes(x = DOY, ymin = low80_ndvi, ymax = hi80_ndvi, fill = traitement),
        alpha = .25, colour = NA
      ) +
      geom_line(
        data = ndvi_stats,
        aes(x = DOY, y = med_ndvi, colour = traitement),
        size = 1
      ) +
      
      # Barre verticale IRGmax
      geom_segment(
        data = irg_peaks,
        aes(
          x     = DOY_peak,
          xend  = DOY_peak,
          y     = baseline_y,
          yend  = irg_max * sf,
          colour = traitement
        ),
        size = 2
      ) +
      
      # Trait horizontal (optionnel) matérialisant IRGmax
      geom_hline(
        data = irg_peaks,
        aes(yintercept = irg_max * sf, colour = traitement),
        linetype = "dashed", size = .8
      ) +
      
      # Axes
      scale_x_continuous(
        breaks       = seq(60, 330, 30),
        minor_breaks = seq(60, 330, 10)
      ) +
      scale_y_continuous(
        name     = "Plant Phenology Index (PPI)",
        limits   = c(baseline_y, y_max_all),
        expand   = c(0, 0),
        sec.axis = sec_axis(~ . / sf, name = "Greenup rate (IRG)")
      ) +
      
      # Couleurs & légendes
      scale_colour_manual(values = cols, labels = labs_trt, name = NULL) +
      scale_fill_manual(  values = cols, labels = labs_trt, name = NULL) +
      
      labs(
        title = hab,
        x     = "Julian date"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor       = element_line(size = .2, linetype = "dotted"),
        panel.grid.major       = element_line(size = .3),
        legend.position.inside = c(.78, .92),
        legend.background      = element_blank()
      )
    
    plots[[hab]] <- p
  }
  
  # ────────────────────────────────────────────────────────────────────────────────
  # 4. Affichage -------------------------------------------------------------------
  for (p in plots) print(p)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ###############################################################################
  # 0. PARAMÈTRES ---------------------------------------------------------------
  ###############################################################################
  fst_file   <- "C:/Users/massocam/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/outputs/10. NDVI Effect/dataset_Cayolle/dataset_pheno_LONG_2022_Cayolle.rds"
  
  grid_tif   <- "maillage_20m.tif"        # GeoTIFF de la trame
  grid_shp   <- "maillage_20m.shp"        # shapefile de polygones
  res_m      <- 20                        # résolution du maillage (m)
  crs_code   <- "EPSG:2154"               # Lambert-93 (à adapter si besoin)
  
  ###############################################################################
  # 1. LIBRAIRIES ---------------------------------------------------------------
  ###############################################################################
  suppressPackageStartupMessages({
    library(data.table)
    library(fst)
    library(terra)
    library(sf)
  })
  
  ###############################################################################
  # 2. LECTURE & EMPRISE --------------------------------------------------------
  ###############################################################################
  dt <- as.data.table(readRDS(fst_file))
  
  # harmonise noms de colonnes → x / y
  for (p in list(c("x","y"), c("X","Y"), c("lon","lat"), c("longitude","latitude")))
    if (all(p %in% names(dt))) { setnames(dt, p, c("x","y")); break }
  
  if (!all(c("x","y") %in% names(dt)))
    stop("Colonnes de coordonnées introuvables dans le .fst")
  
  # emprise (min-max) des points
  ext_dt <- ext(range(dt$x), range(dt$y))
  
  ###############################################################################
  # 3. CRÉATION DU RASTER GÉNÉRIQUE --------------------------------------------
  ###############################################################################
  r <- rast(ext_dt, res = res_m, crs = crs_code)
  values(r) <- seq_len(ncell(r))          # chaque pixel reçoit son indice
  
  writeRaster(r, grid_tif, datatype = "INT4S", overwrite = TRUE)
  message("✓ GeoTIFF écrit : ", grid_tif, "  (", ncell(r), " pixels)")
  
  ###############################################################################
  # 4. POLYGONES DU MAILLAGE (OPTIONNEL) ----------------------------------------
  ###############################################################################
  polys <- as.polygons(r, dissolve = FALSE)
  names(polys) <- "cell"
  
  st_write(st_as_sf(polys), grid_shp, delete_dsn = TRUE, quiet = TRUE)
  message("✓ Shapefile écrit : ", grid_shp, "  (", nrow(polys), " polygones)")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  library(fst)
  library(data.table)
  
  # Charger le fichier .fst
  fst_file <- "C:/Users/massocam/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/outputs/10. NDVI Effect/dataset_Cayolle/dataset_pheno_LONG_2022_Cayolle.fst"
  dt <- as.data.table(read_fst(fst_file, as.data.table = TRUE))
  
  # Afficher les noms de colonnes pour être sûr
  print(names(dt))
  
  # Aperçu des colonnes x et y
  print(head(dt[, .(x, y)]))
  
  # Statistiques de base
  summary(dt[, .(x, y)])
  
  # Vérifier s’il y a des NA dans x ou y
  cat("Nombres de NA:\n")
  cat("x:", sum(is.na(dt$x)), "\n")
  cat("y:", sum(is.na(dt$y)), "\n")
  
  
  
  
  # 1. combien de cellules signalées dépassent la taille du raster ?
  range(sub$cell)               # min et max des indices à affecter
  ncell(gabarit)                # nombre total de cellules du raster
  
  # 2. nombre de coordonnées manquantes
  sum(is.na(xy[, "x"]) | is.na(xy[, "y"]))       # devrait être 0
  
  # 3. quels indices posent problème
  bad <- which(is.na(xy[,1]) | is.na(xy[,2]))
  head(pts[bad])           
  
  
  
  
  library(terra)
  
  gabarit <- rast(gabarit_tif, subds = 1)
  
  # • Emprise et CRS du raster
  print(ext(gabarit))
  
  crs(gabarit)
  
  # • Plage des coordonnées tabulaires
  dt <- data.table::as.data.table(fst::read_fst(in_sft))
  for(nm in list(c("x","y"), c("X","Y")))
    if(all(nm %in% names(dt))) data.table::setnames(dt, nm, c("x","y"))
  
  summary(dt[, .(x, y)])
  
  
  
  ## FONCTION 3 : 
  
  # LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_phenology.R"))
  in_sft  <- file.path(out_dataset,
                       sprintf("dataset_pheno_LONG_%d_%s.fst", YEAR, alpage))
  out_tab <- file.path(out_dataset,
                       sprintf("dataset_filtered_to_stationnelle_%d_%s.fst",
                               YEAR, alpage))
  out_tif <- file.path(out_dataset,
                       sprintf("Area_dataset_%d_%s.tif", YEAR, alpage))
  out_shp <- file.path(out_dataset,
                       sprintf("Area_dataset_%d_%s.shp", YEAR, alpage))
  
  gabarit_tif <- "raster/Phenologie/NDVI_season_Cayolle_2022.tif"
  
  create_zones_topo(
    fst_in    = in_sft,
    up_tif    = gabarit_tif,
    table_out = out_tab,
    ras_out   = out_tif,
    shp_out   = out_shp,
    win_perc  = 1,
    method    = "quantile",
    input_crs = "EPSG:3857"   # ← ajoute ce paramètre !
  )
  
  
  

  ## NEW :  : crée des « zones de traitement topographiques »
  ##  à partir d’un fichier .fst (colonnes : cell, x, y, Charge, DOY,
  ##  dah, fsca, habitat_label …).
  ##  
  ##  → Exporte :
  ##      • un CSV récapitulatif par pixel (zones_topo_habitat_chargement.csv.gz)
  ##      • un raster INT1U (zonesTraitement.tif) : 1 = Faible, 2 = Fort
  ##      • un shapefile ponctuel (zonesTraitement.shp)
  ## -----------------------------------------------------------------
  
  ###############################################################################
  # 0. LIBRAIRIES ---------------------------------------------------------------
  ###############################################################################
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(dbscan)
    library(terra)
    library(sf)
    library(ggplot2)
    library(fst)
    library(patchwork)
  })
  
  ###############################################################################
  # 1. PARAMÈTRES UTILISATEUR ---------------------------------------------------
  ###############################################################################
  fst_in <- "C:/Users/massocam/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/outputs/10. NDVI Effect/dataset_Cayolle/dataset_pheno_LONG_2022_Cayolle.rds"
  

  ras_out <- "outputs/9. Analysis_Phenology/data_joint/zonesTraitement.tif"
  shp_out <- "outputs/9. Analysis_Phenology/data_joint/zonesTraitement.shp"
  
  template_tif <- "raster/Phenologie/NDVI_season_Cayolle_2022.tif"  # raster de référence
  
  thr_high <- 300      # >300  → Fort
  thr_low  <-  50      # ≤50   → Faible
  
  eps_m  <- 1000       # rayon DBSCAN (m)
  minpts <- 1          # pixels min / noyau
  
  hab_exclus <- c(NA, "Formations minérales")
  
  ###############################################################################
  # 2. LECTURE & PRÉ‑TRAITEMENT -------------------------------------------------
  ###############################################################################
  message("Lecture du .fst …")
  
  dt <- as.data.table(readRDS(fst_in))
  
  # harmonisation coordonnées → x / y
  for (p in list(c("x","y"), c("x.x","y.x"), c("x.y","y.y"),
                 c("X","Y"), c("lon","lat"), c("longitude","latitude")))
    if (all(p %in% names(dt))) { setnames(dt, p, c("x","y")); break }
  
  req_cols <- c("cell","x","y","Charge","DOY","dah","fsca","habitat_label")
  if (!all(req_cols %in% names(dt)))
    stop("❌  Colonnes manquantes dans le .fst : ",
         paste(setdiff(req_cols, names(dt)), collapse = ", "))
  
  # variables statiques (premier jour de l'année)
  dt_first <- dt[DOY == min(DOY)]
  
  ###############################################################################
  # 3. HABITATS MAJEURS ---------------------------------------------------------
  ###############################################################################
  message("Détermination des 3 habitats les plus représentés …")
  
  top3 <- dt_first[ !habitat_label %in% hab_exclus,
                    .N, by = habitat_label ][order(-N)][1:3, habitat_label]
  
  ###############################################################################
  # 4. FONCTION : DÉTECTION DES ZONES HOMOGÈNES ---------------------------------
  ###############################################################################
  create_zones <- function(sub_dt) {
    med_dah  <- median(sub_dt$dah , na.rm = TRUE)
    med_fsca <- median(sub_dt$fsca, na.rm = TRUE)
    win_dah  <- 0.10 * abs(med_dah)
    win_fsca <- 0.10 * abs(med_fsca)
    
    cand <- sub_dt[dah  %between% c(med_dah - win_dah , med_dah + win_dah) &
                     fsca %between% c(med_fsca - win_fsca, med_fsca + win_fsca)]
    if (nrow(cand) < minpts) return(NULL)
    
    coords <- as.matrix(cand[, .(x, y)])
    cl <- dbscan(coords, eps = eps_m, minPts = minpts)
    cand$cluster <- cl$cluster
    cand[cluster > 0]
  }
  
  ###############################################################################
  # 5. TRAITEMENT PAR HABITAT ---------------------------------------------------
  ###############################################################################
  message("Recherche des zones de traitement (Faible / Fort) …")
  
  zones_all <- list()
  
  for (hab in top3) {
    
    message("• Habitat : ", hab)
    
    z <- create_zones(dt_first[habitat_label == hab])
    if (is.null(z)) { warning("  aucune zone détectée."); next }
    
    # classes Faible / Fort
    z[ , traitement := fcase(
      Charge >  thr_high, "Fort",
      Charge <= thr_low , "Faible",
      default = NA_character_)]
    z <- z[!is.na(traitement)]
    
    # clusters contenant les deux classes
    tab  <- z[ , .N, by = .(cluster, traitement)]
    wide <- dcast(tab, cluster ~ traitement, value.var = "N", fill = 0)
    if (!all(c("Fort","Faible") %in% names(wide))) next
    keep <- wide[Fort >= minpts & Faible >= minpts, cluster]
    
    z <- z[cluster %in% keep]
    if (!nrow(z)) next
    
    # renumérote 3 plus gros clusters
    top_cl <- z[ , .N, by = cluster][order(-N)][1:3,
                                                .(cluster, zone_id = seq_len(.N))]
    z <- merge(z, top_cl, by = "cluster", all.y = TRUE)
    
    z[ , habitat_label := hab]
    zones_all[[hab]] <- z[, .(cell, x, y, zone_id,
                              habitat_label, traitement,
                              dah, fsca)]
  }
  
  ###############################################################################
  # 6. EXPORT CSV ---------------------------------------------------------------
  ###############################################################################
  if (!length(zones_all)) stop("❌  aucune zone valide.")
  
  final <- rbindlist(zones_all)
  
  fwrite(final, csv_out)
  message("✓ CSV exporté : ", csv_out)
  
  ###############################################################################
  # 7. RASTER & SHAPEFILE -------------------------------------------------------
  ###############################################################################
  template <- rast(template_tif, subds = 1)
  zone_r   <- rast(template); values(zone_r) <- NA_integer_
  
  valid <- !is.na(final$cell)
  code_tab <- c(Faible = 1L, Fort = 2L)
  final[ , code := code_tab[traitement]]
  
  zone_r[ final$cell[valid] ] <- final$code[valid]
  writeRaster(zone_r, ras_out, datatype = "INT1U", overwrite = TRUE)
  message("✓ Raster exporté : ", ras_out, "  (1 Faible, 2 Fort)")
  
  xy  <- terra::xyFromCell(zone_r, final$cell[valid])
  sf_pts <- st_as_sf(data.frame(xy,
                                code = final$code[valid]),
                     coords = c("x","y"),
                     crs = crs(zone_r))
  st_write(sf_pts, shp_out, delete_dsn = TRUE, quiet = TRUE)
  message("✓ Shapefile exporté : ", shp_out, "  (champ 'code' = 1/2)")
  
  ###############################################################################
  # 8. BOXPLOTS & BARPLOTS ------------------------------------------------------
  ###############################################################################
  pix_area_ha <- prod(res(template)) / 1e4  # surface pixel en ha
  
  plot_dt <- final[!is.na(dah) & !is.na(fsca),
                   .(traitement, DAH = dah, FSCA = fsca)]
  
  pal <- c(Faible = "khaki3", Fort = "steelblue")
  base_theme <- theme_bw(base_size = 11) +
    theme(axis.title.x  = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "grey85", linewidth = .3),
          panel.grid.minor.y = element_blank())
  
  p_dah <- ggplot(plot_dt, aes(traitement, DAH, fill = traitement)) +
    geom_boxplot(width = .65, outlier.shape = 21, outlier.size = 1) +
    scale_fill_manual(values = pal) +
    labs(title = "Distribution DAH", y = "DAH") +
    base_theme
  
  p_fsca <- ggplot(plot_dt, aes(traitement, FSCA, fill = traitement)) +
    geom_boxplot(width = .65, outlier.shape = 21, outlier.size = 1) +
    scale_fill_manual(values = pal) +
    labs(title = "Distribution FSCA", y = "FSCA") +
    base_theme
  
  print(p_dah / p_fsca)
  
  message("✓ Terminé.")
  
  
  
  
  
  
  
  
  
  
  
  
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
  
  
  
  
  
  
  
  ###############################################################################
  # 0. LIBRAIRIES ---------------------------------------------------------------
  ###############################################################################
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(dbscan)
    library(terra)
    library(sf)
    library(ggplot2)
  })
  
  ###############################################################################
  # 1. PARAMÈTRES UTILISATEUR ---------------------------------------------------
  ###############################################################################
  csv_in   <- "outputs/9. Analysis_Phenology/data_joint/dataset_pheno_2022_Cayolle.csv.gz"
  
  csv_out  <- "outputs/9. Analysis_Phenology/data_joint/zones_topo_habitat_chargement.csv.gz"
  ras_out  <- "outputs/9. Analysis_Phenology/data_joint/zonesTraitement.tif"
  shp_out  <- "outputs/9. Analysis_Phenology/data_joint/zonesTraitement.shp"
  
  thr_high <- 300      # >300  → Fort
  thr_low  <-  50      # ≤50   → Faible
  
  eps_m  <- 1000       # DBSCAN rayon (m)
  minpts <- 20         # pixels min / noyau
  
  hab_exclus <- c(NA, "Formations minérales")
  
  ###############################################################################
  # 2. LECTURE & PRÉ-TRAITEMENT -------------------------------------------------
  ###############################################################################
  dt <- fread(file = csv_in, showProgress = FALSE)
  
  # harmonisation coordonnées → x y
  for (p in list(c("x","y"), c("x.x","y.x"), c("x.y","y.y"),
                 c("X","Y"), c("lon","lat"), c("longitude","latitude")))
    if (all(p %in% names(dt))) { setnames(dt, p, c("x","y")); break }
  
  dt_first <- dt[DOY == min(DOY)]         # variables statiques
  
  ###############################################################################
  # 3. HABITATS MAJEURS ---------------------------------------------------------
  ###############################################################################
  top3 <- dt_first[ !habitat_label %in% hab_exclus,
                    .N, by = habitat_label ][order(-N)][1:3, habitat_label]
  
  ###############################################################################
  # 4. FONCTION ZONES TOPO HOMOGÈNES -------------------------------------------
  ###############################################################################
  create_zones <- function(sub_dt) {
    med_dah  <- median(sub_dt$dah , na.rm = TRUE)
    med_fsca <- median(sub_dt$fsca, na.rm = TRUE)
    win_dah  <- 0.10 * abs(med_dah)
    win_fsca <- 0.10 * abs(med_fsca)
    
    cand <- sub_dt[dah  %between% c(med_dah - win_dah,  med_dah + win_dah) &
                     fsca %between% c(med_fsca - win_fsca, med_fsca + win_fsca)]
    if (nrow(cand) < minpts) return(NULL)
    
    coords <- as.matrix(cand[, .(x, y)])
    cl <- dbscan(coords, eps = eps_m, minPts = minpts)
    cand$cluster <- cl$cluster
    cand[cluster > 0]
  }
  
  ###############################################################################
  # 5. TRAITEMENT PAR HABITAT ---------------------------------------------------
  ###############################################################################
  zones_all <- list()
  
  for (hab in top3) {
    
    message("• Habitat : ", hab)
    
    z <- create_zones(dt_first[habitat_label == hab])
    if (is.null(z)) { warning("  aucune zone."); next }
    
    # classes Faible / Fort
    z[ , traitement := fcase(
      Charge >  thr_high, "Fort",
      Charge <= thr_low , "Faible",
      default = NA_character_)]
    z <- z[!is.na(traitement)]
    
    # clusters contenant les deux classes
    tab  <- z[ , .N, by = .(cluster, traitement)]
    wide <- dcast(tab, cluster ~ traitement, value.var = "N", fill = 0)
    if (!all(c("Fort","Faible") %in% names(wide))) next
    keep <- wide[Fort >= minpts & Faible >= minpts, cluster]
    
    z <- z[cluster %in% keep]
    if (!nrow(z)) next
    
    # renumérote 3 plus gros clusters
    top_cl <- z[ , .N, by = cluster][order(-N)][1:3,
                                                .(cluster, zone_id = seq_len(.N))]
    z <- merge(z, top_cl, by = "cluster", all.y = TRUE)
    
    z[ , habitat_label := hab]
    zones_all[[hab]] <- z[, .(cell, x, y, zone_id,
                              habitat_label, traitement,
                              dah, fsca)]
  }
  
  ###############################################################################
  # 6. EXPORT CSV ---------------------------------------------------------------
  ###############################################################################
  if (!length(zones_all)) stop("❌  aucune zone valide.")
  final <- rbindlist(zones_all)
  write_csv(final, csv_out)
  message("✓ CSV exporté : ", csv_out)
  
  ###############################################################################
  # 7. RASTER & SHAPEFILE -------------------------------------------------------
  ###############################################################################
  template <- rast("raster/Phenologie/NDVI_season_Cayolle_2022.tif", subds = 1)
  zone_r <- rast(template); values(zone_r) <- NA_integer_
  
  valid <- !is.na(final$cell)
  code_tab <- c(Faible = 1L, Fort = 2L)
  final[ , code := code_tab[traitement]]
  
  zone_r[ final$cell[valid] ] <- final$code[valid]
  writeRaster(zone_r, ras_out, datatype = "INT1U", overwrite = TRUE)
  message("✓ Raster exporté : ", ras_out, "  (1 Faible, 2 Fort)")
  
  xy  <- xyFromCell(zone_r, final$cell[valid])
  sf_pts <- st_as_sf(data.frame(xy,
                                code = final$code[valid]),
                     coords = c("x","y"),
                     crs = crs(zone_r))
  st_write(sf_pts, shp_out, delete_dsn = TRUE, quiet = TRUE)
  message("✓ Shapefile exporté : ", shp_out, "  (champ 'code'=1/2)")
  
  ###############################################################################
  # 8. BOXPLOTS & BARPLOT -------------------------------------------------------
  ###############################################################################
  pix_area_ha <- prod(res(template)) / 1e4      # surface pixel en ha
  
  plot_dt <- final[, .(traitement,
                       DAH  = dah,
                       FSCA = fsca)]
  
  ###############################################################################
  # 9.  BOXPLOTS EMPILÉS  DAH & FSCA  (NA retirés) -------------------------------
  ###############################################################################
  library(ggplot2)
  library(patchwork)   # pour empiler élégamment 2 ggplot
  
  ## 1. Nettoyage : on enlève toute ligne NA sur DAH ou FSCA ----------------------
  plot_dt <- final[!is.na(dah) & !is.na(fsca)]
  
  ## 2. Palette & thème commun ---------------------------------------------------
  pal <- c(Faible = "khaki3", Fort = "steelblue")
  
  base_theme <- theme_bw(base_size = 11) +
    theme(
      axis.title.x  = element_blank(),
      legend.position = "none",
      panel.grid.major.y = element_line(colour = "grey85", linewidth = .3),
      panel.grid.minor.y = element_blank()
    )
  
  ## 3. Boxplot DAH --------------------------------------------------------------
  p_dah <- ggplot(plot_dt, aes(traitement, dah, fill = traitement)) +
    geom_boxplot(width = .65, outlier.shape = 21, outlier.size = 1) +
    scale_fill_manual(values = pal) +
    labs(title = "Distribution DAH", y = "DAH") +
    base_theme
  
  ## 4. Boxplot FSCA -------------------------------------------------------------
  p_fsca <- ggplot(plot_dt, aes(traitement, fsca, fill = traitement)) +
    geom_boxplot(width = .65, outlier.shape = 21, outlier.size = 1) +
    scale_fill_manual(values = pal) +
    labs(title = "Distribution FSCA", y = "FSCA") +
    base_theme
  
  ## 5. Empilement vertical ------------------------------------------------------
  p_dah / p_fsca         # l’opérateur "/" de {patchwork}
  
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



## ET LE PLOT FINAL 
  
  
  # ─────────────────────────────────────────────────────────────────────────────
  # PARAMÈTRES I/O
  # ─────────────────────────────────────────────────────────────────────────────
  csv_zones <- "outputs/9. Analysis_Phenology/data_joint/zones_topo_habitat_chargement.csv.gz"
  csv_full  <- "outputs/9. Analysis_Phenology/data_joint/dataset_pheno_2022_Cayolle.csv.gz"
  
  # ─────────────────────────────────────────────────────────────────────────────
  # LIBRAIRIES
  # ─────────────────────────────────────────────────────────────────────────────
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
  })
  
  # ─────────────────────────────────────────────────────────────────────────────
  # 1. LECTURE DES FICHIERS
  # ─────────────────────────────────────────────────────────────────────────────
  library(readr)
  zones <- read_csv(csv_zones, show_col_types = FALSE)
  full  <- read_csv(csv_full,  show_col_types = FALSE,
                    col_select = c(cell, DOY, NDVI))
  
  
  
  
  
  # ─────────────────────────────────────────────────────────────────────────────
  # 0.  CHARGER data.table  (si ce n’est pas déjà fait plus haut)
  # ─────────────────────────────────────────────────────────────────────────────
  library(data.table)
  
  # convertit en data.table si besoin
  setDT(full)
  setDT(zones)
  
  # ─────────────────────────────────────────────────────────────────────────────
  # 2.  SÉLECTION DES PIXELS D’INTÉRÊT
  # ─────────────────────────────────────────────────────────────────────────────
  sel <- merge(
    full,
    zones[, .(cell, traitement)],   # <- fonctionne car data.table chargé
    by = "cell"
  )
  
  
  
  # ─────────────────────────────────────────────────────────────────────────────
  # 3. RÉSUMÉ JOUR × TRAITEMENT
  # ─────────────────────────────────────────────────────────────────────────────
  plot_df <- sel %>%
    group_by(traitement, DOY) %>%
    summarise(med   = median(NDVI, na.rm = TRUE),
              low95 = quantile(NDVI, .2, na.rm = TRUE),
              hi95  = quantile(NDVI, .8, na.rm = TRUE),
              .groups = "drop")
  
  ###############################################################################
  # 1. COULEURS & ÉTIQUETTES ----------------------------------------------------
  ###############################################################################
  cols <- c(Fort   = "steelblue",
            Faible = "khaki3")
  
  labs_trt <- c(Fort   = "high grazing intensity",
                Faible = "low  grazing intensity")
  
  ###############################################################################
  # 2. LIMITE Y AUTOMATIQUE (±10 % marge) ---------------------------------------
  ###############################################################################
  y_min <- min(plot_df$low95, na.rm = TRUE)
  y_max <- max(plot_df$hi95 , na.rm = TRUE)
  y_pad <- .10 * (y_max - y_min)            # 10 % de marge
  
  ###############################################################################
  # 3. GRILLE : pas majeur = 0.25  |  pas mineur = 0.05 -------------------------
  ###############################################################################
  major_breaks <- seq(floor(y_min), ceiling(y_max), 0.25)
  minor_breaks <- seq(floor(y_min), ceiling(y_max), 0.05)
  
  ###############################################################################
  # 4. PLOT GGLOT ---------------------------------------------------------------
  ###############################################################################
  ggplot(plot_df, aes(DOY, med, colour = traitement, fill = traitement)) +
    geom_ribbon(aes(ymin = low95, ymax = hi95),
                alpha = .25, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = cols, labels = labs_trt, name = NULL) +
    scale_fill_manual(values  = cols, labels = labs_trt, name = NULL) +
    scale_x_continuous(breaks       = seq(60, 330, 30),
                       minor_breaks = seq(60, 330, 10)) +
    scale_y_continuous(limits       = c(y_min - y_pad, y_max + y_pad),
                       breaks       = major_breaks,
                       minor_breaks = minor_breaks,
                       expand       = c(0, 0)) +
    labs(x = "Julian date",
         y = "Plant Phenology Index (PPI)") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_line(size = .2, linetype = "dotted"),
      panel.grid.major = element_line(size = .3),
      legend.position.inside = c(.75, .90),   # nouvelle syntaxe ggplot 3.5
      legend.background = element_blank()
    )



































































































































# PART 3 : Plot 
  
  
  # ENTRE
  
  #Création du dossier de sortie des indicateur pour la visualistaion
  clim_case <- file.path(output_dir, "8. Analysis_Climate")
  data_case <- file.path(clim_case, "Data_Use_Fsca_Alti")
  clim_data_rds_file <- file.path(data_case, paste0("Use_Fsca_Alti_",alpage,".rds"))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Fonction pour le df (en RDS)
  
  
  
  library(raster)
  
  # Chemin du raster FSCA calculé pour l'alpage (par exemple "Cayolle")
  fsca_raster_path <- "C:/Users/masso/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/outputs/8. Analysis_Climate/SMOD/Fsca_Cayolle.tif"
  
  # Chemin du raster DEM couvrant l'alpage
  dem_raster_path <- "C:/Users/masso/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/raster/Alti/Cayolle_MNT.tif"
  
  # Charger les rasters depuis les fichiers
  r_fsca <- raster(fsca_raster_path)
  r_dem  <- raster(dem_raster_path)
  
  # Harmoniser les CRS si nécessaire
  if (!compareCRS(r_fsca, r_dem)) {
    message("Les CRS diffèrent, reprojection de r_fsca...")
    r_fsca <- projectRaster(r_fsca, crs = crs(r_dem))
  }
  
  # Agréger le DEM pour passer de 5 m à 20 m (facteur 4)
  r_dem_agg <- aggregate(r_dem, fact = 4, fun = mean)
  
  # Déterminer l'emprise commune entre les deux rasters
  common_extent <- intersect(extent(r_dem_agg), extent(r_fsca))
  if (is.null(common_extent)) {
    stop("Les extents des rasters ne se recoupent pas. Vérifiez vos données.")
  }
  
  # Recadrer les deux rasters à l'emprise commune
  r_fsca_crop <- crop(r_fsca, common_extent)
  r_dem_crop  <- crop(r_dem_agg, common_extent)
  
  # Resampler le DEM agrégé pour qu'il s'aligne exactement sur la grille du FSCA
  r_dem_resampled <- resample(r_dem_crop, r_fsca_crop, method = "bilinear")
  
  # Maintenant, chaque pixel du FSCA a une altitude associée (dans r_dem_resampled)
  merged_raster <- stack(r_dem_resampled, r_fsca_crop)
  plot(merged_raster)
  
  
  
  
  
  library(raster)
  library(ggplot2)
  library(dplyr)
  
  # Nom de l'alpage (à adapter si besoin)
  alpage <- "Cayolle"
  
  
  
  
  # ENTREE
  
  visu_case <- file.path(output_dir, "5. Indicateurs_visualisation")
  
  polygon_case <- file.path(visu_case, "Utilisation_par_quinzaine")
  
  polygon_use_shp = file.path(polygon_case, paste0("Use_polygon_",YEAR,"_",alpage,".shp"))
  
  
  # Chemin du raster FSCA calculé pour l'alpage
  fsca_raster_path <- "C:/Users/masso/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/outputs/8. Analysis_Climate/SMOD/Fsca_123_342_Cayolle.tif"
  
  # Chemin du raster DEM couvrant l'alpage
  dem_raster_path <- "C:/Users/masso/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/raster/Alti/Cayolle_MNT.tif"
  
  # Charger les rasters depuis les fichiers
  r_fsca <- raster(fsca_raster_path)
  r_dem  <- raster(dem_raster_path)
  
  # Harmoniser le CRS si nécessaire
  if (!compareCRS(r_fsca, r_dem)) {
    message("Les CRS diffèrent, reprojection de r_fsca...")
    r_fsca <- projectRaster(r_fsca, crs = crs(r_dem))
  }
  
  # Agréger le DEM pour passer de 5 m à 20 m (facteur 4)
  r_dem_agg <- aggregate(r_dem, fact = 4, fun = mean)
  
  # Définir l'emprise commune entre les deux rasters
  common_extent <- intersect(extent(r_dem_agg), extent(r_fsca))
  if (is.null(common_extent)) {
    stop("Les extents des rasters ne se recoupent pas. Vérifiez vos données.")
  }
  
  # Recadrer les deux rasters sur l'emprise commune
  r_dem_crop  <- crop(r_dem_agg, common_extent)
  r_fsca_crop <- crop(r_fsca, common_extent)
  
  # Resampler le DEM recadré pour qu'il s'aligne exactement sur la grille du FSCA
  r_dem_resampled <- resample(r_dem_crop, r_fsca_crop, method = "bilinear")
  
  # Créer une pile alignée : couche 1 = altitude, couche 2 = FSCA
  merged_raster <- stack(r_dem_resampled, r_fsca_crop)
  
  # Convertir le raster empilé en data frame (chaque pixel devient une ligne)
  df_merged <- as.data.frame(rasterToPoints(merged_raster))
  # Renommer les colonnes pour faciliter l'utilisation
  names(df_merged)[3:4] <- c("altitude", "FSCA")
  
  # Retirer les éventuels NA
  df_merged <- df_merged %>% filter(!is.na(altitude), !is.na(FSCA))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  library(sf)
  library(dplyr)
  
  add_presence_info <- function(df_merged, polygon_use_shp, crs_raster) {
    # Lire le shapefile qui contient les polygones de présence par quinzaine
    polygon_sf <- st_read(polygon_use_shp, quiet = TRUE)
    
    # S'assurer que le shapefile est dans le même CRS que le raster
    polygon_sf <- st_transform(polygon_sf, crs = crs_raster)
    
    # Convertir le data frame des pixels en objet sf (points)
    # On suppose que df_merged contient les colonnes "x" et "y" (coordonnées)
    df_sf <- st_as_sf(df_merged, coords = c("x", "y"), crs = crs_raster)
    
    # Récupérer les périodes (les niveaux de "month_period")
    periods <- unique(polygon_sf$month_period)
    
    # Pour chaque période, calculer une colonne binaire indiquant la présence du troupeau
    for (p in periods) {
      # Sélectionner le polygone pour la période p
      poly_p <- polygon_sf %>% filter(month_period == p)
      # Pour chaque pixel, déterminer s'il intersecte le polygone de la période p
      inter <- st_intersects(df_sf, poly_p)
      # Créer une nouvelle colonne "presence_<p>" (on remplace espaces et tirets par des underscores)
      colname <- paste0("presence_", gsub("[ -]", "_", tolower(p)))
      df_sf[[colname]] <- sapply(inter, function(x) ifelse(length(x) > 0, 1, 0))
    }
    
    # Optionnel : Conserver les coordonnées en colonnes, en retirant la géométrie
    df_result <- cbind(as.data.frame(st_coordinates(df_sf)), st_drop_geometry(df_sf))
    return(df_result)
  }
  
  # Exemple d'utilisation :
  # df_merged est obtenu à partir de vos rasters (voir votre code de base)
  # polygon_use_shp est le chemin vers le shapefile généré
  # crs_raster peut être obtenu par, par exemple, crs(r_fsca)
  df_merged_with_presence <- add_presence_info(df_merged, polygon_use_shp, crs(r_fsca))
  
  # Vous obtenez ainsi un data frame df_merged_with_presence contenant les colonnes "altitude", "FSCA"
  # et pour chaque quinzaine (ex: "presence_avant_15_jun", "presence_16_-_30_jun", etc.) une valeur 1 ou 0.
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Tracé : nuage de points et courbe lissée
  p <- ggplot(df_merged, aes(x = altitude, y = FSCA)) +
    geom_point(alpha = 0.2, size = 0.5, color = "blue") +
    geom_smooth(method = "loess", color = "red", se = FALSE) +
    scale_y_reverse(limits = c(1, 0)) +
    labs(x = "Altitude (m)", y = "FSCA",
         title = paste("FSCA vs Altitude (lissée) pour", alpage)) +
    theme_minimal()
  
  print(p)
  
  
  
  
  
  
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  
  # -------------------------------
  # Supposons que df_merged existe et contient :
  #   - altitude : valeurs numériques (ex. de 1600 à 3000)
  #   - FSCA     : valeurs entre 0 et 1
  # -------------------------------
  
  # 1) Calcul des densités "inversées"
  dens_alt <- density(df_merged$altitude, na.rm = TRUE)
  dens_fsca <- density(df_merged$FSCA, na.rm = TRUE)
  
  # Data frame pour la densité de l'altitude (multipliée par -1 pour renverser la courbe)
  df_dens_alt <- data.frame(
    altitude = dens_alt$x,
    dens = -dens_alt$y  # négatif => la bosse pointe vers le bas
  )
  
  # Data frame pour la densité du FSCA (multipliée par -1 pour renverser la courbe)
  df_dens_fsca <- data.frame(
    FSCA = dens_fsca$x,
    dens = -dens_fsca$y  # négatif => la bosse pointe vers la gauche
  )
  
  # 2) Nuage de points principal : Altitude vs FSCA
  p_main <- ggplot(df_merged, aes(x = altitude, y = FSCA)) +
    geom_point(alpha = 0.2, size = 0.5, color = "blue") +
    geom_smooth(method = "loess", color = "red", se = FALSE) +
    scale_x_continuous(limits = c(min(df_merged$altitude), max(df_merged$altitude)), expand = c(0, 0)) +
    scale_y_reverse(limits = c(1, 0), expand = c(0, 0)) +  # FSCA : 1 en haut, 0 en bas
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(x = "Altitude (m)", y = "FSCA")
  
  # 3) Densité de l'altitude en bas (p_xdens) : tracée vers le bas (densité négative)
  p_xdens <- ggplot(df_dens_alt, aes(x = altitude, y = dens)) +
    geom_line(color = "black", size = 1) +
    scale_x_continuous(limits = c(min(df_merged$altitude), max(df_merged$altitude)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(df_dens_alt$dens), 0), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  # 4) Densité du FSCA à gauche (p_ydens) : tracée vers la gauche (densité négative)
  p_ydens <- ggplot(df_dens_fsca, aes(x = FSCA, y = dens)) +
    geom_line(color = "black", size = 1) +
    scale_x_reverse(limits = c(1, 0), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(df_dens_fsca$dens), 0), expand = c(0, 0)) +
    coord_flip() +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  # 5) Assemblage final avec patchwork
  # - On crée d'abord une rangée avec la densité FSCA (à gauche) et le nuage de points au centre
  # - Puis on ajoute la densité altitude en dessous
  left_and_main <- wrap_plots(list(p_ydens, p_main), ncol = 2, widths = c(0.2, 1))
  final_plot <- wrap_plots(list(left_and_main, p_xdens), ncol = 1, heights = c(1, 0.25))
  
  # Affichage du graphique final
  print(final_plot)
  
  
  
  
  
  
  
  
  
  
  library(ggplot2)
  library(ggExtra)
  
  # ------------------------------------------------------------------
  # 1) Données d'exemple (remplacez par vos données réelles)
  # ------------------------------------------------------------------
  set.seed(123)
  df_merged <- data.frame(
    altitude = runif(1000, 1600, 3000),
    FSCA     = runif(1000, 0, 1)  # valeurs entre 0 et 1
  )
  
  # ------------------------------------------------------------------
  # 2) Nuage de points + courbe lissée + FSCA inversé
  # ------------------------------------------------------------------
  p <- ggplot(df_merged, aes(x = altitude, y = FSCA)) +
    geom_point(alpha = 0.2, size = 1.5, color = "steelblue") +
    geom_smooth(method = "loess", se = FALSE, color = "red", size = 1.2) +
    # FSCA inversé : 1 en haut, 0 en bas
    scale_y_reverse(limits = c(1, 0), expand = c(0, 0)) +
    # Altitude sans marge
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    labs(x = "Altitude (m)", y = "FSCA")
  
  # ------------------------------------------------------------------
  # 3) Distributions marginales (densités) avec ggMarginal
  #    - On utilise trim=TRUE, from=0, to=1 pour l'axe FSCA
  #      afin que la densité ne dépasse pas le range [0..1].
  #    - size < 2 pour réduire l'espace des marges
  # ------------------------------------------------------------------
  p_marg <- ggMarginal(
    p,
    type       = "density",
    margins    = "both",            # densités en haut et à droite
    fill       = "grey",
    color      = "lightgrey",
    alpha      = 0.5,
    size       = 6,               # plus petit => marges plus étroites
    # Arguments passés à geom_density() pour la marge X (altitude)
    xparams    = list(
      adjust = 1.5,
      trim   = TRUE,
      from   = min(df_merged$altitude),
      to     = max(df_merged$altitude)
    ),
    # Arguments passés à geom_density() pour la marge Y (FSCA)
    yparams    = list(
      adjust = 1.5,
      trim   = TRUE,
      from   = 0,
      to     = 1
    ),
    # Retire axes et ticks sur les marges
    marginalTheme = theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
  )
  
  # ------------------------------------------------------------------
  # 4) Affichage final
  # ------------------------------------------------------------------
  print(p_marg)
  
  
  
  
  
  # Foncyion pour créer les MNT
  
  library(raster)
  
  # 1. Lister les fichiers DEM à fusionner
  dem_files <- c("Cayolle_1.tif", 
                 "Cayolle_2.tif", 
                 "Cayolle_3.tif", 
                 "Cayolle_4.tif")
  
  # 2. Charger chaque fichier en tant qu'objet raster
  dem_list <- lapply(dem_files, raster)
  
  # 3. Fusionner tous les rasters en un seul
  #    do.call(merge, dem_list) applique la fonction 'merge' à tous les éléments de la liste
  mnt_merged <- do.call(merge, dem_list)
  
  # 4. Enregistrer le MNT fusionné sous forme d'un nouveau fichier TIF
  writeRaster(mnt_merged, filename = "Cayolle_MNT.tif", format = "GTiff", overwrite = TRUE)
  
  cat("Fusion terminée. Le fichier 'Cayolle_MNT.tif' a été créé.\n")
  
  