##----------------------------------------------------------------------------##
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



source("config.R")

## PARAMETRE : 

alpage = "Cayolle"
YEAR = "2023"



#### 0. Transformation des données ####
#-------------------------------------#
if(TRUE){
  # Descriprion : préparation du smod de l'année, transformation de DOY en année
  # hydrique (1septembre) en année courante (1 janvier)
  
  ## LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_legacy.R"))
  
  clim_case <- file.path(output_dir, "8. Analysis_Climate")
  SMOD_case <- file.path(clim_case, "SMOD")
  SMOD_tif_file_2 <- file.path(SMOD_case, paste0("SMOD_", alpage, "_", YEAR, ".tif"))
  
  case_UP_file <- file.path(raster_dir, "UP")
  UP_file      <- file.path(case_UP_file, paste0("UP_", alpage, "_bis.shp"))
  
  stack_path <- build_smod_fsca_stack(
    alpage         = alpage,
    YEAR           = YEAR,
    smod_hydro_tif = SMOD_tif_file_2,
    out_dir        = SMOD_case,
    up_shape       = UP_file,
    smooth   = TRUE,
    kernel   = "gauss",   # "gauss" | "median" | "mean" | "circle"
    radius_m = 60,
    round_DOY = TRUE
  )
  


  
  
  
  
  
  
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
  
  
  ##PARAMETRE
  YEAR = 2023
  alpage = "Cayolle"
  
  
  # LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_legacy.R"))
  
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
  
  
  
  #Dossier de sortie des indicateur pour la visualistaion
  output_visu_case <- file.path(output_dir, "5. Indicateurs_visualisation")
  
  #Sous-dossier Indicateur traitée : Chargement
  output_chargement_case <- file.path(output_visu_case, "Taux_chargement")
  
  #Création du sous-sous-dossier alpage et années traitée : Chargement
  output_case_alpage <- file.path(output_chargement_case, paste0(YEAR,"_",alpage))
  
  
  
  ## Un .RDS contenant la date de première utilisation du pixel (seuil a un chargement de 30 brebis.jour.hectare)
  rds_out <- file.path(output_case_alpage,
                       paste0("first_day_use_", YEAR, "_", alpage, ".rds"))
  
  
  # Un .TIF du SMOD Médian sur 10 ans 
  clim_case <- file.path(output_dir, "8. Analysis_Climate")
  SMOD_case <- file.path(clim_case, "SMOD")
  SMOD_tif_file <- file.path(SMOD_case, paste0("Fsca_",alpage,".tif"))
  # Un .TIF du Smod de l'année :
  SMOD_FSCA_stack_tif <- file.path(SMOD_case,paste0("SMOD_FSCA_stack_", alpage, "_", YEAR, "_epsg2154.tif"))
  
  
  
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
  
  
  # A REPRENDRE PROPRE
  total_rds_prev <- file.path(output_dir, "4. Chargements_Calcules",
                              paste0(2022, "_", alpage),
                              paste0("total_", 2022, "_", alpage, ".rds"))
  
  first_use_2023 <- file.path(output_dir, "5. Indicateurs_visualisation",
                              "Taux_chargement",
                              paste0(2023, "_", alpage),
                              paste0("first_day_use_", 2023, "_", alpage, ".rds"))
  
  
  
  
  
  
  
  ## CODE
  
  
  ## FONCTION 1 : 
  # Extraction et agrégation de tout raster, rds, ... pour la création du dataset
  # pour chaque alpage
  
  ## CODE
  build_pheno_dataset(
    alpage             = alpage,
    YEAR               = YEAR,
    irg_tif            = IRG_tif_file,
    ndvi_tif           = NDVI_tif_file,
    habitat_tif        = carto_file,
    dah_tif            = DAH_tif_file,
    smod_fsca_stack_tif= SMOD_FSCA_stack_tif,   # ← le nouveau stack 2 bandes
    extra_tif          = EXTRA_tif,
    charge_prev_rds    = total_rds_prev,
    first_use_rds      = first_use_2023,
    lut_habitat_csv    = lut_csv,
    up_shape           = UP_file,
    output_rds         = out_rds
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
  # Cocaténation des jeu de données pour regrouper les alpages, dans notre cas :
  # agrégation de Viso, Cayolle et Sanguiniere donnent le dataset Alpe-Sud
  alpages <- c("Viso", "Cayolle", "Sanguiniere")
  merge_pheno_dataset(year = 2023, alpages = alpages, combo_name = "Alpe-Sud")







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
  YEAR = 2023
  
  
  # LIBRARY & FUNCTION
  source(file.path(functions_dir, "Functions_legacy.R"))
  library(dplyr)
  
  
  
  ## ENTREE
  out_dir <- file.path(output_dir, "10. NDVI Effect",
                       paste0("dataset_", alpage))
  
  dt_alpesud = file.path(out_dir, paste0("dataset_pheno_LONG_", YEAR, "_", alpage, ".rds"))
  
  ## SORTIE
  dt_filtered = file.path(out_dir, paste0("dataset_legacy_filtered_", YEAR, "_", alpage, ".rds"))
  dt_legacy = file.path(out_dir, paste0("dataset_legacy_", YEAR, "_", alpage, ".rds"))
  
  
  ## CODE 
  
  ## Filtre des pixel dont first_day_use < MAXD
  
  # VERSION 1 : Habitat filtré :
  if (FALSE){
    df_filtre <- readRDS(dt_alpesud) %>% 
      filter(
        is.na(first_day_use) |         # pixel jamais vraiment utilisé"~/STAGE_M2_PERSEE/R_studio/Script_ressource/Climate_synthesis/Inputs/PIX_h18v04b_EAS.csv"
          first_day_use >= MAXD          # pixel dont 1er jour d’usage ≥ date du pic
      ) %>% 
      filter(habitat_label %in% c("P. nivales", 
                                  "P. thermiques écorchées", 
                                  "Queyrellins", 
                                  "Nardaies denses du subalpin",
                                  "Megaphorbiaies et Aulnaies"))
    saveRDS(df_filtre, dt_filtered) %>% 
      cat("Save dataset :", dt_filtered, "\n")
  }
  
  # VERSION 2 : Habitat non filtré :
  if (TRUE){
    df_filtre <- readRDS(dt_alpesud) %>%
      dplyr::filter(
        is.na(first_day_use) | first_day_use >= MAXD
      )
    
    saveRDS(df_filtre, dt_legacy)
    cat("Save dataset :", dt_legacy, "\n")
  }
  
  ## TABLEAU
  
  ## A REPRENDRE 
  
  
  {
  
  
  library(dplyr)
  library(knitr)
  library(kableExtra)
  
  # ─────────────────────────── 1. Jeu réduit + comptages ──────────────────
  habitats_sel <- c("Megaphorbiaies et Aulnaies",
                    "Nardaies denses du subalpin",
                    "P. nivales",
                    "Queyrellins",
                    "P. thermiques écorchées"
                    )
  
  df_init <- readRDS(dt_alpesud) %>% 
    filter(habitat_label %in% habitats_sel) %>% 
    distinct(cell, .keep_all = TRUE)
  
  df_filt <- df_init %>% 
    filter(is.na(first_day_use) | first_day_use >= MAXD)
  
  resume <- df_init %>% 
    count(habitat_label, alpage, name = "Pixels (tot.)") %>% 
    left_join(
      df_filt %>% count(habitat_label, alpage, name = "Pixels retenus"),
      by = c("habitat_label","alpage")
    ) %>% 
    mutate(`Pixels retenus` = replace_na(`Pixels retenus`, 0),
           `Pixels retirés` = `Pixels (tot.)` - `Pixels retenus`) %>% 
    left_join(
      df_filt %>% 
        filter(!is.na(Charge) & Charge > 200) %>% 
        count(habitat_label, alpage, name = "Charge > 200"),
      by = c("habitat_label","alpage")
    ) %>% 
    mutate(`Charge > 200` = replace_na(`Charge > 200`, 0))
  
  totaux <- resume %>% 
    group_by(habitat_label) %>% 
    summarise(`Total retenus` = sum(`Pixels retenus`),
              `Total > 200`  = sum(`Charge > 200`), .groups = "drop")
  
  # ─────────────────────────── 2. Tableau long prêt à styler ───────────────
  tab <- resume %>% 
    left_join(totaux, by = "habitat_label") %>% 
    arrange(match(habitat_label, habitats_sel),
            match(alpage, c("Cayolle","Viso","Sanguinière"))) %>% 
    group_by(habitat_label) %>% 
    mutate(
      Habitat        = if_else(row_number()==1, habitat_label, ""),
      `Total retenus`= if_else(row_number()==1, as.character(`Total retenus`), ""),
      `Total > 200`  = if_else(row_number()==1, as.character(`Total > 200`), "")
    ) %>% 
    ungroup() %>% 
    select(Habitat,
           Alpage = alpage,
           `Pixels (tot.)`,
           `Pixels retenus`,
           `Pixels retirés`,
           `Charge > 200`,
           `Total retenus`,
           `Total > 200`) %>% 
    # index pour styling
    mutate(row_id = row_number(),
           id_hab = cumsum(Habitat != ""))
  
  # lignes grisées : Cayolle, Sanguinière, …
  zebra_rows <- tab %>% 
    filter(Habitat == "") %>% 
    group_by(id_hab) %>% 
    filter(row_number() %% 2 == 1) %>%         # 1re, 3e, …
    pull(row_id)
  
  # ── sur‑styler Habitat + Totaux sur ces lignes (fond blanc, pas de bordure)
  tab <- tab %>% 
    mutate(
      Habitat = if_else(row_id %in% zebra_rows,
                        cell_spec(Habitat, extra_css = "border-bottom:0;background:#ffffff;"),
                        Habitat),
      `Total retenus` = if_else(row_id %in% zebra_rows & `Total retenus`!="",
                                cell_spec(`Total retenus`, bold = TRUE, align = "center",
                                          font_size = 14,
                                          extra_css = "border-bottom:0;background:#ffffff;"),
                                `Total retenus`),
      `Total > 200`   = if_else(row_id %in% zebra_rows & `Total > 200`!="",
                                cell_spec(`Total > 200`, bold = TRUE, align = "center",
                                          font_size = 14,
                                          extra_css = "border-bottom:0;background:#ffffff;"),
                                `Total > 200`)
    )
  
  # ─────────────────────────── 3. Rendu HTML final ─────────────────────────
  kable(tab %>% select(-row_id, -id_hab),
        "html", escape = FALSE,
        align = c("l","l","r","r","r","r","r","r"),
        col.names = c("Habitat","Alpage",
                      "Pixels&nbsp;(tot.)","Pixels retenus","Pixels retirés",
                      "Charge&nbsp;&gt; 200",
                      "Total retenus","Total&nbsp;&gt; 200")) %>% 
    kable_styling(full_width = FALSE,
                  bootstrap_options = c("hover","condensed")) %>% 
    
    # fond gris + bordure sous colonnes 2‑6 uniquement
    row_spec(zebra_rows, background = "#f4f4f4",
             extra_css = "td:nth-child(2),td:nth-child(3),td:nth-child(4),td:nth-child(5),td:nth-child(6){border-bottom:1px solid #ddd;}") %>% 
    
    # fusion verticale cellules identiques (Habitat + Totaux)
    add_header_above(c(" " = 2, "Par alpage" = 4, "Totaux habitat" = 2)) %>% 
    collapse_rows(columns = c(1,7,8), valign = "top")
  
  
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
    YEAR = 2023
    
    # LIBRARY & FUNCTION
    source(file.path(functions_dir, "Functions_legacy.R"))
    library(ggplot2)
    library(dplyr)
    
    ## ENTREE
    out_dir <- file.path(output_dir, "10. NDVI Effect",
                         paste0("dataset_", alpage))
    
    # version habitat filtré
    dt_legacy_filtered = file.path(out_dir, paste0("dataset_legacy_filtered_", YEAR, "_", alpage, ".rds"))
    # version full
    dt_legacy = file.path(out_dir, paste0("dataset_legacy_", YEAR, "_", alpage, ".rds"))
    
    ## SORTIE
    dt_legacy_unique_file = file.path(out_dir, paste0("dataset_legacy_mod_", YEAR, "_", alpage, ".rds"))
    dt_legacy_habitat_unique_file = file.path(out_dir, paste0("dataset_legacy_mod_habitat", YEAR, "_", alpage, ".rds"))
    
    
  
    ## CODE
    
    # Lecture et préparation du dataset
    # Passage au pixel (unique, plus de DOY by pix, donc IRG et NDVI inutilisable
    # , mais non utilie pour cette partie), filtre des NA (7.8% des données : FSCA
    # peut etre a retravaillé) et utlisation de charge log
    
    # VERSION 1 : AVEC 5 Habitats :
    if(FALSE){
    dt_legacy_unique <- readRDS(dt_legacy_filtered) %>% 
      arrange(cell, DOY) %>%        # (optionnel) on range d’abord par pixel puis DOY
      group_by(cell) %>%            # 1 groupe = 1 pixel
      slice_head(n = 1) %>%         # garde la 1ʳᵉ ligne de chaque groupe
      ungroup() %>% 
      filter(!is.na(MAXV), !is.na(dah), !is.na(fsca), !is.na(Charge)) %>% 
      mutate(
        Charge_log = log1p(Charge)          # log(Charge + 1)
      )
    
    saveRDS(dt_legacy_unique, dt_legacy_unique_file) %>% 
    cat("Save dataset :", dt_legacy_unique_file, "\n")
    }
    
    
    # VERSION 2 : COMPLETE
    if(TRUE){
      dt_legacy_unique_habitat <- readRDS(dt_legacy) %>% 
        arrange(cell, DOY) %>%        # (optionnel) on range d’abord par pixel puis DOY
        group_by(cell) %>%            # 1 groupe = 1 pixel
        slice_head(n = 1) %>%         # garde la 1ʳᵉ ligne de chaque groupe
        ungroup() %>% 
        filter(!is.na(MAXV), !is.na(dah), !is.na(SMOD_2023), !is.na(Charge)) %>% 
        mutate(
          Charge_log = log1p(Charge)          # log(Charge + 1)
        )
      
      saveRDS(dt_legacy_unique_habitat, dt_legacy_habitat_unique_file ) %>% 
        cat("Save dataset :", dt_legacy_habitat_unique_file, "\n")
    }
    
    
    
    
    
    
    
    
    
    ## Histogramme de distribution des variables 
    
    # Chargement : 
    ggplot(dt_legacy_unique, aes(x = Charge_log)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~ habitat_label, scales = "free_y") +
      labs(title = "Log du chargement par habitat")
    
    # DAH :
    ggplot(dt_legacy_unique, aes(x = dah)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~ habitat_label, scales = "free_y") +
      labs(title = "DAH par habitat")
    
    # FSCA : 
    ggplot(dt_legacy_unique, aes(x = fsca)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~ habitat_label, scales = "free_y") +
      scale_x_reverse(limits = c(1, 0)) +     # ← 1 à gauche, 0 à droite
      labs(title = "FSCA par habitat",
           x = "fSCA", y = "Nombre de pixels")
    
    # MAXV :
    ggplot(dt_legacy_unique, aes(x = MAXV)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~ habitat_label, scales = "free_y") +
      labs(title = "MAXV par habitat")
    
    
    # MAXV
    ggplot(dt_legacy_unique, aes(x = alpage, y = MAXV)) +
      geom_boxplot() +
      facet_wrap(~ habitat_label, scales = "free_y") +
      labs(title = "MAXV par habitat")
    
    
    
    # ONSET10 :
    ggplot(dt_legacy_unique, aes(x = ONSET10)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~ habitat_label, scales = "free_y") +
      labs(title = "ONSET10 par habitat")
    
    
    # SMOD de l'année
    # par habitat
    ggplot(dt_legacy_unique_habitat , aes(x = SMOD_2023)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~ habitat_label, scales = "free_y") +
      labs(title = "SMOD par habitat")
    
    # global
    ggplot(dt_legacy_unique_habitat , aes(x = SMOD_2023)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      labs(title = "SMOD")
    
    # AUC (intégrale du PPI sur la période de pousse)
    # par habitat
    ggplot(dt_legacy_unique_habitat , aes(x = AUCg_ONSET10_MAXD)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      facet_wrap(~ habitat_label, scales = "free_y") +
      labs(title = "AUC (intégrale du PPI période de pousse) par habitat")
    
    # global
    ggplot(dt_legacy_unique_habitat , aes(x = AUCg_ONSET10_MAXD)) +
      geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
      labs(title = "AUC (intégrale du PPI période de pousse)")
    
    
    
    
   
    ## Relation variable réponse (MAXV) vs variable explicative (FSCA, DAH, Charge)
    
    # MAXV vs Charge
    # non linéaire
    ggplot(dt_legacy_unique, aes(x = Charge_log, y = MAXV)) +
      geom_point(alpha = 0.1) +
      geom_smooth(method = "gam", formula = y ~ s(x), colour = "red") +
      facet_wrap(~ habitat_label, scales = "free") +
      theme_minimal()
    
    # linéaire 
    ggplot(dt_legacy_unique, aes(x = Charge_log, y = MAXV)) +
      geom_point(alpha = 0.1) +
      geom_smooth(method = "lm", formula = y ~ x, colour = "red", se = TRUE) +
      facet_wrap(~ habitat_label, scales = "free") +
      labs(x = "log1p(Charge)", y = "MAXV",
           title = "Relation linéaire brute : MAXV ~ log1p(Charge) par habitat") +
      theme_minimal()
    
    
    # MAXV vs FSCA
    ggplot(dt_legacy_unique, aes(x = fsca, y = MAXV)) +
      geom_point(alpha = 0.1) +
      geom_smooth(method = "gam", formula = y ~ s(x), colour = "red") +
      facet_wrap(~ habitat_label, scales = "free_y") +
      scale_x_reverse(limits = c(1, 0)) +     # ← 1 à gauche, 0 à droite
      labs(title = "FSCA par habitat",
           x = "fSCA", y = "Nombre de pixels")+
      theme_minimal()
    
    # MAXV vs DAH 
    ggplot(dt_legacy_unique, aes(x = dah, y = MAXV)) +
      geom_point(alpha = 0.1) +
      geom_smooth(method = "gam", formula = y ~ s(x), colour = "red") +
      facet_wrap(~ habitat_label, scales = "free") +
      theme_minimal()
    
    dt_plot <- dt_legacy_unique_habitat %>% 
    filter( ! habitat_label %in% "Formations minérales")
    
    
    # ONSET10 vs SMOD_2023
    # par habitat
    ggplot(dt_legacy_unique_habitat, aes(x = ONSET10, y = SMOD_2023)) +
      geom_point(alpha = 0.1) +
      geom_smooth(method = "gam", formula = y ~ s(x), colour = "red") +
      facet_wrap(~ habitat_label, scales = "free") +
      theme_minimal()
    
    # global
    ggplot(dt_plot , aes(x = SMOD_2023, y = ONSET10)) +
      geom_point(alpha = 0.1) +
      geom_smooth(method = "gam", formula = y ~ s(x), colour = "red") +
      theme_minimal()
    
    
    # AUC vs SMOD
    # global
    library(ggplot2)
    library(mgcv)  # pour le GAM
    
    # Ajustements (en enlevant les NA)
    lm_fit  <- lm(AUCg_ONSET10_MAXD ~ SMOD_2023, data = dt_plot)
    gam_fit <- mgcv::gam(AUCg_ONSET10_MAXD ~ s(SMOD_2023), data = dt_plot)
    
    r2_lm  <- summary(lm_fit)$r.squared            # R² du modèle linéaire
    r2_gam <- summary(gam_fit)$r.sq                # "R²" du GAM (pour famille gaussienne)
    
    ggplot(dt_plot, aes(x = SMOD_2023, y = AUCg_ONSET10_MAXD)) +
      geom_point(alpha = 0.1) +
      geom_smooth(method = "lm",  formula = y ~ x,    colour = "blue") +
      geom_smooth(method = "gam", formula = y ~ s(x), colour = "red") +
      annotate("text",
               x = -Inf, y = Inf,
               label = sprintf("LM : R² = %.3f\nGAM : R² = %.3f", r2_lm, r2_gam),
               hjust = -0.1, vjust = 1.1, size = 3.5) +
      theme_minimal()
    
    
    
    
    ## Matrice de corrélation 
    library(corrplot)
    
    vars <- dt_legacy_unique %>% 
      select(DAH = dah, FSCA = fsca, Charge_log) %>% 
      na.omit()                              
    
    mat_cor <- cor(vars, use = "complete.obs", method = "spearman")
    print(mat_cor)
    
    corrplot(mat_cor, method = "color", addCoef.col = "black",
             tl.col = "black", tl.cex = 0.9, number.cex = 0.8)
    
    # Variable OK !
    
    
    
  }
  
  
  
  #### 3.2 Modélisation ###
  ##Description :
  # Modélisation |variable réponse : MAXV (mais aussi : ONSET 10, GREENDUR, IRGMAX)
  #              |variable explicative : FSCA, DAH, CHARGE
  
  
  ## PARAMETRE
  alpage = "Alpe-Sud"
  YEAR = 2023
  
  ## ENTREE
  out_dir <- file.path(output_dir, "10. NDVI Effect",
                       paste0("dataset_", alpage))
  
  ## V1 : FILTRE SUR LES HABITATS
  dt_legacy_unique_file = file.path(out_dir, paste0("dataset_legacy_mod_", YEAR, "_", alpage, ".rds"))
  
  ## V2 : SMOD de l'anné, aucun filtre sur les habitats
  dt_legacy_habitat_unique_file = file.path(out_dir, paste0("dataset_legacy_mod_habitat", YEAR, "_", alpage, ".rds"))
  
  
  
  #### MODELISATION V1 ####
  #-----------------------#
  if (TRUE) {
    
    # LOAD DATASET V1
    df_legacy <-readRDS(dt_legacy_unique_file)
    
    ## Modélisation : Nardaies 
    if (TRUE){
      m0 <- glm(MAXV ~ 1,
                family = Gamma("log"),
                data = df_legacy%>%
                  filter(habitat_label == "Nardaies denses du subalpin"))
      
      
      
      m1 <- glm(MAXV ~ dah + fsca + Charge_log,
                family = Gamma("log"),
                data = df_legacy%>%
                  filter(habitat_label == "Nardaies denses du subalpin"))
      
      
      
      m2 <- glm(MAXV ~ dah * fsca * Charge_log,
                family = Gamma("log"),
                data = df_legacy%>%
                  filter(habitat_label == "Nardaies denses du subalpin"))
      
      
      summary(m1)
      
      
      
      
      
      
      
      
      res_df <- data.frame(
        fitted  = fitted(m1),                          # μ̂
        pearson = residuals(m1, type = "pearson"),     # résidus de Pearson
        dev     = residuals(m1, type = "deviance")     # résidus de déviance
      )
      
      ggplot(res_df, aes(x = fitted, y = pearson)) +
        geom_point(alpha = .3) +
        geom_smooth(method = "loess", colour = "red") +
        labs(title = "Résidus de Pearson vs prédictions",
             x = "Valeurs prédites (μ̂)", y = "Résidus de Pearson") +
        theme_minimal()
      
      
      qqnorm(res_df$dev, main = "QQ‑plot résidus de déviance")
      qqline(res_df$dev)
      
      
      ggplot(res_df, aes(sample = dev)) +
        stat_qq(alpha = .3) +
        stat_qq_line(colour = "red") +
        theme_minimal() +
        labs(title = "QQ‑plot (résidus déviance)")
      
      
      infl <- influence.measures(m1)
      which(apply(infl$is.inf, 1, any))              # indices des points influents
      
      
      plot(cooks.distance(m1), type="h", ylab="Cook's distance")
      abline(h = 4/length(res_df$fitted), col="red", lty = 2)
      
      
      
      
      mu  <- fitted(m1)
      var_pred <- (m1$family$variance(mu)) * summary(m1)$dispersion   # φ * μ²
      
      plot(mu, res_df$pearson^2, pch = 20, col = "grey50",
           xlab = "μ̂", ylab = "Résidu Pearson²", main = "Var ≈ φ μ² ?")
      lines(sort(mu), sort(var_pred), col = "red", lwd = 2)
      
      dispersion_phi <- summary(m1)$dispersion    # déjà calculé par glm
      phi_expected   <- 1                        # pour Gamma théorique
      
      cat("phi =", dispersion_phi, "\n")
      
      library(car)
      vif(m1)     
      
      
    }
    
    ## Modélisation : Mégaphorbaiie et aulnaies
    if (TRUE){
      m1_1 <- glm(MAXV ~ alpage*(dah + fsca + Charge_log),
                  family = Gamma("log"),
                  data = df_legacy%>%
                    filter(habitat_label == "Megaphorbiaies et Aulnaies"))
      summary(m1_1)
      
      
      
      
    }
    
    
    ## Modélisation : P. Thermique
    if (TRUE){
      m2_1 <- glm(MAXV ~ alpage*(dah + fsca + Charge_log),
                  family = Gamma("log"),
                  data = df_legacy%>%
                    filter(habitat_label == "P. thermiques écorchées"))
      summary(m2_1)
      
      
      library(emmeans)
      library(dplyr)
      library(purrr)
      
      # Helper robuste
      tidy_trend <- function(mod, var, by = "alpage") {
        e <- emtrends(mod, specs = as.formula(paste0("~ ", by)), var = var, infer = c(TRUE, TRUE))
        s <- summary(e) %>% as.data.frame()
        
        # Colonnes à détecter automatiquement
        slope_col <- grep("\\.trend$", names(s), value = TRUE)
        se_col    <- intersect(c("SE", "std.error"), names(s))[1]
        stat_col  <- intersect(c("t.ratio", "z.ratio"), names(s))[1]
        has_df    <- "df" %in% names(s)
        
        out <- s %>%
          rename(
            estimate = !!slope_col,
            SE       = !!se_col,
            stat     = !!stat_col
          ) %>%
          mutate(var = var) %>%
          relocate(var)
        
        if (has_df) {
          out <- out %>% select(var, !!by, estimate, SE, df, stat, p.value, lower.CL, upper.CL)
        } else {
          out <- out %>% select(var, !!by, estimate, SE, stat, p.value, lower.CL, upper.CL)
        }
        
        out
      }
      
      # Empiler les trois variables
      vars <- c("dah", "fsca", "Charge_log")
      res_trends <- map_dfr(vars, ~ tidy_trend(m2_1, .x))
      
      res_trends
      
    }
    
    ## Modélisation : P. Nivales
    if (TRUE){
      m3_1<- glmer(MAXV ~ fsca + dah + Charge_log + (1|alpage) ,
                   family = Gamma("log"),
                   data = df_legacy%>%
                     filter(habitat_label == "P. nivales"))
      summary(m3_1)
      
      
      
      
      m_fix <- glm(MAXV ~ alpage + dah + fsca + Charge_log,
                   family = Gamma("log"), data = df_legacy%>%
                     filter(habitat_label == "P. nivales"))
      m_mix <- glmer(MAXV ~ dah + fsca + Charge_log + (1|alpage),
                     family = Gamma("log"), data = df_legacy%>%
                       filter(habitat_label == "P. nivales"))
      
      m_fix_2 <- glm(MAXV ~ alpage * (dah + fsca + Charge_log),
                     family = Gamma("log"), data = df_legacy%>%
                       filter(habitat_label == "P. nivales"))
      AIC(m_fix, m_fix_2 , m_mix)          # souvent AIC ↓ en mixte
      anova(m_fix, m_fix_2, test = "Chisq")
      
      
      
      
      library(emmeans)
      get_trends <- function(var) {
        emtrends(m_fix_2, ~ alpage, var = var) %>%
          summary(infer = TRUE) %>%
          as.data.frame() %>%
          mutate(variable = var)
      }
      
      res_trends <- purrr::map_dfr(c("dah","fsca","Charge_log"), get_trends)
      res_trends
      
      pairs(emtrends(m_fix_2, ~ alpage, var = "dah"))
      pairs(emtrends(m_fix_2, ~ alpage, var = "fsca"))
      pairs(emtrends(m_fix_2, ~ alpage, var = "Charge_log"))
      
      
      
      df_pniv <- df_legacy %>%
        filter(habitat_label == "P. nivales") %>%
        mutate(across(c(dah, fsca, Charge_log), ~ . - mean(., na.rm = TRUE)))
      
      m_fix_2 <- glm(MAXV ~ alpage * (dah + fsca + Charge_log),
                     family = Gamma("log"), data = df_pniv)
      
      
      summary(m_fix_2)
      
      
      
      
      
      ## Modèle retenue et plot : 
      df_pniv <- df_legacy %>%
        filter(habitat_label == "P. nivales") %>%
        mutate(across(c(dah, fsca, Charge_log), ~ . - mean(., na.rm = TRUE)))
      
      m_fix_2 <- glm(
        MAXV ~ alpage * (dah + fsca + Charge_log),
        family = Gamma("log"),
        data = df_pniv
      )
      
      
      
      library(emmeans)
      library(dplyr)
      library(purrr)
      library(tidyr)
      
      library(emmeans)
      library(dplyr)
      
      # 1) On calcule d’abord les emmeans sur l’échelle du lien
      emm <- emmeans(m_fix_2, ~ alpage)
      
      # 2) Puis on demande la transformation sur l’échelle de la réponse dans summary()
      ints_raw <- summary(emm, type = "response", infer = TRUE) %>%
        as.data.frame()
      
      names(ints_raw)   # -> regarde: la colonne s'appelle sûrement "response", pas "emmean"
      
      # 3) Colonne à utiliser pour l’estimate (emmean ou response, selon le cas)
      est_col <- if ("emmean" %in% names(ints_raw)) "emmean" else "response"
      
      ints <- ints_raw %>%
        transmute(
          alpage,
          param     = "Intercept (μ)",
          scale     = "response",
          estimate  = .data[[est_col]],
          SE, df, lower.CL, upper.CL, t.ratio, p.value
        )
      ints
      
      
      
      
      
      
      
      library(purrr)
      
      tidy_trend <- function(var){
        emtrends(m_fix_2, ~ alpage, var = var) %>%
          summary(infer = TRUE) %>%
          as.data.frame() %>%
          transmute(
            alpage,
            param    = paste0("slope(", var, ")"),
            scale    = "link",  # sur l’échelle du lien (log)
            estimate = .data[[paste0(var, ".trend")]],
            SE, df, lower.CL, upper.CL, t.ratio, p.value
          )
      }
      
      res_trends <- map_dfr(c("dah","fsca","Charge_log"), tidy_trend)
      res_trends
      
      tab_final <- bind_rows(ints, res_trends) %>%
        arrange(alpage, desc(param == "Intercept (μ)"))
      
      tab_final
      
      
      
      library(ggeffects)
      
      # Effet de la charge (les autres prédicteurs sont fixés à 0 = leurs moyennes centrées)
      eff_charge <- ggeffect(m_fix_2, terms = c("Charge_log [all]", "alpage"))
      plot(eff_charge) + ggplot2::labs(
        x = "log1p(Charge)", y = "MAXV (prédit)",
        title = "Effet marginal de la charge par alpage (dah=fsca=0)"
      )
      
      # Idem pour dah et fsca
      eff_dah  <- ggeffect(m_fix_2, terms = c("dah [all]", "alpage"))
      eff_fsca <- ggeffect(m_fix_2, terms = c("fsca [all]", "alpage"))
      plot(eff_dah)
      plot(eff_fsca)
      
      
      
      
      
      
      
      library(tidyr)
      library(ggplot2)
      
      ## 1) On part du data.frame utilisé pour fitter m_fix_2
      levels_alp <- levels(df_pniv$alpage)
      
      ## 2) Grille propre : on fait varier Charge_log, on fixe dah = fsca = 0
      grid_charge <- tidyr::expand_grid(
        alpage     = levels_alp,
        Charge_log = seq(min(df_pniv$Charge_log), max(df_pniv$Charge_log), length.out = 200)
      ) %>%
        mutate(
          dah  = 0,
          fsca = 0,
          # très important : s’assurer que 'alpage' est bien un facteur avec les mêmes niveaux
          alpage = factor(alpage, levels = levels_alp)
        )
      
      ## 3) Prédictions
      ## (je préfère prédire sur l’échelle du lien puis revenir à la réponse)
      pred <- predict(m_fix_2, newdata = grid_charge, type = "link", se.fit = TRUE)
      
      grid_charge <- grid_charge %>%
        mutate(
          fit_link = pred$fit,
          se_link  = pred$se.fit,
          fit      = exp(fit_link),                         # retour à l’échelle réponse (Gamma log)
          lwr      = exp(fit_link - 1.96 * se_link),
          upr      = exp(fit_link + 1.96 * se_link)
        )
      
      ## 4) Plot
      ggplot(grid_charge, aes(x = Charge_log, y = fit, colour = alpage, fill = alpage)) +
        geom_line(size = 1) +
        geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .15, colour = NA) +
        labs(x = "log1p(Charge)", y = "MAXV prédit",
             title = "Courbes prédites par alpage (dah = fsca = 0)") +
        theme_minimal()
      
      
      
      
      
      
      
      
      
      
      
      
    }
    
    
    ## Modélisation : Queyrellins
    if (TRUE){
      m4_1 <- glm(MAXV ~ fsca + dah + Charge_log,
                  family = Gamma("log"),
                  data = df_legacy%>%
                    filter(habitat_label == "Queyrellins"))
      summary(m4_1)
      
      
      
      
      
      
    }
    
    
    
    
   
  
  
  
  
}

  
  #### MODELISATION V2 ####
  #-----------------------#
  
  
  
  if (TRUE){
    
    # 0. Préparation des données
    if (TRUE){
      #library:
      library(dplyr)
      
      
    data_mod <- readRDS(dt_legacy_habitat_unique_file) %>% 
      filter(!is.na(AUCg_ONSET10_MAXD),
            !is.na(Charge),      # pour le seuil
            !is.na(Charge_log),  # Gradient
            !is.na(alpage),
            !is.na(dah),
            !is.na(habitat_code),
            !is.na(habitat_label),
            !is.na(SMOD_2023),
            ! habitat_label %in% "Formations minérales") %>% 
      mutate(charge_bin = as.integer(Charge >= 1), # Binérisation de la charge avec un seuil a 1; inf.1 = 0 ; sup.1 = 1
             alpage = factor(alpage),
             AUCg = AUCg_ONSET10_MAXD,
             SMOD_2023_std = as.numeric(scale(SMOD_2023)),
             DAH_std       = as.numeric(scale(dah))) %>% 
      transmute(
        alpage,
        AUCg,
        # topo & neige (brut + standardisé)
        dah,        DAH_std,
        SMOD_2023,  SMOD_2023_std,
        # charge (hurdle)
        Charge, Charge_log, charge_bin, habitat_code, habitat_label
      )
      
    
    
    summary(data_mod)
    sapply(data_mod, function(v) sum(!is.finite(v)))  # NA ?
    table(data_mod$charge_bin) # CHarge global
    table(data_mod$alpage, data_mod$charge_bin) # Charge par alpage

    
    center_scale <- list(
      SMOD_2023 = c(center = mean(data_mod$SMOD_2023, na.rm = TRUE),
                    scale  = sd(data_mod$SMOD_2023,   na.rm = TRUE)),
      dah       = c(center = mean(data_mod$dah,       na.rm = TRUE),
                    scale  = sd(data_mod$dah,         na.rm = TRUE))
    )
    center_scale
    
    }
    
    # 1. Exploration ciblée
    if(TRUE){
      # 1.1 Distribution et valeur extrème :
      
      # Histogramme :
      # Variable réponse : AUCg
      ggplot(data_mod , aes(x = AUCg)) +
        geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
        labs(title = "AUC (intégrale du PPI période de pousse)")
      # Une légère queue ? ks.test pour la suite, connaitre la distribution
      
      # variable explicative : dah
      ggplot(data_mod , aes(x = dah)) +
        geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
        labs(title = "DAH")
      
      # variable explicative : smod
      ggplot(data_mod , aes(x = SMOD_2023_std)) +
        geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
        labs(title = "SMOD 2023 (standardisé))")
      
      # variable explicative d'interet : Charge_log
      ggplot(data_mod, aes(x = Charge_log)) +
        geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
        labs(
          title = "Charge (transformation log)",
          x = "Charge_log = log(Charge + 1)",
          y = "Effectifs"
        )
      
      
      # Variable explicative d'intéret : Charge log - pixel non chargé
      ggplot(dplyr::filter(data_mod, charge_bin == 1), aes(x = Charge_log)) +
        geom_histogram(bins = 40, colour = "white", fill = "steelblue") +
        labs(
          title = "Charge (transformation log)",
          x = "Charge_log = log(Charge + 1)",
          y = "Effectifs"
        )
      
      
      # Quantiles
      q_tbl <- function(df, vars) {
        df %>%
          pivot_longer(all_of(vars), names_to = "variable", values_to = "val") %>%
          group_by(variable) %>%
          summarise(
            p01 = quantile(val, 0.01, na.rm = TRUE),
            p05 = quantile(val, 0.05, na.rm = TRUE),
            p50 = quantile(val, 0.50, na.rm = TRUE),
            p95 = quantile(val, 0.95, na.rm = TRUE),
            p99 = quantile(val, 0.99, na.rm = TRUE),
            .groups = "drop"
          )
      }
      
      qt_global <- q_tbl(data_mod, c("AUCg", "DAH_std", "SMOD_2023_std")) %>%
        mutate(across(where(is.numeric), ~ round(., 3)))
      
      qt_global
      
      
      # 1.2 Structure par alpage 
      ggplot(data_mod, aes(x = alpage, y = AUCg))+
        geom_boxplot()
      
      
      ggplot(data_mod, aes(x = alpage, y = SMOD_2023_std))+
        geom_boxplot()
      
      ggplot(data_mod, aes(x = alpage, y = Charge_log))+
        geom_boxplot()
      
    ggplot(data_mod, aes(x = goup =charge_bin, y = AUCg))+
      geom_boxplot()
      
      
      # 1.3 Matrice de corrélation
      
     library(corrplot)
      
      vars <- data_mod %>% 
        select(DAH = DAH_std, FSCA = SMOD_2023_std, Charge_log, AUCg)                               
      
      mat_cor <- cor(vars, use = "complete.obs", method = "spearman")
      print(mat_cor)
      
      corrplot(mat_cor, method = "color", addCoef.col = "black",
               tl.col = "black", tl.cex = 0.9, number.cex = 0.8)
      
      # 1.4 Structure et portée
      
      library(dplyr)
      
      # GLOBAL
      q_charge_global <- data_mod %>%
        filter(charge_bin == 1) %>%
        summarise(
          n    = n(),
          q05  = quantile(Charge_log, 0.05, na.rm = TRUE),
          q50  = quantile(Charge_log, 0.50, na.rm = TRUE),
          q95  = quantile(Charge_log, 0.95, na.rm = TRUE)
        ) %>%
        mutate(
          # (optionnel) revenir à l’échelle Charge brute pour interpréter
          q05_charge = exp(q05) - 1,
          q50_charge = exp(q50) - 1,
          q95_charge = exp(q95) - 1
        )
      q_charge_global
      
      # PAR ALPAGE
      q_charge_par_alpage <- data_mod %>%
        filter(charge_bin == 1) %>%
        group_by(alpage) %>%
        summarise(
          n    = n(),
          q05  = quantile(Charge_log, 0.05, na.rm = TRUE),
          q50  = quantile(Charge_log, 0.50, na.rm = TRUE),
          q95  = quantile(Charge_log, 0.95, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          # (optionnel) aussi sur l’échelle brute
          q05_charge = exp(q05) - 1,
          q50_charge = exp(q50) - 1,
          q95_charge = exp(q95) - 1
        )
      q_charge_par_alpage
      
      # 1.5 Relation préliminaire : 
      
      library(ggplot2)
      
      library(ggplot2)
      library(mgcv)
      
      yl <- quantile(data_mod$AUCg, c(.01, .99), na.rm = TRUE)
      
      # AUCg ~ DAH_std (global)
      p_sc_da <- ggplot(data_mod, aes(DAH_std, AUCg)) +
        geom_point(alpha = 0.08, size = 0.5) +
        geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = TRUE, color = "black") +
        coord_cartesian(ylim = yl) +
        labs(title = "AUCg ~ DAH_std (global)", x = "DAH (standardisé)", y = "AUCg") +
        theme_bw()
      
      # AUCg ~ SMOD_2023_std (global)
      p_sc_sm <- ggplot(data_mod, aes(SMOD_2023_std, AUCg)) +
        geom_point(alpha = 0.08, size = 0.5) +
        geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = TRUE, color = "black") +
        coord_cartesian(ylim = yl) +
        labs(title = "AUCg ~ SMOD_2023_std (global)", x = "Neige (standardisée)", y = "AUCg") +
        theme_bw()
      
      # Par alpage (optionnel)
      p_sc_da_site <- p_sc_da + facet_wrap(~ alpage, ncol = 3, scales = "free_x")
      p_sc_sm_site <- p_sc_sm + facet_wrap(~ alpage, ncol = 3, scales = "free_x")
      
      p_sc_da_site
      p_sc_sm_site
      
      
      
      
      
      
    }
    
    # 2. Modele de base  
    if(TRUE){
      
      library(mgcv)
      
      m_base_lm <- gam(
        AUCg ~ DAH_std +  SMOD_2023_std ,
        data   = data_mod,
        method = "ML"
      )
      
      
      
      m_base<- gam(
        AUCg ~ s(DAH_std, k = 8) +  s(SMOD_2023_std, k = 8),
        data   = data_mod,
        method = "ML"
      )
      
      
      m_base_re <- gam(
        AUCg ~ s(DAH_std, k = 8) +  s(SMOD_2023_std, k = 8) +  s(alpage, bs = "re"),
        data   = data_mod,
        method = "ML"
      )
      
      
      m_base_re_int <- gam(
        AUCg ~ s(DAH_std, k = 8) +                     # effet lissé de la topo
          s(SMOD_2023_std, k = 8) +               # effet lissé de la neige
          ti(DAH_std, SMOD_2023_std, k = c(6, 6)) +  # interaction lisse
          s(alpage, bs = "re"),                   # intercept aléatoire par alpage
        data   = data_mod,
        method = "ML"
      )
      
      AIC(m_base_lm)
      AIC(m_base)
      summary(m_base)
      AIC(m_base_re)
      AIC(m_base_re_int)
      
      # Le dernier model est retenu : non linéaire; alpage en aléatoire; intéraction neige et topo
      
      summary(m_base_re_int)
      
      gam.check(m_base_re_int)
      
      plot(m_base_re_int)
      
      
      # Réajustement du model final avec affinage des nombre de spline
      
      m_base_re_int2 <- gam(
        AUCg ~ s(DAH_std, k = 10) +
          s(SMOD_2023_std, k = 10) +
          ti(DAH_std, SMOD_2023_std, k = c(10, 10)) +  # ou c(10,10) si besoin
          s(alpage, bs = "re"),
        data = data_mod, method = "REML", select = TRUE
      )
      gam.check(m_base_re_int2)
      
      m_base_re_int3 <- gam(
        AUCg ~ s(DAH_std, k=10, bs="ts") +
          s(SMOD_2023_std, k=10, bs="ts") +
          ti(DAH_std, SMOD_2023_std, k=c(10,10)) +
          s(alpage, bs="re"),
        data=data_mod, method="REML", select=TRUE
      )
      gam.check(m_base_re_int3)
      
      m_base_re_int4 <- gam(
        log(AUCg) ~ s(DAH_std, k=10) +
          s(SMOD_2023_std, k=10) +
          ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
          s(alpage, bs="re"),
        data=data_mod, method="ML", select=TRUE
      )
      gam.check(m_base_re_int4)
      
      concurvity(m_base_re_int4, full = TRUE)
      # FOrte concurvité avec le DAH ! donc un model sans ? 
      
      m_base_re_int4 <- gam(
        log(AUCg) ~ s(DAH_std, k=10) +
          s(SMOD_2023_std, k=10) +
          s(alpage, bs="re"),
        data=data_mod, method="ML", select=TRUE
      )
      gam.check(m_base_re_int4)
      summary(m_base_re_int4)
      
      # FORTE PETRE D'AIC
      
      
      ## COMPARAISON ET CHAINE DE MODEL, une intércation SMOD*DAH, l'autre sans
      
      # Chaîne SANS interaction DAH×SMOD (ML)
      m0 <- bam(log(AUCg) ~ s(DAH_std,k=10,bs="ts") + s(SMOD_2023_std,k=10,bs="ts") +
                  s(alpage, bs="re"),
                data=data_mod, method="ML", select=TRUE)   # <- PAS de discrete ici
      m1 <- update(m0, . ~ . + charge_bin)
      m2 <- update(m1, . ~ . + s(Charge_log, k=10, bs="ts", by=charge_bin))
      m3 <- update(m2, . ~ . + ti(Charge_log, SMOD_2023_std, k=c(15,15), by=charge_bin))
      
      # Chaîne AVEC interaction DAH×SMOD dans la base (ML)
      m0_bis <- bam(log(AUCg) ~ s(DAH_std,k=10,bs="ts") + s(SMOD_2023_std,k=10,bs="ts") +
                      ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
                      s(alpage, bs="re"),
                    data=data_mod, method="ML", select=TRUE)
      m1_bis <- update(m0_bis, . ~ . + charge_bin)
      m2_bis <- update(m1_bis, . ~ . + s(Charge_log, k=10, bs="ts", by=charge_bin))
      m3_bis <- update(m2_bis, . ~ . + ti(Charge_log, SMOD_2023_std, k=c(20,20), by=charge_bin))
      
      AIC(m0, m1, m2, m3, m0_bis, m1_bis, m2_bis, m3_bis)   # comparaisons propres en ML
      
      
      summary(m1_bis)
      summary(m2_bis)
      summary(m3_bis)
      gam.check(m3_bis)
      
      
      
      ## PARAMETRAGE DU MODEL FINAL
      
      
      library(mgcv)
      
      m_final <- gam(
        log(AUCg) ~
          # Base topo + neige
          s(DAH_std,       k=10, bs="ts") +
          s(SMOD_2023_std, k=10, bs="ts") +
          ti(DAH_std, SMOD_2023_std, k=c(20,20)) +   # interaction DAH×neige (base)
          # Aléatoire
          s(alpage, bs="re") +
          # Charge (régime chargé uniquement)
          s(Charge_log,                  k=10,     bs="ts", by=charge_bin) +
          ti(Charge_log, SMOD_2023_std,  k=c(15,15),       by=charge_bin),
        data   = data_mod,
        method = "REML",
        select = TRUE
      )
      
      summary(m_final)
      gam.check(m_final)    # k-index ~1 (sinon k=c(12,12) sur Charge/interaction)
      concurvity(m_final, full=TRUE)
      
      
      
      m_hurdle_final <- gam(
        log(AUCg) ~
          # Base topo + neige (avec interaction)
          s(DAH_std,       k = 10, bs = "ts") +
          s(SMOD_2023_std, k = 10, bs = "ts") +
          ti(DAH_std, SMOD_2023_std, k = c(20, 20)) +
          # Aléatoire par alpage
          s(alpage, bs = "re") +
          # --- HURDLE : saut + gradient/interaction si chargé ---
          charge_bin +                                                # SAUT (0 -> 1)
          s(Charge_log, k = 10, bs = "ts", by = charge_bin) +         # GRADIENT si chargé
          ti(Charge_log, SMOD_2023_std, k = c(15, 15), by = charge_bin),  # INTERACTION si chargé
        data   = data_mod,
        method = "REML",
        select = TRUE
      )
      
      summary(m_hurdle_final)
      gam.check(m_hurdle_final)  # vérifier k-index (~1) ; si besoin, passer k=12 sur DAH/Charge
      
      
      
      
      
      
      
      
      m_final_2 <- gam(
        log(AUCg) ~
          # Base topo + neige
          s(dah,       k=10, bs="ts") +
          s(SMOD_2023, k=10, bs="ts") +
          ti(dah, SMOD_2023, k=c(20,20)) +   # interaction DAH×neige (base)
          # Aléatoire
          s(alpage, bs="re") +
          # Charge (régime chargé uniquement)
          charge_bin +
          s(Charge_log,                  k=10,     bs="ts", by=charge_bin) +
          ti(Charge_log, SMOD_2023,  k=c(15,15),       by=charge_bin),
        data   = data_mod,
        method = "REML",
        select = TRUE
      )
      
      summary(m_final_2)
      gam.check(m_final)    # k-index ~1 (sinon k=c(12,12) sur Charge/interaction)
      concurvity(m_final, full=TRUE)
      

      
      
    ## SELECTION DE MODEL POUR LA VERSION PAR ALPAGE DU GRAPH : 
      
      
      
      # =========================
      # Modèle FULL "by alpage"
      # =========================
      library(mgcv); library(dplyr); library(gt)
      
      # Assure-toi que 'alpage' est un facteur
      if (!is.factor(data_mod$alpage)) data_mod$alpage <- factor(data_mod$alpage)
      
      # Gating propre: lissage de charge uniquement si chargé, et séparé par alpage
      # -> on crée une variable facteur avec une modalité "none" comme référence
      data_mod$alp_when_charged <- ifelse(data_mod$charge_bin == 1,
                                          as.character(data_mod$alpage), "none")
      data_mod$alp_when_charged <- factor(data_mod$alp_when_charged,
                                          levels = c("none", levels(data_mod$alpage)))
      
      # Threads (si pas déjà défini)
      n_th <- if (exists("n_th", inherits = TRUE)) n_th else max(1, parallel::detectCores() - 1)
      
      # ---- Tes modèles existants ----
      m_ref <- bam(
        log(AUCg) ~
          s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
          ti(dah, SMOD_2023, k=c(20,20)) +
          s(alpage, bs="re") +
          charge_bin +
          s(Charge_log, k=10, bs="ts", by=charge_bin) +
          ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin),
        data=data_mod, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th
      )
      
      m_ref_plus <- bam(
        log(AUCg) ~
          s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
          ti(dah, SMOD_2023, k=c(20,20)) +
          s(alpage, bs="re") +
          charge_bin +
          s(Charge_log, k=10, bs="ts", by=charge_bin) +
          ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin) +
          s(alpage, bs="re", by=charge_bin) +
          ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
        data=data_mod, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th
      )
      
      m_sobre <- bam(
        log(AUCg) ~
          s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
          ti(dah, SMOD_2023, k=c(15,15)) +
          s(alpage, bs="re") +
          charge_bin +
          s(Charge_log, k=6, bs="cs", by=charge_bin) +
          ti(Charge_log, SMOD_2023, k=c(6,6), by=charge_bin) +
          s(alpage, bs="re", by=charge_bin) +
          s(alpage, by=I(Charge_log*charge_bin), bs="re"),
        data=data_mod, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th, gamma=1.2
      )
      
      # ---- NOUVEAU : FULL by-alpage (toutes les composantes varient par site) ----
      # k modestes + gamma pour contenir les vagues
      m_full_by_alp <- bam(
        log(AUCg) ~
          # topo + neige: courbes propres à chaque alpage
          s(dah,       by = alpage, bs = "fs", k = 7) +
          s(SMOD_2023, by = alpage, bs = "fs", k = 7) +
          ti(dah, SMOD_2023, by = alpage, bs = c("tp","tp"), k = c(5,5)) +
          
          # aléatoire (niveau) commun
          s(alpage, bs = "re") +
          
          # charge: uniquement après seuil, avec formes par alpage
          charge_bin +
          s(Charge_log,             by = alp_when_charged, bs = "fs", k = 6) +
          ti(Charge_log, SMOD_2023, by = alp_when_charged, bs = c("tp","tp"), k = c(5,5)),
        
        data   = data_mod,
        method = "fREML",
        select = TRUE,
        discrete = TRUE,
        nthreads = n_th,
        gamma = 1.2
      )
      
      # =========================
      # Comparaison compacte
      # =========================
      
      cv_rmse <- function(m, dat=data_mod, K=5, seed=42){
        set.seed(seed); id <- sample(rep(1:K, length.out=nrow(dat)))
        rmse <- function(e) sqrt(mean(e^2, na.rm=TRUE)); y <- log(dat$AUCg)
        mean(sapply(1:K, function(k){
          idx <- id==k
          pr <- try(predict(m, newdata=dat[idx,,drop=FALSE], type="response"), silent=TRUE)
          if(inherits(pr,"try-error")) pr <- rep(NA_real_, sum(idx))
          rmse(pr - y[idx])
        }), na.rm=TRUE)
      }
      
      cv_rmse_by_site <- function(m, dat=data_mod, K=5, seed=42){
        set.seed(seed)
        sites <- levels(dat$alpage)
        sapply(sites, function(s){
          ds <- dat[dat$alpage==s,,drop=FALSE]
          id <- sample(rep(1:K, length.out=nrow(ds)))
          y  <- log(ds$AUCg); rmse <- function(e) sqrt(mean(e^2, na.rm=TRUE))
          mean(sapply(1:K, function(k){
            idx <- id==k
            pr <- try(predict(m, newdata=ds[idx,,drop=FALSE], type="response"), silent=TRUE)
            if(inherits(pr,"try-error")) pr <- rep(NA_real_, sum(idx))
            rmse(pr - y[idx])
          }), na.rm=TRUE)
        })
      }
      
      summ_one <- function(m, name){
        s <- summary(m)
        cv_g <- cv_rmse(m)
        cv_sites <- cv_rmse_by_site(m)
        data.frame(
          Model = name,
          n = nobs(m),
          AIC = AIC(m),
          R2 = 100*s$r.sq,
          RMSE_CV = cv_g,
          RMSE_CV_sites_mean = mean(cv_sites),
          stringsAsFactors = FALSE
        )
      }
      
      mods <- list(ref=m_ref, `ref+`=m_ref_plus, sobre=m_sobre, full_by_alp=m_full_by_alp)
      tab <- do.call(rbind, Map(summ_one, mods, names(mods)))
      tab$DeltaAIC <- tab$AIC - min(tab$AIC, na.rm = TRUE)
      
      tab_disp <- tab[,c("Model","n","AIC","DeltaAIC","R2","RMSE_CV","RMSE_CV_sites_mean")]
      gt(tab_disp) |>
        fmt_number(c("AIC","DeltaAIC"), decimals=0) |>
        fmt_number(c("R2"), decimals=1) |>
        fmt_number(c("RMSE_CV","RMSE_CV_sites_mean"), decimals=3) |>
        cols_label(
          Model="Model",
          n="n",
          AIC="AIC",
          DeltaAIC=html("&Delta;AIC"),
          R2=html("R<sup>2</sup> (%)"),
          RMSE_CV="RMSE (CV, log)",
          RMSE_CV_sites_mean="RMSE (CV) mean by site"
        ) |>
        tab_options(table.font.names=c("Inter","Segoe UI","Helvetica","Arial"))
      
      
      
      
      
      m_ref_plus <- bam(
        log(AUCg) ~
          s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
          ti(dah, SMOD_2023, k=c(15,15)) +
          s(alpage, bs="re") +
          charge_bin +
          s(Charge_log, k=10, bs="ts", by=charge_bin) +
          ti(Charge_log, SMOD_2023, k=c(10,10), by=charge_bin) +
          s(alpage, bs="re", by=charge_bin) +
          ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
        data=data_mod, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th
      )
      
      summary(m_ref_plus)
      
      m_by_alp <- m_ref_plus
      m_final_by_alp <- m_ref_plus
      
      library(mgcv)
      library(dplyr)
      
      ## Assure les bons types
      data_mod <- data_mod %>%
        mutate(
          alpage    = factor(alpage),
          charge_bin = factor(ifelse(Charge_log >= log(2), "post", "pre"),
                              levels = c("pre","post"))  # "pre" = référence
        )
      
      m_final_by_alp <- gam(
        log(AUCg) ~
          # Base topo + neige
          s(dah,       k = 10, bs = "ts") +
          s(SMOD_2023, k = 10, bs = "ts") +
          ti(dah, SMOD_2023, k = c(20,20)) +        # interaction DAH × neige (base)
          
          # Effets aléatoires d'interception par alpage
          s(alpage, bs = "re") +
          
          # Saut au seuil (effet paramétrique)
          charge_bin +
          
          # >>> Courbes de charge DIFFÉRENTES par alpage, uniquement côté "post"
          # (factor-smooth: un lissage par niveau d'alpage, partageant λ ; zéro côté "pre")
          s(Charge_log, alpage, bs = "fs", k = 10, by = charge_bin) +
          
          # Interaction charge × neige (commune), active seulement côté "post"
          ti(Charge_log, SMOD_2023, k = c(15,15), by = charge_bin),
        
        data   = data_mod,
        method = "REML",
        select = TRUE
      )
      
      
      
      # PLot FOcus 2023 by snow class 
      if(TRUE){
      
      ## FINAL COMPLET : 
      
      # =========================================================
      # Libraries
      # =========================================================
      
      
      library(dplyr)
      library(tidyr)
      library(purrr)
      library(tibble)
      library(ggplot2)
      library(conflicted); conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::summarise, dplyr::mutate, dplyr::arrange, dplyr::count)
      # =========================================================
      # Look & feel
      # =========================================================
      base_family <- "Segoe UI"   # ou "Inter"/"Helvetica" selon ta machine# Palette des 3 classes (labels EXACTS)
      pal <- c("Déneigement : Précoce (≤120)" = "#fed976",
               "Déneigement : Moyen (121–150)"= "#41b6c4",
               "Déneigement : Tardif (≥151)"  = "#0c2c84")
      
      # =========================================================
      # Paramètres de support / pointillé
      # =========================================================
      x_thr    <- log(2)     # seuil Charge=1
      width_x  <- 0.25       # fenêtre locale (Charge_log)
      n_min_x_class    <- 100
      prop_min_x_class <- 0.005
      x_dotted_min     <- 5.5 # on n'autorise le pointillé qu'à forte charge
      
      # =========================================================
      # Termes à exclure
      # =========================================================
      term_names <- colnames(predict(m_by_alp, type = "terms"))
      exc_dah    <- term_names[grepl("dah", term_names)]
      
      # =========================================================
      # Préparation des données
      # =========================================================
      rng_by_alp <- data_mod %>%
        group_by(alpage) %>%
        summarise(ch_lo = quantile(Charge_log, .01, na.rm = TRUE),
                  ch_hi = quantile(Charge_log, .99, na.rm = TRUE),
                  .groups = "drop")
      
      used_charged <- data_mod %>%
        inner_join(rng_by_alp, by = "alpage") %>%
        filter(between(Charge_log, ch_lo, ch_hi), charge_bin == 1) %>%
        dplyr::select(alpage, SMOD_2023, Charge_log, AUCg, ch_lo, ch_hi)
      
      # Classes SMOD fixes (+10 j)
      class_bounds <- tibble(
        neige = factor(c("Déneigement : Précoce (≤120)",
                         "Déneigement : Moyen (121–150)",
                         "Déneigement : Tardif (≥151)"),
                       levels = c("Déneigement : Précoce (≤120)",
                                  "Déneigement : Moyen (121–150)",
                                  "Déneigement : Tardif (≥151)")),
        z_min = c(-Inf, 121, 151),
        z_max = c(120, 150,  Inf)
      )
      
      # Support global par classe
      has_support_class <- function(alp, zmin, zmax, nmin = 200, pmin = 0.01){
        d_alp <- dplyr::filter(used_charged, alpage == alp)
        n_tot <- nrow(d_alp)
        d_cls <- dplyr::filter(d_alp, SMOD_2023 >= zmin, SMOD_2023 <= zmax)
        n_cls <- nrow(d_cls)
        (n_cls >= nmin) && (n_cls / pmax(n_tot,1) >= pmin)
      }
      
      # Support local le long de x dans une classe
      support_along_x <- function(d_cls, x_seq){
        n_cls <- nrow(d_cls)
        if (n_cls == 0) return(list(n_loc = rep(0L, length(x_seq)),
                                    prop  = rep(0,   length(x_seq))))
        n_loc <- sapply(x_seq, function(x) sum(abs(d_cls$Charge_log - x) <= width_x, na.rm = TRUE))
        list(n_loc = n_loc, prop = n_loc / n_cls)
      }
      
      # =========================================================
      # Courbes marginalisées (pondérées par la distribution réelle de SMOD)
      # =========================================================
      curves_excl <- tibble(alpage = levels(data_mod$alpage)) %>%
        left_join(rng_by_alp, by = "alpage") %>%
        group_split(alpage) %>%
        map_dfr(function(lim){
          a     <- lim$alpage[[1]]
          x_seq <- seq(lim$ch_lo[[1]], lim$ch_hi[[1]], length.out = 200)
          
          map_dfr(seq_len(nrow(class_bounds)), function(i){
            lab_i <- class_bounds$neige[i]
            zmin  <- class_bounds$z_min[i]
            zmax  <- class_bounds$z_max[i]
            
            d_cls <- dplyr::filter(used_charged, alpage == a,
                                   SMOD_2023 >= zmin, SMOD_2023 <= zmax)
            if (!has_support_class(a, zmin, zmax)) return(tibble())
            
            z_grid <- sort(unique(round(d_cls$SMOD_2023)))
            if (length(z_grid) > 60) z_grid <- z_grid[round(seq(1, length(z_grid), length.out = 60))]
            
            w_tbl <- as.data.frame(table(round(d_cls$SMOD_2023)))
            names(w_tbl) <- c("SMOD_2023","w")
            w_tbl$SMOD_2023 <- as.numeric(as.character(w_tbl$SMOD_2023))
            w_tbl <- dplyr::filter(w_tbl, SMOD_2023 %in% z_grid)
            
            nd <- tidyr::expand_grid(
              alpage     = factor(a, levels = levels(data_mod$alpage)),
              Charge_log = x_seq,
              SMOD_2023  = z_grid
            ) %>%
              left_join(w_tbl, by = "SMOD_2023") %>%
              mutate(
                neige      = lab_i,
                dah        = median(data_mod$dah[data_mod$alpage == a], na.rm = TRUE),
                charge_bin = as.integer(Charge_log >= x_thr),
                side = factor(ifelse(Charge_log < x_thr, "Avant seuil (non pâturé)", "Après seuil (pâturé)"),
                              levels = c("Avant seuil (non pâturé)", "Après seuil (pâturé)"))
              )
            
            pr <- predict(m_by_alp, newdata = nd, type = "link", se.fit = TRUE, exclude = exc_dah)
            out <- tibble(
              alpage = a, neige = lab_i,
              Charge_log = nd$Charge_log, side = nd$side,
              SMOD_2023 = nd$SMOD_2023, w = nd$w,
              AUC    = exp(pr$fit),
              AUC_lo = exp(pr$fit - 1.96*pr$se.fit),
              AUC_hi = exp(pr$fit + 1.96*pr$se.fit)
            )
            
            out_sum <- out %>%
              group_by(alpage, neige, side, Charge_log) %>%
              summarise(
                wsum   = sum(w, na.rm = TRUE),
                AUC    = sum(w*AUC,    na.rm = TRUE) / pmax(wsum, 1),
                AUC_lo = sum(w*AUC_lo, na.rm = TRUE) / pmax(wsum, 1),
                AUC_hi = sum(w*AUC_hi, na.rm = TRUE) / pmax(wsum, 1),
                .groups = "drop"
              )
            
            sup <- support_along_x(d_cls, x_seq)
            out_sum$n_loc <- sup$n_loc
            out_sum$prop  <- sup$prop
            
            # Pointillé uniquement "dans la queue", après seuil et à forte charge
            after_thr <- out_sum$Charge_log >= x_thr & out_sum$Charge_log >= x_dotted_min
            ok_local  <- out_sum$n_loc >= n_min_x_class & out_sum$prop >= prop_min_x_class
            last_ok   <- if (any(after_thr & ok_local)) max(which(after_thr & ok_local)) else -Inf
            
            out_sum$lt <- factor(
              ifelse(out_sum$Charge_log >= x_thr &
                       out_sum$Charge_log >= x_dotted_min &
                       seq_len(nrow(out_sum)) > last_ok, "dotted", "solid"),
              levels = c("solid","dotted")
            )
            out_sum
          })
        }) %>%
        arrange(alpage, neige, side, Charge_log)
      
      # =========================================================
      # Lignes de référence (baseline pré-seuil prolongée)
      # =========================================================
      pre_df <- filter(curves_excl, side == "Avant seuil (non pâturé)", Charge_log <= x_thr)
      post_df <- filter(curves_excl, side == "Après seuil (pâturé)",     Charge_log >  x_thr)
      
      witness <- pre_df %>%
        group_by(alpage, neige) %>%
        summarise(y_thr = approx(Charge_log, AUC, xout = x_thr, rule = 2, ties = "ordered")$y,
                  .groups = "drop") %>%
        left_join(post_df %>% group_by(alpage, neige) %>% summarise(x_end = max(Charge_log), .groups = "drop"),
                  by = c("alpage","neige")) %>%
        mutate(x_start = x_thr) %>%
        filter(!is.na(x_end), x_end > x_start)
      
      # =========================================================
      # Comptes n par (alpage × classe) + positions par défaut
      # =========================================================
      counts_df <- used_charged %>%
        mutate(
          neige = dplyr::case_when(
            SMOD_2023 <= 120 ~ "Déneigement : Précoce (≤120)",
            SMOD_2023 <= 150 ~ "Déneigement : Moyen (121–150)",
            TRUE             ~ "Déneigement : Tardif (≥151)"
          ),
          neige = factor(neige, levels = levels(class_bounds$neige))
        ) %>%
        count(alpage, neige, name = "n_pix") %>%
        left_join(rng_by_alp, by = "alpage") %>%
        mutate(
          x_pos = ch_lo + 0.12*(ch_hi - ch_lo),
          y_pos = c(98, 94, 90)[as.integer(neige)]
        )
      
      # Déplacer les 'n =' du panneau Viso (droite, sous la légende)
      counts_df2 <- counts_df %>%
        mutate(
          span    = ch_hi - ch_lo,
          x_pos2  = if_else(alpage == "Viso", ch_hi - 0.02*span, x_pos),
          y_pos2  = if_else(alpage == "Viso", y_pos - 20,        y_pos),
          hjust2  = if_else(alpage == "Viso", 1,                 0)
        )
      
      # =========================================================
      # Points observés (fond)
      # =========================================================
      obs_curve <- data_mod %>%
        inner_join(rng_by_alp, by = "alpage") %>%
        filter(between(Charge_log, ch_lo, ch_hi)) %>%
        select(alpage, Charge_log, AUCg)
      
      # =========================================================
      # Légende compacte (sans les 3 pointillés colorés)
      # =========================================================
      class_lvls <- levels(class_bounds$neige)
      legend_breaks <- c(
        class_lvls,
        "No-load reference level (pre-threshold baseline)",
        "Loading threshold"
      )
      legend_labels <- c(
        "Déneigement : Précoce (≤120)" = "Snowmelt Early (≤120)",
        "Déneigement : Moyen (121–150)"= "Snowmelt Mid (121–150)",
        "Déneigement : Tardif (≥151)"  = "Snowmelt Late (≥151)",
        "No-load reference level (pre-threshold baseline)" = "No-load reference level",
        "Loading threshold"                                = "Loading threshold"
      )
      legend_cols <- c(
        pal,
        "No-load reference level (pre-threshold baseline)" = "grey40",
        "Loading threshold"                                = "black"
      )
      
      # Pour que la scale de couleur connaisse aussi les 2 niveaux “ref” & “threshold”
      curves_excl$neige <- factor(curves_excl$neige,
                                  levels = c(class_lvls,
                                             "No-load reference level (pre-threshold baseline)",
                                             "Loading threshold"))
      
      legend_ref_key <- data.frame(
        x = c(0,1), y = 0,
        neige = factor("No-load reference level (pre-threshold baseline)", levels = legend_breaks)
      )
      legend_thr_key <- data.frame(
        x = c(0,1), y = 0,
        neige = factor("Loading threshold", levels = legend_breaks)
      )
      
      # Séparation des segments à tracer
      curves_solid  <- subset(curves_excl, lt == "solid")
      curves_dotted <- subset(curves_excl, lt == "dotted")
      
      # =========================================================
      # PLOT
      # =========================================================
      p <- ggplot() +
        geom_point(data = obs_curve, aes(Charge_log, AUCg),
                   colour = "grey60", alpha = 0.10, size = 0.25) +
        geom_ribbon(data = curves_excl,
                    aes(Charge_log, ymin = AUC_lo, ymax = AUC_hi,
                        fill = neige, group = interaction(alpage, neige, side)),
                    alpha = 0.15, colour = NA, show.legend = FALSE) +
        
        geom_line(data = curves_solid,
                  aes(Charge_log, AUC, colour = neige,
                      group = interaction(alpage, neige, side)),
                  linewidth = 1, linetype = "solid", lineend = "round", show.legend = TRUE) +
        geom_line(data = curves_dotted,
                  aes(Charge_log, AUC, colour = neige,
                      group = interaction(alpage, neige, side)),
                  linewidth = 1, linetype = "dotted", lineend = "round", show.legend = FALSE) +
        
        # Référence (gris) — tracée mais hors légende
        geom_segment(data = witness,
                     aes(x = x_start, xend = x_end, y = y_thr, yend = y_thr),
                     colour = "grey40", linetype = "longdash", linewidth = 1, alpha = 0.95,
                     show.legend = FALSE) +
        
        # Seuil vertical
        geom_vline(xintercept = x_thr, linetype = "dotted",
                   colour = "black", linewidth = 1, alpha = 0.98, show.legend = FALSE) +
        
        # n = (plus gros & gras ; Viso déplacé)
        geom_text(data = counts_df2,
                  aes(x = x_pos2, y = y_pos2, label = paste0("n = ", n_pix),
                      colour = neige, hjust = hjust2),
                  size = 4.2, fontface = "bold", show.legend = FALSE) +
        
        # Dummies pour “ref” et “threshold” dans la légende
        geom_line(data = legend_ref_key,
                  aes(x, y, colour = neige),
                  linetype = "longdash", linewidth = 0.8, lineend = "butt",
                  alpha = 0, inherit.aes = FALSE, show.legend = TRUE) +
        geom_line(data = legend_thr_key,
                  aes(x, y, colour = neige),
                  linetype = "dotted", linewidth = 0.7, lineend = "butt",
                  alpha = 0, inherit.aes = FALSE, show.legend = TRUE) +
        
        scale_colour_manual(
          values = legend_cols,
          breaks = legend_breaks,
          labels = legend_labels,
          name   = "Snowmelt class — Line type",
          drop   = FALSE
        ) +
        scale_fill_manual(values = pal, guide = "none") +
        
        guides(
          # 5 clés : 3 solides (classes), ref longdash, seuil dotted
          colour = guide_legend(
            title.position = "top",
            keywidth = unit(42, "pt"),
            override.aes = list(
              alpha     = 1,
              linewidth = c(0.7,0.7,0.7, 0.8,0.7),
              linetype  = c("solid","solid","solid", "longdash","dotted")
            )
          )
        ) +
        
        facet_wrap(~ alpage, ncol = 3) +
        coord_cartesian(ylim = c(0, 100), clip = "off") +
        scale_x_continuous("Stocking rate (log transform)",
                           expand = expansion(mult = c(0.005, 0.02))) +
        scale_y_continuous("Predicted growth production (GPROD)",
                           expand = expansion(mult = c(0.02, 0.08))) +
        
        theme_minimal(base_size = 12, base_family = base_family) +
        theme(
          legend.position      = c(0.985, 0.98),
          legend.justification = c(1, 1),
          legend.direction     = "vertical",
          legend.background    = element_rect(fill = scales::alpha("white", 0.95), colour = "grey80"),
          legend.key.height    = unit(11, "pt"),
          legend.title         = element_text(size = 11, face = "bold"),
          legend.text          = element_text(size = 10.5),
          
          axis.text.x          = element_text(size = 11),
          axis.text.y          = element_text(size = 11),
          axis.title.x         = element_text(size = 14, margin = margin(t = 6)),
          axis.title.y         = element_text(size = 14, margin = margin(r = 6)),
          
          panel.grid.minor     = element_blank(),
          panel.grid.major.x   = element_line(linewidth = 0.25),
          panel.grid.major.y   = element_line(linewidth = 0.25),
          panel.border         = element_rect(colour = "black", fill = NA, linewidth = 0.9),
          
          strip.background     = element_rect(fill = "#f7f7f7", colour = "grey30", linewidth = 0.6),
          strip.text           = element_text(face = "bold", size = 14.5),
          
          plot.margin          = margin(4, 8, 2, 8)
        )
      
      print(p)
      
      # =========================================================
      # Export (PNG)
      # =========================================================
      if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png("Figure_5_Resultat_M2_Predict_model.svg",
                      width = 12.5, height = 6, units = "in", res = 450, scaling = 1)
        print(p); dev.off()
      } else {
        ggsave("Figure_5_Resultat_M2_Predict_model.svg",
               p, width = 12.5, height = 6, dpi = 450, device = cairo_png)
      }
      
      if (requireNamespace("svglite", quietly = TRUE)) {
        svglite::svglite("Figure_5_Resultat_M2_Predict_model.svg",
                         width = 12.5, height = 6)
        print(p); dev.off()
      } else {
        ggsave("Figure_5_Resultat_M2_Predict_model.svg",
               p, width = 12.5, height = 6, dpi = 450, device = cairo_svg)
      }
      
      ## V2 : LEGENDE COMPLETE
      
      
      
      # ===== Counts: move labels for Viso under the legend (right side) =====
      counts_df2 <- counts_df |>
        dplyr::mutate(
          span    = ch_hi - ch_lo,
          # Viso: coller à droite (2% du bord) et descendre sous la légende (~10 unités)
          x_pos2  = dplyr::if_else(alpage == "Viso", ch_hi - 0.02*span, x_pos),
          y_pos2  = dplyr::if_else(alpage == "Viso", y_pos - 30,          y_pos),
          hjust2  = dplyr::if_else(alpage == "Viso", 1,                   0)
        )
      
      # ===== Legend content (single block, perfect patterns) =====
      class_lvls <- levels(class_bounds$neige)
      low_lvls   <- paste0(class_lvls, " — prediction-uncertain (few local points)")
      
      legend_breaks <- c(
        class_lvls, low_lvls,
        "No-load reference level (pre-threshold baseline)",
        "Loading threshold"
      )
      
      legend_labels <- c(
        "Déneigement : Précoce (≤120)" = "Snowmelt Early (≤120)",
        "Déneigement : Moyen (121–150)"= "Snowmelt Mid (121–150)",
        "Déneigement : Tardif (≥151)"  = "Snowmelt Late (≥151)",
        "Déneigement : Précoce (≤120) — prediction-uncertain (few local points)" =
          "Prediction-uncertain",
        "Déneigement : Moyen (121–150) — prediction-uncertain (few local points)"  =
          "Prediction-uncertain",
        "Déneigement : Tardif (≥151) — prediction-uncertain (few local points)"   =
          "Prediction-uncertain",
        "No-load reference level (pre-threshold baseline)" = "No-load reference level",
        "Loading threshold"                                = "Loading threshold"
      )
      
      legend_cols <- c(
        pal,
        stats::setNames(unname(pal), low_lvls),
        "No-load reference level (pre-threshold baseline)" = "grey40",
        "Loading threshold"                                = "black"
      )
      
      # Make sure the scale knows all levels
      curves_excl$neige <- factor(curves_excl$neige, levels = legend_breaks)
      
      # Legend-only dummies (don’t draw on panels)
      legend_low_keys <- data.frame(
        x = rep(c(0,1), 3), y = 0,
        neige = factor(rep(low_lvls, each = 2), levels = legend_breaks)
      )
      legend_ref_key <- data.frame(
        x = c(0,1), y = 0,
        neige = factor("No-load reference level (pre-threshold baseline)", levels = legend_breaks)
      )
      legend_thr_key <- data.frame(
        x = c(0,1), y = 0,
        neige = factor("Loading threshold", levels = legend_breaks)
      )
      
      # Split lines
      curves_solid  <- subset(curves_excl,  lt == "solid")
      curves_dotted <- subset(curves_excl,  lt == "dotted")
      
      # ===== PLOT =====
      p <- ggplot() +
        geom_point(data = obs_curve, aes(Charge_log, AUCg),
                   colour = "grey60", alpha = 0.10, size = 0.25) +
        geom_ribbon(data = curves_excl,
                    aes(Charge_log, ymin = AUC_lo, ymax = AUC_hi,
                        fill = neige, group = interaction(alpage, neige, side)),
                    alpha = 0.15, colour = NA, show.legend = FALSE) +
        
        # Curves (do not change linewidth/linetype)
        geom_line(data = curves_solid,
                  aes(Charge_log, AUC, colour = neige,
                      group = interaction(alpage, neige, side)),
                  linewidth = 1, linetype = "solid", lineend = "round", show.legend = TRUE) +
        geom_line(data = curves_dotted,
                  aes(Charge_log, AUC, colour = neige,
                      group = interaction(alpage, neige, side)),
                  linewidth = 1, linetype = "dotted", lineend = "round", show.legend = FALSE) +
        
        # Reference (drawn)
        geom_segment(data = witness,
                     aes(x = x_start, xend = x_end, y = y_thr, yend = y_thr),
                     colour = "grey40", linetype = "longdash", linewidth = 1, alpha = 0.95) +
        
        # Threshold (drawn)
        geom_vline(xintercept = x_thr, linetype = "dotted",
                   colour = "black", linewidth = 1, alpha = 0.98, show.legend = FALSE) +
        
        # n = labels (bigger + bold; Viso moved under legend at right)
        geom_text(data = counts_df2,
                  aes(x = x_pos2, y = y_pos2, label = paste0("n = ", n_pix),
                      colour = neige, hjust = hjust2),
                  size = 4.2, fontface = "bold", show.legend = FALSE) +
        
        # --- Legend dummies (thin, crisp patterns) ---
        geom_line(
          data = legend_low_keys, aes(x, y, colour = neige),
          linetype = "dotted", linewidth = 0.7, lineend = "butt",
          alpha = 0, inherit.aes = FALSE, show.legend = TRUE
        ) +
        geom_line(
          data = legend_ref_key, aes(x, y, colour = neige),
          linetype = "longdash", linewidth = 0.8, lineend = "butt",
          alpha = 0, inherit.aes = FALSE, show.legend = TRUE
        ) +
        geom_line(
          data = legend_thr_key, aes(x, y, colour = neige),
          linetype = "dotted", linewidth = 0.7, lineend = "butt",
          alpha = 0, inherit.aes = FALSE, show.legend = TRUE
        ) +
        
        scale_colour_manual(
          values = legend_cols,
          breaks = legend_breaks,
          labels = legend_labels,
          name   = "Snowmelt class — Line type",
          drop   = FALSE
        ) +
        scale_fill_manual(values = pal, guide = "none") +
        
        # Override EXACT for the 8 legend keys (solid, 3 dotted, ref longdash, threshold dotted)
        guides(
          colour = guide_legend(
            title.position = "top",
            keywidth = unit(42, "pt"),
            override.aes = list(
              alpha     = 1,
              linewidth = c(0.7,0.7,0.7, 0.7,0.7,0.7, 0.8,0.7),
              linetype  = c("solid","solid","solid", "dotted","dotted","dotted", "longdash","dotted")
            )
          )
        ) +
        
        facet_wrap(~ alpage, ncol = 3) +
        coord_cartesian(ylim = c(0, 100), clip = "off") +
        scale_x_continuous("Stocking rate (log transform)",
                           expand = expansion(mult = c(0.005, 0.02))) +
        scale_y_continuous("Predicted growth production (GPROD)",
                           expand = expansion(mult = c(0.02, 0.08))) +
        
        theme_minimal(base_size = 12, base_family = base_family) +
        theme(
          legend.position      = c(0.985, 0.98),
          legend.justification = c(1, 1),
          legend.direction     = "vertical",
          legend.background    = element_rect(fill = scales::alpha("white", 0.95), colour = "grey80"),
          legend.key.height    = unit(11, "pt"),
          legend.title         = element_text(size = 11, face = "bold"),
          legend.text          = element_text(size = 10.5),
          
          axis.text.x          = element_text(size = 11),
          axis.text.y          = element_text(size = 11),
          axis.title.x         = element_text(size = 14, margin = margin(t = 6)),
          axis.title.y         = element_text(size = 14, margin = margin(r = 6)),
          
          panel.grid.minor     = element_blank(),
          panel.grid.major.x   = element_line(linewidth = 0.25),
          panel.grid.major.y   = element_line(linewidth = 0.25),
          panel.border         = element_rect(colour = "black", fill = NA, linewidth = 0.9),
          
          strip.background     = element_rect(fill = "#f7f7f7", colour = "grey30", linewidth = 0.6),
          strip.text           = element_text(face = "bold", size = 14.5),
          
          plot.margin          = margin(4, 8, 2, 8)
        )
      
      print(p)
      
      # ===== Export (sharp text) =====
      if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png("Figure_5_Resultat_M2_Predict_model_V2.png",
                      width = 12.5, height = 6, units = "in", res = 450, scaling = 1)
        print(p); dev.off()
      } else {
        ggsave("Figure_5_Resultat_M2_Predict_model_V2.png",
               p, width = 12.5, height = 6, dpi = 450, device = cairo_png)
      }
      
      }
      
      
      
      # CARTE DE CHALEUR
       if (F){
      library(dplyr)
      library(tidyr)
      library(purrr)
      library(ggplot2)
      library(mgcv)
      
      x_thr <- log(2)
      
      ## 1) Termes contenant 'dah' à exclure des prédictions
      term_names <- colnames(predict(m_by_alp, type = "terms"))
      exc_dah    <- term_names[grepl("dah", term_names)]  # s(dah:...), ti(dah,SMOD_2023:...), etc.
      
      ## 2) Bornes (SMOD et Charge) par alpage pour tracer sur du réel
      rng_by_alp <- data_mod %>%
        group_by(alpage) %>%
        summarise(
          smod_lo = quantile(SMOD_2023, .01, na.rm = TRUE),
          smod_hi = quantile(SMOD_2023, .99, na.rm = TRUE),
          ch_lo   = quantile(Charge_log, .01, na.rm = TRUE),
          ch_hi   = quantile(Charge_log, .99, na.rm = TRUE),
          .groups = "drop"
        )
      
      ## 3) Grille régulière par alpage
      surf_all <- rng_by_alp %>%
        group_split(alpage) %>%
        map_dfr(function(lim){
          a     <- lim$alpage[[1]]
          x_seq <- seq(lim$ch_lo[[1]],   lim$ch_hi[[1]],   length.out = 120)
          z_seq <- seq(lim$smod_lo[[1]], lim$smod_hi[[1]], length.out = 80)
          
          expand_grid(
            alpage     = factor(a, levels = levels(data_mod$alpage)),
            SMOD_2023  = z_seq,
            Charge_log = x_seq
          ) %>%
            mutate(
              ## valeur quelconque pour 'dah' (les termes seront EXCLUS à la prédiction)
              dah              = 0,
              charge_bin       = as.integer(Charge_log >= x_thr),
              alp_when_charged = factor(ifelse(charge_bin == 1, as.character(alpage), "none"),
                                        levels = c("none", levels(data_mod$alpage)))
            )
        })
      
      ## 4) Prédictions SANS l’effet de dah
      pr <- predict(m_by_alp, newdata = surf_all, type = "link", se.fit = TRUE,
                    exclude = exc_dah)
      
      surf_all <- surf_all %>%
        mutate(
          AUC    = exp(pr$fit),
          AUC_lo = exp(pr$fit - 1.96 * pr$se.fit),
          AUC_hi = exp(pr$fit + 1.96 * pr$se.fit)
        )
      
      ## 5) Points observés (pour contexte visuel)
      obs <- data_mod %>%
        inner_join(rng_by_alp, by = "alpage") %>%
        filter(between(Charge_log, ch_lo, ch_hi),
               between(SMOD_2023,  smod_lo, smod_hi)) %>%
        select(alpage, SMOD_2023, Charge_log)
      
      ## 6) Palette contrôlable par l’utilisateur (remplace 'pal_grad' par ce que tu veux)
      library(ggplot2)
      library(grid)  # pour unit()
      
      # ⇩ Choisis la palette que tu veux (ordre = bas → haut d’AUC)
      pal_grad <- c("#0b1026", "#0e3a5e", "#1976a1", "#46b5b5", "#a8eb12", "#ffca3a")
      
      rng_auc <- range(surf_all$AUC, na.rm = TRUE)
      
      ggplot() +
        # tuiles
        geom_raster(data = surf_all,
                    aes(Charge_log, SMOD_2023, fill = AUC)) +
        # isocontours
        geom_contour(data = surf_all,
                     aes(Charge_log, SMOD_2023, z = AUC),
                     colour = "white", bins = 12, linewidth = 0.28, alpha = .75) +
        # points d’observation (très discrets, au-dessus)
        geom_point(data = obs,
                   aes(Charge_log, SMOD_2023),
                   colour = "white", alpha = 0.06, size = 0.15) +
        # seuil
        geom_vline(xintercept = x_thr, linetype = "22", colour = "white", linewidth = .4, alpha = .7) +
        
        # palette contrôlable
        scale_fill_gradientn(colours = pal_grad, limits = rng_auc,
                             oob = scales::squish, name = "AUC prédite") +
        
        facet_wrap(~ alpage, ncol = 3) +
        coord_cartesian(expand = FALSE) +
        labs(title    = "AUC prédite (Charge × Neige) — effet de « dah » exclu",
             subtitle = "Pointillé : seuil de charge = 1 (log(2)).",
             x = "Charge_log = log(Charge + 1)", y = "SMOD_2023 (jour)") +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid       = element_blank(),
          panel.background = element_rect(fill = "#0b0f1a", colour = NA),
          plot.background  = element_rect(fill = "#0b0f1a", colour = NA),
          axis.text        = element_text(colour = "#e5e7eb"),
          axis.title       = element_text(colour = "#e5e7eb"),
          strip.background = element_rect(fill = "#111827", colour = NA),
          strip.text       = element_text(colour = "white", face = "bold"),
          legend.position  = "right",
          legend.text      = element_text(colour = "#e5e7eb"),
          legend.title     = element_text(colour = "#e5e7eb")
        ) +
        guides(fill = guide_colorbar(title.position = "top",
                                     barwidth = unit(5.5, "cm"),
                                     barheight = unit(0.4, "cm")))
      
      
    
      
    
      
      
      
      
      
       }
      
      
      # Chargement des librairies nécessaires
      library(mgcv)
      library(knitr)
      library(kableExtra)
      library(dplyr)
      library(webshot2)  # Pour l'export PNG
      
      # Votre modèle GAM
      m_ref_plus <- bam(
        log(AUCg) ~ s(dah, k=10, bs="ts") + 
          s(SMOD_2023, k=10, bs="ts") + 
          ti(dah, SMOD_2023, k=c(15,15)) + 
          s(alpage, bs="re") + 
          charge_bin + 
          s(Charge_log, k=10, bs="ts", by=charge_bin) + 
          ti(Charge_log, SMOD_2023, k=c(10,10), by=charge_bin) + 
          s(alpage, bs="re", by=charge_bin) + 
          ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
        data = data_mod, 
        method = "fREML", 
        select = TRUE, 
        discrete = TRUE, 
        nthreads = n_th
      )
      
      # Extraction des résultats du modèle
      model_summary <- summary(m_ref_plus)
      
      # Fonction pour extraire et formatter les résultats
      extract_gam_results <- function(model_summary) {
        # Termes paramétriques
        parametric <- data.frame(
          Term = rownames(model_summary$p.table),
          Estimate = model_summary$p.table[, "Estimate"],
          SE = model_summary$p.table[, "Std. Error"],
          t_value = model_summary$p.table[, "t value"],
          p_value = model_summary$p.table[, "Pr(>|t|)"],
          Type = "Paramétrique"
        )
        
        # Termes smooth
        smooth <- data.frame(
          Term = rownames(model_summary$s.table),
          Estimate = NA,  # Les termes smooth n'ont pas d'estimate direct
          SE = NA,
          t_value = model_summary$s.table[, "F"],  # F-statistic pour les smooth
          p_value = model_summary$s.table[, "p-value"],
          Type = "Smooth"
        )
        
        # Combinaison des résultats
        results <- rbind(parametric, smooth)
        
        # Ajout de la significativité
        results$Significance <- case_when(
          results$p_value < 0.001 ~ "***",
          results$p_value < 0.01 ~ "**",
          results$p_value < 0.05 ~ "*",
          results$p_value < 0.1 ~ ".",
          TRUE ~ ""
        )
        
        return(results)
      }
      
      # Extraction des résultats
      results_table <- extract_gam_results(model_summary)
      
      # Formatage de la table pour présentation
      formatted_table <- results_table %>%
        mutate(
          Term = case_when(
            Term == "(Intercept)" ~ "Intercept",
            Term == "charge_bin1" ~ "Charge binaire (Oui)",
            grepl("^s\\(", Term) ~ gsub("s\\(([^)]+)\\)", "\\1 (smooth)", Term),
            grepl("^ti\\(", Term) ~ gsub("ti\\(([^)]+)\\)", "\\1 (interaction)", Term),
            TRUE ~ Term
          ),
          Estimate = ifelse(is.na(Estimate), "—", sprintf("%.4f", Estimate)),
          SE = ifelse(is.na(SE), "—", sprintf("%.4f", SE)),
          Statistic = ifelse(Type == "Paramétrique", 
                             sprintf("%.3f", t_value), 
                             sprintf("%.3f", t_value)),
          `P-value` = case_when(
            p_value < 0.001 ~ "< 0.001",
            p_value < 0.01 ~ sprintf("%.3f", p_value),
            TRUE ~ sprintf("%.3f", p_value)
          )
        ) %>%
        select(Term, Type, Estimate, SE, Statistic, `P-value`, Significance)
      
      # Création de la table avec kable et kableExtra
      beautiful_table <- formatted_table %>%
        kbl(
          caption = "Résultats du modèle GAM : log(AUCg)",
          col.names = c("Terme", "Type", "Estimation", "Erreur Standard", 
                        "Statistique", "Valeur-p", "Signif."),
          align = c("l", "c", "r", "r", "r", "r", "c")
        ) %>%
        kable_styling(
          bootstrap_options = c("striped", "hover", "condensed", "responsive"),
          full_width = FALSE,
          position = "center",
          font_size = 12
        ) %>%
        row_spec(0, bold = TRUE, background = "#34495e", color = "white") %>%
        column_spec(1, bold = TRUE, width = "3cm") %>%
        column_spec(2, width = "2cm") %>%
        column_spec(c(3:6), width = "1.5cm") %>%
        column_spec(7, width = "1cm") %>%
        pack_rows("Termes paramétriques", 1, sum(results_table$Type == "Paramétrique"), 
                  label_row_css = "background-color: #3498db; color: white;") %>%
        pack_rows("Termes smooth", sum(results_table$Type == "Paramétrique") + 1, nrow(results_table),
                  label_row_css = "background-color: #e74c3c; color: white;") %>%
        footnote(
          general = c("Codes de significativité : 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
                      paste("R² ajusté =", round(model_summary$r.sq, 3)),
                      paste("R² expliqué =", round(model_summary$dev.expl, 3)),
                      paste("AIC =", round(AIC(m_ref_plus), 1)),
                      paste("n =", nobs(m_ref_plus))),
          general_title = "Notes :",
          footnote_as_chunk = TRUE
        )
      
      # Affichage de la table
      print(beautiful_table)
      
      # Export en PNG
      # Sauvegarde temporaire en HTML
      temp_html <- tempfile(fileext = ".html")
      save_kable(beautiful_table, file = temp_html)
      
      # Conversion en PNG
      webshot2::webshot(
        url = temp_html,
        file = "gam_results_table.png",
        vwidth = 1000,
        vheight = 800,
        zoom = 2,
        delay = 1
      )
      
      # Nettoyage
      unlink(temp_html)
      
      cat("Table exportée en PNG : gam_results_table.png\n")
      
      # Informations supplémentaires sur le modèle
      cat("\nInformations supplémentaires du modèle :\n")
      cat("Deviance expliquée :", round(model_summary$dev.expl * 100, 2), "%\n")
      cat("R² ajusté :", round(model_summary$r.sq, 3), "\n")
      cat("AIC :", round(AIC(m_ref_plus), 1), "\n")
      cat("Nombre d'observations :", nobs(m_ref_plus), "\n"
      
      
      
      
      
      
      
      
      
      
      
      
      ## TAILLE d'EFFET : 
      
      # =========================================================
      # Taille d'effet par alpage × régime de neige
      #  -> intégration uniquement sur la zone "solide" (supportée)
      #  -> % = 100 * (∫[x_thr..x1] (AUC(x) - y_thr) dx) / (y_thr * (x1 - x_thr))
      #     > 0 : gain vs témoin ; < 0 : perte ; 0 si gains = pertes
      #
      # ATTENTION : Il faut lancer le script du plot pour certaon objet
      #
      # =========================================================
      
      # Petite intégration trapézoïdale
      trapz <- function(x, y) {
        if (length(x) < 2) return(0)
        sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
      }
      
      # On garde les segments post-seuil "solides" des 3 classes de neige
      post_solid <- curves_excl %>%
        dplyr::filter(
          side == "Après seuil (pâturé)",
          neige %in% levels(class_bounds$neige),
          lt == "solid"
        ) %>%
        dplyr::select(alpage, neige, Charge_log, AUC, AUC_lo, AUC_hi)
      
      # Table de référence témoin (y_thr) et borne x_end
      ref_tbl <- witness %>%
        dplyr::select(alpage, neige, y_thr, x_start, x_end)
      
      # Calculs d'aire + % effet
      eff_size <- post_solid %>%
        dplyr::inner_join(ref_tbl, by = c("alpage","neige")) %>%
        dplyr::group_by(alpage, neige) %>%
        dplyr::arrange(Charge_log, .by_group = TRUE) %>%
        dplyr::mutate(
          x0 = max(x_start[1],
                   suppressWarnings(min(Charge_log[Charge_log >= x_start[1]], na.rm = TRUE))),
          x1 = pmin(max(Charge_log, na.rm = TRUE), x_end[1])
        ) %>%
        dplyr::filter(Charge_log >= x0, Charge_log <= x1) %>%
        dplyr::summarise(
          x0      = dplyr::first(x0),
          x1      = dplyr::first(x1),
          y_thr   = dplyr::first(y_thr),
          n_x     = dplyr::n(),
          area_curve    = trapz(Charge_log, AUC),
          area_curve_lo = trapz(Charge_log, AUC_lo),
          area_curve_hi = trapz(Charge_log, AUC_hi),
          area_ref      = y_thr * pmax(x1 - x0, 0),
          area_diff     = area_curve    - area_ref,
          area_diff_lo  = area_curve_lo - area_ref,
          area_diff_hi  = area_curve_hi - area_ref,
          pos_area = trapz(Charge_log, pmax(AUC - y_thr, 0)),
          neg_area = trapz(Charge_log, pmax(y_thr - AUC, 0)),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          effect_pct = dplyr::if_else(area_ref > 0, 100 * area_diff    / area_ref, NA_real_),
          effect_lo  = dplyr::if_else(area_ref > 0, 100 * area_diff_lo / area_ref, NA_real_),
          effect_hi  = dplyr::if_else(area_ref > 0, 100 * area_diff_hi / area_ref, NA_real_),
          effect_pct = round(effect_pct, 1),
          effect_lo  = round(effect_lo, 1),
          effect_hi  = round(effect_hi, 1),
          effect_label    = sprintf("%+.1f%%", effect_pct),               # <- ajouté
          effect_ci_label = ifelse(is.na(effect_lo), NA_character_,
                                   sprintf("[%0.1f ; %0.1f]%%", effect_lo, effect_hi))  # <- ajouté
        ) %>%
        dplyr::arrange(alpage, neige)
      
      print(eff_size)
      
      # Optionnel : export CSV
      utils::write.csv(eff_size, "effect_size_by_alpage_neige.csv", row.names = FALSE)
      
      
      
      # =========================================================
      # 1) Effet "step" au seuil (par alpage × classe de neige)
      #    % = 100 * (AUC_post_thr - y_thr) / y_thr
      #    -> AUC_post_thr : AUC prédite à x = x_thr côté "pâturé"
      #    -> y_thr        : baseline témoin pré-seuil (déjà calculée dans `witness`)
      # =========================================================
      
      post_at_thr <- dplyr::filter(curves_excl, side == "Après seuil (pâturé)")
      
      post_thr <- post_at_thr %>%
        dplyr::group_by(alpage, neige) %>%
        dplyr::summarise(
          AUC_post_thr = stats::approx(Charge_log, AUC, xout = x_thr, rule = 2, ties = "ordered")$y,
          .groups = "drop"
        )
      
      eff_step_thr <- witness %>%
        dplyr::inner_join(post_thr, by = c("alpage","neige")) %>%
        dplyr::transmute(
          alpage, neige,
          x_thr = x_start,
          y_thr, AUC_post_thr,
          effect_bin_pct = round(100 * (AUC_post_thr - y_thr) / pmax(y_thr, .Machine$double.eps), 1)
        ) %>%
        dplyr::arrange(alpage, neige)
      
      print(eff_step_thr)
      # utils::write.csv(eff_step_thr, "effect_step_threshold_by_alpage_neige.csv", row.names = FALSE)
      
      # =========================================================
      # 2) Effet "toggle pixel-wise"
      #    Pour chaque pixel chargé observé: prédire AUC(chargé=1) et AUC(chargé=0)
      #    -> Taille d'effet pixel = 100 * (AUC1 - AUC0) / AUC0
      #    Résumé par alpage × classe de neige
      # =========================================================
      
      # Reclassement neige (mêmes bornes que `class_bounds`)
      classify_neige <- function(z) {
        dplyr::case_when(
          z <= 120 ~ "Déneigement : Précoce (≤120)",
          z <= 150 ~ "Déneigement : Moyen (121–150)",
          TRUE     ~ "Déneigement : Tardif (≥151)"
        )
      }
      
      px <- used_charged %>%            # on part des pixels "chargés" dans le domaine supporté
        dplyr::mutate(
          neige = factor(classify_neige(SMOD_2023), levels = levels(class_bounds$neige)),
          dah   = stats::median(data_mod$dah[data_mod$alpage == alpage], na.rm = TRUE) # cohérent avec tes prédictions
        )
      
      # Préparations des deux jeux de nouvelles données (seule charge_bin change)
      nd1 <- px %>% dplyr::mutate(charge_bin = 1)  # scénario chargé (observé)
      nd0 <- px %>% dplyr::mutate(charge_bin = 0)  # contre-factuel non chargé
      
      # Prédictions sur l'échelle du lien puis retour sur AUC (exp)
      pr1 <- stats::predict(m_by_alp, newdata = nd1, type = "link", exclude = exc_dah)
      pr0 <- stats::predict(m_by_alp, newdata = nd0, type = "link", exclude = exc_dah)
      
      px$AUC1 <- exp(pr1)
      px$AUC0 <- exp(pr0)
      
      px$effect_bin_pct <- 100 * (px$AUC1 - px$AUC0) / pmax(px$AUC0, .Machine$double.eps)
      
      # Résumé par alpage × classe de neige
      eff_toggle_by_pixel <- px %>%
        dplyr::group_by(alpage, neige) %>%
        dplyr::summarise(
          n         = dplyr::n(),
          mean_pct  = round(mean(effect_bin_pct, na.rm = TRUE), 1),
          median_pct= round(stats::median(effect_bin_pct, na.rm = TRUE), 1),
          q25_pct   = round(stats::quantile(effect_bin_pct, 0.25, na.rm = TRUE), 1),
          q75_pct   = round(stats::quantile(effect_bin_pct, 0.75, na.rm = TRUE), 1),
          .groups = "drop"
        ) %>%
        dplyr::arrange(alpage, neige)
      
      print(eff_toggle_by_pixel)
      # utils::write.csv(eff_toggle_by_pixel, "effect_toggle_pixel_by_alpage_neige.csv", row.names = FALSE)
      
      
      
      
      # =========================================================
      # Effet global (toutes classes & alpages confondus)
      #  -> toggle pixel-wise : % = 100 * (AUC1 - AUC0) / AUC0
      #     où AUC1 = prédiction avec charge_bin=1 ; AUC0 = contre-factuel charge_bin=0
      # =========================================================
      
      # On part des pixels "chargés" déjà filtrés dans used_charged (support en x)
      px <- used_charged %>%
        dplyr::mutate(
          # cohérent avec tes prédictions marginalisées (DAH fixé par alpage)
          dah = stats::median(data_mod$dah[data_mod$alpage == alpage], na.rm = TRUE)
        )
      
      # Deux scénarios : chargé vs non chargé (seule charge_bin change)
      nd1 <- px %>% dplyr::mutate(charge_bin = 1L)
      nd0 <- px %>% dplyr::mutate(charge_bin = 0L)
      
      # Prédictions (échelle du lien -> exponentiation)
      pr1 <- stats::predict(m_by_alp, newdata = nd1, type = "link", exclude = exc_dah)
      pr0 <- stats::predict(m_by_alp, newdata = nd0, type = "link", exclude = exc_dah)
      
      px$AUC1 <- exp(pr1)
      px$AUC0 <- exp(pr0)
      
      # Taille d'effet pixel
      px$effect_bin_pct <- 100 * (px$AUC1 - px$AUC0) / pmax(px$AUC0, .Machine$double.eps)
      
      # Résumé GLOBAL (moyenne sur tous les pixels)
      global_effect <- px %>%
        dplyr::summarise(
          n          = dplyr::n(),
          mean_pct   = round(mean(effect_bin_pct, na.rm = TRUE), 1),
          sd_pct     = stats::sd(effect_bin_pct, na.rm = TRUE),
          se_pct     = sd_pct / sqrt(n),
          ci95_lo    = round(mean(effect_bin_pct, na.rm = TRUE) - 1.96 * se_pct, 1),
          ci95_hi    = round(mean(effect_bin_pct, na.rm = TRUE) + 1.96 * se_pct, 1),
          median_pct = round(stats::median(effect_bin_pct, na.rm = TRUE), 1),
          q25_pct    = round(stats::quantile(effect_bin_pct, 0.25, na.rm = TRUE), 1),
          q75_pct    = round(stats::quantile(effect_bin_pct, 0.75, na.rm = TRUE), 1)
        )
      
      print(global_effect)
      
      # Optionnel : export
      # utils::write.csv(global_effect, "effect_global_toggle_pixel.csv", row.names = FALSE)
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
    ### TABLE PAR MODEL
    library(mgcv)
    
    # Paramètres génériques conseillés
    n_th <- max(1, parallel::detectCores(logical = TRUE) - 1)
    
    # M1 — Topographie seule (DAH) + aléatoire
    m1_bam <- bam(
      log(AUCg) ~
        s(dah,       k = 10, bs = "ts") +
        s(alpage,          bs = "re"),
      data    = data_mod,
      method  = "fREML",
      select  = TRUE,
      discrete = TRUE,
      nthreads = n_th
    )
    
    # M2 — Contexte topo-climatique (DAH + SMOD + DAH×SMOD) + aléatoire
    m2_bam <- bam(
      log(AUCg) ~
        s(dah,       k = 10, bs = "ts") +
        s(SMOD_2023, k = 10, bs = "ts") +
        ti(dah, SMOD_2023, k = c(10, 10)) +
        s(alpage,          bs = "re"),
      data    = data_mod,
      method  = "fREML",
      select  = TRUE,
      discrete = TRUE,
      nthreads = n_th
    )
    
    # M3 — + Pastoralisme (charge) et interaction avec SMOD
    m3_bam <- bam(
      log(AUCg) ~
        s(dah,       k = 10, bs = "ts") +
        s(SMOD_2023, k = 10, bs = "ts") +
        ti(dah, SMOD_2023, k = c(20, 20)) +
        s(alpage, bs = "re") +
        charge_bin +
        s(Charge_log,                 k = 10, bs = "ts", by = charge_bin) +
        ti(Charge_log, SMOD_2023, k = c(15, 15),        by = charge_bin),
      data    = data_mod,
      method  = "fREML",
      select  = TRUE,
      discrete = TRUE,
      nthreads = n_th
    )
    
    
   
    ## ====================== CLEAN JOURNAL-STYLE TABLE ======================
    
    # Packages
    pkgs <- c("mgcv","dplyr","gt","webshot2")
    to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (length(to_install)) install.packages(to_install)
    lapply(pkgs, library, character.only = TRUE)
    
    # ---------- 1) Models (bam > gam) ----------
    get_if <- function(nm) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else NULL
    m1x <- get_if("m1_bam"); if (is.null(m1x)) m1x <- get_if("m1")
    m2x <- get_if("m2_bam"); if (is.null(m2x)) m2x <- get_if("m2")
    m3x <- get_if("m3_bam"); if (is.null(m3x)) m3x <- get_if("m3")
    stopifnot(!is.null(m1x), !is.null(m2x), !is.null(m3x))
    stopifnot(exists("data_mod", inherits = TRUE), "AUCg" %in% names(data_mod))
    
    mods <- list(
      "DAH"                                   = m1x,
      "DAH × SMOD"                            = m2x,
      "DAH × SMOD + SMOD × Pastoralism"       = m3x
    )
    
    # ---------- 2) Choose performance metric ----------
    metric_choice <- "R2"  # set to "Deviance" to display deviance explained instead of R²
    
    # ---------- 3) Fast K-fold CV (RMSE on log(AUCg)) ----------
    cv_rmse <- function(m, data, K = 5, seed = 42) {
      set.seed(seed)
      id <- sample(rep(1:K, length.out = nrow(data)))
      rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
      y <- log(data$AUCg)
      mean(sapply(1:K, function(k) {
        idx <- id == k
        pred <- try(predict(m, newdata = data[idx, , drop = FALSE], type = "response"), silent = TRUE)
        if (inherits(pred, "try-error")) pred <- rep(NA_real_, sum(idx))
        rmse(pred - y[idx])
      }), na.rm = TRUE)
    }
    
    # ---------- 4) Minimal extraction ----------
    generic_of <- function(name){
      if (grepl("Pastoral", name, ignore.case = TRUE)) "Topo-climate + land use"
      else if (grepl("SMOD", name, ignore.case = TRUE)) "Topo-climatic"
      else "Topographic"
    }
    
    extract_simple <- function(m, name) {
      s <- summary(m)
      AICv <- suppressWarnings(tryCatch(AIC(m), error = function(e) NA_real_))
      dplyr::tibble(
        Generic = generic_of(name),
        Model   = name,
        n       = nobs(m),                 # pixels
        AIC     = AICv,
        `ΔAIC`  = NA_real_,
        `R2 (%)` = 100 * s$r.sq,
        `Deviance explained (%)` = 100 * s$dev.expl,
        `RMSE (CV, log)` = cv_rmse(m, data_mod, K = 5)
      )
    }
    
    tab <- dplyr::bind_rows(Map(extract_simple, mods, names(mods))) |>
      dplyr::mutate(`ΔAIC` = AIC - min(AIC, na.rm = TRUE))
    
    if (identical(tolower(metric_choice), "r2")) {
      perf_lab <- "R2 (%)"
      tab <- dplyr::select(tab, Generic, Model, n, AIC, `ΔAIC`, `R2 (%)`, `RMSE (CV, log)`)
    } else {
      perf_lab <- "Deviance explained (%)"
      tab <- dplyr::select(tab, Generic, Model, n, AIC, `ΔAIC`, `Deviance explained (%)`, `RMSE (CV, log)`)
    }
    
    # ---------- 5) GT table: no frame, no title, integers everywhere; RMSE with decimal comma ----------
    gt_tbl <- gt::gt(tab) |>
      # integers (no thousands separators)
      gt::fmt_number(columns = c("n","AIC","ΔAIC"), decimals = 0, use_seps = FALSE) |>
      gt::fmt_number(columns = perf_lab, decimals = 0, use_seps = FALSE) |>
      # RMSE: 3 decimals, decimal comma, no thousand sep
      gt::fmt_number(columns = "RMSE (CV, log)", decimals = 3, use_seps = FALSE, dec_mark = ",") |>
      gt::cols_label(`ΔAIC` = gt::html("&Delta;AIC")) |>
      # very light rules, no outer frame, no title
      gt::tab_style(
        style = gt::cell_borders(sides = "bottom", color = "black", weight = gt::px(0.6)),
        locations = gt::cells_column_labels()
      ) |>
      gt::tab_style(
        style = gt::cell_borders(sides = "bottom", color = "#cccccc", weight = gt::px(0.4)),
        locations = gt::cells_body()
      ) |>
      gt::tab_options(
        table.font.names = c("Times New Roman","Times","Serif"),
        data_row.padding = gt::px(3),
        table.border.top.color = "white",
        table.border.bottom.color = "white"
      ) 
      
    
    # --- Rendre les en-têtes de colonnes en gras ---
    gt_tbl <- gt_tbl |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_column_labels()
      )
    # (optionnel, certaines versions) :
    gt_tbl <- gt_tbl |>
      gt::tab_options(column_labels.font.weight = "bold")
    
    # Export à nouveau
    gt::gtsave(gt_tbl, "gam_model_fit_clean.html", inline_css = TRUE)
    gt::gtsave(gt_tbl, "gam_model_fit_clean.png", vwidth = 1500, vheight = 420)
    
    
    gt_tbl
    
    
    
    # --- Point webshot2 vers un navigateur Chromium (Chrome/Edge/Brave) ---
    set_chromium <- function() {
      cand <- c(
        file.path(Sys.getenv("LOCALAPPDATA"), "Google/Chrome/Application/chrome.exe"),
        "C:/Program Files/Google/Chrome/Application/chrome.exe",
        "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe",
        "C:/Program Files/Microsoft/Edge/Application/msedge.exe",
        "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe",
        "C:/Program Files/BraveSoftware/Brave-Browser/Application/brave.exe"
      )
      cand <- cand[file.exists(cand)]
      if (!length(cand)) stop("Chromium browser not found; set CHROMOTE_CHROME manually.")
      Sys.setenv(CHROMOTE_CHROME = cand[1])
      message("Using Chromium at: ", cand[1])
    }
    set_chromium()
    
    # Export
    gt::gtsave(gt_tbl, "gam_model_fit_clean.html", inline_css = TRUE)
    gt::gtsave(gt_tbl, "gam_model_fit_clean.png", vwidth = 1500, vheight = 420)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # ==================== PSEUDO Δ%DEV PAR GROUPES, PAR ALPAGE + TOTAL ====================
    # Données attendues : data_mod avec colonnes AUCg, dah, SMOD_2023, alpage, charge_bin, Charge_log
    
    # ---- Packages ----
    req_pkgs <- c("mgcv","dplyr","gt")
    to_install <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
    if (length(to_install)) install.packages(to_install)
    invisible(lapply(req_pkgs, library, character.only = TRUE))
    
    # ---- Paramètres ----
    set.seed(42)
    n_th <- max(1, parallel::detectCores(logical = TRUE) - 1)
    clip_neg_to_zero <- FALSE   # <- mets TRUE si tu veux interdire les pourcentages négatifs à l'affichage
    
    # ---- Données ----
    stopifnot(exists("data_mod"), is.data.frame(data_mod))
    needed <- c("AUCg","dah","SMOD_2023","alpage","charge_bin","Charge_log")
    miss <- setdiff(needed, names(data_mod)); if (length(miss)) stop("Colonnes manquantes: ", paste(miss, collapse=", "))
    
    dat <- data_mod |>
      dplyr::select(all_of(needed)) |>
      dplyr::filter(is.finite(AUCg), AUCg > 0,
                    is.finite(dah), is.finite(SMOD_2023),
                    !is.na(alpage), is.finite(Charge_log)) |>
      dplyr::mutate(alpage = factor(alpage),
                    charge_bin = factor(charge_bin))
    if (nrow(dat) < 10) stop("Trop peu d'observations après filtrage: ", nrow(dat))
    
    # ---- Modèles imbriqués M0→M4 ----
    fam <- gaussian()
    M0 <- mgcv::bam(log(AUCg) ~ s(alpage, bs="re"),
                    data=dat, family=fam, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th)
    M1 <- mgcv::bam(log(AUCg) ~ s(alpage, bs="re") + s(dah, k=10, bs="ts"),
                    data=dat, family=fam, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th)
    M2 <- mgcv::bam(log(AUCg) ~ s(alpage, bs="re") + s(dah, k=10, bs="ts") +
                      s(SMOD_2023, k=10, bs="ts"),
                    data=dat, family=fam, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th)
    M3 <- mgcv::bam(log(AUCg) ~ s(alpage, bs="re") + s(dah, k=10, bs="ts") +
                      s(SMOD_2023, k=10, bs="ts") +
                      ti(dah, SMOD_2023, k=c(10,10)),
                    data=dat, family=fam, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th)
    M4 <- mgcv::bam(log(AUCg) ~ s(alpage, bs="re") + s(dah, k=10, bs="ts") +
                      s(SMOD_2023, k=10, bs="ts") +
                      ti(dah, SMOD_2023, k=c(10,10)) +
                      charge_bin +
                      s(Charge_log, k=10, bs="ts", by=charge_bin) +
                      ti(Charge_log, SMOD_2023, k=c(10,10), by=charge_bin),
                    data=dat, family=fam, method="fREML", select=TRUE, discrete=TRUE, nthreads=n_th)
    
    # ---- Prédictions ----
    y_log <- log(dat$AUCg)
    preds <- list(
      M0 = as.numeric(predict(M0, newdata=dat, type="response")),
      M1 = as.numeric(predict(M1, newdata=dat, type="response")),
      M2 = as.numeric(predict(M2, newdata=dat, type="response")),
      M3 = as.numeric(predict(M3, newdata=dat, type="response")),
      M4 = as.numeric(predict(M4, newdata=dat, type="response"))
    )
    
    # ---- Pseudo Δ%Dev séquentiel par alpage ----
    pseudo_dev_by_group <- function(y, yhat_list, group, clip_neg=FALSE){
      g <- factor(group); lev <- levels(g)
      out <- lapply(lev, function(lv){
        idx <- (g == lv); if (!any(idx)) return(NULL)
        yy <- y[idx]
        sse <- sapply(yhat_list, function(ph) sum((yy - ph[idx])^2))
        sse0 <- sse["M0"]
        if (!is.finite(sse0) || sse0 <= .Machine$double.eps) {
          topo <- clim <- inter <- pasto <- tot_dir <- NA_real_
        } else {
          topo   <- 100 * (sse["M0"] - sse["M1"]) / sse0
          clim   <- 100 * (sse["M1"] - sse["M2"]) / sse0
          inter  <- 100 * (sse["M2"] - sse["M3"]) / sse0
          pasto  <- 100 * (sse["M3"] - sse["M4"]) / sse0
          tot_dir<- 100 * (sse["M0"] - sse["M4"]) / sse0  # Total (M4 vs M0), indépendant des arrondis
        }
        if (clip_neg) {
          topo <- pmax(topo, 0); clim <- pmax(clim, 0); inter <- pmax(inter, 0); pasto <- pmax(pasto, 0)
        }
        tot_sum <- topo + clim + inter + pasto
        dplyr::tibble(
          Alpage = lv,
          `Topo (DAH)` = topo,
          `Climat (SMOD)` = clim,
          `Interaction DAH×SMOD` = inter,
          `Pastoralisme` = pasto,
          `Total (somme)` = tot_sum,
          `Total (M4 vs M0)` = tot_dir,
          n = sum(idx)
        )
      })
      dplyr::bind_rows(out)
    }
    
    tab_alpage <- pseudo_dev_by_group(y_log, preds, dat$alpage, clip_neg=clip_neg_to_zero)
    
    # ---- Ligne globale pondérée ----
    tab_global <- tab_alpage |>
      dplyr::summarise(
        Alpage = "Global (pondéré)",
        `Topo (DAH)` = stats::weighted.mean(`Topo (DAH)`, w=n, na.rm=TRUE),
        `Climat (SMOD)` = stats::weighted.mean(`Climat (SMOD)`, w=n, na.rm=TRUE),
        `Interaction DAH×SMOD` = stats::weighted.mean(`Interaction DAH×SMOD`, w=n, na.rm=TRUE),
        `Pastoralisme` = stats::weighted.mean(`Pastoralisme`, w=n, na.rm=TRUE),
        `Total (somme)` = stats::weighted.mean(`Total (somme)`, w=n, na.rm=TRUE),
        `Total (M4 vs M0)` = stats::weighted.mean(`Total (M4 vs M0)`, w=n, na.rm=TRUE),
        n = sum(n)
      )
    
    tab_final <- dplyr::bind_rows(tab_alpage, tab_global)
    
    # ---- Table GT ----
    gt_tbl <- gt::gt(tab_final) |>
      gt::fmt_number(columns = c("Topo (DAH)","Climat (SMOD)","Interaction DAH×SMOD",
                                 "Pastoralisme","Total (somme)","Total (M4 vs M0)"),
                     decimals = 1, use_seps = FALSE, dec_mark = ",") |>
      gt::fmt_number(columns = "n", decimals = 0, use_seps = FALSE) |>
      gt::tab_style(style = gt::cell_borders(sides="bottom", color="black", weight=gt::px(0.6)),
                    locations = gt::cells_column_labels()) |>
      gt::tab_style(style = gt::cell_borders(sides="bottom", color="#cccccc", weight=gt::px(0.4)),
                    locations = gt::cells_body()) |>
      gt::tab_style(style = gt::cell_text(weight="bold"),
                    locations = gt::cells_column_labels()) |>
      gt::tab_options(table.font.names = c("Times New Roman","Times","Serif"),
                      data_row.padding = gt::px(3),
                      table.border.top.color = "white",
                      table.border.bottom.color = "white")
    gt_tbl
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ## SECTION HABITAT TAILLE EFFET FOCUS SUR 2023 : 
    # =========================================================
    # Packages
    # =========================================================
    library(mgcv)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(stringr)
    
    # =========================================================
    # Données d'entrée : data_mod (2023)
    # On suppose que data_mod correspond à 2023 (pas de colonne "YEAR")
    # =========================================================
    # Sanity check des types + colonnes attendues
    stopifnot(all(c("alpage","AUCg","dah","SMOD_2023",
                    "Charge","Charge_log","charge_bin",
                    "habitat_code","habitat_label") %in% names(data_mod)))
    
    data_2023 <- data_mod %>%
      mutate(
        alpage     = factor(alpage),
        # On garde charge_bin en 0/1 NUMÉRIQUE (utile pour "gater" les smooths post-seuil)
        charge_bin = as.integer(charge_bin),
        # Petite sécurité si NA dans charge_bin : on le reconstruit
        charge_bin = ifelse(is.na(charge_bin), as.integer(Charge_log >= log(2)), charge_bin)
      )
    
    # =========================================================
    # Carte des 4 habitats focus
    # =========================================================
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    focus_codes <- c(1, 9, 5, 6)
    focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
    
    # Pour filtrages
    data_2023 <- data_2023 %>%
      mutate(
        habitat_code = as.integer(habitat_code),
        hab_focus = if_else(habitat_code %in% focus_codes,
                            hab_labels_map[as.character(habitat_code)], NA_character_)
      )
    
    # =========================================================
    # Formule modèle (approche “référence” avec courbes par alpage)
    # - base topo + neige
    # - aléatoire d'intercept par alpage
    # - saut paramétrique au seuil (charge_bin)
    # - lissage de Charge_log DIFFÉRENT par alpage (bs="fs"), ACTIVÉ seulement côté post (by=charge_bin)
    # - interaction Charge × neige activée côté post
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
    # Choix du moteur : bam (plus rapide) + threads
    # =========================================================
    n_th <- tryCatch({
      max(1L, parallel::detectCores(logical = TRUE) - 1L)
    }, error = function(e) 1L)
    
    fit_fun <- function(dat) {
      # IMPORTANT : alpage en facteur (subset) et charge_bin NUMÉRIQUE 0/1
      dat <- dat %>%
        mutate(
          alpage     = factor(alpage),              # refactorise sur le subset
          charge_bin = as.integer(Charge_log >= log(2))  # recalcul propre (sécurité)
        )
      bam(
        form_by_alp,
        data     = dat,
        method   = "fREML",
        select   = TRUE,
        discrete = TRUE,
        nthreads = n_th
      )
    }
    
    # =========================================================
    # Fit modèles
    # - Global 2023 (tous habitats)
    # - 4 habitats focus (modèles séparés)
    # =========================================================
    mod_2023_global <- fit_fun(data_2023)
    
    mods_hab_2023 <- map(setNames(as.list(focus_codes), as.character(focus_codes)), function(hc) {
      di <- data_2023 %>% filter(habitat_code == hc)
      if (nrow(di) < 50) return(NULL)  # garde-fou minimal
      fit_fun(di)
    })
    
    # =========================================================
    # Fonction magnitude d'effet Δ% (Q25 -> Q75 de Charge_log)
    # - diff sur l'échelle log : eff = dbar %*% beta
    # - SE via V (delta-method linéaire)
    # - % = 100*(exp(eff)-1)
    # - Crucial : RECOMPUTER charge_bin pour nd_lo / nd_hi
    # =========================================================
    delta_q25_q75 <- function(model, data_in, group_label = "Total") {
      d <- data_in
      # garde uniquement les valeurs > 0 (comme ton brouillon)
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
      
      nd_lo <- d; nd_lo$Charge_log <- ql
      nd_hi <- d; nd_hi$Charge_log <- qh
      # >>> cohérence des termes by=charge_bin :
      nd_lo$charge_bin <- as.integer(nd_lo$Charge_log >= log(2))
      nd_hi$charge_bin <- as.integer(nd_hi$Charge_log >= log(2))
      
      # lpmatrix diff + moyenne sur les lignes (pondération uniforme sur la distribution observée)
      Xhi  <- predict(model, newdata = nd_hi, type = "lpmatrix")
      Xlo  <- predict(model, newdata = nd_lo, type = "lpmatrix")
      dbar <- colMeans(Xhi - Xlo)
      
      b <- coef(model); V <- vcov(model)
      eff <- as.numeric(drop(dbar %*% b))
      se  <- sqrt(as.numeric(drop(dbar %*% V %*% dbar)))
      lwr <- eff - 1.96 * se
      upr <- eff + 1.96 * se
      
      tibble::tibble(
        group = group_label,
        effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
        effect_perc = 100 * (exp(eff) - 1),
        lwr_perc    = 100 * (exp(lwr) - 1),
        upr_perc    = 100 * (exp(upr) - 1),
        sig         = (lwr > 0) | (upr < 0),
        q_low = ql, q_high = qh, n = nrow(d)
      )
    }
    
    # =========================================================
    # Extraction des effets (2023)
    # =========================================================
    eff_total <- delta_q25_q75(mod_2023_global, data_2023, group_label = "Total")
    
    eff_habs <- imap_dfr(mods_hab_2023, function(mod, code_chr) {
      lab <- hab_labels_map[[code_chr]]
      di  <- data_2023 %>% filter(habitat_code == as.integer(code_chr))
      if (is.null(mod)) {
        tibble::tibble(
          group = lab,
          effect_log = NA_real_, se_log = NA_real_, lwr_log = NA_real_, upr_log = NA_real_,
          effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
          sig = FALSE, q_low = NA_real_, q_high = NA_real_, n = nrow(di)
        )
      } else {
        delta_q25_q75(mod, di, group_label = lab)
      }
    })
    
    eff_2023 <- bind_rows(eff_total, eff_habs) %>%
      mutate(
        group = factor(group, levels = c("Total", focus_labels))
      )
    
    # =========================================================
    # Barplot
    # =========================================================
    cols_focus <- c(
      "Total"                        = "grey60",
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    pal <- cols_focus[levels(eff_2023$group)]
    
    p_mag_2023 <- ggplot(eff_2023, aes(x = group, y = effect_perc, fill = group)) +
      geom_col(width = 0.80, aes(alpha = ifelse(sig, 1, 0.4))) +
      geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc), width = 0.20) +
      scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
      scale_alpha_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(
        title = "Effet de la charge sur la production — 2023",
        subtitle = expression(Delta~"% entre "~Q[25]~" et "~Q[75]~" de "~Charge[log]),
        x = NULL, y = "Variation (%) vs Q25"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(size = 12))
    
    print(p_mag_2023)
    
    # (Optionnel) Sauvegarde
    # ggsave("barplot_magnitude_charge_2023.png", p_mag_2023, width = 9, height = 5, dpi = 300)
    
    
    
    
    
    
    
    
    
    # =========================================================
    # Packages
    # =========================================================
    library(mgcv)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(stringr)
    
    # =========================================================
    # Données multi-annuelles : dt_legacy
    #  → on suppose les colonnes suivantes (comme pour m07) :
    #     alpage, GPROD, DAH, SMOD, Charge_log, x, y, year_f, habitat_code
    # =========================================================
    stopifnot(all(c("alpage","GPROD","DAH","SMOD","Charge_log","x","y","year_f","habitat_code") %in% names(dt_legacy)))
    
    # 4 habitats ciblés
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    focus_codes   <- as.integer(names(hab_labels_map))
    focus_labels  <- unname(hab_labels_map)
    
    # Sous-ensemble multi-années
    dt_hab <- dt_legacy %>%
      filter(habitat_code %in% focus_codes) %>%
      mutate(
        habitat_code = as.integer(habitat_code),
        habitat4 = factor(habitat_code, levels = focus_codes, labels = focus_labels),
        alpage   = factor(alpage),
        year_f   = factor(year_f)
      ) %>%
      filter(is.finite(GPROD), is.finite(DAH), is.finite(SMOD),
             is.finite(Charge_log), is.finite(x), is.finite(y)) %>%
      droplevels()
    
    # =========================================================
    # Modèle multi-annuel + habitats (simple comme m07_hab)
    # log(GPROD) ~ base topo + neige + lissage Charge par année + par habitat
    # =========================================================
    
    # k "sûrs" (petits) pour éviter les soucis d’unicité
    min_u_year <- dt_hab %>% group_by(year_f)   %>% summarise(nu = n_distinct(Charge_log), .groups="drop") %>% pull(nu) %>% min(na.rm=TRUE)
    min_u_hab  <- dt_hab %>% group_by(habitat4) %>% summarise(nu = n_distinct(Charge_log), .groups="drop") %>% pull(nu) %>% min(na.rm=TRUE)
    n_xy       <- dt_hab %>% distinct(x,y) %>% nrow()
    
    k_year <- max(4, min(8, min_u_year - 1))
    k_hab  <- max(4, min(8, min_u_hab  - 1))
    k_xy   <- max(10, min(30, n_xy - 1))
    
    m_multi <- bam(
      log(GPROD) ~
        s(DAH,  k = 10, bs = "ts") +
        s(SMOD, k = 10, bs = "ts") +
        ti(DAH, SMOD, k = c(20,20), bs = c("ts","ts")) +
        s(Charge_log, by = year_f,   bs = "ts", k = k_year) +  # effet Charge par année
        s(Charge_log, by = habitat4, bs = "ts", k = k_hab ) +  # effet Charge par habitat
        s(alpage, bs = "re") +
        s(x, y, bs = "tp", k = k_xy) +
        year_f + habitat4,
      data     = dt_hab,
      method   = "fREML",
      select   = TRUE,
      discrete = TRUE
    )
    
    # =========================================================
    # Magnitude d’effet Δ% (Q25 -> Q75) — même logique que ton script
    # (on moyenne sur la distribution observée du sous-ensemble fourni)
    # =========================================================
    delta_q0_q50 <- function(model, data_in, group_label = "Total") {
      d <- data_in
      
      # Comme avant : on ne garde que les valeurs > 0 (cohérent avec ton brouillon).
      # ➜ Si tu veux inclure les zéros ou valeurs ≤0, remplace la ligne suivante par:
      #    pos <- d$Charge_log[is.finite(d$Charge_log)]
      pos <- d$Charge_log[d$Charge_log > 0 & is.finite(d$Charge_log)]
      
      if (length(pos) < 10) {
        return(tibble::tibble(
          group = group_label,
          effect_log = NA_real_, se_log = NA_real_,
          lwr_log = NA_real_, upr_log = NA_real_,
          effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
          sig = FALSE, q_low = NA_real_, q_high = NA_real_, n = nrow(d)
        ))
      }
      
      # Q0 = minimum, Q50 = médiane
      ql <- as.numeric(quantile(pos, 0.00, na.rm = TRUE))  # min
      qh <- as.numeric(quantile(pos, 0.50, na.rm = TRUE))  # médiane
      
      # Fallback de sécurité si degénéré
      if (!is.finite(ql) || !is.finite(qh) || ql >= qh) {
        ql <- as.numeric(quantile(pos, 0.00, na.rm = TRUE))
        qh <- as.numeric(quantile(pos, 0.50, na.rm = TRUE))
      }
      
      nd_lo <- d; nd_lo$Charge_log <- ql
      nd_hi <- d; nd_hi$Charge_log <- qh
      
      Xhi  <- predict(model, newdata = nd_hi, type = "lpmatrix")
      Xlo  <- predict(model, newdata = nd_lo, type = "lpmatrix")
      dbar <- colMeans(Xhi - Xlo)
      
      b <- coef(model); V <- vcov(model)
      eff <- as.numeric(drop(dbar %*% b))
      se  <- sqrt(as.numeric(drop(dbar %*% V %*% dbar)))
      lwr <- eff - 1.96 * se
      upr <- eff + 1.96 * se
      
      tibble::tibble(
        group = group_label,
        effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
        effect_perc = 100 * (exp(eff) - 1),
        lwr_perc    = 100 * (exp(lwr) - 1),
        upr_perc    = 100 * (exp(upr) - 1),
        sig         = (lwr > 0) | (upr < 0),
        q_low = ql, q_high = qh, n = nrow(d)
      )
    }
    
    # =========================================================
    # Extraction des effets — pour CHAQUE année :
    #   - "Total" (tous habitats de l’année)
    #   - 4 habitats (dans l’année)
    # → 5 barres par année, 7 années = 7 facettes
    # =========================================================
    lev_year <- levels(dt_hab$year_f)
    
    eff_all <- map_dfr(lev_year, function(yy) {
      dy <- dt_hab %>% filter(year_f == yy)
      
      # Total
      e_total <- delta_q25_q75(m_multi, dy, group_label = "Total")
      
      # Par habitat de l'année
      e_habs <- map_dfr(focus_codes, function(hc) {
        lab <- hab_labels_map[[as.character(hc)]]
        di  <- dy %>% filter(habitat_code == hc)
        delta_q25_q75(m_multi, di, group_label = lab)
      })
      
      bind_rows(
        e_total %>% mutate(year_f = yy),
        e_habs  %>% mutate(year_f = yy)
      )
    }) %>%
      mutate(
        year_f = factor(year_f, levels = lev_year),
        group  = factor(group, levels = c("Total", focus_labels))
      )
    
    # =========================================================
    # Barplot facetté par année (7 panneaux, 5 barres chacun)
    # =========================================================
    cols_focus <- c(
      "Total"                        = "grey60",
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    pal <- cols_focus[levels(eff_all$group)]
    
    p_mag_multi <- ggplot(eff_all, aes(x = group, y = effect_perc, fill = group)) +
      geom_col(width = 0.80, aes(alpha = ifelse(sig, 1, 0.4))) +
      geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc), width = 0.20) +
      scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
      scale_alpha_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      facet_wrap(~ year_f, ncol = 3, scales = "free_y") +
      labs(
        title = "Effet de la charge sur la production — 7 années",
        subtitle = expression(Delta~"% entre "~Q[25]~" et "~Q[75]~" de "~Charge[log]~" (par année)"),
        x = NULL, y = "Variation (%) vs Q25"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(size = 11, angle = 20, hjust = 1))
    
    print(p_mag_multi)
    
    # (Option) Sauvegarde
    # ggsave("barplot_magnitude_charge_multi_annees.png", p_mag_multi, width = 11, height = 8, dpi = 300)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # =========================================================
    # Packages
    # =========================================================
    library(mgcv)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(stringr)
    
    # =========================================================
    # Données d'entrée : data_mod (2023)
    # =========================================================
    stopifnot(all(c("alpage","AUCg","dah","SMOD_2023",
                    "Charge","Charge_log","charge_bin",
                    "habitat_code","habitat_label") %in% names(data_mod)))
    
    data_2023 <- data_mod %>%
      mutate(
        alpage     = factor(alpage),
        # Sécurise charge_bin en 0/1 (recalcule si NA ou incohérent)
        charge_bin = as.integer(ifelse(is.na(charge_bin), Charge_log >= log(2), charge_bin)),
        habitat_code  = as.integer(habitat_code)
      )
    
    # =========================================================
    # 4 habitats focus
    # =========================================================
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    focus_codes  <- c(1, 9, 5, 6)
    focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
    
    data_2023 <- data_2023 %>%
      mutate(hab_focus = if_else(habitat_code %in% focus_codes,
                                 hab_labels_map[as.character(habitat_code)],
                                 NA_character_))
    
    # =========================================================
    # Modèle BINAIRE (charge seulement)
    # - base topo + neige
    # - aléatoire d'intercept par alpage
    # - effet paramétrique charge_bin (0→1)
    # - random slope charge_bin par alpage (variation d’effet entre alpages)
    # =========================================================
    form_bin <- formula(
      log(AUCg) ~
        s(dah,       k = 10, bs = "ts") +
        s(SMOD_2023, k = 10, bs = "ts") +
        ti(dah, SMOD_2023, k = c(20, 20)) +
        s(alpage, bs = "re") +                 # intercept alpage
        charge_bin +                           # effet moyen 0→1
        s(alpage, bs = "re", by = charge_bin)  # random slope 0→1 par alpage
    )
    
    # =========================================================
    # Fit : bam (rapide) + threads
    # =========================================================
    n_th <- tryCatch(max(1L, parallel::detectCores(logical = TRUE) - 1L), error = function(e) 1L)
    
    fit_fun <- function(dat) {
      dat <- dat %>%
        mutate(
          alpage     = factor(alpage),
          charge_bin = as.integer(Charge_log >= log(2)) # (cohérence)
        )
      # Il faut des 0 ET des 1, sinon on ne peut pas estimer l'effet 0→1
      if (length(unique(dat$charge_bin)) < 2L) return(NULL)
      
      bam(form_bin, data = dat, method = "fREML",
          select = TRUE, discrete = TRUE, nthreads = n_th)
    }
    
    # =========================================================
    # Fit modèles (global + 4 habitats)
    # =========================================================
    mod_2023_global <- fit_fun(data_2023)
    
    mods_hab_2023 <- map(setNames(as.list(focus_codes), as.character(focus_codes)), function(hc) {
      di <- data_2023 %>% filter(habitat_code == hc)
      if (nrow(di) < 50) return(NULL)  # garde-fou minimal
      fit_fun(di)
    })
    
    # =========================================================
    # Effet "0 → 1" (chargé vs non chargé) en % + IC95%
    # - On clone le DF en charge_bin=0 et charge_bin=1
    # - Δ = moyenne_lignes( X1 - X0 ) %*% beta
    # - % = 100 * (exp(Δ) - 1)
    # - Garde-fous : min 30 obs par état & ≥5% de données dans chaque état
    # =========================================================
    effect_bin01 <- function(model, data_in, group_label = "Total",
                             nmin_per_state = 30, prop_min_state = 0.05) {
      if (is.null(model)) {
        return(tibble(
          group = group_label, method = "0→1 (insuffisant)",
          effect_log = NA_real_, se_log = NA_real_, lwr_log = NA_real_, upr_log = NA_real_,
          effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
          sig = FALSE, n = nrow(data_in), n0 = NA_integer_, n1 = NA_integer_
        ))
      }
      d <- data_in %>% mutate(charge_bin = as.integer(Charge_log >= log(2)))
      n0 <- sum(d$charge_bin == 0, na.rm = TRUE)
      n1 <- sum(d$charge_bin == 1, na.rm = TRUE)
      
      if (n0 < nmin_per_state || n1 < nmin_per_state ||
          n0/nrow(d) < prop_min_state || n1/nrow(d) < prop_min_state) {
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
      lwr <- eff - 1.96 * se; upr <- eff + 1.96 * se
      
      tibble(
        group = group_label, method = "0→1",
        effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
        effect_perc = 100 * (exp(eff) - 1),
        lwr_perc    = 100 * (exp(lwr) - 1),
        upr_perc    = 100 * (exp(upr) - 1),
        sig = (lwr > 0) | (upr < 0),
        n = nrow(d), n0 = n0, n1 = n1
      )
    }
    
    # =========================================================
    # Extraction des effets (2023)
    # =========================================================
    eff_total <- effect_bin01(mod_2023_global, data_2023, "Total")
    
    eff_habs <- imap_dfr(mods_hab_2023, function(mod, code_chr) {
      lab <- hab_labels_map[[code_chr]]
      di  <- data_2023 %>% filter(habitat_code == as.integer(code_chr))
      effect_bin01(mod, di, lab)
    })
    
    eff_2023 <- bind_rows(eff_total, eff_habs) %>%
      mutate(group = factor(group, levels = c("Total", focus_labels)),
             alpha_plot = ifelse(sig, 1, 0.4),
             lty_plot   = ifelse(grepl("insuffisant", method, ignore.case = TRUE), "twodash", "solid"))
    
    # =========================================================
    # Barplot
    # =========================================================
    cols_focus <- c(
      "Total"                        = "grey60",
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    pal <- cols_focus[levels(eff_2023$group)]
    
    p_mag_2023 <- ggplot(eff_2023, aes(x = group, y = effect_perc, fill = group)) +
      geom_col(width = 0.80, aes(alpha = alpha_plot)) +
      geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc, linetype = lty_plot), width = 0.20) +
      scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
      scale_alpha_identity() +
      scale_linetype_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(
        title = "Effet binaire de la charge (0 → 1) — 2023",
        subtitle = "Δ% entre états non chargé et chargé (IC95%). 'twodash' = effect non estimable (données insuffisantes).",
        x = NULL, y = "Variation (%)"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(size = 12))
    
    print(p_mag_2023)
    
    
    
    
    
    
    
    
    
    
    # ============================
    # Summaries des modèles 2023
    # ============================
    
    stopifnot(exists("mod_2023_global"),
              exists("mods_hab_2023"),
              exists("hab_labels_map"))
    
    cat("\n================  GLOBAL 2023  ================\n")
    if (is.null(mod_2023_global)) {
      cat("Modèle global introuvable / non estimé.\n")
    } else {
      print(summary(mod_2023_global))
    }
    
    # Habitats focus
    purrr::iwalk(mods_hab_2023, function(mod, code_chr){
      lab <- hab_labels_map[[code_chr]]
      cat("\n================  HABITAT:", lab, "(code", code_chr, ")  ================\n")
      if (is.null(mod)) {
        cat("Modèle non estimé (données insuffisantes : pas assez de 0 ET 1).\n")
      } else {
        print(summary(mod))
      }
    })
    # (Optionnel) sauvegarde
    # ggsave("barplot_charge_binaire_2023.png", p_mag_2023, width = 9, height = 5, dpi = 300)
    
    
    
    
    
    
    
    
    coef_to_pct <- function(mod){
      if (is.null(mod)) return(c(NA,NA,NA))
      b <- coef(mod)["charge_bin"]
      se <- sqrt(vcov(mod)["charge_bin","charge_bin"])
      lwr <- b - 1.96*se; upr <- b + 1.96*se
      c(100*(exp(b)-1), 100*(exp(lwr)-1), 100*(exp(upr)-1))
    }
    
    tab_coef <- tibble::tibble(
      group = c("Total",
                "P. nivales","P. thermiques écorchées",
                "Queyrellins","Nardaies denses du subalpin"),
      pct = c(coef_to_pct(mod_2023_global)[1],
              coef_to_pct(mods_hab_2023[["1"]])[1],
              coef_to_pct(mods_hab_2023[["9"]])[1],
              coef_to_pct(mods_hab_2023[["5"]])[1],
              coef_to_pct(mods_hab_2023[["6"]])[1]),
      lwr = c(coef_to_pct(mod_2023_global)[2],
              coef_to_pct(mods_hab_2023[["1"]])[2],
              coef_to_pct(mods_hab_2023[["9"]])[2],
              coef_to_pct(mods_hab_2023[["5"]])[2],
              coef_to_pct(mods_hab_2023[["6"]])[2]),
      upr = c(coef_to_pct(mod_2023_global)[3],
              coef_to_pct(mods_hab_2023[["1"]])[3],
              coef_to_pct(mods_hab_2023[["9"]])[3],
              coef_to_pct(mods_hab_2023[["5"]])[3],
              coef_to_pct(mods_hab_2023[["6"]])[3])
    )
    
    ggplot(tab_coef, aes(group, pct, fill = group)) +
      geom_col(width=0.8) +
      geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.2) +
      labs(title="Effet conditionnel (coefficient charge_bin)",
           y="Variation (%)", x=NULL) +
      theme_minimal(base_size=14) + theme(legend.position="none")
    
    
    
    
    
    
    ### BIN EFFET ALPAGE EN ALEATOIRE 
    # =========================================================
    # Packages
    # =========================================================
    library(mgcv)
    library(dplyr)
    library(ggplot2)
    library(purrr)
    
    # =========================================================
    # Données 2023
    # =========================================================
    stopifnot(all(c("alpage","AUCg","dah","SMOD_2023",
                    "Charge_log","charge_bin","habitat_code","habitat_label") %in% names(data_mod)))
    
    data_2023 <- data_mod %>%
      mutate(
        alpage        = factor(alpage),
        habitat_code  = as.integer(habitat_code),
        charge_bin    = as.integer(ifelse(is.na(charge_bin), Charge_log >= log(2), charge_bin))
      )
    
    # Habitats focus
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    focus_codes  <- c(1,9,5,6)
    focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
    
    # =========================================================
    # Modèle binaire SANS random slope (effet 0→1 unique)
    # =========================================================
    form_bin_simple <- formula(
      log(AUCg) ~
        s(dah,       k = 10, bs = "ts") +
        s(SMOD_2023, k = 10, bs = "ts") +
        ti(dah, SMOD_2023, k = c(20,20)) +
        s(alpage, bs = "re") +   # intercept aléatoire (différences de niveau par alpage)
        charge_bin               # effet global 0→1 UNIQUE
    )
    
    n_th <- tryCatch(max(1L, parallel::detectCores(TRUE)-1L), error=function(e) 1L)
    fit_fun <- function(dat){
      dat <- dat %>% mutate(
        alpage     = factor(alpage),
        charge_bin = as.integer(Charge_log >= log(2))
      )
      # besoin des deux états 0 et 1
      if (length(unique(dat$charge_bin)) < 2L) return(NULL)
      bam(form_bin_simple, data = dat, method="fREML", select=TRUE,
          discrete=TRUE, nthreads=n_th)
    }
    
    # Fit global + 4 habitats (modèles séparés, même formule)
    mod_2023_global_simple <- fit_fun(data_2023)
    
    mods_hab_2023_simple <- map(setNames(as.list(focus_codes), as.character(focus_codes)), function(hc){
      di <- data_2023 %>% filter(habitat_code == hc)
      if (nrow(di) < 50) return(NULL)
      fit_fun(di)
    })
    
    # =========================================================
    # Effet 0→1 directement depuis le coefficient (conditionnel = marginal ici)
    # =========================================================
    effect_from_coef <- function(model, label){
      if (is.null(model)) return(tibble(
        group = label, effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
        sig = FALSE, method = "coef", n = NA_integer_
      ))
      b  <- coef(model)["charge_bin"]
      se <- sqrt(vcov(model)["charge_bin","charge_bin"])
      lwr <- b - 1.96*se; upr <- b + 1.96*se
      tibble(
        group = label, method = "coef",
        effect_perc = 100*(exp(b)-1),
        lwr_perc    = 100*(exp(lwr)-1),
        upr_perc    = 100*(exp(upr)-1),
        sig = (lwr > 0) | (upr < 0),
        n = length(model$y)
      )
    }
    
    eff_total <- effect_from_coef(mod_2023_global_simple, "Total")
    eff_habs  <- imap_dfr(mods_hab_2023_simple, ~ effect_from_coef(.x, hab_labels_map[[.y]]))
    
    eff_2023 <- bind_rows(eff_total, eff_habs) %>%
      mutate(group = factor(group, levels = c("Total", focus_labels)),
             alpha_plot = ifelse(sig, 1, 0.4))
    
    # =========================================================
    # Barplot 2023
    # =========================================================
    cols_focus <- c(
      "Total"                        = "grey60",
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    pal <- cols_focus[levels(eff_2023$group)]
    
    p_bin_simple_2023 <- ggplot(eff_2023, aes(x = group, y = effect_perc, fill = group)) +
      geom_col(width = 0.80, aes(alpha = alpha_plot)) +
      geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc), width = 0.20) +
      scale_fill_manual(values = pal, guide = "none") +
      scale_alpha_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(
        title = "Effet binaire global de la charge (0 → 1) — 2023",
        subtitle = "Effet unique (alpage = intercept aléatoire). Barres = IC95% issus du coef. charge_bin.",
        x = NULL, y = "Variation (%)"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 12))
    
    print(p_bin_simple_2023)
    
    # Summaries (optionnel)
    cat("\n=== GLOBAL simple ===\n"); if (!is.null(mod_2023_global_simple)) print(summary(mod_2023_global_simple))
    purrr::iwalk(mods_hab_2023_simple, function(m, k){
      cat("\n=== HABITAT simple:", hab_labels_map[[k]], "===\n")
      if (is.null(m)) cat("Modèle non estimé (pas les deux états 0 et 1).\n") else print(summary(m))
    })
    
    # # Sauvegarde si besoin :
    # ggsave("barplot_charge_binaire_global_2023.png", p_bin_simple_2023, width = 9, height = 5, dpi = 300)
    
    
    
    
    
    ### SANS CHARGE BIN AVEC QUE DU CONTINUE : # =========================================================
    # Packages
    # =========================================================
    library(mgcv)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    
    # =========================================================
    # Données 2023
    # =========================================================
    stopifnot(all(c("alpage","AUCg","dah","SMOD_2023",
                    "Charge_log","habitat_code","habitat_label") %in% names(data_mod)))
    
    data_2023 <- data_mod %>%
      mutate(
        alpage        = factor(alpage),
        habitat_code  = as.integer(habitat_code)
      )
    
    # 4 habitats focus
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    focus_codes  <- c(1, 9, 5, 6)
    focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
    
    # =========================================================
    # Modèle GAM continu (pas de charge_bin)
    # - s(Charge_log) + interaction ti(Charge_log, SMOD_2023)
    # - alpage en intercept aléatoire
    # =========================================================
    form_cont <- formula(
      log(AUCg) ~
        s(dah,       k = 10, bs = "ts") +
        s(SMOD_2023, k = 10, bs = "ts") +
        ti(dah, SMOD_2023, k = c(20,20)) +
        s(alpage, bs = "re") +
        s(Charge_log, k = 10, bs = "ts") +
        ti(Charge_log, SMOD_2023, k = c(15,15))
    )
    
    n_th <- tryCatch(max(1L, parallel::detectCores(TRUE)-1L), error=function(e) 1L)
    
    fit_fun <- function(dat){
      dat <- dat %>% mutate(alpage = factor(alpage))
      bam(form_cont, data = dat, method="fREML", select=TRUE,
          discrete=TRUE, nthreads=n_th)
    }
    
    # =========================================================
    # Fit modèles (global + 4 habitats)
    # =========================================================
    mod_2023_global_cont <- fit_fun(data_2023)
    
    mods_hab_2023_cont <- map(setNames(as.list(focus_codes), as.character(focus_codes)), function(hc){
      di <- data_2023 %>% filter(habitat_code == hc)
      if (nrow(di) < 50) return(NULL)
      fit_fun(di)
    })
    
    # =========================================================
    # Taille d'effet Δ% entre Q25 → Q75 de Charge_log
    # (diff de prédiction moyenne sur l’échelle log, IC via V)
    # =========================================================
    delta_q25_q75_cont <- function(model, data_in, group_label="Total"){
      if (is.null(model) || nrow(data_in) < 10) {
        return(tibble(
          group = group_label,
          effect_log = NA_real_, se_log = NA_real_,
          lwr_log = NA_real_, upr_log = NA_real_,
          effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
          sig = FALSE, q_low = NA_real_, q_high = NA_real_, n = nrow(data_in)
        ))
      }
      ql <- as.numeric(quantile(data_in$Charge_log, 0.25, na.rm=TRUE))
      qh <- as.numeric(quantile(data_in$Charge_log, 0.75, na.rm=TRUE))
      if (!is.finite(ql) || !is.finite(qh) || ql >= qh) {
        ql <- as.numeric(quantile(data_in$Charge_log, 0.10, na.rm=TRUE))
        qh <- as.numeric(quantile(data_in$Charge_log, 0.90, na.rm=TRUE))
      }
      
      nd_lo <- data_in; nd_lo$Charge_log <- ql
      nd_hi <- data_in; nd_hi$Charge_log <- qh
      
      Xhi <- predict(model, newdata = nd_hi, type = "lpmatrix")
      Xlo <- predict(model, newdata = nd_lo, type = "lpmatrix")
      dbar <- colMeans(Xhi - Xlo)
      
      b <- coef(model); V <- vcov(model)
      eff <- as.numeric(drop(dbar %*% b))
      se  <- sqrt(as.numeric(drop(dbar %*% V %*% dbar)))
      lwr <- eff - 1.96*se; upr <- eff + 1.96*se
      
      tibble(
        group = group_label,
        effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
        effect_perc = 100*(exp(eff)-1),
        lwr_perc    = 100*(exp(lwr)-1),
        upr_perc    = 100*(exp(upr)-1),
        sig = (lwr > 0) | (upr < 0),
        q_low = ql, q_high = qh, n = nrow(data_in)
      )
    }
    
    # =========================================================
    # Extraction (2023)
    # =========================================================
    eff_total <- delta_q25_q75_cont(mod_2023_global_cont, data_2023, "Total")
    
    eff_habs <- imap_dfr(mods_hab_2023_cont, function(mod, code_chr){
      lab <- hab_labels_map[[code_chr]]
      di  <- data_2023 %>% filter(habitat_code == as.integer(code_chr))
      delta_q25_q75_cont(mod, di, lab)
    })
    
    eff_2023 <- bind_rows(eff_total, eff_habs) %>%
      mutate(group = factor(group, levels = c("Total", focus_labels)),
             alpha_plot = ifelse(sig, 1, 0.4))
    
    # =========================================================
    # Barplot
    # =========================================================
    cols_focus <- c(
      "Total"                        = "grey60",
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    pal <- cols_focus[levels(eff_2023$group)]
    
    p_mag_2023_cont <- ggplot(eff_2023, aes(x = group, y = effect_perc, fill = group)) +
      geom_col(width = 0.80, aes(alpha = alpha_plot)) +
      geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc), width = 0.20) +
      scale_fill_manual(values = pal, guide = "none") +
      scale_alpha_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(
        title = "Effet de la charge continue (Q25 → Q75) — 2023",
        subtitle = expression(Delta~"% entre "~Q[25]~" et "~Q[75]~" de "~Charge[log]),
        x = NULL, y = "Variation (%)"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 12))
    
    print(p_mag_2023_cont)
    
    # (optionnel) summary des modèles
    cat("\n=== GLOBAL (continu) ===\n"); print(summary(mod_2023_global_cont))
    purrr::iwalk(mods_hab_2023_cont, function(m, k){
      cat("\n=== HABITAT (continu):", hab_labels_map[[k]], "===\n")
      if (is.null(m)) cat("Modèle non estimé.\n") else print(summary(m))
    })
    
    # # Sauvegarde si besoin
    # ggsave("barplot_charge_continue_Q25Q75_2023.png", p_mag_2023_cont, width = 9, height = 5, dpi = 300)
    
    
    
    
    
    # =========================================================
    # Packages
    # =========================================================
    library(mgcv)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(ggplot2)
    
    # =========================================================
    # Données 2023 attendues: data_mod avec colonnes utilisées ci-dessous
    # =========================================================
    stopifnot(all(c("alpage","AUCg","dah","SMOD_2023","Charge_log",
                    "habitat_code","habitat_label") %in% names(data_mod)))
    
    data_2023 <- data_mod %>%
      mutate(alpage = factor(alpage),
             habitat_code = as.integer(habitat_code))
    
    # Habitats focus
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    focus_codes  <- c(1,9,5,6)
    focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
    
    # =========================================================
    # Modèle continu (intercept aléatoire alpage)
    # =========================================================
    form_cont <- formula(
      log(AUCg) ~
        s(dah,       k = 10, bs = "ts") +
        s(SMOD_2023, k = 10, bs = "ts") +
        ti(dah, SMOD_2023, k = c(20,20)) +
        s(alpage, bs = "re") +
        s(Charge_log, k = 10, bs = "ts") +
        ti(Charge_log, SMOD_2023, k = c(15,15))
    )
    
    n_th <- tryCatch(max(1L, parallel::detectCores(TRUE)-1L), error=function(e) 1L)
    fit_fun <- function(dat){
      dat <- dat %>% mutate(alpage = factor(alpage))
      bam(form_cont, data = dat, method = "fREML", select = TRUE,
          discrete = TRUE, nthreads = n_th)
    }
    
    # Modèles séparés par habitat (pour avoir une courbe par panneau)
    mods_hab_2023_cont <- map(setNames(as.list(focus_codes), as.character(focus_codes)), function(hc){
      di <- data_2023 %>% filter(habitat_code == hc)
      if (nrow(di) < 50) return(NULL)
      fit_fun(di)
    })
    
    # =========================================================
    # Courbe partielle de Charge_log (s + ti) + IC95% (lpmatrix)
    #  - On fixe dah et SMOD_2023 à leur médiane HABITAT
    #  - On exclut tous les termes sauf ceux qui contiennent Charge_log
    #  - On centre la courbe à 0 au plus petit x (effet relatif)
    #  - On renvoie variation (%) = 100*(exp(delta)-1)
    # =========================================================
    partial_charge_curve <- function(mod, dat, label, n = 200, trim = 0.01){
      if (is.null(mod) || nrow(dat) < 10) return(tibble())
      
      # grille "raisonnable" (1%-99%) du log-charge
      rng <- quantile(dat$Charge_log, probs = c(trim, 1-trim), na.rm = TRUE)
      x_seq <- seq(rng[[1]], rng[[2]], length.out = n)
      
      # covariables représentatives (médianes habitat)
      dah_med  <- median(dat$dah,       na.rm = TRUE)
      smod_med <- median(dat$SMOD_2023, na.rm = TRUE)
      
      # alpage arbitraire (RE intercept retiré ensuite)
      ref_alp <- levels(dat$alpage)[1]
      
      nd <- tibble::tibble(
        Charge_log = x_seq,
        dah        = dah_med,
        SMOD_2023  = smod_med,
        alpage     = factor(ref_alp, levels = levels(dat$alpage))
      )
      
      # Lp-matrix
      X  <- predict(mod, newdata = nd, type = "lpmatrix")
      b  <- coef(mod); V <- vcov(mod)
      
      # Colonnes à garder = coefficients des smooths qui impliquent "Charge_log"
      keep_cols <- rep(FALSE, length(b))
      for (sm in mod$smooth) {
        lab <- sm$label
        if (grepl("Charge_log", lab, fixed = TRUE)) {
          keep_cols[sm$first.para:sm$last.para] <- TRUE
        }
      }
      X_keep <- X
      X_keep[, !keep_cols] <- 0
      
      fit <- as.numeric(X_keep %*% b)
      se  <- sqrt(pmax(0, rowSums((X_keep %*% V) * X_keep)))
      
      # Centrage au x le plus bas (effet relatif)
      fit0 <- fit[1]
      dlt  <- fit  - fit0
      lo   <- (fit - 1.96*se) - fit0
      hi   <- (fit + 1.96*se) - fit0
      
      tibble::tibble(
        habitat = label,
        Charge_log = x_seq,
        eff_pct  = 100*(exp(dlt) - 1),
        lo_pct   = 100*(exp(lo)  - 1),
        hi_pct   = 100*(exp(hi)  - 1)
      ) %>%
        # point d'optimum (max de la courbe)
        mutate(opt = ifelse(eff_pct == max(eff_pct, na.rm = TRUE), TRUE, FALSE))
    }
    
    curves <- imap_dfr(mods_hab_2023_cont, function(mod, code_chr){
      lab <- hab_labels_map[[code_chr]]
      dat <- data_2023 %>% filter(habitat_code == as.integer(code_chr))
      partial_charge_curve(mod, dat, lab)
    })
    
    # Vérifie qu'on a des courbes
    if (nrow(curves) == 0L) stop("Pas assez de données pour tracer les courbes.")
    
    # Position de l'optimum (une ligne verticale par panneau)
    opt_df <- curves %>%
      group_by(habitat) %>%
      slice_max(order_by = eff_pct, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(habitat, x_opt = Charge_log, y_opt = eff_pct)
    
    # Palette
    cols_focus <- c(
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    
    # =========================================================
    # PLOT
    # =========================================================
    p_curves <- ggplot(curves, aes(Charge_log, eff_pct, colour = habitat, fill = habitat)) +
      geom_ribbon(aes(ymin = lo_pct, ymax = hi_pct), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 1) +
      geom_vline(data = opt_df, aes(xintercept = x_opt, colour = habitat),
                 linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
      facet_wrap(~ habitat, scales = "free_x") +
      scale_colour_manual(values = cols_focus, guide = "none") +
      scale_fill_manual(values = cols_focus, guide = "none") +
      geom_hline(yintercept = 0, linewidth = 0.4) +
      labs(
        title    = "Courbe partielle de l'effet de la charge — 2023",
        subtitle = "Contribution de s(Charge_log) + ti(Charge_log, SMOD_2023) (SMOD fixé à sa médiane habitat).\nVariation (%) relative au niveau de charge le plus bas.",
        x = "Charge (log-transform)", y = "Variation (%) vs charge minimale"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"))
    
    print(p_curves)
    
    # (Optionnel) tableau des optima
    opt_df %>%
      mutate(Charge_raw_at_opt = exp(x_opt) - 1) %>%
      arrange(habitat) %>%
      print(n = Inf)
    
    
    
    
    
    
    
    # =========================================================
    # Packages
    # =========================================================
    library(dplyr)
    library(ggplot2)
    library(purrr)
    library(tidyr)
    library(scales)
    
    # =========================================================
    # Données 2023
    # =========================================================
    stopifnot(all(c("alpage","AUCg","dah","SMOD_2023",
                    "Charge_log","habitat_code","habitat_label") %in% names(data_mod)))
    
    data_2023 <- data_mod %>%
      mutate(alpage = factor(alpage),
             habitat_code = as.integer(habitat_code))
    
    # Habitats focus
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    focus_codes  <- c(1, 9, 5, 6)
    focus_levels <- unname(hab_labels_map[as.character(focus_codes)])
    
    df <- data_2023 %>%
      filter(habitat_code %in% focus_codes) %>%
      mutate(habitat = factor(hab_labels_map[as.character(habitat_code)],
                              levels = focus_levels))
    
    # =========================================================
    # Tableau synthèse & détection "queue pauvre" à Q75
    # =========================================================
    thr      <- log(2)     # seuil "chargé"
    width_x  <- 0.25       # fenêtre locale (log1p)
    n_min    <- 30         # nb min d'observations dans la fenêtre locale
    prop_min <- 0.01       # ≥1% de l'échantillon dans la fenêtre
    
    stats_ch <- df %>%
      group_by(habitat) %>%
      summarise(
        n          = n(),
        n_pos      = sum(Charge_log > 0, na.rm = TRUE),
        q10        = quantile(Charge_log, 0.10, na.rm = TRUE),
        q25        = quantile(Charge_log, 0.25, na.rm = TRUE),
        q50        = quantile(Charge_log, 0.50, na.rm = TRUE),
        q75        = quantile(Charge_log, 0.75, na.rm = TRUE),
        q90        = quantile(Charge_log, 0.90, na.rm = TRUE),
        min_x      = min(Charge_log, na.rm = TRUE),
        max_x      = max(Charge_log, na.rm = TRUE),
        n_loc_q25  = sum(abs(Charge_log - q25) <= width_x, na.rm = TRUE),
        n_loc_q75  = sum(abs(Charge_log - q75) <= width_x, na.rm = TRUE),
        prop_q25   = n_loc_q25 / n,
        prop_q75   = n_loc_q75 / n,
        sparse_q25 = (n_loc_q25 < n_min) | (prop_q25 < prop_min),
        sparse_q75 = (n_loc_q75 < n_min) | (prop_q75 < prop_min),
        p_above_thr = mean(Charge_log >= thr),
        .groups = "drop"
      ) %>%
      mutate(
        thr_log   = thr,
        q25_raw   = exp(q25) - 1,
        q75_raw   = exp(q75) - 1,
        thr_raw   = exp(thr) - 1
      )
    
    print(stats_ch %>% arrange(habitat), n = Inf)
    
    # =========================================================
    # Histogrammes facettés + lignes Q25/Q50/Q75 + seuil
    # =========================================================
    p_hist <- ggplot(df, aes(x = Charge_log)) +
      geom_histogram(bins = 40, fill = "grey70", colour = "white") +
      geom_vline(data = stats_ch, aes(xintercept = q25), colour = "#1f77b4",
                 linetype = "dashed", linewidth = 0.7) +
      geom_vline(data = stats_ch, aes(xintercept = q50), colour = "#000000",
                 linetype = "solid",  linewidth = 0.7) +
      geom_vline(data = stats_ch, aes(xintercept = q75), colour = "#d62728",
                 linetype = "dashed", linewidth = 0.7) +
      geom_vline(xintercept = thr, colour = "black", linetype = "dotted", linewidth = 0.7) +
      facet_wrap(~ habitat, scales = "free_y") +
      labs(
        title = "Distribution de la charge (log) — 2023",
        subtitle = "Traits = Q25 (bleu pointillé), Médiane (noir), Q75 (rouge pointillé), Seuil chargé log(2) (noir pointillé fin).",
        x = "Charge_log", y = "Comptes"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"))
    print(p_hist)
    
    # =========================================================
    # Densités + rug (pour visualiser les queues) — même facettage
    # =========================================================
    p_den <- ggplot(df, aes(x = Charge_log)) +
      geom_density(linewidth = 0.9, alpha = 0.2, fill = "grey70") +
      geom_rug(sides = "b", alpha = 0.15) +
      geom_vline(data = stats_ch, aes(xintercept = q25), colour = "#1f77b4",
                 linetype = "dashed", linewidth = 0.7) +
      geom_vline(data = stats_ch, aes(xintercept = q50), colour = "#000000",
                 linetype = "solid",  linewidth = 0.7) +
      geom_vline(data = stats_ch, aes(xintercept = q75), colour = "#d62728",
                 linetype = "dashed", linewidth = 0.7) +
      geom_vline(xintercept = thr, colour = "black", linetype = "dotted", linewidth = 0.7) +
      facet_wrap(~ habitat, scales = "free_y") +
      labs(
        title = "Densité & rug de Charge_log — 2023",
        subtitle = "Vérifie si Q75 tombe dans une queue pauvre (faible densité).",
        x = "Charge_log", y = "Densité"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"))
    print(p_den)
    
    # =========================================================
    # (Optionnel) drapeau visuel "Q75 dans queue pauvre"
    # =========================================================
    lab_sparse <- stats_ch %>%
      transmute(
        habitat, x = q75, y = 0,  # y sera replacé par l'échelle de chaque facet avec annotate
        label = ifelse(sparse_q75,
                       paste0("Q75 sparse\nn≈", n_loc_q75, " (", percent(prop_q75, 0.1), ")"),
                       paste0("Q75 ok\nn≈", n_loc_q75, " (", percent(prop_q75, 0.1), ")"))
      )
    # À utiliser si tu veux annoter manuellement sur p_den/p_hist (les axes 'y' variant selon facet).
    
    
    
    
    
    
    
    
    
    # =========================================================
    # Packages
    # =========================================================
    library(mgcv)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(ggplot2)
    
    # =========================================================
    # Données 2023
    # =========================================================
    stopifnot(all(c("alpage","AUCg","dah","SMOD_2023","Charge_log",
                    "habitat_code","habitat_label") %in% names(data_mod)))
    
    data_2023 <- data_mod %>%
      mutate(alpage = factor(alpage),
             habitat_code = as.integer(habitat_code))
    
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    focus_codes  <- c(1,9,5,6)
    focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
    
    # =========================================================
    # Modèle continu (intercept aléatoire alpage)
    # =========================================================
    form_cont <- formula(
      log(AUCg) ~
        s(dah,       k = 10, bs = "ts") +
        s(SMOD_2023, k = 10, bs = "ts") +
        ti(dah, SMOD_2023, k = c(20,20)) +
        s(alpage, bs = "re") +
        s(Charge_log, k = 10, bs = "ts") +
        ti(Charge_log, SMOD_2023, k = c(15,15))
    )
    
    n_th <- tryCatch(max(1L, parallel::detectCores(TRUE)-1L), error=function(e) 1L)
    fit_fun <- function(dat){
      dat <- dat %>% mutate(alpage = factor(alpage))
      bam(form_cont, data = dat, method="fREML", select=TRUE,
          discrete=TRUE, nthreads=n_th)
    }
    
    # Fit global + 4 habitats (si tu les as déjà, tu peux garder les tiens)
    mod_2023_global_cont <- fit_fun(data_2023)
    mods_hab_2023_cont <- purrr::map(setNames(as.list(focus_codes), as.character(focus_codes)), function(hc){
      di <- data_2023 %>% filter(habitat_code == hc)
      if (nrow(di) < 50) return(NULL)
      fit_fun(di)
    })
    
    # =========================================================
    # Outils : effets entre 2 points et élasticité locale +10%
    # =========================================================
    thr <- log(2)
    log_plus_p <- function(x_log, p = 0.10) log1p((1+p) * pmax(expm1(x_log), 0))
    
    effect_between_points <- function(model, d, x_from, x_to){
      nd_a <- d; nd_a$Charge_log <- x_from
      nd_b <- d; nd_b$Charge_log <- x_to
      Xb <- predict(model, newdata = nd_b, type = "lpmatrix")
      Xa <- predict(model, newdata = nd_a, type = "lpmatrix")
      dbar <- colMeans(Xb - Xa)
      b <- coef(model); V <- vcov(model)
      eff <- as.numeric(drop(dbar %*% b))
      se  <- sqrt(as.numeric(drop(dbar %*% V %*% dbar)))
      lwr <- eff - 1.96*se; upr <- eff + 1.96*se
      list(eff_log = eff, lwr_log = lwr, upr_log = upr)
    }
    
    effect_plus10_at <- function(model, d, x_star, p = 0.10){
      effect_between_points(model, d, x_star, log_plus_p(x_star, p))
    }
    
    # =========================================================
    # Choix ROBUSTE des points (post-seuil) + fallback +10%
    # =========================================================
    robust_effect <- function(model, data_in, label,
                              q_lo = 0.40, q_hi = 0.60,     # central 40→60
                              width_x = 0.25,               # fenêtre locale
                              n_min = 30, prop_min = 0.01,  # support local
                              delta_thr = 0.02,             # un petit cran au-dessus du seuil
                              p_inc = 0.10){                # +10% si fallback
      if (is.null(model) || nrow(data_in) < 10) return(tibble())
      
      d0 <- data_in %>% filter(is.finite(Charge_log))
      post <- d0$Charge_log[d0$Charge_log >= thr + delta_thr]
      if (length(post) < 10) {
        return(tibble(group=label, method="insuffisant", effect_perc=NA_real_,
                      lwr_perc=NA_real_, upr_perc=NA_real_, x_from=NA_real_, x_to=NA_real_,
                      kept="none", n=nrow(d0)))
      }
      
      # 1) points centraux post-seuil
      x1 <- as.numeric(quantile(post, q_lo, na.rm = TRUE))
      x2 <- as.numeric(quantile(post, q_hi, na.rm = TRUE))
      
      n1 <- sum(abs(d0$Charge_log - x1) <= width_x)
      n2 <- sum(abs(d0$Charge_log - x2) <= width_x)
      ok <- (n1 >= n_min && n2 >= n_min &&
               n1/nrow(d0) >= prop_min && n2/nrow(d0) >= prop_min)
      
      if (ok) {
        e <- effect_between_points(model, d0, x1, x2)
        return(tibble(
          group  = label, method = "P40→P60 post",
          effect_perc = 100*(exp(e$eff_log)-1),
          lwr_perc    = 100*(exp(e$lwr_log)-1),
          upr_perc    = 100*(exp(e$upr_log)-1),
          x_from = x1, x_to = x2, kept = sprintf("n1=%d,n2=%d", n1, n2),
          n = nrow(d0)
        ))
      }
      
      # 2) fallback : élasticité +10% au médian post-seuil
      x_star <- as.numeric(quantile(post, 0.50, na.rm = TRUE))
      e <- effect_plus10_at(model, d0, x_star, p = p_inc)
      tibble(
        group  = label, method = "+10% @ P50 post",
        effect_perc = 100*(exp(e$eff_log)-1),
        lwr_perc    = 100*(exp(e$lwr_log)-1),
        upr_perc    = 100*(exp(e$upr_log)-1),
        x_from = x_star, x_to = log_plus_p(x_star, p_inc),
        kept = "fallback", n = nrow(d0)
      )
    }
    
    # =========================================================
    # Extraction (global + 4 habitats)
    # =========================================================
    eff_total <- robust_effect(mod_2023_global_cont, data_2023, "Total")
    
    eff_habs <- purrr::imap_dfr(mods_hab_2023_cont, function(mod, code_chr){
      lab <- hab_labels_map[[code_chr]]
      di  <- data_2023 %>% filter(habitat_code == as.integer(code_chr))
      robust_effect(mod, di, lab)
    })
    
    eff_rb <- bind_rows(eff_total, eff_habs) %>%
      mutate(group = factor(group, levels = c("Total", focus_labels)),
             lty   = ifelse(grepl("fallback", kept), "twodash", "solid"),
             alpha = 1)
    
    # =========================================================
    # Barplot ROBUSTE
    # =========================================================
    cols_focus <- c(
      "Total"                        = "grey60",
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    pal <- cols_focus[levels(eff_rb$group)]
    
    p_rb <- ggplot(eff_rb, aes(x = group, y = effect_perc, fill = group)) +
      geom_col(width = 0.80, aes(alpha = alpha)) +
      geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc, linetype = lty), width = 0.20) +
      scale_fill_manual(values = pal, guide = "none") +
      scale_linetype_identity(name = "Méthode",
                              guide = guide_legend(override.aes = list(color="black"))) +
      scale_alpha_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(
        title = "Effet de la charge — version robuste (2023)",
        subtitle = "Contraste P40→P60 du post-seuil si bien supporté ; sinon fallback = élasticité locale +10% au médian post-seuil.\nBarres = IC95% via lpmatrix.",
        x = NULL, y = "Variation (%)"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "right",
            axis.text.x = element_text(size = 12))
    
    print(p_rb)
    
    # =========================================================
    # (Optionnel) Diagnostic : où sont les points utilisés ?
    # =========================================================
    diag_df <- data_2023 %>%
      filter(habitat_code %in% focus_codes) %>%
      mutate(habitat = factor(hab_labels_map[as.character(habitat_code)], levels = focus_labels))
    
    used_lines <- eff_rb %>%
      filter(!is.na(x_from)) %>%
      mutate(habitat = factor(group, levels = focus_labels)) %>%
      select(habitat, x_from, x_to, method)
    
    p_diag <- ggplot(diag_df, aes(Charge_log)) +
      geom_density(fill="grey85") +
      geom_rug(alpha=.12, sides="b") +
      geom_vline(xintercept = thr, linetype="dotted") +
      geom_vline(data = used_lines, aes(xintercept = x_from, colour = habitat), linetype="dashed") +
      geom_vline(data = used_lines, aes(xintercept = x_to,   colour = habitat), linetype="dashed") +
      facet_wrap(~ habitat, scales = "free_y") +
      scale_colour_manual(values = pal, guide = "none") +
      labs(title = "Points utilisés pour l'estimation robuste — 2023",
           subtitle = "Traits pointillés = bornes du contraste ou du +10% fallback ; pointillé noir = seuil log(2).",
           x = "Charge_log", y = "Densité") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"))
    print(p_diag)
    
    # (Optionnel) tableau récap
    eff_rb %>%
      mutate(x_from_raw = exp(x_from) - 1,
             x_to_raw   = exp(x_to) - 1) %>%
      select(group, method, effect_perc, lwr_perc, upr_perc,
             x_from, x_to, x_from_raw, x_to_raw, kept, n) %>%
      arrange(group) %>% print(n = Inf)
    
    
    
    
    
    
    
    
    
    
    
    # =========================================================
    # Packages
    # =========================================================
    library(lme4)
    library(dplyr)
    library(purrr)
    library(ggplot2)
    
    # =========================================================
    # Données 2023
    # =========================================================
    stopifnot(all(c("alpage","AUCg","dah","SMOD_2023",
                    "Charge_log","habitat_code","habitat_label") %in% names(data_mod)))
    
    data_2023 <- data_mod %>%
      mutate(
        alpage       = factor(alpage),
        habitat_code = as.integer(habitat_code)
      )
    
    # Habitats focus
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    focus_codes  <- c(1,9,5,6)
    focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
    
    # =========================================================
    # Modèle linéaire mixte SANS interaction avec la neige
    #   log(AUCg) ~ dah + SMOD_2023 + Charge_log + (1|alpage)
    # =========================================================
    form_lmm_noSnowInt <- log(AUCg) ~ dah + SMOD_2023 + Charge_log + (1|alpage)
    
    m_lmm_global <- lmer(form_lmm_noSnowInt, data = data_2023, REML = TRUE)
    
    mods_lmm_hab <- purrr::map(setNames(as.list(focus_codes), as.character(focus_codes)), function(hc){
      di <- dplyr::filter(data_2023, habitat_code == hc)
      # Besoin de >= 2 alpages pour estimer un intercept aléatoire
      if (nrow(di) < 50 || length(unique(di$alpage)) < 2) return(NULL)
      lmer(form_lmm_noSnowInt, data = di, REML = TRUE)
    })
    
    cat("\n=== SUMMARY GLOBAL (sans interaction neige) ===\n")
    print(summary(m_lmm_global))
    
    # =========================================================
    # Effets basés UNIQUEMENT sur les Estimates
    #  Ici, sans interaction, l'effet pour un Δ de Charge_log est :
    #     Δ * β_Charge_log
    #  Var = Δ^2 * Var(β_Charge_log)
    # =========================================================
    coef_effect_only <- function(mod, delta){
      b  <- fixef(mod)["Charge_log"]
      se <- sqrt(vcov(mod)["Charge_log","Charge_log"])
      eff_log <- delta * b
      se_log  <- abs(delta) * se
      lwr_log <- eff_log - 1.96*se_log
      upr_log <- eff_log + 1.96*se_log
      c(eff_log = eff_log, lwr_log = lwr_log, upr_log = upr_log)
    }
    
    log_plus_p <- function(x_log, p=0.10) log1p((1+p)*pmax(expm1(x_log),0))
    
    # Q25 → Q75
    eff_Q25_Q75 <- function(mod, d, label){
      ql <- as.numeric(quantile(d$Charge_log, 0.25, na.rm=TRUE))
      qh <- as.numeric(quantile(d$Charge_log, 0.75, na.rm=TRUE))
      if (!is.finite(ql) || !is.finite(qh) || ql >= qh) {
        ql <- as.numeric(quantile(d$Charge_log, 0.10, na.rm=TRUE))
        qh <- as.numeric(quantile(d$Charge_log, 0.90, na.rm=TRUE))
      }
      delta <- qh - ql
      e <- coef_effect_only(mod, delta)
      tibble::tibble(
        group = label, method = "Q25→Q75",
        effect_perc = 100*(exp(e["eff_log"])-1),
        lwr_perc    = 100*(exp(e["lwr_log"])-1),
        upr_perc    = 100*(exp(e["upr_log"])-1),
        x_from = ql, x_to = qh
      )
    }
    
    # +10 % au médian post-seuil (robuste)
    eff_plus10_P50post <- function(mod, d, label, p=0.10, thr=log(2), delta_thr=0.02){
      d_post <- d %>% filter(Charge_log >= thr + delta_thr)
      if (nrow(d_post) < 30) {
        return(tibble::tibble(group=label, method="+10%@P50post (insuff.)",
                              effect_perc=NA_real_, lwr_perc=NA_real_, upr_perc=NA_real_,
                              x_from=NA_real_, x_to=NA_real_))
      }
      x_med <- as.numeric(quantile(d_post$Charge_log, 0.50, na.rm=TRUE))
      delta <- log_plus_p(x_med, p) - x_med
      e <- coef_effect_only(mod, delta)
      tibble::tibble(
        group = label, method = "+10%@P50post",
        effect_perc = 100*(exp(e["eff_log"])-1),
        lwr_perc    = 100*(exp(e["lwr_log"])-1),
        upr_perc    = 100*(exp(e["upr_log"])-1),
        x_from = x_med, x_to = log_plus_p(x_med, p)
      )
    }
    
    # =========================================================
    # Extraction GLOBAL + 4 habitats
    # =========================================================
    # Q25→Q75
    eff_q   <- eff_Q25_Q75(m_lmm_global, data_2023, "Total")
    eff_q_h <- purrr::imap_dfr(mods_lmm_hab, function(mod, code_chr){
      lab <- hab_labels_map[[code_chr]]
      di  <- dplyr::filter(data_2023, habitat_code == as.integer(code_chr))
      if (is.null(mod)) tibble::tibble(group=lab, method="Q25→Q75 (insuff.)",
                                       effect_perc=NA, lwr_perc=NA, upr_perc=NA,
                                       x_from=NA, x_to=NA)
      else eff_Q25_Q75(mod, di, lab)
    })
    eff_q_all <- bind_rows(eff_q, eff_q_h) %>%
      mutate(group = factor(group, levels = c("Total", focus_labels)))
    
    # +10 % local
    eff10   <- eff_plus10_P50post(m_lmm_global, data_2023, "Total")
    eff10_h <- purrr::imap_dfr(mods_lmm_hab, function(mod, code_chr){
      lab <- hab_labels_map[[code_chr]]
      di  <- dplyr::filter(data_2023, habitat_code == as.integer(code_chr))
      if (is.null(mod)) tibble::tibble(group=lab, method="+10%@P50post (insuff.)",
                                       effect_perc=NA, lwr_perc=NA, upr_perc=NA,
                                       x_from=NA, x_to=NA)
      else eff_plus10_P50post(mod, di, lab)
    })
    eff10_all <- bind_rows(eff10, eff10_h) %>%
      mutate(group = factor(group, levels = c("Total", focus_labels)))
    
    # =========================================================
    # Barplots
    # =========================================================
    cols_focus <- c(
      "Total"                        = "grey60",
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    
    plot_eff <- function(df, title_txt){
      pal <- cols_focus[levels(df$group)]
      ggplot(df, aes(x = group, y = effect_perc, fill = group)) +
        geom_col(width = 0.80, na.rm = TRUE) +
        geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc), width = 0.20, na.rm = TRUE) +
        scale_fill_manual(values = pal, guide = "none") +
        geom_hline(yintercept = 0, linewidth = 0.5) +
        labs(title = title_txt, x = NULL, y = "Variation (%)") +
        theme_minimal(base_size = 14) +
        theme(panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 12))
    }
    
    p_q25q75 <- plot_eff(eff_q_all,  "Effet charge (linéaire, sans interaction neige) — Q25→Q75 (2023)")
    p_plus10 <- plot_eff(eff10_all,  "Effet charge (linéaire, sans interaction neige) — +10% @ médian post-seuil (2023)")
    
    print(p_q25q75)
    print(p_plus10)
    
    # (optionnel) tableaux
    # eff_q_all  %>% arrange(group) %>% print(n=Inf)
    # eff10_all  %>% arrange(group) %>% print(n=Inf)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ## MEME  CHOSE AVCE UNE AUGMENTTAION DE 10%
    library(dplyr)
    library(ggplot2)
    library(purrr)
    library(mgcv)
    
    # ---------- Mapping des 4 habitats ----------
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    focus_codes  <- c(1,9,5,6)
    focus_labels <- unname(hab_labels_map[as.character(focus_codes)])
    
    # ---------- Petit utilitaire : +10% sur l'échelle brute mais variable modèle = log1p ----------
    log_plus_p <- function(x_log, p = 0.10) {
      # x_log = log1p(Charge). On convertit -> Charge, on applique +p, puis on revient en log1p
      log1p((1 + p) * pmax(expm1(x_log), 0))
    }
    
    # ---------- Effet moyen d'une augmentation de +10% ----------
    # - on filtre aux charges > 0 (sinon +10% de 0 = 0, ça ”dilue” l’effet)
    # - on RE-CALCULE charge_bin dans nd0/nd1 (crucial pour les termes by=charge_bin)
    effect_plus10 <- function(model, data_in, group_label = "Total", p = 0.10) {
      d <- data_in %>% mutate(alpage = factor(alpage))
      d <- d %>% filter(Charge_log > 0)  # comme pour Q25→Q75
      
      if (nrow(d) < 10) {
        return(tibble::tibble(
          group="(insuffisant)", effect_log=NA_real_, se_log=NA_real_,
          lwr_log=NA_real_, upr_log=NA_real_,
          effect_perc=NA_real_, lwr_perc=NA_real_, upr_perc=NA_real_,
          sig=FALSE, n=nrow(d)
        ))
      }
      
      nd0 <- d
      nd0$charge_bin <- as.integer(nd0$Charge_log >= log(2))
      
      nd1 <- d
      nd1$Charge_log <- log_plus_p(d$Charge_log, p = p)
      nd1$charge_bin <- as.integer(nd1$Charge_log >= log(2))
      
      X1  <- predict(model, newdata = nd1, type = "lpmatrix")
      X0  <- predict(model, newdata = nd0, type = "lpmatrix")
      dbar <- colMeans(X1 - X0)
      
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
        n = nrow(d)
      )
    }
    
    # ---------- Couleurs ----------
    cols_focus <- c(
      "Total"                        = "grey60",
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    
    # ---------- Calcul des effets (+10%) ----------
    # (on suppose que "data_2023", "mod_2023_global" et "mods_hab_2023" existent déjà)
    eff_tot_10 <- effect_plus10(mod_2023_global, data_2023, group_label = "Total", p = 0.10)
    
    eff_hab_10 <- imap_dfr(mods_hab_2023, function(mod, code_chr) {
      di <- data_2023 %>% filter(habitat_code == as.integer(code_chr))
      lab <- hab_labels_map[[code_chr]]
      if (is.null(mod)) {
        tibble::tibble(
          group = lab,
          effect_log = NA_real_, se_log = NA_real_, lwr_log = NA_real_, upr_log = NA_real_,
          effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
          sig = FALSE, n = nrow(di)
        )
      } else {
        effect_plus10(mod, di, group_label = lab, p = 0.10)
      }
    })
    
    eff10_2023 <- bind_rows(eff_tot_10, eff_hab_10) %>%
      mutate(group = factor(group, levels = c("Total", focus_labels)))
    
    # ---------- Barplot ----------
    pal <- cols_focus[levels(eff10_2023$group)]
    
    p_mag_10 <- ggplot(eff10_2023, aes(x = group, y = effect_perc, fill = group)) +
      geom_col(width = 0.80, aes(alpha = ifelse(sig, 1, 0.4))) +
      geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc), width = 0.20) +
      scale_fill_manual(values = pal, limits = levels(eff10_2023$group), drop = FALSE) +
      scale_alpha_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(
        title = "Effet d'une augmentation de +10% de la charge — 2023",
        subtitle = "Δ% moyen (exp(Δlog)−1) ×100, IC95% via matrice de variance des coefficients",
        x = NULL, y = "Variation (%)"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(size = 12))
    
    print(p_mag_10)
    # ggsave("barplot_plus10_charge_2023.png", p_mag_10, width = 9, height = 5, dpi = 300)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # ===== Barplot magnitude de l'effet (Δ% Q25→Q75 de Charge_log) : Total + 4 habitats, par année =====
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(mgcv)
    
    # --- libellés & couleurs ---
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    hab_levels <- unname(hab_labels_map[c("1","9","5","6")])
    cols_focus <- c(
      "Total"                        = "grey60",
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    
    # --- listes modèles & data (doivent exister) ---
    mlist <- list(`2017`=mod_2017, `2018`=mod_2018, `2019`=mod_2019,
                  `2020`=mod_2020, `2021`=mod_2021, `2022`=mod_2022, `2023`=mod_2023)
    dlist <- list(`2017`=d2017, `2018`=d2018, `2019`=d2019,
                  `2020`=d2020, `2021`=d2021, `2022`=d2022, `2023`=d2023)
    
    # --- collecteur résultats ---
    eff_res <- tibble::tibble(
      YEAR = integer(), group = character(),
      effect_log = numeric(), se_log = numeric(),
      lwr_log = numeric(), upr_log = numeric(),
      effect_perc = numeric(), lwr_perc = numeric(), upr_perc = numeric(),
      sig = logical(), q_low = numeric(), q_high = numeric(), n = integer()
    )
    
    # --- extraction (petite boucle OK pour l'agrégation) ---
    for (yy in names(mlist)) {
      mY <- mlist[[yy]]
      dY <- dlist[[yy]]
      
      # ---------- TOTAL (tous habitats, vrai global) ----------
      pos_vals <- dY$Charge_log[dY$Charge_log > 0]
      if (length(pos_vals) >= 10) {
        ql <- as.numeric(quantile(pos_vals, 0.25, na.rm = TRUE))
        qh <- as.numeric(quantile(pos_vals, 0.75, na.rm = TRUE))
        if (!is.finite(ql) || !is.finite(qh) || ql >= qh) {
          ql <- as.numeric(quantile(pos_vals, 0.10, na.rm = TRUE))
          qh <- as.numeric(quantile(pos_vals, 0.90, na.rm = TRUE))
        }
        nd_lo <- dY; nd_lo$Charge_log <- ql
        nd_hi <- dY; nd_hi$Charge_log <- qh
        
        Xhi <- predict(mY, newdata = nd_hi, type = "lpmatrix")
        Xlo <- predict(mY, newdata = nd_lo, type = "lpmatrix")
        dbar <- colMeans(Xhi - Xlo)
        b <- coef(mY); V <- vcov(mY)
        
        eff <- as.numeric(dbar %*% b)
        se  <- sqrt(drop(dbar %*% V %*% dbar))
        lwr <- eff - 1.96*se; upr <- eff + 1.96*se
        
        eff_res <- bind_rows(eff_res, tibble::tibble(
          YEAR = as.integer(yy), group = "Total",
          effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
          effect_perc = 100*(exp(eff)-1),
          lwr_perc = 100*(exp(lwr)-1), upr_perc = 100*(exp(upr)-1),
          sig = (lwr > 0) | (upr < 0),
          q_low = ql, q_high = qh, n = nrow(dY)
        ))
      } else {
        eff_res <- bind_rows(eff_res, tibble::tibble(
          YEAR = as.integer(yy), group = "Total",
          effect_log = NA_real_, se_log = NA_real_, lwr_log = NA_real_, upr_log = NA_real_,
          effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
          sig = FALSE, q_low = NA_real_, q_high = NA_real_, n = nrow(dY)
        ))
      }
      
      # ---------- PAR HABITAT (utilise automatiquement global + déviation) ----------
      # Dans dY, on a les indicateurs h1,h9,h5,h6 posés lors du fit
      hab_specs <- list(
        "P. nivales"                  = dY$h1 == 1,
        "P. thermiques écorchées"     = dY$h9 == 1,
        "Queyrellins"                 = dY$h5 == 1,
        "Nardaies denses du subalpin" = dY$h6 == 1
      )
      
      for (hname in names(hab_specs)) {
        idx <- hab_specs[[hname]]
        dH <- dY[idx, , drop = FALSE]
        pos_h <- dH$Charge_log[dH$Charge_log > 0]
        if (nrow(dH) >= 30 && length(pos_h) >= 10) {
          ql <- as.numeric(quantile(pos_h, 0.25, na.rm = TRUE))
          qh <- as.numeric(quantile(pos_h, 0.75, na.rm = TRUE))
          if (!is.finite(ql) || !is.finite(qh) || ql >= qh) {
            ql <- as.numeric(quantile(pos_h, 0.10, na.rm = TRUE))
            qh <- as.numeric(quantile(pos_h, 0.90, na.rm = TRUE))
          }
          nd_lo <- dH; nd_lo$Charge_log <- ql
          nd_hi <- dH; nd_hi$Charge_log <- qh
          
          Xhi <- predict(mY, newdata = nd_hi, type = "lpmatrix")
          Xlo <- predict(mY, newdata = nd_lo, type = "lpmatrix")
          dbar <- colMeans(Xhi - Xlo)
          b <- coef(mY); V <- vcov(mY)
          
          eff <- as.numeric(dbar %*% b)
          se  <- sqrt(drop(dbar %*% V %*% dbar))
          lwr <- eff - 1.96*se; upr <- eff + 1.96*se
          
          eff_res <- bind_rows(eff_res, tibble::tibble(
            YEAR = as.integer(yy), group = hname,
            effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
            effect_perc = 100*(exp(eff)-1),
            lwr_perc = 100*(exp(lwr)-1), upr_perc = 100*(exp(upr)-1),
            sig = (lwr > 0) | (upr < 0),
            q_low = ql, q_high = qh, n = nrow(dH)
          ))
        } else {
          eff_res <- bind_rows(eff_res, tibble::tibble(
            YEAR = as.integer(yy), group = hname,
            effect_log = NA_real_, se_log = NA_real_, lwr_log = NA_real_, upr_log = NA_real_,
            effect_perc = NA_real_, lwr_perc = NA_real_, upr_perc = NA_real_,
            sig = FALSE, q_low = NA_real_, q_high = NA_real_, n = nrow(dH)
          ))
        }
      }
    }
    
    # --- Barplot : 5 barres par année ---
    eff_res <- eff_res %>%
      mutate(
        YEAR  = factor(YEAR, levels = 2017:2023),
        group = factor(group, levels = c("Total", hab_levels))
      )
    
    pal <- cols_focus[levels(eff_res$group)]
    
    p_eff_charge <- ggplot(eff_res, aes(x = YEAR, y = effect_perc, fill = group)) +
      geom_col(position = position_dodge2(width = 0.9, preserve = "single"),
               width = 0.85, aes(alpha = ifelse(sig, 1, 0.4))) +
      geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc),
                    position = position_dodge2(width = 0.9, preserve = "single"),
                    width = 0.22) +
      scale_fill_manual(values = pal, name = "Groupe", drop = FALSE) +
      scale_alpha_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(
        title = "Effet du gradient pastoral (Δ% entre Q25→Q75 de Charge_log)",
        subtitle = "Barres = IC95%. Transparence = non significatif.",
        x = "Année", y = "Variation de la production (%)"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "right")
    
    p_eff_charge
    
    # (optionnel) tableau résultats
    # eff_res %>% arrange(YEAR, group) %>% print(n = Inf)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ## VERSION AVEC DES MODEL SEPARER : 
    # ===== Modèles annuels : GLOBAL (tous habitats) + 4-habitats (par année) =====
    library(dplyr)
    library(mgcv)
    library(parallel)
    
    nthreads <- max(1, parallel::detectCores() - 1)
    
    # AUCg de secours
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
    
    # ========= 2017 =========
    d2017_all <- dt_legacy %>%
      filter(YEAR == 2017) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2017_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2017_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2017 <- dt_legacy %>%
      filter(YEAR == 2017, habitat_code %in% c(1,9,5,6)) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab4 = factor(as.character(habitat_code),
                      levels=c("1","9","5","6"),
                      labels=unname(hab_labels_map[c("1","9","5","6")]))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2017 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab4 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab4) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab4),
      data=d2017, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2018 =========
    d2018_all <- dt_legacy %>%
      filter(YEAR == 2018) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2018_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2018_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2018 <- dt_legacy %>%
      filter(YEAR == 2018, habitat_code %in% c(1,9,5,6)) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab4 = factor(as.character(habitat_code),
                      levels=c("1","9","5","6"),
                      labels=unname(hab_labels_map[c("1","9","5","6")]))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2018 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab4 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab4) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab4),
      data=d2018, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2019 =========
    d2019_all <- dt_legacy %>%
      filter(YEAR == 2019) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2019_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2019_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2019 <- dt_legacy %>%
      filter(YEAR == 2019, habitat_code %in% c(1,9,5,6)) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab4 = factor(as.character(habitat_code),
                      levels=c("1","9","5","6"),
                      labels=unname(hab_labels_map[c("1","9","5","6")]))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2019 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab4 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab4) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab4),
      data=d2019, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2020 =========
    d2020_all <- dt_legacy %>%
      filter(YEAR == 2020) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2020_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2020_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2020 <- dt_legacy %>%
      filter(YEAR == 2020, habitat_code %in% c(1,9,5,6)) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab4 = factor(as.character(habitat_code),
                      levels=c("1","9","5","6"),
                      labels=unname(hab_labels_map[c("1","9","5","6")]))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2020 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab4 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab4) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab4),
      data=d2020, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2021 =========
    d2021_all <- dt_legacy %>%
      filter(YEAR == 2021) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2021_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2021_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2021 <- dt_legacy %>%
      filter(YEAR == 2021, habitat_code %in% c(1,9,5,6)) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab4 = factor(as.character(habitat_code),
                      levels=c("1","9","5","6"),
                      labels=unname(hab_labels_map[c("1","9","5","6")]))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2021 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab4 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab4) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab4),
      data=d2021, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2022 =========
    d2022_all <- dt_legacy %>%
      filter(YEAR == 2022) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2022_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2022_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2022 <- dt_legacy %>%
      filter(YEAR == 2022, habitat_code %in% c(1,9,5,6)) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab4 = factor(as.character(habitat_code),
                      levels=c("1","9","5","6"),
                      labels=unname(hab_labels_map[c("1","9","5","6")]))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2022 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab4 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab4) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab4),
      data=d2022, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2023 =========
    d2023_all <- dt_legacy %>%
      filter(YEAR == 2023) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2023_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2023_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2023 <- dt_legacy %>%
      filter(YEAR == 2023, habitat_code %in% c(1,9,5,6)) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab4 = factor(as.character(habitat_code),
                      levels=c("1","9","5","6"),
                      labels=unname(hab_labels_map[c("1","9","5","6")]))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2023 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab4 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab4) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab4),
      data=d2023, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # Exemples d’inspection :
    # summary(mod_2017_global); summary(mod_2017)
    # plot(mod_2023_global, pages=1); plot(mod_2023, pages=1)
    # gam.check(mod_2022_global); gam.check(mod_2022)
    
    
    
    
    # ===== Barplot de la magnitude de l'effet (Δ% Q25→Q75 de Charge_log) =====
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(mgcv)
    
    # Libellés/couleurs (les 4 habitats focus)
    hab_labels_map <- c(
      "1" = "P. nivales",
      "9" = "P. thermiques écorchées",
      "5" = "Queyrellins",
      "6" = "Nardaies denses du subalpin"
    )
    hab_levels <- unname(hab_labels_map[c("1","9","5","6")])
    cols_focus <- c(
      "Total"                        = "grey60",
      "P. nivales"                   = "#3B5BDB",
      "P. thermiques écorchées"      = "#E03131",
      "Queyrellins"                  = "#F4A261",
      "Nardaies denses du subalpin"  = "#1B9E77"
    )
    
    # Listes des modèles et jeux de données déjà créés plus haut
    mglob <- list(`2017`=mod_2017_global, `2018`=mod_2018_global, `2019`=mod_2019_global,
                  `2020`=mod_2020_global, `2021`=mod_2021_global, `2022`=mod_2022_global, `2023`=mod_2023_global)
    dglob <- list(`2017`=d2017_all,      `2018`=d2018_all,      `2019`=d2019_all,
                  `2020`=d2020_all,      `2021`=d2021_all,      `2022`=d2022_all,      `2023`=d2023_all)
    
    mhab <- list(`2017`=mod_2017, `2018`=mod_2018, `2019`=mod_2019,
                 `2020`=mod_2020, `2021`=mod_2021, `2022`=mod_2022, `2023`=mod_2023)
    dhab <- list(`2017`=d2017,    `2018`=d2018,    `2019`=d2019,
                 `2020`=d2020,    `2021`=d2021,    `2022`=d2022,    `2023`=d2023)
    
    # Collecteur
    eff_res <- tibble::tibble(
      YEAR = integer(), group = character(),
      effect_log = numeric(), se_log = numeric(),
      lwr_log = numeric(), upr_log = numeric(),
      effect_perc = numeric(), lwr_perc = numeric(), upr_perc = numeric(),
      sig = logical(), q_low = numeric(), q_high = numeric(), n = integer()
    )
    
    # Extraction (boucle agrégative)
    for (yy in names(mglob)) {
      # ----- TOTAL (vrai global) -----
      mg <- mglob[[yy]]; dg <- dglob[[yy]]
      pos <- dg$Charge_log[dg$Charge_log > 0]
      if (length(pos) >= 10) {
        ql <- as.numeric(quantile(pos, 0.25, na.rm=TRUE))
        qh <- as.numeric(quantile(pos, 0.75, na.rm=TRUE))
        if (!is.finite(ql) || !is.finite(qh) || ql >= qh) {
          ql <- as.numeric(quantile(pos, 0.10, na.rm=TRUE))
          qh <- as.numeric(quantile(pos, 0.90, na.rm=TRUE))
        }
        nd_lo <- dg; nd_lo$Charge_log <- ql
        nd_hi <- dg; nd_hi$Charge_log <- qh
        
        Xhi <- predict(mg, newdata=nd_hi, type="lpmatrix")
        Xlo <- predict(mg, newdata=nd_lo, type="lpmatrix")
        dbar <- colMeans(Xhi - Xlo)
        b <- coef(mg); V <- vcov(mg)
        
        eff <- as.numeric(dbar %*% b)
        se  <- sqrt(drop(dbar %*% V %*% dbar))
        lwr <- eff - 1.96*se; upr <- eff + 1.96*se
        
        eff_res <- bind_rows(eff_res, tibble::tibble(
          YEAR = as.integer(yy), group = "Total",
          effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
          effect_perc = 100*(exp(eff)-1),
          lwr_perc = 100*(exp(lwr)-1), upr_perc = 100*(exp(upr)-1),
          sig = (lwr > 0) | (upr < 0),
          q_low = ql, q_high = qh, n = nrow(dg)
        ))
      } else {
        eff_res <- bind_rows(eff_res, tibble::tibble(
          YEAR=as.integer(yy), group="Total",
          effect_log=NA_real_, se_log=NA_real_, lwr_log=NA_real_, upr_log=NA_real_,
          effect_perc=NA_real_, lwr_perc=NA_real_, upr_perc=NA_real_,
          sig=FALSE, q_low=NA_real_, q_high=NA_real_, n=nrow(dg)
        ))
      }
      
      # ----- PAR HABITAT (sur le modèle habitat de la même année) -----
      mh <- mhab[[yy]]; dh <- dhab[[yy]]
      for (hname in hab_levels) {
        dH <- dh %>% filter(hab4 == hname)
        pos_h <- dH$Charge_log[dH$Charge_log > 0]
        if (nrow(dH) >= 30 && length(pos_h) >= 10) {
          ql <- as.numeric(quantile(pos_h, 0.25, na.rm=TRUE))
          qh <- as.numeric(quantile(pos_h, 0.75, na.rm=TRUE))
          if (!is.finite(ql) || !is.finite(qh) || ql >= qh) {
            ql <- as.numeric(quantile(pos_h, 0.10, na.rm=TRUE))
            qh <- as.numeric(quantile(pos_h, 0.90, na.rm=TRUE))
          }
          nd_lo <- dH; nd_lo$Charge_log <- ql
          nd_hi <- dH; nd_hi$Charge_log <- qh
          
          Xhi <- predict(mh, newdata=nd_hi, type="lpmatrix")
          Xlo <- predict(mh, newdata=nd_lo, type="lpmatrix")
          dbar <- colMeans(Xhi - Xlo)
          b <- coef(mh); V <- vcov(mh)
          
          eff <- as.numeric(dbar %*% b)
          se  <- sqrt(drop(dbar %*% V %*% dbar))
          lwr <- eff - 1.96*se; upr <- eff + 1.96*se
          
          eff_res <- bind_rows(eff_res, tibble::tibble(
            YEAR = as.integer(yy), group = hname,
            effect_log = eff, se_log = se, lwr_log = lwr, upr_log = upr,
            effect_perc = 100*(exp(eff)-1),
            lwr_perc = 100*(exp(lwr)-1), upr_perc = 100*(exp(upr)-1),
            sig = (lwr > 0) | (upr < 0),
            q_low = ql, q_high = qh, n = nrow(dH)
          ))
        } else {
          eff_res <- bind_rows(eff_res, tibble::tibble(
            YEAR=as.integer(yy), group=hname,
            effect_log=NA_real_, se_log=NA_real_, lwr_log=NA_real_, upr_log=NA_real_,
            effect_perc=NA_real_, lwr_perc=NA_real_, upr_perc=NA_real_,
            sig=FALSE, q_low=NA_real_, q_high=NA_real_, n=nrow(dH)
          ))
        }
      }
    }
    
    # Barplot : 5 barres par année
    eff_res <- eff_res %>%
      mutate(
        YEAR  = factor(YEAR, levels = 2017:2023),
        group = factor(group, levels = c("Total", hab_levels))
      )
    
    pal <- cols_focus[levels(eff_res$group)]
    
    p_mag <- ggplot(eff_res, aes(x = YEAR, y = effect_perc, fill = group)) +
      geom_col(position = position_dodge2(width = 0.9, preserve = "single"),
               width = 0.85, aes(alpha = ifelse(sig, 1, 0.4))) +
      geom_errorbar(aes(ymin = lwr_perc, ymax = upr_perc),
                    position = position_dodge2(width = 0.9, preserve = "single"),
                    width = 0.22) +
      scale_fill_manual(values = pal, name = "Groupe", drop = FALSE) +
      scale_alpha_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(
        title = "Effet du gradient pastoral sur la production (Δ% entre Q25→Q75)",
        subtitle = "Total = modèle global (tous habitats). Habitats = modèles séparés (1,9,5,6). Barres = IC95%.",
        x = "Année", y = "Variation (%)"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "right")
    
    p_mag
    
    # (optionnel) tableau
    # eff_res %>% arrange(YEAR, group) %>% print(n = Inf)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # ===== Δ% entre Q10→Q90 de Charge_log =====
    library(dplyr); library(tidyr); library(ggplot2); library(mgcv)
    
    # Libellés / couleurs
    hab_labels_map <- c("1"="P. nivales","9"="P. thermiques écorchées",
                        "5"="Queyrellins","6"="Nardaies denses du subalpin")
    hab_levels <- unname(hab_labels_map[c("1","9","5","6")])
    cols_focus <- c("Total"="grey60","P. nivales"="#3B5BDB","P. thermiques écorchées"="#E03131",
                    "Queyrellins"="#F4A261","Nardaies denses du subalpin"="#1B9E77")
    
    # Objets existants
    mglob <- list(`2017`=mod_2017_global, `2018`=mod_2018_global, `2019`=mod_2019_global,
                  `2020`=mod_2020_global, `2021`=mod_2021_global, `2022`=mod_2022_global, `2023`=mod_2023_global)
    dglob <- list(`2017`=d2017_all,      `2018`=d2018_all,      `2019`=d2019_all,
                  `2020`=d2020_all,      `2021`=d2021_all,      `2022`=d2022_all,      `2023`=d2023_all)
    mhab  <- list(`2017`=mod_2017, `2018`=mod_2018, `2019`=mod_2019,
                  `2020`=mod_2020, `2021`=mod_2021, `2022`=mod_2022, `2023`=mod_2023)
    dhab  <- list(`2017`=d2017,    `2018`=d2018,    `2019`=d2019,
                  `2020`=d2020,    `2021`=d2021,    `2022`=d2022,    `2023`=d2023)
    
    # Choix des quantiles (central 80%)
    q_lo <- 0.10; q_hi <- 0.90
    
    eff_res <- tibble::tibble(
      YEAR=integer(), group=character(),
      effect_log=numeric(), se_log=numeric(), lwr_log=numeric(), upr_log=numeric(),
      effect_perc=numeric(), lwr_perc=numeric(), upr_perc=numeric(),
      sig=logical(), q_low=numeric(), q_high=numeric(), n=integer()
    )
    
    for (yy in names(mglob)) {
      # --- TOTAL ---
      mg <- mglob[[yy]]; dg <- dglob[[yy]]
      pos <- dg$Charge_log[dg$Charge_log > 0]
      if (length(pos) >= 10) {
        ql <- as.numeric(quantile(pos, q_lo, na.rm=TRUE))
        qh <- as.numeric(quantile(pos, q_hi, na.rm=TRUE))
        nd_lo <- dg; nd_lo$Charge_log <- ql
        nd_hi <- dg; nd_hi$Charge_log <- qh
        
        Xhi <- predict(mg, newdata=nd_hi, type="lpmatrix")
        Xlo <- predict(mg, newdata=nd_lo, type="lpmatrix")
        dbar <- colMeans(Xhi - Xlo); b <- coef(mg); V <- vcov(mg)
        eff <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
        lwr <- eff - 1.96*se; upr <- eff + 1.96*se
        
        eff_res <- bind_rows(eff_res, tibble::tibble(
          YEAR=as.integer(yy), group="Total",
          effect_log=eff, se_log=se, lwr_log=lwr, upr_log=upr,
          effect_perc=100*(exp(eff)-1),
          lwr_perc=100*(exp(lwr)-1), upr_perc=100*(exp(upr)-1),
          sig=(lwr>0)|(upr<0), q_low=ql, q_high=qh, n=nrow(dg)
        ))
      }
      
      # --- HABITATS ---
      mh <- mhab[[yy]]; dh <- dhab[[yy]]
      for (hname in hab_levels) {
        dH <- dh %>% filter(hab4==hname)
        pos_h <- dH$Charge_log[dH$Charge_log > 0]
        if (nrow(dH)>=30 && length(pos_h)>=10) {
          ql <- as.numeric(quantile(pos_h, q_lo, na.rm=TRUE))
          qh <- as.numeric(quantile(pos_h, q_hi, na.rm=TRUE))
          nd_lo <- dH; nd_lo$Charge_log <- ql
          nd_hi <- dH; nd_hi$Charge_log <- qh
          
          Xhi <- predict(mh, newdata=nd_hi, type="lpmatrix")
          Xlo <- predict(mh, newdata=nd_lo, type="lpmatrix")
          dbar <- colMeans(Xhi - Xlo); b <- coef(mh); V <- vcov(mh)
          eff <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
          lwr <- eff - 1.96*se; upr <- eff + 1.96*se
          
          eff_res <- bind_rows(eff_res, tibble::tibble(
            YEAR=as.integer(yy), group=hname,
            effect_log=eff, se_log=se, lwr_log=lwr, upr_log=upr,
            effect_perc=100*(exp(eff)-1),
            lwr_perc=100*(exp(lwr)-1), upr_perc=100*(exp(upr)-1),
            sig=(lwr>0)|(upr<0), q_low=ql, q_high=qh, n=nrow(dH)
          ))
        }
      }
    }
    
    # Plot
    eff_res <- eff_res %>%
      mutate(YEAR=factor(YEAR, levels=2017:2023),
             group=factor(group, levels=c("Total", hab_levels)))
    pal <- cols_focus[levels(eff_res$group)]
    
    p_q1090 <- ggplot(eff_res, aes(x=YEAR, y=effect_perc, fill=group)) +
      geom_col(position=position_dodge2(width=0.9, preserve="single"),
               width=0.85, aes(alpha=ifelse(sig,1,0.4))) +
      geom_errorbar(aes(ymin=lwr_perc, ymax=upr_perc),
                    position=position_dodge2(width=0.9, preserve="single"), width=0.22) +
      scale_fill_manual(values=pal, name="Groupe", drop=FALSE) +
      scale_alpha_identity() +
      geom_hline(yintercept=0, linewidth=0.5) +
      labs(title="Effet du gradient pastoral (Δ% entre Q10→Q90 de Charge_log)",
           subtitle="Barres = IC95%. Transparence = non significatif.",
           x="Année", y="Variation (%)") +
      theme_minimal(base_size=14) +
      theme(panel.grid.minor=element_blank(), legend.position="right")
    
    p_q1090
    
    
    
    
    
    
    
    
    # ===== Élasticité moyenne : % pour +10% de charge =====
    library(dplyr); library(ggplot2); library(mgcv)
    
    step <- log(1.10)  # +10% sur l'échelle log-charge
    ela_res <- tibble::tibble(
      YEAR=integer(), group=character(),
      elas=numeric(), se=numeric(), lwr=numeric(), upr=numeric(), n=integer()
    )
    
    for (yy in names(mglob)) {
      # --- TOTAL ---
      mg <- mglob[[yy]]; dg <- dglob[[yy]]
      d0 <- dg; d1 <- dg
      d1$Charge_log <- d0$Charge_log + step
      X1 <- predict(mg, newdata=d1, type="lpmatrix")
      X0 <- predict(mg, newdata=d0, type="lpmatrix")
      dbar <- colMeans((X1 - X0) / step); b <- coef(mg); V <- vcov(mg)
      est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
      # Conversion en % pour +10%
      pct10 <- 100*(exp(est*step)-1)
      lwr   <- 100*(exp((est-1.96*se)*step)-1)
      upr   <- 100*(exp((est+1.96*se)*step)-1)
      
      ela_res <- bind_rows(ela_res, tibble::tibble(
        YEAR=as.integer(yy), group="Total",
        elas=pct10, se=se, lwr=lwr, upr=upr, n=nrow(dg)
      ))
      
      # --- HABITATS ---
      mh <- mhab[[yy]]; dh <- dhab[[yy]]
      for (hname in hab_levels) {
        dH0 <- dh %>% filter(hab4==hname)
        if (nrow(dH0) >= 30) {
          dH1 <- dH0; dH1$Charge_log <- dH0$Charge_log + step
          X1 <- predict(mh, newdata=dH1, type="lpmatrix")
          X0 <- predict(mh, newdata=dH0, type="lpmatrix")
          dbar <- colMeans((X1 - X0) / step); b <- coef(mh); V <- vcov(mh)
          est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
          pct10 <- 100*(exp(est*step)-1)
          lwr   <- 100*(exp((est-1.96*se)*step)-1)
          upr   <- 100*(exp((est+1.96*se)*step)-1)
          
          ela_res <- bind_rows(ela_res, tibble::tibble(
            YEAR=as.integer(yy), group=hname,
            elas=pct10, se=se, lwr=lwr, upr=upr, n=nrow(dH0)
          ))
        }
      }
    }
    
    ela_res <- ela_res %>%
      mutate(YEAR=factor(YEAR, levels=2017:2023),
             group=factor(group, levels=c("Total", hab_levels)))
    pal <- c("Total"="grey60","P. nivales"="#3B5BDB","P. thermiques écorchées"="#E03131",
             "Queyrellins"="#F4A261","Nardaies denses du subalpin"="#1B9E77")
    
    p_elast <- ggplot(ela_res, aes(x=YEAR, y=elas, fill=group)) +
      geom_col(position=position_dodge2(width=0.9, preserve="single"),
               width=0.85) +
      geom_errorbar(aes(ymin=lwr, ymax=upr),
                    position=position_dodge2(width=0.9, preserve="single"), width=0.22) +
      scale_fill_manual(values=pal, name="Groupe", drop=FALSE) +
      geom_hline(yintercept=0, linewidth=0.5) +
      labs(title="Élasticité moyenne : effet d’un +10% de charge",
           subtitle="% de variation de la production pour +10% de charge (IC95%)",
           x="Année", y="% pour +10% de charge") +
      theme_minimal(base_size=14) +
      theme(panel.grid.minor=element_blank(), legend.position="right")
    
    p_elast
    
    
    
    
    
    
    
    
    
    # ===== % pour +10% de charge — EXACT, en repartant de la charge brute =====
    library(dplyr); library(ggplot2); library(mgcv); library(tidyr)
    
    # Libellés/couleurs
    hab_labels_map <- c("1"="P. nivales","9"="P. thermiques écorchées",
                        "5"="Queyrellins","6"="Nardaies denses du subalpin")
    hab_levels <- unname(hab_labels_map[c("1","9","5","6")])
    cols_focus <- c("Total"="grey60","P. nivales"="#3B5BDB","P. thermiques écorchées"="#E03131",
                    "Queyrellins"="#F4A261","Nardaies denses du subalpin"="#1B9E77")
    
    # Registres objets existants
    mglob <- list(`2017`=mod_2017_global, `2018`=mod_2018_global, `2019`=mod_2019_global,
                  `2020`=mod_2020_global, `2021`=mod_2021_global, `2022`=mod_2022_global, `2023`=mod_2023_global)
    dglob <- list(`2017`=d2017_all,      `2018`=d2018_all,      `2019`=d2019_all,
                  `2020`=d2020_all,      `2021`=d2021_all,      `2022`=d2022_all,      `2023`=d2023_all)
    mhab  <- list(`2017`=mod_2017, `2018`=mod_2018, `2019`=mod_2019,
                  `2020`=mod_2020, `2021`=mod_2021, `2022`=mod_2022, `2023`=mod_2023)
    dhab  <- list(`2017`=d2017,    `2018`=d2018,    `2019`=d2019,
                  `2020`=d2020,    `2021`=d2021,    `2022`=d2022,    `2023`=d2023)
    
    # Nom de la charge brute (présente dans ton summary)
    raw_col <- "chargement_median_2022-2024"
    
    # helper pour reconstruire Charge_log exactement comme au fit
    rebuild_charge_log <- function(d) {
      stopifnot(raw_col %in% names(d), "Charge_log" %in% names(d))
      raw <- pmax(d[[raw_col]], 0)
      # évalue quel log a été utilisé (log vs log1p)
      err_log  <- mean(abs(d$Charge_log - log(pmax(raw, .Machine$double.eps))), na.rm=TRUE)
      err_log1 <- mean(abs(d$Charge_log - log1p(raw)), na.rm=TRUE)
      if (!is.finite(err_log))  err_log  <- Inf
      if (!is.finite(err_log1)) err_log1 <- Inf
      if (err_log1 < err_log) {
        list(fun = function(x) log1p(pmax(x,0)), type="log1p")
      } else {
        list(fun = function(x) log(pmax(x, .Machine$double.eps)), type="log")
      }
    }
    
    eff10 <- tibble::tibble(
      YEAR=integer(), group=character(),
      pct10=numeric(), lwr=numeric(), upr=numeric(), n=integer(), trans=character()
    )
    
    for (yy in names(mglob)) {
      # ==== GLOBAL ====
      mg <- mglob[[yy]]; dg <- dglob[[yy]]
      if (!raw_col %in% names(dg)) stop("Colonne brute introuvable dans dglob[[",yy,"]]")
      tf <- rebuild_charge_log(dg)
      
      d0 <- dg
      d1 <- dg
      d0$Charge_log <- tf$fun(dg[[raw_col]])
      d1$Charge_log <- tf$fun(dg[[raw_col]] * 1.10)   # +10% exact sur la brute
      
      X1 <- predict(mg, newdata=d1, type="lpmatrix")
      X0 <- predict(mg, newdata=d0, type="lpmatrix")
      dbar <- colMeans(X1 - X0)                       # Δη moyen (log-échelle)
      b <- coef(mg); V <- vcov(mg)
      
      est <- as.numeric(dbar %*% b)                   # moyenne des Δη
      se  <- sqrt(drop(dbar %*% V %*% dbar))
      # passage en % sur GPROD :
      pct  <- 100*(exp(est) - 1)
      lwr  <- 100*(exp(est - 1.96*se) - 1)
      upr  <- 100*(exp(est + 1.96*se) - 1)
      
      eff10 <- bind_rows(eff10, tibble::tibble(
        YEAR=as.integer(yy), group="Total", pct10=pct, lwr=lwr, upr=upr, n=nrow(dg), trans=tf$type
      ))
      
      # ==== HABITATS ====
      mh <- mhab[[yy]]; dh <- dhab[[yy]]
      tfh <- rebuild_charge_log(dh)
      for (hname in hab_levels) {
        dH <- dh %>% filter(hab4 == hname)
        if (nrow(dH) < 30) next
        d0 <- dH; d1 <- dH
        d0$Charge_log <- tfh$fun(dH[[raw_col]])
        d1$Charge_log <- tfh$fun(dH[[raw_col]] * 1.10)
        
        X1 <- predict(mh, newdata=d1, type="lpmatrix")
        X0 <- predict(mh, newdata=d0, type="lpmatrix")
        dbar <- colMeans(X1 - X0); b <- coef(mh); V <- vcov(mh)
        est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
        pct <- 100*(exp(est) - 1)
        lwr <- 100*(exp(est - 1.96*se) - 1)
        upr <- 100*(exp(est + 1.96*se) - 1)
        
        eff10 <- bind_rows(eff10, tibble::tibble(
          YEAR=as.integer(yy), group=hname, pct10=pct, lwr=lwr, upr=upr, n=nrow(dH), trans=tfh$type
        ))
      }
    }
    
    eff10 <- eff10 %>%
      mutate(YEAR=factor(YEAR, levels=2017:2023),
             group=factor(group, levels=c("Total", hab_levels)))
    
    pal <- cols_focus[levels(eff10$group)]
    
    p_elast_exact <- ggplot(eff10, aes(YEAR, pct10, fill=group)) +
      geom_col(position=position_dodge2(width=0.9, preserve="single"), width=0.85) +
      geom_errorbar(aes(ymin=lwr, ymax=upr),
                    position=position_dodge2(width=0.9, preserve="single"), width=0.22) +
      scale_fill_manual(values=pal, name="Groupe", drop=FALSE) +
      geom_hline(yintercept=0, linewidth=0.5) +
      labs(title="Élasticité moyenne : effet d’un +10% de charge (exact, sur la charge brute)",
           subtitle="% de variation de la production attendue (GPROD) pour +10% de charge — IC95%",
           x="Année", y="% pour +10% de charge") +
      theme_minimal(base_size=14) +
      theme(panel.grid.minor=element_blank(), legend.position="right")
    
    p_elast_exact
    
    
    
    
    
    # =============================== #
    # 0) Prépa commune (6 habitats)   #
    # =============================== #
    library(dplyr)
    library(mgcv)
    library(parallel)
    library(tidyr)
    library(ggplot2)
    
    nthreads <- max(1, parallel::detectCores() - 1)
    
    # AUCg de secours
    if (!"AUCg" %in% names(dt_legacy) && "GPROD" %in% names(dt_legacy)) {
      dt_legacy$AUCg <- dt_legacy$GPROD
    }
    dt_legacy$AUCg <- pmax(dt_legacy$AUCg, .Machine$double.eps)
    
    # Liste des habitats à modéliser (1,9,5,6 + 17,21)
    codes_focus <- c("1","9","5","6","17","21")
    
    # Construire des libellés à partir des données (fallback si manquant)
    lab_from_data <- dt_legacy %>%
      filter(as.character(habitat_code) %in% codes_focus) %>%
      group_by(habitat_code, habitat_label) %>%
      tally(sort = TRUE) %>%
      group_by(habitat_code) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(code = as.character(habitat_code),
                label = ifelse(is.na(habitat_label) | habitat_label=="", paste0("Habitat ", habitat_code), habitat_label))
    hab_labels_map <- setNames(lab_from_data$label, lab_from_data$code)
    
    # Valeurs par défaut (si certains codes n'ont pas de label en base)
    defaults <- c("1"="P. nivales",
                  "9"="P. thermiques écorchées",
                  "5"="Queyrellins",
                  "6"="Nardaies denses du subalpin",
                  "17"="Habitat 17",
                  "21"="Habitat 21")
    for (k in names(defaults)) if (is.null(hab_labels_map[k]) || is.na(hab_labels_map[k])) hab_labels_map[k] <- defaults[k]
    
    hab_levels <- unname(hab_labels_map[codes_focus])
    
    # Palette (Total + 6 habitats)
    cols_focus <- setNames(
      c("grey60", "#3B5BDB", "#E03131", "#F4A261", "#1B9E77", "#6F42C1", "#8D6E63"),
      c("Total", hab_labels_map[c("1","9","5","6","17","21")])
    )
    
    # =============================== #
    # 1) Modèles annuels (global + H) #
    # =============================== #
    
    # ========= 2017 =========
    d2017_all <- dt_legacy %>%
      filter(YEAR == 2017) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2017_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2017_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2017 <- dt_legacy %>%
      filter(YEAR == 2017, as.character(habitat_code) %in% codes_focus) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab6 = factor(as.character(habitat_code),
                      levels=codes_focus,
                      labels=hab_levels)
      ) %>% filter(is.finite(Charge_log))
    
    mod_2017 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab6 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab6) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab6),
      data=d2017, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2018 =========
    d2018_all <- dt_legacy %>%
      filter(YEAR == 2018) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2018_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2018_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2018 <- dt_legacy %>%
      filter(YEAR == 2018, as.character(habitat_code) %in% codes_focus) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab6 = factor(as.character(habitat_code),
                      levels=codes_focus,
                      labels=hab_levels)
      ) %>% filter(is.finite(Charge_log))
    
    mod_2018 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab6 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab6) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab6),
      data=d2018, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2019 =========
    d2019_all <- dt_legacy %>%
      filter(YEAR == 2019) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2019_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2019_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2019 <- dt_legacy %>%
      filter(YEAR == 2019, as.character(habitat_code) %in% codes_focus) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab6 = factor(as.character(habitat_code),
                      levels=codes_focus,
                      labels=hab_levels)
      ) %>% filter(is.finite(Charge_log))
    
    mod_2019 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab6 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab6) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab6),
      data=d2019, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2020 =========
    d2020_all <- dt_legacy %>%
      filter(YEAR == 2020) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2020_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2020_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2020 <- dt_legacy %>%
      filter(YEAR == 2020, as.character(habitat_code) %in% codes_focus) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab6 = factor(as.character(habitat_code),
                      levels=codes_focus,
                      labels=hab_levels)
      ) %>% filter(is.finite(Charge_log))
    
    mod_2020 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab6 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab6) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab6),
      data=d2020, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2021 =========
    d2021_all <- dt_legacy %>%
      filter(YEAR == 2021) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2021_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2021_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2021 <- dt_legacy %>%
      filter(YEAR == 2021, as.character(habitat_code) %in% codes_focus) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab6 = factor(as.character(habitat_code),
                      levels=codes_focus,
                      labels=hab_levels)
      ) %>% filter(is.finite(Charge_log))
    
    mod_2021 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab6 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab6) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab6),
      data=d2021, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2022 =========
    d2022_all <- dt_legacy %>%
      filter(YEAR == 2022) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2022_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2022_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2022 <- dt_legacy %>%
      filter(YEAR == 2022, as.character(habitat_code) %in% codes_focus) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab6 = factor(as.character(habitat_code),
                      levels=codes_focus,
                      labels=hab_levels)
      ) %>% filter(is.finite(Charge_log))
    
    mod_2022 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab6 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab6) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab6),
      data=d2022, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # ========= 2023 =========
    d2023_all <- dt_legacy %>%
      filter(YEAR == 2023) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab_all = factor(as.character(habitat_code))
      ) %>% filter(is.finite(Charge_log))
    
    mod_2023_global <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab_all + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts") +
        ti(Charge_log, SMOD_2023_std, k=c(15,15)),
      data=d2023_all, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    d2023 <- dt_legacy %>%
      filter(YEAR == 2023, as.character(habitat_code) %in% codes_focus) %>%
      mutate(
        DAH_std = as.numeric(scale(DAH)),
        SMOD_std = as.numeric(scale(SMOD)),
        SMOD_2023_std = SMOD_std,
        alpage = factor(alpage),
        hab6 = factor(as.character(habitat_code),
                      levels=codes_focus,
                      labels=hab_levels)
      ) %>% filter(is.finite(Charge_log))
    
    mod_2023 <- bam(
      log(AUCg) ~
        s(DAH_std, k=10, bs="ts") +
        s(SMOD_2023_std, k=10, bs="ts") +
        ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
        hab6 + s(alpage, bs="re") +
        s(Charge_log, k=10, bs="ts", by=hab6) +
        ti(Charge_log, SMOD_2023_std, k=c(15,15), by=hab6),
      data=d2023, method="fREML", select=TRUE,
      discrete=TRUE, nthreads=nthreads, na.action=na.exclude
    )
    
    # =============================== #
    # 2) Effet ×2 (médiane -> 2× médiane)
    #    sans fixer SMOD/DAH (marginalisé) #
    # =============================== #
    
    mglob <- list(`2017`=mod_2017_global, `2018`=mod_2018_global, `2019`=mod_2019_global,
                  `2020`=mod_2020_global, `2021`=mod_2021_global, `2022`=mod_2022_global, `2023`=mod_2023_global)
    dglob <- list(`2017`=d2017_all,      `2018`=d2018_all,      `2019`=d2019_all,
                  `2020`=d2020_all,      `2021`=d2021_all,      `2022`=d2022_all,      `2023`=d2023_all)
    
    mhab  <- list(`2017`=mod_2017, `2018`=mod_2018, `2019`=mod_2019,
                  `2020`=mod_2020, `2021`=mod_2021, `2022`=mod_2022, `2023`=mod_2023)
    dhab  <- list(`2017`=d2017,    `2018`=d2018,    `2019`=d2019,
                  `2020`=d2020,    `2021`=d2021,    `2022`=d2022,    `2023`=d2023)
    
    raw_col <- "chargement_median_2022-2024"
    
    res2x <- tibble::tibble(
      YEAR=integer(), group=character(),
      pct2x=numeric(), lwr=numeric(), upr=numeric(), n=integer(),
      c0=numeric(), c1=numeric()
    )
    
    for (yy in names(mglob)) {
      # ----- GLOBAL (tous habitats) -----
      mg <- mglob[[yy]]; dg <- dglob[[yy]]
      lo <- quantile(pmax(dg[[raw_col]],0), 0.01, na.rm=TRUE)
      hi <- quantile(pmax(dg[[raw_col]],0), 0.99, na.rm=TRUE)
      med <- as.numeric(median(pmax(dg[[raw_col]],0), na.rm=TRUE))
      c0 <- max(lo, min(med, hi))
      c1 <- max(lo, min(2*med, hi))
      
      nd0 <- dg; nd1 <- dg
      nd0$Charge_log <- log1p(c0)
      nd1$Charge_log <- log1p(c1)
      
      X1 <- predict(mg, newdata=nd1, type="lpmatrix")
      X0 <- predict(mg, newdata=nd0, type="lpmatrix")
      dbar <- colMeans(X1 - X0); b <- coef(mg); V <- vcov(mg)
      est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
      pct <- 100*(exp(est)-1); lwr <- 100*(exp(est-1.96*se)-1); upr <- 100*(exp(est+1.96*se)-1)
      
      res2x <- bind_rows(res2x, tibble::tibble(
        YEAR=as.integer(yy), group="Total",
        pct2x=pct, lwr=lwr, upr=upr, n=nrow(dg), c0=c0, c1=c1
      ))
      
      # ----- PAR HABITAT (6 groupes) -----
      mh <- mhab[[yy]]; dh <- dhab[[yy]]
      for (hname in hab_levels) {
        dH <- dh %>% filter(hab6 == hname)
        if (nrow(dH) < 30) next
        loH <- quantile(pmax(dH[[raw_col]],0), 0.01, na.rm=TRUE)
        hiH <- quantile(pmax(dH[[raw_col]],0), 0.99, na.rm=TRUE)
        medH <- as.numeric(median(pmax(dH[[raw_col]],0), na.rm=TRUE))
        c0 <- max(loH, min(medH, hiH))
        c1 <- max(loH, min(2*medH, hiH))
        
        nd0 <- dH; nd1 <- dH
        nd0$Charge_log <- log1p(c0)
        nd1$Charge_log <- log1p(c1)
        
        X1 <- predict(mh, newdata=nd1, type="lpmatrix")
        X0 <- predict(mh, newdata=nd0, type="lpmatrix")
        dbar <- colMeans(X1 - X0); b <- coef(mh); V <- vcov(mh)
        est <- as.numeric(dbar %*% b); se <- sqrt(drop(dbar %*% V %*% dbar))
        pct <- 100*(exp(est)-1); lwr <- 100*(exp(est-1.96*se)-1); upr <- 100*(exp(est+1.96*se)-1)
        
        res2x <- bind_rows(res2x, tibble::tibble(
          YEAR=as.integer(yy), group=hname,
          pct2x=pct, lwr=lwr, upr=upr, n=nrow(dH), c0=c0, c1=c1
        ))
      }
    }
    
    # =============================== #
    # 3) Barplot (Total + 6 habitats) #
    # =============================== #
    res2x <- res2x %>%
      mutate(
        YEAR  = factor(YEAR, levels = 2017:2023),
        group = factor(group, levels = c("Total", hab_levels)),
        sig   = (lwr > 0) | (upr < 0)
      )
    
    pal <- cols_focus[levels(res2x$group)]
    
    p_2x <- ggplot(res2x, aes(YEAR, pct2x, fill = group)) +
      geom_col(position = position_dodge2(width = 0.9, preserve = "single"),
               width = 0.85, aes(alpha = ifelse(sig, 1, 0.4))) +
      geom_errorbar(aes(ymin = lwr, ymax = upr),
                    position = position_dodge2(width = 0.9, preserve = "single"),
                    width = 0.22) +
      scale_fill_manual(values = pal, name = "Groupe", drop = FALSE) +
      scale_alpha_identity() +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(
        title = "Effet d’un doublement de la charge (médiane → 2× médiane)",
        subtitle = "Sans fixer SMOD/DAH (marginalisé sur les covariables observées). Barres = IC95%.",
        x = "Année", y = "% de variation de la production (GPROD)"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "right")
    
    p_2x
    
    # (optionnel) Inspecter les valeurs utilisées
    # res2x %>% arrange(YEAR, group) %>% select(YEAR, group, n, c0, c1) %>% print(n=Inf)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # =========================
    # AIC : 2023 global vs ref
    # =========================
    library(mgcv)
    library(dplyr)
    
    # Threads
    if (!exists("nthreads")) nthreads <- 1L
    if (!exists("n_th"))     n_th     <- nthreads
    
    # 1) Modèle 2023 GLOBAL (d2023_all) — fREML
    if (!exists("mod_2023_global")) {
      mod_2023_global <- bam(
        log(AUCg) ~
          s(DAH_std, k=10, bs="ts") +
          s(SMOD_2023_std, k=10, bs="ts") +
          ti(DAH_std, SMOD_2023_std, k=c(20,20)) +
          hab_all + s(alpage, bs="re") +
          s(Charge_log, k=10, bs="ts") +
          ti(Charge_log, SMOD_2023_std, k=c(15,15)),
        data = d2023_all, method = "fREML", select = TRUE,
        discrete = TRUE, nthreads = nthreads, na.action = na.exclude
      )
    }
    
    # 2) Modèle 2023 RÉFÉRENCE (data_mod) — fREML
    m_ref <- bam(
      log(AUCg) ~
        s(dah, k=10, bs="ts") + s(SMOD_2023, k=10, bs="ts") +
        ti(dah, SMOD_2023, k=c(20,20)) +
        s(alpage, bs="re") +
        charge_bin +
        s(Charge_log, k=10, bs="ts", by=charge_bin) +
        ti(Charge_log, SMOD_2023, k=c(15,15), by=charge_bin),
      data = data_mod, method = "fREML", select = TRUE,
      discrete = TRUE, nthreads = n_th, na.action = na.exclude
    )
    
    # --- AIC (REML) ---
    aic_reml <- tibble::tibble(
      model = c("mod_2023_global [REML]", "m_ref [REML]"),
      data  = c("d2023_all", "data_mod"),
      n     = c(nrow(d2023_all), nrow(data_mod)),
      AIC   = c(AIC(mod_2023_global), AIC(m_ref))
    )
    print(aic_reml)
    
    # =========================================
    # (Option recommandé) Comparabilité en ML :
    # refit exactement les mêmes modèles en ML
    # =========================================
    mod_2023_global_ML <- bam(
      formula(mod_2023_global),
      data = d2023_all, method = "ML", select = TRUE,
      discrete = TRUE, nthreads = nthreads, na.action = na.exclude
    )
    
    m_ref_ML <- bam(
      formula(m_ref),
      data = data_mod, method = "ML", select = TRUE,
      discrete = TRUE, nthreads = n_th, na.action = na.exclude
    )
    
    aic_ml <- tibble::tibble(
      model = c("mod_2023_global [ML]", "m_ref [ML]"),
      data  = c("d2023_all", "data_mod"),
      n     = c(nrow(d2023_all), nrow(data_mod)),
      AIC   = c(AIC(mod_2023_global_ML), AIC(m_ref_ML))
    )
    print(aic_ml)
    
    cat("\nNOTE: L’AIC est comparable seulement entre modèles fités sur le même jeu de données (mêmes observations/réponse).\n",
        "Ici, d2023_all vs data_mod diffèrent vraisemblablement : compare plutôt au sein d’un même dataset,\n",
        "ou utilise une évaluation commune (CV, test set) si tu veux trancher entre global et référence.\n", sep = "")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  }
  
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Chargement des librairies nécessaires
  library(mgcv)
  library(knitr)
  library(kableExtra)
  library(dplyr)
  library(webshot2)  # Pour l'export PNG
  
  # Votre modèle GAM
  m_ref_plus <- bam(
    log(AUCg) ~ s(dah, k=10, bs="ts") + 
      s(SMOD_2023, k=10, bs="ts") + 
      ti(dah, SMOD_2023, k=c(15,15)) + 
      s(alpage, bs="re") + 
      charge_bin + 
      s(Charge_log, k=10, bs="ts", by=charge_bin) + 
      ti(Charge_log, SMOD_2023, k=c(10,10), by=charge_bin) + 
      s(alpage, bs="re", by=charge_bin) + 
      ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
    data = data_mod, 
    method = "fREML", 
    select = TRUE, 
    discrete = TRUE, 
    nthreads = n_th
  )
  
  # Extraction des résultats du modèle
  model_summary <- summary(m_ref_plus)
  
  # Fonction pour extraire et formatter les résultats
  extract_gam_results <- function(model_summary) {
    # Termes paramétriques
    parametric <- data.frame(
      Term = rownames(model_summary$p.table),
      Estimate = model_summary$p.table[, "Estimate"],
      SE = model_summary$p.table[, "Std. Error"],
      t_value = model_summary$p.table[, "t value"],
      p_value = model_summary$p.table[, "Pr(>|t|)"],
      Type = "Paramétrique"
    )
    
    # Termes smooth
    smooth <- data.frame(
      Term = rownames(model_summary$s.table),
      Estimate = NA,  # Les termes smooth n'ont pas d'estimate direct
      SE = NA,
      t_value = model_summary$s.table[, "F"],  # F-statistic pour les smooth
      p_value = model_summary$s.table[, "p-value"],
      Type = "Smooth"
    )
    
    # Combinaison des résultats
    results <- rbind(parametric, smooth)
    
    # Ajout de la significativité
    results$Significance <- case_when(
      results$p_value < 0.001 ~ "***",
      results$p_value < 0.01 ~ "**",
      results$p_value < 0.05 ~ "*",
      results$p_value < 0.1 ~ ".",
      TRUE ~ ""
    )
    
    return(results)
  }
  
  # Extraction des résultats
  results_table <- extract_gam_results(model_summary)
  
  # Formatage de la table pour présentation
  formatted_table <- results_table %>%
    mutate(
      Term = case_when(
        Term == "(Intercept)" ~ "Intercept",
        Term == "charge_bin1" ~ "Charge binaire (Oui)",
        grepl("^s\\(", Term) ~ gsub("s\\(([^)]+)\\)", "\\1 (smooth)", Term),
        grepl("^ti\\(", Term) ~ gsub("ti\\(([^)]+)\\)", "\\1 (interaction)", Term),
        TRUE ~ Term
      ),
      Estimate = ifelse(is.na(Estimate), "—", sprintf("%.4f", Estimate)),
      SE = ifelse(is.na(SE), "—", sprintf("%.4f", SE)),
      Statistic = ifelse(Type == "Paramétrique", 
                         sprintf("%.3f", t_value), 
                         sprintf("%.3f", t_value)),
      `P-value` = case_when(
        p_value < 0.001 ~ "< 0.001",
        p_value < 0.01 ~ sprintf("%.3f", p_value),
        TRUE ~ sprintf("%.3f", p_value)
      )
    ) %>%
    select(Term, Type, Estimate, SE, Statistic, `P-value`, Significance)
  
  # Création de la table avec kable et kableExtra
  beautiful_table <- formatted_table %>%
    kbl(
      caption = "Résultats du modèle GAM : log(AUCg)",
      col.names = c("Terme", "Type", "Estimation", "Erreur Standard", 
                    "Statistique", "Valeur-p", "Signif."),
      align = c("l", "c", "r", "r", "r", "r", "c")
    ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE,
      position = "center",
      font_size = 12
    ) %>%
    row_spec(0, bold = TRUE, background = "#34495e", color = "white") %>%
    column_spec(1, bold = TRUE, width = "3cm") %>%
    column_spec(2, width = "2cm") %>%
    column_spec(c(3:6), width = "1.5cm") %>%
    column_spec(7, width = "1cm") %>%
    pack_rows("Termes paramétriques", 1, sum(results_table$Type == "Paramétrique"), 
              label_row_css = "background-color: #3498db; color: white;") %>%
    pack_rows("Termes smooth", sum(results_table$Type == "Paramétrique") + 1, nrow(results_table),
              label_row_css = "background-color: #e74c3c; color: white;") %>%
    footnote(
      general = c("Codes de significativité : 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
                  paste("R² ajusté =", round(model_summary$r.sq, 3)),
                  paste("R² expliqué =", round(model_summary$dev.expl, 3)),
                  paste("AIC =", round(AIC(m_ref_plus), 1)),
                  paste("n =", nobs(m_ref_plus))),
      general_title = "Notes :",
      footnote_as_chunk = TRUE
    )
  
  # Affichage de la table
  print(beautiful_table)
  
  # === MÉTHODES D'EXPORT PNG (plusieurs options) ===
  
  # OPTION 1 : Export avec gt (recommandé)
  gt_table <- formatted_table %>%
    gt(caption = "Table 1. Résultats du modèle additif généralisé pour log(AUCg)") %>%
    tab_header(
      title = "Modèle GAM : log(AUCg)"
    ) %>%
    cols_label(
      Term = "Terme",
      Type = "Type", 
      Estimate = "Estimation",
      SE = "Erreur Standard",
      Statistic = "Statistique",
      `P-value` = "Valeur-p",
      Significance = "Signif."
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    tab_style(
      style = cell_borders(
        sides = c("top", "bottom"),
        color = "black",
        weight = px(2)
      ),
      locations = cells_column_labels()
    ) %>%
    tab_footnote(
      footnote = "Codes de significativité : 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"
    ) %>%
    tab_footnote(
      footnote = paste("R² ajusté =", round(model_summary$r.sq, 3))
    ) %>%
    tab_footnote(
      footnote = paste("AIC =", round(AIC(m_ref_plus), 1))
    ) %>%
    opt_table_font(font = "Times New Roman") %>%
    tab_options(
      table.background.color = "white",
      heading.background.color = "white",
      column_labels.background.color = "white"
    )
  
  # Export avec gt
  tryCatch({
    gtsave(gt_table, "gam_results_table_gt.png", vwidth = 1000, vheight = 800)
    cat("Table exportée avec gt : gam_results_table_gt.png\n")
  }, error = function(e) {
    cat("Erreur avec gt :", e$message, "\n")
  })
  
  # OPTION 2 : Export HTML simple (toujours fonctionnel)
  html_output <- as.character(beautiful_table)
  writeLines(html_output, "gam_results_table.html")
  cat("Table HTML sauvegardée : gam_results_table.html\n")
  
  # OPTION 3 : Essai avec webshot (version classique)
  if(require(webshot, quietly = TRUE)) {
    tryCatch({
      # Installation de phantomjs si nécessaire
      if(!webshot:::find_phantom()) {
        webshot::install_phantomjs()
      }
      
      # Export
      temp_html <- "temp_table.html"
      save_kable(beautiful_table, file = temp_html)
      webshot(temp_html, "gam_results_table_webshot.png", 
              vwidth = 1000, vheight = 800, zoom = 2)
      unlink(temp_html)
      cat("Table exportée avec webshot : gam_results_table_webshot.png\n")
    }, error = function(e) {
      cat("Erreur avec webshot :", e$message, "\n")
    })
  }
  
  # OPTION 4 : Instructions manuelles pour l'export
  cat("\n=== INSTRUCTIONS D'EXPORT MANUEL ===\n")
  cat("Si les exports automatiques ne fonctionnent pas :\n")
  cat("1. Ouvrez le fichier 'gam_results_table.html' dans votre navigateur\n")
  cat("2. Faites clic droit > 'Imprimer' > 'Enregistrer au format PDF'\n")
  cat("3. Ou utilisez l'outil de capture d'écran de votre système\n\n")
  
  # Informations supplémentaires sur le modèle
  cat("\nInformations supplémentaires du modèle :\n")
  cat("Deviance expliquée :", round(model_summary$dev.expl * 100, 2), "%\n")
  cat("R² ajusté :", round(model_summary$r.sq, 3), "\n")
  cat("AIC :", round(AIC(m_ref_plus), 1), "\n")
  cat("Nombre d'observations :", nobs(m_ref_plus), "\n")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Chargement des librairies nécessaires
  library(mgcv)
  library(knitr)
  library(kableExtra)
  library(dplyr)
  library(gt)
  
  # Votre modèle GAM
  m_ref_plus <- bam(
    log(AUCg) ~ s(dah, k=10, bs="ts") + 
      s(SMOD_2023, k=10, bs="ts") + 
      ti(dah, SMOD_2023, k=c(15,15)) + 
      s(alpage, bs="re") + 
      charge_bin + 
      s(Charge_log, k=10, bs="ts", by=charge_bin) + 
      ti(Charge_log, SMOD_2023, k=c(10,10), by=charge_bin) + 
      s(alpage, bs="re", by=charge_bin) + 
      ti(Charge_log, alpage, bs=c("tp","re"), by=charge_bin, k=c(6, NA)),
    data = data_mod, 
    method = "fREML", 
    select = TRUE, 
    discrete = TRUE, 
    nthreads = n_th
  )
  
  # Extraction des résultats du modèle
  model_summary <- summary(m_ref_plus)
  
  # Extraction et formatage simplifiés
  extract_gam_results <- function(model_summary) {
    # Termes paramétriques
    parametric <- data.frame(
      Term = rownames(model_summary$p.table),
      Estimate = model_summary$p.table[, "Estimate"],
      SE = model_summary$p.table[, "Std. Error"],
      p_value = model_summary$p.table[, "Pr(>|t|)"],
      Type = "Paramétrique"
    )
    
    # Termes smooth
    smooth <- data.frame(
      Term = rownames(model_summary$s.table),
      Estimate = NA,
      SE = NA,
      p_value = model_summary$s.table[, "p-value"],
      Type = "Smooth"
    )
    
    # Combinaison
    results <- rbind(parametric, smooth)
    
    # Significativité
    results$Significance <- case_when(
      results$p_value < 0.001 ~ "***",
      results$p_value < 0.01 ~ "**",
      results$p_value < 0.05 ~ "*",
      results$p_value < 0.1 ~ ".",
      TRUE ~ ""
    )
    
    return(results)
  }
  
  # Extraction des résultats
  results_table <- extract_gam_results(model_summary)
  
  # Table ultra-simplifiée type article
  formatted_table <- results_table %>%
    mutate(
      Term = case_when(
        Term == "(Intercept)" ~ "Intercept",
        Term == "charge_bin1" ~ "Charge (binaire)",
        Term == "s(dah)" ~ "DAH",
        Term == "s(SMOD_2023)" ~ "SMOD 2023",
        Term == "ti(dah,SMOD_2023)" ~ "DAH × SMOD 2023",
        Term == "s(alpage)" ~ "Alpage",
        Term == "s(Charge_log):charge_bin" ~ "Charge log × Charge binaire",
        Term == "ti(Charge_log,SMOD_2023):charge_bin" ~ "Charge log × SMOD 2023 × Charge binaire",
        Term == "s(alpage):charge_bin" ~ "Alpage × Charge binaire", 
        Term == "ti(Charge_log,alpage):charge_bin" ~ "Charge log × Alpage × Charge binaire",
        TRUE ~ Term
      ),
      Coefficient = ifelse(is.na(Estimate), "—", sprintf("%.3f", Estimate)),
      SE = ifelse(is.na(SE), "—", sprintf("%.3f", SE)),
      `P-value` = case_when(
        p_value < 0.001 ~ "< 0.001",
        TRUE ~ sprintf("%.3f", p_value)
      )
    ) %>%
    select(Term, Coefficient, SE, `P-value`, Significance)
  
  # Table SANS couleurs, SANS type, SANS statistique
  beautiful_table <- formatted_table %>%
    kbl(
      col.names = c("", "Coefficient", "SE", "p", ""),
      align = c("l", "r", "r", "r", "c"),
      format = "html",
      table.attr = 'style="border-collapse: collapse; font-family: serif;"'
    ) %>%
    kable_styling(
      bootstrap_options = "none",
      full_width = FALSE,
      position = "center",
      font_size = 11
    ) %>%
    # AUCUNE couleur - juste des bordures noires
    row_spec(0, bold = TRUE, 
             extra_css = "border-top: 2px solid black; border-bottom: 1px solid black; background: white;") %>%
    column_spec(1, italic = TRUE, width = "6cm") %>%
    column_spec(2:4, width = "1.5cm") %>%
    column_spec(5, width = "0.8cm", bold = TRUE) %>%
    row_spec(nrow(formatted_table), 
             extra_css = "border-bottom: 2px solid black;") %>%
    footnote(
      general = c("*** p < 0.001, ** p < 0.01, * p < 0.05",
                  paste0("R²adj = ", round(model_summary$r.sq, 3), 
                         ", déviance expliquée = ", round(model_summary$dev.expl * 100, 1), "%"),
                  paste0("AIC = ", round(AIC(m_ref_plus), 0), ", n = ", nobs(m_ref_plus))),
      footnote_as_chunk = TRUE
    )
  
  # Affichage
  print(beautiful_table)
  
  # Export HTML simple
  html_output <- as.character(beautiful_table)
  writeLines(html_output, "table_gam_simple.html")
  cat("Table HTML sauvegardée : table_gam_simple.html\n")
  
  # Export avec gt (version minimaliste)
  tryCatch({
    gt_table <- formatted_table %>%
      gt() %>%
      cols_label(
        Term = "",
        Coefficient = "Coefficient",
        SE = "SE", 
        `P-value` = "p",
        Significance = ""
      ) %>%
      tab_style(
        style = list(
          cell_text(weight = "bold"),
          cell_borders(sides = "top", color = "black", weight = px(2)),
          cell_borders(sides = "bottom", color = "black", weight = px(1))
        ),
        locations = cells_column_labels()
      ) %>%
      tab_style(
        style = cell_text(style = "italic"),
        locations = cells_body(columns = Term)
      ) %>%
      tab_style(
        style = cell_borders(sides = "bottom", color = "black", weight = px(2)),
        locations = cells_body(rows = nrow(formatted_table))
      ) %>%
      tab_source_note("*** p < 0.001, ** p < 0.01, * p < 0.05") %>%
      tab_source_note(paste0("R²adj = ", round(model_summary$r.sq, 3), 
                             ", déviance expliquée = ", round(model_summary$dev.expl * 100, 1), "%")) %>%
      tab_source_note(paste0("AIC = ", round(AIC(m_ref_plus), 0), ", n = ", nobs(m_ref_plus))) %>%
      opt_table_font(font = "serif") %>%
      tab_options(table.background.color = "white", table.font.size = 11)
    
    gtsave(gt_table, "table_gam_gt.png", vwidth = 800, vheight = 600)
    cat("Table GT exportée : table_gam_gt.png\n")
  }, error = function(e) {
    cat("Export GT non disponible\n")
  })
  
  cat("\nTable simplifiée créée : AUCUNE couleur, pas de colonnes Type/Statistique\n")