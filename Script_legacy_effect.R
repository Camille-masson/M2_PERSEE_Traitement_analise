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
        Charge, Charge_log, charge_bin
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
      
      
      m_by_alp <- m_ref_plus
      
      
      
      ## FINAL COMPLET : 
      
      # =========================================================
      # Libraries
      # =========================================================
      library(dplyr)
      library(tidyr)
      library(purrr)
      library(tibble)
      library(ggplot2)
      
      # =========================================================
      # Look & feel
      # =========================================================
      base_family <- "Segoe UI"   # ou "Inter"/"Helvetica" selon ta machine
      
      # Palette des 3 classes (labels EXACTS)
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
        select(alpage, SMOD_2023, Charge_log, AUCg, ch_lo, ch_hi)
      
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
      # Export net (PNG)
      # =========================================================
      if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png("Figure_5_Resultat_M2_Predict_model.png",
                      width = 12.5, height = 6, units = "in", res = 450, scaling = 1)
        print(p); dev.off()
      } else {
        ggsave("Figure_5_Resultat_M2_Predict_model.png",
               p, width = 12.5, height = 6, dpi = 450, device = cairo_png)
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
      
      
      
      
      # CARTE DE CHALEUR
      
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
    
    
  }
  
  