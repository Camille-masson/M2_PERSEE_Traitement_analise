## APPROCHE MONO-ALPAGE

## SORTIE
# Dossier de sortie
out_dataset       <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_",alpage))
dir.create(out_dataset, recursive=TRUE, showWarnings=FALSE)
# Un .fst contetant l'ensemble du dataset 
in_rds <- file.path(out_dataset, paste0("dataset_pheno_LONG_", YEAR, "_", alpage, ".rds"))



## Un alpage : détéction des zones maximisant le nombre de pixel dans les classe de chargement avec les contraintes
## Habitat FSCA et DAH

{

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

## ————————————————————————————————————————————————————————————
## PARAMÈTRES À MODIFIER SEULEMENT ICI
YEAR      <- 2022
alpage    <- "Cayolle"           # ← un seul alpage à la fois
root_dir  <- getwd()
output_dir<- file.path(root_dir, "outputs")

hab_keep  <- c("Nardaies denses du subalpin",
               "Queyrellins",
               "Megaphorbiaies et Aulnaies",
               "P. productives",
               "P. nivales",
               "P. thermiques écorchées")
## ————————————————————————————————————————————————————————————

## 1. Lecture + pré‑filtre complet -------------------------------
df <- readRDS(in_rds) %>% 
  filter(habitat_label %in% hab_keep) %>%
  filter(!is.na(IRG) & !is.na(NDVI) &        # plus de NA IRG/NDVI
           !is.na(dah) & !is.na(fsca))        # plus de NA dah/fsca

## 2. opt_window robuste aux NA ---------------------------------
opt_window <- function(dat, dah_col = "dah", fsca_col = "fsca",
                       win = 0.10, grid_probs = seq(.05, .95, .05)) {
  
  dah_cent  <- quantile(dat[[dah_col]],  grid_probs, na.rm = TRUE)
  fsca_cent <- quantile(dat[[fsca_col]], grid_probs, na.rm = TRUE)
  
  best <- list(total = -1)
  
  for (dc in dah_cent) {
    dah_ok <- dat[[dah_col]] >= (1-win)*dc & dat[[dah_col]] <= (1+win)*dc
    if (!any(dah_ok, na.rm = TRUE)) next
    
    for (fc in fsca_cent) {
      mask <- dah_ok &
        dat[[fsca_col]] >= (1-win)*fc & dat[[fsca_col]] <= (1+win)*fc
      if (!any(mask, na.rm = TRUE)) next        # ⬅️ na.rm = TRUE
      
      counts <- table(dat$load_class[mask])
      if (length(counts) < 3) next
      
      total <- sum(counts)
      if (total > best$total) {
        best <- list(
          total  = total,
          mask   = mask %in% TRUE,              # remplace NA par FALSE
          center = c(dah = dc, fsca = fc),
          counts = counts
        )
      }
    }
  }
  best
}


## 3. Boucle habitat : terciles + fenêtre -------------------------------------
res_list <- vector("list", length(hab_keep))
names(res_list) <- hab_keep
summary_tab <- list()

for (hab in hab_keep) {
  sub <- df %>% filter(habitat_label == hab)
  
  ## 3a. Découpe en terciles
  q1 <- quantile(sub$Charge, .2, na.rm = TRUE)
  q2 <- quantile(sub$Charge, .8, na.rm = TRUE)
  
  sub <- sub %>%
    mutate(load_class = case_when(
      Charge <= q1           ~ "Faible",
      Charge <= q2           ~ "Moyen",
      TRUE                   ~ "Fort"
    ))
  
  ## 3b. Recherche de la fenêtre ±10 %
  best <- opt_window(sub)
  
  if (best$total < 0) {                    # aucun couple valide
    message("⚠️  Habitat «", hab, "» : impossible de trouver une fenêtre ",
            "contenant les 3 classes (terciles).")
    next
  }
  
  res_list[[hab]] <- sub[best$mask, ]      # pixels retenus
  
  summary_tab[[hab]] <- data.frame(
    habitat       = hab,
    Charge_q33    = as.numeric(q1),
    Charge_q66    = as.numeric(q2),
    dah_center    = best$center["dah"],
    fsca_center   = best$center["fsca"],
    n_total       = best$total,
    n_faible      = best$counts["Faible"],
    n_moyen       = best$counts["Moyen"],
    n_fort        = best$counts["Fort"],
    stringsAsFactors = FALSE
  )
}

summary_tab <- bind_rows(summary_tab)

cat("\n=== Seuils de charge (terciles) & effectifs par habitat ===\n")
print(summary_tab, row.names = FALSE)

## 4. Concaténation des pixels retenus + sauvegarde ---------------------------
df_out <- bind_rows(res_list)

out_dir <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_", alpage))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(out_dir,
                     sprintf("dataset_pheno_LONG_%d_%s_TERCILES.rds", YEAR, alpage))
saveRDS(df_out, out_rds)
message("\nSauvegardé : ", out_rds)


}



## La MEME CHOSE sans les classes de chargement juste Habitat DAH et FSCA les centre sont basé sur la median de DAH et FSCA de chaque habitat
alpage = "Alpe-Sud"
YEAR = 2023



## DATASET Alpe-Sud qui regroupe : Cayolle, Viso, Sanguiniere

library(dplyr)
library(purrr)   # map_dfr()

build_ndvi_dataset <- function(year,
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

## Exemple d’appel -------------------------------------------------------------
alpages <- c("Viso", "Cayolle", "Sanguiniere")
build_ndvi_dataset(year = 2023, alpages = alpages)




## REGROUPEMENT DES HABITATS PAR MAXV POUR LES HAB MINORITAIRE 

out_dir <- file.path(output_dir, "10. NDVI Effect",
                     paste0("dataset_", alpage))

dt_alpesud = file.path(out_dir, paste0("dataset_pheno_LONG_", YEAR, "_", alpage, ".rds"))




# ── 5. Boxplot : MAXV par habitat (jeu filtré) ───────────────────────────────


df_plot <- readRDS(dt_alpesud) %>% 
  filter(!is.na(habitat_label),          # habitats renseignés
         !is.na(MAXV)) %>%              # MAXV renseigné
  group_by(cell) %>%                    # dé‑doublonnage pixel
  slice(1) %>%                          #    → garde la 1re ligne de chaque pixel
  ungroup()

ggplot(df_plot,
       aes(x = reorder(habitat_label, MAXV, median, na.rm = TRUE),
           y = MAXV)) +
  geom_boxplot(outlier.alpha = 0.25, na.rm = TRUE) +
  coord_flip() +
  labs(x = NULL, y = "MAXV",
       title = "Distribution de MAXV par habitat – Alpe‑Sud 2022 (un seul pixel par DOY)") +
  theme_bw(base_size = 11)


## On regroupe ceratin habitat avec les habitat faiblement représenté et qui ont un MAXV proche
## On arrive alors a 8 habitats :
## - P. productives
## - Megaphorbiaies et Aulnaies
## - Landes
## - Queyrellins
## - P. nitrophiles + P. humides = P.humides-nitrophiles
## - Nardaies denses du subalpin
## - P. thermiques écorchées
## - P. nivales + P. en bombement de lalpin + P. th. méditerranéo-montagnardes + P. intermédiaires de lalpin = P. alpine

## Habitat supprimé : 
## - Formations minérales
## - Forêts non pastorales




# -----------------------------------------------------------------------------
#  Fenêtre commune (par habitat)  ±10 % FSCA  &  ±0,15 DAH
#  -> conserve un maximum de pixels pour **tous les habitats présents** (≠ NA)
# -----------------------------------------------------------------------------
#  • Lit le jeu complet 2022 – Alpe‑Sud.
#  • Sélectionne, pour chaque habitat, la fenêtre la plus « dense » en pixels
#    autour d’un couple (fsca₀, dah₀) optimisé.
#  • Sauvegarde le sous‑ensemble concaténé + un tableau récapitulatif.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
})

YEAR      <- 2022
alpage    <- "Alpe-Sud"
root_dir  <- getwd()
output_dir<- file.path(root_dir, "outputs")

in_rds <- file.path(output_dir, "10. NDVI Effect",
                    paste0("dataset_", alpage),
                    sprintf("dataset_pheno_LONG_%d_%s.rds", YEAR, alpage))

WIN_FSCA <- 0.1          # ±10 % FSCA
WIN_DAH  <- 0.15          # ±0,15 DAH
GRID_P   <- seq(0.05, 0.95, 0.05)

# ── 1. Lecture, dé‑doublonnage pixel & regroupements explicites ──────────────
library(dplyr)   # assure‑toi qu’il est bien chargé en dernier

# ── 1. Lecture, dé‑doublonnage pixel & regroupements ────────────────────────
df0 <- readRDS(in_rds) %>% 
  filter(!is.na(habitat_label), !is.na(dah), !is.na(fsca)) %>% 
  group_by(cell) %>% slice(1) %>% ungroup() %>% 
  
  mutate(
    habitat_label = case_when(
      # ── Regroupement 1 : P. humides‑nitrophiles ───────────────
      habitat_label %in% c("P. nitrophiles",
                           "P. humides")                                 ~ "P.humides-nitrophiles",
      
      # ── Regroupement 2 : P. alpine ────────────────────────────
      habitat_label %in% c("P. nivales",
                           "P. en bombement de l\u0092alpin",
                           "P. en bombement de l'alpin",
                           "P. th. méditerranéo-montagnardes",
                           "P. th. mediterraneo-montagnardes",
                           "P. intermédiaires de l\u0092alpin",
                           "P. intermédiaires de l'alpin")              ~ "P. alpine",
      
      TRUE ~ habitat_label   # tout le reste inchangé
    )
  ) %>% 
  
  filter(!habitat_label %in% c("Formations minérales",
                               "Forêts non pastorales"))

hab_list <- sort(unique(df0$habitat_label))   # doit en rester 8

# ── 2. Fonction : meilleure fenêtre pour un habitat ──────────────────────────
best_window <- function(sub) {
  med_f <- quantile(sub$fsca, GRID_P, na.rm = TRUE)
  med_d <- quantile(sub$dah , GRID_P, na.rm = TRUE)
  
  best <- list(n = 0)
  for (fc in med_f)
    for (dc in med_d) {
      m <- abs(sub$fsca - fc) <= WIN_FSCA * fc &
        abs(sub$dah  - dc) <= WIN_DAH
      n_pix <- sum(m)
      if (n_pix > best$n)
        best <- list(mask = m, n = n_pix,
                     center = c(dah = dc, fsca = fc))
    }
  best
}

# ── 3. Boucle sur les 8 habitats ─────────────────────────────────────────────
summary_tab <- list()
res_list    <- list()

for (h in hab_list) {
  sub <- filter(df0, habitat_label == h)
  bw  <- best_window(sub)
  res_list[[h]] <- sub[bw$mask, ]
  
  summary_tab[[h]] <- tibble(
    habitat     = h,
    dah_center  = bw$center["dah"],
    fsca_center = bw$center["fsca"],
    n_total     = nrow(sub),
    n_window    = bw$n
  )
}

summary_tab <- bind_rows(summary_tab)

cat("\n=== Pixels retenus (±0,15 DAH & ±10 % FSCA) par habitat ===\n")
print(summary_tab, row.names = FALSE)

# ── 4. Sauvegarde ────────────────────────────────────────────────────────────
out_dir <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_", alpage))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(out_dir,
                     sprintf("dataset_pheno_LONG_%d_%s_WIN10pct.rds", YEAR, alpage))

saveRDS(bind_rows(res_list), out_rds)
write_csv(summary_tab,
          file.path(out_dir, sprintf("summary_window_%d_%s.csv", YEAR, alpage)))

message("✓ Nouveau dataset : ", out_rds)






dt_alpesud_filtered = file.path(out_dir, paste0("dataset_pheno_LONG_", YEAR, "_", alpage, "_WIN10pct.rds"))

## Controle du ONSET10 avec les paramètres DAH et FSCA 
# PARTIE 1 PLOT 

library(hexbin)
df <- readRDS(dt_alpesud_filtered )   
df2 <- df %>%
  group_by(habitat_label) %>%
  mutate(resid = ONSET10 - median(ONSET10, na.rm = TRUE)) %>%
  ungroup()

p_hex <- ggplot(df2, aes(dah, fsca, z = ONSET10))+
  stat_summary_hex(bins = 30, fun = median)+
  scale_fill_distiller(palette = "Spectral", name = "ONSET10")+
  facet_wrap(~ habitat_label, scales = "free")+
  labs(title = "Médiane ONSET10 – Hexbin 2D",
       x = "dah", y = "fSCA")+
  theme_bw(9)

p_violin <- ggplot(df2, aes(habitat_label, resid, fill = habitat_label))+
  geom_violin(trim = FALSE, colour = NA, alpha = .8)+
  geom_hline(yintercept = 0, colour = "grey30")+
  coord_flip()+
  scale_fill_brewer(palette = "Set3", guide = "none")+
  labs(title = "Résidus ONSET10\n(médiane retirée)",
       x = NULL, y = "ONSET10 – médiane")+
  theme_bw(9)

library(patchwork)
p_hex + p_violin + plot_layout(widths = c(3,1))


## FILTRE DES OUTLAYER

library(dplyr)

# ± 2 MAD par habitat  +  comptage des pixels
df_clean <- df2 %>%                                 # df2 = jeu avec resid déjà calculé
  group_by(habitat_label) %>% 
  mutate(
    mad_val = mad(resid, na.rm = TRUE),             # écart médian absolu
    keep    = abs(resid) <= 2 * mad_val             # ← seuil plus strict
  ) %>% 
  ungroup()

# tableau récap : combien de pixels retirés / conservés ?
tbl_out <- df_clean %>% 
  group_by(habitat_label) %>% 
  summarise(
    n_total   = n(),
    n_removed = sum(!keep),
    n_kept    = sum(keep),
    pct_removed = round(100 * n_removed / n_total, 1)
  ) %>% 
  arrange(desc(pct_removed))

print(tbl_out, n = Inf)

# on enlève les outliers pour la suite des graphes
df_clean <- df_clean %>% 
  filter(keep) %>% 
  select(-mad_val, -keep)




df_check <- df_clean %>%                       # calcule à nouveau les résidus
  group_by(habitat_label) %>% 
  mutate(resid = ONSET10 - median(ONSET10, na.rm = TRUE)) %>% 
  ungroup()

p_hex_clean    <- p_hex  %+% df_check          # même objets ggplot, nouvelle data
p_violin_clean <- p_violin %+% df_check

p_hex_clean + p_violin_clean + plot_layout(widths = c(3,1))




saveRDS (df_clean,dt_alpesud_filtered )


## MODELE

dt_alpesud_filtered = file.path(out_dir, paste0("dataset_pheno_LONG_", YEAR, "_", alpage, "_WIN10pct.rds"))
dt_alpesud_nofiltered = file.path(out_dir, paste0("dataset_pheno_LONG_", YEAR, "_", alpage, "_WIN10pct_nofiltered.rds"))

library(dplyr); library(readr)

df_filt <- readRDS(dt_alpesud_filtered ) %>% 
  mutate(version = "filtré")

df_raw  <- readRDS(dt_alpesud_nofiltered) %>% 
  mutate(version = "non filtré")

df_all  <- bind_rows(df_filt, df_raw)

library(dplyr)
library(lme4)
library(broom.mixed)   # tidy()
library(performance)   # r2()

fit_mix <- function(d) lmer(ONSET10 ~ dah + fsca + (1 | alpage), data = d)

mod_mix_tab <- df_all %>% 
  group_by(version, habitat_label) %>% 
  do({
    m  <- fit_mix(.)
    
    ## (a) pentes fixes
    betas <- tidy(m, effects = "fixed")[c("term","estimate")] %>% 
      tidyr::pivot_wider(names_from = term, values_from = estimate)   # β_dah  β_fsca
    
    ## (b) qualité globale
    sig  <- sigma(m)                           # RMSE = σ
    r2   <- performance::r2(m)                 # R2_marg, R2_cond
    
    bind_cols(betas,
              rmse   = sig,
              r2_m   = r2$R2_marginal,
              r2_c   = r2$R2_conditional)
  }) %>% 
  ungroup()

print(mod_mix_tab, n = Inf)




library(dplyr)
library(ggplot2)
library(patchwork)

# On part de ton tableau de résultats mod_mix_tab
# Il contient pour chaque habitat et chaque version :
#   r2_m = R² marginal (dah+fsca)
#   rmse = RMSE intra‑site

# (1) Barplot des R² marginales
p1 <- ggplot(mod_mix_tab, aes(habitat_label, r2_m, fill = version)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Part de variance expliquée par DAH+FSCA (R² marginal)",
    x = NULL, y = "R² marginal"
  ) +
  theme_bw(base_size = 11)

# (2) Barplot des RMSE
p2 <- ggplot(mod_mix_tab, aes(habitat_label, rmse, fill = version)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Précision intra‑site (RMSE en jours)",
    x = NULL, y = "RMSE (jours)"
  ) +
  theme_bw(base_size = 11)

# Composition côte-à-côte
p1 + p2 + 
  plot_layout(ncol = 2) + 
  plot_annotation(
    caption = "Comparaison filtré vs non filtré par habitat"
  )







## JEU DE DONNEE SEMBLE OK 





dt_alpesud_filtered = file.path(out_dir, paste0("dataset_pheno_LONG_", YEAR, "_", alpage, "_WIN10pct.rds"))

library(dplyr)
library(ggplot2)

# calcul du seuil 95e centile global
seuil95 <- quantile(df$Charge, 0.95, na.rm = TRUE)

df_win <- df %>% 
  mutate(charg_win = pmin(Charge, seuil95))  # cap au 95e

ggplot(df_win, aes(x = charg_win, y = MAXV)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = FALSE, colour = "darkblue") +
  facet_wrap(~ habitat_label, scales = "free") +
  labs(
    title = "MAXV vs chargement (winsorisé au 95ᵉ centile)",
    x     = paste0("chargement (capé à ", round(seuil95), ")"),
    y     = "MAXV"
  ) +
  theme_bw(11)



library(dplyr)

df_cat <- df %>% 
  group_by(habitat_label) %>% 
  mutate(
    # calcule les 1ᵉʳ et 2ᵉᵐᵉ tertiles de Charge
    q33 = quantile(Charge, probs = 1/3, na.rm = TRUE),
    q66 = quantile(Charge, probs = 2/3, na.rm = TRUE),
    cat_ch = case_when(
      Charge <= q33 ~ "faible",
      Charge <= q66 ~ "moyen",
      TRUE          ~ "fort"
    )
  ) %>% 
  ungroup() %>% 
  select(-q33, -q66)




library(dplyr)
library(ggplot2)

# 1. Définir q20 et q80 *par habitat*
df_cat2 <- df %>% 
  group_by(habitat_label) %>% 
  mutate(
    q20 = quantile(Charge, 0.20, na.rm = TRUE),
    q80 = quantile(Charge, 0.80, na.rm = TRUE),
    cat_ch = case_when(
      Charge <= q20           ~ "faible",
      Charge >  q80           ~ "fort",
      TRUE                    ~ "moyen"
    ),
    # transformer en facteur pour contrôler l’ordre
    cat_ch = factor(cat_ch, levels = c("faible","moyen","fort"))
  ) %>% 
  ungroup() %>% 
  select(-q20, -q80)

# 2. Calculer les effectifs pour l’annotation
counts2 <- df_cat2 %>% 
  count(habitat_label, cat_ch)

# 3. Tracer les boxplots avec les effectifs
ggplot(df_cat2, aes(x = cat_ch, y = MAXV, fill = cat_ch)) +
  geom_boxplot(outlier.shape = NA, width = .6) +
  geom_jitter(width = .15, alpha = .3, size = .6, colour = "grey20") +
  geom_text(
    data = counts2,
    aes(x = cat_ch, y = Inf, label = n),
    vjust = 1.2, size = 3
  ) +
  scale_fill_brewer(palette = "Set2", name = "Classe") +
  facet_wrap(~ habitat_label, scales = "free_y") +
  labs(
    title    = "MAXV par classe de chargement (20/60/20)",
    subtitle = "Classes définies par percentiles 20 % et 80 %",
    x        = "Classe de chargement",
    y        = "MAXV"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(margin = margin(t = 5))
  )

library(lmerTest)    # pour lmerTest::lmer et p‑values
library(performance) # pour r2()

df <- readRDS(dt_alpesud_filtered)

# 1) Landes
mod_Landes <- lmerTest::lmer(
  MAXV ~ Charge + dah + fsca + (1|alpage),
  data   = df,
  subset = habitat_label == "Landes"
)
summary(mod_Landes)
performance::r2(mod_Landes)

# 2) Megaphorbiaies et Aulnaies
mod_Mega <- lmerTest::lmer(
  MAXV ~ Charge + dah + fsca + (1|alpage),
  data   = df,
  subset = habitat_label == "Megaphorbiaies et Aulnaies"
)
summary(mod_Mega)
performance::r2(mod_Mega)

# 3) Nardaies denses du subalpin
mod_Nardaies <- lmerTest::lmer(
  MAXV ~ Charge + dah + fsca + (1|alpage),
  data   = df,
  subset = habitat_label == "Nardaies denses du subalpin"
)
summary(mod_Nardaies)
performance::r2(mod_Nardaies)

# 4) P. alpine
mod_Palpine <- lmerTest::lmer(
  MAXV ~ Charge + dah + fsca + (1|alpage),
  data   = df,
  subset = habitat_label == "P. alpine"
)
summary(mod_Palpine)
performance::r2(mod_Palpine)

# 5) P. productives
mod_Pprod <- lmerTest::lmer(
  MAXV ~ Charge + dah + fsca + (1|alpage),
  data   = df,
  subset = habitat_label == "P. productives"
)
summary(mod_Pprod)
performance::r2(mod_Pprod)

# 6) P. thermiques écorchées
mod_Ptherm <- lmerTest::lmer(
  MAXV ~ Charge + dah + fsca + (1|alpage),
  data   = df,
  subset = habitat_label == "P. thermiques écorchées"
)
summary(mod_Ptherm)
performance::r2(mod_Ptherm)

# 7) P.humides-nitrophiles
mod_PHum <- lmerTest::lmer(
  MAXV ~ Charge + dah + fsca + (1|alpage),
  data   = df,
  subset = habitat_label == "P.humides-nitrophiles"
)
summary(mod_PHum)
performance::r2(mod_PHum)

# 8) Queyrellins
mod_Quey <- lmerTest::lmer(
  MAXV ~ Charge + dah + fsca + (1|alpage),
  data   = df,
  subset = habitat_label == "Queyrellins"
)
summary(mod_Quey)
performance::r2(mod_Quey)













# ─────────────────────────────────────────────
# 0. Packages
suppressPackageStartupMessages(library(dplyr))

# ─────────────────────────────────────────────
# 1.Fichiers d’E/S
YEAR   <- 2022
alpage <- "Alpe-Sud"

out_dir <- file.path(getwd(), "outputs/10. NDVI Effect", paste0("dataset_", alpage))

dt_alpesud_filtered <- file.path(out_dir,paste0("dataset_pheno_LONG_", YEAR, "_",alpage, "_WIN10pct.rds"))     # jeu filtré (1 DOY/pixel)

in_rds <- file.path(out_dir,paste0("dataset_pheno_LONG_", YEAR, "_",alpage, ".rds"))              # jeu LONG brut

dt_alpesud_full <-  file.path(out_dir,paste0("dataset_pheno_LONG_", YEAR, "_",alpage, "_filtered_full.rds")) 

# ─────────────────────────────────────────────
# 2. Jeu filtré : cellules retenues **et** habitats corrects
df_filt <- readRDS(dt_alpesud_filtered)

cells_ok <- unique(df_filt$cell)                       # liste des pixels utiles
hab_corr <- df_filt %>%                               # table de correspondance (habitat corrigé)
  select(cell, habitat_label) %>% distinct()

# ─────────────────────────────────────────────
# 3. Jeu LONG complet : on garde seulement les pixels voulus
df_full <- readRDS(in_rds) %>% 
  filter(cell %in% cells_ok) %>%                      # 3.1 filtrage sur cellID
  select(-habitat_label) %>%                          # 3.2 on retire l’ancien habitat
  left_join(hab_corr, by = "cell") %>%                # 3.3 on injecte le bon
  relocate(habitat_label, .after = cell)              # (option) remet la colonne au même endroit

# ─────────────────────────────────────────────
# 4. Sauvegarde
saveRDS(df_full, dt_alpesud_full)
message("✓ Jeu complet filtré enregistré dans : ", dt_alpesud_full)




# ─────────────────────────────────────────────
# 0. Packages
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ─────────────────────────────────────────────
# 1. Paramètres
YEAR   <- 2022
alpage <- "Alpe-Sud"

out_dir <- file.path(getwd(), "outputs/10. NDVI Effect", paste0("dataset_", alpage))

# relecture
dt_alpesud_filtered_full <- readRDS(file.path(out_dir,paste0("dataset_pheno_LONG_", YEAR, "_",alpage, "_filtered_full.rds")))

# ─────────────────────────────────────────────
# 2. Classe de chargement (terciles)  – suppression des NA
df <- dt_alpesud_filtered_full %>% 
  mutate(
    load_class = cut(
      Charge,
      breaks = c(-Inf,
                 quantile(Charge, .20, na.rm = TRUE),
                 quantile(Charge, .80, na.rm = TRUE),
                 Inf),
      labels = c("Faible", "Moyen", "Fort")
    ),
    load_class = factor(load_class, levels = c("Faible", "Moyen", "Fort"))
  ) %>% 
  filter(!is.na(load_class))

# ─────────────────────────────────────────────
# 3. Palette & libellés
pal_trt <- c(Faible = "khaki3", Moyen = "seagreen3", Fort = "steelblue")
lab_trt <- c(Faible = "low grazing", Moyen = "medium grazing", Fort = "high grazing")

# ─────────────────────────────────────────────
# 4. Boucle par habitat
habitats <- sort(unique(df$habitat_label))
plots <- vector("list", length(habitats)); names(plots) <- habitats

for (hab in habitats) {
  
  df_h <- filter(df, habitat_label == hab)
  
  # ── 4.1 NDVI : médiane + bande 20–80 %
  ndvi_stats <- df_h %>% 
    group_by(DOY, load_class) %>% 
    summarise(
      med_ndvi = median(NDVI, na.rm = TRUE),
      low80    = quantile(NDVI, .20, na.rm = TRUE),
      hi80     = quantile(NDVI, .80, na.rm = TRUE),
      .groups  = "drop"
    ) %>% 
    group_by(load_class) %>% 
    filter(any(med_ndvi > 0)) %>% 
    ungroup()
  
  # ── 4.2 🔄 IRG : pic (jour + valeur) par classe
  irg_peaks <- df_h %>% 
    group_by(DOY, load_class) %>% 
    summarise(med_irg = median(IRG, na.rm = TRUE), .groups = "drop") %>% 
    filter(load_class %in% ndvi_stats$load_class) %>% 
    group_by(load_class) %>% 
    slice_max(med_irg, n = 1, with_ties = FALSE) %>% 
    ungroup() %>% 
    rename(DOY_peak = DOY, irg_max = med_irg)
  
  # ── 4.3 🔄 Bornes + facteur d’échelle pour l’axe IRG
  y_min <- min(ndvi_stats$low80, na.rm = TRUE)
  y_max <- max(ndvi_stats$hi80 , na.rm = TRUE)
  pad   <- 0.10 * (y_max - y_min)
  
  baseline_y <- y_min
  peak_y     <- baseline_y + 0.90 * (y_max - y_min)      # 90 % de la hauteur NDVI
  sf         <- (peak_y - baseline_y) / max(irg_peaks$irg_max)
  y_max_all  <- max(y_max, peak_y) + 0.30 * pad           # marge en haut
  
  # ── 4.4 Palette réduite aux classes présentes
  cls_here  <- unique(ndvi_stats$load_class)
  cols_here <- pal_trt[cls_here]
  labs_here <- lab_trt[cls_here]
  
  # ── 4.5 Graphique
  p <- ggplot() +
    # NDVI : ruban + médiane
    geom_ribbon(
      data = ndvi_stats,
      aes(DOY, ymin = low80, ymax = hi80,
          fill = load_class), alpha = .25, colour = NA) +
    geom_line(
      data = ndvi_stats,
      aes(DOY, med_ndvi, colour = load_class), size = 1) +
    
    # 🔄 Barre verticale au DOY du pic IRG
    geom_segment(
      data = irg_peaks,
      aes(x = DOY_peak, xend = DOY_peak,
          y = baseline_y, yend = irg_max * sf,
          colour = load_class), size = 2) +
    
    # 🔄 Trait horizontal à la hauteur du pic IRG
    geom_hline(
      data = irg_peaks,
      aes(yintercept = irg_max * sf, colour = load_class),
      linetype = "dashed", size = .8) +
    
    # Axes
    scale_x_continuous(
      breaks = seq(60, 330, 30),
      minor_breaks = seq(60, 330, 10)) +
    scale_y_continuous(
      name     = "Plant Phenology Index (NDVI)",
      limits   = c(baseline_y, y_max_all),
      expand   = c(0, 0),
      sec.axis = sec_axis(~ . / sf, name = "Green‑up rate (IRG)")) +
    
    # Couleurs & légendes
    scale_colour_manual(values = cols_here, labels = labs_here,
                        name = NULL, na.translate = FALSE) +
    scale_fill_manual(values = cols_here, labels = labs_here,
                      name = NULL, na.translate = FALSE) +
    
    labs(title = sprintf("%s – %s – %d", alpage, hab, YEAR),
         x = "Julian date") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_line(size = .2, linetype = "dotted"),
      panel.grid.major = element_line(size = .3),
      legend.position  = c(.80, .92),
      legend.background = element_blank()
    )
  
  plots[[hab]] <- p
}

# ─────────────────────────────────────────────
# 5. Affichage
for (p in plots) print(p)














Hill_file <-"C:/Users/massocam/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/raster/Alti/Hill_1_Viso.tif"
MASK_file <- "C:/Users/massocam/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/raster/Alti/MASK.shp"











# ──────────────────────────────
# 0. Packages
suppressPackageStartupMessages(library(terra))

# ──────────────────────────────
# 1. Chemins des fichiers
hill_path <- "C:/Users/massocam/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/raster/Alti/Hill_1_Viso.tif"
mask_path <- "C:/Users/massocam/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/raster/Alti/MASK.shp"

out_path  <- file.path(dirname(hill_path), "hill_viso_crop.tif")

# ──────────────────────────────
# 2. Lecture des données
hill <- rast(hill_path)   # hillshade
mask <- vect(mask_path)   # polygone(s) de découpe

# ──────────────────────────────
# 3. Harmonisation des projections (si nécessaire)
if (!compareGeom(hill, mask, stopOnError = FALSE)) {
  mask <- project(mask, crs(hill))
}

# ──────────────────────────────
# 4. Découpage + masque
hill_crop <- mask(crop(hill, mask), mask)

# ──────────────────────────────
# 5. Sauvegarde
writeRaster(hill_crop, out_path, overwrite = TRUE)
cat("✓ Raster croppé écrit dans :", out_path, "\n")

























                                                                                    ##Model 

# ── 0. Prérequis ────────────────────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(car)          # pour leveneTest() si besoin
# df_raw  ← jeu non filtré   (dt_alpesud_nofiltered)
# df_filt ← jeu filtré       (dt_alpesud_filtered)

# ── 1. Tableau comparatif des variances d’ONSET10 ───────────────────────────
var_summary <- function(x, tag){
  x %>% 
    group_by(habitat_label) %>% 
    summarise(var = var(ONSET10, na.rm = TRUE), .groups = "drop") %>% 
    mutate(dataset = tag)
}

tab_var <- bind_rows(
  var_summary(df_raw , "Non filtré"),
  var_summary(df_filt, "Filtré")
)

bloc_stats <- tab_var %>%                       # ratio  < 1 ⇒ variance réduite
  tidyr::pivot_wider(names_from = dataset, values_from = var) %>% 
  mutate(ratio = `Filtré` / `Non filtré`)

print(bloc_stats)
#> # A tibble: 6 × 4
#>   habitat_label                `Non filtré`  `Filtré` ratio
#>   <chr>                                <dbl>     <dbl> <dbl>
#> 1 Nardaies denses du subalpin           …         …     …
#> …

# ── 2. Test de réduction de variance (F-test) ───────────────────────────────
test_var <- purrr::map_dfr(unique(df_raw$habitat_label), \(h){
  var.test(df_raw %>% filter(habitat_label==h) %>% pull(ONSET10),
           df_filt %>% filter(habitat_label==h) %>% pull(ONSET10)) |>
    {\(o) tibble(habitat = h, F = o$statistic, pval = o$p.value)}()
})
print(test_var)

# ── 3. Figure barres côte-à-côte ────────────────────────────────────────────
ggplot(tab_var,
       aes(x = habitat_label, y = var, fill = dataset)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = NULL, y = "Variance d’ONSET10",
       fill = NULL,
       title = "Réduction de la variabilité d’ONSET10 après filtrage") +
  theme_bw(base_size = 11)






ggplot(bind_rows(df_raw  %>% mutate(set="Non filtré"),
                 df_filt %>% mutate(set="Filtré")),
       aes(x = ONSET10, fill = set)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ habitat_label, ncol = 2)









## INDIVIDUEL 




# GLM (gaussian) – Effet de FSCA & DAH sur ONSET10
# Version « explicite » : un modèle écrit à la main pour chaque habitat
# (jeu non filtré puis jeu filtré)
# -----------------------------------------------------------------------------
# Prérequis : dplyr + base glm ; aucune boucle ni purrr
# -----------------------------------------------------------------------------

library(dplyr)

YEAR   <- 2022
alpage <- "Alpe-Sud"
root_dir  <- getwd()
output_dir<- file.path(root_dir, "outputs", "10. NDVI Effect",
                       paste0("dataset_", alpage))

dt_alpesud_filtered   <- file.path(output_dir,
                                   sprintf("dataset_pheno_LONG_%d_%s_WIN10pct.rds",
                                           YEAR, alpage))
dt_alpesud_nofiltered <- file.path(output_dir,
                                   sprintf("dataset_pheno_LONG_%d_%s_WIN10pct_nofiltered.rds",
                                           YEAR, alpage))

hab_keep <- c("Nardaies denses du subalpin", "Queyrellins",
              "Megaphorbiaies et Aulnaies", "P. productives",
              "P. nivales", "P. thermiques écorchées")

# --- Chargement --------------------------------------------------------------
df_raw  <- readRDS(dt_alpesud_nofiltered) %>% filter(habitat_label %in% hab_keep)
df_filt <- readRDS(dt_alpesud_filtered)    %>% filter(habitat_label %in% hab_keep)

# =============================================================================
# 1. Modèles GLM – jeu NON FILTRÉ
# =============================================================================

### 1.1 Nardaies denses du subalpin ------------------------------------------
d_nardaies_nf <- df_raw %>% filter(habitat_label == "Nardaies denses du subalpin")
mod_nf_nardaies <- glm(ONSET10 ~ fsca + dah, data = d_nardaies_nf, family = gaussian())
cat("\n===== NON FILTRÉ – Nardaies denses du subalpin =====\n")
print(summary(mod_nf_nardaies))

### 1.2 Queyrellins -----------------------------------------------------------
d_quey_nf <- df_raw %>% filter(habitat_label == "Queyrellins")
mod_nf_quey <- glm(ONSET10 ~ fsca + dah, data = d_quey_nf, family = gaussian())
cat("\n===== NON FILTRÉ – Queyrellins =====\n")
print(summary(mod_nf_quey))

### 1.3 Megaphorbiaies et Aulnaies -------------------------------------------
d_mega_nf <- df_raw %>% filter(habitat_label == "Megaphorbiaies et Aulnaies")
mod_nf_mega <- glm(ONSET10 ~ fsca + dah, data = d_mega_nf, family = gaussian())
cat("\n===== NON FILTRÉ – Megaphorbiaies et Aulnaies =====\n")
print(summary(mod_nf_mega))

### 1.4 P. productives --------------------------------------------------------
d_prod_nf <- df_raw %>% filter(habitat_label == "P. productives")
mod_nf_prod <- glm(ONSET10 ~ fsca + dah, data = d_prod_nf, family = gaussian())
cat("\n===== NON FILTRÉ – P. productives =====\n")
print(summary(mod_nf_prod))

### 1.5 P. nivales ------------------------------------------------------------
d_niv_nf <- df_raw %>% filter(habitat_label == "P. nivales")
mod_nf_niv <- glm(ONSET10 ~ fsca + dah, data = d_niv_nf, family = gaussian())
cat("\n===== NON FILTRÉ – P. nivales =====\n")
print(summary(mod_nf_niv))

### 1.6 P. thermiques écorchées ----------------------------------------------
d_therm_nf <- df_raw %>% filter(habitat_label == "P. thermiques écorchées")
mod_nf_therm <- glm(ONSET10 ~ fsca + dah, data = d_therm_nf, family = gaussian())
cat("\n===== NON FILTRÉ – P. thermiques écorchées =====\n")
print(summary(mod_nf_therm))

# =============================================================================
# 2. Modèles GLM – jeu FILTRÉ
# =============================================================================

### 2.1 Nardaies denses du subalpin ------------------------------------------
d_nardaies_f <- df_filt %>% filter(habitat_label == "Nardaies denses du subalpin")
mod_f_nardaies <- glm(ONSET10 ~ fsca + dah, data = d_nardaies_f, family = gaussian())
cat("\n===== FILTRÉ – Nardaies denses du subalpin =====\n")
print(summary(mod_f_nardaies))

### 2.2 Queyrellins -----------------------------------------------------------
d_quey_f <- df_filt %>% filter(habitat_label == "Queyrellins")
mod_f_quey <- glm(ONSET10 ~ fsca + dah, data = d_quey_f, family = gaussian())
cat("\n===== FILTRÉ – Queyrellins =====\n")
print(summary(mod_f_quey))

### 2.3 Megaphorbiaies et Aulnaies -------------------------------------------
d_mega_f <- df_filt %>% filter(habitat_label == "Megaphorbiaies et Aulnaies")
mod_f_mega <- glm(ONSET10 ~ fsca + dah, data = d_mega_f, family = gaussian())
cat("\n===== FILTRÉ – Megaphorbiaies et Aulnaies =====\n")
print(summary(mod_f_mega))

### 2.4 P. productives --------------------------------------------------------
d_prod_f <- df_filt %>% filter(habitat_label == "P. productives")
mod_f_prod <- glm(ONSET10 ~ fsca + dah, data = d_prod_f, family = gaussian())
cat("\n===== FILTRÉ – P. productives =====\n")
print(summary(mod_f_prod))

### 2.5 P. nivales ------------------------------------------------------------
d_niv_f <- df_filt %>% filter(habitat_label == "P. nivales")
mod_f_niv <- glm(ONSET10 ~ fsca + dah, data = d_niv_f, family = gaussian())
cat("\n===== FILTRÉ – P. nivales =====\n")
print(summary(mod_f_niv))

### 2.6 P. thermiques écorchées ----------------------------------------------
d_therm_f <- df_filt %>% filter(habitat_label == "P. thermiques écorchées")
mod_f_therm <- glm(ONSET10 ~ fsca + dah, data = d_therm_f, family = gaussian())
cat("\n===== FILTRÉ – P. thermiques écorchées =====\n")
print(summary(mod_f_therm))

# -----------------------------------------------------------------------------
# 3. Lecture rapide
#    • Compare visuellement les coefficients fsca & dah entre non filtré / filtré
#      pour voir si l’effet diminue.
# -----------------------------------------------------------------------------





par(mfrow=c(2,2));  plot(mod_nf_nardaies)  # idem pour autres



car::vif(mod_nf_quey)




library(dplyr)

results_iqr <- df_filt %>%                       # jeu filtré
  group_by(habitat_label) %>%
  summarise(
    iqr_fsca = IQR(fsca, na.rm = TRUE),
    iqr_dah  = IQR(dah , na.rm = TRUE),
    coef_fsca = coef(glm(ONSET10 ~ fsca + dah, data = cur_data()))["fsca"],
    coef_dah  = coef(glm(ONSET10 ~ fsca + dah, data = cur_data()))["dah"]
  ) %>%
  mutate(delta_fsca_iqr = coef_fsca * iqr_fsca,
         delta_dah_iqr  = coef_dah  * iqr_dah,
         delta_iqr_max  = pmax(abs(delta_fsca_iqr), abs(delta_dah_iqr))) %>%
  select(habitat_label, starts_with("delta"))

print(results_iqr)























## Controle du ONSET10 avec les paramètres DAH et FSCA 
# PARTIE 1 PLOT 

library(hexbin)
df <- readRDS(dt_alpesud_filtered )   
df2 <- df %>%
  group_by(habitat_label) %>%
  mutate(resid = ONSET10 - median(ONSET10, na.rm = TRUE)) %>%
  ungroup()

p_hex <- ggplot(df2, aes(dah, fsca, z = ONSET10))+
  stat_summary_hex(bins = 30, fun = median)+
  scale_fill_distiller(palette = "Spectral", name = "ONSET10")+
  facet_wrap(~ habitat_label, scales = "free")+
  labs(title = "Médiane ONSET10 – Hexbin 2D",
       x = "dah", y = "fSCA")+
  theme_bw(9)

p_violin <- ggplot(df2, aes(habitat_label, resid, fill = habitat_label))+
  geom_violin(trim = FALSE, colour = NA, alpha = .8)+
  geom_hline(yintercept = 0, colour = "grey30")+
  coord_flip()+
  scale_fill_brewer(palette = "Set3", guide = "none")+
  labs(title = "Résidus ONSET10\n(médiane retirée)",
       x = NULL, y = "ONSET10 – médiane")+
  theme_bw(9)

library(patchwork)
p_hex + p_violin + plot_layout(widths = c(3,1))



# MODEL
# Dsitribution du ONSET 10


library(ggplot2)
library(dplyr)

# df = votre jeu "fenêtre contrôlée" (une ligne = 1 pixel unique)
# ---------------------------------------------------------------

## Histogramme global + densité
ggplot(df, aes(ONSET10)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "steelblue", colour = "white") +
  geom_density(colour = "red", linewidth = .8, adjust = 1.2) +
  labs(title = "Distribution globale de ONSET10",
       x = "Jour de l'année", y = "Densité") +
  theme_bw()

## Histogrammes séparés par habitat
ggplot(df, aes(ONSET10)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "white") +
  facet_wrap(~ habitat_label, scales = "free_y") +
  labs(title = "ONSET10 — distribution par habitat",
       x = "Jour de l'année", y = "Pixels") +
  theme_bw(base_size = 9)




# Q–Q plot global
qqnorm(df$ONSET10); qqline(df$ONSET10, col = "red")

# Test de Shapiro sur un échantillon aléatoire (10 000 px max)
set.seed(1)
shapiro.test(sample(df$ONSET10, size = min(10000, nrow(df))))


library(dplyr)
library(broom)
library(purrr)

## df = pixels retenus (1 ligne = 1 pixel, 1 habitat)
## ............................................................................

analyse_hab <- function(dat) {
  # Standardisation (une unité = 1 écart‑type dans l'habitat)
  dat <- mutate(dat,
                dah_z  = scale(dah )[ ,1],
                fsca_z = scale(fsca)[ ,1])
  
  # modèles imbriqués
  m0 <- lm(ONSET10 ~ 1,                 data = dat)          # nul
  m1 <- lm(ONSET10 ~ dah_z + fsca_z,   data = dat)          # complet
  
  # comparaison F (équiv. anova type I ici)
  comp <- anova(m0, m1)
  p_F  <- comp$`Pr(>F)`[2]
  dAIC <- AIC(m1) - AIC(m0)
  R2   <- summary(m1)$r.squared
  
  # pentes + IC95
  coef_tab <- tidy(m1, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, conf.low, conf.high)
  
  tibble(
    n_pix      = nrow(dat),
    p_F        = p_F,
    delta_AIC  = dAIC,
    R2         = R2,
    dah_eff    = coef_tab$estimate[coef_tab$term=="dah_z"],
    dah_low    = coef_tab$conf.low[coef_tab$term=="dah_z"],
    dah_high   = coef_tab$conf.high[coef_tab$term=="dah_z"],
    fsca_eff   = coef_tab$estimate[coef_tab$term=="fsca_z"],
    fsca_low   = coef_tab$conf.low[coef_tab$term=="fsca_z"],
    fsca_high  = coef_tab$conf.high[coef_tab$term=="fsca_z"]
  )
}

results <- df %>%
  split(.$habitat_label) %>%
  map_df(analyse_hab, .id = "habitat")

print(results, digits = 3)








sd_dah  <- sd(df$dah , na.rm = TRUE)
sd_fsca <- sd(df$fsca, na.rm = TRUE)

eff_dah_real  <- 0.81  * sd_dah     # jours / unité dah
eff_fsca_real <- 1.24 * sd_fsca     # jours / % neige



library(performance)
mods  <- lapply(split(df, df$habitat_label), \(d) 
                lmer(ONSET10 ~ dah_z + fsca_z + (1|cell), data = d))

tbl <- map_dfr(mods, \(m) {
  tibble(
    habitat = sub("\\.\\d+$","",names(m)),
    r2      = r2(m)$R2_marginal,
    eff_dah = fixef(m)["dah_z"],
    eff_fsc = fixef(m)["fsca_z"],
    se_dah  = sqrt(diag(vcov(m)))["dah_z"],
    se_fsc  = sqrt(diag(vcov(m)))["fsca_z"]
  )
) %>% mutate(across(starts_with("eff"),
                    ~sprintf("%.2f ± %.2f", .x, get(paste0("se_", sub("eff_","",cur_column()))))))

print(tbl)









df_filtered <- readRDS(dt_alpesud_filtered ) 
df_nofiltered <- readRDS(dt_alpesud_nofiltered ) 







library(dplyr)
library(broom)

analyse_ds <- function(df) {
  df %>%
    group_by(habitat_label) %>%
    mutate(   # centrage‑réduction par habitat
      dah_z  = scale(dah )[ ,1],
      fsca_z = scale(fsca)[ ,1]
    ) %>%
    group_modify(~{
      m0 <- lm(ONSET10 ~ 1,              data = .x)
      m1 <- lm(ONSET10 ~ dah_z + fsca_z, data = .x)
      comp <- anova(m0, m1)
      tibble(
        n_pix     = nrow(.x),
        p_F       = comp$`Pr(>F)`[2],
        delta_AIC = AIC(m1) - AIC(m0),
        R2        = summary(m1)$r.squared,
        dah_est   = coef(m1)["dah_z"],
        fsca_est  = coef(m1)["fsca_z"]
      )
    }) %>% ungroup()
}






res_nofilt <- analyse_ds(df_nofiltered) %>% mutate(dataset = "nofilter")
res_filt   <- analyse_ds(df_filtered   ) %>% mutate(dataset = "filtered")

res_all <- bind_rows(res_nofilt, res_filt)







library(tidyr)

comp <- res_all %>%
  select(habitat_label, dataset, R2, dah_est, fsca_est) %>%
  pivot_wider(names_from = dataset,
              values_from = c(R2, dah_est, fsca_est),
              names_glue = "{.value}_{dataset}") %>%
  mutate(
    dR2   = R2_filtered   - R2_nofilter,
    dDah  = abs(dah_est_filtered )  - abs(dah_est_nofilter),
    dFsca = abs(fsca_est_filtered)  - abs(fsca_est_nofilter)
  )
print(comp, digits = 3)



wilcox.test(comp$dR2,   alternative = "less")   # R² devrait baisser
wilcox.test(comp$dDah,  alternative = "less")   # pente dah plus petite ?
wilcox.test(comp$dFsca, alternative = "less")   # pente fSCA plus petite ?





library(ggplot2)

# a) barplot des R² avant / après
ggplot(res_all, aes(x = habitat_label, y = R2,
                    fill = dataset)) +
  geom_col(position = "dodge") +
  labs(title = "R² du modèle ONSET10 ~ dah + fSCA",
       subtitle = "Sans fenêtre vs. fenêtre ±10 %",
       y = "R² marginal", x = NULL) +
  theme_bw(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# b) forest‑plot des pentes (dah et fSCA)
res_long <- res_all %>%
  pivot_longer(c(dah_est, fsca_est),
               names_to = "var", values_to = "slope")

ggplot(res_long,
       aes(y = habitat_label, x = slope,
           colour = dataset, shape = dataset)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(size = 2, position = position_dodge(width = .5)) +
  facet_wrap(~ var, scales = "free_x",
             labeller = labeller(var = c(dah_est="dah_z (jours/σ)",
                                         fsca_est="fSCA_z (jours/σ)"))) +
  labs(title = "Pentes avant / après fenêtre",
       x = "Pente (jours par écart‑type)", y = NULL) +
  theme_bw(9)





































readRDS(out_rds)



###############################################################################
#  MAXV ~ Charge (pixels retenus ±10 % dah & fsca)  – par habitat
###############################################################################
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(scales) })

## ---- 1. lecture ------------------------------------------------------------
df <- readRDS(out_rds)                                # chemin _WIN10pct.rds
df <- df %>% filter(!is.na(MAXV), !is.na(Charge), Charge <= 1000)

## ---- 2. graphique ----------------------------------------------------------
p <- ggplot(df, aes(x = Charge, y = MAXV)) +
  geom_point(alpha = 0.12, size = 0.7, colour = "steelblue") +
  geom_smooth(method = "loess", span = 1.1, se = FALSE,
              colour = "red", linewidth = 0.9) +
  scale_x_continuous(limits = c(0, 1000),
                     breaks = seq(0, 1000, 200),
                     labels = comma) +
  facet_wrap(~ habitat_label, scales = "free_y") +
  labs(title    = "Influence du chargement (0–1000 UA) sur MAXV",
       subtitle = sprintf("%s – %d  |  fenêtre ±10 %% dah & fsca", alpage, YEAR),
       x        = "Charge (UA)",
       y        = "MAXV (valeur au pic de saison)",
       caption  = paste("Pixels retenus :", format(nrow(df), big.mark = " "))) +
  theme_bw(base_size = 10)

print(p)


































## ────────────────────────────────────────────────────────────────
##   MAXV ~ Charge  (pixels pâturés 0‑1000 UA)  –  par habitat
## ────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

## 1. paramètres d’E/S --------------------------------------------------------
YEAR    <- 2022
alpage  <- "Cayolle"
out_rds <- file.path("outputs", "10. NDVI Effect",
                     paste0("dataset_", alpage),
                     sprintf("dataset_pheno_LONG_%d_%s_WIN10pct.rds",
                             YEAR, alpage))

## habitats à garder
hab_keep <- c("Nardaies denses du subalpin", "Queyrellins",
              "Megaphorbiaies et Aulnaies", "P. productives",
              "P. nivales", "P. thermiques écorchées")

## 2. lecture + filtre  -------------------------------------------------------
df <- readRDS(out_rds) %>%
  filter(habitat_label %in% hab_keep,
         !is.na(MAXV), !is.na(Charge),
         Charge > 0, Charge <= 1000)          # limite 0–1000 UA

## 3. plot -------------------------------------------------------------------
p <- ggplot(df, aes(Charge, MAXV)) +
  geom_point(alpha = 0.12, size = 0.7, colour = "steelblue") +
  geom_smooth(method = "loess", span = 1.1, se = FALSE,
              colour = "red", linewidth = 0.9, family = "gaussian") +
  scale_x_continuous(limits = c(0, 1000),
                     breaks = seq(0, 1000, 200),
                     labels = comma) +
  facet_wrap(~ habitat_label, scales = "free_y") +
  labs(title    = "Influence du gradient de chargement (0–1000 UA) sur MAXV",
       subtitle = sprintf("%s – %d", alpage, YEAR),
       x        = "Charge (Unité Animale)",
       y        = "MAXV (valeur au pic de saison)",
       caption  = paste("Pixels pâturés avec Charge 0–1000 UA  •  n =",
                        format(nrow(df), big.mark = " "))) +
  theme_bw(base_size = 10)

print(p)  # affiche dans R / RStudio

## 4. sauvegarde optional -----------------------------------------------------
out_png <- file.path("outputs", "10. NDVI Effect",
                     paste0("dataset_", alpage),
                     sprintf("Charge_MAXV_%d_%s_0-1000.png", YEAR, alpage))
ggsave(out_png, plot = p, width = 9, height = 6, dpi = 300)
message("✓ Figure enregistrée : ", out_png)








































































# ──────────────────────────────────────────────────────────────────────────────
# 0. Packages -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# ──────────────────────────────────────────────────────────────────────────────
# 1. Lecture -------------------------------------------------------------------
YEAR   <- 2022
alpage <- "Cayolle"

root_dir   <- getwd()
rds_file   <- file.path(
  root_dir, "outputs/10. NDVI Effect",
  paste0("dataset_", alpage),
  sprintf("dataset_pheno_LONG_%d_%s_TERCILES.rds", YEAR, alpage)
)
df <- readRDS(rds_file)

# ── 1.1 Harmonisation du nom de colonne ---------------------------------------
# 'load_class' → 'traitement' (Faible / Moyen / Fort)
df <- df %>%
  mutate(traitement = factor(load_class, levels = c("Faible", "Moyen", "Fort")))

# ──────────────────────────────────────────────────────────────────────────────
# 2. Palette & libellés --------------------------------------------------------
pal_trt <- c(
  Faible = "khaki3",
  Moyen  = "seagreen3",
  Fort   = "steelblue"
)
lab_trt <- c(
  Faible = "low grazing intensity",
  Moyen  = "medium grazing intensity",
  Fort   = "high grazing intensity"
)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Boucle par habitat --------------------------------------------------------
habitats <- sort(unique(df$habitat_label))
plots <- vector("list", length(habitats)); names(plots) <- habitats

for (hab in habitats) {
  
  df_h <- filter(df, habitat_label == hab)
  
  ## 3.1 NDVI : médiane + bande 20–80 %
  ndvi_stats <- df_h %>%
    group_by(DOY, traitement) %>%
    summarise(
      med_ndvi   = median(NDVI, na.rm = TRUE),
      low80_ndvi = quantile(NDVI, .20, na.rm = TRUE),
      hi80_ndvi  = quantile(NDVI, .80, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    group_by(traitement) %>%                # retire classe vide (NDVI = 0 partout)
    filter(any(med_ndvi > 0)) %>%
    ungroup()
  
  if (nrow(ndvi_stats) == 0) next           # rien à tracer
  
  ## 3.2 Pic IRG par traitement
  irg_peaks <- df_h %>%
    group_by(DOY, traitement) %>%
    summarise(med_irg = median(IRG, na.rm = TRUE), .groups = "drop") %>%
    filter(traitement %in% ndvi_stats$traitement) %>%
    group_by(traitement) %>%
    slice_max(med_irg, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(DOY_peak = DOY, irg_max = med_irg)
  
  ## 3.3 Axes & padding
  y_min <- min(ndvi_stats$low80_ndvi, na.rm = TRUE)
  y_max <- max(ndvi_stats$hi80_ndvi , na.rm = TRUE)
  pad   <- max(1e-6, .10 * (y_max - y_min))
  
  baseline_y <- y_min
  peak_y     <- baseline_y + .90 * (y_max - y_min)
  sf         <- (peak_y - baseline_y) / max(irg_peaks$irg_max)
  y_max_all  <- max(y_max, peak_y) + .30 * pad
  
  ## 3.4 Palette réduite aux classes présentes
  trt_present <- unique(ndvi_stats$traitement)
  cols_here   <- pal_trt[trt_present]
  labs_here   <- lab_trt[trt_present]
  
  ## 3.5 Graphique
  p <- ggplot() +
    geom_ribbon(
      data = ndvi_stats,
      aes(DOY, ymin = low80_ndvi, ymax = hi80_ndvi, fill = traitement),
      alpha = .25, colour = NA
    ) +
    geom_line(
      data = ndvi_stats,
      aes(DOY, med_ndvi, colour = traitement),
      size = 1
    ) +
    geom_segment(
      data = irg_peaks,
      aes(x = DOY_peak, xend = DOY_peak,
          y = baseline_y, yend = irg_max * sf,
          colour = traitement),
      size = 2
    ) +
    geom_hline(
      data = irg_peaks,
      aes(yintercept = irg_max * sf, colour = traitement),
      linetype = "dashed", size = .8
    ) +
    scale_x_continuous(
      breaks = seq(60, 330, 30), minor_breaks = seq(60, 330, 10)
    ) +
    scale_y_continuous(
      name     = "Plant Phenology Index (PPI)",
      limits   = c(baseline_y, y_max_all),
      expand   = c(0, 0),
      sec.axis = sec_axis(~ . / sf, name = "Green‑up rate (IRG)")
    ) +
    scale_colour_manual(values = cols_here, breaks = trt_present,
                        labels = labs_here, name = NULL) +
    scale_fill_manual(values = cols_here, breaks = trt_present,
                      labels = labs_here, name = NULL) +
    labs(title = hab, x = "Julian date") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_line(size = .2, linetype = "dotted"),
      panel.grid.major = element_line(size = .3),
      legend.position  = c(.78, .92),
      legend.background = element_blank()
    )
  
  plots[[hab]] <- p
}

# ──────────────────────────────────────────────────────────────────────────────
# 4. Affichage ----------------------------------------------------------------
for (p in plots) print(p)




























## BLABLA























































































































































































## APPROCHE AVEC PLUSIEUR ALPAGE
library(dplyr)
library(purrr)   # pour map_dfr()

build_ndvi_dataset <- function(year,
                               alpages,
                               root_dir   = getwd(),
                               output_dir = file.path(root_dir, "outputs"),
                               write_fst  = FALSE) {
  ## ————————————————————————————————————————————————————————
  ## 1. Normalisation des alpages
  alpages <- sort(unique(alpages))
  combo   <- paste(alpages, collapse = "_")
  
  ## 2. Fonction de lecture + filtrage
  load_one <- function(ap) {
    in_file <- file.path(output_dir, "10. NDVI Effect",
                         paste0("dataset_", ap),
                         sprintf("dataset_pheno_LONG_%d_%s.rds", year, ap))
    if (!file.exists(in_file))
      stop("Fichier manquant : ", in_file, call. = FALSE)
    
    readRDS(in_file) |>
      filter(!is.na(IRG) & !is.na(NDVI)) |>      # ⬅️ filtrage mis à jour
      mutate(alpage = ap)
  }
  
  ## 3. Concaténation
  dataset_all <- map_dfr(alpages, load_one)
  
  ## 4. Sauvegarde
  out_dir <- file.path(output_dir, "10. NDVI Effect",
                       paste0("dataset_", combo))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_rds <- file.path(out_dir,
                       sprintf("dataset_pheno_LONG_%d_%s.rds", year, combo))
  saveRDS(dataset_all, out_rds)
  message("Sauvegardé : ", out_rds)
  
  if (write_fst) {
    if (!requireNamespace("fst", quietly = TRUE))
      stop("Installez le package {fst} ou mettez write_fst = FALSE.")
    fst::write_fst(dataset_all, sub("\\.rds$", ".fst", out_rds))
  }
  
  invisible(dataset_all)
}

## Exemple d’appel
alpages <- c("Viso", "Cayolle", "Sanguiniere")
build_ndvi_dataset(year = 2022, alpages = alpages)






# # Option .fst si besoin
# write_fst(dataset_all, sub("\\.rds$", ".fst", out_rds))
 readRDS(out_rds)
 
 
 
 
 
 
 

 
 
 # 1) Tu charges / construis ton data.frame concaténé comme avant
 df_all <- readRDS("C:/Users/massocam/Documents/STAGE_M2_PERSEE/R_studio/PERSEE_Traitement_Catlog/outputs/10. NDVI Effect/dataset_Cayolle_Sanguiniere_Viso/dataset_pheno_LONG_2022_Cayolle_Sanguiniere_Viso.rds")
 
 # LIBRARY & FUNCTION
 source(file.path(functions_dir, "Functions_phenology.R"))
 
 ## Exige au moins 100 pixels dans chaque classe
 df_filtre <- optimize_loading_windows2(
   df_all,
   min_pixels = c(faible = 10L, moyen = 10L, fort = 1L)
 )
 
 # 3) Si tu veux ensuite sauvegarder ce sous‑ensemble :
 saveRDS(df_filtre,
         "outputs/10. NDVI Effect/dataset_Cayolle_Sanguiniere_Viso/selection_optimisee.rds")
 
 
 
 
 
 
 
 
 
 # ─────────────────────────────────────────────────────────────────
 # 0. Packages -----------------------------------------------------
 suppressPackageStartupMessages({
   library(dplyr)
   library(ggplot2)
   library(stringr)   # pour str_to_title()
 })
 
 # ─────────────────────────────────────────────────────────────────
 # 1. Lecture ------------------------------------------------------
 rds_file <- file.path(
   "C:/Users/massocam/Documents/STAGE_M2_PERSEE/R_studio",
   "PERSEE_Traitement_Catlog/outputs/10. NDVI Effect",
   "dataset_Cayolle_Sanguiniere_Viso/selection_optimisee.rds"
 )
 df_zones <- readRDS(rds_file)
 
 # ── 1.1 Colonne 'traitement' normalisée --------------------------
 # load_class = « faible » / « moyen » / « fort »  →  Faible / Moyen / Fort
 df_zones <- df_zones %>%
   mutate(
     traitement = str_to_title(load_class),                 # f → F
     traitement = factor(traitement, levels = c("Faible", "Moyen", "Fort"))
   )
 
 # ─────────────────────────────────────────────────────────────────
 # 2. Palette & libellés ------------------------------------------
 pal_trt <- c(
   Faible = "khaki3",
   Moyen  = "seagreen3",
   Fort   = "steelblue"
 )
 lab_trt <- c(
   Faible = "low grazing intensity",
   Moyen  = "medium grazing intensity",
   Fort   = "high grazing intensity"
 )
 
 # ─────────────────────────────────────────────────────────────────
 # 3. Boucle par habitat ------------------------------------------
 habitats <- sort(unique(df_zones$habitat_label))
 plots <- vector("list", length(habitats)); names(plots) <- habitats
 
 for (hab in habitats) {
   
   df_h <- filter(df_zones, habitat_label == hab)
   
   ## 3.1 Statistiques NDVI (médiane + bande 20–80 %)
   ndvi_stats <- df_h %>%
     group_by(DOY, traitement) %>%
     summarise(
       med_ndvi   = median(NDVI, na.rm = TRUE),
       low80_ndvi = quantile(NDVI, .20, na.rm = TRUE),
       hi80_ndvi  = quantile(NDVI, .80, na.rm = TRUE),
       .groups    = "drop"
     )
   
   ## 3.2 Pic IRG par traitement
   irg_peaks <- df_h %>%
     group_by(DOY, traitement) %>%
     summarise(med_irg = median(IRG, na.rm = TRUE), .groups = "drop") %>%
     group_by(traitement) %>%
     slice_max(med_irg, n = 1, with_ties = FALSE) %>%
     ungroup() %>%
     rename(DOY_peak = DOY, irg_max = med_irg)
   
   ## 3.3 Axes & padding
   y_min <- min(ndvi_stats$low80_ndvi, na.rm = TRUE)
   y_max <- max(ndvi_stats$hi80_ndvi , na.rm = TRUE)
   pad   <- .10 * (y_max - y_min)
   
   baseline_y <- y_min
   peak_y     <- baseline_y + .90 * (y_max - y_min)
   sf         <- (peak_y - baseline_y) / max(irg_peaks$irg_max)
   y_max_all  <- max(y_max, peak_y) + .30 * pad
   
   ## 3.4 Palette réduite aux classes présentes
   trt_present <- levels(df_h$traitement)[levels(df_h$traitement) %in% df_h$traitement]
   cols_here   <- pal_trt[trt_present]
   labs_here   <- lab_trt[trt_present]
   
   ## 3.5 Construction du graphique
   p <- ggplot() +
     geom_ribbon(
       data = ndvi_stats,
       aes(DOY, ymin = low80_ndvi, ymax = hi80_ndvi, fill = traitement),
       alpha = .25, colour = NA
     ) +
     geom_line(
       data = ndvi_stats,
       aes(DOY, med_ndvi, colour = traitement),
       size = 1
     ) +
     geom_segment(
       data = irg_peaks,
       aes(x = DOY_peak, xend = DOY_peak,
           y = baseline_y, yend = irg_max * sf,
           colour = traitement),
       size = 2
     ) +
     geom_hline(
       data = irg_peaks,
       aes(yintercept = irg_max * sf, colour = traitement),
       linetype = "dashed", size = .8
     ) +
     scale_x_continuous(
       breaks = seq(60, 330, 30), minor_breaks = seq(60, 330, 10)
     ) +
     scale_y_continuous(
       name     = "Plant Phenology Index (PPI)",
       limits   = c(baseline_y, y_max_all),
       expand   = c(0, 0),
       sec.axis = sec_axis(~ . / sf, name = "Green‑up rate (IRG)")
     ) +
     scale_colour_manual(values = cols_here, breaks = trt_present,
                         labels = labs_here, name = NULL) +
     scale_fill_manual(values = cols_here, breaks = trt_present,
                       labels = labs_here, name = NULL) +
     labs(title = hab, x = "Julian date") +
     theme_bw(base_size = 10) +
     theme(
       panel.grid.minor = element_line(size = .2, linetype = "dotted"),
       panel.grid.major = element_line(size = .3),
       legend.position  = c(.78, .92),
       legend.background = element_blank()
     )
   
   plots[[hab]] <- p
 }
 
 # ─────────────────────────────────────────────────────────────────
 # 4. Affichage ----------------------------------------------------
 for (p in plots) print(p)
 
 
 
 
 
 
 
 
 
 
 
 df_zones %>%                                     # ➊ lignes par classe/habitat
   count(habitat_label, traitement) %>%
   tidyr::pivot_wider(names_from = traitement, values_from = n, values_fill = 0)
 
 df_zones %>%                                     # ➋ NDVI médian par classe
   group_by(habitat_label, traitement) %>%
   summarise(med_NDVI = median(NDVI), .groups = "drop")
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ## NOUVELLE APPROCHE AVEC UNE ACP SUR TOUTE LES VARIABLES : 
 
 
 # APPROCHE MONO-ALPAGE
 
 ## SORTIE
 # Dossier de sortie
 out_dataset       <- file.path(output_dir, "10. NDVI Effect", paste0("dataset_",alpage))
 dir.create(out_dataset, recursive=TRUE, showWarnings=FALSE)
 # Un .fst contetant l'ensemble du dataset 
 out_rds <- file.path(out_dataset, paste0("dataset_pheno_LONG_", YEAR, "_", alpage, ".rds"))
 
 
 
 
 
 
 
 
 
 
 
 
 suppressPackageStartupMessages({
   library(dplyr);  library(readr);  library(ggplot2)
   library(factoextra);  library(caret);  library(scales)
   library(conflicted)
 })
 
 ## préférences anti‑conflits
 conflict_prefer("intersect", "base",  quiet = TRUE)
 conflict_prefer("first",     "dplyr", quiet = TRUE)
 conflict_prefer("filter",    "dplyr", quiet = TRUE)   # ← nouveau
 conflict_prefer("select",    "dplyr", quiet = TRUE)   # (optionnel mais utile)
 conflict_prefer("summarise", "dplyr", quiet = TRUE)   # idem
 
 ## paramètres ---------------------------------------------------------------
 set.seed(42)
 sample_max <- 2e5
 vars_pheno <- c("MINV","MAXV","AMPL","LENGTH","ONSET10","MaxSlope",
                 "GreenUpDur","GreenDownDur","AsymSlope",
                 "LSLOPE","RSLOPE",
                 "SOSD","EOSD","MAXD","OFFSET","GROWTH","SENESC")
 
 ## 1. lecture & dé‑doublonnage par pixel ------------------------------------
 df <- read_rds(out_rds) %>%
   filter(!is.na(Charge) & Charge >= 1) %>%          # pixels pâturés
   select(cell, Charge, all_of(vars_pheno)) %>%
   group_by(cell) %>%                                # 1 ligne / pixel
   summarise(across(everything(),
                    ~ dplyr::first(na.omit(.x))),     # ← explicite dplyr::first
             .groups = "drop")
 
 ## 2. nettoyage Inf / NA -----------------------------------------------------
 df <- df %>%
   mutate(across(all_of(vars_pheno),
                 ~ ifelse(is.infinite(.x), NA, .x))) %>%  # Inf → NA
   drop_na()                                            # retire lignes NA
 
 ## sous‑échantillonnage optionnel
 if (nrow(df) > sample_max) df <- slice_sample(df, n = sample_max)
 
 ## 3. standardisation
 X <- scale(df %>% select(all_of(vars_pheno)))
 y <- df$Charge
 
 ## 4. ACP --------------------------------------------------------------------
 pca <- prcomp(X, center = FALSE, scale. = FALSE)
 fviz_eig(pca, addlabels = TRUE)            # scree‑plot
 fviz_pca_var(pca, col.var = "contrib",
              gradient.cols = c("#00AFBB","#E7B800","#FC4E07"),
              repel = TRUE)
 
 ## 5. corrélation Charge × composantes
 cors <- cor(pca$x, y)
 print(round(cors,3))
 
 ## 6. importance via random‑forest ------------------------------------------
 ctrl <- trainControl(method = "cv", number = 5)
 rf   <- train(x = X, y = y,
               method = "rf", trControl = ctrl, importance = TRUE)
 imp  <- varImp(rf)$importance %>%
   tibble::rownames_to_column("var") %>%
   arrange(desc(Overall))
 
 ggplot(imp, aes(reorder(var, Overall), Overall)) +
   geom_col(fill = "steelblue") +
   coord_flip() +
   labs(x=NULL, y="Importance (RF)",
        title="Importance des variables sur la Charge") +
   theme_minimal()
 
 ## 7. résumé console ---------------------------------------------------------
 cat("\nCorrélations Charge‑PC :\n"); print(round(cors,3))
 cat("\nTop 10 variables importantes (RF) :\n"); print(head(imp,10))
 