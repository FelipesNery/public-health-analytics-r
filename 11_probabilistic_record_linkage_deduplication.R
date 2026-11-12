# =============================================================================
# PROJECT 11 - Probabilistic Record Linkage and Deduplication of Health Data
# =============================================================================
#
# DESCRIPTION
# Identifying duplicate records of the same individual within a notifiable
# disease database (deduplication) and linking two different databases that
# lack a unique identifier (record linkage), using name, mother's name, sex
# and date of birth. An essential step to build cohorts from health
# information systems (master's thesis project).
#
# WHAT WAS DONE
#  - Name standardization: upper case, removal of accents, punctuation,
#    Portuguese name particles (DE, DA, DOS...) and invalid values
#    (e.g., "SEM INFORMACAO" = "no information")
#  - Splitting first name, middle name and last surname
#  - Splitting date of birth into day, month and year
#  - Blocking on multiple keys and phonetic encoding (Soundex)
#  - String comparison (Jaro-Winkler / Levenshtein) and EpiLink weights
#  - Classification of true pairs using a score cut-off
#  - Creation of a unique person identifier (person ID)
#  - Linkage between two different databases (RLBigDataLinkage)
#
# TECHNIQUES / PACKAGES: RecordLinkage, stringr, tidyverse
# NOTE: FICTITIOUS randomly generated names and data; variables anonymized.
# =============================================================================

library(tidyverse)
library(stringr)

set.seed(1123)

# -----------------------------------------------------------------------------
# 1. Fictitious database with duplicates and typos
# -----------------------------------------------------------------------------
pool_n1 <- c("ANA", "MARIA", "JOSE", "JOAO", "PAULO", "LUCAS", "JULIA", "PEDRO",
             "CARLA", "RAFAEL", "BRUNA", "MARCOS", "LUIZA", "TIAGO", "SOFIA")
pool_n2 <- c("CLARA", "LUIZ", "EDUARDO", "BEATRIZ", "HENRIQUE", "VITORIA", "")
pool_sn <- c("SILVA", "SOUZA", "OLIVEIRA", "SANTOS", "LIMA", "PEREIRA", "COSTA",
             "RODRIGUES", "ALMEIDA", "NASCIMENTO", "FERREIRA", "CARVALHO")

make_name <- function(k) {
  str_squish(paste(sample(pool_n1, k, TRUE), sample(pool_n2, k, TRUE),
                   sample(c("DE", "DA", "DOS", ""), k, TRUE),
                   sample(pool_sn, k, TRUE), sample(pool_sn, k, TRUE)))
}

n_ind <- 600
people_lk <- tibble(
  nm_pt  = make_name(n_ind),
  nm_mo  = make_name(n_ind),
  dob_pt = format(as.Date("1950-01-01") + sample(0:25000, n_ind, TRUE), "%d/%m/%Y"),
  sex_pt = sample(c("M", "F"), n_ind, TRUE)
)

# Introduce typos (replace one letter)
add_typo <- function(x) {
  pos <- sample(2:(nchar(x) - 1), 1)
  substr(x, pos, pos) <- sample(LETTERS, 1)
  x
}

dup_lk <- people_lk |>
  slice_sample(n = 120) |>
  mutate(nm_pt = map_chr(nm_pt, add_typo),
         nm_mo = if_else(runif(n()) < 0.3, "SEM INFORMACAO", nm_mo))

notif_lk <- bind_rows(people_lk, dup_lk) |>
  slice_sample(prop = 1) |>
  mutate(dt_nt  = as.Date("2015-01-01") + sample(0:1800, n(), TRUE),
         cod_nt = sprintf("N%06d", row_number()))

# -----------------------------------------------------------------------------
# 2. Name standardization
# -----------------------------------------------------------------------------
# Common placeholders found in Brazilian records meaning "unknown"
invalid_vals <- c("", "NA", "SEM INFO", "SEMINFO", "SEM INFORMACAO", "SEM REGISTRO",
                  "NAO INFORMADO", "NAO DECLARADA", "NAO TEM", "IGNORADO")

std_name <- function(x) {
  x |>
    toupper() |>
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") |>
    str_remove_all("[^[:alpha:][:space:]]") |>
    str_replace_all("\\b(DE|DA|DO|DAS|DOS|E)\\b", " ") |>   # name particles
    str_squish() |>
    (\(v) if_else(v %in% invalid_vals, NA_character_, v))()
}

notif_lk <- notif_lk |>
  mutate(
    nm_pt_s = std_name(nm_pt),
    nm_mo_s = std_name(nm_mo),
    # patient name parts
    p_first  = word(nm_pt_s, 1),
    p_middle = if_else(str_count(nm_pt_s, "\\S+") > 2, word(nm_pt_s, 2), NA_character_),
    p_last   = if_else(str_count(nm_pt_s, "\\S+") > 1, word(nm_pt_s, -1), NA_character_),
    # mother's name parts
    m_first  = word(nm_mo_s, 1),
    m_middle = if_else(str_count(nm_mo_s, "\\S+") > 2, word(nm_mo_s, 2), NA_character_),
    m_last   = if_else(str_count(nm_mo_s, "\\S+") > 1, word(nm_mo_s, -1), NA_character_),
    # date of birth split
    dob_d = as.integer(substr(dob_pt, 1, 2)),
    dob_m = as.integer(substr(dob_pt, 4, 5)),
    dob_y = as.integer(substr(dob_pt, 7, 10))
  )

notif_lk |> select(nm_pt, nm_pt_s, p_first, p_middle, p_last) |> head()

# -----------------------------------------------------------------------------
# 3. Matching variables (the order defines the column indices)
# -----------------------------------------------------------------------------
lk_vars <- c("nm_pt_s",   # 1
             "nm_mo_s",   # 2
             "p_first",   # 3
             "p_middle",  # 4
             "p_last",    # 5
             "m_first",   # 6
             "m_middle",  # 7
             "m_last",    # 8
             "dob_d",     # 9
             "dob_m",     # 10
             "dob_y",     # 11
             "sex_pt",    # 12
             "cod_nt")    # 13 (identifier only - excluded from comparison)

df_lk <- as.data.frame(notif_lk[, lk_vars])

# -----------------------------------------------------------------------------
# 4. Probabilistic deduplication
# -----------------------------------------------------------------------------
if (requireNamespace("RecordLinkage", quietly = TRUE)) {
  library(RecordLinkage)

  rl_dd <- RLBigDataDedup(
    df_lk,
    blockfld  = list(c(3, 5), c(9, 10, 11), c(6, 8), c(1, 2, 12)),  # blocks
    exclude   = 13,
    strcmp    = c(1, 2, 4, 7),
    strcmpfun = "jarowinkler",
    phonetic  = c(3, 5, 6, 8)                                        # Soundex
  )

  rl_dd <- epiWeights(rl_dd)

  # Inspect the weight distribution to choose the cut-off
  all_pairs <- getPairs(rl_dd, single.rows = TRUE)
  hist(all_pairs$Weight, breaks = 50, main = "Distribution of EpiLink weights",
       xlab = "Weight")

  cutoff <- 0.62          # chosen after manual review of the grey zone

  res_dd   <- epiClassify(rl_dd, threshold.upper = cutoff)
  true_prs <- getPairs(res_dd, filter.link = "link", single.rows = TRUE)

  true_prs |>
    select(id.1, id.2, nm_pt_s.1, nm_pt_s.2, nm_mo_s.1, nm_mo_s.2, Weight) |>
    arrange(Weight) |>
    head(20)

  # ---------------------------------------------------------------------------
  # 5. Unique person identifier
  # ---------------------------------------------------------------------------
  id_map <- true_prs |>
    select(id.1, id.2) |>
    arrange(id.1) |>
    filter(!duplicated(id.2))

  notif_lk <- notif_lk |>
    mutate(row_id = row_number()) |>
    left_join(id_map, by = c("row_id" = "id.2")) |>
    mutate(person_id = if_else(is.na(id.1), row_id, id.1)) |>
    select(-id.1) |>
    arrange(person_id, dt_nt)

  # how many individuals have more than one notification
  notif_lk |> count(person_id) |> filter(n > 1) |> nrow()

  # keep only the first notification of each person
  unique_db <- notif_lk |> group_by(person_id) |> slice_min(dt_nt, n = 1) |> ungroup()

  # ---------------------------------------------------------------------------
  # 6. Linkage between two different databases
  # ---------------------------------------------------------------------------
  db_a <- df_lk[1:400, ]
  db_b <- df_lk[sample(nrow(df_lk), 300), ]

  rl_lk <- RLBigDataLinkage(
    db_a, db_b,
    blockfld  = list(c(3, 5), c(9, 10, 11), c(6, 8)),
    exclude   = 13,
    strcmp    = c(1, 2),
    strcmpfun = "levenshtein",
    phonetic  = c(3, 6)
  )
  rl_lk <- epiWeights(rl_lk)
  pairs_ab <- getPairs(rl_lk, min.weight = cutoff, single.rows = TRUE)
  pairs_ab |> select(nm_pt_s.1, nm_pt_s.2, nm_mo_s.1, nm_mo_s.2, Weight) |> head()
}
