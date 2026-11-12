# =============================================================================
# PROJECT 14 - ETL Pipeline for Large Health Databases: DATASUS + DuckDB +
#              Parquet
# =============================================================================
#
# DESCRIPTION
# Pipeline to extract, standardize and classify large volumes of health data
# (millions of records) that do not fit comfortably in memory:
# (1) automated download of microdata from DATASUS (Brazilian Ministry of
# Health) via FTP, (2) reading .dbc files and converting them to Parquet and
# (3) standardization and classification with SQL in DuckDB (e.g., lab test
# results from different providers grouped by pathogen and subtype).
#
# WHAT WAS DONE
#  - Listing files on the public FTP server and batch download with tryCatch
#  - Reading .dbc files (read.dbc) and stacking several years
#  - Writing to the columnar Parquet format (arrow)
#  - DuckDB connection with a memory limit and chained VIEWs
#  - Text standardization in SQL: upper, strip_accents, trim
#  - Classification with CASE WHEN (test -> pathogen group / subtype)
#  - Positive/negative result and aggregation by week
#  - Exporting the result with COPY ... TO ... (FORMAT PARQUET)
#
# TECHNIQUES / PACKAGES: RCurl, curl, read.dbc, arrow, duckdb, DBI, tidyverse
# NOTE: simulated data; variable, table and object names anonymized.
# =============================================================================

library(tidyverse)
library(arrow)
library(duckdb)
library(DBI)

set.seed(4096)

# -----------------------------------------------------------------------------
# 1. EXTRACT - batch download from the DATASUS FTP server
# -----------------------------------------------------------------------------
download_ftp <- function(ftp_url, file_pattern, dest_dir) {
  file_list <- RCurl::getURL(ftp_url, dirlistonly = TRUE)
  file_list <- strsplit(file_list, "\r?\n")[[1]]
  targets   <- file_list[grepl(file_pattern, file_list)]

  dir.create(dest_dir, showWarnings = FALSE)

  for (f in targets) {
    tryCatch({
      curl::curl_download(paste0(ftp_url, f), file.path(dest_dir, f))
      message("Downloaded: ", f)
    }, error = function(e) {
      message("Error in ", f, " -> ", e$message)
    })
  }
  invisible(targets)
}

# Example (SINAN notifiable diseases, final files):
# download_ftp("ftp://ftp.datasus.gov.br/dissemin/publicos/SINAN/DADOS/FINAIS/",
#              file_pattern = "^XXXXBR.*\\.dbc$", dest_dir = "raw_files")

# Read and stack the .dbc files -> Parquet
stack_dbc <- function(folder, output_file) {
  files <- list.files(folder, pattern = "\\.dbc$", full.names = TRUE)
  map(files, \(f) read.dbc::read.dbc(f) |> mutate(across(everything(), as.character))) |>
    bind_rows() |>
    write_parquet(output_file)
}
# stack_dbc("raw_files", "stacked_data.parquet")

# -----------------------------------------------------------------------------
# 2. Simulated lab test database (replaces the extraction step)
# -----------------------------------------------------------------------------
test_pool <- c("Respiratory molecular panel", "Rapid test influenza A and B",
               "SARS-CoV-2 PCR", "RSV molecular test", "Mycoplasma IgM",
               "Dengue IgM serology", "Chikungunya PCR", "Complete blood count")
detail_pool <- c("FLUA", "FLUAH3", "H1PDM09", "FLUB", "RSV", "COV2", "ADV", "RHINO", "")
result_pool <- c("Detected", "Not detected", "detected", "NOT DETECTED",
                 "Reactive", "Non-reactive", "Inconclusive")

n_ex <- 50000
tibble(
  sample_id = sprintf("A%07d", 1:n_ex),
  dt_coll   = as.Date("2022-01-01") + sample(0:1000, n_ex, TRUE),
  location  = sample(c(" são paulo", "Manaus ", "BELÉM", "recife", "Curitiba"), n_ex, TRUE),
  test      = sample(test_pool, n_ex, TRUE),
  test_det  = sample(detail_pool, n_ex, TRUE),
  result    = sample(result_pool, n_ex, TRUE)
) |>
  write_parquet("raw_tests.parquet")

# -----------------------------------------------------------------------------
# 3. TRANSFORM with DuckDB
# -----------------------------------------------------------------------------
con_db <- dbConnect(duckdb())
dbExecute(con_db, "SET memory_limit = '4GB'")

dbExecute(con_db, "
  CREATE OR REPLACE VIEW v0_raw AS
  SELECT * FROM read_parquet('raw_tests.parquet')
")

dbGetQuery(con_db, "DESCRIBE v0_raw")
dbGetQuery(con_db, "SELECT COUNT(*) AS n FROM v0_raw")

# 3.1 Text standardization (upper case, no accents, trimmed)
dbExecute(con_db, "
  CREATE OR REPLACE VIEW v1_std AS
  SELECT * REPLACE (
    trim(strip_accents(upper(location))) AS location,
    trim(strip_accents(upper(test)))     AS test,
    trim(strip_accents(upper(result)))   AS result
  )
  FROM v0_raw
")

dbGetQuery(con_db, "SELECT DISTINCT location FROM v1_std")

# 3.2 Classification of the test into a pathogen group
dbExecute(con_db, "
  CREATE OR REPLACE VIEW v2_group AS
  SELECT *,
    CASE
      WHEN test IN ('RESPIRATORY MOLECULAR PANEL', 'RAPID TEST INFLUENZA A AND B',
                    'SARS-COV-2 PCR', 'RSV MOLECULAR TEST', 'MYCOPLASMA IGM')
        THEN 'respiratory'
      WHEN test IN ('DENGUE IGM SEROLOGY', 'CHIKUNGUNYA PCR')
        THEN 'arbovirus'
      ELSE 'other'
    END AS pathogen_grp
  FROM v1_std
")

# 3.3 Subtype (for respiratory panels) and binary result
dbExecute(con_db, "
  CREATE OR REPLACE VIEW v3_subtype AS
  SELECT *,
    CASE
      WHEN pathogen_grp <> 'respiratory' THEN NULL
      WHEN test_det IN ('FLUA')     THEN 'influenza_a'
      WHEN test_det IN ('FLUAH3')   THEN 'influenza_a_h3n2'
      WHEN test_det IN ('H1PDM09')  THEN 'influenza_a_h1n1pdm09'
      WHEN test_det IN ('FLUB')     THEN 'influenza_b'
      WHEN test_det IN ('RSV')      THEN 'rsv'
      WHEN test_det IN ('COV2')     THEN 'sars_cov_2'
      WHEN test_det IN ('ADV')      THEN 'adenovirus'
      WHEN test_det IN ('RHINO')    THEN 'rhinovirus'
      ELSE 'unspecified'
    END AS subtype,
    CASE
      WHEN result IN ('DETECTED', 'REACTIVE')         THEN 1
      WHEN result IN ('NOT DETECTED', 'NON-REACTIVE') THEN 0
      ELSE NULL
    END AS positive
  FROM v2_group
")

dbGetQuery(con_db, "
  SELECT pathogen_grp, subtype, COUNT(*) AS n, AVG(positive) AS pos_rate
  FROM v3_subtype
  GROUP BY pathogen_grp, subtype
  ORDER BY pathogen_grp, n DESC
")

# 3.4 Weekly positivity
pos_week <- dbGetQuery(con_db, "
  SELECT pathogen_grp,
         date_trunc('week', dt_coll) AS week,
         COUNT(*)                    AS tests,
         SUM(positive)               AS positives,
         ROUND(100.0 * SUM(positive) / COUNT(positive), 2) AS pct_pos
  FROM v3_subtype
  WHERE pathogen_grp <> 'other'
  GROUP BY pathogen_grp, week
  ORDER BY week
")

ggplot(pos_week, aes(x = as.Date(week), y = pct_pos, color = pathogen_grp)) +
  geom_line() +
  labs(x = "Week", y = "% positivity", color = "Group") +
  theme_minimal()

# -----------------------------------------------------------------------------
# 4. LOAD - materialize and export to Parquet
# -----------------------------------------------------------------------------
dbExecute(con_db, "CREATE OR REPLACE TABLE tb_final AS SELECT * FROM v3_subtype")
dbGetQuery(con_db, "SELECT COUNT(*) AS n FROM tb_final")
dbExecute(con_db, "COPY tb_final TO 'standardized_tests.parquet' (FORMAT PARQUET)")

# Lazy reading of the result with arrow + dplyr
open_dataset("standardized_tests.parquet") |>
  filter(pathogen_grp == "respiratory") |>
  count(subtype) |>
  collect()

dbDisconnect(con_db, shutdown = TRUE)
