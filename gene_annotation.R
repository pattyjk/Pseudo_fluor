#!/usr/bin/env Rscript
# =============================================================================
# BLASTx Best-Hit Selection + GO Term Annotation via UniProt REST API
# =============================================================================
# Input : BLASTx vs SwissProt results in tabular format (-outfmt 6)
#         Significant genes CSV (with "Gene" column of group_XXXX IDs)
# Output: (1) best_hits.tsv             — one best hit per query
#         (2) go_results.tsv            — GO term annotations per query
#         (3) significant_annotated.tsv — sig genes with GO annotations
# =============================================================================

# ── 0. Install / load packages ───────────────────────────────────────────────
packages <- c("dplyr", "stringr", "tidyr", "readr", "progress", "httr", "jsonlite")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# ── 1. User settings ─────────────────────────────────────────────────────────
BLAST_FILE   <- "~/Documents/GitHub/Pseudo_fluor/blastx_swissprot_results.tsv"
SIG_FILE     <- "~/Documents/GitHub/Pseudo_fluor/Figures/09_conness_fisher_significant.csv"
OUT_HITS     <- "~/Documents/GitHub/Pseudo_fluor/best_hits.tsv"
OUT_GO       <- "~/Documents/GitHub/Pseudo_fluor/go_results.tsv"
OUT_SIG      <- "~/Documents/GitHub/Pseudo_fluor/significant_annotated.tsv"

MIN_PIDENT   <- 30    # minimum % identity
MIN_LENGTH   <- 50    # minimum alignment length
MAX_EVALUE   <- 1e-5  # maximum e-value
BATCH_SIZE   <- 50    # IDs per UniProt API request

# ── 2. Load BLAST results ────────────────────────────────────────────────────
col_names <- c(
  "qseqid", "sseqid", "stitle",
  "pident", "length", "mismatch",
  "gapopen", "qstart", "qend",
  "sstart", "send", "evalue", "bitscore"
)

blast <- read_tsv(
  BLAST_FILE,
  col_names  = col_names,
  col_types  = cols(
    qseqid   = col_character(),
    sseqid   = col_character(),
    stitle   = col_character(),
    pident   = col_double(),
    length   = col_integer(),
    mismatch = col_integer(),
    gapopen  = col_integer(),
    qstart   = col_integer(),
    qend     = col_integer(),
    sstart   = col_integer(),
    send     = col_integer(),
    evalue   = col_double(),
    bitscore = col_double()
  ),
  show_col_types = FALSE
)

message(sprintf("  Loaded %d hits for %d queries.", nrow(blast), n_distinct(blast$qseqid)))

# ── 3. Filter and select best hit ────────────────────────────────────────────
best_hits <- blast %>%
  filter(pident >= MIN_PIDENT, length >= MIN_LENGTH, evalue <= MAX_EVALUE) %>%
  arrange(qseqid, evalue, desc(bitscore), desc(pident)) %>%
  group_by(qseqid) %>%
  slice(1) %>%
  ungroup()

message(sprintf(
  "  %d queries kept after filtering (from %d total).",
  nrow(best_hits), n_distinct(blast$qseqid)
))

# ── 4. Parse UniProt accessions ──────────────────────────────────────────────
best_hits <- best_hits %>%
  mutate(
    uniprot_acc  = str_extract(sseqid, "(?<=sp\\|)[A-Z0-9]+"),
    uniprot_full = str_extract(sseqid, "(?<=sp\\|)[A-Z0-9.]+")
  )

message(sprintf(
  "  Parsed %d UniProt accessions (%d unique).",
  sum(!is.na(best_hits$uniprot_acc)),
  n_distinct(best_hits$uniprot_acc, na.rm = TRUE)
))

write_tsv(best_hits, OUT_HITS)
message(sprintf("  Best hits written to '%s'.", OUT_HITS))

# ── 5. Fetch GO terms from UniProt REST API ───────────────────────────────────
# Uses GET with fixed size=500 to avoid 400 errors from dynamic size parameter
fetch_uniprot_go <- function(ids) {
  query <- paste0("accession:", ids, collapse = " OR ")
  
  res <- GET(
    "https://rest.uniprot.org/uniprotkb/search",
    query = list(
      query  = query,
      fields = "accession,gene_names,protein_name,go_id,go,organism_name",
      format = "json",
      size   = 500        # fixed — must not equal length(ids)
    ),
    add_headers(Accept = "application/json")
  )
  
  if (http_error(res)) {
    stop(sprintf("UniProt API error: %s", http_status(res)$message))
  }
  
  content(res, as = "text", encoding = "UTF-8") %>%
    fromJSON(flatten = TRUE) %>%
    .[["results"]]
}

uniprot_ids <- unique(na.omit(best_hits$uniprot_acc))
batches     <- split(uniprot_ids, ceiling(seq_along(uniprot_ids) / BATCH_SIZE))

message(sprintf(
  "  Fetching GO terms for %d UniProt IDs in %d batches...",
  length(uniprot_ids), length(batches)
))

pb <- progress_bar$new(
  format = "  UniProt [:bar] :current/:total batches  ETA: :eta",
  total  = length(batches), clear = FALSE, width = 60
)

go_raw_list <- vector("list", length(batches))

for (i in seq_along(batches)) {
  pb$tick()
  tryCatch({
    go_raw_list[[i]] <- fetch_uniprot_go(batches[[i]])
  }, error = function(e) {
    message(sprintf("\n  Warning: batch %d failed — %s", i, conditionMessage(e)))
  })
  Sys.sleep(0.5)
}

# ── 6. Parse GO terms from API response ──────────────────────────────────────
parse_go_terms <- function(result_row) {
  acc  <- result_row$primaryAccession
  refs <- result_row$uniProtKBCrossReferences[[1]]
  
  if (is.null(refs) || nrow(refs) == 0) {
    return(tibble(
      uniprot_acc = acc,
      go_id       = NA_character_,
      go_name     = NA_character_,
      go_aspect   = NA_character_
    ))
  }
  
  go_refs <- refs %>% filter(database == "GO")
  
  if (nrow(go_refs) == 0) {
    return(tibble(
      uniprot_acc = acc,
      go_id       = NA_character_,
      go_name     = NA_character_,
      go_aspect   = NA_character_
    ))
  }
  
  tibble(
    uniprot_acc = acc,
    go_id       = go_refs$id,
    go_name     = vapply(go_refs$properties, function(p) {
      p$value[p$key == "GoTerm"][1]
    }, character(1)),
    go_aspect   = vapply(go_refs$properties, function(p) {
      p$value[p$key == "GoEvidenceType"][1]
    }, character(1))
  )
}

go_df <- bind_rows(go_raw_list) %>%
  split(seq_len(nrow(.))) %>%
  lapply(parse_go_terms) %>%
  bind_rows()



# ── 7. Assemble final annotated table ────────────────────────────────────────
final <- best_hits %>%
  dplyr::select(qseqid, sseqid, uniprot_acc, stitle, pident, length, evalue, bitscore) %>%
  left_join(go_df, by = "uniprot_acc")

message(sprintf(
  "  %d queries have at least one GO term; %d have no GO annotation.",
  n_distinct(filter(final, !is.na(go_id))$qseqid),
  n_distinct(filter(final,  is.na(go_id))$qseqid)
))
#4308 queries have at least one GO term; 194 have no GO annotation.

write_tsv(final, OUT_GO)

#Map back to significant genes
sig_genes <- read_csv(SIG_FILE, show_col_types = FALSE)

final_sig <- sig_genes %>%
  left_join(final, by = c("Gene" = "qseqid"))

write_tsv(final_sig, OUT_SIG)

