# HGVS -> genomic VCF fields, using R and the VariantValidator REST API.
# Assembly selected for this project: GRCh38.
#
# Usage in R:
#   install.packages("httr2")
#   source("hgvs_para_vcf.R")
#   triagem <- hgvs_para_vcf(hgvs_exemplo, consultar = FALSE)
#   resultado <- hgvs_para_vcf(hgvs_exemplo, assembly = "GRCh38")
#   resultado[, c("HGVS_INPUT", "CHROM", "POS", "REF", "ALT", "STATUS")]
#   salvar_tsv_novo(resultado, "hgvs_GRCh38.tsv")
#
# Optional capture kits (any number; all BEDs must use the selected assembly):
#   beds <- c(Twist = "twist.bed", Agilent = "agilent.bed.gz")
#   resultado <- hgvs_para_vcf(hgvs_exemplo, assembly = "GRCh38", beds = beds)
# Or annotate a previously converted table without making API requests:
#   resultado <- anotar_beds_captura(resultado, beds, assembly = "GRCh38")
# Each kit adds BED_<name>: IN, OUT, PARTIAL, BORDER, or NA (not evaluated).
# See anotar_beds_captura below for exact overlap rules.
#
# Alternatively: hgvs <- readLines("hgvs.txt", warn = FALSE, encoding = "UTF-8")
# One variant description per line. Input order and duplicates are retained.
# Sourcing this file does not send requests, read VCFs, or write any files.
#
# Scope:
# - Validates and maps versioned RefSeq HGVS descriptions to primary assembly loci.
# - Preserves the stated transcript version; never substitutes a default transcript.
# - Trims protein annotations for the DNA query; keeps the complete input string.
#   The protein annotation supplied by the user is NOT validated against the DNA.
# - Missing reference accessions and genomic deletions >=50 bp are flagged for
#   review. No transcript is inferred from a gene symbol or another input row.
# - Structural variants require an explicit representation and assembly-mapping
#   strategy. Four VCF fields alone may be insufficient (END/SVTYPE are needed).
# - API warnings and corrected HGVS descriptions produce REVIEW, not OK.
# - Results describe an equivalent VCF representation, not necessarily the exact
#   representation used by the original caller in a repetitive sequence.
# - Only the cleaned HGVS and requested assembly are sent to the public service.
# - Raw API responses and database metadata are retained in attributes(result).
#   saveRDS(result, "a_new_filename.rds") can preserve these attributes locally.
# - Live service availability and access requirements may change. If required,
#   set an authorized token with Sys.setenv(VV_TOKEN = "...") before running.
#
# API: https://rest.variantvalidator.org/
# Docs: https://openvar.github.io/variantValidator/rest-vv/rest_VariantValidator.html
# HGVS: https://hgvs-nomenclature.org/stable/recommendations/general/

hgvs_exemplo <- c(
  "NM_001042492.3:c.4279_4280dup;p.(Leu1427PhefsTer2)",
  "NM_001042492.3:c.3827G>A;p.(Arg1276Gln)",
  "NM_001042492.3:c.4677G>A;p.(Trp1559Ter)",
  "NM_001042492.3:c.2072T>C;p.(Leu691Pro)",
  "NM_006914.4:c.224T>G;p.(Met75Arg)",
  "NM_000170.3:c.2186del;p.(Ala729GlufsTer3)",
  "NC_000009.11:g.6642359_6701885del",
  "NM_003722.5:c.109C>T;p.(Arg37Ter)",
  "NM_001042492.3:c.3739_3742del;p.(Phe1247IlefsTer18)",
  "NM_001278116.2:c.2686C>T;p.(Gln896Ter)",
  "NM_004586.3:c.710C>T;p.(Pro237Leu)",
  "NM_001042492.3:c.2409+2T>G;p.(?)",
  "NM_003718.5:c.2579G>A;p.(Arg860Gln)",
  "NM_001042492.3:c.4382T>C;p.(Met1461Thr)",
  "NM_000271.5:c.3104C>T;p.(Ala1035Val)",
  "NM_000271.5:c.3213_3216dup;p.(Gly1073Ter)",
  "NM_005957.5:c.3G>C; p.(Met1?)",
  "NM_005957.5:c.1753-16G>A;p.(?)",
  "NC_000017.11:g.31350792_31355029del",
  "NM_001042492.3:c.7316del;p.(Leu2439Ter)",
  "NM_000038.6:c.1480dup;p.(Ser494LysfsTer43)",
  "NM_015557.3:c.5742G>C;p.(Gln1914His)",
  "NM_001042492.3:c.4532T>C;p.(Leu1511Pro)",
  "NM_000702.4:c.587G>A;p.(Arg196His)",
  "NM_005052.3:c.47A>C;p.(Lys16Thr)",
  "NM_013275.6:c.1491del;p.(Ser498AlafsTer12)",
  "c.2072T>C p.Leu691Pro",
  "NM_004168.4:c.1753C>T,NP_004159.2:p.Arg585Trp",
  "NM_031844.3:c.643_652del",
  "NM_003545.4:c.155A>G:p.Tyr52Cys",
  "NM_001042492.3 c.910C>T",
  "NM_001042492.3 c.1110del",
  "PTPN11:c.236A>G chr12-112450416 A>G p.Gln79Arg",
  "NM_001032221.6:c.751G>C",
  "FBN1:c.4271dup",
  "NM_001042492.3:c.4332+1G>T,",
  "NM_005559.4:c.4993dup,NP_005550.2:p.Gln1665ProfsTer19",
  "NM_031263.4:c.571A>G:p.Arg191Gly",
  "NM_006180.6 c.526T>C p.(Cys176Arg)"
)

.scalar <- function(x) {
  if (is.null(x) || !is.atomic(x) || length(x) == 0L) return(NA_character_)
  as.character(x[[1L]])
}

.clean_hgvs <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) {
    return(list(status = "EMPTY_INPUT", note = "Descricao vazia."))
  }
  x <- gsub("\\", "", x, fixed = TRUE)
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  x <- gsub("<[^>]*>", " ", x)
  x <- trimws(gsub("^\\s*\\||\\|\\s*$", "", x, perl = TRUE))
  dna_mentions <- regmatches(x, gregexpr("[cgnm]\\.[0-9*(?\\[]", x, perl = TRUE))[[1L]]
  if (length(dna_mentions) > 1L) {
    return(list(status = "MULTIPLE_DESCRIPTIONS", note = "Use uma variante de DNA por linha."))
  }
  match <- regexec(
    "((?:NM|NR|NC|NG|NT|NW)_[0-9]+\\.[0-9]+)\\s*:?\\s*([cgnm]\\.\\S+)",
    x, perl = TRUE
  )
  pieces <- regmatches(x, match)[[1L]]
  if (length(pieces) != 3L) {
    return(list(status = "MISSING_REFERENCE", note = "Informe accession.version, como NM_001042492.3. Nao inferido automaticamente."))
  }
  accession <- pieces[[2L]]
  change <- pieces[[3L]]
  if (grepl("[\\[\\]]", change, perl = TRUE)) {
    return(list(status = "COMPLEX_HGVS_REVIEW", note = "Alelos ou descricoes complexas exigem tratamento especifico."))
  }
  change <- sub("[;,].*$", "", change)
  change <- sub(":p\\..*$", "", change)
  hgvs <- paste0(accession, ":", change)
  info <- list(status = "READY", note = NA_character_, query = hgvs,
               accession = accession, change = change,
               source_start = NA_integer_, source_end = NA_integer_, source_build = NA_character_)
  deletion <- regmatches(change, regexec("^g\\.([0-9]+)(?:_([0-9]+))?del(?:[ACGT]*)$", change, perl = TRUE))[[1L]]
  if (length(deletion) == 3L) {
    start <- as.integer(deletion[[2L]])
    end <- if (nzchar(deletion[[3L]])) as.integer(deletion[[3L]]) else start
    if (is.na(start) || is.na(end) || start < 1L || end < start) {
      info$status <- "INVALID_INTERVAL"
      info$note <- "Intervalo de delecao invalido."
      return(info)
    }
    info$source_start <- start
    info$source_end <- end
    # Explicitly documented reference accessions occurring in the supplied list.
    known_builds <- c("NC_000009.11" = "GRCh37", "NC_000017.11" = "GRCh38")
    if (accession %in% names(known_builds)) info$source_build <- unname(known_builds[[accession]])
    if (end - start + 1L >= 50L) {
      info$status <- "SV_REVIEW"
      info$note <- paste0("Delecao genomica de ", end - start + 1L,
                          " bp. Definir representacao SV/END e validar remapeamento para a montagem alvo.")
    }
  }
  info
}

.request_vv <- function(hgvs, assembly, base_url, token) {
  accession <- sub(":.*$", "", hgvs)
  select <- if (grepl("^(NM|NR)_", accession)) accession else "mane_select"
  # Encode each path component, including '+' and '>', to preserve intronic HGVS.
  components <- c("VariantValidator", "variantvalidator", assembly, hgvs, select)
  encoded <- vapply(components, utils::URLencode, character(1L), reserved = TRUE)
  url <- paste0(sub("/+$", "", base_url), "/", paste(encoded, collapse = "/"))
  req <- httr2::request(url)
  req <- httr2::req_url_query(req, `content-type` = "application/json")
  req <- httr2::req_headers(req, Accept = "application/json")
  req <- httr2::req_user_agent(req, "hgvs-para-vcf-R/1.0")
  req <- httr2::req_timeout(req, 90)
  req <- httr2::req_retry(req, max_tries = 3L)
  if (nzchar(token)) req <- httr2::req_auth_bearer_token(req, token)
  httr2::resp_body_json(httr2::req_perform(req), simplifyVector = FALSE)
}

.parse_vv <- function(data, hgvs, assembly, use_chr) {
  nodes <- Filter(is.list, data[setdiff(names(data), c("metadata", "flag"))])
  warning_text <- unique(unlist(lapply(nodes, function(node) node$validation_warnings), use.names = FALSE))
  warning_text <- warning_text[!is.na(warning_text) & nzchar(warning_text)]
  warnings <- if (length(warning_text)) paste(warning_text, collapse = " | ") else NA_character_
  failure <- function(status, note) list(STATUS = status, WARNINGS = paste(na.omit(c(note, warnings)), collapse = " | "))
  # Never take an unrelated first result from a multi-record JSON response.
  candidates <- Filter(function(node) identical(.scalar(node$submitted_variant), hgvs) &&
                         is.list(node$primary_assembly_loci[[tolower(assembly)]]), nodes)
  if (!length(candidates)) {
    return(failure("NO_TARGET_MAPPING", "Sem mapeamento explicito para a montagem alvo; revisar resposta e avisos da API."))
  }
  accession <- sub(":.*$", "", hgvs)
  if (grepl("^(NM|NR)_", accession)) {
    candidates <- Filter(function(node) {
      returned <- .scalar(node$hgvs_transcript_variant)
      !is.na(returned) && identical(sub(":.*$", "", returned), accession)
    }, candidates)
    if (!length(candidates)) return(failure("TRANSCRIPT_VERSION_MISMATCH", "Transcrito ou versao original nao preservado na resposta."))
  }
  signatures <- vapply(candidates, function(node) {
    v <- node$primary_assembly_loci[[tolower(assembly)]]$vcf
    paste(v$chr, v$pos, v$ref, v$alt, sep = ":")
  }, character(1L))
  if (length(unique(signatures)) != 1L) return(failure("AMBIGUOUS_MAPPING", "A API retornou representacoes genomicas diferentes."))
  node <- candidates[[1L]]
  locus <- node$primary_assembly_loci[[tolower(assembly)]]
  vcf <- locus$vcf
  chromosome <- .scalar(vcf$chr)
  pos <- suppressWarnings(as.integer(.scalar(vcf$pos)))
  ref <- .scalar(vcf$ref)
  alt <- .scalar(vcf$alt)
  if (anyNA(c(chromosome, pos, ref, alt)) || pos < 1L ||
      !grepl("^[ACGTN]+$", ref) || !grepl("^[ACGTN]+$", alt)) {
    return(failure("UNSUPPORTED_VCF_ALLELES", "Campos VCF ausentes ou alelos nao explicitos."))
  }
  chromosome <- sub("^chr", "", chromosome)
  if (!chromosome %in% c(as.character(1:22), "X", "Y", "M", "MT")) {
    return(failure("NON_PRIMARY_CONTIG", "Contig fora dos cromossomos primarios esperados."))
  }
  if (use_chr) chromosome <- paste0("chr", if (chromosome %in% c("M", "MT")) "M" else chromosome)
  normalized <- .scalar(node$hgvs_transcript_variant)
  genomic <- .scalar(locus$hgvs_genomic_description)
  changed <- if (grepl("^(NM|NR)_", accession)) !identical(normalized, hgvs) else !identical(genomic, hgvs)
  if (changed) warnings <- paste(na.omit(c(warnings, "HGVS retornado difere da entrada; revisar normalizacao/remapeamento.")), collapse = " | ")
  if (grepl("N", paste0(ref, alt))) warnings <- paste(na.omit(c(warnings, "Alelo contem base N.")), collapse = " | ")
  list(CHROM = chromosome, POS = pos, REF = ref, ALT = alt,
       GENE = .scalar(node$gene_symbol), HGVS_C = normalized, HGVS_G = genomic,
       HGVS_P_RETURNED = .scalar(node$hgvs_predicted_protein_consequence$tlr),
       STATUS = if (!is.na(warnings) && nzchar(warnings)) "REVIEW" else "OK",
       WARNINGS = warnings)
}

hgvs_para_vcf <- function(hgvs, assembly = c("GRCh38", "GRCh37"), use_chr = TRUE,
                          consultar = TRUE, base_url = "https://rest.variantvalidator.org",
                          token = Sys.getenv("VV_TOKEN", unset = ""),
                          beds = NULL, bed_mode = c("variant", "pos"), harmonize_chr = TRUE) {
  assembly <- match.arg(assembly)
  bed_mode <- match.arg(bed_mode)
  stopifnot(is.character(hgvs), length(hgvs) > 0L, length(use_chr) == 1L,
            !is.na(use_chr), is.logical(use_chr))
  if (consultar && !requireNamespace("httr2", quietly = TRUE)) {
    stop("Instale o pacote httr2: install.packages('httr2').", call. = FALSE)
  }
  if (length(beds)) beds <- preparar_beds_captura(beds, assembly, harmonize_chr)
  rows <- vector("list", length(hgvs))
  responses <- list()
  for (i in seq_along(hgvs)) {
    info <- .clean_hgvs(hgvs[[i]])
    row <- data.frame(
      ROW_ID = i, HGVS_INPUT = hgvs[[i]], HGVS_QUERY = .scalar(info$query),
      ASSEMBLY = assembly, CHROM = NA_character_, POS = NA_integer_,
      REF = NA_character_, ALT = NA_character_, GENE = NA_character_,
      HGVS_C = NA_character_, HGVS_G = NA_character_, HGVS_P_RETURNED = NA_character_,
      STATUS = info$status, WARNINGS = .scalar(info$note),
      SOURCE_REFERENCE = .scalar(info$accession),
      SOURCE_START = if (is.null(info$source_start)) NA_integer_ else info$source_start,
      SOURCE_END = if (is.null(info$source_end)) NA_integer_ else info$source_end,
      SOURCE_BUILD = .scalar(info$source_build), stringsAsFactors = FALSE
    )
    if (identical(info$status, "READY") && consultar) {
      query <- info$query
      if (is.null(responses[[query]])) {
        message(sprintf("Consultando %d/%d: %s", i, length(hgvs), query))
        Sys.sleep(0.4) # Below the public API's recommended 3 requests/second.
        responses[[query]] <- tryCatch(.request_vv(query, assembly, base_url, token), error = identity)
      }
      data <- responses[[query]]
      if (inherits(data, "error")) {
        if (inherits(data, c("httr2_http_401", "httr2_http_403"))) {
          stop("A API recusou acesso. Verifique as credenciais autorizadas e os requisitos do servico.", call. = FALSE)
        }
        row$STATUS <- "API_ERROR"
        row$WARNINGS <- conditionMessage(data)
      } else {
        parsed <- tryCatch(.parse_vv(data, query, assembly, use_chr), error = identity)
        if (inherits(parsed, "error")) {
          row$STATUS <- "RESPONSE_ERROR"
          row$WARNINGS <- conditionMessage(parsed)
        } else {
          for (name in names(parsed)) row[[name]] <- parsed[[name]]
        }
      }
    }
    rows[[i]] <- row
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  attr(result, "api_responses") <- responses
  attr(result, "assembly") <- assembly
  attr(result, "retrieval_time_utc") <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  if (length(beds)) {
    result <- anotar_beds_captura(result, beds, assembly = assembly,
                                  mode = bed_mode, harmonize_chr = harmonize_chr)
  }
  result
}

# Capture BED support uses base R only. No VCFs or BEDs are modified.
# BED coordinates are 0-based, half-open; POS is 1-based.
.capture_chr <- function(x, harmonize_chr = TRUE) {
  x <- trimws(as.character(x))
  if (!harmonize_chr) return(x)
  primary <- !is.na(x) & grepl("^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", x)
  x[primary] <- sub("^chr", "", x[primary])
  x[!is.na(x) & x %in% c("M", "MT")] <- "MT"
  x
}

.read_capture_bed <- function(path, harmonize_chr) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection), add = TRUE)
  lines <- trimws(sub("^\ufeff", "", readLines(connection, warn = FALSE)))
  lines <- lines[nzchar(lines) & !grepl("^(#|track([[:space:]]|$)|browser([[:space:]]|$))", lines)]
  if (!length(lines)) stop("BED sem intervalos: ", path, call. = FALSE)
  fields <- strsplit(lines, "[[:space:]]+", perl = TRUE)
  if (any(lengths(fields) < 3L)) stop("BED com menos de tres colunas: ", path, call. = FALSE)
  if (any(lengths(fields) >= 12L)) {
    stop("BED12 nao suportado. Forneca os alvos individuais em BED3, sem intervalos que unam blocos: ", path, call. = FALSE)
  }
  chrom <- .capture_chr(vapply(fields, `[[`, character(1L), 1L), harmonize_chr)
  start <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1L), 2L)))
  end <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1L), 3L)))
  if (any(!is.finite(start) | !is.finite(end) | start < 0 | end <= start |
          start != floor(start) | end != floor(end)) || any(chrom == ".")) {
    stop("BED com coordenadas invalidas (START >= 0, END > START, inteiros): ", path, call. = FALSE)
  }
  # Sort and merge overlapping/adjacent targets so coverage cannot be double-counted.
  lapply(split(seq_along(chrom), chrom), function(index) {
    index <- index[order(start[index], end[index])]
    left <- start[index]
    right <- cummax(end[index])
    group <- cumsum(left > c(-Inf, head(right, -1L)))
    list(start = left[!duplicated(group)], end = right[!duplicated(group, fromLast = TRUE)])
  })
}

preparar_beds_captura <- function(beds, assembly = c("GRCh38", "GRCh37"), harmonize_chr = TRUE) {
  assembly <- match.arg(assembly)
  stopifnot(is.logical(harmonize_chr), length(harmonize_chr) == 1L, !is.na(harmonize_chr))
  if (inherits(beds, "capture_beds")) {
    if (!identical(attr(beds, "assembly"), assembly) ||
        !identical(attr(beds, "harmonize_chr"), harmonize_chr)) {
      stop("Os BEDs preparados usam outra montagem ou regra de cromossomos.", call. = FALSE)
    }
    return(beds)
  }
  if (!is.character(beds) || !length(beds) || anyNA(beds) || any(!nzchar(beds))) {
    stop("beds deve ser um vetor de caminhos, por exemplo c(Twist = 'twist.bed').", call. = FALSE)
  }
  labels <- names(beds)
  if (is.null(labels)) labels <- rep("", length(beds))
  missing_label <- is.na(labels) | !nzchar(trimws(labels))
  labels[missing_label] <- sub("\\.bed(\\.gz)?$", "", basename(beds[missing_label]), ignore.case = TRUE)
  labels <- paste0("BED_", make.names(labels, unique = FALSE))
  if (anyDuplicated(labels)) stop("Nomes de kits repetidos; use nomes distintos no vetor beds.", call. = FALSE)
  paths <- path.expand(unname(beds))
  if (any(!file.exists(paths) | dir.exists(paths))) {
    stop("BED inexistente ou caminho de diretorio: ", paste(paths[!file.exists(paths) | dir.exists(paths)], collapse = ", "), call. = FALSE)
  }
  paths <- normalizePath(paths, mustWork = TRUE)
  regions <- lapply(paths, .read_capture_bed, harmonize_chr = harmonize_chr)
  names(regions) <- labels
  structure(regions, class = c("capture_beds", "list"), assembly = assembly,
            harmonize_chr = harmonize_chr, files = setNames(paths, labels))
}

.capture_covered_bases <- function(regions, start, end) {
  if (is.null(regions) || end <= start) return(0)
  i <- findInterval(start, regions$end) + 1L
  covered <- 0
  while (i <= length(regions$start) && regions$start[[i]] < end) {
    covered <- covered + max(0, min(end, regions$end[[i]]) - max(start, regions$start[[i]]))
    i <- i + 1L
  }
  covered
}

.capture_span <- function(pos, ref, alt) {
  ref <- toupper(as.character(ref))
  alt <- toupper(as.character(alt))
  if (is.na(ref) || is.na(alt) || !grepl("^[ACGT]+$", ref) ||
      !grepl("^[ACGT]+$", alt) || ref == alt) return(NULL)
  prefix <- 0L
  limit <- min(nchar(ref), nchar(alt))
  while (prefix < limit && substr(ref, prefix + 1L, prefix + 1L) == substr(alt, prefix + 1L, prefix + 1L)) {
    prefix <- prefix + 1L
  }
  ref <- substring(ref, prefix + 1L)
  alt <- substring(alt, prefix + 1L)
  while (nzchar(ref) && nzchar(alt) &&
         substr(ref, nchar(ref), nchar(ref)) == substr(alt, nchar(alt), nchar(alt))) {
    ref <- substr(ref, 1L, nchar(ref) - 1L)
    alt <- substr(alt, 1L, nchar(alt) - 1L)
  }
  start <- pos - 1 + prefix
  c(start = start, end = start + nchar(ref))
}

# Adds one BED_<kit> column per file, preserving input rows and conversion warnings.
# - mode='variant' (default): excludes unchanged VCF prefix/suffix anchors.
#   IN: entire remaining REF span is inside the union of targets; OUT: zero overlap;
#   PARTIAL: only part of the span overlaps. For a pure insertion, the breakpoint
#   has zero width: IN requires both adjacent reference bases in the targets;
#   BORDER means only one flank is in a target; OUT means neither is.
# - mode='pos': checks only the VCF POS base; outputs IN/OUT/NA. This can differ
#   from affected-sequence overlap for indels, because POS may be an anchor.
# - NA: missing/invalid coordinates, unusable STATUS, or unsupported alleles.
#   If STATUS is present, only OK/REVIEW rows are eligible. REVIEW stays REVIEW.
# - Literal A/C/G/T alleles only in variant mode; symbolic SVs are not evaluated.
# - Primary contig naming is harmonized in memory: chr1/1, chrX/X, chrM/M/MT.
#   A valid chromosome with no intervals in a kit is OUT. Other contig names are
#   matched exactly. Harmonizing names does NOT convert genome assemblies.
# - BEDs must contain actual capture intervals (BED3 through BED9); BED12 blocks
#   are rejected. Extra columns are ignored. No padding is added.
# - assembly DECLARES the build of all BEDs; it cannot be inferred from BED3.
#   Existing ASSEMBLY values in the table are checked. No liftover is performed.
anotar_beds_captura <- function(resultado, beds, assembly = c("GRCh38", "GRCh37"),
                                mode = c("variant", "pos"), harmonize_chr = TRUE) {
  assembly <- match.arg(assembly)
  mode <- match.arg(mode)
  required <- if (mode == "variant") c("CHROM", "POS", "REF", "ALT") else c("CHROM", "POS")
  if (!is.data.frame(resultado) || !all(required %in% names(resultado))) {
    stop("A tabela deve conter: ", paste(required, collapse = ", "), call. = FALSE)
  }
  prepared <- preparar_beds_captura(beds, assembly, harmonize_chr)
  if (any(names(prepared) %in% names(resultado))) {
    stop("A tabela ja possui coluna(s) para esse(s) kit(s). Use nomes novos ou remova explicitamente essas colunas.", call. = FALSE)
  }
  if ("ASSEMBLY" %in% names(resultado)) {
    declared <- as.character(resultado$ASSEMBLY)
    if (any(!is.na(declared) & nzchar(declared) & declared != assembly)) {
      stop("A montagem da tabela difere da montagem declarada para os BEDs.", call. = FALSE)
    }
  }
  chrom <- .capture_chr(resultado$CHROM, harmonize_chr)
  pos <- suppressWarnings(as.numeric(as.character(resultado$POS)))
  valid <- !is.na(chrom) & nzchar(chrom) & chrom != "." & is.finite(pos) & pos >= 1 & pos == floor(pos)
  if ("STATUS" %in% names(resultado)) valid <- valid & resultado$STATUS %in% c("OK", "REVIEW")
  if ("ASSEMBLY" %in% names(resultado)) valid <- valid & !is.na(declared) & declared == assembly
  valid[is.na(valid)] <- FALSE
  spans <- lapply(seq_len(nrow(resultado)), function(i) {
    if (!valid[[i]]) return(NULL)
    if (mode == "pos") return(c(start = pos[[i]] - 1, end = pos[[i]]))
    .capture_span(pos[[i]], resultado$REF[[i]], resultado$ALT[[i]])
  })
  for (kit in names(prepared)) {
    regions <- prepared[[kit]]
    result <- rep(NA_character_, nrow(resultado))
    for (i in which(lengths(spans) == 2L)) {
      span <- spans[[i]]
      start <- unname(span[[1L]])
      end <- unname(span[[2L]])
      intervals <- regions[[chrom[[i]]]]
      if (start == end) {
        left <- start > 0 && .capture_covered_bases(intervals, start - 1, start) > 0
        right <- .capture_covered_bases(intervals, start, start + 1) > 0
        result[[i]] <- if (left && right) "IN" else if (left || right) "BORDER" else "OUT"
      } else {
        overlap <- .capture_covered_bases(intervals, start, end)
        result[[i]] <- if (overlap == end - start) "IN" else if (overlap > 0) "PARTIAL" else "OUT"
      }
    }
    resultado[[kit]] <- result
  }
  attr(resultado, "capture_beds") <- list(files = attr(prepared, "files"), assembly = assembly,
                                          mode = mode, harmonize_chr = harmonize_chr)
  resultado
}

salvar_tsv_novo <- function(x, path) {
  if (!grepl("\\.tsv$", path, ignore.case = TRUE)) stop("Use um NOVO arquivo .tsv.")
  link <- Sys.readlink(path)
  if (file.exists(path) || (!is.na(link) && nzchar(link))) stop("Destino ja existe; sobrescrita recusada.")
  # TSV is for the flat table; saveRDS on a new path preserves raw-response attributes.
  utils::write.table(x, file = path, sep = "\t", quote = TRUE,
                     row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}
