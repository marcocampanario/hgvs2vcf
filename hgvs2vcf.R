#!/usr/bin/env Rscript

# Terminal runner for hgvs_para_vcf.R.
# Keep both scripts in the same directory.
#
# Input:
#   Tab-separated table with exactly two columns:
#
#   SAMPLE    HGVS
#   sample1   NM_000059.4:c.7007G>A;p.(Arg2336His)
#   sample2   NM_007294.4:c.5266dup;p.(Gln1756ProfsTer74)
#
# Usage:
#   Rscript hgvs2vcf.R variantes.tsv resultado.tsv [GRCh38|GRCh37]
#
# Optional BED files:
#   Rscript hgvs2vcf.R variantes.tsv resultado.tsv GRCh38 \
#     --bed Twist=twist.bed \
#     --bed Agilent=agilent.bed.gz
#
# Reuse an already converted TSV without querying the API:
#   Rscript hgvs2vcf.R resultado.tsv anotado.tsv GRCh38 \
#     --table \
#     --bed Twist=twist.bed
#
# Optional:
#   --bed-mode variant|pos
#
# BED_<NAME> columns contain IN/OUT/PARTIAL/BORDER.
# NA is saved as an empty cell.
#
# All supplied BEDs must use the selected assembly.
# No liftover is performed.
#
# Install httr2 once before running.
# No packages are installed automatically.


main <- function(args = commandArgs(trailingOnly = TRUE),
                 script_path = NULL) {

  usage <- paste(
    "Uso:",
    "  Rscript hgvs2vcf.R entrada.tsv saida.tsv [GRCh38|GRCh37] [opcoes]",
    "",
    "Entrada:",
    "  TSV com duas colunas obrigatorias:",
    "    SAMPLE",
    "    HGVS",
    "",
    "Opcoes:",
    "  --bed NOME=arquivo.bed",
    "      Adiciona um kit; pode repetir quantas vezes quiser.",
    "",
    "  --bed arquivo.bed.gz",
    "      Usa o nome do arquivo como nome do kit.",
    "",
    "  --table",
    "      Entrada e um TSV ja convertido; nao consulta a API.",
    "",
    "  --bed-mode variant|pos",
    "      Regra de sobreposicao; padrao: variant.",
    "",
    "BEDs e coordenadas devem estar na montagem selecionada",
    "(padrao: GRCh38).",
    sep = "\n"
  )


  # ---------------------------------------------------------------------------
  # Parse command-line arguments
  # ---------------------------------------------------------------------------

  if (any(args %in% c("-h", "--help"))) {
    cat(usage, "\n")
    return(invisible(NULL))
  }

  positional <- character()
  beds <- character()

  table_mode <- FALSE
  bed_mode <- "variant"

  i <- 1L

  while (i <= length(args)) {

    arg <- args[[i]]

    if (arg == "--table") {

      table_mode <- TRUE

    } else if (arg %in% c("--bed", "--bed-mode")) {

      if (i == length(args) ||
          startsWith(args[[i + 1L]], "--")) {

        stop(
          "Falta o valor de ",
          arg,
          call. = FALSE
        )
      }

      i <- i + 1L
      value <- args[[i]]

      if (arg == "--bed-mode") {

        bed_mode <- match.arg(
          value,
          c("variant", "pos")
        )

      } else {

        equals <- regexpr(
          "=",
          value,
          fixed = TRUE
        )[[1L]]

        label <- if (equals > 0L) {
          substr(value, 1L, equals - 1L)
        } else {
          ""
        }

        path <- if (equals > 0L) {
          substring(value, equals + 1L)
        } else {
          value
        }

        if (!nzchar(path) ||
            (equals > 0L && !nzchar(trimws(label)))) {

          stop(
            "Use --bed NOME=arquivo.bed ou --bed arquivo.bed.",
            call. = FALSE
          )
        }

        beds <- c(
          beds,
          setNames(path, label)
        )
      }

    } else if (startsWith(arg, "--")) {

      stop(
        "Opcao desconhecida: ",
        arg,
        call. = FALSE
      )

    } else {

      positional <- c(
        positional,
        arg
      )
    }

    i <- i + 1L
  }


  # ---------------------------------------------------------------------------
  # Positional arguments
  # ---------------------------------------------------------------------------

  if (!length(positional) %in% c(2L, 3L)) {
    stop(usage, call. = FALSE)
  }

  input <- path.expand(
    positional[[1L]]
  )

  output <- path.expand(
    positional[[2L]]
  )

  assembly <- if (length(positional) == 3L) {
    positional[[3L]]
  } else {
    "GRCh38"
  }

  if (!assembly %in% c("GRCh38", "GRCh37")) {
    stop(
      "A montagem deve ser GRCh38 ou GRCh37.",
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Validate paths
  # ---------------------------------------------------------------------------

  if (!file.exists(input) ||
      dir.exists(input)) {

    stop(
      "Arquivo de entrada nao encontrado: ",
      input,
      call. = FALSE
    )
  }

  if (!grepl(
    "\\.tsv$",
    output,
    ignore.case = TRUE
  )) {

    stop(
      "O arquivo de saida deve ter extensao .tsv.",
      call. = FALSE
    )
  }

  link <- Sys.readlink(output)

  if (file.exists(output) ||
      (!is.na(link) && nzchar(link))) {

    stop(
      "A saida ja existe; sobrescrita recusada: ",
      output,
      call. = FALSE
    )
  }

  if (!dir.exists(dirname(output))) {

    stop(
      "A pasta de saida nao existe: ",
      dirname(output),
      call. = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Find and source hgvs_para_vcf.R
  # ---------------------------------------------------------------------------

  if (is.null(script_path)) {

    script_arg <- grep(
      "^--file=",
      commandArgs(trailingOnly = FALSE),
      value = TRUE
    )

    if (length(script_arg) != 1L) {
      stop(
        "Execute este script com Rscript.",
        call. = FALSE
      )
    }

    script_path <- sub(
      "^--file=",
      "",
      script_arg
    )
  }

  script <- normalizePath(
    script_path,
    mustWork = TRUE
  )

  helper <- file.path(
    dirname(script),
    "hgvs_para_vcf.R"
  )

  if (!file.exists(helper)) {

    stop(
      paste(
        "Coloque hgvs_para_vcf.R na mesma pasta",
        "de hgvs2vcf.R."
      ),
      call. = FALSE
    )
  }

  source(
    helper,
    local = TRUE,
    encoding = "UTF-8"
  )


  # ---------------------------------------------------------------------------
  # Prepare capture BED files
  # ---------------------------------------------------------------------------

  if (length(beds)) {

    beds <- preparar_beds_captura(
      beds,
      assembly = assembly
    )
  }


  # ---------------------------------------------------------------------------
  # Already-converted table mode
  # ---------------------------------------------------------------------------

  if (table_mode) {

    resultado <- utils::read.delim(
      input,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = "",
      fileEncoding = "UTF-8-BOM",
      na.strings = c("", "NA")
    )

    if (!nrow(resultado)) {

      stop(
        "A tabela de entrada esta vazia.",
        call. = FALSE
      )
    }


    # ---------------------------------------------------------------------------
    # SAMPLE + HGVS input mode
    # ---------------------------------------------------------------------------

  } else {

    entrada <- utils::read.delim(
      input,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = "",
      quote = "\"",
      fileEncoding = "UTF-8-BOM",
      na.strings = c("", "NA")
    )


    # -------------------------------------------------------------------------
    # Validate input table
    # -------------------------------------------------------------------------

    if (!nrow(entrada)) {

      stop(
        "A tabela de entrada esta vazia.",
        call. = FALSE
      )
    }

    required_columns <- c(
      "SAMPLE",
      "HGVS"
    )

    if (!identical(
      names(entrada),
      required_columns
    )) {

      stop(
        paste0(
          "A entrada deve conter exatamente duas colunas, nesta ordem: ",
          "SAMPLE e HGVS.\n",
          "Colunas encontradas: ",
          paste(names(entrada), collapse = ", ")
        ),
        call. = FALSE
      )
    }


    # SAMPLE is essential to trace the diagnostic variant back to the patient.
    sample_invalid <- (
      is.na(entrada$SAMPLE) |
        !nzchar(trimws(entrada$SAMPLE))
    )

    if (any(sample_invalid)) {

      stop(
        paste0(
          "Existem ",
          sum(sample_invalid),
          " linha(s) sem SAMPLE."
        ),
        call. = FALSE
      )
    }


    # Do not remove duplicated HGVS descriptions:
    # different patients may carry the same variant.
    samples <- as.character(
      entrada$SAMPLE
    )

    hgvs <- as.character(
      entrada$HGVS
    )


    # -------------------------------------------------------------------------
    # Initial HGVS screening
    # -------------------------------------------------------------------------

    triagem <- hgvs_para_vcf(
      hgvs,
      assembly = assembly,
      consultar = FALSE
    )

    message("Triagem das entradas:")

    print(
      table(
        triagem$STATUS,
        useNA = "ifany"
      )
    )


    # -------------------------------------------------------------------------
    # VariantValidator conversion
    # -------------------------------------------------------------------------

    resultado <- hgvs_para_vcf(
      hgvs,
      assembly = assembly
    )


    # -------------------------------------------------------------------------
    # Sanity check: one result row must correspond to one input row
    # -------------------------------------------------------------------------

    if (nrow(resultado) != length(samples)) {

      stop(
        paste0(
          "Erro interno: numero de resultados (",
          nrow(resultado),
          ") difere do numero de entradas (",
          length(samples),
          ")."
        ),
        call. = FALSE
      )
    }


    # -------------------------------------------------------------------------
    # Restore patient/sample information
    # -------------------------------------------------------------------------

    resultado <- data.frame(
      SAMPLE = samples,
      resultado,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }


  # ---------------------------------------------------------------------------
  # Capture BED annotation
  # ---------------------------------------------------------------------------

  if (length(beds)) {

    resultado <- anotar_beds_captura(
      resultado,
      beds,
      assembly = assembly,
      mode = bed_mode
    )
  }


  # ---------------------------------------------------------------------------
  # Save output
  # ---------------------------------------------------------------------------

  salvar_tsv_novo(
    resultado,
    output
  )


  # ---------------------------------------------------------------------------
  # Terminal summary
  # ---------------------------------------------------------------------------

  columns <- intersect(
    c(
      "SAMPLE",
      "HGVS_INPUT",
      "CHROM",
      "POS",
      "REF",
      "ALT",
      "STATUS",
      names(beds)
    ),
    names(resultado)
  )

  print(
    resultado[
      ,
      columns,
      drop = FALSE
    ],
    row.names = FALSE
  )

  message(
    "Tabela completa salva em: ",
    output
  )

  invisible(resultado)
}


if (sys.nframe() == 0L) {

  tryCatch(

    main(),

    error = function(e) {

      message(
        "ERRO: ",
        conditionMessage(e)
      )

      quit(
        save = "no",
        status = 1L
      )
    }
  )
}
