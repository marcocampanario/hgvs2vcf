# hgvs2vcf

Convert and validate HGVS variant descriptions into genomic VCF coordinates while preserving sample-level information.

`hgvs2vcf` is an R-based workflow designed to facilitate the validation of reported diagnostic or clinically relevant variants against genomic VCF data. It converts versioned HGVS descriptions into genomic coordinates (`CHROM`, `POS`, `REF`, and `ALT`) using VariantValidator and can optionally determine whether each variant is covered by one or more exome capture BED files.

The workflow supports **GRCh38** and **GRCh37**.

## Features

- Converts transcript-level HGVS descriptions into genomic VCF coordinates.
- Preserves the original `SAMPLE` associated with each variant.
- Preserves transcript versions provided by the user.
- Supports GRCh38 and GRCh37.
- Returns chromosome, position, reference allele, alternate allele, gene, and genomic/transcript HGVS representations.
- Flags ambiguous, unsupported, or potentially problematic HGVS descriptions for review.
- Supports annotation against **multiple exome capture BED files**.
- Distinguishes variants that are completely, partially, or not covered by a capture kit.
- Supports both SNVs and small INDELs.
- Can annotate an already converted table without repeating API requests.
- Preserves input order and duplicated variants occurring in different samples.

## Repository structure

```text
hgvs2vcf/
├── hgvs2vcf.R
├── hgvs_para_vcf.R
├── hg38_Twist_Bioscience_for_Illumina_Exome_2_5.bed
└── README.md
```

- `hgvs2vcf.R` — command-line interface and main workflow.
- `hgvs_para_vcf.R` — HGVS parsing, VariantValidator querying, VCF conversion, and BED annotation functions.
- `hg38_Twist_Bioscience_for_Illumina_Exome_2_5.bed` — example GRCh38 exome capture BED file.

## Requirements

The workflow requires R and the `httr2` package.

Install the dependency once with:

```r
install.packages("httr2")
```

No packages are installed automatically by the scripts.

## Installation

Clone the repository:

```bash
git clone git@github.com:marcocampanario/hgvs2vcf.git
cd hgvs2vcf
```

The scripts `hgvs2vcf.R` and `hgvs_para_vcf.R` must remain in the same directory.

## Input

The standard input is a tab-separated file containing **exactly two columns**, in this order:

```text
SAMPLE	HGVS
sample1	NM_000059.4:c.7007G>A;p.(Arg2336His)
sample2	NM_007294.4:c.5266dup;p.(Gln1756ProfsTer74)
sample3	NM_000546.6:c.1010G>A;p.(Arg337His)
```

The `SAMPLE` column allows each reported variant to remain linked to the individual in whom it was identified.

The HGVS description should contain a **versioned reference accession**, for example:

```text
NM_000059.4:c.7007G>A
```

Protein annotations may also be included:

```text
NM_000059.4:c.7007G>A;p.(Arg2336His)
```

The DNA-level HGVS is used for conversion. The protein annotation supplied in the input is retained as part of the original description but is **not independently validated against the DNA change**.

## Basic usage

The general command is:

```bash
Rscript hgvs2vcf.R input.tsv output.tsv [GRCh38|GRCh37]
```

For example:

```bash
Rscript hgvs2vcf.R variants.tsv variants_GRCh38.tsv GRCh38
```

If no assembly is provided, the default is:

```text
GRCh38
```

The script first screens the HGVS descriptions and then queries VariantValidator for variants that can be processed.

## Output

The output is a tab-separated table containing the original sample and HGVS information together with the converted genomic representation.

Main columns include:

| Column | Description |
|---|---|
| `SAMPLE` | Sample or patient identifier from the input |
| `ROW_ID` | Original row number |
| `HGVS_INPUT` | Original HGVS description |
| `HGVS_QUERY` | Cleaned DNA-level HGVS submitted for validation |
| `ASSEMBLY` | Target genome assembly |
| `CHROM` | Genomic chromosome |
| `POS` | 1-based VCF position |
| `REF` | VCF reference allele |
| `ALT` | VCF alternate allele |
| `GENE` | Gene symbol returned by VariantValidator |
| `HGVS_C` | Transcript-level HGVS returned by VariantValidator |
| `HGVS_G` | Genomic HGVS returned by VariantValidator |
| `HGVS_P_RETURNED` | Predicted protein consequence returned by VariantValidator |
| `STATUS` | Conversion/validation status |
| `WARNINGS` | Warnings or additional information |

An example successful result may look like:

```text
SAMPLE	HGVS_INPUT	CHROM	POS	REF	ALT	STATUS
sample1	NM_000059.4:c.7007G>A;p.(Arg2336His)	chr13	...	G	A	OK
```

The complete output is written to the requested TSV file.

For safety, the workflow **does not overwrite an existing output file**.

## Status values

The `STATUS` column provides information about the conversion and should always be inspected before downstream validation.

The most important statuses are:

- `OK` — HGVS was successfully mapped without warnings.
- `REVIEW` — a genomic representation was obtained, but VariantValidator returned warnings or normalization/remapping differences that should be reviewed.
- `MISSING_REFERENCE` — no valid versioned reference accession was detected.
- `MULTIPLE_DESCRIPTIONS` — more than one DNA-level variant description was detected in the same input entry.
- `COMPLEX_HGVS_REVIEW` — complex HGVS representation requiring specific handling.
- `SV_REVIEW` — large genomic deletion requiring explicit structural-variant handling.
- `NO_TARGET_MAPPING` — no explicit mapping to the requested genome assembly was obtained.
- `TRANSCRIPT_VERSION_MISMATCH` — the requested transcript/version was not preserved in the returned mapping.
- `AMBIGUOUS_MAPPING` — multiple genomic representations were returned.
- `UNSUPPORTED_VCF_ALLELES` — the result could not be represented using supported explicit VCF alleles.
- `NON_PRIMARY_CONTIG` — the variant maps outside the expected primary chromosomes.
- `API_ERROR` — the external API request failed.
- `RESPONSE_ERROR` — an error occurred while interpreting the API response.

Variants marked as `REVIEW` should not automatically be considered incorrect; the status indicates that the returned representation requires inspection.

## Exome capture BED annotation

One or more capture kits can be supplied using `--bed`.

For example:

```bash
Rscript hgvs2vcf.R variants.tsv output.tsv GRCh38 \
  --bed Twist=hg38_Twist_Bioscience_for_Illumina_Exome_2_5.bed
```

Multiple BED files can be supplied:

```bash
Rscript hgvs2vcf.R variants.tsv output.tsv GRCh38 \
  --bed Twist=twist.bed \
  --bed Agilent=agilent.bed \
  --bed IDT=idt.bed.gz
```

Each BED file adds a new column to the output:

```text
BED_Twist
BED_Agilent
BED_IDT
```

BED files may be compressed with gzip.

All BED files must correspond to the same genome assembly selected for the analysis. **No liftover is performed.**

### BED annotation values

With the default `variant` mode, each variant can receive:

| Value | Meaning |
|---|---|
| `IN` | The affected sequence is completely inside the capture intervals |
| `OUT` | There is no overlap with the capture intervals |
| `PARTIAL` | Only part of the affected sequence overlaps the capture intervals |
| `BORDER` | For an insertion, only one side of the insertion breakpoint is captured |
| `NA` | The variant could not be evaluated |

For SNVs, `IN` or `OUT` will normally be sufficient.

For INDELs, the workflow evaluates the actual affected sequence rather than relying exclusively on the VCF anchor position.

## BED overlap modes

Two overlap strategies are available.

### `variant` — default

```bash
--bed-mode variant
```

Evaluates the genomic span affected by the variant.

For deletions and substitutions, unchanged VCF prefix/suffix anchor bases are excluded before overlap is calculated.

For insertions, the insertion breakpoint is evaluated using the two adjacent reference bases.

This is the recommended mode when determining whether the actual variant is expected to fall within a capture target.

### `pos`

```bash
--bed-mode pos
```

Evaluates only the VCF `POS` coordinate.

Possible results are:

```text
IN
OUT
NA
```

Because VCF positions for INDELs may correspond to an anchor base, `pos` and `variant` modes can produce different results.

## Annotating an existing converted table

An existing hgvs2vcf output can be annotated with additional BED files without querying VariantValidator again.

Use:

```bash
Rscript hgvs2vcf.R previous_output.tsv annotated_output.tsv GRCh38 \
  --table \
  --bed Twist=twist.bed
```

This is useful when adding or comparing additional capture kits after the original HGVS conversion.

## Using the functions directly in R

The conversion functions can also be used interactively.

```r
source("hgvs_para_vcf.R")

hgvs <- c(
  "NM_000059.4:c.7007G>A;p.(Arg2336His)",
  "NM_007294.4:c.5266dup;p.(Gln1756ProfsTer74)"
)

result <- hgvs_para_vcf(
  hgvs,
  assembly = "GRCh38"
)

result[, c(
  "HGVS_INPUT",
  "CHROM",
  "POS",
  "REF",
  "ALT",
  "STATUS"
)]
```

A preliminary HGVS screening can be performed without querying the API:

```r
screening <- hgvs_para_vcf(
  hgvs,
  assembly = "GRCh38",
  consultar = FALSE
)

table(screening$STATUS)
```

Capture BEDs can also be supplied directly:

```r
beds <- c(
  Twist = "twist.bed",
  Agilent = "agilent.bed.gz"
)

result <- hgvs_para_vcf(
  hgvs,
  assembly = "GRCh38",
  beds = beds
)
```

## Data handling

Only the cleaned HGVS variant description and selected genome assembly are required for the external VariantValidator request.

The `SAMPLE` identifier is maintained locally by the workflow and is not part of the HGVS query.

Nevertheless, users working with clinical or identifiable data should follow their institution's data protection and ethical requirements.

## Important considerations

This workflow is intended to support variant validation and coordinate harmonization. It does **not** replace manual review of clinically relevant variants.

In particular:

- transcript accessions should include their version;
- transcript versions are not automatically substituted;
- protein annotations supplied by the user are not independently validated against the DNA description;
- API warnings result in `REVIEW` rather than `OK`;
- structural variants require additional representation and validation;
- the returned VCF representation may be equivalent to, but not necessarily identical to, the representation produced by the original variant caller in repetitive regions;
- BED files must use the same genome assembly as the converted variants;
- chromosome-name harmonization (`chr1` versus `1`, for example) does not perform genome-assembly conversion;
- no liftover is performed.

For diagnostic or research validation, the final genomic coordinates should therefore be checked against the corresponding VCF and, when appropriate, the original sequencing data.

## External service

HGVS validation and genomic mapping are performed using the **VariantValidator REST API**.

Service availability and authentication requirements are controlled by the external provider and may change.

If an authorized token is required, it can be provided through:

```bash
export VV_TOKEN="your_token"
```

before running the workflow.

## Author

**Marco Antonio Campanário**

Human genomics and bioinformatics.