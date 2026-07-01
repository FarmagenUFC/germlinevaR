# Modular, individually-toggleable filtering of a read.gvr table

Applies a set of independent variant filters to the MAF-like table
produced by
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md).
Each distinct filter is its own argument; setting an argument to `NULL`
disables that filter entirely (no rows removed by it). With all
defaults, `gvr_filter(gvr)` reproduces the canonical rare /
clinically-relevant / called-genotype pipeline (AF filters + CLIN_SIG +
GT exclusion).

## Usage

``` r
gvr_filter(
  gvr,
  gnomADe_AF = 0.01,
  AF = 0.01,
  ABraOM_AF = 0.01,
  gnomADe_AF_keep_missing = TRUE,
  AF_keep_missing = TRUE,
  ABraOM_AF_keep_missing = TRUE,
  clin_sig_terms = c("likely_pathogenic", "pathogenic", "uncertain_significance"),
  clin_sig_keep_missing = TRUE,
  remove_benign = FALSE,
  biotype_keep = NULL,
  gt_exclude = c("0", "0/0"),
  vc_nonSyn = FALSE,
  missense_only = FALSE,
  genes = NULL,
  panel = NULL,
  hpo = NULL,
  save_excel = FALSE,
  out_dir = NULL,
  file_prefix = "gvr_filter",
  verbose = TRUE
)
```

## Arguments

- gvr:

  data.table / data.frame from
  [`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
  (or compatible). Filtered in a copy; the input object is not modified.

- gnomADe_AF:

  Numeric upper threshold for the gnomAD exome AF column `gnomADe_AF`
  (keep rows with AF \< threshold). `NULL` disables this filter. Default
  0.01.

- AF:

  Numeric upper threshold for the `AF` column (gnomAD genome / VCF
  allele frequency). `NULL` disables this filter. Default 0.01.

- ABraOM_AF:

  Numeric upper threshold for the Brazilian-cohort (ABraOM SABE 609)
  column `ABraOM_AF`. `NULL` disables this filter. Default 0.01.

- gnomADe_AF_keep_missing:

  Logical; if TRUE (default), keep rows whose `gnomADe_AF` is missing
  (NA or ""); if FALSE drop them. Ignored when `gnomADe_AF` is NULL.

- AF_keep_missing:

  Logical; missing-value handling for the `AF` filter. TRUE (default)
  keeps missing. Ignored when `AF` is NULL.

- ABraOM_AF_keep_missing:

  Logical; missing-value handling for the ABraOM filter. TRUE (default)
  retains variants absent from the Brazilian cohort (where absence often
  means "not catalogued", not "common"). Ignored when `ABraOM_AF` is
  NULL.

- clin_sig_terms:

  Character vector of clinical-significance terms to keep (substring,
  case-insensitive, OR-combined). `NULL` disables the CLIN_SIG filter.
  Default: c("likely_pathogenic","pathogenic","uncertain_significance").

- clin_sig_keep_missing:

  Logical; if TRUE (default) rows with missing CLIN_SIG (NA/"") are
  kept. Only relevant when `clin_sig_terms` is non-NULL.

- remove_benign:

  Logical; if TRUE, remove rows whose `CLIN_SIG` contains "benign"
  (substring, case-insensitive). This catches `benign`, `likely_benign`,
  and compound annotations like
  `"uncertain_significance&likely_benign"`. Applied AFTER the
  `clin_sig_terms` keep-filter, so a row that matched a wanted term but
  also contains "benign" is still removed. `FALSE` (default) does not
  remove benign rows.

- biotype_keep:

  Character vector of BIOTYPE values to keep (exact match via %in%).
  `NULL` (default) disables the biotype filter — all biotypes are kept.
  Pass e.g. `c("protein_coding", "protein_coding_LoF")` to restrict to
  protein-coding transcripts.

- gt_exclude:

  Character vector of GT values to remove (exact match). `NULL` disables
  the genotype filter. Default: c("0","0/0").

- vc_nonSyn:

  Logical or character vector. Controls which `Variant_Classification`
  values are retained. `FALSE` (default) keeps all. `TRUE` keeps only
  the 9 protein-altering classes (Frame_Shift_Del, Frame_Shift_Ins,
  Splice_Site, Translation_Start_Site, Nonsense_Mutation,
  Nonstop_Mutation, In_Frame_Del, In_Frame_Ins, Missense_Mutation). A
  custom character vector keeps only those classifications. Rows with
  missing/blank `Variant_Classification` are removed when this filter is
  active.

- missense_only:

  Logical; if `TRUE`, keep only rows whose `Variant_Classification`
  equals `"Missense_Mutation"` (added in vN+5). Default `FALSE`
  preserves prior behaviour byte-for-byte. Combines non-contradictorily
  with `vc_nonSyn`: `vc_nonSyn` runs first (keeping 9 protein-altering
  classes), then `missense_only` narrows to the missense subset. Errors
  with a clear message if `Variant_Classification` is missing.

- genes:

  Character vector of `Hugo_Symbol`s to keep (exact, case-insensitive),
  or `NULL` (default) to keep all genes.

- panel:

  Character vector of curated disease panel name(s) (e.g.
  `"breast cancer"`, `"hereditary prostate cancer"`, `"gist"`). Resolved
  to gene vectors via
  [`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md)
  and unioned with `genes` (dedup, uppercased) before the "Gene subset"
  filter runs. See
  [`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md)
  for the full registry. `NULL` (default) disables panel filtering;
  behaviour is then byte-identical to omitting the argument.

- hpo:

  Character vector of Human Phenotype Ontology (HPO) term id(s), e.g.
  `"HP:0003002"` (Breast carcinoma). Resolved to gene vectors via
  [`gvr_hpo_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_hpo_genes.md)
  and unioned with any genes from `genes` / `panel` before the "Gene
  subset" filter runs. Lenient input accepted: `"HP:0003002"`,
  `"hp:0003002"`, `"3002"`, and `"0003002"` all normalise to canonical
  `"HP:0003002"`. Only exact-term associations are used (no
  ontology-descendant expansion). The HPO table is downloaded and cached
  under `tools::R_user_dir("germlinevaR", "cache")`; see
  [`gvr_hpo_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_hpo_genes.md)
  for offline / air-gapped use via `hpo_path=`. `NULL` (default)
  disables HPO filtering; behaviour is then byte-identical to omitting
  the argument.

- save_excel:

  Logical; if TRUE, also write the FILTERED table to an `.xlsx` workbook
  (single `"Filtered"` sheet) at `<out_dir>/<file_prefix>.xlsx`.
  Requires the openxlsx package (a `Suggests` dependency); if it is not
  installed the export is skipped with a warning. Default FALSE. The
  write is a side effect only: the returned `data.table` is identical
  whether or not `save_excel` is TRUE.

- out_dir:

  Output directory for the Excel file. `NULL` (default) uses the current
  working directory. Created if it does not exist. Only used when
  `save_excel = TRUE`.

- file_prefix:

  Filename prefix (without extension) for the Excel file. Default
  `"gvr_filter"` -\> `gvr_filter.xlsx`. Only used when
  `save_excel = TRUE`.

- verbose:

  Logical; if TRUE (default) print a per-filter breakdown (rows in -\>
  out and rows removed by each active step).

## Value

A `data.table` of the surviving rows, with the same columns as the
input. A plain `data.frame` input is returned as a `data.table`. The
input object is not modified. With `verbose = TRUE`, a per-filter
breakdown (rows in -\> out, and rows removed by each active step) is
printed as it runs.

## Details

Filters are applied in a fixed order; each step operates on the
survivors of the previous one:

1.  gnomAD exome AF - `gnomADe_AF` (+ `gnomADe_AF_keep_missing`)

2.  gnomAD genome / VCF AF - `AF` (+ `AF_keep_missing`)

3.  ABraOM AF - `ABraOM_AF` (+ `ABraOM_AF_keep_missing`)

4.  Clinical significance - `clin_sig_terms` (+ `clin_sig_keep_missing`)

5.  Remove benign - `remove_benign`

6.  Biotype - `biotype_keep`

7.  Genotype exclusion - `gt_exclude`

8.  Variant classification - `vc_nonSyn`

9.  Gene subset - `genes` (unioned with `panel` and `hpo`)

Important data notes (true of
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
output):

- AF columns are CHARACTER (e.g. `"0.8781"`), so they are coerced with
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) before
  comparison.

- "Missing" means EITHER `NA` OR empty string `""`
  ([`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
  uses `""` for absent values). Both are treated as missing everywhere
  in this function.

- `CLIN_SIG` matching is SUBSTRING + case-insensitive. A compound
  annotation such as `"pathogenic&benign"` or
  `"uncertain_significance&likely_benign&benign"` is KEPT because it
  CONTAINS a wanted term. This matches the dplyr `str_detect()` /
  [`grepl()`](https://rdrr.io/r/base/grep.html) convention. Use
  exact-token matching only if you split `CLIN_SIG` yourself.

- The default `gt_exclude = c("0", "0/0")` is a no-op on data whose `GT`
  column only contains called alt genotypes (e.g. `0/1`, `1/1`, `1/2`);
  it is retained for portability to data that does carry `"0"`/`"0/0"`.

The input is never modified: filtering operates on an internal
`data.table` copy.

## See also

[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
to build the table,
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
to summarise the filtered variants.

Other germlinevaR:
[`gvr_hpo_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_hpo_genes.md),
[`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md),
[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md),
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md),
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md),
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)

## Author

germlinevaR authors

## Examples

``` r
## Load the shipped example table
gvr <- readRDS(system.file("extdata", "example_gvr.rds",
    package = "germlinevaR"))
## Default filter (rare + clinically relevant + called genotypes).
## The example table was built without ABraOM annotation, so we
## disable the ABraOM_AF filter to avoid a "column not found" warning.
filt <- gvr_filter(gvr, ABraOM_AF = NULL, verbose = FALSE)
dim(filt)
#> [1]   7 115

## Variations on the default pipeline using the shipped example table

## Add protein-coding biotype filter:
gvr_filter(gvr, ABraOM_AF = NULL,
    biotype_keep = c("protein_coding", "protein_coding_LoF"),
    verbose = FALSE)
#>    Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position
#>         <char>         <char> <char>     <char>     <char>          <num>
#> 1:      ATAD3B              0      .     GRCh38       chr1        1479267
#> 2:        ESPN              0      .     GRCh38       chr1        6445777
#> 3:    PRAMEF33              0      .     GRCh38       chr1       13306220
#> 4:    C11orf21              0      .     GRCh38      chr11        2301755
#> 5:       MUC19              0      .     GRCh38      chr12       40485636
#> 6:       BRCA2              0      .     GRCh38      chr13       32338482
#> 7:      ZNF747              0      .     GRCh38      chr16       30534677
#>    End_Position Strand Variant_Classification Variant_Type Reference_Allele
#>           <num> <char>                 <char>       <char>           <char>
#> 1:      1479267      +                 Intron          SNP                C
#> 2:      6445781      +        Frame_Shift_Del          DEL            AGCTT
#> 3:     13306221      +        Frame_Shift_Ins          INS                -
#> 4:      2301755      +            Splice_Site          SNP                C
#> 5:     40485636      +       Nonstop_Mutation          SNP                A
#> 6:     32338482      +      Missense_Mutation          SNP                G
#> 7:     30534677      + Translation_Start_Site          SNP                C
#>    Tumor_Seq_Allele1 Tumor_Seq_Allele2         Genotype     dbSNP_RS
#>               <char>            <char>           <char>       <char>
#> 1:                 C                 T              C/T             
#> 2:             AGCTT                 -          AGCTT/-  rs753994746
#> 3:                 -    GGCCCAGAAGGTTC -/GGCCCAGAAGGTTC rs1553122356
#> 4:                 C                 G              C/G             
#> 5:                 A                 T              A/T  rs192078109
#> 6:                 G                 C              G/C             
#> 7:                 C                 T              C/T  rs933944167
#>    Tumor_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2
#>                  <char>                 <char>                 <char>
#> 1:            Sample_01                                              
#> 2:            Sample_01                                              
#> 3:            Sample_01                                              
#> 4:            Sample_01                                              
#> 5:            Sample_01                                              
#> 6:            Sample_01                                              
#> 7:            Sample_01                                              
#>             HGVSc                 HGVSp   HGVSp_Short   Transcript_ID
#>            <char>                <char>        <char>          <char>
#> 1:   c.444+159C>T                                     ENST00000673477
#> 2: c.1306_1310del  p.(Ser436ProfsTer31)  p.S436Pfs*31 ENST00000645284
#> 3:   c.268_281dup   p.(Trp98PhefsTer37)   p.W98Ffs*37 ENST00000437300
#> 4:      c.53+1G>C                                     ENST00000381153
#> 5:     c.12684A>T p.(Ter4228CysextTer?) p.*4228Cext*? ENST00000454784
#> 6:      c.4127G>C        p.(Gly1376Ala)      p.G1376A ENST00000380152
#> 7:         c.3G>A             p.(Met1?)         p.M1? ENST00000693075
#>             Consequence t_depth t_ref_count t_alt_count         Allele   IMPACT
#>                  <char>  <char>      <char>      <char>         <char>   <char>
#> 1:       intron_variant      13           9           4              T MODIFIER
#> 2:   frameshift_variant      16          13           3              -     HIGH
#> 3:   frameshift_variant      16          12           4 GGCCCAGAAGGTTC     HIGH
#> 4: splice_donor_variant     101          81          20              G     HIGH
#> 5:            stop_lost      37          30           7              T     HIGH
#> 6:     missense_variant     165         136          29              C MODERATE
#> 7:           start_lost      70          34          36              T     HIGH
#>      SYMBOL            Gene Feature_type         Feature            BIOTYPE
#>      <char>          <char>       <char>          <char>             <char>
#> 1:   ATAD3B ENSG00000160072   Transcript ENST00000673477     protein_coding
#> 2:     ESPN ENSG00000187017   Transcript ENST00000645284     protein_coding
#> 3: PRAMEF33 ENSG00000237700   Transcript ENST00000437300     protein_coding
#> 4: C11orf21 ENSG00000110665   Transcript ENST00000381153     protein_coding
#> 5:    MUC19 ENSG00000205592   Transcript ENST00000454784 protein_coding_LoF
#> 6:    BRCA2 ENSG00000139618   Transcript ENST00000380152     protein_coding
#> 7:   ZNF747 ENSG00000169955   Transcript ENST00000693075     protein_coding
#>      EXON INTRON cDNA_position CDS_position Protein_position Amino_acids
#>    <char> <char>        <char>       <char>           <char>      <char>
#> 1:          4/15                                                        
#> 2:   7/13            1486-1490    1306-1310          436-437        SF/X
#> 3:    2/4              339-340      266-267               89    V/VAQKVX
#> 4:           1/3                                                        
#> 5: 56/173                12684        12684             4228         */C
#> 6:  11/27                 4326         4127             1376         G/A
#> 7:    1/3                  182            3                1         M/I
#>                   Codons     Existing_variation DISTANCE STRAND  FLAGS
#>                   <char>                 <char>   <char> <char> <char>
#> 1:                                                            1       
#> 2:              AGCTTc/c            rs753994746               1       
#> 3: gtg/gtGGCCCAGAAGGTTCg           rs1553122356               1       
#> 4:                                                           -1       
#> 5:               tgA/tgT            rs192078109               1       
#> 6:               gGa/gCa CD1413320&COSV66459147               1       
#> 7:               atG/atA            rs933944167              -1       
#>    VARIANT_CLASS SYMBOL_SOURCE    HGNC_ID CANONICAL        MANE    MANE_SELECT
#>           <char>        <char>     <char>    <char>      <char>         <char>
#> 1:           SNV          HGNC HGNC:24007       YES MANE_Select    NM_031921.6
#> 2:      deletion          HGNC HGNC:13281       YES MANE_Select    NM_031475.3
#> 3:     insertion          HGNC HGNC:49193       YES MANE_Select NM_001291381.1
#> 4:           SNV          HGNC HGNC:13231       YES MANE_Select NM_001329958.2
#> 5:           SNV          HGNC HGNC:14362       YES                           
#> 6:           SNV          HGNC  HGNC:1101       YES MANE_Select    NM_000059.4
#> 7:           SNV          HGNC HGNC:28350       YES MANE_Select NM_001305018.2
#>    MANE_PLUS_CLINICAL    TSL APPRIS        CCDS            ENSP     SWISSPROT
#>                <char> <char> <char>      <char>          <char>        <char>
#> 1:                               P1    CCDS30.1 ENSP00000500094    Q5T9A4.167
#> 2:                               P1    CCDS70.1 ENSP00000496593    B1AK53.128
#> 3:                         1     P1 CCDS85928.1 ENSP00000492439 A0A0G2JMD5.51
#> 4:                         1     A2 CCDS86168.1 ENSP00000370545    Q9P2W6.113
#> 5:                         5     P1             ENSP00000508949              
#> 6:                         5     A2  CCDS9344.1 ENSP00000369497    P51587.242
#> 7:                               A2 CCDS92140.1 ENSP00000509633              
#>           TREMBL       UNIPARC UNIPROT_ISOFORM GENE_PHENO
#>           <char>        <char>          <char>     <char>
#> 1:               UPI000013E044        Q5T9A4-1           
#> 2:               UPI000013D2B6        B1AK53-1          1
#> 3:               UPI000442CEFE                           
#> 4:               UPI0000127A63                           
#> 5:                                                       
#> 6:               UPI00001FCBCC                          1
#> 7: A0A8I5KWK6.11 UPI000004CC0E                           
#>                              SIFT     PolyPhen
#>                            <char>       <char>
#> 1:                                            
#> 2:                                            
#> 3:                                            
#> 4:                                            
#> 5:                                            
#> 6:                tolerated(0.18) benign(0.29)
#> 7: tolerated_low_confidence(0.31)   unknown(0)
#>                                                                                                                                      DOMAINS
#>                                                                                                                                       <char>
#> 1:                                                                                                                                          
#> 2: PANTHER:PTHR24153&Low_complexity_(Seg):seg&Prints:PR01217&MobiDB_lite:mobidb-lite&MobiDB_lite:mobidb-lite&AFDB-ENSP_mappings:AF-B1AK53-F1
#> 3:                                                                   PANTHER:PTHR14224&PIRSF:PIRSF038286&AFDB-ENSP_mappings:AF-A0A0G2JMD5-F1
#> 4:                                                                                                                                          
#> 5:                                                                                                                                          
#> 6:                                                                                                       PANTHER:PTHR11289&PIRSF:PIRSF002397
#> 7:                                                                                                                   MobiDB_lite:mobidb-lite
#>     miRNA HGVS_OFFSET     AF AFR_AF AMR_AF EAS_AF EUR_AF SAS_AF gnomADe_AF
#>    <char>      <char> <char> <char> <char> <char> <char> <char>     <char>
#> 1:                                                                        
#> 2:                                                                0.002621
#> 3:                 15                                            0.0003576
#> 4:                                                                        
#> 5:                                                                0.004832
#> 6:                                                                        
#> 7:                                                               1.443e-06
#>    gnomADe_AFR_AF gnomADe_AMR_AF gnomADe_ASJ_AF gnomADe_EAS_AF gnomADe_FIN_AF
#>            <char>         <char>         <char>         <char>         <char>
#> 1:                                                                           
#> 2:       0.002249       0.001438        0.00167       0.003474       0.001335
#> 3:        0.01328      0.0005464              0              0              0
#> 4:                                                                           
#> 5:       0.003774        0.01013       0.007724        0.00188              0
#> 6:                                                                           
#> 7:              0      2.759e-05              0              0              0
#>    gnomADe_MID_AF gnomADe_NFE_AF gnomADe_REMAINING_AF gnomADe_SAS_AF gnomADg_AF
#>            <char>         <char>               <char>         <char>     <char>
#> 1:                                                                             
#> 2:       0.002595       0.002841             0.003058      0.0008685    0.01296
#> 3:       0.000831      8.943e-06            0.0008676      1.429e-05   0.003988
#> 4:                                                                             
#> 5:       0.005085       0.004811             0.005435       0.005452    0.03926
#> 6:                                                                             
#> 7:              0      9.265e-07                    0              0  6.571e-06
#>    gnomADg_AFR_AF gnomADg_AMI_AF gnomADg_AMR_AF gnomADg_ASJ_AF gnomADg_EAS_AF
#>            <char>         <char>         <char>         <char>         <char>
#> 1:                                                                           
#> 2:        0.01518        0.01202        0.01318        0.01669       0.009633
#> 3:        0.01408              0       0.001057              0              0
#> 4:                                                                           
#> 5:        0.03665        0.03018        0.04658        0.06261        0.02314
#> 6:                                                                           
#> 7:              0              0              0              0              0
#>    gnomADg_FIN_AF gnomADg_MID_AF gnomADg_NFE_AF gnomADg_REMAINING_AF
#>            <char>         <char>         <char>               <char>
#> 1:                                                                  
#> 2:        0.01315       0.004202        0.01189              0.01238
#> 3:              0              0      4.445e-05             0.002388
#> 4:                                                                  
#> 5:        0.02334        0.04444        0.04173               0.0444
#> 6:                                                                  
#> 7:              0              0       1.47e-05                    0
#>    gnomADg_SAS_AF    MAX_AF MAX_AF_POPS CLIN_SIG SOMATIC  PHENO   PUBMED
#>            <char>    <char>      <char>   <char>  <char> <char>   <char>
#> 1:                                                                      
#> 2:       0.009645   0.01669 gnomADg_ASJ                         33968136
#> 3:              0   0.01408 gnomADg_AFR                                 
#> 4:                                                                      
#> 5:        0.04893   0.06261 gnomADg_ASJ                                 
#> 6:                                                   1&1    1&1         
#> 7:              0 2.759e-05 gnomADe_AMR                                 
#>    MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE TRANSCRIPTION_FACTORS
#>        <char>    <char>       <char>             <char>                <char>
#> 1:                                                                       <NA>
#> 2:                                                                       <NA>
#> 3:                                                                       <NA>
#> 4:                                                                       <NA>
#> 5:                                                                       <NA>
#> 6:                                                                       <NA>
#> 7:                                                                       <NA>
#>    FILTER   QUAL INFO_DP INFO_AC INFO_AF INFO_MQ INFO_QD  CNN_1D     GT     AD
#>    <char> <char>  <char>  <char>  <char>  <char>  <char>  <char> <char> <char>
#> 1:   PASS  98.64      13       1   0.500   45.94    7.59  -0.887    0/1    9,4
#> 2:   PASS  79.60      17       1   0.500   55.44    4.98  -8.136    0/1   13,3
#> 3:   PASS 123.60      22       1   0.500   27.42    7.73  -5.337    0/1   12,4
#> 4:   PASS 151.64     107       1   0.500   60.00    1.50  -8.904    0/1  81,20
#> 5:   PASS 196.64      37       1   0.500   57.92    5.31  -5.636    0/1   30,7
#> 6:   PASS 142.64     169       1   0.500   60.00    0.86 -15.246    0/1 136,29
#> 7:   PASS 741.64      70       1   0.500   60.00   10.59   0.399    0/1  34,36
#>    sample_DP     GQ
#>       <char> <char>
#> 1:        13     99
#> 2:        16     87
#> 3:        16     99
#> 4:       101     99
#> 5:        37     99
#> 6:       165     99
#> 7:        70     99

## Only the rarity filter on gnomAD exome AF, nothing else:
gvr_filter(gvr, gnomADe_AF = 0.001, AF = NULL, ABraOM_AF = NULL,
    clin_sig_terms = NULL, gt_exclude = NULL,
    vc_nonSyn = FALSE, genes = NULL, verbose = FALSE)
#>     Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position
#>          <char>         <char> <char>     <char>     <char>          <num>
#>  1:      KLHL17              0      .     GRCh38       chr1         965299
#>  2:     PLEKHN1              0      .     GRCh38       chr1         966179
#>  3:      ATAD3B              0      .     GRCh38       chr1        1479267
#>  4:    PRAMEF33              0      .     GRCh38       chr1       13306220
#>  5:        CD1C              0      .     GRCh38       chr1      158288048
#>  6:    C11orf21              0      .     GRCh38      chr11        2301755
#>  7:       BRCA2              0      .     GRCh38      chr13       32338482
#>  8:     ARHGEF7              0      .     GRCh38      chr13      111245537
#>  9:     Unknown              0      .     GRCh38      chr14       35405648
#> 10:      ZNF747              0      .     GRCh38      chr16       30534677
#> 11:        RFFL              0      .     GRCh38      chr17       35015003
#> 12:         MFF              0      .     GRCh38       chr2      227343060
#> 13:     FAM246C              0      .     GRCh38      chr22       19029801
#>     End_Position Strand Variant_Classification Variant_Type Reference_Allele
#>            <num> <char>                 <char>       <char>           <char>
#>  1:       965299      +                  3'UTR          SNP                A
#>  2:       966179      +                5'Flank          SNP                G
#>  3:      1479267      +                 Intron          SNP                C
#>  4:     13306221      +        Frame_Shift_Ins          INS                -
#>  5:    158288048      +                5'Flank          SNP                A
#>  6:      2301755      +            Splice_Site          SNP                C
#>  7:     32338482      +      Missense_Mutation          SNP                G
#>  8:    111245537      +                 Intron          SNP                T
#>  9:     35405648      +                    RNA          SNP                T
#> 10:     30534677      + Translation_Start_Site          SNP                C
#> 11:     35015003      +                 Intron          SNP                C
#> 12:    227343060      +                 Intron          SNP                C
#> 13:     19029801      +        Frame_Shift_Del          DEL                G
#>     Tumor_Seq_Allele1 Tumor_Seq_Allele2         Genotype     dbSNP_RS
#>                <char>            <char>           <char>       <char>
#>  1:                 A                 G              A/G   rs61531461
#>  2:                 G                 A              G/A   rs13303160
#>  3:                 C                 T              C/T             
#>  4:                 -    GGCCCAGAAGGTTC -/GGCCCAGAAGGTTC rs1553122356
#>  5:                 A                 C              A/C   rs10797006
#>  6:                 C                 G              C/G             
#>  7:                 G                 C              G/C             
#>  8:                 C                 C              C/C    rs9522171
#>  9:                 T                 C              T/C    rs3138053
#> 10:                 C                 T              C/T  rs933944167
#> 11:                 C                 G              C/G    rs2269857
#> 12:                 A                 A              A/A    rs6740870
#> 13:                 -                 -              -/-   rs11356224
#>     Tumor_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2
#>                   <char>                 <char>                 <char>
#>  1:            Sample_01                                              
#>  2:            Sample_01                                              
#>  3:            Sample_01                                              
#>  4:            Sample_01                                              
#>  5:            Sample_01                                              
#>  6:            Sample_01                                              
#>  7:            Sample_01                                              
#>  8:            Sample_01                                              
#>  9:            Sample_01                                              
#> 10:            Sample_01                                              
#> 11:            Sample_01                                              
#> 12:            Sample_01                                              
#> 13:            Sample_01                                              
#>             HGVSc               HGVSp HGVSp_Short   Transcript_ID
#>            <char>              <char>      <char>          <char>
#>  1:     c.*108A>G                                 ENST00000338591
#>  2:                                               ENST00000379410
#>  3:  c.444+159C>T                                 ENST00000673477
#>  4:  c.268_281dup p.(Trp98PhefsTer37) p.W98Ffs*37 ENST00000437300
#>  5:                                               ENST00000368170
#>  6:     c.53+1G>C                                 ENST00000381153
#>  7:     c.4127G>C      p.(Gly1376Ala)    p.G1376A ENST00000380152
#>  8: c.950+1243T>C                                 ENST00000646102
#>  9:      n.953T>C                                 ENST00000848851
#> 10:        c.3G>A           p.(Met1?)       p.M1? ENST00000693075
#> 11:  c.887-240G>C                                 ENST00000394597
#> 12: c.440+2680C>A                                 ENST00000304593
#> 13:      c.280del p.(Ala94LeufsTer43) p.A94Lfs*43 ENST00000652053
#>                            Consequence t_depth t_ref_count t_alt_count
#>                                 <char>  <char>      <char>      <char>
#>  1:                3_prime_UTR_variant      13           5           8
#>  2:              upstream_gene_variant      39          20          19
#>  3:                     intron_variant      13           9           4
#>  4:                 frameshift_variant      16          12           4
#>  5:              upstream_gene_variant      90          47          43
#>  6:               splice_donor_variant     101          81          20
#>  7:                   missense_variant     165         136          29
#>  8:                     intron_variant      38           0          38
#>  9: non_coding_transcript_exon_variant     100          50          50
#> 10:                         start_lost      70          34          36
#> 11:                     intron_variant      13           8           5
#> 12:                     intron_variant      12           0          12
#> 13:                 frameshift_variant      28           0          28
#>             Allele   IMPACT   SYMBOL            Gene Feature_type
#>             <char>   <char>   <char>          <char>       <char>
#>  1:              G MODIFIER   KLHL17 ENSG00000187961   Transcript
#>  2:              A MODIFIER  PLEKHN1 ENSG00000187583   Transcript
#>  3:              T MODIFIER   ATAD3B ENSG00000160072   Transcript
#>  4: GGCCCAGAAGGTTC     HIGH PRAMEF33 ENSG00000237700   Transcript
#>  5:              C MODIFIER     CD1C ENSG00000158481   Transcript
#>  6:              G     HIGH C11orf21 ENSG00000110665   Transcript
#>  7:              C MODERATE    BRCA2 ENSG00000139618   Transcript
#>  8:              C MODIFIER  ARHGEF7 ENSG00000102606   Transcript
#>  9:              C MODIFIER          ENSG00000310289   Transcript
#> 10:              T     HIGH   ZNF747 ENSG00000169955   Transcript
#> 11:              G MODIFIER     RFFL ENSG00000092871   Transcript
#> 12:              A MODIFIER      MFF ENSG00000168958   Transcript
#> 13:              -     HIGH  FAM246C ENSG00000286025   Transcript
#>             Feature            BIOTYPE   EXON INTRON cDNA_position CDS_position
#>              <char>             <char> <char> <char>        <char>       <char>
#>  1: ENST00000338591     protein_coding  12/12                 2147             
#>  2: ENST00000379410     protein_coding                                         
#>  3: ENST00000673477     protein_coding          4/15                           
#>  4: ENST00000437300     protein_coding    2/4              339-340      266-267
#>  5: ENST00000368170     protein_coding                                         
#>  6: ENST00000381153     protein_coding           1/3                           
#>  7: ENST00000380152     protein_coding  11/27                 4326         4127
#>  8: ENST00000646102     protein_coding          8/21                           
#>  9: ENST00000848851             lncRNA    1/1                  953             
#> 10: ENST00000693075     protein_coding    1/3                  182            3
#> 11: ENST00000394597     protein_coding           5/6                           
#> 12: ENST00000304593     protein_coding           5/8                           
#> 13: ENST00000652053 protein_coding_LoF    1/1                  278          278
#>     Protein_position Amino_acids                Codons     Existing_variation
#>               <char>      <char>                <char>                 <char>
#>  1:                                                                rs61531461
#>  2:                                                                rs13303160
#>  3:                                                                          
#>  4:               89    V/VAQKVX gtg/gtGGCCCAGAAGGTTCg           rs1553122356
#>  5:                                                                rs10797006
#>  6:                                                                          
#>  7:             1376         G/A               gGa/gCa CD1413320&COSV66459147
#>  8:                                                                 rs9522171
#>  9:                                                        rs3138053&CR078726
#> 10:                1         M/I               atG/atA            rs933944167
#> 11:                                                                 rs2269857
#> 12:                                                                 rs6740870
#> 13:               93         R/X                cGg/cg             rs11356224
#>     DISTANCE STRAND  FLAGS VARIANT_CLASS SYMBOL_SOURCE    HGNC_ID CANONICAL
#>       <char> <char> <char>        <char>        <char>     <char>    <char>
#>  1:               1                  SNV          HGNC HGNC:24023       YES
#>  2:      303      1                  SNV          HGNC HGNC:25284       YES
#>  3:               1                  SNV          HGNC HGNC:24007       YES
#>  4:               1            insertion          HGNC HGNC:49193       YES
#>  5:     1875      1                  SNV          HGNC  HGNC:1636       YES
#>  6:              -1                  SNV          HGNC HGNC:13231       YES
#>  7:               1                  SNV          HGNC  HGNC:1101       YES
#>  8:               1                  SNV          HGNC HGNC:15607       YES
#>  9:               1                  SNV                                YES
#> 10:              -1                  SNV          HGNC HGNC:28350       YES
#> 11:              -1                  SNV          HGNC HGNC:24821       YES
#> 12:               1                  SNV          HGNC HGNC:24858       YES
#> 13:               1             deletion          HGNC HGNC:54842       YES
#>            MANE    MANE_SELECT MANE_PLUS_CLINICAL    TSL APPRIS        CCDS
#>          <char>         <char>             <char> <char> <char>      <char>
#>  1: MANE_Select    NM_198317.3                         1     P1 CCDS30550.1
#>  2: MANE_Select    NM_032129.3                         1     P2     CCDS4.1
#>  3: MANE_Select    NM_031921.6                               P1    CCDS30.1
#>  4: MANE_Select NM_001291381.1                         1     P1 CCDS85928.1
#>  5: MANE_Select    NM_001765.3                         1     P1  CCDS1175.1
#>  6: MANE_Select NM_001329958.2                         1     A2 CCDS86168.1
#>  7: MANE_Select    NM_000059.4                         5     A2  CCDS9344.1
#>  8: MANE_Select NM_001354046.2                                  CCDS86360.1
#>  9:                                                                        
#> 10: MANE_Select NM_001305018.2                               A2 CCDS92140.1
#> 11: MANE_Select NM_001017368.2                         1     P4 CCDS11286.1
#> 12: MANE_Select NM_001277062.2                         2     P1 CCDS63140.1
#> 13:                                                          P1            
#>                ENSP     SWISSPROT        TREMBL       UNIPARC UNIPROT_ISOFORM
#>              <char>        <char>        <char>        <char>          <char>
#>  1: ENSP00000343930    Q6TDP4.154               UPI00001DFBF0                
#>  2: ENSP00000368720    Q494U1.138               UPI00001416D8        Q494U1-1
#>  3: ENSP00000500094    Q5T9A4.167               UPI000013E044        Q5T9A4-1
#>  4: ENSP00000492439 A0A0G2JMD5.51               UPI000442CEFE                
#>  5: ENSP00000357152    P29017.197               UPI000013DF78                
#>  6: ENSP00000370545    Q9P2W6.113               UPI0000127A63                
#>  7: ENSP00000369497    P51587.242               UPI00001FCBCC                
#>  8: ENSP00000495631               A0A2R8YG42.32 UPI000387C911                
#>  9:                                                                          
#> 10: ENSP00000509633               A0A8I5KWK6.11 UPI000004CC0E                
#> 11: ENSP00000378096    Q8WZ73.184               UPI000006D6B9        Q8WZ73-1
#> 12: ENSP00000304898    Q9GZY8.150               UPI00000720C1        Q9GZY8-2
#> 13: ENSP00000498832                                                          
#>     GENE_PHENO                           SIFT     PolyPhen
#>         <char>                         <char>       <char>
#>  1:                                                       
#>  2:                                                       
#>  3:                                                       
#>  4:                                                       
#>  5:                                                       
#>  6:                                                       
#>  7:          1                tolerated(0.18) benign(0.29)
#>  8:                                                       
#>  9:                                                       
#> 10:            tolerated_low_confidence(0.31)   unknown(0)
#> 11:                                                       
#> 12:          1                                            
#> 13:                                                       
#>                                                                     DOMAINS
#>                                                                      <char>
#>  1:                                                                        
#>  2:                                                                        
#>  3:                                                                        
#>  4: PANTHER:PTHR14224&PIRSF:PIRSF038286&AFDB-ENSP_mappings:AF-A0A0G2JMD5-F1
#>  5:                                                                        
#>  6:                                                                        
#>  7:                                     PANTHER:PTHR11289&PIRSF:PIRSF002397
#>  8:                                                                        
#>  9:                                                                        
#> 10:                                                 MobiDB_lite:mobidb-lite
#> 11:                                                                        
#> 12:                                                                        
#> 13:                                                                        
#>      miRNA HGVS_OFFSET     AF AFR_AF AMR_AF EAS_AF EUR_AF SAS_AF gnomADe_AF
#>     <char>      <char> <char> <char> <char> <char> <char> <char>     <char>
#>  1:                    0.0100 0.0363 0.0029      0      0      0  0.0008317
#>  2:                    0.6516 0.2194 0.7032 0.7123 0.9056 0.8753           
#>  3:                                                                        
#>  4:                 15                                            0.0003576
#>  5:                    0.4189 0.6778 0.4251 0.3552 0.2087 0.3466           
#>  6:                                                                        
#>  7:                                                                        
#>  8:                    0.5393 0.3699 0.4409 0.7083 0.5815 0.6207           
#>  9:                    0.2368 0.2307 0.2493   0.13 0.2883 0.2935           
#> 10:                                                               1.443e-06
#> 11:                    0.2778 0.6831 0.0994 0.1508 0.0915 0.1789           
#> 12:                    0.9557 0.8729 0.9899  0.999  0.992 0.9611           
#> 13:                  2 1.0000      1      1      1      1      1           
#>     gnomADe_AFR_AF gnomADe_AMR_AF gnomADe_ASJ_AF gnomADe_EAS_AF gnomADe_FIN_AF
#>             <char>         <char>         <char>         <char>         <char>
#>  1:         0.0256       0.002306              0              0              0
#>  2:                                                                           
#>  3:                                                                           
#>  4:        0.01328      0.0005464              0              0              0
#>  5:                                                                           
#>  6:                                                                           
#>  7:                                                                           
#>  8:                                                                           
#>  9:                                                                           
#> 10:              0      2.759e-05              0              0              0
#> 11:                                                                           
#> 12:                                                                           
#> 13:                                                                           
#>     gnomADe_MID_AF gnomADe_NFE_AF gnomADe_REMAINING_AF gnomADe_SAS_AF
#>             <char>         <char>               <char>         <char>
#>  1:      0.0007645      4.317e-05             0.002097      7.476e-05
#>  2:                                                                  
#>  3:                                                                  
#>  4:       0.000831      8.943e-06            0.0008676      1.429e-05
#>  5:                                                                  
#>  6:                                                                  
#>  7:                                                                  
#>  8:                                                                  
#>  9:                                                                  
#> 10:              0      9.265e-07                    0              0
#> 11:                                                                  
#> 12:                                                                  
#> 13:                                                                  
#>     gnomADg_AF gnomADg_AFR_AF gnomADg_AMI_AF gnomADg_AMR_AF gnomADg_ASJ_AF
#>         <char>         <char>         <char>         <char>         <char>
#>  1:    0.00804         0.0278              0       0.002851              0
#>  2:     0.7174         0.3169          0.867         0.7396         0.8764
#>  3:                                                                       
#>  4:   0.003988        0.01408              0       0.001057              0
#>  5:     0.3701         0.6226         0.3114         0.3835         0.2638
#>  6:                                                                       
#>  7:                                                                       
#>  8:     0.5312         0.4173         0.3695         0.4954         0.6205
#>  9:     0.2619         0.2479         0.2258         0.2313         0.2641
#> 10:  6.571e-06              0              0              0              0
#> 11:     0.2553         0.6091        0.07237         0.1393         0.1928
#> 12:     0.9586         0.8769         0.9967         0.9797         0.9718
#> 13:          1              1              1              1              1
#>     gnomADg_EAS_AF gnomADg_FIN_AF gnomADg_MID_AF gnomADg_NFE_AF
#>             <char>         <char>         <char>         <char>
#>  1:              0              0              0      0.0001185
#>  2:         0.7359         0.9126         0.8129         0.9026
#>  3:                                                            
#>  4:              0              0              0      4.445e-05
#>  5:         0.3664         0.2951         0.2551         0.2347
#>  6:                                                            
#>  7:                                                            
#>  8:         0.7248         0.5964          0.517         0.5748
#>  9:         0.1256         0.2503         0.2789         0.2895
#> 10:              0              0              0       1.47e-05
#> 11:         0.1575         0.1362         0.1531         0.1042
#> 12:          0.999         0.9987         0.9592         0.9929
#> 13:              1         0.9999              1              1
#>     gnomADg_REMAINING_AF gnomADg_SAS_AF    MAX_AF
#>                   <char>         <char>    <char>
#>  1:             0.004327      0.0002127    0.0363
#>  2:               0.7547         0.8649    0.9126
#>  3:                                              
#>  4:             0.002388              0   0.01408
#>  5:               0.3452         0.3383    0.6778
#>  6:                                              
#>  7:                                              
#>  8:                0.548         0.6204    0.7248
#>  9:               0.2457         0.2723    0.2935
#> 10:                    0              0 2.759e-05
#> 11:               0.2171          0.182    0.6831
#> 12:               0.9677           0.96     0.999
#> 13:                    1              1         1
#>                                                                                                                               MAX_AF_POPS
#>                                                                                                                                    <char>
#>  1:                                                                                                                                   AFR
#>  2:                                                                                                                           gnomADg_FIN
#>  3:                                                                                                                                      
#>  4:                                                                                                                           gnomADg_AFR
#>  5:                                                                                                                                   AFR
#>  6:                                                                                                                                      
#>  7:                                                                                                                                      
#>  8:                                                                                                                           gnomADg_EAS
#>  9:                                                                                                                                   SAS
#> 10:                                                                                                                           gnomADe_AMR
#> 11:                                                                                                                                   AFR
#> 12:                                                                                                                       EAS&gnomADg_EAS
#> 13: AFR&AMR&EAS&EUR&SAS&gnomADg_AFR&gnomADg_AMI&gnomADg_AMR&gnomADg_ASJ&gnomADg_EAS&gnomADg_MID&gnomADg_NFE&gnomADg_REMAINING&gnomADg_SAS
#>     CLIN_SIG SOMATIC  PHENO
#>       <char>  <char> <char>
#>  1:                        
#>  2:                        
#>  3:                        
#>  4:                        
#>  5:                        
#>  6:                        
#>  7:              1&1    1&1
#>  8:                       1
#>  9:   benign     0&1    1&1
#> 10:                        
#> 11:                        
#> 12:   benign              1
#> 13:                        
#>                                                                                                                                                                                                                                                                                                                PUBMED
#>                                                                                                                                                                                                                                                                                                                <char>
#>  1:                                                                                                                                                                                                                                                                                                                  
#>  2:                                                                                                                                                                                                                                                                                                          29422604
#>  3:                                                                                                                                                                                                                                                                                                                  
#>  4:                                                                                                                                                                                                                                                                                                                  
#>  5:                                                                                                                                                                                                                                                                                                                  
#>  6:                                                                                                                                                                                                                                                                                                                  
#>  7:                                                                                                                                                                                                                                                                                                                  
#>  8:                                                                                                                                                                                                                                                                                                          35835914
#>  9: 22295056&25541970&24176007&22742663&26870821&21424379&31475028&23867959&26371589&26834482&21274447&25223483&28389768&19500386&19797428&19798070&26488500&31226330&17463416&23487427&23996241&24368589&30943907&26252270&31519222&29754263&26060603&24066085&36353205&26885097&24578542&16540234&19036641&38143267
#> 10:                                                                                                                                                                                                                                                                                                                  
#> 11:                                                                                                                                                                                                                                                                                                                  
#> 12:                                                                                                                                                                                                                                                                                                                  
#> 13:                                                                                                                                                                                                                                                                                                                  
#>     MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE TRANSCRIPTION_FACTORS
#>         <char>    <char>       <char>             <char>                <char>
#>  1:                                                                       <NA>
#>  2:                                                                       <NA>
#>  3:                                                                       <NA>
#>  4:                                                                       <NA>
#>  5:                                                                       <NA>
#>  6:                                                                       <NA>
#>  7:                                                                       <NA>
#>  8:                                                                       <NA>
#>  9:                                                                       <NA>
#> 10:                                                                       <NA>
#> 11:                                                                       <NA>
#> 12:                                                                       <NA>
#> 13:                                                                       <NA>
#>     FILTER    QUAL INFO_DP INFO_AC INFO_AF INFO_MQ INFO_QD  CNN_1D     GT
#>     <char>  <char>  <char>  <char>  <char>  <char>  <char>  <char> <char>
#>  1:   PASS  283.64      14       1   0.500   60.00   21.82   1.847    0/1
#>  2:   PASS  632.64      39       1   0.500   60.00   16.22   1.167    0/1
#>  3:   PASS   98.64      13       1   0.500   45.94    7.59  -0.887    0/1
#>  4:   PASS  123.60      22       1   0.500   27.42    7.73  -5.337    0/1
#>  5:   PASS 1208.64      97       1   0.500   60.00   13.43  -0.462    0/1
#>  6:   PASS  151.64     107       1   0.500   60.00    1.50  -8.904    0/1
#>  7:   PASS  142.64     169       1   0.500   60.00    0.86 -15.246    0/1
#>  8:   PASS 1632.06      38       2    1.00   60.00   28.65  -0.698    1/1
#>  9:   PASS 1299.64     102       1   0.500   60.00   13.00  -1.026    0/1
#> 10:   PASS  741.64      70       1   0.500   60.00   10.59   0.399    0/1
#> 11:   PASS  178.64      16       1   0.500   60.00   13.74   2.530    0/1
#> 12:   PASS  511.06      13       2    1.00   60.00   33.52   2.529    1/1
#> 13:   PASS  986.03      36       2    1.00   57.07   29.13   0.606    1/1
#>         AD sample_DP     GQ
#>     <char>    <char> <char>
#>  1:    5,8        13     99
#>  2:  20,19        39     99
#>  3:    9,4        13     99
#>  4:   12,4        16     99
#>  5:  47,43        90     99
#>  6:  81,20       101     99
#>  7: 136,29       165     99
#>  8:   0,38        38     99
#>  9:  50,50       100     99
#> 10:  34,36        70     99
#> 11:    8,5        13     99
#> 12:   0,12        12     36
#> 13:   0,28        28     84

## Pathogenic-only, protein-coding:
gvr_filter(gvr, ABraOM_AF = NULL,
    clin_sig_terms = c("pathogenic", "likely_pathogenic"),
    biotype_keep = "protein_coding", verbose = FALSE)
#>    Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position
#>         <char>         <char> <char>     <char>     <char>          <num>
#> 1:      ATAD3B              0      .     GRCh38       chr1        1479267
#> 2:        ESPN              0      .     GRCh38       chr1        6445777
#> 3:    PRAMEF33              0      .     GRCh38       chr1       13306220
#> 4:    C11orf21              0      .     GRCh38      chr11        2301755
#> 5:       BRCA2              0      .     GRCh38      chr13       32338482
#> 6:      ZNF747              0      .     GRCh38      chr16       30534677
#>    End_Position Strand Variant_Classification Variant_Type Reference_Allele
#>           <num> <char>                 <char>       <char>           <char>
#> 1:      1479267      +                 Intron          SNP                C
#> 2:      6445781      +        Frame_Shift_Del          DEL            AGCTT
#> 3:     13306221      +        Frame_Shift_Ins          INS                -
#> 4:      2301755      +            Splice_Site          SNP                C
#> 5:     32338482      +      Missense_Mutation          SNP                G
#> 6:     30534677      + Translation_Start_Site          SNP                C
#>    Tumor_Seq_Allele1 Tumor_Seq_Allele2         Genotype     dbSNP_RS
#>               <char>            <char>           <char>       <char>
#> 1:                 C                 T              C/T             
#> 2:             AGCTT                 -          AGCTT/-  rs753994746
#> 3:                 -    GGCCCAGAAGGTTC -/GGCCCAGAAGGTTC rs1553122356
#> 4:                 C                 G              C/G             
#> 5:                 G                 C              G/C             
#> 6:                 C                 T              C/T  rs933944167
#>    Tumor_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2
#>                  <char>                 <char>                 <char>
#> 1:            Sample_01                                              
#> 2:            Sample_01                                              
#> 3:            Sample_01                                              
#> 4:            Sample_01                                              
#> 5:            Sample_01                                              
#> 6:            Sample_01                                              
#>             HGVSc                HGVSp  HGVSp_Short   Transcript_ID
#>            <char>               <char>       <char>          <char>
#> 1:   c.444+159C>T                                   ENST00000673477
#> 2: c.1306_1310del p.(Ser436ProfsTer31) p.S436Pfs*31 ENST00000645284
#> 3:   c.268_281dup  p.(Trp98PhefsTer37)  p.W98Ffs*37 ENST00000437300
#> 4:      c.53+1G>C                                   ENST00000381153
#> 5:      c.4127G>C       p.(Gly1376Ala)     p.G1376A ENST00000380152
#> 6:         c.3G>A            p.(Met1?)        p.M1? ENST00000693075
#>             Consequence t_depth t_ref_count t_alt_count         Allele   IMPACT
#>                  <char>  <char>      <char>      <char>         <char>   <char>
#> 1:       intron_variant      13           9           4              T MODIFIER
#> 2:   frameshift_variant      16          13           3              -     HIGH
#> 3:   frameshift_variant      16          12           4 GGCCCAGAAGGTTC     HIGH
#> 4: splice_donor_variant     101          81          20              G     HIGH
#> 5:     missense_variant     165         136          29              C MODERATE
#> 6:           start_lost      70          34          36              T     HIGH
#>      SYMBOL            Gene Feature_type         Feature        BIOTYPE   EXON
#>      <char>          <char>       <char>          <char>         <char> <char>
#> 1:   ATAD3B ENSG00000160072   Transcript ENST00000673477 protein_coding       
#> 2:     ESPN ENSG00000187017   Transcript ENST00000645284 protein_coding   7/13
#> 3: PRAMEF33 ENSG00000237700   Transcript ENST00000437300 protein_coding    2/4
#> 4: C11orf21 ENSG00000110665   Transcript ENST00000381153 protein_coding       
#> 5:    BRCA2 ENSG00000139618   Transcript ENST00000380152 protein_coding  11/27
#> 6:   ZNF747 ENSG00000169955   Transcript ENST00000693075 protein_coding    1/3
#>    INTRON cDNA_position CDS_position Protein_position Amino_acids
#>    <char>        <char>       <char>           <char>      <char>
#> 1:   4/15                                                        
#> 2:            1486-1490    1306-1310          436-437        SF/X
#> 3:              339-340      266-267               89    V/VAQKVX
#> 4:    1/3                                                        
#> 5:                 4326         4127             1376         G/A
#> 6:                  182            3                1         M/I
#>                   Codons     Existing_variation DISTANCE STRAND  FLAGS
#>                   <char>                 <char>   <char> <char> <char>
#> 1:                                                            1       
#> 2:              AGCTTc/c            rs753994746               1       
#> 3: gtg/gtGGCCCAGAAGGTTCg           rs1553122356               1       
#> 4:                                                           -1       
#> 5:               gGa/gCa CD1413320&COSV66459147               1       
#> 6:               atG/atA            rs933944167              -1       
#>    VARIANT_CLASS SYMBOL_SOURCE    HGNC_ID CANONICAL        MANE    MANE_SELECT
#>           <char>        <char>     <char>    <char>      <char>         <char>
#> 1:           SNV          HGNC HGNC:24007       YES MANE_Select    NM_031921.6
#> 2:      deletion          HGNC HGNC:13281       YES MANE_Select    NM_031475.3
#> 3:     insertion          HGNC HGNC:49193       YES MANE_Select NM_001291381.1
#> 4:           SNV          HGNC HGNC:13231       YES MANE_Select NM_001329958.2
#> 5:           SNV          HGNC  HGNC:1101       YES MANE_Select    NM_000059.4
#> 6:           SNV          HGNC HGNC:28350       YES MANE_Select NM_001305018.2
#>    MANE_PLUS_CLINICAL    TSL APPRIS        CCDS            ENSP     SWISSPROT
#>                <char> <char> <char>      <char>          <char>        <char>
#> 1:                               P1    CCDS30.1 ENSP00000500094    Q5T9A4.167
#> 2:                               P1    CCDS70.1 ENSP00000496593    B1AK53.128
#> 3:                         1     P1 CCDS85928.1 ENSP00000492439 A0A0G2JMD5.51
#> 4:                         1     A2 CCDS86168.1 ENSP00000370545    Q9P2W6.113
#> 5:                         5     A2  CCDS9344.1 ENSP00000369497    P51587.242
#> 6:                               A2 CCDS92140.1 ENSP00000509633              
#>           TREMBL       UNIPARC UNIPROT_ISOFORM GENE_PHENO
#>           <char>        <char>          <char>     <char>
#> 1:               UPI000013E044        Q5T9A4-1           
#> 2:               UPI000013D2B6        B1AK53-1          1
#> 3:               UPI000442CEFE                           
#> 4:               UPI0000127A63                           
#> 5:               UPI00001FCBCC                          1
#> 6: A0A8I5KWK6.11 UPI000004CC0E                           
#>                              SIFT     PolyPhen
#>                            <char>       <char>
#> 1:                                            
#> 2:                                            
#> 3:                                            
#> 4:                                            
#> 5:                tolerated(0.18) benign(0.29)
#> 6: tolerated_low_confidence(0.31)   unknown(0)
#>                                                                                                                                      DOMAINS
#>                                                                                                                                       <char>
#> 1:                                                                                                                                          
#> 2: PANTHER:PTHR24153&Low_complexity_(Seg):seg&Prints:PR01217&MobiDB_lite:mobidb-lite&MobiDB_lite:mobidb-lite&AFDB-ENSP_mappings:AF-B1AK53-F1
#> 3:                                                                   PANTHER:PTHR14224&PIRSF:PIRSF038286&AFDB-ENSP_mappings:AF-A0A0G2JMD5-F1
#> 4:                                                                                                                                          
#> 5:                                                                                                       PANTHER:PTHR11289&PIRSF:PIRSF002397
#> 6:                                                                                                                   MobiDB_lite:mobidb-lite
#>     miRNA HGVS_OFFSET     AF AFR_AF AMR_AF EAS_AF EUR_AF SAS_AF gnomADe_AF
#>    <char>      <char> <char> <char> <char> <char> <char> <char>     <char>
#> 1:                                                                        
#> 2:                                                                0.002621
#> 3:                 15                                            0.0003576
#> 4:                                                                        
#> 5:                                                                        
#> 6:                                                               1.443e-06
#>    gnomADe_AFR_AF gnomADe_AMR_AF gnomADe_ASJ_AF gnomADe_EAS_AF gnomADe_FIN_AF
#>            <char>         <char>         <char>         <char>         <char>
#> 1:                                                                           
#> 2:       0.002249       0.001438        0.00167       0.003474       0.001335
#> 3:        0.01328      0.0005464              0              0              0
#> 4:                                                                           
#> 5:                                                                           
#> 6:              0      2.759e-05              0              0              0
#>    gnomADe_MID_AF gnomADe_NFE_AF gnomADe_REMAINING_AF gnomADe_SAS_AF gnomADg_AF
#>            <char>         <char>               <char>         <char>     <char>
#> 1:                                                                             
#> 2:       0.002595       0.002841             0.003058      0.0008685    0.01296
#> 3:       0.000831      8.943e-06            0.0008676      1.429e-05   0.003988
#> 4:                                                                             
#> 5:                                                                             
#> 6:              0      9.265e-07                    0              0  6.571e-06
#>    gnomADg_AFR_AF gnomADg_AMI_AF gnomADg_AMR_AF gnomADg_ASJ_AF gnomADg_EAS_AF
#>            <char>         <char>         <char>         <char>         <char>
#> 1:                                                                           
#> 2:        0.01518        0.01202        0.01318        0.01669       0.009633
#> 3:        0.01408              0       0.001057              0              0
#> 4:                                                                           
#> 5:                                                                           
#> 6:              0              0              0              0              0
#>    gnomADg_FIN_AF gnomADg_MID_AF gnomADg_NFE_AF gnomADg_REMAINING_AF
#>            <char>         <char>         <char>               <char>
#> 1:                                                                  
#> 2:        0.01315       0.004202        0.01189              0.01238
#> 3:              0              0      4.445e-05             0.002388
#> 4:                                                                  
#> 5:                                                                  
#> 6:              0              0       1.47e-05                    0
#>    gnomADg_SAS_AF    MAX_AF MAX_AF_POPS CLIN_SIG SOMATIC  PHENO   PUBMED
#>            <char>    <char>      <char>   <char>  <char> <char>   <char>
#> 1:                                                                      
#> 2:       0.009645   0.01669 gnomADg_ASJ                         33968136
#> 3:              0   0.01408 gnomADg_AFR                                 
#> 4:                                                                      
#> 5:                                                   1&1    1&1         
#> 6:              0 2.759e-05 gnomADe_AMR                                 
#>    MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE TRANSCRIPTION_FACTORS
#>        <char>    <char>       <char>             <char>                <char>
#> 1:                                                                       <NA>
#> 2:                                                                       <NA>
#> 3:                                                                       <NA>
#> 4:                                                                       <NA>
#> 5:                                                                       <NA>
#> 6:                                                                       <NA>
#>    FILTER   QUAL INFO_DP INFO_AC INFO_AF INFO_MQ INFO_QD  CNN_1D     GT     AD
#>    <char> <char>  <char>  <char>  <char>  <char>  <char>  <char> <char> <char>
#> 1:   PASS  98.64      13       1   0.500   45.94    7.59  -0.887    0/1    9,4
#> 2:   PASS  79.60      17       1   0.500   55.44    4.98  -8.136    0/1   13,3
#> 3:   PASS 123.60      22       1   0.500   27.42    7.73  -5.337    0/1   12,4
#> 4:   PASS 151.64     107       1   0.500   60.00    1.50  -8.904    0/1  81,20
#> 5:   PASS 142.64     169       1   0.500   60.00    0.86 -15.246    0/1 136,29
#> 6:   PASS 741.64      70       1   0.500   60.00   10.59   0.399    0/1  34,36
#>    sample_DP     GQ
#>       <char> <char>
#> 1:        13     99
#> 2:        16     87
#> 3:        16     99
#> 4:       101     99
#> 5:       165     99
#> 6:        70     99

## Remove benign annotations (including likely_benign and compound entries):
gvr_filter(gvr, ABraOM_AF = NULL, remove_benign = TRUE, verbose = FALSE)
#>    Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position
#>         <char>         <char> <char>     <char>     <char>          <num>
#> 1:      ATAD3B              0      .     GRCh38       chr1        1479267
#> 2:        ESPN              0      .     GRCh38       chr1        6445777
#> 3:    PRAMEF33              0      .     GRCh38       chr1       13306220
#> 4:    C11orf21              0      .     GRCh38      chr11        2301755
#> 5:       MUC19              0      .     GRCh38      chr12       40485636
#> 6:       BRCA2              0      .     GRCh38      chr13       32338482
#> 7:      ZNF747              0      .     GRCh38      chr16       30534677
#>    End_Position Strand Variant_Classification Variant_Type Reference_Allele
#>           <num> <char>                 <char>       <char>           <char>
#> 1:      1479267      +                 Intron          SNP                C
#> 2:      6445781      +        Frame_Shift_Del          DEL            AGCTT
#> 3:     13306221      +        Frame_Shift_Ins          INS                -
#> 4:      2301755      +            Splice_Site          SNP                C
#> 5:     40485636      +       Nonstop_Mutation          SNP                A
#> 6:     32338482      +      Missense_Mutation          SNP                G
#> 7:     30534677      + Translation_Start_Site          SNP                C
#>    Tumor_Seq_Allele1 Tumor_Seq_Allele2         Genotype     dbSNP_RS
#>               <char>            <char>           <char>       <char>
#> 1:                 C                 T              C/T             
#> 2:             AGCTT                 -          AGCTT/-  rs753994746
#> 3:                 -    GGCCCAGAAGGTTC -/GGCCCAGAAGGTTC rs1553122356
#> 4:                 C                 G              C/G             
#> 5:                 A                 T              A/T  rs192078109
#> 6:                 G                 C              G/C             
#> 7:                 C                 T              C/T  rs933944167
#>    Tumor_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2
#>                  <char>                 <char>                 <char>
#> 1:            Sample_01                                              
#> 2:            Sample_01                                              
#> 3:            Sample_01                                              
#> 4:            Sample_01                                              
#> 5:            Sample_01                                              
#> 6:            Sample_01                                              
#> 7:            Sample_01                                              
#>             HGVSc                 HGVSp   HGVSp_Short   Transcript_ID
#>            <char>                <char>        <char>          <char>
#> 1:   c.444+159C>T                                     ENST00000673477
#> 2: c.1306_1310del  p.(Ser436ProfsTer31)  p.S436Pfs*31 ENST00000645284
#> 3:   c.268_281dup   p.(Trp98PhefsTer37)   p.W98Ffs*37 ENST00000437300
#> 4:      c.53+1G>C                                     ENST00000381153
#> 5:     c.12684A>T p.(Ter4228CysextTer?) p.*4228Cext*? ENST00000454784
#> 6:      c.4127G>C        p.(Gly1376Ala)      p.G1376A ENST00000380152
#> 7:         c.3G>A             p.(Met1?)         p.M1? ENST00000693075
#>             Consequence t_depth t_ref_count t_alt_count         Allele   IMPACT
#>                  <char>  <char>      <char>      <char>         <char>   <char>
#> 1:       intron_variant      13           9           4              T MODIFIER
#> 2:   frameshift_variant      16          13           3              -     HIGH
#> 3:   frameshift_variant      16          12           4 GGCCCAGAAGGTTC     HIGH
#> 4: splice_donor_variant     101          81          20              G     HIGH
#> 5:            stop_lost      37          30           7              T     HIGH
#> 6:     missense_variant     165         136          29              C MODERATE
#> 7:           start_lost      70          34          36              T     HIGH
#>      SYMBOL            Gene Feature_type         Feature            BIOTYPE
#>      <char>          <char>       <char>          <char>             <char>
#> 1:   ATAD3B ENSG00000160072   Transcript ENST00000673477     protein_coding
#> 2:     ESPN ENSG00000187017   Transcript ENST00000645284     protein_coding
#> 3: PRAMEF33 ENSG00000237700   Transcript ENST00000437300     protein_coding
#> 4: C11orf21 ENSG00000110665   Transcript ENST00000381153     protein_coding
#> 5:    MUC19 ENSG00000205592   Transcript ENST00000454784 protein_coding_LoF
#> 6:    BRCA2 ENSG00000139618   Transcript ENST00000380152     protein_coding
#> 7:   ZNF747 ENSG00000169955   Transcript ENST00000693075     protein_coding
#>      EXON INTRON cDNA_position CDS_position Protein_position Amino_acids
#>    <char> <char>        <char>       <char>           <char>      <char>
#> 1:          4/15                                                        
#> 2:   7/13            1486-1490    1306-1310          436-437        SF/X
#> 3:    2/4              339-340      266-267               89    V/VAQKVX
#> 4:           1/3                                                        
#> 5: 56/173                12684        12684             4228         */C
#> 6:  11/27                 4326         4127             1376         G/A
#> 7:    1/3                  182            3                1         M/I
#>                   Codons     Existing_variation DISTANCE STRAND  FLAGS
#>                   <char>                 <char>   <char> <char> <char>
#> 1:                                                            1       
#> 2:              AGCTTc/c            rs753994746               1       
#> 3: gtg/gtGGCCCAGAAGGTTCg           rs1553122356               1       
#> 4:                                                           -1       
#> 5:               tgA/tgT            rs192078109               1       
#> 6:               gGa/gCa CD1413320&COSV66459147               1       
#> 7:               atG/atA            rs933944167              -1       
#>    VARIANT_CLASS SYMBOL_SOURCE    HGNC_ID CANONICAL        MANE    MANE_SELECT
#>           <char>        <char>     <char>    <char>      <char>         <char>
#> 1:           SNV          HGNC HGNC:24007       YES MANE_Select    NM_031921.6
#> 2:      deletion          HGNC HGNC:13281       YES MANE_Select    NM_031475.3
#> 3:     insertion          HGNC HGNC:49193       YES MANE_Select NM_001291381.1
#> 4:           SNV          HGNC HGNC:13231       YES MANE_Select NM_001329958.2
#> 5:           SNV          HGNC HGNC:14362       YES                           
#> 6:           SNV          HGNC  HGNC:1101       YES MANE_Select    NM_000059.4
#> 7:           SNV          HGNC HGNC:28350       YES MANE_Select NM_001305018.2
#>    MANE_PLUS_CLINICAL    TSL APPRIS        CCDS            ENSP     SWISSPROT
#>                <char> <char> <char>      <char>          <char>        <char>
#> 1:                               P1    CCDS30.1 ENSP00000500094    Q5T9A4.167
#> 2:                               P1    CCDS70.1 ENSP00000496593    B1AK53.128
#> 3:                         1     P1 CCDS85928.1 ENSP00000492439 A0A0G2JMD5.51
#> 4:                         1     A2 CCDS86168.1 ENSP00000370545    Q9P2W6.113
#> 5:                         5     P1             ENSP00000508949              
#> 6:                         5     A2  CCDS9344.1 ENSP00000369497    P51587.242
#> 7:                               A2 CCDS92140.1 ENSP00000509633              
#>           TREMBL       UNIPARC UNIPROT_ISOFORM GENE_PHENO
#>           <char>        <char>          <char>     <char>
#> 1:               UPI000013E044        Q5T9A4-1           
#> 2:               UPI000013D2B6        B1AK53-1          1
#> 3:               UPI000442CEFE                           
#> 4:               UPI0000127A63                           
#> 5:                                                       
#> 6:               UPI00001FCBCC                          1
#> 7: A0A8I5KWK6.11 UPI000004CC0E                           
#>                              SIFT     PolyPhen
#>                            <char>       <char>
#> 1:                                            
#> 2:                                            
#> 3:                                            
#> 4:                                            
#> 5:                                            
#> 6:                tolerated(0.18) benign(0.29)
#> 7: tolerated_low_confidence(0.31)   unknown(0)
#>                                                                                                                                      DOMAINS
#>                                                                                                                                       <char>
#> 1:                                                                                                                                          
#> 2: PANTHER:PTHR24153&Low_complexity_(Seg):seg&Prints:PR01217&MobiDB_lite:mobidb-lite&MobiDB_lite:mobidb-lite&AFDB-ENSP_mappings:AF-B1AK53-F1
#> 3:                                                                   PANTHER:PTHR14224&PIRSF:PIRSF038286&AFDB-ENSP_mappings:AF-A0A0G2JMD5-F1
#> 4:                                                                                                                                          
#> 5:                                                                                                                                          
#> 6:                                                                                                       PANTHER:PTHR11289&PIRSF:PIRSF002397
#> 7:                                                                                                                   MobiDB_lite:mobidb-lite
#>     miRNA HGVS_OFFSET     AF AFR_AF AMR_AF EAS_AF EUR_AF SAS_AF gnomADe_AF
#>    <char>      <char> <char> <char> <char> <char> <char> <char>     <char>
#> 1:                                                                        
#> 2:                                                                0.002621
#> 3:                 15                                            0.0003576
#> 4:                                                                        
#> 5:                                                                0.004832
#> 6:                                                                        
#> 7:                                                               1.443e-06
#>    gnomADe_AFR_AF gnomADe_AMR_AF gnomADe_ASJ_AF gnomADe_EAS_AF gnomADe_FIN_AF
#>            <char>         <char>         <char>         <char>         <char>
#> 1:                                                                           
#> 2:       0.002249       0.001438        0.00167       0.003474       0.001335
#> 3:        0.01328      0.0005464              0              0              0
#> 4:                                                                           
#> 5:       0.003774        0.01013       0.007724        0.00188              0
#> 6:                                                                           
#> 7:              0      2.759e-05              0              0              0
#>    gnomADe_MID_AF gnomADe_NFE_AF gnomADe_REMAINING_AF gnomADe_SAS_AF gnomADg_AF
#>            <char>         <char>               <char>         <char>     <char>
#> 1:                                                                             
#> 2:       0.002595       0.002841             0.003058      0.0008685    0.01296
#> 3:       0.000831      8.943e-06            0.0008676      1.429e-05   0.003988
#> 4:                                                                             
#> 5:       0.005085       0.004811             0.005435       0.005452    0.03926
#> 6:                                                                             
#> 7:              0      9.265e-07                    0              0  6.571e-06
#>    gnomADg_AFR_AF gnomADg_AMI_AF gnomADg_AMR_AF gnomADg_ASJ_AF gnomADg_EAS_AF
#>            <char>         <char>         <char>         <char>         <char>
#> 1:                                                                           
#> 2:        0.01518        0.01202        0.01318        0.01669       0.009633
#> 3:        0.01408              0       0.001057              0              0
#> 4:                                                                           
#> 5:        0.03665        0.03018        0.04658        0.06261        0.02314
#> 6:                                                                           
#> 7:              0              0              0              0              0
#>    gnomADg_FIN_AF gnomADg_MID_AF gnomADg_NFE_AF gnomADg_REMAINING_AF
#>            <char>         <char>         <char>               <char>
#> 1:                                                                  
#> 2:        0.01315       0.004202        0.01189              0.01238
#> 3:              0              0      4.445e-05             0.002388
#> 4:                                                                  
#> 5:        0.02334        0.04444        0.04173               0.0444
#> 6:                                                                  
#> 7:              0              0       1.47e-05                    0
#>    gnomADg_SAS_AF    MAX_AF MAX_AF_POPS CLIN_SIG SOMATIC  PHENO   PUBMED
#>            <char>    <char>      <char>   <char>  <char> <char>   <char>
#> 1:                                                                      
#> 2:       0.009645   0.01669 gnomADg_ASJ                         33968136
#> 3:              0   0.01408 gnomADg_AFR                                 
#> 4:                                                                      
#> 5:        0.04893   0.06261 gnomADg_ASJ                                 
#> 6:                                                   1&1    1&1         
#> 7:              0 2.759e-05 gnomADe_AMR                                 
#>    MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE TRANSCRIPTION_FACTORS
#>        <char>    <char>       <char>             <char>                <char>
#> 1:                                                                       <NA>
#> 2:                                                                       <NA>
#> 3:                                                                       <NA>
#> 4:                                                                       <NA>
#> 5:                                                                       <NA>
#> 6:                                                                       <NA>
#> 7:                                                                       <NA>
#>    FILTER   QUAL INFO_DP INFO_AC INFO_AF INFO_MQ INFO_QD  CNN_1D     GT     AD
#>    <char> <char>  <char>  <char>  <char>  <char>  <char>  <char> <char> <char>
#> 1:   PASS  98.64      13       1   0.500   45.94    7.59  -0.887    0/1    9,4
#> 2:   PASS  79.60      17       1   0.500   55.44    4.98  -8.136    0/1   13,3
#> 3:   PASS 123.60      22       1   0.500   27.42    7.73  -5.337    0/1   12,4
#> 4:   PASS 151.64     107       1   0.500   60.00    1.50  -8.904    0/1  81,20
#> 5:   PASS 196.64      37       1   0.500   57.92    5.31  -5.636    0/1   30,7
#> 6:   PASS 142.64     169       1   0.500   60.00    0.86 -15.246    0/1 136,29
#> 7:   PASS 741.64      70       1   0.500   60.00   10.59   0.399    0/1  34,36
#>    sample_DP     GQ
#>       <char> <char>
#> 1:        13     99
#> 2:        16     87
#> 3:        16     99
#> 4:       101     99
#> 5:        37     99
#> 6:       165     99
#> 7:        70     99

## Keep only protein-altering variants and a gene panel:
gvr_filter(gvr, ABraOM_AF = NULL, vc_nonSyn = TRUE,
    genes = c("TP53", "BRCA1", "BRCA2"), verbose = FALSE)
#>    Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position
#>         <char>         <char> <char>     <char>     <char>          <num>
#> 1:       BRCA2              0      .     GRCh38      chr13       32338482
#>    End_Position Strand Variant_Classification Variant_Type Reference_Allele
#>           <num> <char>                 <char>       <char>           <char>
#> 1:     32338482      +      Missense_Mutation          SNP                G
#>    Tumor_Seq_Allele1 Tumor_Seq_Allele2 Genotype dbSNP_RS Tumor_Sample_Barcode
#>               <char>            <char>   <char>   <char>               <char>
#> 1:                 G                 C      G/C                     Sample_01
#>    Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2     HGVSc          HGVSp
#>                    <char>                 <char>    <char>         <char>
#> 1:                                               c.4127G>C p.(Gly1376Ala)
#>    HGVSp_Short   Transcript_ID      Consequence t_depth t_ref_count t_alt_count
#>         <char>          <char>           <char>  <char>      <char>      <char>
#> 1:    p.G1376A ENST00000380152 missense_variant     165         136          29
#>    Allele   IMPACT SYMBOL            Gene Feature_type         Feature
#>    <char>   <char> <char>          <char>       <char>          <char>
#> 1:      C MODERATE  BRCA2 ENSG00000139618   Transcript ENST00000380152
#>           BIOTYPE   EXON INTRON cDNA_position CDS_position Protein_position
#>            <char> <char> <char>        <char>       <char>           <char>
#> 1: protein_coding  11/27                 4326         4127             1376
#>    Amino_acids  Codons     Existing_variation DISTANCE STRAND  FLAGS
#>         <char>  <char>                 <char>   <char> <char> <char>
#> 1:         G/A gGa/gCa CD1413320&COSV66459147               1       
#>    VARIANT_CLASS SYMBOL_SOURCE   HGNC_ID CANONICAL        MANE MANE_SELECT
#>           <char>        <char>    <char>    <char>      <char>      <char>
#> 1:           SNV          HGNC HGNC:1101       YES MANE_Select NM_000059.4
#>    MANE_PLUS_CLINICAL    TSL APPRIS       CCDS            ENSP  SWISSPROT
#>                <char> <char> <char>     <char>          <char>     <char>
#> 1:                         5     A2 CCDS9344.1 ENSP00000369497 P51587.242
#>    TREMBL       UNIPARC UNIPROT_ISOFORM GENE_PHENO            SIFT     PolyPhen
#>    <char>        <char>          <char>     <char>          <char>       <char>
#> 1:        UPI00001FCBCC                          1 tolerated(0.18) benign(0.29)
#>                                DOMAINS  miRNA HGVS_OFFSET     AF AFR_AF AMR_AF
#>                                 <char> <char>      <char> <char> <char> <char>
#> 1: PANTHER:PTHR11289&PIRSF:PIRSF002397                                        
#>    EAS_AF EUR_AF SAS_AF gnomADe_AF gnomADe_AFR_AF gnomADe_AMR_AF gnomADe_ASJ_AF
#>    <char> <char> <char>     <char>         <char>         <char>         <char>
#> 1:                                                                             
#>    gnomADe_EAS_AF gnomADe_FIN_AF gnomADe_MID_AF gnomADe_NFE_AF
#>            <char>         <char>         <char>         <char>
#> 1:                                                            
#>    gnomADe_REMAINING_AF gnomADe_SAS_AF gnomADg_AF gnomADg_AFR_AF gnomADg_AMI_AF
#>                  <char>         <char>     <char>         <char>         <char>
#> 1:                                                                             
#>    gnomADg_AMR_AF gnomADg_ASJ_AF gnomADg_EAS_AF gnomADg_FIN_AF gnomADg_MID_AF
#>            <char>         <char>         <char>         <char>         <char>
#> 1:                                                                           
#>    gnomADg_NFE_AF gnomADg_REMAINING_AF gnomADg_SAS_AF MAX_AF MAX_AF_POPS
#>            <char>               <char>         <char> <char>      <char>
#> 1:                                                                      
#>    CLIN_SIG SOMATIC  PHENO PUBMED MOTIF_NAME MOTIF_POS HIGH_INF_POS
#>      <char>  <char> <char> <char>     <char>    <char>       <char>
#> 1:              1&1    1&1                                         
#>    MOTIF_SCORE_CHANGE TRANSCRIPTION_FACTORS FILTER   QUAL INFO_DP INFO_AC
#>                <char>                <char> <char> <char>  <char>  <char>
#> 1:                                     <NA>   PASS 142.64     169       1
#>    INFO_AF INFO_MQ INFO_QD  CNN_1D     GT     AD sample_DP     GQ
#>     <char>  <char>  <char>  <char> <char> <char>    <char> <char>
#> 1:   0.500   60.00    0.86 -15.246    0/1 136,29       165     99
```
