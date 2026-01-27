# IQ Genetics — MC1R (E-locus) Sanger Automation (Gallus gallus)

Windows offline desktop app to batch-analyze Sanger `.ab1` chromatograms against an MC1R reference, call genotypes, flag dirty/mixed templates, and export CSV.

## MVP (implemented)
- Primary input: `.ab1` (ABIF) trace files
- Reference input: FASTA / FNA / FA / TXT (no forced renaming)
- Automatic:
  - AB1 parsing (basecalls, qualities, peak locations, raw traces)
  - Local alignment to reference (Smith–Waterman)
  - Genotype calling at **c.274** and **c.644**
  - Dirty/mixed-template flagging (peak purity + quality metrics)
- UI:
  - Reference file path
  - File list
  - Status/progress
  - Results table
  - Export CSV button (PDF button stubbed for later)

## Build an installer (MSI)
This repo includes a GitHub Actions workflow that builds a **self-contained** Windows app and a single **MSI installer** (runtime bundled).

- Workflow: `.github/workflows/build-installer.yml`
- Output artifact: `IQGenetics_MC1R_Setup.msi`

## Notes
- cDNA coordinates assume the reference sequence is **CDS starting at ATG = c.1**.
  - If you provide a longer NCBI genomic `.fna`, the app automatically attempts to extract the **945 bp ORF** starting at the first ATG that yields a 945 bp ORF (typical MC1R CDS).
  - For maximum robustness, supply the CDS-only reference.

