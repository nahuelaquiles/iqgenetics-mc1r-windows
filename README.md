# IQ Genetics — MC1R E-locus Haplotype Analysis (Gallus gallus)

Offline Windows desktop software for analysis of Sanger `.ab1` chromatograms from the IQ Genetics/GECORP MC1R assay.

## Scientific model (v1.0, 2026 update)

The original software classified E status mainly from c.274 (E92K). Version 1.0 replaces that logic with sample-level consensus and the 18 common MC1R haplotypes reported by Ma et al. (PNAS, 2026; doi:10.1073/pnas.2605288123).

Haplotype sites: c.212 T>C (M71T), c.274 G>A (E92K), c.376 G>A (V126I), c.398 T>A/C (L133Q/P), c.409 G>A (A137T), c.427 A>G (T143A), c.637 C>T (R213C), and c.644 A>C (H215P).

The software loads all reads for a batch, groups reads belonging to the same bird, excludes globally dirty chromatograms, builds per-site consensus, and enumerates all compatible diplotypes among the 18 common haplotypes.

## Breeder-facing categories

- `E/E-compatible`: both compatible haplotypes are strict Extended Black E1/E2.
- `E/ALT`: one E1/E2 haplotype plus a different common MC1R haplotype.
- `ALT/ALT`: no strict E1/E2 haplotype is present among the compatible diplotype(s).
- `AMBIGUOUS`: more than one breeder-facing category remains compatible with the unphased Sanger data.
- `UNRESOLVED`: the genotype is not represented by the 18 common haplotypes.
- `NoCall`: no usable sample-level result.

`E/R1` and `E/R2` retain the names used in Ma et al. 2026 but are not used to certify a line as fixed `E/E` in this conservative breeding implementation.

## Input and consensus

Load all `.ab1` reads for all birds together. Common suffixes such as `-F`, `-R`, `-2-P-F`, `-2-P-R`, and `-2-Ext-F` are normalized so reads from one bird are combined before interpretation. A validated MC1R reference is bundled in the installer and selected automatically.

Persistent secondary peaks trigger global QC exclusion. A low Phred score at a truly heterozygous site is not accepted merely by lowering a threshold: the site is rescued only when the whole read passed QC and a strong localized dual peak is present. Discordant reads are resolved only by a 2:1-or-stronger sample-level majority; otherwise that site becomes `NoCall`.

## Output

One row per bird/sample with read count, QC, genotype at all eight 2026 haplotype sites, call status, breeder-facing category, compatible diplotype(s), compatible phenotype class(es), and interpretation.

## Scope

This is a targeted MC1R/E-locus breeding assay. It is not a whole-pigmentation-genome test and does not certify total black phenotype independently of EDN3/fibromelanosis, PMEL, GJA5/Melanotic, phenotype records, or other genetic background effects.

## Build

GitHub Actions builds a self-contained Windows x64 application and WiX MSI. The workflow runs an MC1R haplotype self-test before publishing the installer artifact.
