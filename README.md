# HAT-FLEX
### Tychele N. Turner, Ph.D.
### Washington University in St. Louis

**HAT-FLEX**: Flexible Trio DNV detection on existing VCFs.

HAT-FLEX is a caller-agnostic, drop-in trio DNV detection tool that operates directly on existing VCFs. It introduces allele-level intersection, sex/PAR-aware logic, clustering, comprehensive audit outputs, and streamlined operations. HAT-FLEX supports both trio-level and large multi-sample VCFs, producing tidy, per-child outputs with full provenance. Developed in response to user feedback from HAT, HAT-FLEX enables use of existing VCFs, can extend to non-human species with diploid genomes, improves performance, increases configurability, and provides robust handling of sex chromosomes (X/Y).

> [!NOTE]
> Supports VCFs generated on highly accurate **short-read** and **long-read** sequencing datasets, respectively.

For reference, HAT is our original de novo variant calling tool and it begins with the CRAM or BAM files for each trio members. It is available here: https://github.com/tnturnerlab/HAT. If using **HAT** or **HAT-FLEX**, please cite our HAT paper:

Ng JK, Turner TN. HAT: de novo variant calling for highly accurate short-read and long-read sequencing data. Bioinformatics. 2024 Jan 2;40(1):btad775. doi: 10.1093/bioinformatics/btad775. PMID: 38175776; PMCID: PMC10777354.

Already have DNVs and looking to QC them? Check out our tool acorn at: https://github.com/TNTurnerLab/acorn. If using **acorn**, please cite our paper:

Turner TN. Acorn: an R package for de novo variant analysis. BMC Bioinformatics. 2023 Sep 2;24(1):330. doi: 10.1186/s12859-023-05457-z. PMID: 37660114; PMCID: PMC10475174.

## HAT-FLEX vs HAT (cool differentiators)

## Input & orchestration
- **Multi-sample VCF aware (input):** HAT-FLEX can read multi-sample cohort VCFs and extract each trio on the fly.
- **Single-caller mode:** Works with just `--caller1-vcf`; logs caller1-only per child.
- **No heavy workflow:** Pure Python CLI (optional `pysam`); easy to drop into any stack.

## Intersection & multiallelic precision
- **Allele-level intersection (default):** Splits multiallelic sites and keeps **only the ALT the child carries**; toggleable to locus/site mode.
- **VCF-spec normalization:** Trims REF/ALT common prefix/suffix before keying to avoid false mismatches.

## Sex/PAR & haploid intelligence
- **Built-in hg38 PARs** (override via `--par-bed`).
- **Haploid-aware thresholds:** Separate DP/AB minima for male non-PAR X/Y.
- **Irrelevant parent logic:** Ignores the non-contributing parent at haploid loci and **annotates** which one (`IRRELPARENT`).
- **Unusual haploid flag:** Marks het-like GT or suspicious mid-range AB (`HAPUNUSUAL`) for the male non-PAR X/Y, with smart exceptions for partial haploid GTs.

## Filtering finesse
- **Parent ALT-leak detection:** By AD counts and/or by GT (`./1`, `1/.`, etc.).
- **Strand-bias checks:** Require ALT on both strands; optional Fisher exact test.
- **Homopolymer screen:** Configurable A/T run filter.
- **Min child ALT reads:** Simple knob to tame noisy calls.
- **Flexible FORMAT enforcement:** Require keys (`GT,DP,AD,GQ` by default) with optional **strict** mode and concise error logging.

## Region masking (plug-and-play)
- **Any BEDs:** Point at a **folder** or a **tar(.gz/.bz2)**; auto-loads `.bed`/`.bed.gz`.
- **chr/chrless bridging:** Handles `chrX`/`X` naming seamlessly.

## Clustering context
- **PASS-only sliding-window clustering:** Adds `INFO/CLUSTER=yes|no`.
- **Retro-upgrade logic:** When a new PASS creates a cluster, prior PASS in-window records get upgraded to `CLUSTER=yes` for consistency.

## Output quality & auditability
- **Rich INFO tags:** `DENOVO`, `SEXCTX`, `IRRELPARENT`, `HAPUNUSUAL`, `CLUSTER`.
- **Per-site TSV (optional):** Trio GT/DP/GQ/AB (irrelevant parent blanked).
- **Metrics (TSV/JSON):** Counts by stage, DP/GQ/AB medians, %X/%Y/%PAR, parent-leak sites.
- **Run manifest JSON:** Full provenance (timestamp, args, inputs, outputs) per child.
- **Deterministic sorting:** Honors header `##contig` order for stable diffs.
- **Updates headers:** Auto-adds missing INFO/FILTER declarations once.
- **Indexed outputs:** bgzip+tabix when pysam available; graceful gzip fallback.

## What HAT focuses on instead
- **End-to-end pipeline from BAM/CRAM:** Parabricks (GPU) DeepVariant + GATK HC → GLnexus → simple filters.
- **Workflow tooling:** Snakemake/WDL + Docker; many intermediates (often not retained).

---

### TL;DR
If you want a **turnkey, GPU-accelerated calling pipeline from reads**, HAT is the stack.  
If you already have VCFs and want **sharper, sex-aware, allele-level DNV filtering with strong auditability and minimal ops**, **HAT-FLEX** brings the cool stuff.


## Usage:

### Installation
```
git clone https://github.com/tnturnerlab/HAT-FLEX.git
cd HAT-FLEX/
pip3 install .
```

### Usage to get output similar to HAT, examples below assume build 38 PAR, if using a different genome use the `--par-bed` flag 
(note: please see HAT GitHub for recommended regions files including centromeres, LCR, repeats).

Example input data is available at https://zenodo.org/records/17602491
```
hat-flex --family-file HG03732_family_file.txt \
--caller1-vcf HG03732.trio.glnexus.dv.vcf.gz \
--caller2-vcf HG03732.trio.glnexus.hc.vcf.gz \
--out-dir . \
--verbose \
--require-format-keys GT,DP,AD,GQ \
--normalize-alleles \
--regions regions/ \
--metrics HG03732.metrics.txt \
--default-child-sex male
```

### Usage to keep the failing sites and tag them in the VCF
```
hat-flex --family-file HG03732_family_file.txt \
--caller1-vcf HG03732.trio.glnexus.dv.vcf.gz \
--caller2-vcf HG03732.trio.glnexus.hc.vcf.gz \
--out-dir . \
--verbose \
--require-format-keys GT,DP,AD,GQ \
--normalize-alleles \
--regions regions/ \
--metrics HG03732.metrics.txt \
--emit-sites-tsv \
--keep-failures \
--default-child-sex male
```

### Usage with only one caller
```
hat-flex --family-file HG03732_family_file.txt \
--caller1-vcf HG03732.trio.glnexus.dv.vcf.gz \
--out-dir . \
--verbose \
--require-format-keys GT,DP,AD,GQ \
--normalize-alleles \
--regions regions/ \
--metrics HG03732.metrics.txt \
--emit-sites-tsv \
--keep-failures \
--default-child-sex male
```

### Usage with no region masking
```
hat-flex --family-file HG03732_family_file.txt \
--caller1-vcf HG03732.trio.glnexus.dv.vcf.gz \
--out-dir . \
--verbose \
--require-format-keys GT,DP,AD,GQ \
--normalize-alleles \
--metrics HG03732.metrics.txt \
--emit-sites-tsv \
--keep-failures \
--default-child-sex male
```

## Check out more by running `hat-flex -h`. It will give the message below:

```
Welcome to HAT-FLEX for calling de novo variants. Check out TNTurnerLab GitHub to learn more.

options:
  -h, --help            show this help message and exit
  --family-file FAMILY_FILE
                        This file contains the fatherName,motherName,childName and is optional. If not used the pedigree file will be used instead. (default: None)
  --ped PED             Standard pedigree file of fid, iid, pid, mid, sex, pheno (default: None)
  --prefer-source {ped,family}
                        If both exist, which one is preferred by the user the family-file or the pedigree file (default: ped)
  --caller1-vcf CALLER1_VCF
                        Path to Caller 1 VCF (default: None)
  --caller2-vcf CALLER2_VCF
                        Path to Caller 2 VCF (default: None)
  --intersect-mode {locus,allele}
                        Type of intersection; exact=allele or position=locus (default: allele)
  --normalize-alleles   option to normalize alleles (default: False)
  --regions REGIONS     Directory of regions (as compressed or uncompressed bed files) you want masked in the output (e.g., segmental duplications) (default: None)
  --par-bed PAR_BED     bed file containing coordinates of pseudoautosomal regions on the X and Y chromosome. If this is not provided, the default is b38 of the human genome. (default: None)
  --gq-value GQ_VALUE   Minimum GQ desired for genotypes. Default is 20. (default: 20)
  --depth-value DEPTH_VALUE
                        Minimum DP desired for site. Default is 10. (default: 10)
  --haploid-depth-value HAPLOID_DEPTH_VALUE
                        DP minimum when the child is haploid (male non-PAR X/Y). Also applied to the relevant parent at haploid loci. Default is 5. (default: 5)
  --ab-min AB_MIN       Minimum AB desired for alternate allele. Default is 0.25 (default: 0.25)
  --haploid-ab-min HAPLOID_AB_MIN
                        Minimum AB desired for alternate allele in the haploid state. Default is 0.85 (default: 0.85)
  --haploid-suspicious-ab-low HAPLOID_SUSPICIOUS_AB_LOW
                        Suspicious low AB for the haploid state. Default is 0.35 (default: 0.35)
  --haploid-suspicious-ab-high HAPLOID_SUSPICIOUS_AB_HIGH
                        Suspicious high AB for the haploid state. Default is 0.65 (default: 0.65)
  --default-child-sex {male,female}
                        Child Sex. If not provided, only the autosomes will be run (default: None)
  --require-format-keys REQUIRE_FORMAT_KEYS
                        Format keys required for consideration in the filtering scheme. (default: GT,DP,AD,GQ)
  --strict-format       If this is provided than all the format keys in --require-format-keys must be present in the variant line for it to be considered. (default: False)
  --sb-min-each SB_MIN_EACH
                        SB minimum optional. (default: 1)
  --sb-fisher-pmax SB_FISHER_PMAX
                        Right-tail Fisher exact p-value maximum for strand bias on ALT reads (set to, e.g., 0.001). (default: None)
  --min-child-alt-reads MIN_CHILD_ALT_READS
                        minimum number of reads the alternate allele must be seen in the child. (default: 0)
  --max-parent-alt-reads MAX_PARENT_ALT_READS
                        Maximum alternate allele reads allowed in each parent (AD-based). Overridden for the relevant parent at haploid loci by --haploid-max-parent-alt-reads if provided. (default: 0)
  --haploid-max-parent-alt-reads HAPLOID_MAX_PARENT_ALT_READS
                        If set, overrides --max-parent-alt-reads for the relevant parent at haploid loci (MaleX→mother; MaleY→father). (default: None)
  --max-homopolymer MAX_HOMOPOLYMER
                        Maximum homopolymer length of A or T. Default is 10. (default: 10)
  --parent-gt-alt-triggers-parentalt
                        If set, relevant parent's GT containing alternate allele triggers ParentAlt even if AD missing. (default: True)
  --no-parent-gt-alt-triggers-parentalt
                        Disable GT-based ParentAlt for the relevant parent. (default: True)
  --haploid-suspicious-skip-missing-allele
                        Skip HaploidUnusual when child haploid GT is './ALT' or 'ALT/.' (partial-missing) at haploid loci. (default: True)
  --no-haploid-suspicious-skip-missing-allele
                        Do not skip HaploidUnusual for partially missing haploid GTs. (default: True)
  --cluster-window-bp CLUSTER_WINDOW_BP
                        Clustering window (bp) for PASS variants; sets INFO/CLUSTER on PASS records to yes/no based on window content. (default: 100)
  --cluster-min-count CLUSTER_MIN_COUNT
                        Minimum PASS de novos within the window (including this site) to call CLUSTER=yes (default 2 = requires ≥1 other PASS site). (default: 2)
  --out-dir OUT_DIR     output directory for the results. (default: /path/to/outdir)
  --metrics METRICS     optional metrics file to store critical metrics. (default: None)
  --emit-sites-tsv      optional tsv file of sites (default: False)
  --keep-failures       option to keep the failing sites in the VCF with the listing of why they failed to be a true de novo (default: False)
  --write-index {none,tbi,csi}
                        option to write a tabix-index for the output vcf (default: none)
  --dry-run             dry run of the analysis (default: False)
  -v, --verbose         verbose messaging of processing steps (default: 0)
```

### HAT versus HAT-FLEX (at a glance)

| Aspect | **HAT** | **HAT-FLEX** |
|---|---|---|
| **Input flexibility** | Expects CRAM/BAM files to start and per-family/trio joint VCFs produced by its own pipeline | Accepts **multi-sample VCFs** (or per-trio files), pulls just Father/Mother/Child for each trio |
| **Start point** | From **BAM/CRAM** → callers → GLnexus → filter | From **existing VCFs** → intersect + **rich filtering** (Pure Python, no workflow required) |
| **Intersection** | Site/locus-level intersection of callers. HAT requires two callers. | **Allele-level** (default) **or** locus-level; VCF-spec allele normalization. The intersection is optional and HAT-FLEX can perform DNV calling on only one caller. |
| **Multiallelic handling** | Site-level, no explicit per-ALT logic | **True per-ALT**: splits multiallelics; keeps **only the ALT** the child carries; AD/GT checks keyed to that ALT |
| **Single-caller mode** | Designed around two callers (DeepVariant + GATK HaplotypeCaller) | **Yes**, HAT-FLEX can optionally run with just `--caller1-vcf`; logs **caller1-only** mode per child |
| **Sex/PAR awareness** | Uniform thresholds; no explicit haploid logic | **Male X/Y haploid logic**, PAR detection, **haploid DP/AB** thresholds, `SEXCTX` tag |
| **Irrelevant parent handling** | Not explicit | Ignores irrelevant parent at haploid loci; annotates `IRRELPARENT=F|M` |
| **Parent ALT “leak” checks** | Not explicit | **AD-based** limits and optional **GT-based** triggers (`./1`, `1/.`, etc.) |
| **Strand bias** | Not considered | Requires ALT on both strands if selected as an option; optional **Fisher exact** p-value |
| **Homopolymer screen** | Only checks for 10 As or 10 Ts| **Max homopolymer** filter (configurable) |
| **Region masks** | Prescribed RepeatMasker/LCR/repeats/CpG sets | **Any BEDs** from a **folder or tar(.gz/.bz2)**; easy to swap resources |
| **Clustering** | None | **PASS-only sliding-window** → `INFO/CLUSTER=yes|no` |
| **Annotations** | Standard | Adds `DENOVO`, `SEXCTX`, `IRRELPARENT`, `HAPUNUSUAL`, `CLUSTER` |
| **Auditability** | Filenames encode thresholds | **Metrics** (DP/GQ/AB medians, counts, PAR/X/Y %), **sites.tsv**, **run_manifest.json** |
| **Compression/index** | External tools in workflow | Writes **bgzip+tabix** when `pysam` present (or gzip) |
| **Dependencies** | Parabricks/GLnexus/bcftools/bedtools + Snakemake/WDL/Docker | **Pure Python 3** (+ optional `pysam`) |
| **Output granularity** | Per-child **trio VCF** | Per-child **trio VCF** (even when input was multi-sample) |
| **Number & type of output files (per child)** | **1 “final”:** filtered **DNV VCF**; plus many pipeline intermediates (usually not retained) | **1 to 4 deliverables:** `<child>.final.vcf.gz` (+ optional `.tbi`/`.csi`), optional `sites.tsv`, optional **metrics** (TSV/JSON), and `run_manifest.json` |
