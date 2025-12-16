# Maize Bioinformatics Scripts

This repository contains scripting assets for maize (and related projects) covering resequencing pipelines, GWAS, population partitioning and discriminant analysis, genomic selection (GS), and web backend management. This document consolidates all `.r`, `.sh`, and `.py` scripts, and provides an overview, file structure, dependencies, usage, examples, and caveats consistent with the actual implementations.

## Overview
- R scripts: GWAS workflows and plotting (GAPIT, CMplot), genetic distance computation and export, population clustering and discriminant analysis, GBLUP/rrBLUP GS workflows.
- Shell scripts: resequencing end-to-end (FASTQ → BAM/CRAM → gVCF → cohort VCF), QC (fastp/FastQC/MultiQC), mapping and variant calling (BWA/GATK), VCF QC (vcftools/plink), TASSEL association analysis and LD, batch/parallel execution (GNU parallel, Slurm).
- Python scripts: Flask backend services and test runners; genomic selection example with ridge regression.

## File Structure
Only author-maintained scripts are listed; third-party dependency folders such as `renv/library` and `.venv/Lib/site-packages` are excluded. Absolute paths are grouped by extension:

- .r (R scripts)
  - `d:\maize1512\01_Eventshorizon\GS\GBLUP\GBLUP.R`
  - `d:\maize1512\01_Eventshorizon\GWAS\CMplot.r`
  - `d:\maize1512\01_Eventshorizon\GWAS\GAPIT\gapit_pip.r`
  - `d:\maize1512\01_Eventshorizon\Genetic_distance\distance.R`
  - `d:\maize1512\01_Eventshorizon\Genetic_distance\final_excel.R`
  - `d:\maize1512\bioweb\k_save\scripts\gblup_pipeline.R`
  - `d:\maize1512\bioweb\k_save\scripts\rrblup_pipeline.R`
  - `d:\maize1512\02_RealProjects\群体划分与判别分析\群体划分.R`
  - `d:\maize1512\02_RealProjects\群体划分与判别分析\判别分析.R`
  - Archived/experimental (not recommended for production):
    - `d:\maize1512\trash\ancestor_plot.R`
    - `d:\maize1512\trash\distance.R`
    - `d:\maize1512\trash\filtersamples_from_hmp.R`
    - `d:\maize1512\trash\get_major_allele_rigorous.R`
    - `d:\maize1512\trash\get_ref_info.R`
    - `d:\maize1512\trash\refAlleles_cal.R`
    - `d:\maize1512\trash\training_and_pre.R`

- .sh (Shell scripts)
  - Resequencing and pipelines:
    - `d:\maize1512\reseq\fq_processing\generate_manifest.sh`
    - `d:\maize1512\reseq\maize_reseq_pipeline_v4\maize_reseq_pipeline.sh`
  - Web deployment:
    - `d:\maize1512\bioweb\bio_web\deploy.sh`
  - TASSEL/GWAS and tooling:
    - `d:\maize1512\02_RealProjects\gapit_project\my-r-project\sh_Script_file\mirror_tutorial.sh`
    - `d:\maize1512\02_RealProjects\gapit_project\my-r-project\sh_Script_file\tassel_LD.sh`
    - `d:\maize1512\02_RealProjects\gapit_project\my-r-project\sh_Script_file\tassel_pipline_v1.sh`
    - `d:\maize1512\02_RealProjects\gapit_project\my-r-project\sh_Script_file\conda镜像配置.sh`
  - QC and GATK:
    - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\admix.sh`
    - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\map.sh`
    - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\map2.sh`
    - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\qc.sh`
    - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\qc_array.sh`
    - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\snp.sh`
    - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\snp_noBQSR.sh`
    - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\vcf_qc.sh`
    - Parallel/chunking:
      - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\fromDesktop\HaplotypeCaller_array.sh`
      - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\fromDesktop\map.sh`
      - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\fromDesktop\parallel_creat_allTasks.sh`
      - `d:\maize1512\02_RealProjects\LH_gapit\lh_gapit\02_scripts\QC流程脚本\fromDesktop\parallel_HaplotypeCaller.sh`

- .py (Python scripts)
  - Backend services and testing:
    - `d:\maize1512\bioweb\bio_web\backend\app.py`
    - `d:\maize1512\bioweb\bio_web\backend\run_tests.py`
    - `d:\maize1512\bioweb\bio_web\run_all_tests.py`
    - Test modules (pytest):
      - `d:\maize1512\bioweb\bio_web\backend\tests\conftest.py`
      - `d:\maize1512\bioweb\bio_web\backend\tests\test_auth.py`
      - `d:\maize1512\bioweb\bio_web\backend\tests\test_gwas.py`
      - `d:\maize1512\bioweb\bio_web\backend\tests\test_scripts.py`
      - `d:\maize1512\bioweb\bio_web\backend\tests\test_traits.py`
  - Genomic selection (Python):
    - `d:\maize1512\02_RealProjects\py_gs\ml\rig.r.py`

## Dependencies and Environment Requirements
- General:
  - `.sh` scripts target Linux/WSL; ensure executable permissions and a correct `PATH`. Some scripts require `GNU parallel` and Slurm (`sbatch`).
  - `.r` scripts require R ≥ 4.1; `renv` or CRAN mirrors recommended. Common packages: `data.table`, `tidyverse`, `snpReady`, `rrBLUP`, `qqman`, `CMplot`, `caret`, `randomForest`, `e1071`, `openxlsx`, `ComplexHeatmap`, `circlize`, `optparse`, `doParallel`, `foreach`.
  - `.py` backend targets Python 3.10+; typical deps: `Flask`, `flask-cors`, `flask-sqlalchemy`, `flask-migrate`, `flask-jwt-extended`, `python-dotenv`. Test requirements are specified in `backend/requirements-test.txt`.
- Resequencing pipeline (`maize_reseq_pipeline.sh`):
  - Tools: `samtools`, `bcftools`, `bgzip`, `tabix`, `fastp`, `bwa` or `bwa-mem2`, `gatk`, `java`, optional `pigz`, optional `GNU parallel`.
  - Inputs: reference `FASTA`, sample `manifest.tsv` (columns `sample lane r1 r2`, etc.).
- Variant calling and VCF QC:
  - `map*.sh`: `bwa`, `samtools`
  - `snp.sh`/`snp_noBQSR.sh`: `gatk`, `java`, `samtools`, optional indexed dbSNP VCF
  - `vcf_qc.sh`: `vcftools`, `plink`, `bcftools`, `bgzip`, `tabix`
- GAPIT/GWAS:
  - R package `GAPIT3` or `GAPIT`; plotting with `qqman`, `CMplot`
- Population partitioning and discriminant analysis:
  - `adegenet`, `poppr`, `ape`, `ggtree`, `treeio`, `caret`, `randomForest`, `e1071`
- Web deployment:
  - `docker`, `docker-compose`; `.env` configuration file

## Usage (Per Script)

- `reseq/fq_processing/generate_manifest.sh`
  - Purpose: Generate a `manifest.tsv` template from a FASTQ folder (paired R1/R2 detection)
  - Run: `bash generate_manifest.sh <input_dir> [output_file] [extension]`
  - Output: `manifest_template.tsv` with `sample lane r1 r2`; supports Illumina naming patterns

- `reseq/maize_reseq_pipeline_v4/maize_reseq_pipeline.sh`
  - Purpose: End-to-end pipeline (`check` / `prep-ref` / `sample(s)` / `joint` / `all`)
  - Key options: `--manifest --ref --outdir --threads --jobs --java-mem-gb [--intervals] [--known-sites] [--cram]`
  - Examples:
    - Check: `./maize_reseq_pipeline.sh check --manifest manifest.tsv --ref ref.fa --outdir out`
    - Prepare reference: `./maize_reseq_pipeline.sh prep-ref --ref ref.fa --outdir out --threads 8`
    - Run all samples in parallel: `./maize_reseq_pipeline.sh samples --manifest manifest.tsv --ref ref.fa --outdir out --jobs 12 --threads 8`
    - Joint genotyping: `./maize_reseq_pipeline.sh joint --manifest manifest.tsv --ref ref.fa --outdir out --jobs 4 --threads 8`
    - One-shot flow: `./maize_reseq_pipeline.sh all --manifest manifest.tsv --ref ref.fa --outdir out --jobs 8 --threads 8`

- `bioweb/bio_web/deploy.sh`
  - Purpose: Local/production Docker deployment, testing, backup/restore
  - Subcommands: `dev | prod | test | stop | restart | logs | clean | backup | restore | help`
  - Example: `bash deploy.sh dev`

- `01_Eventshorizon/GWAS/GAPIT/gapit_pip.r`
  - Purpose: Multi-trait/multi-model GAPIT workflow (optional parallel), aggregation and plotting
  - Run: `Rscript gapit_pip.r --pheno pheno.csv --geno_hmp chr1.hmp.txt,chr2.hmp.txt --outdir GAPIT_out --models GLM,MLM --traits Height`
  - Common options: `--GD/--GM` as an alternative to HapMap; `--covar --kinship --n_pcs --maf --parallel --cores --alpha_fdr --plot_extra`
  - Outputs: per-trait results and combined table `GWAS_AllTraits_AllModels_Combined.csv`, significant hits `result/GWAS_Significant_*.csv`

- `01_Eventshorizon/GWAS/CMplot.r`
  - Purpose: Manhattan/QQ/circular manhattan/SNP density plots; input must contain `Trait/Marker/Chr/Pos/p`
  - Run: load and execute in R; or `Rscript CMplot.r`
  - Key step: standardize column names to `chr/BP/P/SNP` before calling `CMplot`

- `01_Eventshorizon/Genetic_distance/distance.R`
  - Purpose: Compute pairwise genetic distances from HapMap, export Excel and heatmap
  - Run: `Rscript distance.R` (default input `365geno.hmp.txt` inside the script)
  - Outputs: `GeneticDistance_full_365lines.xlsx`, nearest/farthest neighbors `*_nearest20.xlsx` / `*_farthest20.xlsx`, heatmap `GeneticDistance_heatmap.pdf`

- `01_Eventshorizon/Genetic_distance/final_excel.R`
  - Purpose: Format nearest/farthest neighbor lists into a consolidated Excel sheet
  - Run: `Rscript final_excel.R` (defaults to reading files under `GD_output`)
  - Output: `Formatted_Genotype_Distance.xlsx`

- `bioweb/k_save/scripts/gblup_pipeline.R`
  - Purpose: rrBLUP GBLUP CV and independent prediction from CSV inputs
  - Run: `Rscript gblup_pipeline.R` (reads `./dataset/train_data.csv`, `trait_imputed.csv`, `pred_data.csv`)
  - Outputs: `./result/gblup_cv_pred.csv`, `gblup_independent_pred.csv` and metric files

- `bioweb/k_save/scripts/rrblup_pipeline.R`
  - Purpose: HapMap IUPAC recode → numeric → SNPReady QC → rrBLUP GBLUP
  - Run: `Rscript rrblup_pipeline.R` (reads `./dataset/train_data.hmp.txt`, `pred_data.hmp.txt`, phenotype file)
  - Outputs: validation/prediction CSVs and metrics, `GBLUP_SNPReady.RData`

- `01_Eventshorizon/GS/GBLUP/GBLUP.R`
  - Purpose: SNPReady + rrBLUP GBLUP (80/20 validation + independent prediction)
  - Run: `Rscript GBLUP.R` (reads `dataset/train_data.hmp.txt`, `pred_data.hmp.txt`, `trait_imputed.txt`)
  - Outputs: validation/prediction CSVs and metrics

- `02_RealProjects/群体划分与判别分析/群体划分.R`
  - Purpose: Hierarchical clustering with Rogers distance, automatic `k` selection, exports assignments and tree plots
  - Run: `Rscript 群体划分.R` (default input `365geno.hmp.txt`)
  - Outputs: `cluster_assignment_and_k.xlsx`, `01_k_selection_silhouette.*`, `02_hclust_tree_colored.*`

- `02_RealProjects/群体划分与判别分析/判别分析.R`
  - Purpose: RF/SVM discriminant analysis on QC’d genotypes, exports evaluation and plots
  - Run: `Rscript 判别分析.R` (defaults: `365geno.hmp.txt`, `YJ-524.hmp.txt`, `s_to_d.xlsx`)
  - Outputs: confusion matrices, per-class metrics plots, Excel report, predictions `YJ_pre_result.xlsx`

- `02_RealProjects/LH_gapit/lh_gapit/02_scripts/admix.sh`
  - Purpose: Run `admixture --cv` on PLINK `bed/bim/fam` over a range of K in parallel and collect CV error
  - Run: edit header variables then `bash admix.sh`
  - Output: `admixture_logs/logK.out`; the script prints CV errors by K at the end

- `02_RealProjects/LH_gapit/lh_gapit/02_scripts/QC流程脚本/map*.sh`
  - Purpose: BWA MEM mapping → SAM→BAM conversion → sorting → indexing
  - Run: set `THREADS`, `REF_GENOME`, `INPUT_FQ1/2`, `RG_STRING` then `bash map.sh` / `map2.sh`
  - Output: `sample.sorted.bam` and `.bai`

- `02_RealProjects/LH_gapit/lh_gapit/02_scripts/QC流程脚本/qc*.sh`
  - `qc.sh` (standalone): fastp filtering → FastQC → MultiQC aggregation; edit paths/switches then run
  - `qc_array.sh` (Slurm): array jobs driven by `samples.tsv` for fastp/FastQC/MultiQC

- `02_RealProjects/LH_gapit/lh_gapit/02_scripts/QC流程脚本/snp.sh`
  - Purpose: GATK `MarkDuplicates → BQSR → HaplotypeCaller → GenotypeGVCFs → hard-filter → merge → PASS`
  - Run: edit `PROJECT/REF/DBSNP/GATK/THREADS` at the top then execute
  - Output: `*_pass_variants.vcf.gz` (with index `*.tbi`)

- `02_RealProjects/LH_gapit/lh_gapit/02_scripts/QC流程脚本/snp_noBQSR.sh`
  - Purpose: Run `HaplotypeCaller` directly from a MarkDuplicates BAM (skip BQSR), supports `-i/-r/-o`
  - Run: `bash snp_noBQSR.sh -i sample.sorted.markdup.bam -r ref.fa -o outdir`

- `02_RealProjects/LH_gapit/lh_gapit/02_scripts/QC流程脚本/vcf_qc.sh`
  - Purpose: SNP-only extraction, sample/site QC, PLINK binary output
  - Run: edit thresholds at the top then `bash vcf_qc.sh`
  - Output: `OUTPUT_PREFIX.{bed,bim,fam}`

- Parallel and chunking (fromDesktop)
  - `HaplotypeCaller_array.sh`: Slurm array to run GATK HC per sample×interval
  - `parallel_creat_allTasks.sh`: build `all_tasks.txt` (sample/interval pairs)
  - `parallel_HaplotypeCaller.sh`: GNU parallel to run HC tasks

- TASSEL and environment
  - `tassel_LD.sh`: per-chromosome sliding-window LD; merges outputs
  - `tassel_pipline_v1.sh`: filter → kinship → MLM association workflow
  - `mirror_tutorial.sh`: Docker/TASSEL usage and terminal theming; includes TASSEL MLM examples
  - `conda镜像配置.sh`: configure Tsinghua mirror channels for conda

- Python backend and tests
  - `bioweb/bio_web/backend/app.py`: Flask API (auth, traits, GWAS, scripts). Run `python app.py`
    - Config: `.env` with `SECRET_KEY`, `JWT_SECRET_KEY`, `DATABASE_URL` (defaults to `sqlite:///maize_bio.db`)
  - `backend/run_tests.py`: install test requirements and run pytest with coverage; `python backend/run_tests.py`
  - `run_all_tests.py`: combined frontend/backend test runner and coverage report; `python run_all_tests.py`
  - `backend/tests/*`: pytest modules used by the runners, not executed directly

- Python GS example
  - `02_RealProjects/py_gs/ml/rig.r.py`: ridge regression GS with Monte Carlo CV; inputs are HapMap and phenotype TSV; run `python rig.r.py`

## Example Usage
- Build manifest and run full resequencing pipeline:
  - `bash reseq/fq_processing/generate_manifest.sh data/fastq manifest.tsv`
  - `bash reseq/maize_reseq_pipeline_v4/maize_reseq_pipeline.sh all --manifest manifest.tsv --ref ref.fa --outdir out --jobs 8 --threads 8`
- Run GAPIT end-to-end:
  - `Rscript 01_Eventshorizon/GWAS/GAPIT/gapit_pip.r --pheno traits.csv --geno_hmp chr1.hmp.txt,chr2.hmp.txt --outdir GAPIT_out --models GLM,MLM --traits Height`
- VCF QC and PLINK output:
  - `bash 02_RealProjects/LH_gapit/lh_gapit/02_scripts/QC流程脚本/vcf_qc.sh`
- Web backend (dev):
  - `bash bioweb/bio_web/deploy.sh dev` → frontend `http://localhost:3000`, backend `http://localhost:5000`
- Python GS:
  - `python 02_RealProjects/py_gs/ml/rig.r.py`

## Notes
- Platform and paths:
  - `.sh` scripts should run on Linux/WSL; Windows path examples must be adapted to Linux mount points (e.g., `/mnt/d/...`).
  - Some R scripts include `setwd(...)` and hard-coded filenames; adjust or remove as needed.
- Resources and parallelism:
  - GATK/HaplotypeCaller, joint genotyping, and parallel tasks require sufficient CPU/RAM and temp storage; tune `GNU parallel`/Slurm job settings to your hardware.
- Indexes and references:
  - Reference FASTA must have `*.fai` and `*.dict`; build `bwa`/`bwa-mem2` indexes beforehand.
  - `vcf_qc.sh` expects compressed/indexed VCF (`bgzip` + `tabix`); it auto-handles plain `.vcf`.
- Dependencies:
  - Use reliable CRAN/Bioc mirrors for R packages; `impute` requires `BiocManager` first.
  - Backend test deps: install via `pip install -r backend/requirements-test.txt` before running tests.
- Data consistency:
  - HapMap must include `rs# / alleles / chrom / pos`; scripts contain strict checks and recoding (IUPAC → diploid bases → 0/1/2).
  - GAPIT scripts switching between HapMap and GD/GM require aligned column names and sample sets.
- Known issues/limitations:
  - `fromDesktop/map.sh` contains duplicated fragments and stray backticks near the end; prefer `QC流程脚本/map.sh` or clean it up.
  - `GBLUP.R`/`rrblup_pipeline.R` include hard-coded paths and trait column names (`Trait02`, `fenzhi_BLUP`); update to your dataset.
  - `snp.sh` requires `--known-sites` dbSNP VCF with index.
  - `admix.sh` requires `GNU parallel` and `admixture`; PLINK files must share the same prefix.

## Version & Maintenance
- This README reflects the state of scripts as of 2025-12-16. Please update it when scripts or parameters change.
- Scripts under `trash/` are historical or experimental and may not be stable.

