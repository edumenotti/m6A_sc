# Nextflow Pipeline Completion (Scripts 09–12) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire scripts 09–12 into the Nextflow pipeline, completing the full orchestrated workflow from raw counts through final annotation figures.

**Architecture:** Four new `.nf` module files following the existing pattern in `pipeline/modules/`. Scripts 09+10 form a human-in-the-loop pair (same pattern as scripts 13+14): MANUAL_MARKER_REVIEW produces diagnostics, a human-filled TSV config triggers APPLY_MANUAL_ANNOTATION. Scripts 12 and 11 follow unconditionally. All new processes are wired into `pipeline/main.nf` after RECONCILE_ANNOTATIONS (08). Since the annotation map `pipeline/config/manual_annotation_level1_map_leiden_r2.0.tsv` is already committed, APPLY_MANUAL_ANNOTATION runs unconditionally in this dataset.

**Tech Stack:** Nextflow DSL2, pixi, Python

---

## Data Flow

```
RECONCILE_ANNOTATIONS (08) → adata_reconciled.h5ad
  └─► MANUAL_MARKER_REVIEW (09) → review plots + cluster_review_checkpoint.json
      [HUMAN: fill pipeline/config/manual_annotation_level1_map_leiden_r2.0.tsv]
  └─► APPLY_MANUAL_ANNOTATION (10) → adata_manual_level1.h5ad
        └─► SUBSET_RECLUSTER (12) → adata_annotated_final.h5ad
              └─► FINAL_FIGURES (11) → diagnostic figures
```

Note: 09 and 10 both take `adata_reconciled.h5ad` as input independently. The annotation map is a config file (committed), not a process output.

## File Map

| File | Action | Purpose |
|------|--------|---------|
| `pipeline/modules/manual_marker_review.nf` | Create | Nextflow process for script 09 |
| `pipeline/modules/apply_manual_annotation.nf` | Create | Nextflow process for script 10 |
| `pipeline/modules/subset_recluster.nf` | Create | Nextflow process for script 12 |
| `pipeline/modules/final_figures.nf` | Create | Nextflow process for script 11 |
| `pipeline/nextflow.config` | Modify | Add params: `manual_annotation_map`, `annotation_version`, `leiden_resolution_manual` |
| `pipeline/params.yaml` | Modify | Same 3 new params |
| `pipeline/main.nf` | Modify | Add includes + workflow steps for new modules |
| `pipeline/conf/slurm.config` | Modify | Add resource allocations for new processes |

---

### Task 1: Create manual_marker_review.nf

**Files:**
- Create: `pipeline/modules/manual_marker_review.nf`

- [ ] **Step 1: Read an existing medium-weight module to follow the pattern**

Run:
```bash
cat pipeline/modules/reconcile_annotations.nf
```

- [ ] **Step 2: Create the module file**

Create `pipeline/modules/manual_marker_review.nf` with this content:
```nextflow
process MANUAL_MARKER_REVIEW {
    tag "manual_marker_review"
    publishDir "${params.outdir}/09_manual_marker_review", mode: 'copy'
    memory '24 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "*.png"
    path "*.csv"
    path "*.json"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/09_manual_marker_review.py \
        --input ${h5ad} \
        --out . \
        --marker-db ${projectDir}/config/manual_annotation_markers_skull_immune.tsv \
        --cluster-key leiden_r${params.leiden_resolution_manual}
    """
}
```

- [ ] **Step 3: Verify Nextflow can parse the file**

Run:
```bash
nextflow -version 2>/dev/null || echo "nextflow not in PATH — syntax check skipped"
```

If Nextflow is available:
```bash
nextflow inspect pipeline/modules/manual_marker_review.nf 2>&1 | grep -i error || echo "Parsed OK"
```

- [ ] **Step 4: Commit**

```bash
git add pipeline/modules/manual_marker_review.nf
git commit -m "feat: add Nextflow module for script 09 (manual marker review)"
```

---

### Task 2: Create apply_manual_annotation.nf

**Files:**
- Create: `pipeline/modules/apply_manual_annotation.nf`

- [ ] **Step 1: Create the module file**

Create `pipeline/modules/apply_manual_annotation.nf` with this content:
```nextflow
process APPLY_MANUAL_ANNOTATION {
    tag "apply_manual_annotation"
    publishDir "${params.outdir}/10_manual_level1", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    path h5ad
    path annotation_map

    output:
    path "adata_manual_level1.h5ad", emit: h5ad
    path "*.csv"
    path "*.png"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/10_apply_manual_level1_annotation.py \
        --input ${h5ad} \
        --map ${annotation_map} \
        --out . \
        --cluster-key leiden_r${params.leiden_resolution_manual} \
        --annotation-version ${params.annotation_version}
    """
}
```

- [ ] **Step 2: Commit**

```bash
git add pipeline/modules/apply_manual_annotation.nf
git commit -m "feat: add Nextflow module for script 10 (apply manual annotation)"
```

---

### Task 3: Create subset_recluster.nf

**Files:**
- Create: `pipeline/modules/subset_recluster.nf`

- [ ] **Step 1: Create the module file**

Create `pipeline/modules/subset_recluster.nf` with this content:
```nextflow
process SUBSET_RECLUSTER {
    tag "subset_recluster"
    publishDir "${params.outdir}/12_subset_recluster", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "adata_annotated_final.h5ad", emit: h5ad
    path "*.png"
    path "*.csv"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/12_subset_recluster.py \
        --input ${h5ad} \
        --out .
    """
}
```

- [ ] **Step 2: Commit**

```bash
git add pipeline/modules/subset_recluster.nf
git commit -m "feat: add Nextflow module for script 12 (subset recluster)"
```

---

### Task 4: Create final_figures.nf

**Files:**
- Create: `pipeline/modules/final_figures.nf`

- [ ] **Step 1: Create the module file**

Create `pipeline/modules/final_figures.nf` with this content:
```nextflow
process FINAL_FIGURES {
    tag "final_figures"
    publishDir "${params.outdir}/11_final_figures", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    path h5ad

    output:
    path "*.png"
    path "*.csv"
    path "*.md"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/11_final_figures.py \
        --input ${h5ad} \
        --out .
    """
}
```

- [ ] **Step 2: Commit**

```bash
git add pipeline/modules/final_figures.nf
git commit -m "feat: add Nextflow module for script 11 (final figures)"
```

---

### Task 5: Add new params to nextflow.config and params.yaml

**Files:**
- Modify: `pipeline/nextflow.config`
- Modify: `pipeline/params.yaml`

- [ ] **Step 1: Read current nextflow.config params block**

Run:
```bash
grep -n "leiden_resolution\|progenitor" pipeline/nextflow.config | head -10
```

- [ ] **Step 2: Add new params to nextflow.config**

In `pipeline/nextflow.config`, inside the `params { }` block, add after the existing `leiden_resolution` line:

```
    leiden_resolution_manual = 2.0
    annotation_version = "manual_level1_r2.0_2026-05-05"
    manual_annotation_map = "${projectDir}/config/manual_annotation_level1_map_leiden_r2.0.tsv"
```

- [ ] **Step 3: Add the same params to params.yaml**

In `pipeline/params.yaml`, add:
```yaml
leiden_resolution_manual: 2.0
annotation_version: "manual_level1_r2.0_2026-05-05"
manual_annotation_map: "pipeline/config/manual_annotation_level1_map_leiden_r2.0.tsv"
```

- [ ] **Step 4: Verify nextflow.config parses (no unclosed braces)**

Run:
```bash
grep -c '{' pipeline/nextflow.config && grep -c '}' pipeline/nextflow.config
```
Expected: both counts are equal.

- [ ] **Step 5: Commit**

```bash
git add pipeline/nextflow.config pipeline/params.yaml
git commit -m "feat: add leiden_resolution_manual, annotation_version, manual_annotation_map params"
```

---

### Task 6: Wire new modules into main.nf

**Files:**
- Modify: `pipeline/main.nf`

- [ ] **Step 1: Read current main.nf includes and workflow block**

Run:
```bash
cat pipeline/main.nf
```

- [ ] **Step 2: Add includes for the four new modules**

After the existing `include { MARKERS }` line, add:
```nextflow
include { MANUAL_MARKER_REVIEW }    from './modules/manual_marker_review'
include { APPLY_MANUAL_ANNOTATION } from './modules/apply_manual_annotation'
include { SUBSET_RECLUSTER }        from './modules/subset_recluster'
include { FINAL_FIGURES }           from './modules/final_figures'
```

- [ ] **Step 3: Add workflow steps after MARKERS**

In the `workflow { }` block, after the `MARKERS(RECONCILE_ANNOTATIONS.out.h5ad)` line, add:

```nextflow
    /*
     * Post-reconciliation annotation sub-workflow (scripts 09–12):
     *   MANUAL_MARKER_REVIEW  — produces diagnostics for human review
     *   [HUMAN: fill pipeline/config/manual_annotation_level1_map_leiden_r2.0.tsv]
     *   APPLY_MANUAL_ANNOTATION — applies filled map; map is committed so this always runs
     *   SUBSET_RECLUSTER      — B cell / erythroid level2 sub-clustering
     *   FINAL_FIGURES         — curated summary figures for the annotated dataset
     */
    annotation_map_ch = Channel.fromPath(params.manual_annotation_map, checkIfExists: true)

    MANUAL_MARKER_REVIEW(RECONCILE_ANNOTATIONS.out.h5ad)
    APPLY_MANUAL_ANNOTATION(RECONCILE_ANNOTATIONS.out.h5ad, annotation_map_ch)
    SUBSET_RECLUSTER(APPLY_MANUAL_ANNOTATION.out.h5ad)
    FINAL_FIGURES(SUBSET_RECLUSTER.out.h5ad)
```

- [ ] **Step 4: Verify main.nf has balanced braces and expected includes**

Run:
```bash
grep -c 'include {' pipeline/main.nf
```
Expected: `14` (10 original + 4 new)

Run:
```bash
grep -c '{' pipeline/main.nf && grep -c '}' pipeline/main.nf
```
Expected: equal counts (or equal when accounting for string literals).

- [ ] **Step 5: Commit**

```bash
git add pipeline/main.nf
git commit -m "feat: wire scripts 09-12 into main.nf workflow"
```

---

### Task 7: Update slurm.config with resource allocations for new processes

**Files:**
- Modify: `pipeline/conf/slurm.config`

- [ ] **Step 1: Read current slurm.config**

Run:
```bash
cat pipeline/conf/slurm.config
```

- [ ] **Step 2: Add CPU-only resource allocations for new processes**

In `pipeline/conf/slurm.config`, inside the `withName:` block that currently lists `'QC|DOUBLETS|...'`, extend the pipe-separated list to include the new CPU-only processes:

```groovy
    withName: 'QC|DOUBLETS|NORMALIZE|CLUSTER|RECONCILE_ANNOTATIONS|MARKERS|MANUAL_MARKER_REVIEW|APPLY_MANUAL_ANNOTATION|SUBSET_RECLUSTER|FINAL_FIGURES' {
        queue = 'general'
        clusterOptions = '--partition=general'
    }
```

- [ ] **Step 3: Commit**

```bash
git add pipeline/conf/slurm.config
git commit -m "feat: add slurm resource allocations for new processes 09-12"
```

---

## Self-Review

**Spec coverage:**
- ✅ All four scripts (09, 10, 11, 12) have Nextflow modules
- ✅ Human-in-loop pattern documented in main.nf comments
- ✅ New params added to both nextflow.config and params.yaml
- ✅ Slurm profile updated
- ✅ Data flow: 08 → 09/10 → 12 → 11 is correct
- ✅ Annotation map path uses the committed config file

**Type consistency:** All `emit: h5ad` names match what downstream processes receive as `path h5ad`. ✅

**Placeholder scan:** No TBD. All process commands use exact script paths and param names. ✅
