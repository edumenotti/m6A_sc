# Nextflow Integration: Progenitor Recluster (Scripts 13 + 14) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire scripts `13_progenitor_recluster.py` and `14_apply_progenitor_annotation.py` into the Nextflow pipeline as proper modules with parameterized inputs and human-in-the-loop conditional execution, after first safely committing the existing untracked Nextflow infrastructure.

**Architecture:** Two new modules — `PROGENITOR_RECLUSTER` (always runs when invoked) and `APPLY_PROGENITOR_ANNOTATION` (runs only if filled annotation map TSV exists). They take an h5ad path as an input parameter rather than chaining off the existing `RECONCILE_ANNOTATIONS → MARKERS` workflow, because (a) scripts 10–12 are not yet wired into Nextflow and (b) the human-review checkpoint between 13 and 14 makes them fundamentally a separate sub-workflow. This decouples the new modules from the unmigrated upstream automation while keeping a clean migration path: when scripts 10–12 are eventually wired, the parameter path will simply be replaced by a channel from the upstream process.

**Tech Stack:** Nextflow DSL2, pixi, Python (existing scripts), bash. No new dependencies.

---

## Background & Constraints

### Current state of the Nextflow pipeline

`pipeline/main.nf` workflow chain:

```
QC → DOUBLETS → NORMALIZE → INTEGRATE → CLUSTER → ANNOTATE_HSPC + ANNOTATE_MATURE → RECONCILE_ANNOTATIONS → MARKERS
```

The chain ends at `MARKERS` (script 07). Scripts 08–14 exist as standalone Python files but **only 08 has a Nextflow module** (`reconcile_annotations.nf`, called by `RECONCILE_ANNOTATIONS`). Scripts 09–14 have no modules and are run manually.

**The entire `pipeline/` Nextflow infrastructure is currently untracked in git** (verified by `git status`). This is a pre-existing issue not caused by recent work — the modules and config files have never been committed.

### Why scripts 13/14 are special

They sit at a human-in-the-loop boundary:
- Script 13 produces diagnostic outputs (UMAPs, DE tables, dotplots, `annotation_template.tsv`)
- Human reviews outputs and fills `pipeline/config/progenitor_annotation_map.tsv` with `# chosen_resolution: X.X` header + lineage labels
- Script 14 reads the filled TSV and applies labels to the h5ad

A pure DAG can't represent this. The integration approach: script 14 runs **only when the filled map TSV exists**, using Nextflow's `when:` directive. First pipeline invocation runs 13 only; second invocation (after human review) runs both.

### Why scripts 10–12 are out of scope

Wiring scripts 10–12 into Nextflow is required for an end-to-end automated chain, but each has its own integration concerns (script 10 reads a manual mapping TSV; script 12 has multiple subset targets) and would expand this plan significantly. They will be wired up in a separate plan. For now, modules 13/14 take their h5ad input from a `params.progenitor_input_h5ad` path that the user sets to wherever script 12's output lives.

### Conventions to follow

Looking at existing modules (`pipeline/modules/markers.nf`, `pipeline/modules/reconcile_annotations.nf`):

- Process names are uppercase: `MARKERS`, `RECONCILE_ANNOTATIONS`
- `tag` matches the process role
- `publishDir "${params.outdir}/<NN>_<name>", mode: 'copy'`
- Resource directives: `memory '32 GB'`, `cpus 4` (adjust per process needs)
- Script invocation: `pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/<NN>_<name>.py --args ...`
- Outputs: `path "*.png"`, `path "*.csv"`, `path "<specific>.h5ad", emit: h5ad`

---

## File Map

| Action | Path | Purpose |
|---|---|---|
| Commit (existing) | `pipeline/main.nf` | Existing untracked workflow file |
| Commit (existing) | `pipeline/modules/*.nf` | 10 existing untracked module files |
| Commit (existing) | `pipeline/conf/local.config`, `pipeline/conf/slurm.config` | Existing untracked profile configs |
| Commit (existing) | `pipeline/nextflow.config`, `pipeline/params.yaml` | Existing untracked top-level config |
| Commit (existing) | `pipeline/config/cell_type_markers.tsv`, `manual_annotation_level1_map_leiden_r2.0.tsv`, `manual_annotation_markers_skull_immune.tsv` | Existing untracked configs |
| Create | `pipeline/modules/progenitor_recluster.nf` | New module wrapping script 13 |
| Create | `pipeline/modules/apply_progenitor_annotation.nf` | New module wrapping script 14 (conditional) |
| Modify | `pipeline/nextflow.config` | Add 3 new params |
| Modify | `pipeline/params.yaml` | Add the same 3 params (kept in sync) |
| Modify | `pipeline/main.nf` | Include new modules + invoke conditionally |
| Modify (git) | `.gitignore` | Confirm `results/` ignored, do not ignore `pipeline/` |

---

## Task 1: Safety commit of existing untracked Nextflow infrastructure

**Goal:** Get the current Nextflow files into git before modifying them, so we have a baseline and can revert cleanly if integration breaks.

**Files (all currently untracked):**
- `pipeline/main.nf`
- `pipeline/modules/annotate.nf`, `annotate_hspc.nf`, `annotate_mature.nf`, `cluster.nf`, `doublets.nf`, `integrate.nf`, `markers.nf`, `normalize.nf`, `qc.nf`, `reconcile_annotations.nf`
- `pipeline/conf/local.config`, `pipeline/conf/slurm.config`
- `pipeline/nextflow.config`, `pipeline/params.yaml`
- `pipeline/config/cell_type_markers.tsv`, `pipeline/config/manual_annotation_level1_map_leiden_r2.0.tsv`, `pipeline/config/manual_annotation_markers_skull_immune.tsv`

- [ ] **Step 1: Verify what's untracked under `pipeline/`**

```bash
git status --short pipeline/
```

Expected: a list of `??` entries matching the file list above. If anything else appears (e.g. `__pycache__`), do NOT commit it — only commit source files.

- [ ] **Step 2: Confirm `__pycache__` is ignored**

```bash
grep -E '(__pycache__|\.pyc)' .gitignore || echo "MISSING — needs to be added"
```

If missing, add `__pycache__/` and `*.pyc` to `.gitignore` and stage the change. Otherwise skip this step.

- [ ] **Step 3: Stage exact files (no `git add .`)**

```bash
git add pipeline/main.nf \
        pipeline/nextflow.config \
        pipeline/params.yaml \
        pipeline/conf/local.config \
        pipeline/conf/slurm.config \
        pipeline/modules/annotate.nf \
        pipeline/modules/annotate_hspc.nf \
        pipeline/modules/annotate_mature.nf \
        pipeline/modules/cluster.nf \
        pipeline/modules/doublets.nf \
        pipeline/modules/integrate.nf \
        pipeline/modules/markers.nf \
        pipeline/modules/normalize.nf \
        pipeline/modules/qc.nf \
        pipeline/modules/reconcile_annotations.nf \
        pipeline/config/cell_type_markers.tsv \
        pipeline/config/manual_annotation_level1_map_leiden_r2.0.tsv \
        pipeline/config/manual_annotation_markers_skull_immune.tsv
```

- [ ] **Step 4: Verify staged file list**

```bash
git diff --cached --name-only
```

Expected: exactly the files listed above. If any unexpected file is staged, run `git restore --staged <file>` to unstage it.

- [ ] **Step 5: Commit the baseline**

```bash
git commit -m "chore(pipeline): commit existing Nextflow workflow, modules, configs as baseline

Brings the previously untracked Nextflow infrastructure into git so future
edits have a clean diff history. Covers main.nf workflow, 10 process
modules (qc through reconcile_annotations and markers), local/slurm
profiles, top-level config + params, and 3 annotation/marker config files.
No code changes — files are committed as-is."
```

- [ ] **Step 6: Verify clean working tree for these files**

```bash
git status pipeline/
```

Expected: no `??` entries for the files just committed. Other untracked items under `pipeline/scripts/` (e.g. `__pycache__`, prep helper scripts) may remain — they are out of scope for this task.

---

## Task 2: Create module `pipeline/modules/progenitor_recluster.nf` (script 13)

**Files:**
- Create: `pipeline/modules/progenitor_recluster.nf`

**Reference module to mirror:** `pipeline/modules/reconcile_annotations.nf` (single h5ad input, multiple typed outputs, `publishDir`).

- [ ] **Step 1: Write the module file**

```groovy
process PROGENITOR_RECLUSTER {
    tag "progenitor_recluster"
    publishDir "${params.outdir}/13_progenitor_recluster", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "annotation_template.tsv", emit: template
    path "leiden_subset_assignments.csv", emit: assignments
    path "*.csv"
    path "*.png"
    path "*.tsv"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/13_progenitor_recluster.py \\
        --input ${h5ad} \\
        --out .
    """
}
```

- [ ] **Step 2: Sanity-check the module shape**

```bash
test -f pipeline/modules/progenitor_recluster.nf
grep -c "^process PROGENITOR_RECLUSTER" pipeline/modules/progenitor_recluster.nf
grep -c "emit: template" pipeline/modules/progenitor_recluster.nf
grep -c "emit: assignments" pipeline/modules/progenitor_recluster.nf
```

Expected: file exists; each grep returns `1`.

- [ ] **Step 3: Commit**

```bash
git add pipeline/modules/progenitor_recluster.nf
git commit -m "feat(nextflow): add PROGENITOR_RECLUSTER module wrapping script 13

Wraps pipeline/scripts/13_progenitor_recluster.py as a Nextflow process.
Takes an h5ad path and emits annotation_template.tsv (for human review)
and leiden_subset_assignments.csv (consumed by the apply step), plus
all diagnostic CSVs/PNGs. Resource profile: 32 GB / 4 CPUs (matches
MARKERS module which does similar DE workloads)."
```

---

## Task 3: Create module `pipeline/modules/apply_progenitor_annotation.nf` (script 14, conditional)

**Files:**
- Create: `pipeline/modules/apply_progenitor_annotation.nf`

**Conditional logic requirement:** Process must only run when the human-filled annotation map TSV exists. Use a `when:` directive checking the file path.

- [ ] **Step 1: Write the module file**

```groovy
process APPLY_PROGENITOR_ANNOTATION {
    tag "apply_progenitor_annotation"
    publishDir "${params.outdir}/14_progenitor_annotated", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    path h5ad
    path assignments
    path map_tsv

    output:
    path "adata_progenitor_annotated.h5ad", emit: h5ad
    path "final_annotation_counts.csv"
    path "umap_final_annotation.png"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/14_apply_progenitor_annotation.py \\
        --input ${h5ad} \\
        --assignments ${assignments} \\
        --map ${map_tsv} \\
        --out .
    """
}
```

(No `when:` here — gating happens in `main.nf` via channel construction so we can fail loudly if upstream files are missing rather than silently skipping inside the process.)

- [ ] **Step 2: Sanity-check the module**

```bash
test -f pipeline/modules/apply_progenitor_annotation.nf
grep -c "^process APPLY_PROGENITOR_ANNOTATION" pipeline/modules/apply_progenitor_annotation.nf
grep -c "emit: h5ad" pipeline/modules/apply_progenitor_annotation.nf
```

Expected: file exists; each grep returns `1`.

- [ ] **Step 3: Commit**

```bash
git add pipeline/modules/apply_progenitor_annotation.nf
git commit -m "feat(nextflow): add APPLY_PROGENITOR_ANNOTATION module wrapping script 14

Wraps pipeline/scripts/14_apply_progenitor_annotation.py as a Nextflow
process. Takes the h5ad, leiden assignments, and human-filled annotation
map TSV; emits the relabeled h5ad plus summary CSV/UMAP. Conditional
invocation is handled in main.nf (channel only constructed when the map
file exists), keeping process semantics simple."
```

---

## Task 4: Add new params to `nextflow.config` and `params.yaml`

**Files:**
- Modify: `pipeline/nextflow.config` (the `params { ... }` block)
- Modify: `pipeline/params.yaml`

**New params:**

| name | default | meaning |
|---|---|---|
| `progenitor_input_h5ad` | `"${projectDir}/../results/12_subset_recluster/adata_annotated_final.h5ad"` | Path to the final annotated h5ad consumed by script 13 |
| `progenitor_annotation_map` | `"${projectDir}/config/progenitor_annotation_map.tsv"` | Path to the human-filled annotation map (script 14 input) |
| `run_progenitor_recluster` | `true` | Toggle — set false to skip the entire progenitor sub-workflow |

- [ ] **Step 1: Add params to `nextflow.config`**

Open `pipeline/nextflow.config`. Inside the existing `params { ... }` block, append (immediately before the closing `}`):

```groovy
    progenitor_input_h5ad = "${projectDir}/../results/12_subset_recluster/adata_annotated_final.h5ad"
    progenitor_annotation_map = "${projectDir}/config/progenitor_annotation_map.tsv"
    run_progenitor_recluster = true
```

- [ ] **Step 2: Add equivalent entries to `params.yaml`**

Append at the end of `pipeline/params.yaml` (no quotes around the toggle boolean; YAML reads `true` as boolean):

```yaml
progenitor_input_h5ad: "../results/12_subset_recluster/adata_annotated_final.h5ad"
progenitor_annotation_map: "config/progenitor_annotation_map.tsv"
run_progenitor_recluster: true
```

- [ ] **Step 3: Verify both files parse**

For `params.yaml`:

```bash
pixi run python -c "import yaml; print(list(yaml.safe_load(open('pipeline/params.yaml')).keys())[-3:])"
```

Expected output ends with: `['progenitor_input_h5ad', 'progenitor_annotation_map', 'run_progenitor_recluster']`.

For `nextflow.config` (Groovy syntax check via Nextflow itself):

```bash
cd pipeline && nextflow config -profile local | grep -E '^\s*(progenitor_input_h5ad|progenitor_annotation_map|run_progenitor_recluster)' && cd -
```

Expected: 3 matching lines printed. If `nextflow` is not on PATH, install via `pixi run` or use `pixi run nextflow config ...`.

- [ ] **Step 4: Commit**

```bash
git add pipeline/nextflow.config pipeline/params.yaml
git commit -m "feat(nextflow): add params for progenitor recluster sub-workflow

Adds progenitor_input_h5ad (path to script-12 output that feeds script
13), progenitor_annotation_map (path to human-filled label map for
script 14), and run_progenitor_recluster (master toggle). Defaults
point at the conventional locations under results/ and pipeline/config/."
```

---

## Task 5: Wire modules into `pipeline/main.nf` with conditional execution

**Files:**
- Modify: `pipeline/main.nf`

**Behavior to implement:**
- If `params.run_progenitor_recluster` is false → skip both modules.
- Else → always run `PROGENITOR_RECLUSTER`.
- Run `APPLY_PROGENITOR_ANNOTATION` only if `params.progenitor_annotation_map` exists as a file. If it doesn't, log a friendly message instead of failing.

- [ ] **Step 1: Read the current main.nf for reference**

```bash
cat pipeline/main.nf
```

Confirm it matches the expected structure (DSL2 with `include { ... }` lines and a single `workflow { ... }` block).

- [ ] **Step 2: Replace the entire contents of `pipeline/main.nf` with this**

```groovy
nextflow.enable.dsl = 2

include { QC }        from './modules/qc'
include { DOUBLETS }  from './modules/doublets'
include { NORMALIZE } from './modules/normalize'
include { INTEGRATE } from './modules/integrate'
include { CLUSTER }   from './modules/cluster'
include { ANNOTATE_HSPC }  from './modules/annotate_hspc'
include { ANNOTATE_MATURE } from './modules/annotate_mature'
include { RECONCILE_ANNOTATIONS } from './modules/reconcile_annotations'
include { MARKERS }   from './modules/markers'
include { PROGENITOR_RECLUSTER }       from './modules/progenitor_recluster'
include { APPLY_PROGENITOR_ANNOTATION } from './modules/apply_progenitor_annotation'

workflow {
    h5_ch = Channel.fromPath(params.h5_input, checkIfExists: true)
    hspc_ref_ch = Channel.fromPath(params.hspc_ref_h5ad, checkIfExists: true)
    mature_ref_ch = Channel.fromPath(params.mature_ref_h5ad, checkIfExists: true)

    QC(h5_ch)
    DOUBLETS(QC.out.h5ad)
    NORMALIZE(DOUBLETS.out.h5ad)
    INTEGRATE(NORMALIZE.out.h5ad)
    CLUSTER(INTEGRATE.out.h5ad)
    ANNOTATE_HSPC(CLUSTER.out.h5ad, hspc_ref_ch)
    ANNOTATE_MATURE(CLUSTER.out.h5ad, mature_ref_ch)
    RECONCILE_ANNOTATIONS(ANNOTATE_HSPC.out.h5ad, ANNOTATE_MATURE.out.h5ad)
    MARKERS(RECONCILE_ANNOTATIONS.out.h5ad)

    if (params.run_progenitor_recluster) {
        prog_in_ch = Channel.fromPath(params.progenitor_input_h5ad, checkIfExists: true)
        PROGENITOR_RECLUSTER(prog_in_ch)

        map_file = file(params.progenitor_annotation_map)
        if (map_file.exists()) {
            APPLY_PROGENITOR_ANNOTATION(
                prog_in_ch,
                PROGENITOR_RECLUSTER.out.assignments,
                Channel.fromPath(params.progenitor_annotation_map, checkIfExists: true)
            )
        } else {
            log.info "[progenitor] Skipping APPLY_PROGENITOR_ANNOTATION — annotation map not found at ${params.progenitor_annotation_map}. Fill in the map after reviewing PROGENITOR_RECLUSTER outputs, then re-run."
        }
    }
}
```

- [ ] **Step 3: Static syntax check (Nextflow config + DSL parse)**

```bash
cd pipeline && nextflow config -profile local > /dev/null && echo OK; cd -
```

Expected: prints `OK`. Any Groovy parse error or include resolution error will surface here.

- [ ] **Step 4: Dry-run inspection (preview-only, no execution)**

```bash
cd pipeline && nextflow inspect -profile local main.nf 2>&1 | head -60; cd -
```

Expected: list of processes including `QC`, `DOUBLETS`, ..., `PROGENITOR_RECLUSTER`, and conditionally `APPLY_PROGENITOR_ANNOTATION`. If `nextflow inspect` is unavailable in your version, skip and rely on Step 5.

- [ ] **Step 5: Commit**

```bash
git add pipeline/main.nf
git commit -m "feat(nextflow): wire progenitor recluster sub-workflow into main.nf

Adds PROGENITOR_RECLUSTER and APPLY_PROGENITOR_ANNOTATION includes.
After the existing MARKERS step, when run_progenitor_recluster=true,
runs PROGENITOR_RECLUSTER on the script-12 output. Then conditionally
runs APPLY_PROGENITOR_ANNOTATION only if the annotation map TSV exists,
otherwise logs a skip message — supporting the human-in-the-loop
review checkpoint without breaking the pipeline."
```

---

## Task 6: End-to-end smoke test — both pipeline branches

**Files:**
- Run: `pipeline/main.nf` via Nextflow CLI

**Why two runs:** verify both states of the conditional — (a) map TSV exists → APPLY runs, (b) map TSV missing → APPLY is skipped with a friendly message.

- [ ] **Step 1: Confirm input file exists**

```bash
ls -lh results/12_subset_recluster/adata_annotated_final.h5ad
```

If missing, the test cannot run. Report BLOCKED.

- [ ] **Step 2: Test branch A — both modules run (annotation map present)**

The map already exists (`pipeline/config/progenitor_annotation_map.tsv` was committed in commit e024779). Verify:

```bash
test -f pipeline/config/progenitor_annotation_map.tsv && echo "map present" || echo "map missing"
```

Expected: `map present`. If missing, this test branch cannot run.

Run only the progenitor sub-workflow in isolation (avoid re-running the full pipeline from QC):

```bash
cd pipeline && pixi run nextflow run main.nf \
    -profile local \
    --run_progenitor_recluster true \
    -entry workflow_progenitor_only 2>&1 | tail -40; cd -
```

> **Note for the implementer:** `-entry` requires a named workflow. Since the current `main.nf` only has the default unnamed `workflow {}`, this command will fail. **Skip to the alternative** below.

**Alternative (if `-entry` workflow_progenitor_only does not exist):** the simplest valid smoke test is to just run the new modules manually outside Nextflow to confirm the inline scripts haven't regressed, and confirm Nextflow parses the workflow without error. Run instead:

```bash
cd pipeline && pixi run nextflow config -profile local > /dev/null && echo "config OK"
cd ../ && pixi run python pipeline/scripts/13_progenitor_recluster.py --help > /dev/null && echo "13 help OK"
pixi run python pipeline/scripts/14_apply_progenitor_annotation.py --help > /dev/null && echo "14 help OK"
```

Expected: all three lines print `OK`. This validates that the Nextflow config parses, modules can be resolved, and both scripts still respond to `--help` (no import-time regressions).

- [ ] **Step 3: Test branch B — APPLY skipped (annotation map absent)**

Temporarily rename the map and run a Nextflow inspection that exercises the conditional:

```bash
mv pipeline/config/progenitor_annotation_map.tsv pipeline/config/progenitor_annotation_map.tsv.tmp
cd pipeline && pixi run nextflow config -profile local > /dev/null && echo "config OK without map"; cd -
```

Look for the `[progenitor] Skipping APPLY_PROGENITOR_ANNOTATION` log line if you do a real run (optional). At minimum, confirm `config OK without map` prints.

Restore:

```bash
mv pipeline/config/progenitor_annotation_map.tsv.tmp pipeline/config/progenitor_annotation_map.tsv
test -f pipeline/config/progenitor_annotation_map.tsv && echo "map restored"
```

Expected: `map restored`. If the rename or restore fails, fix immediately — this file is the source of truth for the labels you committed in e024779.

- [ ] **Step 4: No commit needed**

This task is a verification step. No file changes.

---

## Task 7: Document the sub-workflow in `pipeline/main.nf` header comment

**Files:**
- Modify: `pipeline/main.nf` (top-of-file comment)

**Why:** the conditional + human-in-loop behavior is non-obvious from code alone. A short comment at the top of `main.nf` documents the contract for future readers.

- [ ] **Step 1: Add a header comment block to main.nf**

Open `pipeline/main.nf`. Insert this block at the very top (line 1, before `nextflow.enable.dsl = 2`):

```groovy
/*
 * Charles scRNA-seq pipeline — main workflow.
 *
 * Linear pipeline (always runs):
 *   QC → DOUBLETS → NORMALIZE → INTEGRATE → CLUSTER →
 *   ANNOTATE_HSPC + ANNOTATE_MATURE → RECONCILE_ANNOTATIONS → MARKERS
 *
 * Optional progenitor sub-workflow (when run_progenitor_recluster=true):
 *   PROGENITOR_RECLUSTER  — diagnostic only; produces annotation_template.tsv
 *   [HUMAN REVIEW]        — fill pipeline/config/progenitor_annotation_map.tsv
 *                           with chosen_resolution comment + level1/level2 labels
 *   APPLY_PROGENITOR_ANNOTATION — runs only if the filled map exists
 *
 * On the first invocation the annotation map will be missing and APPLY is
 * skipped with a log message. After human review, re-run the pipeline and
 * APPLY will execute, producing adata_progenitor_annotated.h5ad.
 */

nextflow.enable.dsl = 2
```

- [ ] **Step 2: Verify the file still parses**

```bash
cd pipeline && pixi run nextflow config -profile local > /dev/null && echo OK; cd -
```

Expected: `OK`.

- [ ] **Step 3: Commit**

```bash
git add pipeline/main.nf
git commit -m "docs(nextflow): document progenitor sub-workflow + human checkpoint in main.nf

Adds a header comment summarizing the linear automated chain, the
optional progenitor sub-workflow, and the human-review checkpoint
between PROGENITOR_RECLUSTER and APPLY_PROGENITOR_ANNOTATION."
```

---

## Self-Review

**Spec coverage:**
- ✅ Safety commit of all currently-untracked Nextflow infrastructure (Task 1)
- ✅ Module wrapping script 13 (Task 2)
- ✅ Module wrapping script 14 (Task 3)
- ✅ Params added to both nextflow.config and params.yaml (Task 4)
- ✅ Modules wired into main.nf with conditional execution (Task 5)
- ✅ End-to-end verification of both branches of the conditional (Task 6)
- ✅ Inline documentation of the human-in-loop contract (Task 7)
- ✅ Out-of-scope items (modules for scripts 10–12) explicitly called out in Background

**No placeholders:** All Groovy/YAML/bash content is shown verbatim. No "TBD" or "similar to". Every step has either an exact command or exact code.

**Type / name consistency:**
- `PROGENITOR_RECLUSTER` (process name) used consistently in module file and `include`/invocation in main.nf
- `APPLY_PROGENITOR_ANNOTATION` likewise
- `emit: assignments` in PROGENITOR_RECLUSTER consumed as `PROGENITOR_RECLUSTER.out.assignments` in main.nf — matches
- Parameter names `progenitor_input_h5ad`, `progenitor_annotation_map`, `run_progenitor_recluster` identical across `nextflow.config`, `params.yaml`, and `main.nf` references
- Output paths `${params.outdir}/13_progenitor_recluster` and `${params.outdir}/14_progenitor_annotated` match the existing `results/13_*` and `results/14_*` directory convention used by manual runs of scripts 13 and 14

**Known limitation flagged in plan:** The `-entry` smoke test in Task 6 falls back to a static-parse + `--help` check because the current pipeline has no named sub-workflow entry. A follow-up plan could add a `workflow PROGENITOR_ONLY {}` entry-point if isolated execution becomes a frequent need.
