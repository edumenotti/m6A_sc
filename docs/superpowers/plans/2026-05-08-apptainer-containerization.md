# Apptainer Containerization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Package the full pipeline environment into an Apptainer `.sif` container so the pipeline can run reproducibly on Yale HPC without requiring pixi to be installed on the cluster.

**Architecture:** Build an Apptainer image from a `.def` file that installs pixi, copies `pixi.toml`+`pixi.lock` into the container, and runs `pixi install --locked` to bake the exact locked environment. The Nextflow modules currently call `pixi run -m ${projectDir}/../pixi.toml python ...`; inside the container pixi is installed and the project manifest lives at `/project/pixi.toml`, so a new `params.pixi_manifest` parameter lets the apptainer profile override the manifest path without changing individual module scripts. A new `apptainer` profile in `nextflow.config` enables the container and sets the correct manifest path.

**Tech Stack:** Apptainer ≥1.0, Nextflow DSL2, pixi, Ubuntu 22.04

---

## Key Design Decisions

- **Why pixi inside the container, not a plain conda env?** pixi.lock is the authoritative pin of all 300+ packages. Re-creating the env via conda's own solver risks drift. Using `pixi install --locked` inside the container guarantees bit-for-bit reproduction.
- **Why `params.pixi_manifest`?** All 14 module scripts use `pixi run -m <path> python`. A single parameter override in the apptainer profile is less fragile than patching 14 files.
- **GPU support:** Apptainer passes through host NVIDIA drivers with `--nv`. The container only needs CUDA runtime libraries, which pixi installs from conda-forge (already in `pixi.lock` via `cuda = "12.8"` system-requirement). The Yale HPC job submission in `slurm.config` already requests `--gres=gpu:1`.
- **`.sif` file location:** Built once locally, then copied to HPC. It is large (several GB) and is NOT committed to git. Add it to `.gitignore`.

---

## File Map

| File | Action | Purpose |
|------|--------|---------|
| `charles-scrna.def` | Create | Apptainer definition file — builds the container |
| `pipeline/conf/apptainer.config` | Create | Nextflow apptainer profile config |
| `pipeline/nextflow.config` | Modify | Add `params.pixi_manifest` default + `apptainer` profile include |
| `pipeline/params.yaml` | Modify | Add `pixi_manifest` default |
| `pipeline/modules/*.nf` (all 14) | Modify | Replace hardcoded `${projectDir}/../pixi.toml` with `${params.pixi_manifest}` |
| `.gitignore` | Modify | Add `*.sif` to ignore built container images |

---

### Task 1: Write the Apptainer definition file

**Files:**
- Create: `charles-scrna.def`

- [ ] **Step 1: Check pixi and Apptainer availability**

Run:
```bash
which apptainer && apptainer --version
which pixi && pixi --version
```
Note the versions. If Apptainer is not installed locally, the build step (Task 5) must be done on HPC or a machine with Apptainer. The `.def` file can still be written and committed.

- [ ] **Step 2: Create the definition file**

Create `charles-scrna.def` with this content:
```
Bootstrap: docker
From: ubuntu:22.04

%labels
    Maintainer Eduardo Menotti
    Project charles-scrna
    Description scRNA-seq pipeline environment (pixi-locked)

%files
    pixi.toml /project/pixi.toml
    pixi.lock /project/pixi.lock

%post
    export DEBIAN_FRONTEND=noninteractive
    apt-get update && apt-get install -y \
        curl \
        ca-certificates \
        build-essential \
        libssl-dev \
        libffi-dev \
        git \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

    # Install pixi to /opt/pixi
    export PIXI_HOME=/opt/pixi
    curl -fsSL https://pixi.sh/install.sh | bash
    export PATH="/opt/pixi/bin:$PATH"

    # Install the locked environment into the container
    # This reads pixi.lock exactly — no solver, no drift
    pixi install --locked --manifest-path /project/pixi.toml

    # Smoke-test: Python and R must be importable
    /project/.pixi/envs/default/bin/python -c "import scanpy; print('scanpy', scanpy.__version__)"
    /project/.pixi/envs/default/bin/Rscript -e "cat('R OK\n')"

%environment
    export PIXI_HOME=/opt/pixi
    export PATH="/opt/pixi/bin:/project/.pixi/envs/default/bin:$PATH"
    export PIXI_PROJECT_DIR=/project
    export R_LIBS_USER=/project/.pixi/envs/default/lib/R/library

%runscript
    exec "$@"
```

- [ ] **Step 3: Verify the file was created**

Run:
```bash
head -5 charles-scrna.def
```
Expected: starts with `Bootstrap: docker`

- [ ] **Step 4: Commit the definition file**

```bash
git add charles-scrna.def
git commit -m "feat: add Apptainer definition file for charles-scrna pipeline"
```

---

### Task 2: Add *.sif to .gitignore

**Files:**
- Modify: `.gitignore`

- [ ] **Step 1: Read current .gitignore**

Run:
```bash
cat .gitignore
```

- [ ] **Step 2: Add container image patterns**

Append to `.gitignore`:
```
# Apptainer/Singularity container images (large binary, not versioned in git)
*.sif
*.simg
```

- [ ] **Step 3: Commit**

```bash
git add .gitignore
git commit -m "chore: ignore Apptainer container image files (*.sif, *.simg)"
```

---

### Task 3: Add params.pixi_manifest to nextflow.config and all modules

**Files:**
- Modify: `pipeline/nextflow.config`
- Modify: `pipeline/params.yaml`
- Modify: all 14 files in `pipeline/modules/*.nf`

- [ ] **Step 1: Verify all modules use the hardcoded pixi path**

Run:
```bash
grep -l 'projectDir.*pixi.toml' pipeline/modules/*.nf
```
Expected: all 14 `.nf` files are listed.

- [ ] **Step 2: Add pixi_manifest param to nextflow.config**

In `pipeline/nextflow.config`, inside `params { }`, add after the existing params:
```
    pixi_manifest = "${projectDir}/../pixi.toml"
```

- [ ] **Step 3: Add pixi_manifest param to params.yaml**

In `pipeline/params.yaml`, add:
```yaml
pixi_manifest: "pipeline/../pixi.toml"
```

- [ ] **Step 4: Replace hardcoded pixi.toml path in all 14 modules**

Run this sed command to replace in all module files at once:
```bash
sed -i 's|${projectDir}/../pixi.toml|${params.pixi_manifest}|g' pipeline/modules/*.nf
```

- [ ] **Step 5: Verify the replacement worked in every module**

Run:
```bash
grep -r 'projectDir.*pixi.toml' pipeline/modules/
```
Expected: no output (all replaced).

Run:
```bash
grep -c 'params.pixi_manifest' pipeline/modules/*.nf | grep -v ':0'
```
Expected: all 14 modules show at least 1 match.

- [ ] **Step 6: Commit**

```bash
git add pipeline/nextflow.config pipeline/params.yaml pipeline/modules/*.nf
git commit -m "feat: parameterize pixi manifest path via params.pixi_manifest for container support"
```

---

### Task 4: Create the Nextflow apptainer profile

**Files:**
- Create: `pipeline/conf/apptainer.config`
- Modify: `pipeline/nextflow.config`

- [ ] **Step 1: Create pipeline/conf/apptainer.config**

Create `pipeline/conf/apptainer.config` with this content:
```groovy
apptainer {
    enabled    = true
    autoMounts = true
}

params {
    // Inside the container, pixi and the locked env are at /project/
    pixi_manifest   = '/project/pixi.toml'
    container_image = "${projectDir}/../charles-scrna.sif"
}

process {
    container = params.container_image

    // GPU processes still need --nv to pass through host NVIDIA drivers
    withName: 'INTEGRATE|ANNOTATE_HSPC|ANNOTATE_MATURE' {
        containerOptions = '--nv'
    }
}
```

- [ ] **Step 2: Register the apptainer profile in nextflow.config**

In `pipeline/nextflow.config`, inside the `profiles { }` block, add:
```
    apptainer {
        includeConfig 'conf/apptainer.config'
    }
```

- [ ] **Step 3: Verify the profiles block now has three entries**

Run:
```bash
grep -A 1 'profiles {' pipeline/nextflow.config
grep 'includeConfig' pipeline/nextflow.config
```
Expected: three `includeConfig` lines (local, slurm, apptainer).

- [ ] **Step 4: Commit**

```bash
git add pipeline/conf/apptainer.config pipeline/nextflow.config
git commit -m "feat: add Nextflow apptainer profile with container path and pixi manifest override"
```

---

### Task 5: Build the container image

> **Note:** This step requires Apptainer installed and root or fakeroot privileges. On Yale HPC, use `apptainer build --fakeroot`. If building locally, `sudo apptainer build` works. The `.sif` file is large (3–6 GB) and is not committed to git.

**Files:**
- Output: `charles-scrna.sif` (not committed)

- [ ] **Step 1: Confirm you are in the project root**

Run:
```bash
ls charles-scrna.def pixi.toml pixi.lock
```
Expected: all three files present.

- [ ] **Step 2: Build the image (locally with sudo)**

Run (takes 15–40 minutes due to pixi install):
```bash
sudo apptainer build charles-scrna.sif charles-scrna.def 2>&1 | tee build.log
```

Or on Yale HPC (fakeroot, no sudo):
```bash
apptainer build --fakeroot charles-scrna.sif charles-scrna.def 2>&1 | tee build.log
```

Expected final lines in build.log:
```
scanpy 1.9.x
R OK
INFO:    Adding environment to container
INFO:    Creating SIF file...
INFO:    Build complete: charles-scrna.sif
```

- [ ] **Step 3: Verify scanpy and R work inside the container**

Run:
```bash
apptainer exec charles-scrna.sif python -c "import scanpy; print(scanpy.__version__)"
apptainer exec charles-scrna.sif Rscript -e "library(Seurat); cat('Seurat OK\n')"
```
Expected: version strings printed without errors.

- [ ] **Step 4: Verify pixi run works inside the container**

Run:
```bash
apptainer exec charles-scrna.sif pixi run -m /project/pixi.toml python -c "import scvi; print('scvi OK')"
```
Expected: `scvi OK`

- [ ] **Step 5: Clean up build log (not needed in git)**

Run:
```bash
rm -f build.log
```

---

### Task 6: Test the full pipeline with the apptainer profile (dry-run)

**Files:** None modified.

- [ ] **Step 1: Run a Nextflow dry-run with the apptainer profile**

Run from the project root:
```bash
cd pipeline && nextflow run main.nf -profile apptainer -params-file params.yaml -stub 2>&1 | tail -30
```

If `-stub` is not supported on this Nextflow version, use `-preview`:
```bash
nextflow run main.nf -profile apptainer -params-file params.yaml -preview 2>&1 | tail -30
```

Expected: all processes listed, no `WARN` about missing container, no `ERROR`.

- [ ] **Step 2: Verify the container path is resolved correctly in the dry-run output**

Run:
```bash
cd pipeline && nextflow run main.nf -profile apptainer -params-file params.yaml -stub 2>&1 | grep -i "container\|sif"
```
Expected: references to `charles-scrna.sif`.

- [ ] **Step 3: Document how to run on Yale HPC in a README comment**

Add a comment block to `pipeline/main.nf` at the very top (before `nextflow.enable.dsl = 2`):

```nextflow
/*
 * Usage:
 *   Local (pixi installed):  nextflow run main.nf -profile local  -params-file params.yaml
 *   Yale HPC (Apptainer):    nextflow run main.nf -profile apptainer,slurm -params-file params.yaml
 *
 * Build container (once, from project root):
 *   sudo apptainer build charles-scrna.sif charles-scrna.def        # local
 *   apptainer build --fakeroot charles-scrna.sif charles-scrna.def  # HPC (no sudo)
 *   scp charles-scrna.sif <netid>@grace.hpc.yale.edu:<project-dir>/
 */
```

- [ ] **Step 4: Commit**

```bash
git add pipeline/main.nf
git commit -m "docs: add HPC usage instructions and apptainer build notes to main.nf header"
```

---

## Self-Review

**Spec coverage:**
- ✅ Apptainer `.def` file builds environment from `pixi.lock` (locked, reproducible)
- ✅ GPU pass-through handled via `--nv` in apptainer.config for GPU processes
- ✅ No module scripts need changing beyond the `pixi_manifest` param substitution
- ✅ Local and slurm profiles unaffected — `params.pixi_manifest` defaults to project-local pixi.toml
- ✅ Container image excluded from git via `.gitignore`
- ✅ Usage instructions added to main.nf header
- ✅ Yale HPC specifics: fakeroot build documented, `apptainer,slurm` combined profile documented

**Placeholder scan:** No TBD. Build command, expected output, and verification steps are all explicit. ✅

**Risk:** `pixi install --locked` inside the container requires network access during build. The `--locked` flag prevents solver drift but still downloads packages from conda-forge/bioconda. On air-gapped HPC systems, the `.sif` must be built externally and copied over (documented in Task 5). ✅
