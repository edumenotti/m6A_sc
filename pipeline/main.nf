/*
 * Usage:
 *   Local (pixi installed):  nextflow run main.nf -profile local  -params-file params.yaml
 *   Yale HPC (Apptainer):    nextflow run main.nf -profile apptainer,slurm -params-file params.yaml
 *
 * Build container (once, from project root):
 *   sudo apptainer build charles-scrna.sif charles-scrna.def        # local with sudo
 *   apptainer build --fakeroot charles-scrna.sif charles-scrna.def  # HPC (no sudo needed)
 *   scp charles-scrna.sif <netid>@grace.hpc.yale.edu:<project-dir>/
 */

/*
 * Charles scRNA-seq pipeline — main workflow.
 *
 * Always runs:
 *   QC → DOUBLETS → NORMALIZE → INTEGRATE → CLUSTER → MARKERS →
 *   MANUAL_MARKER_REVIEW → APPLY_MANUAL_ANNOTATION → SUBSET_RECLUSTER → FINAL_FIGURES
 *
 * Optional reference-based annotation (run_label_transfer=true):
 *   CLUSTER → ANNOTATE_HSPC + ANNOTATE_MATURE → RECONCILE_ANNOTATIONS
 *   Produces popv_prediction / scanvi_label / final_cell_type as context for
 *   MANUAL_MARKER_REVIEW. The manual_level1/level2 labels in the final h5ad
 *   come from the human-filled annotation map (script 10), not from popV, so
 *   this step is informational. Default off to save GPU time.
 *
 * Optional progenitor sub-workflow (run_progenitor_recluster=true):
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

include { QC }        from './modules/qc'
include { DOUBLETS }  from './modules/doublets'
include { NORMALIZE } from './modules/normalize'
include { INTEGRATE } from './modules/integrate'
include { CLUSTER }   from './modules/cluster'
include { ANNOTATE_HSPC }  from './modules/annotate_hspc'
include { ANNOTATE_MATURE } from './modules/annotate_mature'
include { RECONCILE_ANNOTATIONS } from './modules/reconcile_annotations'
include { MARKERS }   from './modules/markers'
include { MANUAL_MARKER_REVIEW }    from './modules/manual_marker_review'
include { APPLY_MANUAL_ANNOTATION } from './modules/apply_manual_annotation'
include { SUBSET_RECLUSTER }        from './modules/subset_recluster'
include { FINAL_FIGURES }           from './modules/final_figures'
include { PROGENITOR_RECLUSTER }       from './modules/progenitor_recluster'
include { APPLY_PROGENITOR_ANNOTATION } from './modules/apply_progenitor_annotation'
include { MACROPHAGE_STATES }  from './modules/macrophage_states'
include { COMPOSITION_SCCODA } from './modules/composition'
include { PSEUDOBULK_DEG }     from './modules/pseudobulk_deg'
include { PATHWAY_ACTIVITY }   from './modules/pathway_activity'

workflow {
    h5_ch = Channel.fromPath(params.h5_input, checkIfExists: true)

    QC(h5_ch)
    DOUBLETS(QC.out.h5ad)
    NORMALIZE(DOUBLETS.out.h5ad)
    INTEGRATE(NORMALIZE.out.h5ad)
    CLUSTER(INTEGRATE.out.h5ad)

    /*
     * Optional reference-based annotation (run_label_transfer):
     * If enabled, runs popV against the Nestorowa (HSPC) and Kucinski (mature)
     * references and reconciles them. The reconciled h5ad carries popV/scanvi
     * columns + a 'final_cell_type' that the manual review (script 09) uses as
     * context. Manual labels (script 10) override these regardless.
     * If disabled, the cluster output is fed directly into MARKERS / manual review;
     * scripts 09 and 10 already guard their popV/final_cell_type lookups with
     * `if col in obs:`, so the missing columns are handled gracefully.
     */
    if (params.run_label_transfer) {
        hspc_ref_ch = Channel.fromPath(params.hspc_ref_h5ad, checkIfExists: true)
        mature_ref_ch = Channel.fromPath(params.mature_ref_h5ad, checkIfExists: true)
        ANNOTATE_HSPC(CLUSTER.out.h5ad, hspc_ref_ch)
        ANNOTATE_MATURE(CLUSTER.out.h5ad, mature_ref_ch)
        RECONCILE_ANNOTATIONS(ANNOTATE_HSPC.out.h5ad, ANNOTATE_MATURE.out.h5ad)
        annotated_h5ad_ch = RECONCILE_ANNOTATIONS.out.h5ad
    } else {
        annotated_h5ad_ch = CLUSTER.out.h5ad
    }

    MARKERS(annotated_h5ad_ch)

    /*
     * Post-cluster annotation sub-workflow (scripts 09–12):
     *   MANUAL_MARKER_REVIEW  — produces diagnostics for human review
     *   [HUMAN: fill pipeline/config/manual_annotation_level1_map_leiden_r2.0.tsv]
     *   APPLY_MANUAL_ANNOTATION — applies filled map; map is committed so this always runs
     *   SUBSET_RECLUSTER      — B cell / erythroid level2 sub-clustering
     *   FINAL_FIGURES         — curated summary figures for the annotated dataset
     */
    annotation_map_ch = Channel.fromPath(params.manual_annotation_map, checkIfExists: true)

    MANUAL_MARKER_REVIEW(annotated_h5ad_ch)
    APPLY_MANUAL_ANNOTATION(annotated_h5ad_ch, annotation_map_ch)
    SUBSET_RECLUSTER(APPLY_MANUAL_ANNOTATION.out.h5ad)
    FINAL_FIGURES(SUBSET_RECLUSTER.out.h5ad)

    if (params.run_progenitor_recluster) {
        PROGENITOR_RECLUSTER(SUBSET_RECLUSTER.out.h5ad)

        map_file = file(params.progenitor_annotation_map)
        if (map_file.exists()) {
            APPLY_PROGENITOR_ANNOTATION(
                SUBSET_RECLUSTER.out.h5ad,
                PROGENITOR_RECLUSTER.out.assignments,
                Channel.fromPath(params.progenitor_annotation_map, checkIfExists: true)
            )
        } else {
            log.info "[progenitor] Skipping APPLY_PROGENITOR_ANNOTATION — annotation map not found at ${params.progenitor_annotation_map}. Fill in the map after reviewing PROGENITOR_RECLUSTER outputs, then re-run."
        }
    }

    /*
     * Optional downstream analysis block (run_downstream_analysis=true):
     *   Requires APPLY_PROGENITOR_ANNOTATION to have produced the annotated h5ad
     *   with genotype column (replicate 1=WT, 2=Mutant).
     *
     *   - MACROPHAGE_STATES: per-cell pathway scoring (decoupler ULM, not affected
     *     by composition confound — see plan 2026-05-12).
     *   - COMPOSITION_SCCODA: Bayesian compositional analysis (scCODA, CPU-only).
     *   - PSEUDOBULK_DEG: per-celltype factorial DEG with PyDESeq2 (donor as replicate).
     *   - PATHWAY_ACTIVITY: pathway scoring per celltype × contrast from DEG stats.
     *
     *   The old CellChat/NicheNet block was deprecated 2026-05-12 (composition
     *   confound + n=2 → exploratory only). Outputs were moved to
     *   results/_archive_exploratory/; see plan 2026-05-12-composition-and-pseudobulk-analysis.md.
     */
    if (params.run_downstream_analysis) {
        downstream_map_file = file(params.progenitor_annotation_map)
        if (!downstream_map_file.exists()) {
            error "run_downstream_analysis=true requires progenitor_annotation_map to exist. Run the progenitor sub-workflow first."
        }
        // Use file() (value channel) so the same path can feed multiple processes
        // without queue-channel exhaustion. Recreate Channel.of for each .map()
        // because operators consume queue channels once.
        prog_annotated = file(
            "${params.outdir}/14_progenitor_annotated/adata_progenitor_annotated.h5ad",
            checkIfExists: true
        )

        MACROPHAGE_STATES(Channel.of(prog_annotated))

        COMPOSITION_SCCODA(
            Channel.of('manual_level1', 'manual_level2').map { lvl -> tuple(prog_annotated, lvl) }
        )

        PSEUDOBULK_DEG(
            Channel.of('manual_level1', 'manual_level2').map { lvl -> tuple(prog_annotated, lvl) }
        )

        PATHWAY_ACTIVITY(PSEUDOBULK_DEG.out.deg_dir)
    }
}
