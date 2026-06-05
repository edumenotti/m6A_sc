/*
 * Usage:
 *   Local workstation (pixi installed):
 *     nextflow run main.nf -params-file params.yaml --h5_input /path/to/filtered_feature_bc_matrix.h5
 *   Yale Bouchet (Apptainer + SLURM):
 *     nextflow run main.nf -profile bouchet -params-file params.yaml \
 *       --h5_input /path/to/filtered_feature_bc_matrix.h5
 *   Other SLURM cluster: copy the `bouchet` profile in nextflow.config, set your
 *     CPU/GPU partition names, and run with that profile.
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
include { FREEZE_MANUAL_ANNOTATION }    from './modules/freeze_manual_annotation'
include { MACROPHAGE_STATES }  from './modules/macrophage_states'
include { COMPOSITION_SCCODA } from './modules/composition'
include { PSEUDOBULK_DEG }     from './modules/pseudobulk_deg'
include { PATHWAY_ACTIVITY }   from './modules/pathway_activity'
include { EXPORT_FOR_HEMASCRIBE } from './modules/export_for_hemascribe'
include { HEMASCRIBE }            from './modules/hemascribe'
include { MERGE_HEMASCRIBE }      from './modules/merge_hemascribe'
include { FINALIZE_ANNOTATION }   from './modules/finalize_annotation'
include { CLEAN_COLUMNS }         from './modules/clean_columns'

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
     * Manual annotation → adata_progenitor_annotated.h5ad (feeds the HemaScribe chain).
     *
     * Default (rederive_manual_annotation=false): FREEZE_MANUAL_ANNOTATION transfers
     * the canonical per-cell manual_level1/level2 labels onto the fresh object BY
     * BARCODE. This is reproducible on any hardware. Cluster-ID-based maps are not:
     * scVI integration + Leiden drift across machines (different cluster count/IDs),
     * so a `cluster_id -> label` map cannot be reapplied to a re-clustered object.
     *
     * rederive_manual_annotation=true: re-run the original cluster-based chain
     * (MANUAL_MARKER_REVIEW → APPLY_MANUAL_ANNOTATION → SUBSET_RECLUSTER →
     * FINAL_FIGURES → PROGENITOR_RECLUSTER → APPLY_PROGENITOR_ANNOTATION). Only
     * reproducible on the exact machine/library versions the maps were built on.
     */
    prog_annotated_ch = Channel.empty()
    if (params.rederive_manual_annotation) {
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
                prog_annotated_ch = APPLY_PROGENITOR_ANNOTATION.out.h5ad.first()
            } else {
                log.info "[progenitor] Skipping APPLY_PROGENITOR_ANNOTATION — annotation map not found at ${params.progenitor_annotation_map}. Fill in the map after reviewing PROGENITOR_RECLUSTER outputs, then re-run."
            }
        }
    } else {
        frozen_ann_ch = Channel.fromPath(params.frozen_per_cell_labels, checkIfExists: true)
        FREEZE_MANUAL_ANNOTATION(annotated_h5ad_ch, frozen_ann_ch)
        prog_annotated_ch = FREEZE_MANUAL_ANNOTATION.out.h5ad.first()
    }

    /*
     * Optional downstream analysis block (run_downstream_analysis=true):
     *   Reproduces the canonical annotation and the analyses reported to the lab.
     *
     *   Canonical-annotation chain (scripts 21a/b/c → 21z → 22):
     *     EXPORT_FOR_HEMASCRIBE → HEMASCRIBE → MERGE_HEMASCRIBE
     *       → FINALIZE_ANNOTATION (adds cell_type) → CLEAN_COLUMNS
     *     produces results/14_progenitor_annotated/adata_hemascribe.h5ad
     *     (canonical `cell_type`, 32 types).
     *
     *   Analyses on the canonical `cell_type` level (params.downstream_levels):
     *     - MACROPHAGE_STATES: per-cell pathway scoring (decoupler ULM) on the manual object.
     *     - COMPOSITION_SCCODA: Bayesian compositional analysis (scCODA, CPU-only).
     *     - PSEUDOBULK_DEG: per-celltype factorial DEG with PyDESeq2 (donor as replicate).
     *     - PATHWAY_ACTIVITY: pathway scoring per celltype × contrast from DEG stats.
     *
     *   The old CellChat/NicheNet block was deprecated 2026-05-12 (composition
     *   confound + n=2 → exploratory only); see results/_archive_exploratory/.
     */
    if (params.run_downstream_analysis) {
        if (params.rederive_manual_annotation) {
            if (!params.run_progenitor_recluster || !file(params.progenitor_annotation_map).exists()) {
                error "rederive_manual_annotation=true + run_downstream_analysis=true requires run_progenitor_recluster=true and a filled ${params.progenitor_annotation_map}, so APPLY_PROGENITOR_ANNOTATION can produce adata_progenitor_annotated.h5ad first."
            }
        } else if (!file(params.frozen_per_cell_labels).exists()) {
            error "run_downstream_analysis=true requires the frozen per-cell labels at ${params.frozen_per_cell_labels} (or set rederive_manual_annotation=true to re-derive them from clustering)."
        }

        // Per-cell macrophage states run on the manual (pre-HemaScribe) object.
        MACROPHAGE_STATES(prog_annotated_ch)

        // Canonical annotation chain → adata_hemascribe.h5ad (cell_type).
        EXPORT_FOR_HEMASCRIBE(prog_annotated_ch)
        HEMASCRIBE(EXPORT_FOR_HEMASCRIBE.out.export_dir)
        MERGE_HEMASCRIBE(prog_annotated_ch.combine(HEMASCRIBE.out.labels))
        FINALIZE_ANNOTATION(MERGE_HEMASCRIBE.out.h5ad)
        CLEAN_COLUMNS(FINALIZE_ANNOTATION.out.h5ad)
        canonical_ch = CLEAN_COLUMNS.out.h5ad.first()

        // Downstream analyses on the canonical cell_type level(s).
        levels_ch = Channel.fromList(params.downstream_levels.tokenize(','))
        COMPOSITION_SCCODA(canonical_ch.combine(levels_ch))
        PSEUDOBULK_DEG(canonical_ch.combine(levels_ch))
        PATHWAY_ACTIVITY(PSEUDOBULK_DEG.out.deg_dir)
    }
}
