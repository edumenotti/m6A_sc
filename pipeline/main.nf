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
