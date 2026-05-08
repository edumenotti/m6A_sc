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
}
