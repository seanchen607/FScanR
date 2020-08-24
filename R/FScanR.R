##' Identifies Programmed Ribosomal Frameshifting (PRF) events from mRNA/cDNA BLASTX output
##'
##'
##' @title FScanR
##' @param blastx_output Input file with 14 columns in tab-delimited format, output from BLASTX using parameters: 
##' -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe'  -max_target_seqs 1
##' @param evalue_cutoff Threshold of E-value for BLASTX hits, default 1e-5
##' @param frameDist_cutoff Threshold for gap size (bp) to detect frameshifting between BLASTX hits of same mRNA/cDNA sequence, default 10
##' @return dataframe
##' @importFrom stats complete.cases
##' @export 
##' @author Xiao Chen
##' @references 1. X Chen, Y Jiang, F Gao, W Zheng, TJ Krock, NA Stover, C Lu, LA Katz & W Song (2019). 
##' Genome analyses of the new model protist Euplotes vannus focusing on genome rearrangement and resistance 
##' to environmental stressors. Molecular Ecology Resources, 19(5):1292-1308. 
##' <https://doi.org/10.1111/1755-0998.13023>
##' @examples
##' FScanR(FScanR:::test_data)

## Main
FScanR <- function(blastx_output    = FScanR:::test_data, 
                   evalue_cutoff    = 1e-5, 
                   frameDist_cutoff = 10
                   ) {

    blastx <- blastx_output
    colnames(blastx) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
    blastx <- blastx[complete.cases(blastx) & blastx$evalue <= evalue_cutoff,,drop=FALSE]

    blastx_freq <- table(blastx$qseqid)
    blastx_freq_cutoff <- blastx_freq[blastx_freq > 1]

    blastx_cutoff <- blastx[blastx$qseqid %in% names(blastx_freq_cutoff),,drop=FALSE]

    blastx_cutoff_sort <- blastx_cutoff[order(blastx_cutoff$qseqid, blastx_cutoff$sseqid, blastx_cutoff$qstart),,drop=FALSE]

    prf <- data.frame()
    if (nrow(blastx_cutoff_sort) > 0) {
        for (i in 2:nrow(blastx_cutoff_sort)) {
            qseqid <- blastx_cutoff_sort[i,1]
            sseqid <- blastx_cutoff_sort[i,2]
            qstart <- blastx_cutoff_sort[i,7]
            qend <- blastx_cutoff_sort[i,8]
            sstart <- blastx_cutoff_sort[i,9]
            send <- blastx_cutoff_sort[i,10]
            qframe <- blastx_cutoff_sort[i,13]
            qseqid_last <- blastx_cutoff_sort[i-1,1]
            sseqid_last <- blastx_cutoff_sort[i-1,2]
            qstart_last <- blastx_cutoff_sort[i-1,7]
            qend_last <- blastx_cutoff_sort[i-1,8]
            sstart_last <- blastx_cutoff_sort[i-1,9]
            send_last <- blastx_cutoff_sort[i-1,10]
            qframe_last <- blastx_cutoff_sort[i-1,13]

            if (qseqid == qseqid_last & sseqid == sseqid_last & qframe != qframe_last & qframe * qframe_last > 0) {
                if (qframe > 0 & qframe_last > 0) {
                    frameStart <- qend_last
                    frameEnd <- qstart
                } else if (qframe < 0 & qframe_last < 0) {
                    frameStart <- qstart_last
                    frameEnd <- qend
                }
                qDist <- frameEnd - frameStart - 1
                sDist <- sstart - send_last
                qframe_last <- ifelse(qframe_last %in% c(-3, 3), 0, qframe_last)
                qframe <- ifelse(qframe %in% c(-3, 3), 0, qframe)
                FS_type <- qDist + (1 - sDist) * 3
                if (qDist >= frameDist_cutoff * (-1) & qDist <= frameDist_cutoff & sDist <= floor(frameDist_cutoff/3) & sDist >= floor(frameDist_cutoff/3) * (-1)) {
                    prf_sub <- data.frame(as.character(qseqid), frameStart, frameEnd, as.character(sseqid), send_last + 1, sstart, FS_type)
                    prf <- rbind(prf, prf_sub)
                }
                prf <- prf[prf$FS_type < 3 & prf$FS_type > -3,,drop=FALSE]
            }
        }
        colnames(prf) <- c("DNA_seqid", "FS_start", "FS_end", "Pep_seqid", "Pep_FS_start", "Pep_FS_end", "FS_type")
    } else {
        message("No PRF events detected!")
    }
    
    return(prf)
}
