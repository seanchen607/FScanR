#' BLASTX output
#'
#' A dataset containing the output of homolog sequence alignment by BLASTX.
#'
#' @format A data frame with 38549 rows and 14 variables:
#' \describe{
#'   \item{qseqid}{query sequence ID}
#'   \item{sseqid}{subject sequence ID}
#'   \item{pident}{sequence alignment identity}
#'   \item{length}{sequence alignment length}
#'   \item{mismatch}{sequence alignment mismatch numbers}
#'   \item{gapopen}{gap size within sequence alignment}
#'   \item{qstart}{start location of alignment on query sequence}
#'   \item{qend}{end location of alignment on query sequence}
#'   \item{sstart}{start location of alignment on subject sequence}
#'   \item{send}{end location of alignment on subject sequence}
#'   \item{evalue}{sequence alignment E-value}
#'   \item{bitscore}{sequence alignment Bitscore}
#'   \item{qframe}{frame No. of query sequence within the alignment}
#'   \item{sframe}{frame No. of subject sequence within the alignment}
#' }
#' @source \url{http://www.diamondse.info/}
"test"
