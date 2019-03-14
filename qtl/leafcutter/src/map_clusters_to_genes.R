# library(leafcutter, quietly=TRUE)
suppressMessages(library(argparser, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))
suppressMessages(library(foreach, quietly=TRUE))

# leafcutter functions:

#' Make a data.frame of meta data about the introns
#' @param introns Names of the introns
#' @return Data.frame with chr, start, end, cluster id
#' @export
get_intron_meta <- function(introns) {
  intron_meta <- do.call(rbind, strsplit(introns,":"))
  colnames(intron_meta) <- c("chr","start","end","clu")
  intron_meta <- as.data.frame(intron_meta, stringsAsFactors=FALSE)
  intron_meta$start <- as.numeric(intron_meta$start)
  intron_meta$end <- as.numeric(intron_meta$end)
  intron_meta
}

#' Work out which gene each cluster belongs to. Note the chromosome names used in the two inputs must match.
#' @param intron_meta Data frame describing the introns, usually from get_intron_meta
#' @param exons_table Table of exons, see e.g. /data/gencode19_exons.txt.gz
#' @return Data.frame with cluster ids and genes separated by commas
#' @import dplyr
#' @export
map_clusters_to_genes <- function(intron_meta, exons_table) {
  gene_df <- foreach (chr=sort(unique(intron_meta$chr)), .combine=rbind) %dopar% {

    intron_chr <- intron_meta[ intron_meta$chr==chr, ]
    exons_chr <- exons_table[exons_table$chr==chr, ]

    exons_chr$temp <- exons_chr$start
    intron_chr$temp <- intron_chr$end
    three_prime_matches <- inner_join( intron_chr, exons_chr, by="temp")

    exons_chr$temp <- exons_chr$end
    intron_chr$temp <- intron_chr$start
    five_prime_matches <- inner_join( intron_chr, exons_chr, by="temp")

    all_matches <- rbind(three_prime_matches, five_prime_matches)[ , c("clu", "gene_name")]

    all_matches <- all_matches[!duplicated(all_matches),]

    if (nrow(all_matches)==0) return(NULL)
    all_matches$clu <- paste(chr,all_matches$clu,sep=':')
    all_matches
  }

  clu_df <- gene_df %>% group_by(clu) %>% summarize(genes=paste(gene_name, collapse = ","))
  class(clu_df) <- "data.frame"
  clu_df
}


p <- arg_parser("LeafCutter: map clusters to genes")
p <- add_argument(p, "intron_counts_file", help="Intron counts file from LeafCutter, typically <prefix>_perind.counts.gz")
p <- add_argument(p, "exon_file", help="File listing all unique exons in annotation. Must have columns: chr, start, end, strand, gene_id[, gene_name].")
p <- add_argument(p, "output_name", help="Output file name")
p <- add_argument(p, "--output_dir", short="-o", help="Output directory", default=".")
argv <- parse_args(p)

cat("LeafCutter: mapping clusters to genes\n")
intron_counts <- read.table(argv$intron_counts_file, header=TRUE, check.names=FALSE, row.names=1)
intron_meta <- get_intron_meta(rownames(intron_counts))

exon_table <- read.table(argv$exon_file, header=TRUE, stringsAsFactors=FALSE)
stopifnot(is.element('gene_id', colnames(exon_table)))
exon_table[, 'gene_name'] <- exon_table[, 'gene_id']

m <- map_clusters_to_genes(intron_meta, exon_table)
write.table(m, file.path(argv$output_dir, argv$output_name), sep = "\t", quote=FALSE, row.names=FALSE)
