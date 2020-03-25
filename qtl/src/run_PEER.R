# Author: Francois Aguet

library(peer, quietly=TRUE)  # https://github.com/PMBio/peer
library(argparser, quietly=TRUE)

WriteTable <- function(data, filename, index.name) {
    datafile <- file(filename, open = "wt")
    on.exit(close(datafile))
    header <- c(index.name, colnames(data))
    writeLines(paste0(header, collapse="\t"), con=datafile, sep="\n")
    write.table(data, datafile, sep="\t", col.names=FALSE, quote=FALSE)
}

p <- arg_parser("Run PEER factor estimation")
p <- add_argument(p, "expr.file", help="")
p <- add_argument(p, "prefix", help="")
p <- add_argument(p, "n", help="Number of hidden confounders to estimate")
p <- add_argument(p, "--covariates", help="Observed covariates")
p <- add_argument(p, "--alphaprior_a", help="", default=0.001)
p <- add_argument(p, "--alphaprior_b", help="", default=0.01)
p <- add_argument(p, "--epsprior_a", help="", default=0.1)
p <- add_argument(p, "--epsprior_b", help="", default=10)
p <- add_argument(p, "--max_iter", help="", default=1000)
p <- add_argument(p, "--output_dir", short="-o", help="Output directory", default=".")
argv <- parse_args(p)

cat("PEER: loading expression data ... ")
if (grepl('.gz$', argv$expr.file)) {
    nrows <- as.integer(system(paste0("zcat ", argv$expr.file, " | wc -l | cut -d' ' -f1 "), intern=TRUE, wait=TRUE))
} else {
    nrows <- as.integer(system(paste0("wc -l ", argv$expr.file, " | cut -d' ' -f1 "), intern=TRUE, wait=TRUE))
}
if (grepl('.bed$', argv$expr.file) || grepl('.bed.gz$', argv$expr.file)) {
    df <- read.table(argv$expr.file, sep="\t", nrows=nrows, header=TRUE, check.names=FALSE, comment.char="")
    row.names(df) <- df[, 4]
    df <- df[, 5:ncol(df)]
} else {
    df <- read.table(argv$expr.file, sep="\t", nrows=nrows, header=TRUE, check.names=FALSE, comment.char="", row.names=1)
}
M <- t(as.matrix(df))
cat("done.\n")

# run PEER
cat(paste0("PEER: estimating hidden confounders (", argv$n, ")\n"))
model <- PEER()
invisible(PEER_setNk(model, argv$n))
invisible(PEER_setPhenoMean(model, M))
invisible(PEER_setPriorAlpha(model, argv$alphaprior_a, argv$alphaprior_b))
invisible(PEER_setPriorEps(model, argv$epsprior_a, argv$epsprior_b))
invisible(PEER_setNmax_iterations(model, argv$max_iter))
if (!is.null(argv$covariates) && !is.na(argv$covariates)) {
    has.cov <- TRUE
    covar.df <- read.table(argv$covariates, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
    covar.df[] <- sapply(covar.df, as.numeric)
    cat(paste0("  * including ", dim(covar.df)[2], " covariates", "\n"))
    invisible(PEER_setCovariates(model, as.matrix(covar.df[rownames(M), ])))  # samples x covariates
} else {
    has.cov <- FALSE
}
time <- system.time(PEER_update(model))

X <- PEER_getX(model)  # samples x PEER factors
A <- PEER_getAlpha(model)  # PEER factors x 1
R <- t(PEER_getResiduals(model))  # genes x samples

# add relevant row/column names
if (has.cov) {
    cols <- c(colnames(covar.df), paste0("InferredCov",1:(ncol(X)-dim(covar.df)[2])))
} else {
    cols <- paste0("InferredCov",1:ncol(X))
}
rownames(X) <- rownames(M)
colnames(X) <- cols
rownames(A) <- cols
colnames(A) <- "Alpha"
A <- as.data.frame(A)
A$Relevance <- 1.0 / A$Alpha
rownames(R) <- colnames(M)
colnames(R) <- rownames(M)

# write results
cat("PEER: writing results ... ")
WriteTable(t(X), file.path(argv$output_dir, paste0(argv$prefix, ".PEER_covariates.txt")), "ID")  # format(X, digits=6)
WriteTable(A, file.path(argv$output_dir, paste0(argv$prefix, ".PEER_alpha.txt")), "ID")
WriteTable(R, file.path(argv$output_dir, paste0(argv$prefix, ".PEER_residuals.txt")), "ID")
cat("done.\n")
