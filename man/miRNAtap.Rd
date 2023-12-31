% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miRNAtap.R
\docType{package}
\name{miRNAtap}
\alias{miRNAtap}
\alias{miRNAtap-package}
\title{miRNAtap: microRNA Targets - Aggregated Predictions.}
\description{
It is a package with tools to facilitate implementation of
workflows
requiring miRNA prediction through access to multiple prediction results 
(DIANA, Targetscan, PicTar, Miranda, and miRDB) and their aggregation. 
Three aggregation methods are available: minimum, maximum and geometric mean,
additional parameters provide further tuning of the results. 
Predictions are available for Homo sapiens, Mus musculus 
and Rattus norvegicus (the last one through homology translation).
}
\examples{
#direct targets in mouse aggregated from all sources:
targets_mouse <- getPredictedTargets('let-7a',species='mmu', method='geom') 
#homology-translated targets in rat aggregated from all sources
targets_rat <- getPredictedTargets('let-7a',species='rno', method='geom') 
}
\author{
Maciej Pajak \email{m.pajak@sms.ed.ac.uk}, Ian Simpson
}

