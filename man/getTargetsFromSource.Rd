% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbDataAccess.R
\name{getTargetsFromSource}
\alias{getTargetsFromSource}
\title{Get target list from a single source}
\usage{
getTargetsFromSource(mirna, species = "mmu", source = "diana",
  synonyms = TRUE, both_strands = FALSE)
}
\arguments{
\item{mirna}{miRNA in a standard format}

\item{species}{species in a standard three-letter acronym, default 
\code{'mmu'}}

\item{source}{a source target prediction algorithm table to query, default 
\code{'diana'}, other possible values are \code{'miranda'}, 
\code{'targetscan'}, and \code{'pictar'}.}

\item{synonyms}{when searching for -3p miRNA automatically also searches for
miRNA with the same name but ending with * (some databases list -3p miRNA
this way) and other way around, similarly for -5p miRNA, default TRUE}

\item{both_strands}{overrides \code{synonyms} and searches for targets of
both -5p and -3p strands together}
}
\value{
\code{data.frame} object with entrez IDs of target genes and their
scores, if there are no targets found for a given miRNA in a given 
table then an empty
}
\description{
This function queries precompiled annotation SQLite database
which contains miRNA - target gene associations with their respective scores.
}
\examples{
targets <- getTargetsFromSource('let-7a', species='hsa', source='targetscan')
head(targets) 
#top of the listof human targets of let-7a from TargetScan only
}
\author{
Maciej Pajak \email{m.pajak@sms.ed.ac.uk}
}
\references{
Friedman, R. C., Farh, K. K.-H., Burge, C. B., and Bartel, D. P. (2009). 
Most mammalian mRNAs are conserved targets of microRNAs. Genome research, 
19(1):92-105.

Griffiths-Jones, S., Saini, H. K., van Dongen, S., and Enright, A. J. 
(2008). miRBase: tools for microRNA genomics. Nucleic acids research, 
36(Database issue):D154-8.

Lall, S., Grun, D., Krek, A., Chen, K., Wang, Y.-L., Dewey, C. N., ... 
Rajewsky, N. (2006). A genome-wide map of conserved microRNA targets in 
C. elegans. Current biology : CB, 16(5):460-71.

Maragkakis, M., Vergoulis, T., Alexiou, P., Reczko, M., Plomaritou, K., 
Gousis, M., ... Hatzigeorgiou, A. G. (2011). DIANA-microT Web server 
upgrade supports Fly and Worm miRNA target prediction and bibliographic 
miRNA to disease association. Nucleic Acids Research, 39(Web Server issue), 
W145-8.
}

