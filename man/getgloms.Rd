% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{getgloms}
\alias{getgloms}
\title{Get Gloms}
\usage{
getgloms(ps, ranks = c("Phylum", "Class", "Order", "Family", "Genus"))
}
\arguments{
\item{ps}{A phyloseq object.}

\item{ranks}{A character vector specifying the taxonomic ranks at which agglomeration should be performed. Default is c("Phylum", "Class", "Order", "Family", "Genus").}
}
\value{
A named list of phyloseq objects, each representing the agglomeration at a specified taxonomic rank.
}
\description{
Pre-calculate and return a named list of phyloseq objects agglomerated at specified taxonomic ranks.
}
\note{
The resulting list is named based on the specified taxonomic ranks.
}
\examples{
# Example usage:
gloms_list <- getgloms(physeq_obj, ranks = c("Phylum", "Class", "Order", "Family", "Genus"))

}
\seealso{
\code{\link{phyloseq}}, \code{\link{speedyseq::tax_glom}}
}
