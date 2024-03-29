% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{do_microshades_plot}
\alias{do_microshades_plot}
\title{Create Microshades Plot}
\usage{
do_microshades_plot(
  ps,
  rankhigh = "Phylum",
  ranklow = "Genus",
  x,
  facet = NULL
)
}
\arguments{
\item{ps}{An object representing the microbiome data. It should be preprocessed with the microshades package.}

\item{rankhigh}{The higher taxonomic rank to consider for grouping taxa. Default is "Phylum".}

\item{ranklow}{The lower taxonomic rank to consider for grouping taxa. Default is "Genus".}

\item{x}{The variable to be plotted on the x-axis.}

\item{facet}{An optional variable to facet the plot by. Default is NULL.}
}
\value{
A Microshades plot visualizing taxonomic abundance across samples.
}
\description{
This function generates a Microshades plot using data prepared with the microshades package. Microshades is a package designed for analyzing microbiome data, particularly for visualizing taxonomic abundance across samples.
}
\examples{
# Load necessary libraries and data
data(GlobalPatterns)

# Create Microshades plot
do_microshades_plot(ps = GlobalPatterns, rankhigh = "Phylum", ranklow = "Genus", x = "SampleID", facet = SampleType)

}
\references{
Microshades package: https://github.com/KarstensLab/microshades

Karstens, L., Asquith, M., Davin, S., Stauffer, P., Fair, D., & Gregory, W. T. (2018). Microbiome analyses and
endophenotypic associations in colorectal cancer. Gut Pathogens, 10(1), 1-14. https://doi.org/10.1186/s13099-018-0275-6
}
\seealso{
\code{\link{prep_mdf}}, \code{\link{melt_to_top_n}}, \code{\link{create_color_dfs}}, \code{\link{custom_legend}}, \code{\link{plot_microshades}}
}
\concept{visualization}
\keyword{microbiome}
\keyword{microshades}
\keyword{plot}
