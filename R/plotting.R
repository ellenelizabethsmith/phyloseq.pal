#various functions for plotting phyloseq objects

library(ggplot2)

#Nice high contrast 20 colours
nice20 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',   '#000075', '#808080')
#larger set from rcolor brewer - a bit ugly but necessary
bigcolset<-c(RColorBrewer::brewer.pal(8,"Set2"), RColorBrewer::brewer.pal(12,"Set3"), RColorBrewer::brewer.pal(9,"Set1"), RColorBrewer::brewer.pal(12,"Paired"))
#a last resort!!!!!!
ridiculouslybigcolset <- c(rep(c(nice20,bigcolset),10))

theme_pp <- function (base_size = 12, base_family = "", xlabels = TRUE)
{
  res <- ggplot2::theme_bw(base_size = base_size,base_family = base_family) +
    theme(strip.background =element_rect(fill="white"))+
    theme(axis.text.y = element_blank())+
    theme(axis.ticks = element_blank())+
    scale_x_discrete(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme(panel.spacing.x = unit(0,"lines"))
  if (!xlabels) {
    res <- res +
    theme(axis.text.x = element_blank())+
    theme(axis.title.x = element_blank())
  }
  res
}



#' theme_pp
#'
#' @param base_size
#' @param base_family
#' @param xlabels
#'
#' @return
#' @export
#'
#' @examples
theme_pp <- function (base_size = 12, base_family = "", xlabels = TRUE)
{
  res <- list(
    theme_bw(),
    theme(strip.background =element_rect(fill="white")),
    theme(axis.text.y = element_blank()),
    theme(axis.ticks = element_blank()),
    scale_x_discrete(expand = c(0,0)),
    scale_y_continuous(expand = c(0,0)),
    theme(axis.title.x = element_blank()),
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
    theme(legend.position = "bottom"),
    theme(panel.spacing = unit(0, "lines")),  # Remove gaps between panels
    theme(text = element_text(size = base_size))
  )
  if (!xlabels) {
    res <- c(res,list(
      theme(axis.text.x = element_blank()),
      theme(axis.title.x = element_blank())
    )
    )
  }
  res
}

#' DADA2 output to physloseq object
#'
#' @param directory
#' @param metadata
#' @param ranks
#'
#' @return
#' @export
#'
#' @examples
ps_from_ampliseq <- function(directory,metadata=NULL,ranks){
  #asvs
  asv <- phyloseq::otu_table(as.matrix(read.table(paste0(directory,"/dada2/ASV_table.tsv"),header = TRUE,row.names = 1)),taxa_are_rows = TRUE)

  #taxonomy
  tax <- read.csv(paste0(directory,"/dada2/ASV_tax.tsv"),header = TRUE,sep = "\t",row.names = 1)
  tax <- tax[,1:(ncol(tax)-2)]
  tax <- phyloseq::tax_table(as.matrix(tax))

  if(!is.null(metadata)){
    #metadata - it's on you to format this properly.
    metadata <- phyloseq::sample_data(metadata)
    ps <- phyloseq::phyloseq(asv,tax,metadata)
  }else{
    ps <- phyloseq::phyloseq(asv,tax)
  }

  ps@tax_table@.Data[is.na(ps@tax_table@.Data)] <- "Unassigned"
  ps@tax_table@.Data[ps@tax_table@.Data == ""] <- "Unassigned"

  return(ps)
}


#Strongly suspect this is bad practise
#' Return 20 pretty colours
#'
#' @return
#' @export
#'
#' @examples
get_nice_20 <- function(){
  return(c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',   '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',   '#000075', '#808080'))
}

#' Return a lot of ugly colours
#'
#' @return
#' @export
#'
#' @examples
get_big_colours <- function(){
  bigcolset<-c(RColorBrewer::brewer.pal(8,"Set2"), RColorBrewer::brewer.pal(12,"Set3"), RColorBrewer::brewer.pal(9,"Set1"), RColorBrewer::brewer.pal(12,"Paired"))
  #a last resort!!!!!!
  ridiculouslybigcolset <- c(rep(c(nice20,bigcolset),10))
  return(ridiculouslybigcolset)
  }

#' Get Gloms
#'
#' Pre-calculate and return a named list of phyloseq objects agglomerated at specified taxonomic ranks.
#'
#' @param ps A phyloseq object.
#' @param ranks A character vector specifying the taxonomic ranks at which agglomeration should be performed. Default is c("Phylum", "Class", "Order", "Family", "Genus").
#'
#' @return A named list of phyloseq objects, each representing the agglomeration at a specified taxonomic rank.
#' @export
#'
#' @examples
#' # Example usage:
#' gloms_list <- getgloms(physeq_obj, ranks = c("Phylum", "Class", "Order", "Family", "Genus"))
#'
#' @note The resulting list is named based on the specified taxonomic ranks.
#'
#' @seealso
#' \code{\link{phyloseq}}, \code{\link{speedyseq::tax_glom}}
#'
getgloms <- function(ps,ranks= c("Phylum" ,"Class","Order","Family","Genus")){
  ps@tax_table@.Data[is.na(ps@tax_table@.Data)] <- "Unassigned"
  gloms <- lapply(ranks, function(x) speedyseq::tax_glom(ps,taxrank=x,NArm=FALSE)) #agglomerate counts per taxon at each of the specied ranks
  names(gloms) <- ranks
  return(gloms)
}

#' Melt to top n taxa
#' Takes a phyloseq object and returns a melted table with only the top n taxa labelled.
#' Taxa with assignments at the specified rank that are not in the top n will be labelled "Other".
#' Those without assignments will be labelled "Unassigned"
#'
#' @param ps
#' @param n
#' @param rank
#'
#' @return
#' @export
#'
#' @examples
melt_to_top_n <- function(ps,n,rank){
  #check if phyloseq agglomerated at appropriate rank
  #and if not, perform agglomeration
  if(!length(unique(as.data.frame(ps@tax_table@.Data)[,rank])) == length(unique(rownames((as.data.frame(ps@tax_table@.Data)))))){
    print("Agglomerating at specified rank...")
    ps@tax_table@.Data[is.na(ps@tax_table@.Data)] <- "Unassigned"
    glom <-speedyseq::tax_glom(ps, taxrank = rank, NArm = FALSE)
    print("Done.")
  }  else if(rank == "OTU"){
    print("Using OTUs")
    glom <- ps
    tt <- data.frame(glom@tax_table@.Data)
    tt$OTU <- rownames(tt)
    tt <- as.matrix(tt)
    glom@tax_table <- tax_table(tt)
  }  else {
    print("Using provided agglomerated object.")
    glom <- ps
  }

  tt <- data.frame(glom@tax_table@.Data)
  keep <- rownames(tt[tt[[rank]]!="Unassigned",])
  ass_only <- prune_taxa(keep, glom)
  topntaxa <- tt[names(sort(taxa_sums(ass_only), TRUE)[1:min(nrow(ass_only@tax_table),
                                                             n)]),][[rank]]

  #create a taxonomy table containing these ASVs with everything else set to "Other"
  tt <- data.frame(ps@tax_table@.Data)
  tt$taxon <- "Other"
  tt[tt[[rank]] == "Unassigned",]$taxon <- "Unassigned"
  tt[tt[[rank]] %in% topntaxa,]$taxon <- tt[tt[[rank]] %in% topntaxa,][[rank]]


  #taxtabn = data.frame(cbind(tax_table(glom), taxon = "Other"))
  #taxtabn[topnotus, "taxon"] <- as(tax_table(glom)[topnotus, rank],"character")

  tax_table(ps) <- tax_table(as.matrix(tt))
  melt <- speedyseq::psmelt(ps)

  #get names for reordering factors
  labels <- (unique(melt$taxon))
  nlabs <- length(labels)

  #move "Unassigned" to end
  if ("Unassigned" %in% labels){
    labels <- labels[!labels %in% "Unassigned"]
    labels[nlabs] = "Unassigned"
  }
  #if > n taxa, move "Other" to end
  if("Other" %in% melt$taxon){
    labels <- labels[!labels %in% "Other"]
    labels[nlabs] <- "Other"
  }
  for(i in seq(1:nlabs)){ #making other and unassigned have consistent colours
    if(labels[i] == "Unassigned"){
    } else if(labels[i] == "Other"){
    }
  }

  #labels by abundance w/ unassigned and otehr at end
  melt$taxon <- forcats::fct_relevel(melt$taxon,labels)

  return(melt)

}




#' Plot Taxa Abundance
#'
#' This function generates stacked bar plots to visualize the abundance of taxa in a phyloseq object.
#'
#' @param ps A phyloseq object.
#' @param rank (default: Genus) The taxonomic rank at which agglomeration should be performed if necessary.
#' @param x (default: SampleID) The variable to be plotted on the x-axis.
#' @param wrap (Optional) A variable to be used for faceting the plot.
#' @param n (Optional, defualt: 20) The number of top taxa to display in the plot.
#' @param byabundance If TRUE, taxa will be ordered by decreasing abundance; if FALSE, they will be ordered by taxonomic rank.
#' @param abs If TRUE, plots will use absolute rather than relative counts per sample.
#' @param size The font size for plot labels.
#'
#' @return A ggplot object representing the taxa abundance plot.
#' @export
#'
#' @examples
#' # Example usage:
#' data("GlobalPatterns")
#` ps <- subset_taxa(GlobalPatterns, Phylum %in% c("Proteobacteria", "Bacteroidetes"))
#` plot_taxa_abundance(psx, "Genus", "Sample", wrap = "SampleType", n = 15, abs = FALSE)
#'
#' @note If the phyloseq object is not agglomerated at the specified rank, the function will perform agglomeration.

#' @seealso
#' \code{\link{phyloseq}}, \code{\link{psmelt}}, \code{\link{tax_glom}}, \code{\link{ggplot2}}
#'
#' @importFrom ggplot2 aes_string geom_bar scale_fill_manual labs theme scale_y_continuous scale_x_discrete facet_grid
#'
plot_taxa_abundance <- function(ps,rank,x, wrap = NULL, n=20, byabundance=TRUE,abs=FALSE,size=10){

  #set cols
  cols.n <- c(c(nice20, ridiculouslybigcolset)[1:n],"lightgrey")

  #check if phyloseq agglomerated at appropriate rank
  #and if not, perform agglomeration
  if(!length(unique(as.data.frame(ps@tax_table@.Data)[,rank])) == length(unique(rownames((as.data.frame(ps@tax_table@.Data)))))){
    print("Agglomerating at specified rank...")
    ps@tax_table@.Data[is.na(ps@tax_table@.Data)] <- "Unassigned"
    glom <-speedyseq::tax_glom(ps, taxrank = rank, NArm = FALSE)
    print("Done.")
  }  else if(rank == "OTU"){
    print("Using OTUs")
    glom <- ps
    tt <- data.frame(glom@tax_table@.Data)
    tt$OTU <- rownames(tt)
    tt <- as.matrix(tt)
    glom@tax_table <- tax_table(tt)
  }  else {
    print("Using provided agglomerated object.")
    glom <- ps
  }

  tt <- data.frame(glom@tax_table@.Data)
  keep <- rownames(tt[tt[[rank]]!="Unassigned",])
  ass_only <- prune_taxa(keep, glom)
  topnotus <- names(sort(taxa_sums(ass_only), TRUE)[1:min(nrow(ass_only@tax_table),
                                                          n)])

  #create a taxonomy table containing these ASVs with everything else set to "Other"
  tt <- data.frame(glom@tax_table@.Data)
  tt$taxon <- "Other"
  tt[tt[[rank]] == "Unassigned",]$taxon <- "Unassigned"
  tt[rownames(tt) %in% topnotus,]$taxon <- tt[rownames(tt) %in% topnotus,][[rank]]


  #taxtabn = data.frame(cbind(tax_table(glom), taxon = "Other"))
  #taxtabn[topnotus, "taxon"] <- as(tax_table(glom)[topnotus, rank],"character")

  tax_table(glom) <- tax_table(as.matrix(tt))
  melt <- speedyseq::psmelt(glom)

  #get names for reordering factors
  labels <- (unique(melt$taxon))
  nlabs <- length(labels)

  #move "Unassigned" to end
  if ("Unassigned" %in% labels){
    labels <- labels[!labels %in% "Unassigned"]
    labels[nlabs] = "Unassigned"
  }
  #if > n taxa, move "Other" to end
  if("Other" %in% melt$taxon){
    labels <- labels[!labels %in% "Other"]
    labels[nlabs] <- "Other"
  }
  for(i in seq(1:nlabs)){ #making other and unassigned have consistent colours
    if(labels[i] == "Unassigned"){
      cols.n[i] <- "grey"
    } else if(labels[i] == "Other"){
      cols.n[i] <- "lightblue"
    }
  }

  #labels by abundance w/ unassigned and otehr at end
  melt$taxon <- forcats::fct_relevel(melt$taxon,labels)

  #make the stacked bar chart
  i <- ggplot(melt, aes_string(x = x, y = "Abundance", fill = "taxon")) +
    geom_bar(stat = "identity", width = 1, position = position_fill()) +
    scale_fill_manual(values=cols.n, na.value = "grey")+
    theme(axis.title.x = element_blank())+
    labs(fill=rank)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "bottom")+
    theme(text = element_text(size = size))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))

  if(abs){
    i <- i+ geom_bar(stat = "identity", width = 1, position = "stack")
  }

  if (is.character(wrap)){
    i <- i + facet_grid(as.formula(paste("~", wrap)), scales = "free_x", space = "free")
  }

  return(i)
}

#' Merge Two Differently oriented DADA2 runs from the same samples
#'
#' @param dada_fwd
#' @param dada_rvs
#'
#' @return
#' @export
#'
#' @examples
merge_orientations <- function(dada_fwd, dada_rvs){

  rvs_seqs <- Biostrings::DNAStringSet(rownames(dada_rvs))
  rc <- reverseComplement(rvs_seqs)

  rownames(dada_rvs) <- as.character(rc)


  length(intersect(rownames(dada_rvs),rownames(dada_fwd)))
  length(rownames(dada_rvs))
  length(rownames(dada_fwd))

  combined <- dada2::mergeSequenceTables(dada_fwd,dada_rvs,repeats = "sum")
  length(rownames(combined))

  combined.derep <- dada2::collapseNoMismatch(combined,verbose = TRUE)
  length(rownames(combined.derep))

  summary <- data.frame(fwd=length(rownames(dada_fwd)),
                        rev=length(rownames(dada_rvs)),
                        identical.after.rc=length(intersect(rownames(dada_rvs),rownames(dada_fwd))),
                        expected=(length(rownames(dada_fwd))+ length(rownames(dada_rvs)))-length(intersect(rownames(dada_rvs),rownames(dada_fwd))),
                        after.merge=length(rownames(combined)),
                        after.collapse=length(rownames(combined.derep))
  )
  rownames(summary) <- "ASVS"
  print(summary)

  return(combined.derep)

}

#' Create Microshades Plot
#'
#' This function generates a Microshades plot using data prepared with the microshades package. Microshades is a package designed for analyzing microbiome data, particularly for visualizing taxonomic abundance across samples.
#'
#' @param ps An object representing the microbiome data. It should be preprocessed with the microshades package.
#' @param rankhigh The higher taxonomic rank to consider for grouping taxa. Default is "Phylum".
#' @param ranklow The lower taxonomic rank to consider for grouping taxa. Default is "Genus".
#' @param x The variable to be plotted on the x-axis.
#' @param facet An optional variable to facet the plot by. Default is NULL.
#'
#' @return A Microshades plot visualizing taxonomic abundance across samples.
#' @export
#'
#' @examples
#' # Load necessary libraries and data
#' data(GlobalPatterns)
#'
#' # Create Microshades plot
#' do_microshades_plot(ps = GlobalPatterns, rankhigh = "Phylum", ranklow = "Genus", x = "SampleID", facet = SampleType)
#'
#' @importFrom microshades prep_mdf create_color_dfs custom_legend plot_microshades
#' @importFrom ggplot2 scale_y_discrete theme element_text facet_grid
#' @importFrom ggpubr ggarrange
#'
#' @seealso
#' \code{\link{prep_mdf}}, \code{\link{melt_to_top_n}}, \code{\link{create_color_dfs}}, \code{\link{custom_legend}}, \code{\link{plot_microshades}}
#'
#' @references
#' Microshades package: https://github.com/KarstensLab/microshades
#'
#' Karstens, L., Asquith, M., Davin, S., Stauffer, P., Fair, D., & Gregory, W. T. (2018). Microbiome analyses and
#' endophenotypic associations in colorectal cancer. Gut Pathogens, 10(1), 1-14. https://doi.org/10.1186/s13099-018-0275-6
#'
#' @keywords microshades plot microbiome
#' @family visualization
do_microshades_plot <- function(ps, rankhigh = "Phylum", ranklow = "Genus", x, facet = NULL) {

  # Prepare microbiome data
  mdf <- microshades::prep_mdf(ps)

  # Extract top taxa
  top <- melt_to_top_n(ps, n = 5, rank = rankhigh)

  # Arrange taxa
  fam <- rev(as.character(unique(top$taxon[!top$taxon %in% c("Other", "Unassigned")])))

  # Create color palette
  color_obj <- microshades::create_color_dfs(mdf = mdf, selected_groups = fam, group_level = rankhigh, subgroup_level = ranklow, top_n_subgroups = 4, cvd = TRUE)

  # Extract color dataframes
  mdf_group <- color_obj$mdf
  cdf <- color_obj$cdf

  # Create custom legend
  leg <- microshades::custom_legend(mdf = color_obj$mdf, cdf = color_obj$cdf, group_level = rankhigh, subgroup_level = ranklow, x = "sample_alias", legend_text_size = 20)

  # Generate Microshades plot
  p <- microshades::plot_microshades(mdf_group, cdf, group_label = paste0(rankhigh, " ", ranklow), x = x, y = "Abundance") +
    scale_y_discrete(expand = c(0, 0)) +
    theme_pp() +
    theme(text = element_text(size = 15), legend.position = "none") +
    ylab(paste(rankhigh, ranklow, sep = " "))

  # Facet the plot if specified
  if (is.character(facet)) {
    p <- p + facet_grid(as.formula(paste("~", facet)), scales = "free_x", space = "free")
  }

  # Arrange plot and legend
  ggpubr::ggarrange(p, leg, widths = c(0.6, 0.1))
}


#Bray curtis plots. PCOA is table from "ordinate(ps)".
#' Bray Plot
#'
#' @param ps
#' @param pcoa
#' @param colour
#' @param shape
#'
#' @return
#' @export
#'
#' @examples
brayplot_pp <- function(ps, pcoa, colour, shape=NULL){
  p <- plot_ordination(ps, pcoa, color = colour, axes = c(1,2),shape = shape) +
    geom_point(size = 2) +
    labs(title = colour, color = colour)+
    scale_color_manual(values = nice20)


  return(p)
}

#' Calculate the percentage of assigned taxa at each taxonomic rank
#'
#' This function calculates the percentage of assigned taxa at each taxonomic rank
#' based on a given phyloseq object.
#'
#' @param ps A phyloseq object.
#'
#' @return A data frame containing the percentage of assigned taxa at each taxonomic rank.
#'
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' get_percentage_assigned(GlobalPatterns)
#'
#' @import phyloseq
#' @export
get_percentage_assigned <- function(ps) {

  ps@tax_table@.Data[is.na(ps@tax_table@.Data)] <- "Unassigned"
  ps@tax_table@.Data[ps@tax_table@.Data == ""] <- "Unassigned"

  table <- as.data.frame(ps@tax_table)
  pct_un <- list()
  for (rank in colnames(table)) {
    l <- length(table[[rank]])
    u <- length(table[[rank]][table[[rank]] != "Unassigned"])
    pct_un[[rank]] <- round(u / l * 100, 1)
  }
  pct_tab <- data.frame(Rank = names(pct_un), Pct.Assigned = unlist(pct_un))
  return(pct_tab)
}


#alpha diverstiy plots. richtable is from estimate_richness(ps).
#' Alpha Diversity Plot
#'
#' @param richtable
#' @param xvar
#' @param yvar
#' @param wrapvar
#' @param colno
#'
#' @return
#' @export
#'
#' @examples
alphaplot_pp <- function(richtable, xvar, yvar, wrapvar = NULL, colno = NULL){
  p <- ggplot(richtable, aes_string(x = xvar, y = yvar, fill = xvar)) +
    geom_boxplot() +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = "none")

  if (is.character(wrapvar)){
    p <- p + facet_grid(as.formula(paste("~", wrapvar)),scales = "free_x")
  }
  return(p)
}

#for generating a table from kruskall significance tests. richtable is generated by estimate_richness(ps)
#' Do Kruskall Tests
#'
#' @param richtable
#' @param testvars
#' @param metrics
#'
#' @return
#' @export
#'
#' @examples
makekruskalltests <- function(richtable,testvars,metrics){
  #create df
  df <- data.frame()
  for (k in metrics) df[[k]] <- as.character()

  for (v in testvars){
    results <- vector()
    for(i in metrics){
      results[[i]] <- kruskal.test(formula(paste(i, "~ ", v)), data = richtable)$p.value
    }
    newrow <- nrow(df) + 1
    df[newrow,] <- results
    rownames(df)[newrow] <- v
  }
  return(df)
}

#plot expression scatter from DESEq2
deseq2_scatter <- function(sigtab,dds,x,wrap = NULL){
  plots <- list()
  for(t in rownames(sigtab)){
    title <- paste(sigtab[t,]$Family, sigtab[t,]$Genus,sep = "_")
    data <- plotCounts(dds = dds,gene = t,intgroup = colnames(dds@colData), main= title, returnData = TRUE)

    p <- ggplot(data, aes_string(x=x, y="count", color=x)) +
      scale_y_log10() +
      theme_classic() +
      geom_point(position=position_jitter(width=.1,height=0), size=5)+
      ggtitle(title) +
      theme(legend.position = "none") +
      theme(text = element_text(size = 10))

    if (is.character(wrap)){
      p <- p + facet_grid(as.formula(paste("~", wrap)), scales = "free_x", space="free")
    }
    plots[[t]] <- p
  }
  return(plots)
}

tax_to_box <- function(x,wrap,data){
  p <- ggplot(data, aes_string(x=x, y="count", fill=x)) +
    scale_y_log10() +
    theme_classic() +
    geom_boxplot()+
    geom_point(position=position_jitter(width=.1,height=0), size=5)+
    ggtitle(title) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 15))

  if (is.character(wrap)){
    p <- p + facet_grid(as.formula(paste("~", wrap)), scales = "free_x", space="free")
  }
  return(p)
}

#plot expression box from DESEq2
deseq2_box <- function(sigtab,dds,x,wrap = NULL){
  plots <- list()
  for(t in rownames(sigtab)){
    title <- paste(sigtab[t,]$Order,sigtab[t,]$Family, sigtab[t,]$Genus,sep = "_")
    data <- plotCounts(dds = dds,gene = t,intgroup = colnames(dds@colData), main= title, returnData = TRUE)

    p <- ggplot(data, aes_string(x=x, y="count", fill=x)) +
      scale_y_log10() +
      theme_classic() +
      geom_boxplot()+
      geom_point(position=position_jitter(width=.1,height=0), size=5)+
      ggtitle(title) +
      theme(legend.position = "none") +
      theme(text = element_text(size = 15))

    if (is.character(wrap)){
      p <- p + facet_grid(as.formula(paste("~", wrap)), scales = "free_x", space="free")
    }
    plots[[t]] <- p
  }
  return(plots)
}


sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

#function for plotting DA_bars from a sig table generated in DESeq2
da_bars <- function(sigtab, title,cols = NULL, x_limits = NULL,taxa_order = NULL,setcolor = FALSE,orderbyphylum = FALSE){
  bar_width =.8
  sigtabgen = subset(sigtab, !is.na(Genus) & Genus != "Unassigned" & Genus != "")

  # Specify custom order for Genus
  if (!is.null(taxa_order)) {
    sigtabgen <- sigtabgen[sigtabgen$Genus %in% taxa_order,]
    sigtabgen$Genus = factor(sigtabgen$Genus, (levels = rev(taxa_order)))

  } else if (orderbyphylum) {

    # Order the data by Phylum and then by log2FoldChange within each Phylum
    sigtabgen <- sigtabgen[rev(order(sigtabgen$Phylum, -sigtabgen$log2FoldChange)), ]

    # Create a custom factor order based on the sorted data
    genus_order <- sigtabgen$Genus

    # Reorder Genus factor levels based on the custom order
    sigtabgen$Genus <- factor(sigtabgen$Genus, levels = genus_order)

  } else {
    x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
    x = sort(x, TRUE)
    sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = rev(names(x)))
  }

  num_taxa = length(unique(sigtabgen$Genus))
  plot_height = max(2,num_taxa * 2.5)  # Minimum height of 3 units

  p <- ggplot(sigtabgen, aes(y = Genus, x = log2FoldChange, color = Phylum, fill = Phylum)) +
    theme_pubr() +
    geom_vline(xintercept = 0.0, color = "black", size = 1) +
    geom_bar(stat = "identity", position = position_dodge(width = bar_width), width = bar_width) +
    theme(axis.text.x = element_text(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 15)) +
    ggtitle(title) +
    guides(fill = guide_legend(override.aes = list(size = 5))) #+
  #coord_fixed(ratio = plot_height / num_taxa,expand = FALSE)

  if(! is.null(x_limits)){
    p <- p +
      xlim(x_limits)
  }

  if(!is.null(cols)){
    p <- p +
      scale_fill_manual(values = cols) +
      scale_color_manual(values = cols)
  }

  #this allows consistent colours across plots
  if(setcolor){
    p <- p +
      scale_fill_manual(values = c(Firmicutes = "orange",Actinobacteria="blue",Proteobacteria="purple",Bacteroidetes= "forestgreen")) +
      scale_color_manual(values = c(Firmicutes = "orange",Actinobacteria="blue",Proteobacteria="purple",Bacteroidetes= "forestgreen"))
  }

  return(p)
}

