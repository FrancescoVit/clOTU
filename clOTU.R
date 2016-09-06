#### clOTU, a function for plotting word clouds to represent abundance data in microbial ecology. 

#### Author and mantainer: Francesco Vitali; francesco.vitali.bio@gmail.com
#### Credits: Claudio Ruggeri

#### Basic idea for this function is to substitute the classical barplot representation of microbial
#### abundance when low numbers (1-2) samples are considered. It was designed for the analysis of microbiome
#### data obtained after NGS sequencing, which were analyzed through the QIIME pipeline, and which are handled 
#### in R through phyloseq package.

#### clOTU is still very beta, it is quite slow and plotting is still far from optimal. 

#### ARGUMENTS ####
#### phyloseq.obj = name of the phyloseq object, if data have already been imported in R session 
#### biom.file = if data are not already been imported in R session, path to input biom table obtained by qiime
#### mapping.file = if data are not already imported in R session, path to input mapping file 
#### (e.g. the one created for analysis with qiime)
#### tax.level = taxonomic level at wich the OTU table should be agglomerated for plot (using tax_glom function), 
#### use "OTU" to avoid agglomeration
#### category (optional) = name of the category in mapping file to use for sample agglomeration/grouping
#### to.plot = item to be plotted. It must indicate the column of the otu table to be plotted (i.e. the sample); 
#### can be a number, indicating the number of column, or can be the name of the column.  Can be "OTU" to obtain 
#### a word cloud with otu name


clOTU <- function(phyloseq.obj, biom.file, mapping.file, tax.level, category, to.plot, max.word) {
  
  #loading required packages
  require(phyloseq)
  require(wordcloud)
  require(plyr)
  require(shape)
  require(tm)
  require(NLP)
  
  #Reading in the biom table, the mapping file and creating a phyloseq object.
  if (missing(phyloseq.obj))  {
    biom.file <- import_biom(BIOMfilename=biom.file)
    env <- import_qiime_sample_data(mapfilename = mapping.file)
    phyloseq_obj <- merge_phyloseq(biom.file, env)
  }  else  {
    phyloseq_obj = phyloseq.obj
  }
  #Renaming phylogenetic rank
  colnames(tax_table(phyloseq_obj)) <- c(Rank1 = "Kingdom", Rank2 = "Phylum", Rank3 = "Class", Rank4 = "Order", Rank5 = "Family", Rank6 = "Genus", Rank7 = "Species")
  if (tax.level == "OTU")  {
    otus <- as.data.frame(otu_table(phyloseq_obj))
    taxon <- as.data.frame(tax_table(phyloseq_obj))
    
    #Extracting otus and taxonomy information for the plot 
    j <- data.frame(otus[to.plot],row.names(taxon))
    colnames(j) <- c("Counts","Tax_lvl")
  }  else  {
    #Modifying phyloseq object with provided taxonomy level and sample category
    if (missing(category))    {
      #Using phyloseq agglomerate function to agglomerate OTU table at the chosen phylogenetic rank
      choosen_tax <- tax_glom(physeq = phyloseq_obj, taxrank = tax.level)
      #extracting OTU table and taxon table
      otus <- as.data.frame(otu_table(choosen_tax))
      taxon <- as.data.frame(tax_table(choosen_tax))
      
    } else {  
      #Using phylose agglomerate function to agglomerate samples on a category in mapping file
      agglom_sample <- merge_samples(x = phyloseq_obj, group = as.character(category))
      #Using phyloseq agglomerate function to agglomerate OTU table at the chosen phylogenetic rank
      choosen_tax <- tax_glom(physeq = agglom_sample, taxrank = tax.level)
      #extracting OTU table and taxon table
      otus <- as.data.frame(t(otu_table(choosen_tax)))
      taxon <- as.data.frame(tax_table(choosen_tax))
    }
    
    #Extracting otus and taxonomy information for the plot 
    j <- data.frame(otus[to.plot],taxon[tax.level])
    colnames(j) <- c("Counts","Tax_lvl")
  } 
  
  if (missing(max.word)) { 
    mw = 500
  } else {
    mw=max.word
  }
  
  vec = NULL
  for (i in 1:nrow(j))  {
    count <- as.vector(j[i,1])
    word <- as.vector(j[i,2])
    vec <- c(vec, rep(word,count))
  }
  
  # remove inizial string of taxonomy (i.e. p__ or c__ and so on)
  vec <- gsub("^.__","",vec)
  
  # building legend components
  fact <- as.factor(vec)
  freq_per_level <- count(fact)
  #extracting max, min and median number of observation per level
  max_legend <- max(freq_per_level$freq)
  min_legend <- min(freq_per_level$freq)
  median_legend <- median(freq_per_level$freq)
  
  #plot
  #setting background color and removing margins
  par(bg="grey90",mar=c(0,0,0,0), par(oma=c(0,0,0,0)))
  #setting up division of the screen
  close.screen(all = T)
  split.screen(matrix(c(0,0,1,1,0.97,0,1,0.97), ncol = 4))
  split.screen(matrix(c(0,0.85,0.85,1,0,0,1,1), ncol=4), screen = 2)
  #plot the cloud
  screen(3)
  wordcloud(vec, scale=c(4,0.4), max.words=mw, random.order=FALSE, random.color=FALSE, rot.per=0, use.r.layout=FALSE, colors=brewer.pal(11, "RdYlGn"))
  #add text reference to taxonomic level, sample
  screen(1)
  text(x = 0.2, y = 0.2 , paste("Taxonomic level:", tax.level), cex = 0.7, font=3)
  text(x = 0.8, y = 0.2 , paste("Sample:", to.plot), cex = 0.7, font=3)
  #plot legends. There are 2 calls to have the 50% value on the left
  screen(4)
  colorlegend(posx = c(0.35, 0.4), posy = c(0.25,0.75), col = brewer.pal(11, "RdYlGn"), zlim = c(min_legend, max_legend), zval = c(min_legend, max_legend), lab.col = "black")
  colorlegend(posx = c(0.35, 0.4), posy = c(0.25,0.75), col = brewer.pal(11, "RdYlGn"), zlim = c(min_legend, max_legend), zval = c(median_legend), lab.col = "Black", left= TRUE)
  #add reference letter "a" for dimension of charahcter in word cloud
  #text(x = 0.4, y= 0.85, "a", col = "black", cex= 4)
  #text(x = 0.4, y= 0.2, "a", col = "black", cex= 0.4)
  #Removed as appear redundant; number of observation is also mapped with color.
}
