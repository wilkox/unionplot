#Load OTU table and prepare
OTUTable <- read.tidy("../residences_project/picked_OTUs_winter_forward/final_otu_table_mc2_w_taxonomy.clean.tidy.txt")
OTUTable <- OTUTable[which(OTUTable$OTU %in% levels(OTUTable$OTU)[1:1000]), ]
Samples <- read.tidy("../residences_project/samples/samples.txt")[c("Sample", "Type")]
OTUTable <- merge(OTUTable, Samples, by = "Sample", all.x = TRUE)

GroupFactor <- "Type"
ColourFactor <- "Phylum"
SizeFactor <- "RelativeAbundance"

Data <- OTUTable[which(OTUTable$Sample == levels(OTUTable$Sample)[1]), ]
Centre <- c(50,50)
Max.radius <- 20

#Return the capacity of ring r
ring.capacity <- function(r) {
  if (r == 0) {
    return(1) 
  } else if (r == 1) {
    return(6) 
  } else {
    return(2 * ring.capacity(r - 1))  
  }
}

#Return the cartesian coordinates for datum n
place.datum <- function(n) {

  #Get ring and index within ring
  Ring <- floor(log(((n + 4) / 6), 2)) + 1
  Index <- n - abs(((2 ^ (Ring - 1) - 1) * 6) + 2)

  #Convert to cartesian coordinates
  r <- Ring
  phi <- (2 * pi * Index) / ring.capacity(Ring)
  x <- r * cos(phi)
  y <- r * sin(phi)

  return(data.frame(x = x, y = y))
}

Data <- ldply(1:500, place.datum)
Data$n <- factor(1:500)
Plot <- ggplot(Data, aes(x = x, y = y, colour = n)) + geom_point()
Plot

draw.overlap.dotplot <- function(OTUTable, GroupFactor = "Sample", ColourFactor = "Phylum", SizeFactor = "RelativeAbundance") {
  
  #For each group, generate list of OTUs in that group
  GroupOTUs <- dlply(OTUTable, c(GroupFactor), function(x) as.character(x$OTU))

  #Routine to determine intersection of >2 vectors
  deep.intersect <- function(List) {
    if (length(List) == 2) {
      return(intersect(unlist(List[[1]]), unlist(List[[2]])))
    } else {
      return(intersect(unlist(List[[1]]), deep.intersect(List[c(2:length(List))])))
    }
  }

  #Routine to identify overlap between list of groups
  identify.overlap <- function(Groups) {

    #Get set intersection for groups, unless there is only one
    if (length(Groups) == 1) {
      Intersection <- unlist(GroupOTUs[Groups[1]])
    } else {
      Intersection <- unlist(deep.intersect(GroupOTUs[Groups]))
    }

    #Get set difference with other groups, unless all groups are in
    if (length(Groups) == length(GroupOTUs)) {
      Diff <- Intersection
    } else {
      Diff <- setdiff(Intersection, unlist(GroupOTUs[-(which(names(GroupOTUs) %in% Groups))]))
    }

    #Return list
    return(Diff)
  }

  #Generate all combinations for the group factor
  Groups <- levels(OTUTable[[GroupFactor]])
  Combinations <- Reduce(c, llply(1:length(Groups), function(m) combn(Groups, m = m, simplify = FALSE)))

  #Get lists of overlaps for the group factor
  Overlaps <- llply(Combinations, identify.overlap, .progress = "time")

}
