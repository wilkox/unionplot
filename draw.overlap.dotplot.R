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

#Return coordinates for i circles, with hexagonal pack layout
place.circles <- function(i) {

  #Determine 'degree' of pack i.e. how many rings
  # Example:
  #
  #      O O O
  #     O O O O 
  #    O O O O O
  #     O O O O 
  #      O O O
  #
  # (19 circles) has degree = 3 (three rings, three circles on each side)
  #
  # The number of circles in each ring (excluding the innermost) proceeds
  #  as an arithmetic progression: 6, 12, 18...
  # The cumulative number of circles is the sum of the progression
  #  i.e. an arithmetic series: 6, 18, 30 ...
  # So, the ith circle is the sum of the arithmetic progression, plus change.
  # The sum of the first n terms in an arithmetic progression is given by:
  #
  #         n(a_1 + a_n)
  # S_n =  --------------
  #               2
  #
  # So we rearrange and solve for n with the quadratic equation.

  if (i == 1) {
    Degree = 1  
  } else {
    Degree <- floor((-3 + sqrt(9 + (12 * (i - 2)))) / 6) + 2
  }

  return(Degree)

}


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
