#Return coordinates for i points, with hexagonal pack layout
place.points <- function(i) {

  #Determine 'degree' of pack i.e. how many rings
  # Example:
  #
  #      O O O
  #     . O O O 
  #    . O O O O
  #     . O O O 
  #      . O O
  #
  # (15 points) has degree = 3 (three rings, three points on each side; 
  #  outermost ring is not full, with dots marking the empty spaces)
  #
  # The number of points in each ring (excluding the innermost) proceeds
  #  as an arithmetic progression: 6, 12, 18...
  # The cumulative number of points is the sum of the progression
  #  i.e. an arithmetic series: 6, 18, 30 ...
  # So, the ith point is the sum of the arithmetic progression, plus change,
  #  where the number of terms is the number of rings less one
  # The sum of the first n terms in an arithmetic progression is given by:
  #
  #                   n(a_1 + a_n)
  #           S_n =  --------------
  #                         2
  #
  #  where a_j is the value of the jth term.
  # So we rearrange and solve for n with the quadratic equation(with a few
  #  trimmings to account for the innermost ring and 1-indexing)

  if (i == 1) {
    Degree <- 1
  } else {
    Degree <- floor((-3 + sqrt(9 + (12 * (i - 2)))) / 6) + 2
  }

  #Place the points in a rastering fashing from the middle outwards
  # Bottom row has Degree points, middle has 2Degree + 1, there are 
  # 2Degree + 1 rows in total
  points.in.row <- function(Row, Degree) {
    if (Row > Degree) {
      return(points.in.row((2 * Degree) - Row, Degree))
    } else {
      return(Degree + Row - 1)
    }
  }
  Points <- data.frame(x = numeric(), y = numeric()) 
  Row <- Degree
  Position <- 1
  for (j in 1:i) {

    x <- Position - Degree + (points.in.row(Degree, Degree) - points.in.row(Row, Degree)) / 2
    y <- Row - Degree
    Points <- rbind(Points, data.frame(x = x, y = y))

    if (Position == points.in.row(Row, Degree)) {
      if (Row == Degree) {
        Row <- Row + 1 
      } else if (Row > Degree) {
        Row <- (2 * Degree) - Row
      } else {
        Row <- (2 * Degree) - Row + 1
      }
      Position <- 1
    } else {
      Position <- Position + 1
    }
  }

  return(Points)

}

#Load OTU table and prepare
OTUTable <- read.tidy("../residences_project/picked_OTUs_winter_forward/final_otu_table_mc2_w_taxonomy.clean.tidy.txt")
OTUTable <- OTUTable[which(OTUTable$OTU %in% levels(OTUTable$OTU)[1:1000]), ]
Samples <- read.tidy("../residences_project/samples/samples.txt")[c("Sample", "Type")]
OTUTable <- merge(OTUTable, Samples, by = "Sample", all.x = TRUE)

GroupFactor <- "Type"
ColourFactor <- "Phylum"

draw.overlap.dotplot <- function(OTUTable, GroupFactor = "Sample", ColourFactor = "Phylum") {
  
  #For each group, generate list of OTUs in that group
  GroupOTUs <- dlply(OTUTable, c(GroupFactor), function(x) as.character(x$OTU), .progress = "time")

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

    #Add colour factor
    Overlap <- data.frame(OTU = Diff)
    Overlap <- unique(merge(Overlap, OTUTable[c("OTU", ColourFactor)], by = "OTU", all.x = TRUE))

    #Sort by colour factor
    Overlap <- Overlap[order(Overlap[ColourFactor]), ]

    #Return list
    return(Overlap)
  }

  #Generate all combinations for the group factor
  Groups <- levels(OTUTable[[GroupFactor]])
  Combinations <- Reduce(c, llply(1:length(Groups), function(m) combn(Groups, m = m, simplify = FALSE), .progress = "time"))

  #Get lists of overlaps for the group factor
  Overlaps <- llply(Combinations, identify.overlap, .progress = "time")
  names(Overlaps) <- unlist(llply(Combinations, function(x) paste(x, collapse = ", "), .progress = "time"))

  #Add point coordinates to overlaps
  add.coordinates <- function(Overlap) {
    Points <- place.points(nrow(Overlap))
    Overlap <- cbind(Overlap, Points)
    return(Overlap)
  }
  Overlaps <- llply(Overlaps, add.coordinates, .progress = "time")

  #Coordinates for the centres of each region in the Venn diagram
  # (A, B, C) are the three group factor levels
  # Triangle has A bottom left, B bottom right, C apex
  # Triangle starts equilateral but is warped by the relative sizes
  #  of the different groups

  L <- 20 # Base length for one side of the triangle
  M <- 100 # Base count for points in a group

  RegionsX <- numeric(length = 7)
  names(RegionsX) <- c("A", "B", "C", "AB", "AC", "BC", "ABC")
  RegionsY <- numeric(length = 7)
  names(RegionsY) <- c("A", "B", "C", "AB", "AC", "BC", "ABC")
  CirclesX <- numeric(length = 7)
  names(CirclesX) <- c("A", "B", "C")
  CirclesY <- numeric(length = 7)
  names(CirclesY) <- c("A", "B", "C")

  #A
  AWeight <- (nrow(Overlaps[[1]]) + nrow(Overlaps[[4]]) + nrow(Overlaps[[5]])) / M
  RegionsX["A"] <- 0 - ((L / 2) * AWeight)
  RegionsY["A"] <- 0 - ((L / 2) * AWeight)
  CirclesX["A"] <- 0 - ((L / 3) * AWeight)
  CirclesY["A"] <- 0 - ((L / 3) * AWeight)

  #B
  BWeight <- (nrow(Overlaps[[2]]) + nrow(Overlaps[[4]]) + nrow(Overlaps[[6]])) / M
  RegionsX["B"] <- L + ((L / 2) * BWeight)
  RegionsY["B"] <- L - ((L / 2) * BWeight)
  CirclesX["B"] <- L + ((L / 3) * BWeight)
  CirclesY["B"] <- L - ((L / 3) * BWeight)

  #C
  CWeight <- (nrow(Overlaps[[3]]) + nrow(Overlaps[[5]]) + nrow(Overlaps[[6]])) / M
  RegionsX["C"] <- L / 2
  RegionsY["C"] <- sqrt((L^2) - ((L / 2)^2)) + ((L / 2) * CWeight)
  CirclesX["C"] <- L / 3
  CirclesY["C"] <- sqrt((L^2) - ((L / 2)^2)) + ((L / 3) * CWeight)

  #AB
  RegionsX["AB"] <- mean(c(RegionsX["A"], RegionsX["B"]))
  RegionsY["AB"] <- mean(c(RegionsY["A"], RegionsY["B"]))
  
  #AC
  RegionsX["AC"] <- mean(c(RegionsX["A"], RegionsX["C"]))
  RegionsY["AC"] <- mean(c(RegionsY["A"], RegionsY["C"]))

  #BC
  RegionsX["BC"] <- mean(c(RegionsX["B"], RegionsX["C"]))
  RegionsY["BC"] <- mean(c(RegionsY["B"], RegionsY["C"]))

  #ABC
  RegionsX["ABC"] <- mean(c(RegionsX["A"], RegionsX["B"], RegionsX["C"]))
  RegionsY["ABC"] <- mean(c(RegionsY["A"], RegionsY["B"], RegionsY["C"]))

  names(RegionsX) <- names(Overlaps)
  names(RegionsY) <- names(Overlaps)

  Radii <- c(AWeight, BWeight, CWeight)
  Radii <- Radii * L * 0.6
  names(Radii) <- c("A", "B", "C")

  #Routine to transform point coordinates for each set of overlap points
  # to the correct location
  transform.points <- function(OverlapName) {

    Overlap <- Overlaps[[OverlapName]]

    x <- RegionsX[OverlapName]
    y <- RegionsY[OverlapName]

    Overlap$x <- Overlap$x + x
    Overlap$y <- Overlap$y + y
    Overlap$Overlap <- rep(OverlapName, nrow(Overlap))

    return(Overlap)
  }
  Overlaps <- ldply(names(Overlaps), transform.points, .progress = "time")

  #Generate geom_path circles for the three group factor levels
  # From http://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
  make.region.circle <- function(Region) {

    t <- seq(0, 2 * pi, length.out = 100)
    x <- CirclesX[Region] + (Radii[Region] * cos(t))
    y <- CirclesY[Region] + (Radii[Region] * sin(t))
    Circle <- data.frame(x = x, y = y)
    return(geom_path(data = Circle, mapping = aes(x = x, y = y), colour = "black"))
  }

  Plot <- ggplot(Overlaps, aes_string(x = "x", y = "y", colour = ColourFactor))
  Plot <- Plot + geom_point(size = 1)
#  Plot <- Plot + make.region.circle(1) + make.region.circle(2) + make.region.circle(3)
  Plot <- Plot + theme_minimal()
  print(Plot)

  return(Plot)

}
