#Return coordinates for n points, in trapezoidal layout with p points
# in bottom row
#Indexing is hex, odd-r
place.points.trapezoid <- function(n, p) {

  Points <- data.frame(q = numeric(), r = numeric()) 
  Row <- 1
  Position <- 1

  for (j in 1:n) {

    q <- Position - floor(Row / 2) - 1
    r <- Row - 1

    Points <- rbind(Points, data.frame(r = r, q = q))

    if (Position == (p + Row - 1)) {
      Row <- Row + 1
      Position <- 1
    } else {
      Position <- Position + 1
    }
    
  }

  return(Points)
  
}

#Convert a hex grid in odd-r indexing to cartesian coordinates
convert.odd.r.to.cartesian <- function(Points) {
  Points$x <- ifelse(Points$r %% 2 == 0, Points$q, Points$q + 0.5)
  Points$y <- Points$r
  Points$q <- NULL
  Points$r <- NULL
  return(Points)
}

#Rotate a hex grid in odd-r indexing
rotate.odd.r <- function(Points, Rotations = 1) {

  if (!Rotations == 1) {
    Points <- rotate.odd.r(Points, Rotations = Rotations - 1)
  }
  
  #Convert to cubic
  x <- Points$q - (Points$r - (Points$r %% 2)) / 2
  z <- Points$r
  y <- -x-z

  #Rotate
  xx <- -z
  yy <- -x
  zz <- -y

  #Convert back to odd-r
  Points$q <- xx + (zz - (zz %% 2)) / 2
  Points$r <- zz

  return(Points)
    
}

#Return cartesian coordinates for i points, with hexagonal pack layout
place.points.hexagon <- function(i) {

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

  Degree <- degree(i)

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

  #We need the degree of the centre hexagon to know where to place the
  # surrounding trapazoids
  if (nrow(Overlaps[[7]]) == 1) {
    Degree <- 1
  } else {
    Degree <- floor((-3 + sqrt(9 + (12 * (nrow(Overlaps[[7]]) - 2)))) / 6) + 2
  }

  #Routine to generate a trapezoid
  make.trapezoid <- function(OverlapIndex, BaseRow, Rotation, Offset = c(0,0)) {
    Trap <- place.points.trapezoid(nrow(Overlaps[[OverlapIndex]]), BaseRow)
    Trap <- rotate.odd.r(Trap, Rotation)
    Trap <- convert.odd.r.to.cartesian(Trap)
    Trap$x <- Trap$x + Offset[1]
    Trap$y <- Trap$y + Offset[2]
    Trap <- cbind(Overlaps[[OverlapIndex]], Trap)
    return(Trap)
  }

  #Place the centre hexagon
  Points <- cbind(Overlaps[[7]], place.points.hexagon(nrow(Overlaps[[7]])))

  #Place the top trapezoid
  Points <- rbind(Points, make.trapezoid(1, Degree, 6, c(0.5 - (0.5 * Degree), Degree + 1)))

  #Place the bottom right trapezoid
  Points <- rbind(Points, make.trapezoid(2, Degree, 4, c(Degree + 0.5, -1)))
  
  #Place the bottom left trapezoid
  Points <- rbind(Points, make.trapezoid(3, Degree, 2, c(-(Degree / 2) - 1, -Degree)))

  #Place the top right trapezoid
  Points <- rbind(Points, make.trapezoid(4, Degree, 5, c((Degree / 2) + 1, Degree)))

  #Place the bottom trapezoid
  Points <- rbind(Points, make.trapezoid(6, Degree, 3, c(- 0.5 + (Degree / 2), - Degree - 1)))

  #Place the top left trapezoid
  Points <- rbind(Points, make.trapezoid(5, Degree, 1, c(- Degree - 1, 1)))

  #Draw dividing lines
  #Main hex
  Hex <- data.frame(x = c(-Degree / 2, Degree / 2, Degree, Degree / 2, -Degree / 2, -Degree, -Degree / 2), y = c(Degree, Degree, 0, -Degree, -Degree, 0, Degree))

  Plot <- ggplot(Points, aes(x = x, y = y))
  Plot <- Plot + geom_point(aes(colour = Phylum))
  Plot <- Plot + geom_path(data = Hex)
  Plot <- Plot + theme_classic()
  Plot

}
