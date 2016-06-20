#' @title Place points in a trapezoid
#' @description 
#' Return coordinates for n points, in trapezoidal layout with p points
#' in bottom row.
#' 
#' Coordinate system is hex, odd-r.
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

#' @title Convert odd-r to cartesian coordinates
convert.odd.r.to.cartesian <- function(Points) {
  Points$x <- ifelse(Points$r %% 2 == 0, Points$q, Points$q + 0.5)
  Points$y <- Points$r
  Points$q <- NULL
  Points$r <- NULL
  return(Points)
}

#' @title Rotate odd-r points
#'
#' @description
#' Given a set of points on a hex grid in odd-r coordinates,
#' rotate 60º widdershins about (0, 0)
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

#' @title Determine 'degree' of hexagonal pack i.e. how many rings
#' @description
#'
#'      O O O
#'     . O O O 
#'    . O O O O
#'     . O O O 
#'      . O O
#'
#' (15 points) has degree = 3 (three rings, three points on each side; 
#' outermost ring is not full, with dots marking the empty spaces)
#'
#' The number of points in each ring (excluding the innermost) proceeds
#' as an arithmetic progression: 6, 12, 18...
#'
#' The cumulative number of points is the sum of the progression
#' i.e. an arithmetic series: 6, 18, 30 ...
#'
#' So, the ith point is the sum of the arithmetic progression, plus change,
#' where the number of terms is the number of rings less one
#'
#' The sum of the first n terms in an arithmetic progression is given by:
#'
#'                   n(a_1 + a_n)
#'           S_n =  --------------
#'                         2
#'
#' where a_j is the value of the jth term.
#'
#' So we rearrange and solve for n with the quadratic equation(with a few
#' trimmings to account for the innermost ring and 1-indexing)
degree.of.hexagon <- function(i) {
  if (i == 1) {
    Degree <- 1
  } else {
    Degree <- floor((-3 + sqrt(9 + (12 * (i - 2)))) / 6) + 2
  }
  return(Degree)
}

#' @title Return cartesian coordinates for i points in hexagonal pack layout
#' @description
#' Points are placed in a spiral fashion, from centre outwards
place.points.hexagon <- function(i) {

  Degree <- degree.of.hexagon(i)

  Points <- data.frame(x = numeric(), y = numeric())
  Points <- rbind(Points, c(0, 0))
  Degree <- degree.of.hexagon(i)
  Directions <- list(c(1, 0), c(0.5, -1), c(-0.5, -1), c(-1, 0), c(-0.5, 1), c(0.5, 1))

  for (Ring in 1:Degree) {
    Point <- Directions[[5]] * Ring
    for (Side in 1:6) {
      for (SidePos in 1:Ring) {
        Points <- rbind(Points, Point)
        Point <- Point + Directions[[Side]]
      }
    }
    Ring <- Ring + 1
  }
  names(Points) <- c("x", "y")

  #Select only requested number of points, and reverse order
  # for "outside-in" effect
  Points <- Points[i:1, ]

  return(Points)

}

#' @title Draw a union plot
#' @export
#'
#' @description
#'
#' A union plot is like a Venn diagram, in that it shows the overlap between three
#' groups. Unlike a Venn diagram, the number of overlapping units (in this case, OTUs)
#' is not indicated by a number but by a point drawn for each OTU, which can then be
#' coloured by a factor of interest e.g. Phylum. Also unlike a Venn diagram, the regions
#' are not indicated by overlapping circles but by hexagons and trapezoids; the points
#' are arranged on a hexagonal grid. This allows the regions to be resized according
#' to the number of shared OTUs.
#' 
#' There's probably already a name for this kind of plot, but I couldn't find it after
#' some cursory googling so just named it "union plot". Email me if you know the actual
#' name and I'll fix it.
#' 
#' @param OTUTable an OTU table in tidy format, i.e. a data frame with at least an "OTU" 
#' column as well as columns for the group and colour factors.
#' @param GroupFactor a string corresponding to the name of a factor column in OTUTable, 
#' which determines the three groups that will form the three cardinal regions in the plot.
#' Defaults to "Sample".
#' @param ColourFactor a string corresponding to the name of a factor column in OTUTable,
#' which determines the colours of the points. ColourFactor is mandatory; if you don't want
#' the points to be coloured, you should be drawing a Venn diagram instead. Defaults to
#' "Phylum".
#' @param Pointsize (optional) number to be passed to geom_point() as the point size. Defaults to 1.
#' @param Collapse (optional) integer. If > 1, each point in the plot will not represent
#' a single OTU but this number of OTUs, with the true number rounded up to the nearest
#' multiple of this value (there are no fractional points).
#' @param silent suppress messages
#'
#' @return
#' Returns a ggplot2 grob containing the union plot, which can then be viewed (e.g. with
#' print()) or saved (e.g. with ggsave()) at leisure. Note that the placement of text labels
#' in the plot isn't the greatest; you'll probably need to go in and move them using
#' an image editing suite.
#'
#' @references Amit Patel's page on hexagonal grids was invaluable in creating this
#' function, and many of the routines in here are derived from algorithms found there
#' \url{http://www.redblobgames.com/grids/hexagons}
draw.union.plot <- function(OTUTable, GroupFactor = "Sample", ColourFactor = "Phylum", Pointsize = 1, Collapse = 1, silent = FALSE) {

  Progress <- ifelse(silent, "none", "time")
  
  #For each group, generate list of OTUs in that group
  if (! silent) {
    message(paste0("Generate list of OTUs in each ", GroupFactor, "..."))
  }
  GroupOTUs <- dlply(OTUTable, c(GroupFactor), function(x) as.character(x$OTU), .progress = Progress)

  #Routine to determine intersection of >2 vectors
  deep.intersect <- function(List) {
    if (length(List) == 2) {
      return(intersect(unlist(List[[1]]), unlist(List[[2]])))
    } else {
      return(intersect(unlist(List[[1]]), deep.intersect(List[c(2:length(List))])))
    }
  }

  #Routine to identify union between list of groups
  identify.union <- function(Groups) {

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
  if (! silent) {
    message(paste0("Generate all combinations of ", GroupFactor, "s..."))
  }
  Groups <- levels(OTUTable[[GroupFactor]])
  Combinations <- Reduce(c, llply(1:length(Groups), function(m) combn(Groups, m = m, simplify = FALSE), .progress = Progress))

  #Get lists of unions for the group factor
  if (! silent) {
    message(paste0("Generate lists of OTUs shared between each combination of ", GroupFactor, "s..."))
  }
  Overlaps <- llply(Combinations, identify.union, .progress = Progress)
  names(Overlaps) <- unlist(llply(Combinations, function(x) paste(x, collapse = ", "), .progress = Progress))

  #If a collapse was requested, subsample as appropriate
  OriginalOverlaps <- Overlaps
  if (Collapse > 1) {
    subsample.group <- function(Group) {
      m <- ceiling(nrow(Group) / Collapse)
      return(Group[sample(1:nrow(Group), m), ])
    }
    Overlaps <- llply(Overlaps, function(Overlap) ddply(Overlap, c(ColourFactor), subsample.group))
  }

  #We need the degree of the centre hexagon to know where to place the
  # surrounding trapazoids
  if (! silent) {
    message("Determine degree of centre hexagon...")
  }
  Degree <- degree.of.hexagon(nrow(Overlaps[[7]]))

  #Routine to generate a trapezoid
  make.trapezoid <- function(OverlapIndex, BaseRow, Rotation, Offset = c(0,0)) {

    # Return nothing if there are no points in this overlap
    if (nrow(Overlaps[[OverlapIndex]]) == 0) {
      return()  
    }

    Trap <- place.points.trapezoid(nrow(Overlaps[[OverlapIndex]]), BaseRow)
    Trap <- rotate.odd.r(Trap, Rotation)
    Trap <- convert.odd.r.to.cartesian(Trap)
    Trap$x <- Trap$x + Offset[1]
    Trap$y <- Trap$y + Offset[2]
    Trap <- cbind(Overlaps[[OverlapIndex]], Trap)
    return(Trap)
  }

  if (! silent) {
    message("Generate hexagons and trapezoids...")
  }

  #Place the centre hexagon
  if (! silent) {
    message(paste0(names(Overlaps)[7], "..."))
  }
  Points <- cbind(Overlaps[[7]], place.points.hexagon(nrow(Overlaps[[7]])))

  #Place the top trapezoid
  if (! silent) {
    message(paste0(names(Overlaps)[1], "..."))
  }
  Points <- rbind(Points, make.trapezoid(1, Degree, 6, c(0.5 - (0.5 * Degree), Degree + 1)))

  #Place the bottom right trapezoid
  if (! silent) {
    message(paste0(names(Overlaps)[2], "..."))
  }
  Points <- rbind(Points, make.trapezoid(2, Degree, 4, c(Degree + 0.5, -1)))
  
  #Place the bottom left trapezoid
  if (! silent) {
    message(paste0(names(Overlaps)[3], "..."))
  }
  Points <- rbind(Points, make.trapezoid(3, Degree, 2, c(-(Degree / 2) - 1, -Degree)))

  #Place the top right trapezoid
  if (! silent) {
    message(paste0(names(Overlaps)[4], "..."))
  }
  Points <- rbind(Points, make.trapezoid(4, Degree, 5, c((Degree / 2) + 1, Degree)))

  #Place the bottom trapezoid
  if (! silent) {
    message(paste0(names(Overlaps)[6], "..."))
  }
  Points <- rbind(Points, make.trapezoid(6, Degree, 3, c(- 0.5 + (Degree / 2), - Degree - 1)))

  #Place the top left trapezoid
  if (! silent) {
    message(paste0(names(Overlaps)[5], "..."))
  }
  Points <- rbind(Points, make.trapezoid(5, Degree, 1, c(- Degree - 0.5, 1)))

  if (! silent) {
    message("Drawing lines and labels...")
  }

  #Draw dividing lines
  #Main hex
  Hex <- data.frame(x = c(-Degree / 2, Degree / 2, Degree, Degree / 2, -Degree / 2, -Degree, -Degree / 2), y = c(Degree, Degree, 0, -Degree, -Degree, 0, Degree))

  #To draw the lines between trapezoids, we need to know how many rows are in each trapezoid
  # Blah blah arithmetic progression quadratic solve for blah
  # Where n is number of points in trapezoid, p is number in base row
  # We use ceiling to account for unfilled rows
  trapezoid.degree <- function(n, p) {
    return(ceiling((1 - (2 * p) + sqrt((((2 * p) - 1) ^ 2) + (8 * n))) / 2))
  }
  trapezoid.degree(nrow(Overlaps[[1]]), Degree)

  #Routine to add dividing line
  add.dividing.line <- function(Start = c(0,0), OverlapIndices = c(0,0), xDir = 1, yDir = 1) {
    x <- Start[1]
    y <- Start[2]
    MaxTrapDegree <- max(trapezoid.degree(nrow(Overlaps[[OverlapIndices[1]]]), Degree), trapezoid.degree(nrow(Overlaps[[OverlapIndices[2]]]), Degree))
    xDist <- ifelse(yDir == 0, 1, 0.5)
    xend <- x + (xDir * MaxTrapDegree * xDist) + (xDir * xDist)
    yend <- y + (yDir * (MaxTrapDegree + 1))
    return(data.frame(x = x, xend = xend, y = y, yend = yend))
  }

  DividingLines <- data.frame(x = numeric(), xend = numeric(), y = numeric(), yend = numeric)

  #Add top left dividing line
  DividingLines <- rbind(DividingLines, add.dividing.line(c(-Degree / 2, Degree), c(1, 5), -1, 1))

  #Add top right dividing line
  DividingLines <- rbind(DividingLines, add.dividing.line(c(Degree / 2, Degree), c(1, 4), 1, 1))

  #Add mid right dividing line
  DividingLines <- rbind(DividingLines, add.dividing.line(c(Degree, 0), c(2, 4), 1, 0))

  #Add bottom right dividing line
  DividingLines <- rbind(DividingLines, add.dividing.line(c(Degree / 2, -Degree), c(2, 6), 1, -1))

  #Add bottom left dividing line
  DividingLines <- rbind(DividingLines, add.dividing.line(c(-Degree / 2, -Degree), c(3, 6), -1, -1))

  #Add mid left dividing line
  DividingLines <- rbind(DividingLines, add.dividing.line(c(-Degree, 0), c(3, 5), -1, 0))

  #Routine to add text label
  add.label <- function(OverlapIndex, xdir, ydir, hjust, vjust) {
    x <- (1 + Degree) * xdir
    y <- Degree * ydir
    OTUs <- nrow(OriginalOverlaps[[OverlapIndex]])
    Name <- names(Overlaps)[OverlapIndex]
    label <- paste0(Name, " (", OTUs, ")")
    return(data.frame(label = label, x = x, y = y, hjust = hjust, vjust = vjust))
  }

  ##Text labels for each segment
  Labels <- data.frame(label = character(), x = numeric(), y = numeric(), hjust = numeric(), vjust = numeric())

  #Add text label for centre
  Labels <- rbind(Labels, add.label(7, 0, 1, 0.5, 1))

  #Add text label for top segment
  Labels <- rbind(Labels, add.label(1, 0, 1, 0.5, 0))

  #Add text label for top right segment
  Labels <- rbind(Labels, add.label(4, 1, 0, 0, 0))

  #Add text label for bottom right segment
  Labels <- rbind(Labels, add.label(2, 1, 0, 0, 1))

  #Add text label for bottom segment
  Labels <- rbind(Labels, add.label(6, 0, -1, 0.5, 1))

  #Add text label for bottom left segment
  Labels <- rbind(Labels, add.label(3, -1, 0, 1, 1))

  #Add text label for top left segment
  Labels <- rbind(Labels, add.label(5, -1, 0, 1, 0))

  #If a collapse was requested, add an annotation
  if (Collapse > 1) {
    ColourFactorDescription <- paste0(ColourFactor, " (1-", Collapse, " OTUs)")
  } else {
    ColourFactorDescription <- ColourFactor
  }

  Plot <- ggplot(Points, aes(x = x, y = y))
  Plot <- Plot + geom_point(aes_string(colour = ColourFactor), size = Pointsize)
  Plot <- Plot + geom_path(data = Hex)
  Plot <- Plot + geom_segment(data = DividingLines, aes(x = x, y = y, xend = xend, yend = yend))
  Plot <- Plot + geom_text(data = Labels, aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust), family = "Helvetica", fontface = "bold", size = 3)
  if (length(levels(Points[ColourFactor])) <= 8) {
    Plot <- Plot + scale_colour_brewer(palette = "Set2", name = ColourFactorDescription)
  } else {
    Plot <- Plot + scale_colour_discrete(name = ColourFactorDescription)
  }
  Plot <- Plot + theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.key = element_blank()
  
  )
  Plot <- Plot + coord_fixed()

  return(Plot)

}
