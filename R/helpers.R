# Return coordinates for n points, in trapezoidal layout with p points in
# bottom row.
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

# Convert odd-r to cartesian coordinates
convert.odd.r.to.cartesian <- function(Points) {
  Points$x <- ifelse(Points$r %% 2 == 0, Points$q, Points$q + 0.5)
  Points$y <- Points$r
  Points$q <- NULL
  Points$r <- NULL
  return(Points)
}

# Given a set of points on a hex grid in odd-r coordinates, rotate 60º
# widdershins about (0, 0)
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

# Determine 'degree' of hexagonal pack i.e. how many rings
#
#      O O O
#     . O O O 
#    . O O O O
#     . O O O 
#      . O O
#
# (15 points) has degree = 3 (three rings, three points on each side; 
# outermost ring is not full, with dots marking the empty spaces)
#
# The number of points in each ring (excluding the innermost) proceeds
# as an arithmetic progression: 6, 12, 18...
#
# The cumulative number of points is the sum of the progression
# i.e. an arithmetic series: 6, 18, 30 ...
#
# So, the ith point is the sum of the arithmetic progression, plus change,
# where the number of terms is the number of rings less one
#
# The sum of the first n terms in an arithmetic progression is given by:
#
#                   n(a_1 + a_n)
#           S_n =  --------------
#                         2
#
# where a_j is the value of the jth term.
#
# So we rearrange and solve for n with the quadratic equation(with a few
# trimmings to account for the innermost ring and 1-indexing)
degree.of.hexagon <- function(i) {
  if (i <= 1) {
    Degree <- 1
  } else {
    Degree <- floor((-3 + sqrt(9 + (12 * (i - 2)))) / 6) + 2
  }
  return(Degree)
}

# Return cartesian coordinates for i points in hexagonal pack layout
# Points are placed in a spiral fashion, from centre outwards
place.points.hexagon <- function(i) {

  Degree <- degree.of.hexagon(i)

  Points <- data.frame(x = numeric(), y = numeric())
  Points <- rbind(Points, c(0, 0))
  Degree <- degree.of.hexagon(i)
  Directions <- list(c(1, 0), c(0.5, -1), c(-0.5, -1), c(-1, 0), c(-0.5, 1), c(0.5, 1))

  for (Ring in seq_len(Degree)) {
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
