#' @title Draw a union plot
#' @export
#'
#' @description
#'
#' A union plot is like a Venn diagram, in that it shows the union between three
#' groups. Unlike a Venn diagram, the number of shared units (in this case, OTUs)
#' is not indicated by a number but by a point drawn for each OTU, which can then be
#' coloured by a factor of interest e.g. Phylum. Also unlike a Venn diagram, the regions
#' are not indicated by shared circles but by hexagons and trapezoids; the points
#' are arranged on a hexagonal grid. This allows the regions to be resized according
#' to the number of shared OTUs.
#' 
#' There's probably already a name for this kind of plot, but I couldn't find it after
#' some cursory googling so just named it "union plot". Email me if you know the actual
#' name and I'll fix it.
#' 
#' @param OTUTable an OTU table in tidy format, i.e. a data frame with at least an "OTU" 
#' column as well as columns for the group and colour factors.
#' @param group a string corresponding to the name of a factor column in OTUTable, 
#' which determines the three groups that will form the three cardinal regions in the plot.
#' Defaults to "Sample".
#' @param colour a string corresponding to the name of a factor column in OTUTable,
#' which determines the colours of the points. colour is mandatory; if you don't want
#' the points to be coloured, you should be drawing a Venn diagram instead. Defaults to
#' "Phylum".
#' @param point_size (optional) number to be passed to geom_point() as the point size. Defaults to 1.
#' @param collapse (optional) integer. If > 1, each point in the plot will not represent
#' a single OTU but this number of OTUs, with the true number rounded up to the nearest
#' multiple of this value (there are no fractional points).
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
unionplot <- function(OTUTable, group = "Sample", colour = "Phylum", point_size = 1, collapse = 1) {

  # For each group, generate list of OTUs in that group
  group_OTUs <- split(OTUTable, OTUTable[[group]])
  group_OTUs <- purrr::map(group_OTUs, ~ unique(.x[["OTU"]]))

  # Generate all seven combinations of groups
  groups <- as.character(unique(OTUTable[[group]]))
  combinations <- purrr::map(
    seq_len(length(groups)),
    ~ combn(groups, m = ., simplify = FALSE)
  ) 
  combinations <- Reduce(c, combinations)

  # Get lists of unions for the group factor
  get_intersect <- function(x) {
    purrr::reduce(
      purrr::map(x, ~ OTUTable[which(OTUTable[[group]] %in% .), ][["OTU"]]),
      intersect
    )
  }
  unions <- purrr::map(combinations, get_intersect)
  names(unions) <- purrr::map(combinations, ~ paste(., collapse = ", "))

  # Generate table of individual OTUs for each union
  unions <- purrr::map(unions, ~ OTUTable[OTUTable[["OTU"]] %in% ., ])
  unions <- purrr::map(unions, ~ .[, c("OTU", colour, "Count")])

  # Collapse counts if a collapse was requested
  original_unions <- unions
  if (collapse > 1) {
    unions <- purrr::map(
      unions,
      ~ dplyr::mutate(., Count = ceiling(Count / collapse))
    )
  }

  # Convert OTU table into one-row-per-point
  unions <- purrr::map(unions, ~ dplyr::slice(., rep(seq_len(nrow(.)), times = .[["Count"]])))
  unions <- purrr::map(unions, ~ .[, c("OTU", colour)])

  # Sort by colour
  unions <- purrr::map(unions, ~ dplyr::arrange_(., colour))

  # We need the degree of the centre hexagon to know where to place the
  # surrounding trapazoids
  centre_hex_degree <- degree.of.hexagon(nrow(unions[[7]]))

  # Routine to generate a trapezoid
  make.trapezoid <- function(union, base_row_n, rotation, offset = c(0,0)) {

    # Return nothing if there are no points in this union
    if (nrow(union) == 0) {
      message("Skipping, no points in this union")
      return(union)  
    }

    # Place points in trapezoid
    trapezoid <- place.points.trapezoid(nrow(union), base_row_n)

    # Rotate to fit desired rotation
    trapezoid <- rotate.odd.r(trapezoid, rotation)

    # Convert coordinates from odd-r to Cartesian
    trapezoid <- convert.odd.r.to.cartesian(trapezoid)

    # Adjust for offset
    trapezoid$x <- trapezoid$x + offset[1]
    trapezoid$y <- trapezoid$y + offset[2]

    # Bind to union data and return
    trapezoid <- dplyr::bind_cols(union, trapezoid)
    trapezoid
  }

  # Place the centre hexagon
  if (nrow(unions[[7]]) > 0) {
    point_coords <- dplyr::bind_cols(unions[[7]], place.points.hexagon(nrow(unions[[7]])))
  } else {
    point_coords <- dplyr::bind_cols(unions[[7]], data.frame(x = double(), y = double()))
  }

  # Place the top trapezoid
  point_coords <- dplyr::bind_rows(
    point_coords,
    make.trapezoid(
      unions[[1]],
      centre_hex_degree,
      6,
      c(0.5 - (0.5 * centre_hex_degree), centre_hex_degree + 1)
    )
  )

  # Place the bottom right trapezoid
  point_coords <- dplyr::bind_rows(
    point_coords,
    make.trapezoid(unions[[2]], centre_hex_degree, 4, c(centre_hex_degree + 0.5, -1))
  )
  
  # Place the bottom left trapezoid
  point_coords <- dplyr::bind_rows(
    point_coords,
    make.trapezoid(
      unions[[3]],
      centre_hex_degree,
      2,
      c(-(centre_hex_degree / 2) - 1, -centre_hex_degree)
    )
  )

  # Place the top right trapezoid
  point_coords <- dplyr::bind_rows(
    point_coords,
    make.trapezoid(
      unions[[4]],
      centre_hex_degree,
      5,
      c((centre_hex_degree / 2) + 1, centre_hex_degree)
    )
  )

  # Place the bottom trapezoid
  point_coords <- dplyr::bind_rows(
    point_coords,
    make.trapezoid(
      unions[[6]],
      centre_hex_degree,
      3,
      c(- 0.5 + (centre_hex_degree / 2), - centre_hex_degree - 1)
    )
  )

  # Place the top left trapezoid
  point_coords <- dplyr::bind_rows(
    point_coords,
    make.trapezoid(unions[[5]], centre_hex_degree, 1, c(- centre_hex_degree - 0.5, 1))
  )

  # Draw dividing lines
  # Main hex
  hex_lines <- data.frame(
    x = c(-centre_hex_degree / 2, centre_hex_degree / 2, centre_hex_degree, centre_hex_degree / 2, -centre_hex_degree / 2, -centre_hex_degree, -centre_hex_degree / 2),
    y = c(centre_hex_degree, centre_hex_degree, 0, -centre_hex_degree, -centre_hex_degree, 0, centre_hex_degree)
  )

  # To draw the lines between trapezoids, we need to know how many rows are in each trapezoid
  # Blah blah arithmetic progression quadratic solve for blah
  # Where n is number of points in trapezoid, p is number in base row
  # We use ceiling to account for unfilled rows
  trapezoid_degree <- function(n, p) {
    return(ceiling((1 - (2 * p) + sqrt((((2 * p) - 1) ^ 2) + (8 * n))) / 2))
  }

  # Routine to create coordinates for inter-trapezoid dividing lines
  dividing_line <- function(start = c(0, 0), overlap_indices = c(0, 0), x_dir = 1, y_dir = 1) {
    x <- start[1]
    y <- start[2]
    max_trap_degree <- max(
      trapezoid_degree(nrow(unions[[overlap_indices[1]]]), centre_hex_degree),
      trapezoid_degree(nrow(unions[[overlap_indices[2]]]), centre_hex_degree)
    )
    x_dist <- ifelse(y_dir == 0, 1, 0.5)
    xend <- x + (x_dir * max_trap_degree * x_dist) + (x_dir * x_dist)
    yend <- y + (y_dir * (max_trap_degree + 1))
    return(data.frame(x = x, xend = xend, y = y, yend = yend))
  }

  dividing_lines <- data.frame(
    x = double(),
    xend = double(),
    y = double(),
    yend = double()
  )

  # Add top left dividing line
  dividing_lines <- dplyr::bind_rows(
    dividing_lines,
    dividing_line(c(-centre_hex_degree / 2, centre_hex_degree), c(1, 5), -1, 1)
  )

  # Add top right dividing line
  dividing_lines <- dplyr::bind_rows(
    dividing_lines,
    dividing_line(c(centre_hex_degree / 2, centre_hex_degree), c(1, 4), 1, 1)
  )

  # Add mid right dividing line
  dividing_lines <- dplyr::bind_rows(
    dividing_lines,
    dividing_line(c(centre_hex_degree, 0), c(2, 4), 1, 0)
  )

  # Add bottom right dividing line
  dividing_lines <- dplyr::bind_rows(
    dividing_lines,
    dividing_line(c(centre_hex_degree / 2, -centre_hex_degree), c(2, 6), 1, -1)
  )

  # Add bottom left dividing line
  dividing_lines <- dplyr::bind_rows(
    dividing_lines,
    dividing_line(c(-centre_hex_degree / 2, -centre_hex_degree), c(3, 6), -1, -1)
  )

  # Add mid left dividing line
  dividing_lines <- dplyr::bind_rows(
    dividing_lines,
    dividing_line(c(-centre_hex_degree, 0), c(3, 5), -1, 0)
  )

  # Routine to add text label
  add_label <- function(union_i, xdir, ydir, hjust, vjust) {
    x <- (1 + centre_hex_degree) * xdir
    y <- centre_hex_degree * ydir
    OTUs <- nrow(original_unions[[union_i]])
    name <- names(unions)[union_i]
    label <- stringr::str_c(name, " (", OTUs, ")")
    return(data.frame(
      label = label,
      x = x,
      y = y,
      hjust = hjust,
      vjust = vjust,
      stringsAsFactors = FALSE
    ))
  }

  # Text labels for each segment
  segment_labels <- data.frame(
    label = character(),
    x = double(),
    y = double(),
    hjust = double(),
    vjust = double(),
    stringsAsFactors = FALSE
  )

  # Add text label for centre
  segment_labels <- dplyr::bind_rows(segment_labels, add_label(7, 0, 1, 0.5, 1))

  # Add text label for top segment
  segment_labels <- dplyr::bind_rows(segment_labels, add_label(1, 0, 1, 0.5, 0))

  # Add text label for top right segment
  segment_labels <- dplyr::bind_rows(segment_labels, add_label(4, 1, 0, 0, 0))

  # Add text label for bottom right segment
  segment_labels <- dplyr::bind_rows(segment_labels, add_label(2, 1, 0, 0, 1))

  # Add text label for bottom segment
  segment_labels <- dplyr::bind_rows(segment_labels, add_label(6, 0, -1, 0.5, 1))

  # Add text label for bottom left segment
  segment_labels <- dplyr::bind_rows(segment_labels, add_label(3, -1, 0, 1, 1))

  # Add text label for top left segment
  segment_labels <- dplyr::bind_rows(segment_labels, add_label(5, -1, 0, 1, 0))

  # If a collapse was requested, add an annotation
  if (collapse > 1) {
    colour_legend_name <- stringr::str_c(colour, " (1-", collapse, " OTUs)")
  } else {
    colour_legend_name <- colour
  }

  # Fix for fussy R CMD check
  x <- y <- xend <- yend <- label <- hjust <- vjust <- 0

  # Build basic plot
  p <- ggplot2::ggplot(point_coords, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes_string(colour = colour), size = point_size) +
    ggplot2::geom_path(data = hex_lines) +
    ggplot2::geom_segment(
      data = dividing_lines,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    ggplot2::geom_text(
      data = segment_labels,
      ggplot2::aes(label = label, x = x, y = y, hjust = hjust, vjust = vjust),
      family = "Helvetica",
      fontface = "bold",
      size = 3
    )

  # Set colour scale
  if (nrow(unique(point_coords[colour])) <= 8) {
    p <- p + ggplot2::scale_colour_brewer(palette = "Set2", name = colour_legend_name)
  } else {
    p <- p + ggplot2::scale_colour_discrete(name = colour_legend_name)
  }

  # Set theme
  p <- p + ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank()
  
  )

  # Set coordinate system
  p <- p + ggplot2::coord_fixed()

  # Return plot
  p

}
