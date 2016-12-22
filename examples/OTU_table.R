# Prepare example OTU table
library(stringr)
library(devtools)
library(dplyr)

OTUs <- str_c("OTU", 1:1500)
Families <- c(
  rep("Entreebacteriaceae", 500),
  rep("Mainbacteriaceae", 500),
  rep("Dessertbacteriaceae", 500)
)
OTUFamily <- data.frame(OTU = OTUs, Family = Families)
Samples <- c("Restaurant A", "Restaurant B", "Restaurant C")
OTUTable <- data.frame(
    Sample = rep(Samples, 500)
  ) %>%
  arrange(Sample) %>%
  mutate(OTU = c(
    sample(
      OTUs,
      500,
      prob = c(rep(1.5, 500), rep(0.5, 1000)),
      replace = TRUE
    ),
    sample(
      OTUs,
      500,
      prob = c(rep(0.5, 500), rep(1.5, 500), rep(0.5, 500)),
      replace = TRUE
    ),
    sample(
      OTUs,
      500,
      prob = c(rep(0.5, 1000), rep(1.5, 500)),
      replace = TRUE
    )
  )) %>%
  left_join(OTUFamily) %>%
  as.tbl %>%
  group_by(Sample, OTU, Family) %>%
  summarise(Count = n()) %>%
  ungroup %>%
  as.data.frame
devtools::use_data(OTUTable, overwrite = TRUE)
