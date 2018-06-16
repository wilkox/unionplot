context("unionplot")

test_that("unionplot works", {
  data(OTUTable)
  expect_silent({
    unionplot(OTUTable, group = "Sample", colour = "Family")
    })
})
