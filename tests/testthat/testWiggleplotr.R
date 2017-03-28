context("Rescale introns")

test_that("joinExons produces a single GRanges object", {
  new_exons = joinExons(ncoa7_exons)
  expect_equal(length(new_exons),21) #21 elements
  expect_equal(class(new_exons)[[1]], "GRanges")
})

test_that("rescaleIntrons truncates introns correctly",{
  new_exons = joinExons(ncoa7_exons)
  rescale_introns = rescaleIntrons(ncoa7_exons, ncoa7_cdss, new_exons, 50, c(50,50))
  expect_equal(min(IRanges::end(rescale_introns$new_introns)), 50) #Correct last intron
  expect_equal(max(IRanges::end(rescale_introns$new_introns)), 8376) #Correct last intron
  expect_equal(length(rescale_introns$old_introns), length(rescale_introns$new_introns))
})

test_that("Total length of exons does not change with translation",{
  new_exons = joinExons(ncoa7_exons)
  rescale_introns = rescaleIntrons(ncoa7_exons, ncoa7_cdss, new_exons, 50, c(50,50))
  trans_exons = translateExonCoordinates(new_exons, rescale_introns$old_introns, rescale_introns$new_introns)
  expect_equal(sum(IRanges::width(trans_exons)), sum(IRanges::width(new_exons)))
})