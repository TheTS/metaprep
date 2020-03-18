test_that("Testing normal workflow", {

  data(velocity, package = 'metaprep')

  x <- pool_studies(data = velocity,
                    n.e = n.e,
                    mean.e = mean.e,
                    sd.e = sd.e,
                    n.c = n.c,
                    mean.c = mean.c,
                    sd.c = sd.c,
                    type = type,
                    studlab = study)


  expect_equal(nrow(x), 21)
  expect_equal(ncol(x), 11)
  expect_equal(sum(x$se_g), 7.218245)

})
