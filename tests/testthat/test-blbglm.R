test_that("Functions return same results", {
  fit <- blbglm(Species ~ Sepal.Length * Sepal.Width, data = iris, family = binomial(), m = 3, B = 100)
  expect_identical(coef(fit), coef.blbglm(fit))
  expect_identical(confint(fit), confint.blbglm(fit))
  expect_identical(sigma(fit), sigma.blbglm(fit))
})

test_that("Parallelization runs faster", {
  expect_gt(
    system.time(blbglm(Species ~ Sepal.Length * Sepal.Width, data = iris, family = binomial(), m = 3, B = 100))[1],
    system.time(blbglm(Species ~ Sepal.Length * Sepal.Width, data = iris, family = binomial(), m = 3, B = 100, nthreads = 4))[1]
  )
})
