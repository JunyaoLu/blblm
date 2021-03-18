test_that("Functions return same results", {
  fit <- blbglm(Species ~ Sepal.Length * Sepal.Width, data = iris, family = binomial(), m = 3, B = 100)
  expect_identical(coef(fit), coef.blbglm(fit))
  expect_identical(confint(fit, c("wt", "hp")), confint.blbglm(fit, c("wt", "hp")))
  expect_identical(sigma(fit), sigma.blbglm(fit))
  expect_identical(
    predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170))),
    predict.blbglm(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))
  )
})

test_that("Parallelization works and takes less time to run", {
  expect_lt(
    system.time(blbglm(Species ~ Sepal.Length * Sepal.Width, data = iris, family = binomial(), m = 30, B = 10000, nthreads = 4))[1],
    system.time(blbglm(Species ~ Sepal.Length * Sepal.Width, data = iris, family = binomial(), m = 30, B = 10000))[1]
  )
})
