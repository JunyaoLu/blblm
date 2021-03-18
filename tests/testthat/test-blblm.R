test_that("Functions return same results", {
  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
  expect_identical(coef(fit), coef.blblm(fit))
  expect_identical(confint(fit, c("wt", "hp")), confint.blblm(fit, c("wt", "hp")))
  expect_identical(sigma(fit), sigma.blblm(fit))
  expect_identical(
    predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170))),
    predict.blblm(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))
  )
})

test_that("Parallelization works and takes less time to run", {
  rand_data <- data.frame(replicate(10, sample(0:100, 1000, rep = TRUE)))
  expect_lt(
    system.time(blblm(X1 ~ X3 * X5, data = rand_data, m = 10, B = 10000, nthreads = 4))[1],
    system.time(blblm(X1 ~ X3 * X5, data = rand_data, m = 10, B = 10000))[1]
  )
})
