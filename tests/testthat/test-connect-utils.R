test_that("mrf3 normalizes block names and basic accessors work", {
  fit <- structure(
    list(
      data = list(expr = data.frame(a = 1:3), meth = data.frame(b = 4:6)),
      models = list(expr_meth = list(ntree = 10)),
      weights = list(expr = c(a = 0.9)),
      selected_vars = list(expr = "a"),
      vs_summary = data.frame(block = "expr", p_total = 1, p_selected = 1),
      connection = list(c("expr", "meth"))
    ),
    class = "mrf3_fit"
  )

  expect_equal(get_data(fit)$expr$a, 1:3)
  expect_equal(get_models(fit)$expr_meth$ntree, 10)
  expect_equal(get_weights(fit)$expr, c(a = 0.9))
  expect_equal(get_selected_vars(fit)$expr, "a")
  expect_equal(get_vs_summary(fit)$block, "expr")
  expect_equal(get_connection(fit), list(c("expr", "meth")))
})
