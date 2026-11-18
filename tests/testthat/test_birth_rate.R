test_that("density intercept is accurate at both extremes", {
    intercept <- test_model$functions$compute_density_dependence_intercept

    tiny <- 1e-20
    expect_equal(
        intercept(tiny, tiny),
        -log(tiny + tiny * (1 - tiny)),
        tolerance = 1e-14
    )

    near_one <- 1 - 1e-10
    complement <- (1 - near_one)^2
    # Check relative accuracy: an absolute tolerance would also accept zero.
    expect_equal(
        intercept(near_one, near_one) / (-log1p(-complement)),
        1,
        tolerance = 1e-14
    )
    expect_gt(intercept(near_one, near_one), 0)
    expect_equal(intercept(0.5, 0.5), -log(0.75))
})

test_that("density intercept agrees across the branch boundary", {
    intercept <- test_model$functions$compute_density_dependence_intercept
    # With one input zero, complement = 1 - the other input.
    for (b in c(0.5 - 1e-10, 0.5, 0.5 + 1e-10)) {
        expect_equal(intercept(b, 0), -log(b), tolerance = 1e-14)
        expect_equal(intercept(0, b), -log(b), tolerance = 1e-14)
    }
})

test_that("density intercept retains its exact boundary behavior", {
    intercept <- test_model$functions$compute_density_dependence_intercept
    expect_equal(intercept(0, 0), Inf)
    expect_equal(intercept(1, 0.5), 0)
    expect_equal(intercept(0.5, 1), 0)
    expect_equal(intercept(0, 0.5), -log(0.5))
    expect_equal(intercept(0.5, 0), -log(0.5))
})
