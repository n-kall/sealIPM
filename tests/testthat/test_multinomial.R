test_that("multinomial allocation approximates rmultinom", {
    set.seed(123)

    B <- 10000
    N <- 100
    expected_counts <- c(25, 25, 25, 25)
    p <- expected_counts / sum(expected_counts)

    allocation <- replicate(
        B,
        test_model$functions$multinomial_allocation(
            expected_counts = expected_counts,
            z = rnorm(length(expected_counts) - 1),
            N = N
        )
    )

    multinomial <- replicate(
        B,
        as.numeric(
            rmultinom(
                n = 1,
                size = N,
                prob = p
            )
        )
    )

    # The continuous allocation must preserve the total population.
    expect_equal(
        colSums(allocation),
        rep(N, B),
        tolerance = 1e-8
    )

    # Both methods should have approximately the same means.
    allocation_mean <- rowMeans(allocation)
    multinomial_mean <- rowMeans(multinomial)

    expect_equal(
        allocation_mean,
        N * p,
        tolerance = 0.5
    )

    expect_equal(
        multinomial_mean,
        N * p,
        tolerance = 0.5
    )

    expect_equal(
        allocation_mean,
        multinomial_mean,
        tolerance = 0.5
    )

    # use a relative tolerance rather than expect_identical().
    allocation_var <- apply(allocation, 1, var)
    multinomial_var <- apply(multinomial, 1, var)

    multinomial_var_theory <- N * p * (1 - p)

    expect_equal(
        multinomial_var,
        multinomial_var_theory,
        tolerance = 0.75
    )

    expect_equal(
        allocation_var,
        multinomial_var_theory,
        tolerance = 3
    )

    # Check the negative dependence between categories.
    allocation_cov <- cov(t(allocation))
    multinomial_cov <- cov(t(multinomial))

    off_diagonal <- row(allocation_cov) != col(allocation_cov)

    expect_true(all(allocation_cov[off_diagonal] < 0))
    expect_true(all(multinomial_cov[off_diagonal] < 0))

    # The covariance matrices should be broadly similar.
    expect_equal(
        allocation_cov,
        multinomial_cov,
        tolerance = 3
    )
})


test_that("allocation approximation works with unequal probabilities", {
    set.seed(456)

    B <- 10000
    N <- 100
    expected_counts <- c(10, 20, 30, 40)
    p <- expected_counts / sum(expected_counts)

    allocation <- replicate(
        B,
        test_model$functions$multinomial_allocation(
            expected_counts = expected_counts,
            z = rnorm(length(expected_counts) - 1),
            N = N
        )
    )

    multinomial <- replicate(
        B,
        as.numeric(rmultinom(1, size = N, prob = p))
    )

    # Total is preserved exactly by the continuous approximation.
    expect_equal(
        colSums(allocation),
        rep(N, B),
        tolerance = 1e-8
    )

    # The means should be close to N * p.
    expect_equal(
        rowMeans(allocation),
        N * p,
        tolerance = 1
    )

    expect_equal(
        rowMeans(multinomial),
        N * p,
        tolerance = 1
    )

    # Compare the variance structure approximately.
    allocation_var <- apply(allocation, 1, var)
    multinomial_var <- apply(multinomial, 1, var)

    expect_equal(
        allocation_var,
        multinomial_var,
        tolerance = 5
    )
})
