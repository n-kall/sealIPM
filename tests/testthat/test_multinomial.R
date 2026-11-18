test_that("allocation remains finite when the survivor reference is rare", {
    for (expected in list(
        c(1e-20, 100, 100, 100),
        c(1e-20, 100, 0, 200),
        c(1e-20, 0, 100, 200),
        c(1e-20, 0, 0, 100),
        c(100, 1e-20, 100, 100)
    )) {
        for (noise in c(-2, 0, 2)) {
            counts <- as.numeric(test_model$functions$allocate_fates(
                expected, rep(noise, 3), sum(expected),
                as.integer(expected[2] > 0), as.integer(expected[3] > 0)
            ))
            expect_true(all(is.finite(counts)))
            expect_true(all(counts >= 0))
            expect_equal(sum(counts), sum(expected), tolerance = 1e-12)
        }
    }
})

test_that("structured factor preserves the original allocation mapping", {
    # At ordinary scales, compare against the original dense Cholesky factor.
    for (expected in list(c(10, 20, 30, 40), c(40, 20, 30, 10),
                         c(10, 40, 0, 20), c(10, 0, 40, 20), c(10, 0, 0, 20))) {
        active <- which(expected > 0)
        K <- length(active)
        a <- expected[active] * (1 + 1 / min(expected[active]))
        old_mean <- digamma(a[-1]) - digamma(a[1])
        old_cov <- diag(trigamma(a[-1]), nrow = K - 1) + trigamma(a[1])
        z <- seq_len(K - 1) / 3
        old_logits <- c(0, old_mean + t(chol(old_cov)) %*% z)
        weights <- exp(old_logits - max(old_logits))
        expected_allocation <- numeric(4)
        expected_allocation[active] <- sum(expected) * weights / sum(weights)
        actual <- as.numeric(test_model$functions$allocate_fates(
            expected, z, sum(expected),
            as.integer(expected[2] > 0), as.integer(expected[3] > 0)
        ))
        expect_equal(actual, expected_allocation, tolerance = 1e-12)
    }
})

test_that("multinomial allocation approximates rmultinom", {
    set.seed(123)

    B <- 10000
    N <- 100
    expected_counts <- c(25, 25, 25, 25)
    p <- expected_counts / sum(expected_counts)

    allocation <- replicate(
        B,
        as.numeric(test_model$functions$allocate_fates(
            expected_fate_counts = expected_counts,
            fate_allocation_noise = rnorm(3),
            demographic_population_count = N,
            hunting_active_sweden = 1L,
            hunting_active_finland = 1L
        ))
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
        as.numeric(test_model$functions$allocate_fates(
            expected_fate_counts = expected_counts,
            fate_allocation_noise = rnorm(3),
            demographic_population_count = N,
            hunting_active_sweden = 1L,
            hunting_active_finland = 1L
        ))
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
