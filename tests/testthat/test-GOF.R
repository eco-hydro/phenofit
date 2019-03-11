context("GOF")

test_that("GOF works", {
    Y_obs <- c(1, 3, 4, 2, NA, 6)
    Y_sim <- c(2, 2, 4, 1, NA, 5)
    w     <- c(1, 1, 1, 0.5, 1, 1)

    expect_silent(GOF(Y_obs, Y_sim, w, include.cv = T))

    I <- 3:5
    expect_message(GOF(Y_obs[I], Y_sim[I], w[I], include.cv = T, include.r = T),
        "not enough finite observations")

    I <- 5 # NA value test
    expect_true(all(is.na(GOF(Y_obs[I], Y_sim[I], w[I], include.cv = T))))
})

