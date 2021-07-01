
#' @keywords internal
# ' @importFrom spam diag.spam
v_point = function(y, w = 0 * y + 1, lambda = 100, d = 2) {
    # Compute the value of the normalized V-curve for one value of lambda
    # Prepare for smoothing
    n = length(y)
    E = spam::diag.spam(n)
    D = diff(E, diff = d)
    P = t(D) %*% D

    # Smooth for  log-lambdas to the left and to the right
    z     = whit2(y, lambda, w)
    pz    = P %*% z #D' * D * z
    zgrad = lambda * log(10) * whit2(- pz/ w, lambda, w) #whit2(- pz * lambda, lambda, 1)
    # zgrad1 = whit2(-lambda * pz, lambda, w)

    fit   = sum(w * (y - z) ^ 2)
    dlfit = 2 * sum(-zgrad * w * (y - z)) / fit
    pen   = sum(z * pz)
    dlpen = 2 * sum(pz * zgrad) / pen

    # Take distance
    v = sqrt(dlfit ^ 2 + dlpen ^ 2)
    return(v)
}

# sometimes not converge
# v_opt = function(y, w = 0 * y + 1, d = 2, lambdas = c(0, 4), tol = 0.01) {
#     # Locate the optimal value of log10(lambda) with optimizer
#     # Specify bounds of search range for log10(lambda) in paramter 'lambdas'

#     v_fun = function(lla, y, w, d) v_point(y, w, 10 ^ lla, d)
#     op = optimize(v_fun, lambdas, y, w, d, tol = tol)
#     return(op$minimum)
# }
