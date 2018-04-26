smooth.spline <- function (x,
              y = NULL,
              w = NULL,
              df,
              spar = NULL,
              lambda = NULL,
              cv = FALSE,
              all.knots = FALSE,
              nknots = .nknots.smspl,
              keep.data = TRUE,
              df.offset = 0,
              penalty = 1,
              control.spar = list(),
              tol = 1e-06 *
                  IQR(x),
              keep.stuff = FALSE)
    {
        contr.sp <- list(
            low = -1.5,
            high = 1.5,
            tol = 1e-04,
            eps = 2e-08,
            maxit = 500,
            trace = getOption("verbose")
        )
        contr.sp[names(control.spar)] <- control.spar
        ctrl.Num <- contr.sp[1:4]
        if (!all(vapply(ctrl.Num, is.numeric, NA)) || contr.sp$tol <
            0 || contr.sp$eps <= 0 || contr.sp$maxit <= 0)
            stop("invalid 'control.spar'")
        xy <- xy.coords(x, y, setLab = FALSE)
        y <- xy$y
        x <- xy$x
        if (!all(is.finite(c(x, y))))
            stop("missing or infinite values in inputs are not allowed")
        n <- length(x)
        if (is.na(n))
            stop("invalid number of points")
        no.wgts <- is.null(w)
        w <- if (no.wgts)
            1
        else {
            if (n != length(w))
                stop("lengths of 'x' and 'w' must match")
            if (any(w < 0))
                stop("all weights should be non-negative")
            if (all(w == 0))
                stop("some weights should be positive")
            (w * sum(w > 0)) / sum(w)
        }
        if (!is.finite(tol) || tol <= 0)
            stop("'tol' must be strictly positive and finite")
        if (!match(keep.stuff, c(FALSE, TRUE)))
            stop("invalid 'keep.stuff'")
        xx <- round((x - mean(x)) / tol)
        nd <- !duplicated(xx)
        ux <- sort(x[nd])
        uxx <- sort(xx[nd])
        nx <- length(ux)
        if (nx <= 3L)
            stop("need at least four unique 'x' values")
        if (nx == n) {
            ox <- TRUE
            tmp <- cbind(w, w * y, w * y ^ 2)[order(x), ]
        }
        else {
            ox <- match(xx, uxx)
            tapply1 <-
                function(X,
                         INDEX,
                         FUN = NULL,
                         ...,
                         simplify = TRUE) {
                    sapply(
                        X = unname(split(X, INDEX)),
                        FUN = FUN,
                        ...,
                        simplify = simplify,
                        USE.NAMES = FALSE
                    )
                }
            tmp <-
                matrix(unlist(
                    tapply1(
                        seq_len(n),
                        ox,
                        if (length(w) ==
                            1L)
                            function(i)
                                c(length(i), sum(y[i]), sum(y[i] ^ 2))
                        else
                            function(i)
                                c(sum(w[i]), sum(w[i] * y[i]), sum(w[i] *
                                                                       y[i] ^
                                                                       2))
                    ),
                    use.names = FALSE
                ),
                ncol = 3,
                byrow = TRUE)
        }
        wbar <- tmp[, 1L]
        ybar <- tmp[, 2L] / ifelse(wbar > 0, wbar, 1)
        yssw <- sum(tmp[, 3L] - wbar * ybar ^ 2)
        if (is.na(cv) && !missing(df))
            stop("'cv' must not be NA when 'df' is specified")
        CV <- !is.na(cv) && cv
        if (CV && nx < n)
            warning("cross-validation with non-unique 'x' values seems doubtful")
        r.ux <- ux[nx] - ux[1L]
        xbar <- (ux - ux[1L]) / r.ux
        if (is.numeric(all.knots)) {
            if (is.unsorted(all.knots, strictly = TRUE))
                stop("Numeric 'all.knots' must be strictly increasing")
            if (!missing(nknots) && !is.null(nknots))
                warning("'all.knots' is vector of knots; 'nknots' specification is disregarded")
            nknots <- length(all.knots)
            if (0 < all.knots[1] || all.knots[nknots] < 1)
                stop("numeric 'all.knots' must cover [0,1] (= the transformed data-range)")
            knot <-
                c(rep(all.knots[1], 3), all.knots, rep(all.knots[nknots],
                                                       3))
        }
        else {
            if (all.knots) {
                if (!missing(nknots) && !is.null(nknots))
                    warning("'all.knots' is TRUE; 'nknots' specification is disregarded")
                nknots <- nx
            }
            else if (is.null(nknots))
                nknots <- .nknots.smspl(nx)
            else {
                if (is.function(nknots))
                    nknots <- nknots(nx)
                else if (!is.numeric(nknots))
                    stop("'nknots' must be numeric (in {1,..,n})")
                if (nknots < 1)
                    stop("'nknots' must be at least 1")
                else if (nknots > nx)
                    stop("cannot use more inner knots than unique 'x' values")
            }
            knot <-
                c(rep(xbar[1], 3), if (all.knots)
                    xbar
                  else
                      xbar[seq.int(1,
                                   nx, length.out = nknots)], rep(xbar[nx], 3))
        }
        nk <- nknots + 2L
        spar.is.lambda <- !missing(lambda)
        if (spar.is.lambda <- !missing(lambda)) {
            if (!missing(spar))
                stop("must not specify both 'spar' and 'lambda'")
            ispar <- 1L
        }
        else
            ispar <- if (is.null(spar) || missing(spar)) {
                if (contr.sp$trace)
                    - 1L
                else
                    0L
            }
        else
            1L
        spar <- if (spar.is.lambda)
            as.double(lambda)
        else if (ispar == 1L)
            as.double(spar)
        else
            double(1)
        if (length(spar) != 1)
            stop("'spar' must be of length 1")
        icrit <- if (is.na(cv))
            0L
        else if (cv)
            2L
        else
            1L
        dofoff <- df.offset
        if (!missing(df)) {
            if (df > 1 && df <= nx) {
                icrit <- 3L
                dofoff <- df
            }
            else
                warning("not using invalid df; must have 1 < df <= n := #{unique x} = ",
                        nx)
        }
        iparms <-
            c(
                icrit = icrit,
                ispar = ispar,
                iter = as.integer(contr.sp$maxit),
                spar.is.lambda
            )
        ans.names <- c("coef",
                       "ty",
                       "lev",
                       "spar",
                       "parms",
                       "crit",
                       "iparms",
                       "ier",
                       if (keep.stuff)
                           "scratch")
        fit <- .Fortran(
            C_rbart,
            as.double(penalty),
            as.double(dofoff),
            x = as.double(xbar),
            y = as.double(ybar),
            w = as.double(wbar),
            ssw = as.double(yssw),
            as.integer(nx),
            as.double(knot),
            as.integer(nk),
            coef = double(nk),
            ty = double(nx),
            lev = double(if (is.na(cv))
                1L
                else
                    nx),
            crit = double(1),
            iparms = iparms,
            spar = spar,
            parms = c(unlist(ctrl.Num),
                      ratio = -1),
            scratch = double((17L + 1L) * nk + 1L),
            ld4 = 4L,
            ldnk = 1L,
            ier = integer(1L)
        )[ans.names]
        if (is.na(cv))
            lev <- df <- NA
        else {
            lev <- fit$lev
            df <- sum(lev)
            if (is.na(df))
                stop("NA lev[]; probably smoothing parameter 'spar' way too large!")
        }
        if (fit$ier > 0L) {
            offKind <- if (spar.is.lambda)
                "extreme"
            else if (sml <- fit$spar < 0.5)
                "small"
            else
                "large"
            wtxt <- paste("smoothing parameter value too", offKind)
            if (spar.is.lambda || sml) {
                stop(wtxt)
            }
            else {
                fit$ty <- rep(mean(y), nx)
                df <- 1
                warning(wtxt, "\nsetting df = 1  __use with care!__")
            }
        }
        cv.crit <- if (is.na(cv))
            NA
        else {
            r <- y - fit$ty[ox]
            if (cv) {
                ww <- wbar
                ww[ww == 0] <- 1
                r <- r / (1 - (lev[ox] * w) / ww[ox])
                if (no.wgts)
                    mean(r ^ 2)
                else
                    weighted.mean(r ^ 2, w)
            }
            else
                (if (no.wgts)
                    mean(r ^ 2)
                 else
                     weighted.mean(r ^ 2, w))
            / (1 - (df.offset + penalty *
                        df) / n) ^ 2
        }
        structure(
            list(
                x = ux,
                y = fit$ty,
                w = wbar,
                yin = ybar,
                tol = tol,
                data = if (keep.data)
                    list(x = x, y = y, w = w),
                no.weights = no.wgts,
                lev = lev,
                cv.crit = cv.crit,
                pen.crit = sum(wbar *
                                   (ybar - fit$ty) ^
                                   2),
                crit = fit$crit,
                df = df,
                spar = if (spar.is.lambda)
                    NA
                else
                    fit$spar,
                ratio = if (spar.is.lambda)
                    NA
                else
                    fit$parms[["ratio"]],
                lambda = fit$parms[["low"]],
                iparms = c(fit$iparms, errorI = if (fit$ier)
                    fit$ier
                    else
                        NA),
                auxM = if (keep.stuff)
                    list(
                        XWy = fit$scratch[seq_len(nk)],
                        XWX = fit$scratch[nk + seq_len(4 * nk)],
                        Sigma = fit$scratch[5 *
                                                nk + seq_len(4 * nk)],
                        R = fit$scratch[9 * nk +
                                            seq_len(4 * nk)]
                    ),
                fit = structure(
                    list(
                        knot = knot,
                        nk = nk,
                        min = ux[1L],
                        range = r.ux,
                        coef = fit$coef
                    ),
                    class = "smooth.spline.fit"
                ),
                call = match.call()
            ),
            class = "smooth.spline"
        )
    }
