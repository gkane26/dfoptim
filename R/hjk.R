##
##  h o o k e j e e v e s . R  Hooke-Jeeves Minimization Algorithm
##




#' Hooke-Jeeves derivative-free minimization algorithm
#' 
#' An implementation of the Hooke-Jeeves algorithm for derivative-free
#' optimization.  A bounded and an unbounded version are provided.
#' 
#' Argument \code{control} is a list specifing changes to default values of
#' algorithm control parameters.  Note that parameter names may be abbreviated
#' as long as they are unique.
#' 
#' The list items are as follows: \describe{ \item{list("tol")}{Convergence
#' tolerance. Iteration is terminated when the step length of the main loop
#' becomes smaller than \code{tol}. This does \emph{not} imply that the optimum
#' is found with the same accuracy.  Default is 1.e-06.}
#' 
#' \item{list("maxfeval")}{Maximum number of objective function evaluations
#' allowed. Default is Inf, that is no restriction at all.}
#' 
#' \item{list("maximize")}{A logical indicating whether the objective function
#' is to be maximized (TRUE) or minimized (FALSE). Default is FALSE.}
#' 
#' \item{list("target")}{A real number restricting the absolute function value.
#' The procedure stops if this value is exceeded.  Default is Inf, that is no
#' restriction.}
#' 
#' \item{list("info")}{A logical variable indicating whether the step number,
#' number of function calls, best function value, and the first component of
#' the solution vector will be printed to the console. Default is FALSE.} } If
#' the minimization process threatens to go into an infinite loop, set either
#' \code{maxfeval} or \code{target}.
#' 
#' @aliases hjk hjkb
#' @param par Starting vector of parameter values.  The initial vector may lie
#' on the boundary. If \code{lower[i]=upper[i]} for some \code{i}, the
#' \code{i}-th component of the solution vector will simply be kept fixed.
#' @param fn Nonlinear objective function that is to be optimized.  A scalar
#' function that takes a real vector as argument and returns a scalar that is
#' the value of the function at that point.
#' @param lower,upper Lower and upper bounds on the parameters.  A vector of
#' the same length as the parameters.  If a single value is specified, it is
#' assumed that the same bound applies to all parameters. The starting
#' parameter values must lie within the bounds.
#' @param control A list of control parameters.  See \bold{Details} for more
#' information.
#' @param \dots Additional arguments passed to \code{fn}.
#' @return A list with the following components: \item{par}{Best estimate of
#' the parameter vector found by the algorithm.}
#' 
#' \item{value}{value of the objective function at termination.}
#' 
#' \item{convergence}{indicates convergence (\code{=0}) or not (\code{=1}).}
#' 
#' \item{feval}{number of times the objective \code{fn} was evaluated.}
#' 
#' \item{niter}{number of iterations in the main loop.}
#' @note This algorithm is based on the Matlab code of Prof. C. T. Kelley,
#' given in his book ``Iterative methods for optimization".  It is implemented
#' here with the permission of Prof. Kelley.
#' 
#' This version does not (yet) implement a cache for storing function values
#' that have already been computed as searching the cache makes it slower.
#' @author Hans W Borchers <hwborchers@@googlemail.com>
#' @seealso \code{\link{optim}}, \code{\link{nmk}}
#' @references C.T. Kelley (1999), Iterative Methods for Optimization, SIAM.
#' 
#' Quarteroni, Sacco, and Saleri (2007), Numerical Mathematics, Springer.
#' @keywords optimize
#' @examples
#' 
#' ##  Hooke-Jeeves solves high-dim. Rosenbrock function
#'   rosenbrock <- function(x){
#'     n <- length(x)
#'     sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
#'   }
#' par0 <- rep(0, 10)
#' hjk(par0, rosenbrock)
#' 
#' hjkb(c(0, 0, 0), rosenbrock, upper = 0.5)
#' # $par
#' # [1] 0.50000000 0.25742722 0.06626892
#' 
#' 
#' ##  Hooke-Jeeves does not work well on non-smooth functions
#'   nsf <- function(x) {
#' 	f1 <- x[1]^2 + x[2]^2
#' 	f2 <- x[1]^2 + x[2]^2 + 10 * (-4*x[1] - x[2] + 4)
#' 	f3 <- x[1]^2 + x[2]^2 + 10 * (-x[1] - 2*x[2] + 6)
#' 	max(f1, f2, f3)
#'   }
#' par0 <- c(1, 1)                                 # true min 7.2 at (1.2, 2.4)
#' hjk(par0, nsf) # fmin=8 at xmin=(2,2)
#' 
hjk <- function(par, fn, control = list(), ...) {
    if (!is.numeric(par))
        stop("Argument 'par' must be a numeric vector.", call. = FALSE)
    n <- length(par)
    if (n == 1)
        stop("For univariate functions use some different method.", call. = FALSE)
   
    #-- Control list handling ----------
    cntrl <- list(tol      = 1.e-06,
                  maxfeval = Inf,       # set to Inf if no limit wanted
                  maximize = FALSE,     # set to TRUE  for maximization
                  target   = Inf,       # set to Inf for no restriction
                  info     = FALSE)     # for printing interim information
    nmsCo <- match.arg(names(control), choices = names(cntrl), several.ok = TRUE)
    if (!is.null(names(control))) cntrl[nmsCo] <- control

    tol      <- cntrl$tol;      
    maxfeval <- cntrl$maxfeval
    maximize <- cntrl$maximize
    target   <- cntrl$target
    info     <- cntrl$info

	scale <- if (maximize) -1 else 1
    fun <- match.fun(fn)
    f <- function(x) scale * fun(x, ...)

    #-- Setting steps and stepsize -----
    nsteps <- floor(log2(1/tol))        # number of steps
    steps  <- 2^c(-(0:(nsteps-1)))      # decreasing step size
    dir <- diag(1, n, n)                # orthogonal directions

    x <- par                            # start point
    fx <- fbest <- f(x)                 # smallest value so far
    fcount <- 1                         # counts number of function calls

    if (info) cat("step\tnofc\tfmin\txpar\n")   # info header

    #-- Start the main loop ------------
    ns <- 0
    while (ns < nsteps && fcount < maxfeval && abs(fx) < target) {
        ns <- ns + 1
        hjs    <- .hjsearch(x, f, steps[ns], dir, fcount, maxfeval, target)
        x      <- hjs$x
        fx     <- hjs$fx
        sf     <- hjs$sf
        fcount <- fcount + hjs$finc

        if (info)
            cat(ns, "\t",  fcount, "\t", fx/scale, "\t", x[1], "...\n")
    }

    if (fcount > maxfeval) {
        warning("Function evaluation limit exceeded -- may not converge.")
        conv <- 1
    } else if (abs(fx) > target) {
        warning("Function exceeds min/max value -- may not converge.")
        conv <- 1
    } else {
        conv <- 0
    }

    fx <- fx / scale                    # undo scaling
    return(list(par = x, value = fx,
                convergence = conv, feval = fcount, niter = ns))
}

##  Search with a single scale -----------------------------
.hjsearch <- function(xb, f, h, dir, fcount, maxfeval, target) {
    x  <- xb
    xc <- x
    sf <- 0
    finc <- 0
    hje  <- .hjexplore(xb, xc, f, h, dir)
    x    <- hje$x
    fx   <- hje$fx
    sf   <- hje$sf
    finc <- finc + hje$numf

    # Pattern move
    while (sf == 1) {
        d  <- x-xb
        xb <- x
        xc <- x+d
        fb <- fx
        hje  <- .hjexplore(xb, xc, f, h, dir, fb)
        x    <- hje$x
        fx   <- hje$fx
        sf   <- hje$sf
        finc <- finc + hje$numf

        if (sf == 0) {  # pattern move failed
           hje  <- .hjexplore(xb, xb, f, h, dir, fb)
           x    <- hje$x
           fx   <- hje$fx
           sf   <- hje$sf
           finc <- finc + hje$numf
        }
        if (fcount + finc > maxfeval || abs(fx) > target) break
    }

    return(list(x = x, fx = fx, sf = sf, finc = finc))
}

##  Exploratory move ---------------------------------------
.hjexplore <- function(xb, xc, f, h, dir, fbold) {
    n <- length(xb)
    x <- xb

    if (missing(fbold)) {
        fb <- f(x)
        numf <- 1
    } else {
        fb <- fbold
        numf <- 0
    }

    fx <- fb
    xt <- xc
    sf <- 0                             # do we find a better point ?
    dirh <- h * dir
    fbold <- fx
    for (k in sample.int(n, n)) {       # resample orthogonal directions
        p1 <- xt + dirh[, k]
        ft1 <- f(p1)
        numf <- numf + 1

        p2 <- xt - dirh[, k]
        ft2 <- f(p2)
        numf <- numf + 1

        if (min(ft1, ft2) < fb) {
            sf <- 1
            if (ft1 < ft2) {
                xt <- p1
                fb <- ft1
            } else {
                xt <- p2
                fb <- ft2
            }
        }
    }

    if (sf == 1) {
        x <- xt
        fx <- fb
    }
    return(list(x = x, fx = fx, sf = sf, numf = numf))
}
