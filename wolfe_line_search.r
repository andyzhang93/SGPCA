
############
wolfe <- function(problem, d, x0, c1, c2){
#function [alphas,fs,gs] = wolfe(myFx,d,x0,fx0,gx0)
# Function wolfe performs Line search satisfying Wolfe conditions
# Input
#   myFx:   the optimized function handle
#   d:      the direction we want to search
#   x0:     vector of initial start
#   fx0:    the function value at x0
#   gx0:    the gradient value at x0
# Output
#   alphas: step size
#   fs:     the function value at x0+alphas*d
#   gs:     the gradient value at x0+alphas*d
#   
# Notice
#   I use f and g to save caculation. This funcion wolfe is called by LBFGS_opt.m.

    fx0 <- problem.cost(x0)
    gx0 <- problem.full_grad(x0)

    maxIter <- 3
    alpham <- 20
    alphap <- 0
    alphax <- 1
    gx0 <- t(d)%*%gx0
    fxp <- fx0
    gxp <- gx0
    i <- 1
    # sparsity constriant
    s<-0.01
    # Line search algorithm satisfying Wolfe conditions
    # alphap is alpha_{i-1}
    # alphax is alpha_i
    # alphas is what we want.
    while (true){
        xx <- x0 + alphax*d

        fxx <- problem.cost(xx)
        gxx <- problem.full_grad(xx)
        fs <- fxx
        gs <- gxx
        gxx <- t(d)%*%gxx
        if ((fxx > fx0 + c1*alphax*gx0) || ((i > 1) && (fxx >= fxp))){
            alphas <- Zoom(problem,x0,d,alphap,alphax,fx0,gx0)
        }
        if (gxx >= c2*gx0){
            alphas <- alphax
        }
        if (gxx >= 0){
            alphas <- Zoom(problem,x0,d,alphax,alphap,fx0,gx0)
        }

        alphap <- alphax
        fxp <- fxx
        gxp <- gxx

        if (i > maxIter){
          alphas <- alphax
        }
        #randomly choose alphax from interval (alphax,alpham)
        r <- runif(1)
        alphax <- alphax + (alpham-alphax)*r
        i <- i+1
    }
}


#######
Zoom <- function(problem,x0,d,alphal,alphah,fx0,gx0){
    # This function is called by wolfe
    c1 <- 1e-4
    c2 <- 0.9
    i  <- 0
    maxIter <- 5

    while (true){
        # bisection
        alphax <- 0.5*(alphal+alphah)
        alphas <- alphax
        xx <- x0 + alphax*d

        fxx <- problem.cost(xx)
        gxx <- problem.full_grad(xx)
        fs <- fxx
        gs <- gxx
        gxx <- t(d)%*%gxx
        xl <- x0 + alphal*d

        fxl <- problem.cost(xl)

        if ((fxx > fx0 + c1*alphax*gx0) || (fxx >= fxl)){
            alphah <- alphax
        } else {
            if (gxx >= c2*gx0){
                alphas <- alphax
            }
            if (gxx*(alphah-alphal) >= 0){
                alphah <- alphal
            }
            alphal <- alphax
        }
        i <- i+1
        if (i > maxIter){
            alphas <- alphax
        }
    }
    return(alphas)
}


