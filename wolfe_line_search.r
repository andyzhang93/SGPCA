
###########################       wolfe   ###################
wolfe <- function(x0, c1, c2){
# Function wolfe performs Line search satisfying Wolfe conditions
# Input
#   myFx:   the optimized function handle
#   d:      the direction we want to search
#   x0:     vector of initial start
#   fx0:    the function value at alpha=0
#   gx0:    the gradient value at alpha=0
# Output
#   alphas: step size
#   fs:     the function value at alpha=alphas
#   gs:     the gradient value at alpha=alphas
#   
    # loss and gradient for initial start
    fx0 <- log_like(x0)
    gx0 <- deriv(x0)

    maxIter <- 3
    alpham <- 20
    alphap <- 0
    alphax <- 1
    # gradient about alpha when alpha=0
    gx0 <- -(1/2)*norm(A, type = c("F"))
    # loss and gradient for previous iteration
    fxp <- fx0
    gxp <- gx0
    i <- 1
    # sparsity constriant
    s<-0.1
    # Line search algorithm satisfying Wolfe conditions
    # alphap is alpha_{i-1}
    # alphax is alpha_i
    # alphas is what we want.
    while (true){
        # update for next iteration
        xx <- solve(diag(n)+(alphax/2)*W)%*%(diag(n)-(alphax/2)*W)%*%x0 
        # loss and gradient for next iteration
        fxx <- log_like(xx)
        gxx <- deriv(xx)
        fs <- fxx
        gs <- gxx
        # gradient about alpha
        gxx<-sum(diag(t(gxx)%*%(-(1/2)*solve(diag(n)+(alphax/2)*W)%*%(x0+xx))))
        
        if ((fxx > fx0 + c1*alphax*gx0) || ((i > 1) && (fxx >= fxp))){
            alphas <- Zoom(problem,x0,d,alphap,alphax,fx0,gx0)
        }
        if (gxx >= c2*gx0){
            alphas <- alphax
        }
        if (gxx >= 0){
            alphas <- Zoom(problem,x0,d,alphax,alphap,fx0,gx0)
        }
        if (sum(sapply(xx,abs)) <=s){
          
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


###########################  Zoom   #######################
Zoom <- function(x0,alphal,alphah,fx0,gx0){
    # This function is called by wolfe
    c1 <- 1e-4
    c2 <- 0.9
    i  <- 0
    maxIter <- 5

    while (true){
        # bisection
        alphax <- 0.5*(alphal+alphah)
        alphas <- alphax
        # update for next iteration
        xx <- solve(diag(n)+(alphax/2)*W)%*%(diag(n)-(alphax/2)*W)%*%x0 
        # loss and gradient for next iteration
        fxx <- log_like(xx)
        gxx <- deriv(xx)
        fs <- fxx
        gs <- gxx
        # gradient about alpha
        gxx<-sum(diag(t(gxx)%*%(-(1/2)*solve(diag(n)+(alphax/2)*W)%*%(x0+xx))))
        
        # update for next iteration
        xl <- solve(diag(n)+(alphal/2)*W)%*%(diag(n)-(alphal/2)*W)%*%x0
        fxl <- log_like(xl)

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


