generateXs = function(X0, params, algo=NULL){
    Xs = list()
    for( setNo in 1:length(params) ){
        set = params[[setNo]]
        cat(paste("generating data for ", paste( names(set), set, sep="=", collapse=", " ), "\n", sep=""))
        if(algo=="VEnKF"){
            X = VEnKF::Lorenz04M2SimCpp(X0, Force, K, set$dt, set$M, iter = set$Tmax, burn = 0, newAlgo = FALSE)
        } else {
            X = Lorenz04M2SimCpp(X0, Force, K, set$dt, set$M, set$Tmax, 0, FALSE)
        }
        Xs[[ setNo ]] = X[,set$Tmax]
    }
    return(Xs)
}


plotComparison = function(Xs, params, save=FALSE) {

    colors = c("red", "blue", "black", "forestgreen")
    ylim = 0
    for( X in Xs ){
        ylim = range(c(ylim, Xs))
    }

    legendStrings = c()

    if( save ){
        pdf(file = "compare.pdf", width=8, height=4)
    }
        
    plot(0, ylim=ylim, xlab="grid points", xlim = c(1, 960), ylab="", type="n", main="Lorenz model under several differencing schemes")
    for( Xno in 1:length(Xs) ){
        set = params[[Xno]]
        X = Xs[[Xno]]
        lines( X, col = colors[[Xno]], type="l")
        legendStrings = c(legendStrings, paste( names(set), set, sep="=", collapse=", " ))
    }
    
    legend("topright",
           legend=legendStrings,
           col=colors,
           lty=rep(1, length(params)),
           lwd = rep(2, length(params)))
    if( save ){
        dev.off()
    }
}    

N = 960
Force = 10
Tmax = 5
K = 32
M = 40
dt = 0.005
X0 = getX0(N, Force, K, dt)


params = list(
    list(dt = 0.005, M = 1, Tmax = Tmax)
)


X = generateXs(X0, params, "new")
plotComparison(X, params, save=FALSE)
