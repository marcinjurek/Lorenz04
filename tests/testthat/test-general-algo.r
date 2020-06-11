N = 960
Force = 10
Tmax = 1
K = 31
M = 5
dt = 0.005
Tmax = 5


params = list(
    list(dt = dt/5, M = 5*M, K = K,  Tmax = Tmax),
    list(dt = dt,   M = M,   K = K,  Tmax = Tmax),
    list(dt = dt*5, M = M/5, K = K,  Tmax = Tmax),
    list(dt = dt/5, M = 5*M, K = K+1, Tmax = Tmax),
    list(dt = dt,   M = M,   K = K+1, Tmax = Tmax),
    list(dt = dt*5, M = M/5, K = K+1, Tmax = Tmax)
)

for( setNo in 1:length(params)){

    set = params[[setNo]]
        
    X0 = getX0(N, Force, set$K, dt)
    
    XV   = VEnKF::Lorenz04M2SimCpp(X0, Force, set$K, set$dt, set$M, set$Tmax, burn = 0, newAlgo = FALSE)
    Xnew =        Lorenz04M2SimCpp(X0, Force, set$K, set$dt, set$M, set$Tmax, 0, FALSE)
    test_name = paste("new algorithm produces the same results as Will's algorithm, set ", setNo, sep="")
    test_that(test_name, {
        expect_equal( sum(abs(Xnew[,set$Tmax] - XV[,set$Tmax])), 0 )
    })
}






                                        
