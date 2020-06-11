N = 960
Force = 15
dt = 0.005
M = 40
K = 32
Tmax = 50

set.seed(1988)
Z0 = runif(N)
print(Z0[1:3])

t0 = proc.time()
Z = Lorenz04M2SimCpp(Z0, Force, K, dt, M, iter = Tmax, burn = 0, vectorAlgo=FALSE)
t1 = proc.time() - t0
cat(paste("Simulation took", t1[3], "\n"))
plot(Z[,Tmax], type="l")

