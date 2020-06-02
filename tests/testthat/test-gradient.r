roll = function(x, n){
  if ( n == 0 )
    return( x )
  c( tail(x,n), head(x,-n) )
}



N = 20
Force = 10
dt = 0.005
M = 40
Tmax = 1
K = 6

X0 = matrix(rep(1,N))

evolFun = function(X) Lorenz04M2SimCpp(X, Force, K, dt, M, iter = 1, burn = 0, newAlgo=FALSE)

 
test_that("Correct gradient for K=1", {
  G = exactGrad(X0, 1)
  expect_gt(N, 5)
  expect_equal(G[1,1], -1)
  expect_equal(G[1,2], 1)

  expect_equal(G[1,N - 1], -1)
  expect_equal(G[1,N], 0)
  for(i in 2:N){
    expect_true( all(G[i,] == roll(G[1,], i-1)) )
  }
})

test_that("Correct gradient for K=2", {
    G = (2**2) * exactGrad(X0, 2)
    expect_gt(10, 5)
    expect_equal(G[1,1], -0.75)
    expect_equal(G[1,2], 1)
    expect_equal(G[1,3], 2)
    expect_equal(G[1,4], 1)

    expect_equal(G[1,N - 4], -1)
    expect_equal(G[1,N - 3], -1.75)
    expect_equal(G[1,N - 2], -1)
    expect_equal(G[1,N - 1], -0.5)
    expect_equal(G[1,N], 0)
    for (i in 2:N) {
      expect_true( all(G[i,] == roll(G[1,], i-1)) )
    }
})



test_that("Correct gradient for K=3", {
    G = (3**2) * exactGrad(X0, 3)
    expect_gt(N, 15)
    expect_equal(G[1,1], -1)
    expect_equal(G[1,2],  0)
    expect_equal(G[1,3],  3)
    expect_equal(G[1,4],  3)
    expect_equal(G[1,5],  3)

    expect_equal(G[1,N - 6], -3)
    expect_equal(G[1,N - 5], -3)
    expect_equal(G[1,N - 4], -2)
    expect_equal(G[1,N - 3], -1)
    expect_equal(G[1,N - 2],  0)
    expect_equal(G[1,N - 1], -1)
    expect_equal(G[1,N],      1)
    for(i in 2:N){
      expect_true( all(G[i,] == roll(G[1,], i-1)) )
    }
})
