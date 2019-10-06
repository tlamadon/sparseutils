context("Sparse methods")
library(sparseutils)
library(testthat)
library(Matrix)

test_that("testing the sparse functions", {

  # create a dense matrix
  M = array(runif(5*5),c(5,5))
  A = Matrix::t(1.0*as(array( (runif(10*5)>0.7)*rnorm(10*5) ,c(10,5)), "sparseMatrix"))
  B = Matrix::t(1.0*as(array( (runif(10*5)>0.7)*rnorm(10*5) ,c(10,5)), "sparseMatrix"))
  r0 = sparseDiagCross(A,M,B,dim(A)[2])
  r1 = Matrix::diag(Matrix::t(A) %*% M %*% B)
  expect_true( mean( (r0-r1)^2 ) < 1e-6)

  # testing a second type of function
  A = Matrix::t(1.0*as(array( (runif(10*5)>0.7)*rnorm(10*5) ,c(10,5)), "sparseMatrix"))
  B = Matrix::t(1.0*as(array( (runif(10*5)>0.7)*rnorm(10*5) ,c(10,5)), "sparseMatrix"))
  Q = Matrix::t(1.0*as(array( (runif(5*5)>0.7)*rnorm(5*5) ,c(5,5)), "sparseMatrix"))
  M = array(runif(5*5),c(5,5))
  M = .5* (M + t(M))
  
  r0 = Matrix::diag(Matrix::t(A) %*% M %*% Q %*% M %*% B)
  r1 = sparseSandwichTrace(A, M, Q, B, 10, 5)
  expect_true( mean( (r0-r1)^2 ) < 1e-6)
  
  
  A = array(runif(10*5),c(10,5))
  B = array(runif(5*10),c(5,10))

  r0 = sum(diag( A %*% B))  
  r1 = denseTraceProd(A,B)  
  expect_true( mean( (r0-r1)^2 ) < 1e-6)
  
})