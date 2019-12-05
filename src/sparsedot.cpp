#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sparseDiagCross(S4 A, NumericMatrix M, S4 B, int n) {
  IntegerVector i1 = A.slot("i");
  IntegerVector p1 = A.slot("p");
  NumericVector x1 = A.slot("x");
  IntegerVector i2 = B.slot("i");
  IntegerVector p2 = B.slot("p");
  NumericVector x2 = B.slot("x");

  NumericVector R(n);
  int n1,s1,n2,s2,row1,row2;
  double total;

  for (int j=0; j<n; j++) {
    n1 = p1[j+1] - p1[j]; // number of elements in column j in A
    s1 = p1[j];           // start  of column j in A
    n2 = p2[j+1] - p2[j]; // number of elements in column j in B
    s2 = p2[j];           // start  of column j in B
    total = 0;

    // we skip if there are no elements
    if (n1*n2==0) continue;

    // we compute the cross-product
    for (int k1=0; k1<n1; k1++) {
      for (int k2=0; k2<n2; k2++) {
        row1 = i1[s1 + k1];
        row2 = i2[s2 + k2];
        // get the row of current element
        total = total + M(row1,row2) * x1[s1 + k1] * x2[s2 + k2];
      }
    }
    R(j) = total;
  }
  return R;


//
//   double total = 0;
//   NumericVector::iterator it;
//   for(it = x.begin(); it != x.end(); ++it) {
//     total += *it;
//   }


  // NumericVector::iterator it;
  // IntegerVector dims = mat.slot("Dim");
  // arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
  // arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
  // arma::vec x = Rcpp::as<arma::vec>(mat.slot("x"));
  // int nrow = dims[0], ncol = dims[1];
  // arma::sp_mat res(i, p, x, nrow, ncol);
  // if (show) Rcpp::Rcout << res << std::endl;
  // return res;
}


// [[Rcpp::export]]
NumericVector sparseDiagDot(S4 A, NumericVector D, S4 B, int n) {
  IntegerVector i1 = A.slot("i");
  IntegerVector p1 = A.slot("p");
  NumericVector x1 = A.slot("x");
  IntegerVector i2 = B.slot("i");
  IntegerVector p2 = B.slot("p");
  NumericVector x2 = B.slot("x");

  NumericVector R(n);
  int n1,s1,n2,s2,row1,row2;
  double total;

  for (int j=0; j<n; j++) {
    n1 = p1[j+1] - p1[j]; // number of elements in column j in A
    s1 = p1[j];           // start  of column j in A
    n2 = p2[j+1] - p2[j]; // number of elements in column j in B
    s2 = p2[j];           // start  of column j in B
    total = 0;

    // we skip if there are no elements
    if (n1*n2==0) continue;

    // we compute the cross-product
    for (int k1=0; k1<n1; k1++) {
      for (int k2=0; k2<n2; k2++) {
        row1 = i1[s1 + k1];
        row2 = i2[s2 + k2];
        if (row1 != row2) continue;

        // get the row of current element
        total = total + D(row1) * x1[s1 + k1] * x2[s2 + k2];
      }
    }
    R(j) = total;
  }
  return R;


  //
  //   double total = 0;
  //   NumericVector::iterator it;
  //   for(it = x.begin(); it != x.end(); ++it) {
  //     total += *it;
  //   }


  // NumericVector::iterator it;
  // IntegerVector dims = mat.slot("Dim");
  // arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
  // arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
  // arma::vec x = Rcpp::as<arma::vec>(mat.slot("x"));
  // int nrow = dims[0], ncol = dims[1];
  // arma::sp_mat res(i, p, x, nrow, ncol);
  // if (show) Rcpp::Rcout << res << std::endl;
  // return res;
}


// [[Rcpp::export]]
NumericVector sparseSandwichTrace(S4 L, NumericMatrix M, S4 Q, S4 R, int n, int n0) {
  IntegerVector il = L.slot("i");
  IntegerVector pl = L.slot("p");
  NumericVector xl = L.slot("x");

  IntegerVector ir = R.slot("i");
  IntegerVector pr = R.slot("p");
  NumericVector xr = R.slot("x");

  IntegerVector iq = Q.slot("i");
  IntegerVector pq = Q.slot("p");
  NumericVector xq = Q.slot("x");

  int nr,sr,nl,sl,row,col;
  double total;

  NumericVector D(n);
  NumericVector L0(n0);
  NumericVector R0(n0);

  // we loop over the outer dimension, should be given by n
  for (int j=0; j<n; j++) {

    // we start by pre-computing the within vector L0 and R0
    nl = pl[j+1] - pl[j]; // number of elements in column j in A
    sl = pl[j];           // start  of column j in A
    nr = pr[j+1] - pr[j]; // number of elements in column j in B
    sr = pr[j];           // start  of column j in B
    for (int i=0; i<n0; i++) {
      L0[i] = 0;
      R0[i] = 0;

      for (int k=0; k<nl; k++) {L0(i) += xl(sl + k)*M(i,il(sl + k)) ;} // *M(j,il(sl + k))
      for (int k=0; k<nr; k++) {R0(i) += xr(sr + k)*M(i,ir(sr + k)) ;} // *M(j,il(sr + k))
    }

    // we loop over the elements of Q and construct the sum
    // we can use an iterator, and we need to keep track
    // of row and columns
    total = 0;
    col   = 0;
    row   = 0;
    for (int i=0; i<iq.size(); i++) {
      row = iq(i);
      while (pq(col+1) <= i) col+=1;
      total += L0(row) * R0(col) * xq(i);
      //if (j==0) Rcout << "j=" << j << " row=" << row << " col=" << col << " q=" << xq(i) << "L0=" << L0(row) << "R0= "<< R0(row) << " total= " << total << "\n";
    }

    D(j) = total;
  }
  return D;
}

// [[Rcpp::export]]
NumericVector sparseSandwichDiag(S4 L, NumericMatrix M, S4 Q, S4 R, int n, int n0) {
  IntegerVector il = L.slot("i");
  IntegerVector pl = L.slot("p");
  NumericVector xl = L.slot("x");

  IntegerVector ir = R.slot("i");
  IntegerVector pr = R.slot("p");
  NumericVector xr = R.slot("x");

  IntegerVector iq = Q.slot("i");
  IntegerVector pq = Q.slot("p");
  NumericVector xq = Q.slot("x");

  int nr,sr,nl,sl,row,col;
  double total;

  NumericVector D(n);
  NumericVector L0(n0);
  NumericVector R0(n0);

  // we loop over the outer dimension, should be given by n
  for (int j=0; j<n; j++) {

    // we start by pre-computing the within vector L0 and R0
    nl = pl[j+1] - pl[j]; // number of elements in column j in A
    sl = pl[j];           // start  of column j in A
    nr = pr[j+1] - pr[j]; // number of elements in column j in B
    sr = pr[j];           // start  of column j in B
    for (int i=0; i<n0; i++) {
      L0[i] = 0;
      R0[i] = 0;

      for (int k=0; k<nl; k++) {L0(i) += xl(sl + k)*M(i,il(sl + k)) ;} // *M(j,il(sl + k))
      for (int k=0; k<nr; k++) {R0(i) += xr(sr + k)*M(i,ir(sr + k)) ;} // *M(j,il(sr + k))
    }

    // we loop over the elements of Q and construct the sum
    // we can use an iterator, and we need to keep track
    // of row and columns
    total = 0;
    col   = 0;
    row   = 0;
    for (int i=0; i<iq.size(); i++) {
      row = iq(i);
      while (pq(col+1) <= i) col+=1;
      total += L0(row) * R0(col) * xq(i);
      //if (j==0) Rcout << "j=" << j << " row=" << row << " col=" << col << " q=" << xq(i) << "L0=" << L0(row) << "R0= "<< R0(row) << " total= " << total << "\n";
    }

    D(j) = total;
  }
  return D;
}

// [[Rcpp::export]]
double denseTraceProd(NumericMatrix A, NumericMatrix B) {
  double total =0;
  int ni = A.nrow();
  int nj = A.ncol();
  
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      total += A(i,j) * B(j,i);
    }
  }
  
  return total;
}


// accu.moms.pairs <-function(y1,y2,accu) {
//   if (length(y1)==0) return(accu)
//     if (length(y2)==0) return(accu)
//       
//       if (any(is.na(accu))) {
//         accu=list(n=0,v12=0,s1=0,s2=0)
//       }
//       
//       for ( i in 1:length(y1)) {
//         accu$v12  = accu$v12 + sum(y1[i]*y2)
//         accu$s1   = accu$s1  + y1[i]*length(y2)
//         accu$s2   = accu$s2  + sum(y2)
//         accu$n    = accu$n   + length(y2)
//       }
//       
//       return(accu)
// }

// [[Rcpp::export]]
List paircov(NumericVector A, NumericVector B) {
  double prod =0;
  double m1   =0;
  double m2   =0;
  int    nn   =0;
  
  int ni = A.size();
  int nj = B.size();
  
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      prod += A(i) * B(j);
      m1   += A(i);
      m2   += B(j);
      nn++;
    }
  }
  
  List L = List::create(Named("v12") = prod , Named("m1") = m1, Named("m2") = m2, Named("nn") = nn);
  return L;
}

// [[Rcpp::export]]
List paircov_exclude(NumericVector A, NumericVector B, IntegerVector CA, IntegerVector CB) {
  double prod =0;
  double m1   =0;
  double m2   =0;
  int    nn   =0;
  
  int ni = A.size();
  int nj = B.size();
  
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
      if (CA(i) != CB(j)) {
        prod += A(i) * B(j);
        m1   += A(i);
        m2   += B(j);
        nn++;
      }
    }
  }
  
  List L = List::create(Named("v12") = prod , Named("m1") = m1, Named("m2") = m2, Named("nn") = nn);
  return L;
}


