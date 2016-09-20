FastLoop <- cxxfunction(signature(N = "numeric", M = "numeric", F = "numeric", R = "numeric"), 
      plugin = "RcppArmadillo" , body = '
arma::mat M1 = Rcpp::as < arma::mat >(N);
arma::mat M2 = Rcpp::as < arma::mat >(M);
arma::mat M3 = Rcpp::as < arma::mat >(F);
arma::mat M4 = Rcpp::as < arma::mat >(R);

int i=0,j=0;
for(i = 0; i < M1.n_rows - 1; ++i){

M1(i,0) = M1(i,0) + M4(i,0); // Add the recruitment
//std::cout << "Adding " << i << " is " << M4.at(i,0) << std::endl;
for(j = 0; j < M1.n_cols - 1; j++){
M1(i+1,j+1) = M1(i,j) * exp(- (M2(i,j) + M3(i,j))); 
}
}
return Rcpp::List::create(Rcpp::Named("result") = M1);
')

# Usage
#tmp <- FastLoop(as.matrix(Nb.at.age), as.matrix(Nat.mortality), as.matrix(Fish.mortality), as.matrix(Recruitment))

 



               
