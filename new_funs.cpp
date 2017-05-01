// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

NumericVector varFun(IntegerVector toB_d, 
                     NumericVector toB_s, 
                     IntegerVector B_d, 
                     NumericVector B_s, 
                     double theta, 
                     double dT, 
                     double sT) {
  
  int n = toB_d.size();
  int m = B_d.size();
  
  
  double duv = 0;
  double suv = 0;
  
  NumericVector var_vec(n);
  
  for(int u = 0; u < n; ++u) {
  
    double total_var_u = 0;
  
      for(int v = 0; v < m; ++v) {
  
        duv = toB_d[u] * B_d[v] / dT;
        if(duv > 1){
          duv = 1;
        }
        suv = toB_s[u] * B_s[v] / sT;
  
        total_var_u += suv * suv * (1 - duv + theta) / duv;
  
      }
  
    var_vec[u] = total_var_u;
  
  }
  
  return var_vec;
}

// [[Rcpp::export]]

double get_set_pvalue0(double stat,
                       NumericVector sB, 
                       IntegerVector dB, 
                       double theta, 
                       double dT, 
                       double sT) {
  // set size
  int m = sB.size();
  
  // starting calculations
  double total_mean = pow(sum(sB), 2) / sT;
  NumericVector var_vec = varFun(dB, sB, dB, sB, theta, dT, sT);
  double total_var = 2 * sum(var_vec);
  
  // getting sum of diagonal variance
  double diag_sum = 0;
  double suu = 0;
  double duu = 0;
  for(int u = 0; u < m; u++) {
    suu = sB[u] * sB[u] / sT;
    duu = dB[u] * dB[u] / dT;
    diag_sum += pow(suu, 2) * (1 - duu + theta) / duu;
  }
  
  total_var -= diag_sum;
  double pval = 1 - R::pnorm(stat, total_mean, sqrt(total_var), 1, 0);
  return pval;
}

// [[Rcpp::export]]

double get_set_pvalue(double stat,
                      NumericVector B_s, 
                      IntegerVector B_d, 
                      double theta, 
                      double dT, 
                      double sT) {
  // set size
  int m = B_s.size();
  
  // getting mean
  double total_mean = (pow(sum(B_s), 2) - sum(B_s * B_s)) / sT;

  //Rcout << total_mean << std::endl;
  
  // getting half-variance sum
  double half_var_sum = 0;
  double duv = 0;
  double suv = 0;
  
  for(int u = 0; u < m; ++u) {
    
    for(int v = 0; v < u; ++v) {
      
      duv = B_d[u] * B_d[v] / dT;
      if(duv > 1){
        duv = 1;
      }
      suv = B_s[u] * B_s[v] / sT;
      
      half_var_sum += suv * suv * (1 - duv + theta) / duv;
      
    }
    
  }

  //Rcout << half_var_sum << std::endl;
  //Rcout << stat << std::endl;
  //Rcout << R::pnorm(stat, total_mean, sqrt(4 * half_var_sum), 1, 0) << std::endl;

  double pval = R::pnorm(stat, total_mean, sqrt(4 * half_var_sum), 0, 0);
  return pval;
}

// [[Rcpp::export]]

double get_set_z(double stat,
                 NumericVector B_s, 
                 IntegerVector B_d, 
                 double theta, 
                 double dT, 
                 double sT) {
  // set size
  int m = B_s.size();
  
  // getting mean
  double total_mean = (pow(sum(B_s), 2) - sum(B_s * B_s)) / sT;
  
  //Rcout << total_mean << std::endl;
  
  // getting half-variance sum
  double half_var_sum = 0;
  double duv = 0;
  double suv = 0;
  
  for(int u = 0; u < m; ++u) {
    
    for(int v = 0; v < u; ++v) {
      
      duv = B_d[u] * B_d[v] / dT;
      if(duv > 1){
        duv = 1;
      }
      suv = B_s[u] * B_s[v] / sT;
      
      half_var_sum += suv * suv * (1 - duv + theta) / duv;
      
    }
    
  }
  
  //Rcout << half_var_sum << std::endl;
  //Rcout << stat << std::endl;
  //Rcout << R::pnorm(stat, total_mean, sqrt(4 * half_var_sum), 1, 0) << std::endl;
  
  double z = (stat - total_mean) / sqrt(4 * half_var_sum);
  return z;
}

// [[Rcpp::export]]

IntegerVector indx_seq(int start, int end) {
  IntegerVector ret_range(end - start);
  ret_range = seq_len(end - start + 1) + start - 1;
  return ret_range;
}

// [[Rcpp::export]]
IntegerVector csample_int(IntegerVector x, 
                          int size,
                          bool replace, 
                          NumericVector prob = NumericVector::create()) {
  IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}

// [[Rcpp::export]]
IntegerVector simple_match_int(IntegerVector x,
                               IntegerVector y){
  return match(x, y);
}

// [[Rcpp::export]]
IntegerVector which_in_int(IntegerVector x,
                           IntegerVector y){

  IntegerVector z = match(x, y);
  LogicalVector res0 = !is_na(z);
  IntegerVector indx = seq_along(x);
  indx = indx[res0];
  return indx;
}

// [[Rcpp::export]]
LogicalVector in_int(IntegerVector x,
                     IntegerVector y){
  
  IntegerVector z = match(x, y);
  LogicalVector res0 = !is_na(z);
  return res0;
}


// [[Rcpp::export]]
NumericVector get_set_stats(IntegerVector b_indx,
                            IntegerVector e_indx,
                            IntegerVector set_indx,
                            IntegerVector i,
                            IntegerVector j,
                            NumericVector w){
  
  int nSets = b_indx.size();
  NumericVector set_stats(nSets);
  double increment = nSets / 10;
  int inc_pos = 0;
  Rcout << "--0%";
  
  for(int k = 0; k < nSets; k++){
    
    if (floor(k / increment) > inc_pos) {
      inc_pos += 1;
      Rcout << "--" << inc_pos * 10 << "%";
    }

    //Rcout << "b_indx[k] - 1 = " << b_indx[k] - 1 << std::endl;
    //Rcout << "e_indx[k] - 1 = " << e_indx[k] - 1 << std::endl;
    //Rcout << "indx_seq(14, 18) = " << indx_seq(14, 18) << std::endl;
    
    IntegerVector indx_seq_k = indx_seq(b_indx[k] - 1, e_indx[k] - 1);

    //Rcout << "indx_seq_k = " << indx_seq_k << std::endl;

    IntegerVector B = set_indx[indx_seq_k];
    IntegerVector set_indx_i = which_in_int(i, B);
    IntegerVector j_indxs = j[set_indx_i - 1];
    IntegerVector set_indx_ij = which_in_int(j_indxs, B);
    set_indx_ij = set_indx_i[set_indx_ij - 1];
    NumericVector weights_k = w[set_indx_ij - 1];
    set_stats[k] = 2 * sum(weights_k);
    
  }
  
  //Rcout << "--100%" << std::endl;
  
  return set_stats;
  
}

