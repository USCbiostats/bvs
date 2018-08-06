#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// z <- eta[good] + (y - mu)[good] / mu.eta.val[good]

// [[Rcpp::export]]
NumericVector compute_response(const LogicalVector & good,
                               const NumericVector & eta,
                               const NumericVector & y,
                               const NumericVector & mu,
                               const NumericVector & mu_eta_val) {
    NumericVector z(sum(good));
    for (int i = 0; i < good.length(); ++i) {
        if (good[i]) {
            z[i] = eta[i] + (y[i] - mu[i]) / mu_eta_val[i];
        }
    }
    return z;
}

// sqrt((weights[good] * mu.eta.val[good]^2) / variance(mu)[good])

// [[Rcpp::export]]
NumericVector compute_weights(const LogicalVector & good,
                              const NumericVector & weights,
                              const NumericVector & mu_eta_val,
                              const NumericVector & mu,
                              const NumericVector & varmu) {
    NumericVector wgt(sum(good));
    for (int i = 0; i < good.length(); ++i) {
        if (good[i]) {
            wgt[i] = sqrt((weights[i] * mu_eta_val[i] * mu_eta_val[i]) / varmu[i]);
        }
    }
    return wgt;
}

// [[Rcpp::export]]
void identity_mu_eta(const NumericVector & eta, NumericVector & eta_val) {
    std::fill(eta_val.begin(), eta_val.end(), 1);
}

// [[Rcpp::export]]
void logit_mu_eta(const NumericVector & eta, NumericVector & eta_val) {

    for (int i = 0; i < eta.length(); ++i) {
        double etai = eta[i];
        if (etai > 30 || etai < -30) {
            eta_val[i] = DOUBLE_EPS;
        } else {
            etai = exp(etai);
            eta_val[i] = etai / ((1 + etai) * (1 + etai));
        }
    }
}

// [[Rcpp::export]]
double logit_deviance(const NumericVector & y,
                      const NumericVector & mu,
                      const NumericVector & wt) {

    double dev = 0.0;
    for (int i = 0; i < y.length(); ++i) {
        dev += wt[i] * ((y[i] > 0) ? log(mu[i]) : log(1 - mu[i]));
    }
    return -2 * dev;
}

// [[Rcpp::export]]
void logit_linkinv(const NumericVector & eta, NumericVector & mu) {
    for (int i = 0; i < eta.length(); ++i) {
        double etai = eta[i], tmp;
        tmp = (etai > 30) ? DOUBLE_EPS : (etai < -30) ? 1 / DOUBLE_EPS : exp(etai);
        mu[i] = tmp / (1 + tmp);
    }
}

// [[Rcpp::export]]
Rcpp::XPtr< std::unordered_map< std::string, int > > new_table(int i) {

    // Creating an unordered_map pointer
    std::unordered_map< std::string, int >* x = new std::unordered_map< std::string, int >;
    // Creating the R external pointer
    Rcpp::XPtr< std::unordered_map< std::string , int > > p(x, true);
    return p;
}

// [[Rcpp::export]]
void set_element_in_table(std::string key, int value, Rcpp::XPtr< std::unordered_map< std::string , int > > table) {
    (*table)[key] = value;
}


// [[Rcpp::export]]
int get_element_from_table(std::string key, Rcpp::XPtr< std::unordered_map< std::string , int > > table) {
    std::unordered_map< std::string, int >::const_iterator it = (*table).find(key);
    if (it != (*table).end()) {
        return it->second;
    }
    return 0;
}

