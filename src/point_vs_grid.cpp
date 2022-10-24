#include <Rcpp.h>
#include "windowMean.cpp"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector vector_to_bin(NumericVector indat, float threshold)
{

  int ni = indat.length();
  NumericVector result(ni);

  for (int i = 0; i < ni; i++)
  {
    result(i) = (indat(i) >= threshold);
  }
  return result;
}

// [[Rcpp::export]]
NumericVector window_sum_from_cumsum_for_ij(NumericMatrix indat, int rad, NumericMatrix indices)
{
  // windowed average
  // input matrix is output from cumsum2d[_bin]
  // rad is an integer (>=0), window size is 2*rad+1
  // zero-padding
  // TODO : other boundary options:
  //    - periodic or mirror
  //    - reduce rad close to border
  int i, j, ni = indat.nrow(), nj = indat.ncol();
  int no = indices.ncol();

  int imax, jmax;
  NumericVector result(no);
  for (int k = 0; k < no; k++)
  {
    i = (int)indices(0, k);
    j = (int)indices(1, k);

    imax = std::min(i + rad, ni - 1);
    jmax = std::min(j + rad, nj - 1);
    result(i, j) = indat(imax, jmax);
    if (i > rad)
    {
      result(i, j) -= indat(i - rad - 1, jmax);
      if (j > rad)
        result(i, j) += indat(i - rad - 1, j - rad - 1) - indat(imax, j - rad - 1);
    }
    else if (j > rad)
      result(i, j) -= indat(imax, j - rad - 1);
  }

  return result;
}

// [[Rcpp::export]]
DataFrame harpSpatial_point_vs_grid_scores(NumericVector obfield, NumericMatrix indices, NumericMatrix fcfield,
                                           NumericVector thresholds, NumericVector scales, int startegy)
{
  // indices is Nx2 matrix the location of the grid point in the model passed grid.
  // indices(0,:) is the first index and indices(1,:) is the second index of fcfield which represents the model grid.

  int i, j, k, th, sc;
  double a, b, c, dd;
  int n_thresholds = thresholds.length();
  int n_scales = scales.length();
  int ni = fcfield.nrow(), nj = fcfield.ncol();
  int no = obfield.length();

  NumericVector sum_fc(no), bin_ob(no);
  NumericMatrix cum_fc(ni, nj);

  // numeric vectors for the result
  NumericVector res_thresh(n_thresholds * n_scales);
  NumericVector res_size(n_thresholds * n_scales);

  for (th = 0; th < n_thresholds; th++)
  {
    // calculate cumsum matrices for given threshold

    switch (startegy)
    {
    case 0:
      // Pragmatic stratigy
      cum_fc = cumsum2d_bin(fcfield, thresholds[th]);
      break;
    case 1:
      // code block
      break;
    default:
      // code block
    }

    bin_ob = vector_to_bin(obfield, thresholds[th]);
    for (sc = 0; sc < n_scales; sc++)
    {
      k = th * n_scales + sc;
      res_thresh(k) = thresholds(th);
      res_size(k) = scales(sc);
      // fraction matrices

      switch (startegy)
      {
      case 0:
        // Pragmatic stratigy
        sum_fc = window_sum_from_cumsum_for_ij(cum_fc, (int)scales[sc], indices);
        break;
      case 1:
        // code block
        break;
      default:
        // code block
      }

      res_fss[k] = (fss2 < 1.0E-3) ? 0. : 1. - fss1 / fss2;
      res_a[k] = a / (ni * nj);
      res_b[k] = b / (ni * nj);
      res_c[k] = c / (ni * nj);
      res_d[k] = 1. - res_a[k] - res_b[k] - res_c[k];
    } // sc
  }   // th

  return Rcpp::DataFrame::create(Named("threshold") = res_thresh,
                                 Named("scale") = res_size,
                                 Named("fss") = res_fss,
                                 Named("a") = res_a,
                                 Named("b") = res_b,
                                 Named("c") = res_c,
                                 Named("d") = res_d);
}
