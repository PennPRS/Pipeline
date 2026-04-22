#ifndef SFBM_H
#define SFBM_H

/******************************************************************************/

#ifndef STRICT_R_HEADERS
#define STRICT_R_HEADERS
#endif

#include <mmap.hpp>
#include <system_error> // for std::error_code
#include <Rcpp.h>

/******************************************************************************/

class SFBM {
public:
  SFBM(std::string path, int n, int m, std::vector<size_t> p) : n(n), m(m), p(p) {

    std::error_code error;
    this->ro_mmap.map(path, error);
    if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());

    this->data = reinterpret_cast<const double*>(ro_mmap.data());
  }

  size_t nrow() const { return n; }
  size_t ncol() const { return m; }

  template<class C>
  C prod(const C& x) {

    C res(n);
    for (int i = 0; i < n; i++) res[i] = 0;

    for (int j = 0; j < m; j++) {

      size_t lo = 2 * p[j];
      size_t up = 2 * p[j + 1];

      for (size_t k = lo; k < up; k += 2) {
        int    ind = data[k];
        double val = data[k + 1];
        res[ind] += val * x[j];
      }
    }

    return res;
  }

  template<class C>
  C cprod(const C& x) {

    C res(m);

    for (int j = 0; j < m; j++)
      res[j] = dot_col(j, x);

    return res;
  }

  template<class C>
  double dot_col(int j, const C& x) {

    size_t lo = 2 * p[j];
    size_t up = 2 * p[j + 1];

    double cp = 0;
    for (size_t k = lo; k < up; k += 2) {
      int    ind = data[k];
      double val = data[k + 1];
      cp += val * x[ind];
    }

    return cp;
  }

private:
  mio::mmap_source ro_mmap;
  const double * data;
  int n;
  int m;
  std::vector<size_t> p;
};

/******************************************************************************/

#endif // SFBM_H
