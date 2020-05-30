#include "trigo-basis.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <fstream>

using DoubleVector = std::vector<double>;
using Term = std::array<int, 4>;
using Poly = std::vector<Term>;
using Row = std::vector<Poly>;

namespace {

  constexpr size_t DERIVATIVES = 2;
  size_t table_rows;

  std::vector<Row> triangles[DERIVATIVES+1];

  double evalPoly(const Poly &poly, double u) {
    double S = std::sin(M_PI * u / 2), C = std::cos(M_PI * u / 2);
    double a = 1 - S, b = S + C - 1, c = 1 - C;
    double result = 0.0;
    for (const auto &term : poly)
      result += term[3] * std::pow(a, term[0]) * std::pow(b, term[1]) * std::pow(c, term[2]);
    return result;
  }

}

// Assumes that all triangles are empty
void trigoinit(std::string filename) {
  std::ifstream f(filename.c_str());
  f.exceptions(std::ios::failbit | std::ios::badbit);
  size_t polys, terms;
  int a, b, c, k;
  f >> table_rows;
  for (size_t r = 0; r < table_rows; ++r) {
    Row row[DERIVATIVES+1];
    f >> polys;
    for (size_t p = 0; p < polys; ++p) {
      for (size_t d = 0; d <= DERIVATIVES; ++d) {
        Poly poly;
        f >> terms;
        for (size_t t = 0; t < terms; ++t) {
          f >> a >> b >> c >> k;
          poly.push_back({ a, b, c, k });
        }
        row[d].push_back(poly);
      }
    }
    for (size_t d = 0; d <= DERIVATIVES; ++d)
      triangles[d].push_back(row[d]);
  }
}

void trigobasis(size_t n, double u, size_t derivatives, std::vector<DoubleVector> &coeffs) {
  if (derivatives > DERIVATIVES)
    throw std::runtime_error(std::string("The table only has ") + std::to_string(DERIVATIVES) +
                             " derivatives");
  if (n < 2 || n > table_rows + 1)
    throw std::runtime_error(std::string("The table only has rows for 3 to ") +
                             std::to_string(table_rows + 2) + " control points");
  coeffs.resize(derivatives + 1);
  for (size_t d = 0; d <= derivatives; ++d) {
    const auto &row = triangles[d][n-2];
    coeffs[d].clear();
    std::transform(row.begin(), row.end(), std::back_inserter(coeffs[d]),
                   [u](const Poly &p) { return evalPoly(p, u); });
  }
}
