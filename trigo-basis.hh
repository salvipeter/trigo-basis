#pragma once

#include <string>
#include <vector>

void trigoinit(std::string filename);

void trigobasis(size_t n, double u, size_t derivatives,
                std::vector<std::vector<double>> &coeffs);
