#pragma once
#include <chrono>
#include <cmath>
#include <geogram/mesh/mesh.h>
#include <omp.h>

namespace untangler
{

inline std::vector<int> range(int n)
{
    std::vector<int> r(n);
    for (int i = 0; i < n; ++i)
        r[i] = i;
    return r;
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct UntangleResult
{
    bool has_inversions;
    double energy;
    UntangleResult(bool inv = 0, double E = 0) : has_inversions(inv), energy(E)
    {
    }
};

template <typename T> constexpr T square(T x)
{
    return x * x;
}

inline double chi(double eps, double det)
{
    if (det > 0.0)
    {
        return 0.5 * (det + std::sqrt(eps * eps + det * det));
    }
    return 0.5 * eps * eps / (std::sqrt(eps * eps + det * det) - det);
}

inline double chi_deriv(double eps, double det)
{
    return 0.5 + det / (2.0 * std::sqrt(eps * eps + det * det));
}

struct Timer
{
    const char *name;
    double t0;
    Timer(const char *n) : name(n), t0(omp_get_wtime())
    {
    }
    ~Timer()
    {
        std::cout << name << ": " << omp_get_wtime() - t0 << " s\n";
    }
};

} // namespace untangler
