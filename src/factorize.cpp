#include <iostream>
#include <cstdint>
#include <cmath>
#include <cassert>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>

namespace mp = boost::multiprecision;
using Number = mp::cpp_int;

bool is_perfect_square(Number n)
{
  auto root = mp::sqrt(n);
  return root * root == n;
}

bool is_perfect_square(unsigned long n, unsigned long p)
{
  if (n == 0)
    return true;

  // This problem is solved by using a construct `(n|p)` called the Legendre symbol,
  // which is a function that gives an indication of whether a number is a square modulo p.
  // (n|p) = 1 iff n is a quadratic residue modulo p (i.e., if there is an integer k such that k^2 = n [mod p]).
  //
  // The Legendre symbol (n|p) can be computed using Euler's Criterion, which states that
  // (n|p) = n^((p-1)/2) [mod p]. If the result is 1, then a is a quadratic residue modulo p, otherwise it is not.

  // Calculate the exponent.
  unsigned long exp = (p - 1) / 2;      // For example, exp = 1 + 4 + 32, then (n|p) = n^(1 + 4 + 32) = n^1 * n^4 * n^32.
  unsigned np = 1;                      // Legendre symbol (n|p).
  unsigned long nq = n % p;             // n^q [mod p], starting with q = 1.
  for (unsigned long p2 = 1; p2 != 0; p2 <<= 1)
  {
    if ((exp & p2))
      np = (np * nq) % p;
    nq = (nq * nq) % p;                 // Update n^q [mod p], as if q is multiplied with 2.
  }

  return np == 1;
}

mp::cpp_int ceil_sqrt(mp::cpp_int const& n)
{
  mp::cpp_int root = mp::sqrt(n);
  return (root * root == n) ? root : root + 1;
}

std::vector<unsigned long> make_primes(unsigned long upperlimit)
{
  std::vector<unsigned long> primes{2};
  std::vector<bool> sieve(upperlimit / 2, true); // prime = 2 * index + 3.

  for (unsigned long i = 0;; ++i)
    if (sieve[i])
    {
      unsigned long prime = 2 * i + 3;
      primes.push_back(prime);

      unsigned long j = (prime * prime - 3) / 2;
      if (j >= sieve.size())
      {
        for (unsigned long k = i + 1; k < sieve.size(); ++k)
          if (sieve[k])
            primes.push_back(2 * k + 3);
        break;
      }

      while (j < sieve.size())
      {
        sieve[j] = false;
        j += prime;
      }
    }

  return primes;
}

unsigned long sqrt_mod(unsigned long n, unsigned long p)
{
  assert(n < p);
  for (unsigned long r = 0; r < p; ++r)
    if ((r * r) % p == n)
      return r;
  // Don't call sqrt_mod if n isn't a perfect square mod p.
  assert(false);
}

int main()
{
  unsigned long const prime_upper_limit = 1000; //1160479;
  int const mw = 3;
  auto primes = make_primes(prime_upper_limit);

  Number N;
  Number M;

  // Fabricate a large composit number N that exists of two prime factors.
  Number x = 273869673293;
  Number y = 4839572399;
  N = x * y;

  M = ceil_sqrt(N);
  std::cout << "M = " << M << '\n';

  Number L = (x + y) / 2 - M;                 // 102948413048 = 2^3 · 13 · 853 · 1160479
  std::cout << "L = " << L << '\n';

  Number k = (x - y) / 2;
  std::cout << "k = " << k << '\n';

  // N = x * y
  //
  // x = (M + L + k)
  // y = (M + L - k)
  //
  // k = (x - y)/2                      --> 4 k^2 = (x - y)^2 = x^2 - 2xy + y^2
  // L = (x + y)/2 - M                  --> 4 (L + M)^2 = (x + y)^2 = x^2 + 2xy + y^2
  //
  // 4 (L + M)^2 - 4 k^2 = x^2 + 2xy + y^2 - (x^2 - 2xy + y^2) = 4 xy --> (L + M)^2 - k^2 = N --> k^2 = (L + M)^2 - N
  //
  // x + y = 2M + 2L, L >= 0
  //
  // s(x) = x + y = x + N/x >= 2 * sqrt(N) >= 2M.

  // Check that N = x * y, and what the associated values of L and k are for all calculated primes.
  for (int i = 1; i < primes.size(); ++i)       // Start with prime = 3.
  {
    unsigned long p = primes[i];
    std::cout << "p = " << std::setw(mw) << p;

    unsigned long N_mod_p = (N % p).convert_to<unsigned long>();
    unsigned long M_mod_p = (M % p).convert_to<unsigned long>();
    unsigned long x_mod_p = (x % p).convert_to<unsigned long>();
    unsigned long y_mod_p = (y % p).convert_to<unsigned long>();
    unsigned long L_mod_p = (L % p).convert_to<unsigned long>();
    unsigned long k_mod_p = (k % p).convert_to<unsigned long>();

    unsigned long xy_mod_p = (x_mod_p * y_mod_p) % p;
    assert(N_mod_p == xy_mod_p);

    unsigned long L_plus_M_mod_p = (L_mod_p + M_mod_p) % p;
    unsigned long LplusM_squared_mod_p = (L_plus_M_mod_p * L_plus_M_mod_p) % p;
    unsigned long LplusMsquared_minus_N_mod_p = (LplusM_squared_mod_p - N_mod_p + p) % p;       // This is `e` below.
    unsigned long k_squared_mod_p = (k_mod_p * k_mod_p) % p;
    assert(k_squared_mod_p == LplusMsquared_minus_N_mod_p);

    std::cout << "; N = " << std::setw(mw) << N_mod_p;
    std::cout << "; M = " << std::setw(mw) << M_mod_p;
    std::cout << "; x = " << std::setw(mw) << x_mod_p;
    std::cout << "; y = " << std::setw(mw) << y_mod_p;
    std::cout << "; L = " << std::setw(mw) << L_mod_p;
    std::cout << "; k = " << std::setw(mw) << k_mod_p;
    std::cout << "; e = " << std::setw(mw) << LplusMsquared_minus_N_mod_p;
    std::cout << " [mod " << p << "]\n";
  }

  std::cout << "Trying to factorize " << N << "...\n";

  for (int i = 0; i < primes.size(); ++i)
  {
    unsigned long p = primes[i];
    std::cout << "Start p = " << p << '\n';

    std::vector<int> count_k(p);
    std::vector<unsigned long> sol_x(p);
    std::vector<unsigned long> sol_y(p);

    unsigned long Nmp = (N % p).convert_to<unsigned long>();
    unsigned long Mmp = (M % p).convert_to<unsigned long>();

    // For debugging, use the known value of L, x and y.
    unsigned L_mod_p = (L % p).convert_to<unsigned long>();
    unsigned x_mod_p = (x % p).convert_to<unsigned long>();
    unsigned y_mod_p = (y % p).convert_to<unsigned long>();

    //std::cout << "N mod " << p << " = " << Nmp << '\n';
    //std::cout << "M mod " << p << " = " << Mmp << '\n';
    unsigned long emp = (Mmp * Mmp - Nmp + p) % p;                      // e = (L + M)^2 - N, thus with L=0, e = M^2 - N.
    for (unsigned long Lmp = 0; Lmp < p; ++Lmp)
    {
      std::cout << "  Trying L = " << Lmp << " [mod " << p << "]\n";
      std::cout << "    --> e = " << emp << " [mod " << p << "]";
      // The question is: for which L is e = (M + L)^2 - N a perfect square?
      if (is_perfect_square(emp, p))
      {
        std::cout << " which is a perfect square mod p!\n";
        // Note that if (L1 + M)^2 - N is a perfect square, then so is (-L1 - M)^2 - N, because that is the same value.
        // Therefore, as long as L1 + M != 0 [mod p], there will always exist an another solution L2 != L1 such that
        // L2 + M = -L1 - M --> L2 = -L1 - 2M corresponding to x2 = -y1 and y2 = -x1.
        // Prove:
        // L2 negates 'L + M'.
        // If x1 = M + L1 + k and y1 = M + L1 - k then negating M + L1 gives:
        //    x2 = M + L2 + k = -(M + L1) + k = -(M + L1 - k) = -y1.
        //    y2 = M + L2 - k = -(M + L1) - k = -(M + L1 + k) = -x1.
        // Thus x2 * y2  = (-y1) * (-x1) = x1 * y1 = N.
        //
        // Also note that L + M = 0 then there is only a single such value L that gives rise this perfect square.
        // In that case, since L = (x + y)/2 - M, x + y = 2(L + M) = 0 and thus x = -y [mod p].
        // Then from k = (x - y)/2 follows that k = x.
        //
        // For a given N = x * y [mod p], it seems that there should be a solution for
        // every value of x, except 0; giving rise to p - 1 solutions.
        //
        // x * x^-1 = 1, where 1 is its own inverse and one other element not 0 or 1 is its own inverse.
        // The remaining p - 3 elements get paired up.
        //
        // Since N != 0, it has an inverse N^-1 (possibly equal to N).
        //
        // Hence, just pick an x, not 0, and then calculate y = x^-1 N [mod p].

        unsigned long kmp = sqrt_mod(emp, p);
        unsigned long xmp = (Mmp + Lmp + kmp) % p;
        unsigned long ymp = (Mmp + Lmp - kmp + p) % p;
        unsigned long product_mp = (xmp * ymp) % p;
        assert(product_mp == Nmp);

        std::cout << "  k = " << kmp << " [mod " << p << "]; " <<
          product_mp << " = " << xmp << " * " << ymp << " [mod " << p << "]\n";

        if (L_mod_p != Lmp)
        {
          std::cout << "    but this is not a real solution! The real solution is: ";
          std::cout << Nmp << " = " << x_mod_p << " * " << y_mod_p << " [mod " << p << "]\n";
        }
        else
          std::cout << "    Success for p = " << p << '\n';
      }
      else
      {
        std::cout << " which is not a square.\n";
      }

      // e = (L + M)^2 - N (for the current L).
      // The next value of e should be (L+1 + M)^2 - N.
      // This means we need to add the difference to the current e.
      // The difference is: ((L+1 + M)^2 - N) - ((L + M)^2 - N) = (L+1)^2 + 2(L+1)M + M^2 - N - (L^2 + 2LM + M^2 - N) =
      // (L+1)^2 + 2(L+1)M - L^2 - 2LM = 2L + 1 + 2M = 2(L + M) + 1.
      emp = (emp + 2 * (Lmp + Mmp) + 1) % p;
    }
  }
}
