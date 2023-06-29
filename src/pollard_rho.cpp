#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>

using namespace boost::multiprecision;

cpp_int powerMod(cpp_int base, cpp_int exp, cpp_int mod)
{
  cpp_int result = 1;
  base = base % mod;
  while (exp > 0)
  {
    if (exp % 2 == 1)
      result = (result * base) % mod;
    exp = exp >> 1;
    base = (base * base) % mod;
  }
  return result;
}

#if 0
cpp_int gcd(cpp_int a, cpp_int b)
{
    if (b == 0)
        return a;
    return gcd(b, a % b);
}
#endif

int main()
{
    cpp_int N("9292514988101982959252099254231932340096077522980419051");
    cpp_int F = 1;                      // b!
    cpp_int twoF = 2;                   // 2^F
    cpp_int GCD = 1;                    // gcd(twoF - 1, N)
    for (int b = 1;; b++)
    {
      if (GCD != 1 && GCD != N)
      {
        std::cout << "b: " << b << ", F mod N: " << F << ", 2^F mod N: " << twoF << ", gcd(2^F - 1, N): " << GCD << std::endl;
        break;
      }
      F *= b + 1;
      F %= N;
      twoF = powerMod(2, F, N);
      GCD = gcd(twoF - 1, N);
    }
}

