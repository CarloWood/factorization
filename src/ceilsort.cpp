#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>

namespace mp = boost::multiprecision;
using Number = mp::cpp_int;

Number ceil_sqrt(Number const& n)
{
  Number root = mp::sqrt(n);
  return (root * root == n) ? root : root + 1;
}

Number square(Number const& n)
{
  return n * n;
}

int main()
{
  Number N = 12389236789;
  Number M = ceil_sqrt(N);
  Number Q = M;

  Number Qi_plus_1;
  unsigned int iterations = 0;
  for (Number Qi = Q;; Qi = Qi_plus_1)
  {
    //std::cout << "Q = " << Qi << '\n';
    Qi_plus_1 = ceil_sqrt(square(ceil_sqrt(square(Qi) - N)) + N);
    ++iterations;

    if (Qi_plus_1 == Qi)
    {
      Q = Qi_plus_1;
      break;
    }
    if (Qi_plus_1 - Qi > 1)
      std::cout << (Qi_plus_1 - Qi) << '\n';
  }

  Number k = mp::sqrt(square(Q) - N);
  Number x = Q - k;
  Number y = Q + k;
  std::cout << "N = " << N << " = " << x << " * " << y << std::endl;
  std::cout << "iterations = " << iterations << std::endl;
}
