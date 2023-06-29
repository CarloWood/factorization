#include "sys.h"
#include "matplotlibcpp.h"
#include <vector>
#include <iostream>
#include <tuple>
#include <cmath>
#include "debug.h"

using Number = unsigned long;

// Returns true if n is a perfect square; the root is put in root.
bool is_perfect_square(Number n, Number& root)
{
  root = std::sqrt(n);
  return root * root == n;
}

std::tuple<Number, Number, int> factorize_and_count(Number N)
{
  Number M = std::ceil(std::sqrt(N));
  Number Q = M;
  int iterations = 1;
  Number upper_limit_Q = (N + 1) / 2;
  while (Q < upper_limit_Q)
  {
    Number D = Q*Q - N;
    Number k;
    if (is_perfect_square(D, k))
    {
      Number x = Q + k;
      Number y = Q - k;
      return {x, y, iterations};
    }
    ++Q;
    ++iterations;
  }
  return {N, 1, iterations};
}

// If N is prime then we don't find a factorization and will leave the loop
// after number of iterations: I = upper_limit_Q - M = (N + 1) / 2 - ceil(sqrt(N)) --> I ~ N/2
//
// If N = 3 * p, where p is prime, then k = (N/3 - 3)/2 --> Q = sqrt(((N/3 - 3)/2)^2 + N) ~ N/6
// In general, N = f * p, where p is prime and f << p, then Q ~ N/(2 * f).
//
// However, if f approaches sqrt(N) ... say, f = M - R where R is (much?) smaller than M, then
// k = (N / f - f)/2 = (N / (M - R) - (M - R)) / 2 --> 2k(M - R) = N - (M - R)^2 = N - M^2 + 2MR - R^2 -->
// 2k(1 - R/M) = N/M - M + 2R - R^2/M ~ R (2 - R/M) --> k ~ R --> iterations is Q - M = sqrt(N + k^2) - M ~
// sqrt(M^2 + R^2) - M = sqrt(M^2 + 2 M (R^2 / 2M) + (R^2 / 2M)^2) - M = sqrt((M + (R^2 / 2M))^2) - M = R^2 / 2M.
// This is much much smaller than N/(2f).

namespace plt = matplotlibcpp;

class DataPoint
{
 private:
  double x_;
  double y_;
  bool undefined_;

 public:
  DataPoint() : undefined_(true) { }
  DataPoint(double x, double y) : x_(x), y_(y), undefined_(false) { }

  double x() const { ASSERT(!undefined_); return x_; }
  double y() const { ASSERT(!undefined_); return y_; }
};

class Plot
{
 private:
  DataPoint xy_min_;            // The most bottom-left data point in the graph.
  DataPoint xy_max_;            // The most upper-right data point in the graph.

  std::vector<DataPoint> data_points_;

 public:
  void add(DataPoint data_point);
  void draw();
};

void Plot::add(DataPoint data_point)
{
  data_points_.push_back(data_point);
}

void Plot::draw()
{
  // Sort the data by x coordinate.
  std::sort(data_points_.begin(), data_points_.end(),
      [](DataPoint const& lhs, DataPoint const& rhs) {
        return lhs.x() < rhs.x();
      });

  // Split the data into two vectors.
  std::vector<double> plot_x;
  std::vector<double> plot_y;
  for (DataPoint const& data_point : data_points_)
  {
    plot_x.push_back(data_point.x());
    plot_y.push_back(data_point.y());
  }

  // Show the plot.
  plt::scatter(plot_x, plot_y);
  plt::show();
}

/*
Number const Nmax = 1000000;
Number const ymax = 0.0;
Number const xmax = std::log10(Nmax);

void draw_horizontal_line_at(double y)
{
  std::vector<double> x_line = {ymax, xmax};
  std::vector<double> y_line(x_line.size(), y);

  plt::plot(x_line, y_line, {{"color", "lightgray"}, {"linewidth", "0.5"}});
}
*/

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  Plot plot;

  plot.add({3, 9});
  plot.add({2, 4});
  plot.add({4, 16});

  plot.draw();
}
