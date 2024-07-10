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
  bool undefined_y_;

 public:
  DataPoint() : undefined_y_(true) { }
  DataPoint(double x, double y) : x_(x), y_(y), undefined_y_(false) { }

  double x() const { return x_; }
  double y() const { ASSERT(!undefined_y_); return y_; }

  double& x() { return x_; }
  double& y() { return y_; }

  void update_min(double y);
  void update_max(double y);
};

void DataPoint::update_min(double y)
{
  if (undefined_y_ || y < y_)
  {
    y_ = y;
    undefined_y_ = false;
  }
}

void DataPoint::update_max(double y)
{
  if (undefined_y_ || y > y_)
  {
    y_ = y;
    undefined_y_ = false;
  }
}

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

void draw_horizontal_line_at(double y, double xmin, double xmax)
{
  std::vector<double> x_line = {xmin, xmax};
  std::vector<double> y_line(x_line.size(), y);

  plt::plot(x_line, y_line, {{"color", "lightgray"}, {"linewidth", "0.5"}});
}

void draw_vertical_line_at(double x, double ymin, double ymax)
{
  std::vector<double> y_line = {ymin, ymax};
  std::vector<double> x_line(y_line.size(), x);

  plt::plot(x_line, y_line, {{"color", "lightgray"}, {"linewidth", "0.5"}});
}

void Plot::draw()
{
  if (!data_points_.empty())
  {
#if 0
    // Sort the data by x coordinate.
    std::sort(data_points_.begin(), data_points_.end(),
        [](DataPoint const& lhs, DataPoint const& rhs) {
          return lhs.x() < rhs.x();
        });
#endif

    xy_min_.x() = data_points_[0].x();
    xy_max_.x() = data_points_[data_points_.size() - 1].x();
  }

  // Split the data into two vectors.
  std::vector<double> plot_x;
  std::vector<double> plot_y;
  for (DataPoint const& data_point : data_points_)
  {
    plot_x.push_back(data_point.x());
    plot_y.push_back(data_point.y());

    xy_min_.update_min(data_point.y());
    xy_max_.update_max(data_point.y());
  }

  double xmin = xy_min_.x();
  double xmax = xy_max_.x();
  double ymin = xy_min_.y();
  double ymax = xy_max_.y();


  for (double x = std::floor(xmin); x <= std::ceil(xmax); ++x)
    draw_vertical_line_at(x, ymin, ymax);
  for (double y = std::floor(ymin); y <= std::ceil(ymax); ++y)
    draw_horizontal_line_at(y, xmin, xmax);

  // Show the plot.
  plt::plot(plot_x, plot_y, {{"color", "green"}, {"linewidth", "0.5"}});
  plt::show();
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  Plot plot;

  double N = 5777;

  for (double x = 48; x <= 58; x += 1)
  {
    double y = N / x;
    double yr = 109.0 - (x - 53.0) * 2.0 + (x - 53.0) * (x - 53.0) * 0.055;

    for (double dx = -0.5; dx <= 0.5; dx += 0.125)
    {
      y = N / (x + dx);
      plot.add({x + dx, y - yr});
    }
  }

  plot.draw();
}
