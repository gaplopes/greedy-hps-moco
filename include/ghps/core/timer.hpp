#ifndef TIME_HPP_
#define TIME_HPP_

#include <chrono>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace moco::ghps {


class TimeLimitReached : public std::runtime_error {
 public:
  TimeLimitReached() : std::runtime_error("Computation time limit reached") {}
};

class Timer {
 public:
  Timer() : Timer(std::numeric_limits<double>::infinity()) {}

  Timer(double timelimit) : m_start(Clock::now()), m_deadline(Clock::time_point::max()) {
    if (timelimit < std::chrono::duration<double>(Clock::time_point::max() - m_start).count())
      m_deadline = m_start + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(std::max(0.0, timelimit)));
  }

  Timer(double timelimit, const Timer& deadline) : Timer(timelimit) {
    m_deadline = std::min(m_deadline, deadline.m_deadline);
  }

  double elapsed() const {
    auto now = Clock::now();
    auto dur = std::chrono::duration<double>(now - m_start);
    return dur.count();
  }

  bool finished() const {
    return Clock::now() >= m_deadline;
  }

  double left() const {
    if (m_deadline == Clock::time_point::max()) return std::numeric_limits<double>::infinity();
    return std::max(0.0, std::chrono::duration<double>(m_deadline - Clock::now()).count());
  }

  void check() const { if (finished()) throw TimeLimitReached(); }

 private:
  using Clock = std::chrono::steady_clock;
  Clock::time_point m_start;
  Clock::time_point m_deadline;
};

} // namespace moco::ghps

#endif  // TIME_HPP_
