#ifndef TIME_HPP_
#define TIME_HPP_

#include <chrono>
#include <limits>

class Timer {
 public:
  Timer() : m_timelimit(std::numeric_limits<double>::infinity()),
            m_start(std::chrono::high_resolution_clock::now()) {}

  Timer(double timelimit) : m_timelimit(timelimit),
                            m_start(std::chrono::high_resolution_clock::now()) {}

  double elapsed() const {
    auto now = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration<double>(now - m_start);
    return dur.count();
  }

  bool finished() const {
    return elapsed() > m_timelimit;
  }

 private:
  double m_timelimit;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
};

#endif  // TIME_HPP_
