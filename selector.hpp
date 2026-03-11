#pragma once
#include <algorithm>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cut.hpp"

enum class EvalMode {
  EarlyStopping,  // return false on first failed cut (default, fast)
  EvaluateAll     // always run every cut, return false if any failed
};

template <typename P>
class Selector {
 public:
  explicit Selector(std::string name, EvalMode mode = EvalMode::EarlyStopping) : name_(std::move(name)), mode_(mode) {}

  // Add a cut (chainable).
  Selector& add(Cut<P> cut) {
    cuts_.push_back(std::move(cut));
    return *this;
  }

  Selector& set_mode(EvalMode m) {
    mode_ = m;
    return *this;
  }

  // Core evaluation — call this once per particle.
  bool operator()(P p) {
    ++n_calls_;

    if (mode_ == EvalMode::EarlyStopping) {
      for (auto& cut : cuts_) {
        if (!cut(p)) return false;  // cut already increments its own counters
      }
      return true;
    } else {
      // EvaluateAll: run every cut even if one fails
      bool all_pass = true;
      for (auto& cut : cuts_) {
        if (!cut(p)) all_pass = false;
      }
      return all_pass;
    }
  }

  // --- Statistics ---
  long long calls() const { return n_calls_; }

  void print_cutflow(std::ostream& os = std::cout) const {
    os << "\n=== Cutflow: " << name_ << " (invocations: " << n_calls_ << ") ===\n";
    os << std::left << std::setw(30) << "Cut" << std::setw(12) << "Calls" << std::setw(12) << "Passed" << std::setw(12)
       << "Failed"
       << "Eff [%]\n";
    os << std::string(70, '-') << '\n';
    for (const auto& c : cuts_) {
      double eff = c.calls() > 0 ? 100.0 * c.passed() / c.calls() : 0.0;
      os << std::setw(30) << c.name() << std::setw(12) << c.calls() << std::setw(12) << c.passed() << std::setw(12)
         << c.failed() << std::fixed << std::setprecision(1) << eff << '\n';
    }
  }

  void reset() {
    n_calls_ = 0;
    for (auto& c : cuts_) c.reset();
  }

  const std::string& name() const { return name_; }

 private:
  std::string name_;
  EvalMode mode_;
  std::vector<Cut<P>> cuts_;
  long long n_calls_{0};
};