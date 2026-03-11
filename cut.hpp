#pragma once
#include <atomic>
#include <functional>
#include <string>

// A named, self-counting boolean predicate.
// Wraps any callable: free function, lambda, functor, std::bind expression.
// P is the "primary" particle pointer type; additional args are bound at
// construction time via a capturing lambda so the stored signature is always
// bool(P).
template <typename P>
class Cut {
 public:
  using Fn = std::function<bool(P)>;

  Cut(std::string name, Fn fn) : name_(std::move(name)), fn_(std::move(fn)) {}

  // Evaluate and increment counters.
  bool operator()(P p) {
    ++n_calls_;
    bool result = fn_(p);
    if (result) ++n_pass_;
    return result;
  }

  // --- Accessors ---
  const std::string& name() const { return name_; }
  long long calls() const { return n_calls_; }
  long long passed() const { return n_pass_; }
  long long failed() const { return n_calls_ - n_pass_; }

  void reset() {
    n_calls_ = 0;
    n_pass_ = 0;
  }

 private:
  std::string name_;
  Fn fn_;
  long long n_calls_{0};
  long long n_pass_{0};
};