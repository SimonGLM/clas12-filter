#pragma once

#include <iostream>
#include <string>
#include <vector>

enum class EvaluationMode { EarlyReturn, CompleteTrace };

struct CutResult {
  std::string name;
  bool passed;
};

class ParticleContext {
 private:
  EvaluationMode mode;
  std::vector<CutResult> cut_history;
  bool overall_passed = true;

 public:
  ParticleContext(EvaluationMode mode = EvaluationMode::EarlyReturn) : mode(mode) {}

  template <typename Func, typename... Args>
  bool apply_cut(const std::string& cut_name, Func&& cut_func, Args&&... args) {
    bool result = cut_func(std::forward<Args>(args)...);
    cut_history.push_back({cut_name, result});

    if (!result) {
      overall_passed = false;
      return mode == EvaluationMode::CompleteTrace;  // continue if tracing, stop if early-return
    }
    return true;  // continue regardless
  }

  bool passed() const { return overall_passed; }

  void print_trace() const {
    std::cout << "Cut history:\n";
    for (const auto& result : cut_history) {
      std::cout << "  " << (result.passed ? "✓" : "✗") << " " << result.name << "\n";
    }
  }

  const std::vector<CutResult>& get_cut_history() const { return cut_history; }
};
