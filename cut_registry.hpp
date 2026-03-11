#pragma once
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "clas12defs.h"
#include "cut.hpp"

struct CutStatistics {
  std::string name;
  uint64_t invocations = 0;
  uint64_t passes = 0;
  uint64_t failures = 0;

  double pass_rate() const { return invocations > 0 ? (100.0 * passes / invocations) : 0.0; }
};

class CutRegistry {
 private:
  std::vector<std::unique_ptr<Cut>> cuts_;
  std::map<std::string, CutStatistics> stats_;

 public:
  void add_cut(std::unique_ptr<Cut> cut) {
    cuts_.push_back(std::move(cut));
    stats_[cuts_.back()->name()] = {};
  }

  const std::vector<std::unique_ptr<Cut>>& get_cuts() const { return cuts_; }

  CutStatistics& get_stats(const std::string& cut_name) { return stats_[cut_name]; }

  const std::map<std::string, CutStatistics>& get_all_stats() const { return stats_; }

  void print_summary() const {
    for (const auto& [name, stat] : stats_) {
      printf("%s: %lu/%lu (%.1f%%)\n", name.c_str(), stat.passes, stat.invocations, stat.pass_rate());
    }
  }
};