#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

class StatisticsCollector {
 private:
  struct CutStats {
    int invocations = 0;
    int rejections = 0;

    int getInvocations() const { return invocations; }
    int getRejections() const { return rejections; }
    int getAcceptance() const { return invocations - rejections; }
    double getRejectionRate() const { return rejections * 1. / invocations; }
  };

  struct SelectorStats {
    int invocations = 0;
    int rejections = 0;
    std::unordered_map<std::string, CutStats> cuts;
  };

  static std::map<std::string, CutStats>& registry() {
    static std::map<std::string, CutStats> reg;
    return reg;
  }

  static std::map<std::string, SelectorStats>& selector_registry() {
    static std::map<std::string, SelectorStats> reg;
    return reg;
  }

 public:
  static void record_cut(const std::string& cut_name, bool result) {
    auto& stats = registry()[cut_name];
    stats.invocations++;
    if (!result) {
      stats.rejections++;
    }
  }

  static void record_cut_in_selector(const std::string& selector_name, const std::string& cut_name, bool result) {
    auto& selector_stats = selector_registry()[selector_name];
    auto& cut_stats = selector_stats.cuts[cut_name];
    cut_stats.invocations++;
    if (!result) {
      cut_stats.rejections++;
    }
  }

  static void record_selector_invocation(const std::string& selector_name, bool passed) {
    auto& selector_stats = selector_registry()[selector_name];
    selector_stats.invocations++;
    if (!passed) {
      selector_stats.rejections++;
    }
  }

  static void print_hierarchical_statistics() {
    std::cout << "\n" << std::string(120, '=') << std::endl;
    std::cout << std::format("{:<47} | {:>15} | {:>15} | {:>15} | {:>15} ",
                             "SELECTOR STATISTICS (with child cut usage)", "Invokations", "Accepted", "Rejected",
                             "Rejection rate")
              << std::endl;
    std::cout << std::string(120, '=') << std::endl;

    // Convert selector registry to a vector and sort by invocations
    std::vector<std::pair<std::string, SelectorStats>> sorted_selectors(selector_registry().begin(),
                                                                        selector_registry().end());
    std::unordered_map<std::string, int> explicit_order = {{"electron", 0}, {"proton", 1},  {"neutron", 2},
                                                           {"piplus", 3},   {"piminus", 4}, {"Kplus", 5},
                                                           {"Kminus", 6},   {"photon", 7}};
    std::sort(sorted_selectors.begin(), sorted_selectors.end(), [&](const auto& a, const auto& b) {
      int idx1, idx2;
      try {
        idx1 = explicit_order.at(a.first);
        idx2 = explicit_order.at(b.first);
      } catch (std::out_of_range& e) {
        // If either selector is not in the explicit order map, fall back to sorting by invocations
        return a.second.invocations > b.second.invocations;
      }
      return idx1 < idx2;
    });

    for (const auto& [selector_name, selector_stats] : sorted_selectors) {
      if (!selector_stats.cuts.empty()) {
        // The selector's invocation count is tracked separately
        int selector_invocations = selector_stats.invocations;
        int selector_rejections = selector_stats.rejections;
        int selector_acceptances = selector_invocations - selector_rejections;

        std::cout << std::format("\033[1m{:<47} | {:>15} | {:>15} | {:>15} | {:>13.2f}%\033[0m", selector_name,
                                 selector_invocations, selector_acceptances, selector_rejections,
                                 selector_invocations > 0 ? selector_rejections * 1. / selector_invocations * 100 : 0)
                  << std::endl;

        // Sort child cuts by invocations
        std::vector<std::pair<std::string, CutStats>> sorted_cuts(selector_stats.cuts.begin(),
                                                                  selector_stats.cuts.end());
        std::sort(sorted_cuts.begin(), sorted_cuts.end(),
                  [](const auto& a, const auto& b) { return a.second.invocations > b.second.invocations; });

        for (const auto& [cut_name, stats] : sorted_cuts) {
          std::cout << std::format("    {:<43} | {:>15} | {:>15} | {:>15} | {:>13.2f}%", cut_name,
                                   stats.getInvocations(), stats.getAcceptance(), stats.getRejections(),
                                   stats.getInvocations() > 0 ? stats.getRejectionRate() * 100 : 0)
                    << std::endl;
        }
      }
    }
    std::cout << "\033[0m" << std::string(120, '=') << std::endl;
  }
};
