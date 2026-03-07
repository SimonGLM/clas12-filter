#pragma once

#include <fmt/format.h>
#include <iostream>
#include <map>
#include <string>

class StatisticsCollector {
 private:
  struct CutStats {
    int invocations = 0;
    int rejections = 0;

    int getAcceptance() const { return invocations - rejections; }
    double getRejectionRate() const { return rejections * 1. / invocations; }
  };

  static std::map<std::string, CutStats>& registry() {
    static std::map<std::string, CutStats> reg;
    return reg;
  }

  static std::map<std::string, std::map<std::string, CutStats>>& selector_registry() {
    static std::map<std::string, std::map<std::string, CutStats>> reg;
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
    auto& stats = selector_registry()[selector_name][cut_name];
    stats.invocations++;
    if (!result) {
      stats.rejections++;
    }
  }

  static void print_summary() {
    std::cout << "\n" << std::string(120, '=') << std::endl;
    std::cout << "CUT STATISTICS SUMMARY" << std::endl;
    std::cout << std::string(120, '=') << std::endl;

    bool odd_row = false;
    for (const auto& [name, stats] : registry()) {
      if (odd_row) std::cout << "\033[0;30m";
      std::cout << fmt::format("{:<47} | {:>15} | {:>15} | {:>15} | {:>13.2f}%", name, stats.invocations,
                              stats.getAcceptance(), stats.rejections, stats.getRejectionRate() * 100)
                << std::endl;
      if (odd_row) std::cout << "\033[0m";
      odd_row = !odd_row;
    }
    std::cout << "\033[0m" << std::string(120, '=') << std::endl;
  }

  static void print_hierarchical_statistics() {
    std::cout << "\n" << std::string(120, '=') << std::endl;
    std::cout << "SELECTOR STATISTICS (with child cut usage)" << std::endl;
    std::cout << std::string(120, '=') << std::endl;

    for (const auto& [selector_name, cuts] : selector_registry()) {
      if (!cuts.empty()) {
        int total_invocations = 0;
        int total_acceptances = 0;
        int total_rejections = 0;

        for (const auto& [cut_name, stats] : cuts) {
          total_invocations += stats.invocations;
          total_rejections += stats.rejections;
        }
        total_acceptances = total_invocations - total_rejections;

        std::cout << fmt::format("\033[1m{:<47} | {:>15} | {:>15} | {:>15} | {:>13.2f}%\033[0m", selector_name,
                                total_invocations, total_acceptances, total_rejections,
                                total_rejections * 1. / total_invocations * 100)
                  << std::endl;

        for (const auto& [cut_name, stats] : cuts) {
          std::cout << fmt::format("    {:<43} | {:>15} | {:>15} | {:>15} | {:>13.2f}%", cut_name,
                                  stats.invocations, stats.getAcceptance(), stats.rejections,
                                  stats.getRejectionRate() * 100)
                    << std::endl;
        }
      }
    }
    std::cout << "\033[0m" << std::string(120, '=') << std::endl;
  }
};
