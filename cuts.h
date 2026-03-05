#pragma once

#include <atomic>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "region_particle.h"

class StatisticsDecorator {
 private:
  struct ChildStatistics {
    int invocations = 0;
    int rejections = 0;
    int getAcceptance() const { return invocations - rejections; }
    double getRejectionRate() const { return rejections * 1. / invocations; }
  };

  std::string name;
  std::atomic<int> invocations{0};
  std::atomic<int> rejections{0};
  std::unordered_map<StatisticsDecorator*, ChildStatistics> child_stats;

  static std::map<std::string, StatisticsDecorator*>& registry() {
    static std::map<std::string, StatisticsDecorator*> reg;
    return reg;
  }

 protected:
  StatisticsDecorator(const std::string& cut_name) : name(cut_name) { registry()[name] = this; }

  void track_call(bool result) {
    invocations++;
    if (!result) {
      rejections++;
    }
  }

  static std::vector<StatisticsDecorator*>& call_stack() {
    static std::vector<StatisticsDecorator*> stack;
    return stack;
  }

 public:
  virtual ~StatisticsDecorator() { registry().erase(name); }

  void on_child_called(StatisticsDecorator* child, bool result) {
    child_stats[child].invocations++;
    if (!result) {
      child_stats[child].rejections++;
    }
  }

  int getInvocations() const { return invocations; }
  int getRejections() const { return rejections; }
  int getAcceptance() const { return invocations - rejections; }
  double getRejectionRate() const { return rejections * 1. / invocations; }
  const std::string& get_name() const { return name; }
  const std::unordered_map<StatisticsDecorator*, ChildStatistics>& getChildStats() const { return child_stats; }

  // Static method to print all statistics
  static void printAllStatistics() {
    std::cout << "\n" << std::string(120, '=') << std::endl;
    std::cout << "CUT STATISTICS SUMMARY" << std::endl;
    std::cout << std::format("{:<47} | {:>15} | {:>15} | {:>15} | {:>15} ", "Cut Name", "Invocations", "Accepted",
                             "Rejections", "Rejection Rate")
              << std::endl;
    std::cout << std::string(120, '=') << std::endl;
    // std::sort(registry().begin(), registry().end(),
    //           [](std::pair<std::string, StatisticsDecorator*> a, std::pair<std::string, StatisticsDecorator*> b) {
    //             return a.first < b.first;
    //           });
    bool odd_row = false;
    for (const auto& [name, decorator] : registry()) {
      if (odd_row) std::cout << "\033[0;30m";
      std::cout << std::format("{:<47} | {:>15} | {:>15} | {:>15} | {:>13.2f}%", name, decorator->getInvocations(),
                               decorator->getAcceptance(), decorator->getRejections(),
                               decorator->getRejectionRate() * 100)
                << std::endl;
      if (odd_row) std::cout << "\033[0m";
      odd_row = !odd_row;
    }
    std::cout << "\033[0m" << std::string(120, '=') << std::endl;
  }

  // Static method to print hierarchical statistics
  static void printHierarchicalStatistics() {
    std::cout << "\n" << std::string(120, '=') << std::endl;
    std::cout << "SELECTOR STATISTICS (with child cut usage)" << std::endl;
    std::cout << std::string(120, '=') << std::endl;

    // Print selectors and their observed children
    for (const auto& [name, decorator] : registry()) {
      if (name.find("selector::") == 0) {
        std::cout << "\n" << name << ":" << std::endl;
        std::cout << std::format("  Invocations: {:>10}  |  Accepted: {:>10}  |  Rejected: {:>10}  |  Rejection Rate: {:>6.2f}%",
                                 decorator->getInvocations(),
                                 decorator->getAcceptance(),
                                 decorator->getRejections(),
                                 decorator->getRejectionRate() * 100)
                  << std::endl;
        
        if (!decorator->child_stats.empty()) {
          std::cout << "  Used cuts (per-selector statistics):" << std::endl;
          for (const auto& [child, stats] : decorator->child_stats) {
            std::cout << std::format("    {:<45} | Invocations: {:>10}  |  Accepted: {:>10}  |  Rejected: {:>10}  |  Rej Rate: {:>6.2f}%",
                                     child->get_name(),
                                     stats.invocations,
                                     stats.getAcceptance(),
                                     stats.rejections,
                                     stats.getRejectionRate() * 100)
                      << std::endl;
          }
        }
      }
    }

    // Print flat list of all cuts
    std::cout << "\n" << std::string(120, '=') << std::endl;
    std::cout << "ALL CUTS STATISTICS SUMMARY" << std::endl;
    std::cout << std::format("{:<47} | {:>15} | {:>15} | {:>15} | {:>15} ", "Cut Name", "Invocations", "Accepted",
                             "Rejections", "Rejection Rate")
              << std::endl;
    std::cout << std::string(120, '=') << std::endl;
    bool odd_row = false;
    for (const auto& [name, decorator] : registry()) {
      if (odd_row) std::cout << "\033[0;30m";
      std::cout << std::format("{:<47} | {:>15} | {:>15} | {:>15} | {:>13.2f}%", name, decorator->getInvocations(),
                               decorator->getAcceptance(), decorator->getRejections(),
                               decorator->getRejectionRate() * 100)
                << std::endl;
      if (odd_row) std::cout << "\033[0m";
      odd_row = !odd_row;
    }
    std::cout << "\033[0m" << std::string(120, '=') << std::endl;
  }
};

template <typename Func>
class DecoratedCut : public StatisticsDecorator {
 private:
  Func func;

 public:
  DecoratedCut(const std::string& cut_name, Func f) : StatisticsDecorator(cut_name), func(f) {}

  template <typename... Args>
  bool operator()(Args&&... args) {
    // Push ourselves onto the call stack
    StatisticsDecorator::call_stack().push_back(this);
    
    // Execute the wrapped function
    bool result = func(std::forward<Args>(args)...);
    
    // Pop ourselves from the call stack
    StatisticsDecorator::call_stack().pop_back();
    
    // Track the result in our own statistics
    track_call(result);
    
    // Notify parent (if any) with our result
    if (!StatisticsDecorator::call_stack().empty()) {
      StatisticsDecorator::call_stack().back()->on_child_called(this, result);
    }
    
    return result;
  }
};

#define DECLARE_CUT(name, func_ptr) DecoratedCut<decltype(func_ptr)> name(#name, func_ptr)

namespace cuts {
  namespace generic {
    bool _forward_detector_cut(clas12::region_particle*);
    bool _central_detector_cut(clas12::region_particle*);
    bool _forward_tagger_cut(clas12::region_particle*);
    bool _PID_cut(clas12::region_particle*, int pid);
    bool _charge_cut(clas12::region_particle*, int charge);

    static DecoratedCut<bool (*)(clas12::region_particle*)> forward_detector_cut("forward_detector_cut",
                                                                                 cuts::generic::_forward_detector_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*)> central_detector_cut("central_detector_cut",
                                                                                 cuts::generic::_central_detector_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*)> forward_tagger_cut("forward_tagger_cut",
                                                                               cuts::generic::_forward_tagger_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*, int)> PID_cut("PID_cut", cuts::generic::_PID_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*, int)> charge_cut("charge_cut", cuts::generic::_charge_cut);
  }  // namespace generic

  namespace FD {
    bool _HTCC_nphe_cut(clas12::region_particle*);
    bool _EC_sampling_fraction_cut(clas12::region_particle*, bool inbending, bool simulation, bool spring2019);
    bool _EC_hit_position_fiducial_cut_homogeneous(clas12::region_particle*, int tightness, bool inbending);
    bool _EC_outer_vs_EC_inner_cut(clas12::region_particle*, int tightness);
    bool _DC_fiducial_cut_edge(clas12::region_particle*, int region, bool inbending);
    bool _DC_z_vertex_cut(clas12::region_particle*);
    bool _phot_EC_sampling_fraction_cut(clas12::region_particle*);
    bool _phot_EC_outer_vs_EC_inner_cut(clas12::region_particle*);

    static DecoratedCut<bool (*)(clas12::region_particle*)> HTCC_nphe_cut("HTCC_nphe_cut", cuts::FD::_HTCC_nphe_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*, bool, bool, bool)> EC_sampling_fraction_cut(
        "EC_sampling_fraction_cut", cuts::FD::_EC_sampling_fraction_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*, int, bool)> EC_hit_position_fiducial_cut_homogeneous(
        "EC_hit_position_fiducial_cut_homogeneous", cuts::FD::_EC_hit_position_fiducial_cut_homogeneous);
    static DecoratedCut<bool (*)(clas12::region_particle*, int)> EC_outer_vs_EC_inner_cut(
        "EC_outer_vs_EC_inner_cut", cuts::FD::_EC_outer_vs_EC_inner_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*, int, bool)> DC_fiducial_cut_edge(
        "DC_fiducial_cut_edge", cuts::FD::_DC_fiducial_cut_edge);
    static DecoratedCut<bool (*)(clas12::region_particle*)> DC_z_vertex_cut("DC_z_vertex_cut",
                                                                            cuts::FD::_DC_z_vertex_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*)> phot_EC_sampling_fraction_cut(
        "phot_EC_sampling_fraction_cut", cuts::FD::_phot_EC_sampling_fraction_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*)> phot_EC_outer_vs_EC_inner_cut(
        "phot_EC_outer_vs_EC_inner_cut", cuts::FD::_phot_EC_outer_vs_EC_inner_cut);
  }  // namespace FD
  namespace FT {
    bool _FT_eid_FTCAL_fiducial_cut(clas12::region_particle*);
    bool _FT_eid_FTTRK_fiducial_cut(clas12::region_particle*);
    bool _FT_eid_FTHODO_fiducial_cut(clas12::region_particle*);
    bool _FT_photid_FTCAL_fiducial_cut(clas12::region_particle*);
    bool _FT_eid_energy_vs_radius_cut(clas12::region_particle*);
    bool _FT_photid_beta_cut(clas12::region_particle*);

    static DecoratedCut<bool (*)(clas12::region_particle*)> FT_eid_FTCAL_fiducial_cut(
        "FT_eid_FTCAL_fiducial_cut", cuts::FT::_FT_eid_FTCAL_fiducial_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*)> FT_eid_FTTRK_fiducial_cut(
        "FT_eid_FTTRK_fiducial_cut", cuts::FT::_FT_eid_FTTRK_fiducial_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*)> FT_eid_FTHODO_fiducial_cut(
        "FT_eid_FTHODO_fiducial_cut", cuts::FT::_FT_eid_FTHODO_fiducial_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*)> FT_photid_FTCAL_fiducial_cut(
        "FT_photid_FTCAL_fiducial_cut", cuts::FT::_FT_photid_FTCAL_fiducial_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*)> FT_eid_energy_vs_radius_cut(
        "FT_eid_energy_vs_radius_cut", cuts::FT::_FT_eid_energy_vs_radius_cut);
    static DecoratedCut<bool (*)(clas12::region_particle*)> FT_photid_beta_cut("FT_photid_beta_cut",
                                                                               cuts::FT::_FT_photid_beta_cut);
  }  // namespace FT
  namespace CD {}
  namespace vertex {
    bool _delta_vz_cut(clas12::region_particle*, double reference_vertex_z);
    static DecoratedCut<bool (*)(clas12::region_particle*, double)> delta_vz_cut("delta_vz_cut",
                                                                                 cuts::vertex::_delta_vz_cut);
  }  // namespace vertex
  bool _phot_beta_cut(clas12::region_particle*, int tightness);
  bool _neutr_beta_cut(clas12::region_particle*, int run);
  bool _CD_neutr_beta_cut(clas12::region_particle*, int run);
  bool _basic_FTOF_cut(clas12::region_particle*);

  static DecoratedCut<bool (*)(clas12::region_particle*, int)> phot_beta_cut("phot_beta_cut", cuts::_phot_beta_cut);
  static DecoratedCut<bool (*)(clas12::region_particle*, int)> neutr_beta_cut("neutr_beta_cut", cuts::_neutr_beta_cut);
  static DecoratedCut<bool (*)(clas12::region_particle*, int)> CD_neutr_beta_cut("CD_neutr_beta_cut",
                                                                                 cuts::_CD_neutr_beta_cut);
  static DecoratedCut<bool (*)(clas12::region_particle*)> basic_FTOF_cut("basic_FTOF_cut", cuts::_basic_FTOF_cut);
}  // namespace cuts

// forward_detector_cut
// central_detector_cut
// forward_tagger_cut
// PID_cut
// charge_cut
// HTCC_nphe_cut
// EC_sampling_fraction_cut
// EC_hit_position_fiducial_cut_homogeneous
// EC_outer_vs_EC_inner_cut
// DC_fiducial_cut_edge
// DC_z_vertex_cut
// phot_EC_sampling_fraction_cut
// phot_EC_outer_vs_EC_inner_cut
// FT_eid_FTCAL_fiducial_cut
// FT_eid_FTTRK_fiducial_cut
// FT_eid_FTHODO_fiducial_cut
// FT_photid_FTCAL_fiducial_cut
// FT_eid_energy_vs_radius_cut
// FT_photid_beta_cut
// delta_vz_cut
// phot_beta_cut
// neutr_beta_cut
// CD_neutr_beta_cut
// basic_FTOF_cut