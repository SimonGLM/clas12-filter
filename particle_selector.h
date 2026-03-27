#pragma once
#include "cuts.h"
#include "particle_context.h"
#include "statistics_collector.h"

// Helper macro to check cut and continue/stop based on mode
#define CHECK_CUT(ctx, name, func, ...) \
  if (!ctx.apply_cut(name, func, __VA_ARGS__)) return ctx.passed()

namespace selectors {
  EvaluationMode evaluation_mode = EvaluationMode::EarlyReturn;
  std::map<int, bool> detector_flags = {{clas12::FT, true}, {clas12::FD, true}, {clas12::CD, true}};

  namespace impl {
    bool electron(ParticleContext& ctx, clas12::region_particle* p, bool inbending, cuts::tightness tightness) {
      // CHECK_CUT(ctx, "PID_cut(11)", cuts::generic::impl::PID_cut, p, 11);
      // CHECK_CUT(ctx, "charge_cut(-1)", cuts::generic::impl::charge_cut, p, -1);

      if (p->getRegion() == clas12::FD && detector_flags[clas12::FD]) {
        CHECK_CUT(ctx, "HTCC_nphe_cut", cuts::FD::impl::HTCC_nphe_cut, p);
        CHECK_CUT(ctx, "EC_outer_vs_EC_inner_cut", cuts::FD::impl::EC_outer_vs_EC_inner_cut, p, tightness);
        CHECK_CUT(ctx, "EC_sampling_fraction_cut", cuts::FD::impl::EC_sampling_fraction_cut, p, inbending, false,
                  false);
        CHECK_CUT(ctx, "EC_hit_position_fiducial_cut_homogeneous",
                  cuts::FD::impl::EC_hit_position_fiducial_cut_homogeneous, p, cuts::tightness::loose, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::DC_fiducial_cut_edge_reg3, p, inbending);
        CHECK_CUT(ctx, "DC_z_vertex_cut", cuts::FD::impl::DC_z_vertex_cut, p, inbending);
        // CHECK_CUT(ctx, "is_in_FD_check", []([[maybe_unused]] auto p) { return true; }, true); // just collect
        // statistics on FD counts
      }
      if (p->getRegion() == clas12::CD && detector_flags[clas12::CD]) {
        // no cuts specified
      }
      if (p->getRegion() == clas12::FT && detector_flags[clas12::FT]) {
        CHECK_CUT(ctx, "FT_eid_FTCAL_fiducial_cut", cuts::FT::impl::FT_eid_FTCAL_fiducial_cut, p);
        CHECK_CUT(ctx, "FT_eid_FTTRK_fiducial_cut", cuts::FT::impl::FT_eid_FTTRK_fiducial_cut, p);
        CHECK_CUT(ctx, "FT_eid_FTHODO_fiducial_cut", cuts::FT::impl::FT_eid_FTHODO_fiducial_cut, p);
        CHECK_CUT(ctx, "FT_eid_energy_vs_radius_cut", cuts::FT::impl::FT_eid_energy_vs_radius_cut, p);
      }
      return ctx.passed();
    }

    bool proton(ParticleContext& ctx, clas12::region_particle* p, bool inbending, cuts::tightness tightness,
                double reference_vertex_z) {
      // CHECK_CUT(ctx, "PID_cut(2212)", cuts::generic::impl::PID_cut, p, 2212);
      // CHECK_CUT(ctx, "charge_cut(+1)", cuts::generic::impl::charge_cut, p, +1);

      if (p->getRegion() == clas12::FD && detector_flags[clas12::FD]) {
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::DC_fiducial_cut_edge_reg3, p, inbending);
        CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::delta_vz_cut, p, reference_vertex_z);
      }
      if (p->getRegion() == clas12::CD && detector_flags[clas12::CD]) {
        // no cuts specified
      }
      if (p->getRegion() == clas12::FT && detector_flags[clas12::FT]) {
        // no cuts specified
      }
      return ctx.passed();
    }

    bool neutron(ParticleContext& ctx, clas12::region_particle* p, bool inbending, cuts::tightness tightness,
                 double reference_vertex_z) {
      // CHECK_CUT(ctx, "PID_cut(2112)", cuts::generic::impl::PID_cut, p, 2112);
      // CHECK_CUT(ctx, "charge_cut(0)", cuts::generic::impl::charge_cut, p, 0);
      CHECK_CUT(ctx, "momentum_cut", cuts::generic::impl::momentum_cut, p);

      if (p->getRegion() == clas12::FD && detector_flags[clas12::FD]) {
        CHECK_CUT(ctx, "neutr_beta_cut", cuts::impl::neutr_beta_cut, p);
        CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::delta_vz_cut, p, reference_vertex_z);
      }
      if (p->getRegion() == clas12::CD && detector_flags[clas12::CD]) {
        CHECK_CUT(ctx, "neutr_beta_cut", cuts::impl::neutr_beta_cut, p);
      }
      if (p->getRegion() == clas12::FT && detector_flags[clas12::FT]) {
        // no cuts specified
      }
      return ctx.passed();
    }

    bool piplus(ParticleContext& ctx, clas12::region_particle* p, bool inbending, cuts::tightness tightness,
                double reference_vertex_z) {
      // if (!ctx.apply_cut("PID_cut(211)", cuts::generic::impl::PID_cut, p, 211)) return ctx.passed();
      // CHECK_CUT(ctx, "charge_cut(+1)", cuts::generic::impl::charge_cut, p, +1);

      if (p->getRegion() == clas12::FD && detector_flags[clas12::FD]) {
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::DC_fiducial_cut_edge_reg3, p, inbending);
        CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::delta_vz_cut, p, reference_vertex_z);
      }
      if (p->getRegion() == clas12::CD && detector_flags[clas12::CD]) {
        // no cuts specified
      }
      if (p->getRegion() == clas12::FT && detector_flags[clas12::FT]) {
        // no cuts specified
      }
      return ctx.passed();
    }

    bool piminus(ParticleContext& ctx, clas12::region_particle* p, bool inbending, cuts::tightness tightness,
                 double reference_vertex_z) {
      // if (!ctx.apply_cut("PID_cut(-211)", cuts::generic::impl::PID_cut, p, -211)) return ctx.passed();
      // CHECK_CUT(ctx, "charge_cut(-1)", cuts::generic::impl::charge_cut, p, -1);

      if (p->getRegion() == clas12::FD && detector_flags[clas12::FD]) {
        // CHECK_CUT(ctx, "EC_outer_vs_EC_inner_cut", cuts::FD::impl::EC_outer_vs_EC_inner_cut, p, tightness); //
        // ignored in old filter
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::DC_fiducial_cut_edge_reg3, p, inbending);
        CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::delta_vz_cut, p, reference_vertex_z);
      }
      if (p->getRegion() == clas12::CD && detector_flags[clas12::CD]) {
        // no cuts specified
      }
      if (p->getRegion() == clas12::FT && detector_flags[clas12::FT]) {
        // no cuts specified
      }
      return ctx.passed();
    }

    bool Kplus(ParticleContext& ctx, clas12::region_particle* p, bool inbending, cuts::tightness tightness,
               double reference_vertex_z) {
      // CHECK_CUT(ctx, "PID_cut(321)", cuts::generic::impl::PID_cut, p, 321);
      // CHECK_CUT(ctx, "charge_cut(+1)", cuts::generic::impl::charge_cut, p, +1);

      if (p->getRegion() == clas12::FD && detector_flags[clas12::FD]) {
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::DC_fiducial_cut_edge_reg3, p, inbending);
        CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::delta_vz_cut, p, reference_vertex_z);
      }
      if (p->getRegion() == clas12::CD && detector_flags[clas12::CD]) {
        // no cuts specified
      }
      if (p->getRegion() == clas12::FT && detector_flags[clas12::FT]) {
        // no cuts specified
      }
      return ctx.passed();
    }

    bool Kminus(ParticleContext& ctx, clas12::region_particle* p, bool inbending, cuts::tightness tightness,
                double reference_vertex_z) {
      // CHECK_CUT(ctx, "PID_cut(-321)", cuts::generic::impl::PID_cut, p, -321);
      // CHECK_CUT(ctx, "charge_cut(-1)", cuts::generic::impl::charge_cut, p, -1);

      if (p->getRegion() == clas12::FD && detector_flags[clas12::FD]) {
        // CHECK_CUT(ctx, "EC_outer_vs_EC_inner_cut", cuts::FD::impl::EC_outer_vs_EC_inner_cut, p, tightness); //
        // ignored in old filter
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::DC_fiducial_cut_edge_reg3, p, inbending);
        CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::delta_vz_cut, p, reference_vertex_z);
      }
      if (p->getRegion() == clas12::CD && detector_flags[clas12::CD]) {
        // no cuts specified
      }
      if (p->getRegion() == clas12::FT && detector_flags[clas12::FT]) {
        // no cuts specified
      }
      return ctx.passed();
    }

    bool photon(ParticleContext& ctx, clas12::region_particle* p, bool inbending, cuts::tightness tightness) {
      // CHECK_CUT(ctx, "PID_cut(22)", cuts::generic::impl::PID_cut, p, 22);
      // CHECK_CUT(ctx, "charge_cut(0)", cuts::generic::impl::charge_cut, p, 0);

      if (p->getRegion() == clas12::FD && detector_flags[clas12::FD]) {
        CHECK_CUT(ctx, "phot_beta_cut", cuts::impl::phot_beta_cut, p, cuts::tightness::medium);
        CHECK_CUT(ctx, "EC_hit_position_fiducial_cut_homogeneous",
                  cuts::FD::impl::EC_hit_position_fiducial_cut_homogeneous, p, cuts::tightness::medium, inbending);
        // rest of cuts are computed but ignored in old filter
        // CHECK_CUT(ctx, "phot_EC_sampling_fraction_cut", cuts::FD::impl::phot_EC_sampling_fraction_cut, p);
        // CHECK_CUT(ctx, "phot_EC_outer_vs_EC_inner_cut", cuts::FD::impl::phot_EC_outer_vs_EC_inner_cut, p);
      }
      if (p->getRegion() == clas12::CD && detector_flags[clas12::CD]) {
        // no cuts specified
      }
      if (p->getRegion() == clas12::FT && detector_flags[clas12::FT]) {
        CHECK_CUT(ctx, "FT_photid_FTCAL_fiducial_cut", cuts::FT::impl::FT_photid_FTCAL_fiducial_cut, p);
        CHECK_CUT(ctx, "FT_photid_beta_cut", cuts::FT::impl::FT_photid_beta_cut, p);
      }
      return ctx.passed();
    }
  }  // namespace impl

  // Helper function to record cut history to statistics
  inline void record_selector_cuts(const std::string& selector_name, const ParticleContext& ctx, bool result) {
    StatisticsCollector::record_selector_invocation(selector_name, result);
    for (const auto& cut : ctx.get_cut_history()) {
      StatisticsCollector::record_cut_in_selector(selector_name, cut.name, cut.passed);
    }
  }

  // Public interface functions that create ParticleContext internally
  // These match the calling interface in new_filter.cxx
  inline bool electron(clas12::region_particle* p, bool inbending, cuts::tightness tightness) {
    ParticleContext ctx(evaluation_mode);
    bool result = impl::electron(ctx, p, inbending, tightness);
    record_selector_cuts("electron", ctx, result);
    return result;
  }

  inline bool proton(clas12::region_particle* p, bool inbending, cuts::tightness tightness, double reference_vertex_z) {
    ParticleContext ctx(evaluation_mode);
    bool result = impl::proton(ctx, p, inbending, tightness, reference_vertex_z);
    record_selector_cuts("proton", ctx, result);
    return result;
  }

  inline bool neutron(clas12::region_particle* p, bool inbending, cuts::tightness tightness,
                      double reference_vertex_z) {
    ParticleContext ctx(evaluation_mode);
    bool result = impl::neutron(ctx, p, inbending, tightness, reference_vertex_z);
    record_selector_cuts("neutron", ctx, result);
    return result;
  }

  inline bool piplus(clas12::region_particle* p, bool inbending, cuts::tightness tightness, double reference_vertex_z) {
    ParticleContext ctx(evaluation_mode);
    bool result = impl::piplus(ctx, p, inbending, tightness, reference_vertex_z);
    record_selector_cuts("piplus", ctx, result);
    return result;
  }

  inline bool piminus(clas12::region_particle* p, bool inbending, cuts::tightness tightness,
                      double reference_vertex_z) {
    ParticleContext ctx(evaluation_mode);
    bool result = impl::piminus(ctx, p, inbending, tightness, reference_vertex_z);
    record_selector_cuts("piminus", ctx, result);
    return result;
  }

  inline bool Kplus(clas12::region_particle* p, bool inbending, cuts::tightness tightness, double reference_vertex_z) {
    ParticleContext ctx(evaluation_mode);
    bool result = impl::Kplus(ctx, p, inbending, tightness, reference_vertex_z);
    record_selector_cuts("Kplus", ctx, result);
    return result;
  }

  inline bool Kminus(clas12::region_particle* p, bool inbending, cuts::tightness tightness, double reference_vertex_z) {
    ParticleContext ctx(evaluation_mode);
    bool result = impl::Kminus(ctx, p, inbending, tightness, reference_vertex_z);
    record_selector_cuts("Kminus", ctx, result);
    return result;
  }

  inline bool photon(clas12::region_particle* p, bool inbending, cuts::tightness tightness, double reference_vertex_z) {
    ParticleContext ctx(evaluation_mode);
    bool result = impl::photon(ctx, p, inbending, tightness);
    record_selector_cuts("photon", ctx, result);
    return result;
  }
}  // namespace selectors
