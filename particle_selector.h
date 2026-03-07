#pragma once
#include "cuts.h"
#include "particle_context.h"

// Helper macro to check cut and continue/stop based on mode
#define CHECK_CUT(ctx, name, func, ...) \
  if (!ctx.apply_cut(name, func, __VA_ARGS__)) return ctx.passed()

namespace selectors {
  namespace impl {
    bool _electron(ParticleContext& ctx, clas12::region_particle* p, bool inbending, int tightness) {
      CHECK_CUT(ctx, "PID_cut(11)", cuts::generic::impl::_PID_cut, p, 11);
      CHECK_CUT(ctx, "charge_cut(-1)", cuts::generic::impl::_charge_cut, p, -1);

      if (p->getRegion() == clas12::FD) {
        CHECK_CUT(ctx, "HTCC_nphe_cut", cuts::FD::impl::_HTCC_nphe_cut, p);
        CHECK_CUT(ctx, "EC_outer_vs_EC_inner_cut", cuts::FD::impl::_EC_outer_vs_EC_inner_cut, p, tightness);
        CHECK_CUT(ctx, "EC_sampling_fraction_cut", cuts::FD::impl::_EC_sampling_fraction_cut, p, inbending, false,
                  false);
        CHECK_CUT(ctx, "EC_hit_position_fiducial_cut_homogeneous",
                  cuts::FD::impl::_EC_hit_position_fiducial_cut_homogeneous, p, tightness, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::_DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::_DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::_DC_fiducial_cut_edge_reg3, p, inbending);
        CHECK_CUT(ctx, "DC_z_vertex_cut", cuts::FD::impl::_DC_z_vertex_cut, p);
      }
      if (p->getRegion() == clas12::FT) {
        CHECK_CUT(ctx, "FT_eid_FTCAL_fiducial_cut", cuts::FT::impl::_FT_eid_FTCAL_fiducial_cut, p);
        CHECK_CUT(ctx, "FT_eid_FTTRK_fiducial_cut", cuts::FT::impl::_FT_eid_FTTRK_fiducial_cut, p);
        CHECK_CUT(ctx, "FT_eid_FTHODO_fiducial_cut", cuts::FT::impl::_FT_eid_FTHODO_fiducial_cut, p);
        CHECK_CUT(ctx, "FT_eid_energy_vs_radius_cut", cuts::FT::impl::_FT_eid_energy_vs_radius_cut, p);
      }
      return ctx.passed();
    }

    bool _proton(ParticleContext& ctx, clas12::region_particle* p, bool inbending, int tightness,
                 double reference_vertex_z) {
      CHECK_CUT(ctx, "PID_cut(2212)", cuts::generic::impl::_PID_cut, p, 2212);
      CHECK_CUT(ctx, "charge_cut(+1)", cuts::generic::impl::_charge_cut, p, +1);

      if (p->getRegion() == clas12::FD) {
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::_DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::_DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::_DC_fiducial_cut_edge_reg3, p, inbending);
        CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::_delta_vz_cut, p, reference_vertex_z);
      }
      return ctx.passed();
    }

    bool _neutron(ParticleContext& ctx, clas12::region_particle* p, bool inbending, int tightness,
                  double reference_vertex_z) {
      CHECK_CUT(ctx, "PID_cut(2112)", cuts::generic::impl::_PID_cut, p, 2112);
      CHECK_CUT(ctx, "charge_cut(0)", cuts::generic::impl::_charge_cut, p, 0);

      if (p->getRegion() == clas12::FD) {
        CHECK_CUT(ctx, "neutr_beta_cut", cuts::impl::_neutr_beta_cut, p, 0);
        CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::_delta_vz_cut, p, reference_vertex_z);
      }
      if (p->getRegion() == clas12::CD) {
        CHECK_CUT(ctx, "neutr_beta_cut", cuts::impl::_neutr_beta_cut, p, 0);
      }
      return ctx.passed();
    }

    bool _piplus(ParticleContext& ctx, clas12::region_particle* p, bool inbending, int tightness,
                 double reference_vertex_z) {
      if (!ctx.apply_cut("PID_cut(221)", cuts::generic::impl::_PID_cut, p, 221)) return ctx.passed();
      CHECK_CUT(ctx, "charge_cut(+1)", cuts::generic::impl::_charge_cut, p, +1);

      if (p->getRegion() == clas12::FD) {
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::_DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::_DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::_DC_fiducial_cut_edge_reg3, p, inbending);
        CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::_delta_vz_cut, p, reference_vertex_z);
      }
      return ctx.passed();
    }

    bool _piminus(ParticleContext& ctx, clas12::region_particle* p, bool inbending, int tightness,
                  double reference_vertex_z) {
      if (!ctx.apply_cut("PID_cut(221)", cuts::generic::impl::_PID_cut, p, 221)) return ctx.passed();
      CHECK_CUT(ctx, "charge_cut(-1)", cuts::generic::impl::_charge_cut, p, -1);

      if (p->getRegion() == clas12::FD) {
        CHECK_CUT(ctx, "EC_outer_vs_EC_inner_cut", cuts::FD::impl::_EC_outer_vs_EC_inner_cut, p, tightness);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::_DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::_DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::_DC_fiducial_cut_edge_reg3, p, inbending);
        CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::_delta_vz_cut, p, reference_vertex_z);
      }
      return ctx.passed();
    }

    bool _Kplus(ParticleContext& ctx, clas12::region_particle* p, bool inbending, int tightness,
                double reference_vertex_z) {
      CHECK_CUT(ctx, "PID_cut(321)", cuts::generic::impl::_PID_cut, p, 321);
      CHECK_CUT(ctx, "charge_cut(+1)", cuts::generic::impl::_charge_cut, p, +1);
      CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::_delta_vz_cut, p, reference_vertex_z);

      if (p->getRegion() == clas12::FD) {
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::_DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::_DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::_DC_fiducial_cut_edge_reg3, p, inbending);
      }
      return ctx.passed();
    }

    bool _Kminus(ParticleContext& ctx, clas12::region_particle* p, bool inbending, int tightness,
                 double reference_vertex_z) {
      CHECK_CUT(ctx, "PID_cut(321)", cuts::generic::impl::_PID_cut, p, 321);
      CHECK_CUT(ctx, "charge_cut(-1)", cuts::generic::impl::_charge_cut, p, -1);
      CHECK_CUT(ctx, "delta_vz_cut", cuts::vertex::impl::_delta_vz_cut, p, reference_vertex_z);

      if (p->getRegion() == clas12::FD) {
        CHECK_CUT(ctx, "EC_outer_vs_EC_inner_cut", cuts::FD::impl::_EC_outer_vs_EC_inner_cut, p, tightness);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region1", cuts::FD::impl::_DC_fiducial_cut_edge_reg1, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region2", cuts::FD::impl::_DC_fiducial_cut_edge_reg2, p, inbending);
        CHECK_CUT(ctx, "DC_fiducial_cut_edge_region3", cuts::FD::impl::_DC_fiducial_cut_edge_reg3, p, inbending);
      }
      return ctx.passed();
    }

    bool _photon(ParticleContext& ctx, clas12::region_particle* p, bool inbending, int tightness) {
      CHECK_CUT(ctx, "PID_cut(22)", cuts::generic::impl::_PID_cut, p, 22);
      CHECK_CUT(ctx, "charge_cut(0)", cuts::generic::impl::_charge_cut, p, 0);

      if (p->getRegion() == clas12::FD) {
        CHECK_CUT(ctx, "phot_EC_sampling_fraction_cut", cuts::FD::impl::_phot_EC_sampling_fraction_cut, p);
        CHECK_CUT(ctx, "phot_EC_outer_vs_EC_inner_cut", cuts::FD::impl::_phot_EC_outer_vs_EC_inner_cut, p);
      }
      return ctx.passed();
    }

    bool _electron_forward(ParticleContext& ctx, clas12::region_particle* p, bool inbending, int tightness) {
      return _electron(ctx, p, inbending, tightness);
    }
  }  // namespace impl

  // Public interface functions that create ParticleContext internally
  // These match the calling interface in new_filter.cxx
  inline bool electron(clas12::region_particle* p, bool inbending, int tightness) {
    ParticleContext ctx(EvaluationMode::EarlyReturn);
    return impl::_electron(ctx, p, inbending, tightness);
  }

  inline bool proton(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
    ParticleContext ctx(EvaluationMode::EarlyReturn);
    return impl::_proton(ctx, p, inbending, tightness, reference_vertex_z);
  }

  inline bool neutron(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
    ParticleContext ctx(EvaluationMode::EarlyReturn);
    return impl::_neutron(ctx, p, inbending, tightness, reference_vertex_z);
  }

  inline bool piplus(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
    ParticleContext ctx(EvaluationMode::EarlyReturn);
    return impl::_piplus(ctx, p, inbending, tightness, reference_vertex_z);
  }

  inline bool piminus(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
    ParticleContext ctx(EvaluationMode::EarlyReturn);
    return impl::_piminus(ctx, p, inbending, tightness, reference_vertex_z);
  }

  inline bool Kplus(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
    ParticleContext ctx(EvaluationMode::EarlyReturn);
    return impl::_Kplus(ctx, p, inbending, tightness, reference_vertex_z);
  }

  inline bool Kminus(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
    ParticleContext ctx(EvaluationMode::EarlyReturn);
    return impl::_Kminus(ctx, p, inbending, tightness, reference_vertex_z);
  }

  inline bool photon(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
    ParticleContext ctx(EvaluationMode::EarlyReturn);
    return impl::_photon(ctx, p, inbending, tightness);
  }
}  // namespace selectors
