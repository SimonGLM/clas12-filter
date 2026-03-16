#pragma once

#include "region_particle.h"

namespace cuts {

  enum tightness { loose = 1, medium, tight };

  namespace generic {
    namespace impl {
      bool _forward_detector_cut(clas12::region_particle*);
      bool _central_detector_cut(clas12::region_particle*);
      bool _forward_tagger_cut(clas12::region_particle*);
      bool _PID_cut(clas12::region_particle*, int pid);
      bool _charge_cut(clas12::region_particle*, int charge);
    }  // namespace impl

    inline bool (*forward_detector_cut)(clas12::region_particle*) = cuts::generic::impl::_forward_detector_cut;
    inline bool (*central_detector_cut)(clas12::region_particle*) = cuts::generic::impl::_central_detector_cut;
    inline bool (*forward_tagger_cut)(clas12::region_particle*) = cuts::generic::impl::_forward_tagger_cut;
    inline bool (*PID_cut)(clas12::region_particle*, int) = cuts::generic::impl::_PID_cut;
    inline bool (*charge_cut)(clas12::region_particle*, int) = cuts::generic::impl::_charge_cut;
  }  // namespace generic

  namespace FD {
    namespace impl {
      bool _HTCC_nphe_cut(clas12::region_particle*);
      bool _EC_sampling_fraction_cut(clas12::region_particle*, bool inbending, bool simulation, bool spring2019);
      bool _EC_hit_position_fiducial_cut_homogeneous(clas12::region_particle*, tightness tightness, bool inbending);
      bool _EC_outer_vs_EC_inner_cut(clas12::region_particle*, tightness tightness);
      bool _DC_fiducial_cut_edge(clas12::region_particle*, int region, bool inbending);
      bool _DC_z_vertex_cut(clas12::region_particle*, bool);
      bool _phot_EC_sampling_fraction_cut(clas12::region_particle*);
      bool _phot_EC_outer_vs_EC_inner_cut(clas12::region_particle*);

      bool _DC_fiducial_cut_edge_reg1(clas12::region_particle* p, bool inbending) {
        return _DC_fiducial_cut_edge(p, 1, inbending);
      };
      bool _DC_fiducial_cut_edge_reg2(clas12::region_particle* p, bool inbending) {
        return _DC_fiducial_cut_edge(p, 2, inbending);
      };
      bool _DC_fiducial_cut_edge_reg3(clas12::region_particle* p, bool inbending) {
        return _DC_fiducial_cut_edge(p, 3, inbending);
      };
    }  // namespace impl

    inline bool (*HTCC_nphe_cut)(clas12::region_particle*) = cuts::FD::impl::_HTCC_nphe_cut;
    inline bool (*EC_sampling_fraction_cut)(clas12::region_particle*, bool, bool,
                                            bool) = cuts::FD::impl::_EC_sampling_fraction_cut;
    inline bool (*EC_hit_position_fiducial_cut_homogeneous)(clas12::region_particle*, tightness, bool) =
        cuts::FD::impl::_EC_hit_position_fiducial_cut_homogeneous;
    inline bool (*EC_outer_vs_EC_inner_cut)(clas12::region_particle*,
                                            tightness) = cuts::FD::impl::_EC_outer_vs_EC_inner_cut;
    inline bool (*DC_fiducial_cut_edge_region1)(clas12::region_particle*,
                                                bool) = cuts::FD::impl::_DC_fiducial_cut_edge_reg1;
    inline bool (*DC_fiducial_cut_edge_region2)(clas12::region_particle*,
                                                bool) = cuts::FD::impl::_DC_fiducial_cut_edge_reg2;
    inline bool (*DC_fiducial_cut_edge_region3)(clas12::region_particle*,
                                                bool) = cuts::FD::impl::_DC_fiducial_cut_edge_reg3;
    inline bool (*DC_z_vertex_cut)(clas12::region_particle*, bool) = cuts::FD::impl::_DC_z_vertex_cut;
    inline bool (*phot_EC_sampling_fraction_cut)(clas12::region_particle*) =
        cuts::FD::impl::_phot_EC_sampling_fraction_cut;
    inline bool (*phot_EC_outer_vs_EC_inner_cut)(clas12::region_particle*) =
        cuts::FD::impl::_phot_EC_outer_vs_EC_inner_cut;
  }  // namespace FD

  namespace FT {
    namespace impl {
      bool _FT_photid_FTCAL_check_hole(double X, double Y, double holeX, double holeY, double holeR) {
        pow(X - holeX, 2) + pow(Y - holeY, 2) - pow(holeR, 2) < 1;
      };
      bool _FT_eid_FTCAL_fiducial_cut(clas12::region_particle*);
      bool _FT_eid_FTTRK_fiducial_cut(clas12::region_particle*);
      bool _FT_eid_FTHODO_fiducial_cut(clas12::region_particle*);
      bool _FT_photid_FTCAL_fiducial_cut(clas12::region_particle*);
      bool _FT_eid_energy_vs_radius_cut(clas12::region_particle*);
      bool _FT_photid_beta_cut(clas12::region_particle*);
    }  // namespace impl

    inline bool (*FT_eid_FTCAL_fiducial_cut)(clas12::region_particle*) = cuts::FT::impl::_FT_eid_FTCAL_fiducial_cut;
    inline bool (*FT_eid_FTTRK_fiducial_cut)(clas12::region_particle*) = cuts::FT::impl::_FT_eid_FTTRK_fiducial_cut;
    inline bool (*FT_eid_FTHODO_fiducial_cut)(clas12::region_particle*) = cuts::FT::impl::_FT_eid_FTHODO_fiducial_cut;
    inline bool (*FT_photid_FTCAL_fiducial_cut)(clas12::region_particle*) =
        cuts::FT::impl::_FT_photid_FTCAL_fiducial_cut;
    inline bool (*FT_eid_energy_vs_radius_cut)(clas12::region_particle*) = cuts::FT::impl::_FT_eid_energy_vs_radius_cut;
    inline bool (*FT_photid_beta_cut)(clas12::region_particle*) = cuts::FT::impl::_FT_photid_beta_cut;
  }  // namespace FT

  namespace CD {
    namespace impl {
      [[maybe_unused]] bool _CD_neutr_beta_cut(clas12::region_particle*);
    }
    [[maybe_unused]] inline bool (*CD_neutr_beta_cut)(clas12::region_particle*) = cuts::CD::impl::_CD_neutr_beta_cut;
  }  // namespace CD

  namespace vertex {
    namespace impl {
      bool _delta_vz_cut(clas12::region_particle*, double reference_vertex_z);
    }  // namespace impl
    inline bool (*delta_vz_cut)(clas12::region_particle*, double) = cuts::vertex::impl::_delta_vz_cut;
  }  // namespace vertex

  namespace impl {
    bool _phot_beta_cut(clas12::region_particle*, tightness tightness);
    bool _neutr_beta_cut(clas12::region_particle*);
    bool _basic_FTOF_cut(clas12::region_particle*);
  }  // namespace impl

  inline bool (*phot_beta_cut)(clas12::region_particle*, tightness) = cuts::impl::_phot_beta_cut;
  inline bool (*neutr_beta_cut)(clas12::region_particle*) = cuts::impl::_neutr_beta_cut;
  inline bool (*basic_FTOF_cut)(clas12::region_particle*) = cuts::impl::_basic_FTOF_cut;
}  // namespace cuts
