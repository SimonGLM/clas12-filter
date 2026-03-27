#pragma once

#include "region_particle.h"

namespace cuts {

  enum tightness { loose = 1, medium, tight };

  namespace generic {
    namespace impl {
      bool forward_detector_cut(clas12::region_particle*);
      bool central_detector_cut(clas12::region_particle*);
      bool forward_tagger_cut(clas12::region_particle*);
      bool PID_cut(clas12::region_particle*, int pid);
      bool charge_cut(clas12::region_particle*, int charge);
      bool momentum_cut(clas12::region_particle*);
    }  // namespace impl

    inline bool (*forward_detector_cut)(clas12::region_particle*) = cuts::generic::impl::forward_detector_cut;
    inline bool (*central_detector_cut)(clas12::region_particle*) = cuts::generic::impl::central_detector_cut;
    inline bool (*forward_tagger_cut)(clas12::region_particle*) = cuts::generic::impl::forward_tagger_cut;
    inline bool (*PID_cut)(clas12::region_particle*, int) = cuts::generic::impl::PID_cut;
    inline bool (*charge_cut)(clas12::region_particle*, int) = cuts::generic::impl::charge_cut;
    inline bool (*momentum_cut)(clas12::region_particle*) = cuts::generic::impl::momentum_cut;
  }  // namespace generic

  namespace FD {
    namespace impl {
      bool HTCC_nphe_cut(clas12::region_particle*);
      bool EC_sampling_fraction_cut(clas12::region_particle*, bool inbending, bool simulation, bool spring2019);
      bool EC_hit_position_fiducial_cut_homogeneous(clas12::region_particle*, tightness tightness, bool inbending);
      bool EC_outer_vs_EC_inner_cut(clas12::region_particle*, tightness tightness);
      bool DC_fiducial_cut_edge(clas12::region_particle*, int region, bool inbending);
      bool DC_z_vertex_cut(clas12::region_particle*, bool);
      bool phot_EC_sampling_fraction_cut(clas12::region_particle*);
      bool phot_EC_outer_vs_EC_inner_cut(clas12::region_particle*);

      bool DC_fiducial_cut_edge_reg1(clas12::region_particle* p, bool inbending) {
        return impl::DC_fiducial_cut_edge(p, 1, inbending);
      };
      bool DC_fiducial_cut_edge_reg2(clas12::region_particle* p, bool inbending) {
        return impl::DC_fiducial_cut_edge(p, 2, inbending);
      };
      bool DC_fiducial_cut_edge_reg3(clas12::region_particle* p, bool inbending) {
        return impl::DC_fiducial_cut_edge(p, 3, inbending);
      };
    }  // namespace impl

    inline bool (*HTCC_nphe_cut)(clas12::region_particle*) = cuts::FD::impl::HTCC_nphe_cut;
    inline bool (*EC_sampling_fraction_cut)(clas12::region_particle*, bool, bool,
                                            bool) = cuts::FD::impl::EC_sampling_fraction_cut;
    inline bool (*EC_hit_position_fiducial_cut_homogeneous)(clas12::region_particle*, tightness, bool) =
        cuts::FD::impl::EC_hit_position_fiducial_cut_homogeneous;
    inline bool (*EC_outer_vs_EC_inner_cut)(clas12::region_particle*,
                                            tightness) = cuts::FD::impl::EC_outer_vs_EC_inner_cut;
    inline bool (*DC_fiducial_cut_edge_region1)(clas12::region_particle*,
                                                bool) = cuts::FD::impl::DC_fiducial_cut_edge_reg1;
    inline bool (*DC_fiducial_cut_edge_region2)(clas12::region_particle*,
                                                bool) = cuts::FD::impl::DC_fiducial_cut_edge_reg2;
    inline bool (*DC_fiducial_cut_edge_region3)(clas12::region_particle*,
                                                bool) = cuts::FD::impl::DC_fiducial_cut_edge_reg3;
    inline bool (*DC_z_vertex_cut)(clas12::region_particle*, bool) = cuts::FD::impl::DC_z_vertex_cut;
    inline bool (*phot_EC_sampling_fraction_cut)(clas12::region_particle*) =
        cuts::FD::impl::phot_EC_sampling_fraction_cut;
    inline bool (*phot_EC_outer_vs_EC_inner_cut)(clas12::region_particle*) =
        cuts::FD::impl::phot_EC_outer_vs_EC_inner_cut;
  }  // namespace FD

  namespace FT {
    namespace impl {
      bool FT_photid_FTCAL_check_hole(double X, double Y, double holeX, double holeY, double holeR) {
        return pow(X - holeX, 2) + pow(Y - holeY, 2) - pow(holeR, 2) < 1;
      };
      bool FT_eid_FTCAL_fiducial_cut(clas12::region_particle*);
      bool FT_eid_FTTRK_fiducial_cut(clas12::region_particle*);
      bool FT_eid_FTHODO_fiducial_cut(clas12::region_particle*);
      bool FT_photid_FTCAL_fiducial_cut(clas12::region_particle*);
      bool FT_eid_energy_vs_radius_cut(clas12::region_particle*);
      bool FT_photid_beta_cut(clas12::region_particle*);
    }  // namespace impl

    inline bool (*FT_eid_FTCAL_fiducial_cut)(clas12::region_particle*) = cuts::FT::impl::FT_eid_FTCAL_fiducial_cut;
    inline bool (*FT_eid_FTTRK_fiducial_cut)(clas12::region_particle*) = cuts::FT::impl::FT_eid_FTTRK_fiducial_cut;
    inline bool (*FT_eid_FTHODO_fiducial_cut)(clas12::region_particle*) = cuts::FT::impl::FT_eid_FTHODO_fiducial_cut;
    inline bool (*FT_photid_FTCAL_fiducial_cut)(clas12::region_particle*) =
        cuts::FT::impl::FT_photid_FTCAL_fiducial_cut;
    inline bool (*FT_eid_energy_vs_radius_cut)(clas12::region_particle*) = cuts::FT::impl::FT_eid_energy_vs_radius_cut;
    inline bool (*FT_photid_beta_cut)(clas12::region_particle*) = cuts::FT::impl::FT_photid_beta_cut;
  }  // namespace FT

  namespace CD {
    namespace impl {
      [[maybe_unused]] bool CD_neutr_beta_cut(clas12::region_particle*);
    }
    [[maybe_unused]] inline bool (*CD_neutr_beta_cut)(clas12::region_particle*) = cuts::CD::impl::CD_neutr_beta_cut;
  }  // namespace CD

  namespace vertex {
    namespace impl {
      bool delta_vz_cut(clas12::region_particle*, double reference_vertex_z);
    }  // namespace impl
    inline bool (*delta_vz_cut)(clas12::region_particle*, double) = cuts::vertex::impl::delta_vz_cut;
  }  // namespace vertex

  namespace impl {
    bool phot_beta_cut(clas12::region_particle*, tightness tightness);
    bool neutr_beta_cut(clas12::region_particle*);
    bool basic_FTOF_cut(clas12::region_particle*);
  }  // namespace impl

  inline bool (*phot_beta_cut)(clas12::region_particle*, tightness) = cuts::impl::phot_beta_cut;
  inline bool (*neutr_beta_cut)(clas12::region_particle*) = cuts::impl::neutr_beta_cut;
  inline bool (*basic_FTOF_cut)(clas12::region_particle*) = cuts::impl::basic_FTOF_cut;
}  // namespace cuts
