#pragma once
#include "cuts.h"

namespace selectors {
  namespace impl {
    bool _electron(clas12::region_particle* p, bool inbending, int tightness) {
      // std::cout << std::format("Checking electron...") << std::endl;
      // basic electron selection

      if (!cuts::generic::PID_cut(p, 11)) return false;
      if (!cuts::generic::charge_cut(p, -1)) return false;

      if (p->getRegion() == clas12::FD) {
        if (!cuts::FD::HTCC_nphe_cut(p)) return false;
        if (!cuts::FD::EC_outer_vs_EC_inner_cut(p, tightness)) return false;
        if (!cuts::FD::EC_sampling_fraction_cut(p, inbending, false, false)) return false;
        if (!cuts::FD::EC_hit_position_fiducial_cut_homogeneous(p, tightness, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region1(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region2(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region3(p, inbending)) return false;
        if (!cuts::FD::DC_z_vertex_cut(p)) return false;
      }
      if (p->getRegion() == clas12::FT) {
        if (!cuts::FT::FT_eid_FTCAL_fiducial_cut(p)) return false;
        if (!cuts::FT::FT_eid_FTTRK_fiducial_cut(p)) return false;
        if (!cuts::FT::FT_eid_FTHODO_fiducial_cut(p)) return false;
        if (!cuts::FT::FT_eid_energy_vs_radius_cut(p)) return false;
      }
      if (p->getRegion() == clas12::CD) {
        // no filters specified
      }
      // std::cout << std::format("Valid electron.") << std::endl;
      return true;
    }

    bool _proton(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
      // std::cout << std::format("Checking proton...") << std::endl;
      if (!cuts::generic::PID_cut(p, 2212)) return false;
      if (!cuts::generic::charge_cut(p, +1)) return false;

      if (p->getRegion() == clas12::FD) {
        if (!cuts::FD::DC_fiducial_cut_edge_region1(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region2(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region3(p, inbending)) return false;
        if (!cuts::vertex::delta_vz_cut(p, reference_vertex_z)) return false;
      }
      if (p->getRegion() == clas12::FT) {
        // no filters specified
      }
      if (p->getRegion() == clas12::CD) {
        // no filters specified
      }

      // std::cout << std::format("Valid proton.") << std::endl;
      return true;
    }

    bool _neutron(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
      // std::cout << std::format("Checking neutron...") << std::endl;
      if (!cuts::generic::PID_cut(p, 2112)) return false;
      if (!cuts::generic::charge_cut(p, 0)) return false;

      if (p->getRegion() == clas12::FD) {
        if (!cuts::neutr_beta_cut(p, 0 /*run*/)) return false;
        if (!cuts::vertex::delta_vz_cut(p, reference_vertex_z)) return false;
      }
      if (p->getRegion() == clas12::FT) {
        // no filters specified
      }
      if (p->getRegion() == clas12::CD) {
        if (!cuts::neutr_beta_cut(p, 0 /*run*/)) return false;
      }

      // std::cout << std::format("Valid neutron.") << std::endl;
      return true;
    }

    bool _piplus(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
      // std::cout << std::format("Checking piplus...") << std::endl;
      // basic pion selection
      if (!cuts::generic::PID_cut(p, 221 && !cuts::generic::PID_cut(p, -221))) return false;
      if (!cuts::generic::charge_cut(p, +1)) return false;

      if (p->getRegion() == clas12::FD) {
        if (!cuts::FD::DC_fiducial_cut_edge_region1(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region2(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region3(p, inbending)) return false;
        if (!cuts::vertex::delta_vz_cut(p, reference_vertex_z)) return false;
      }
      if (p->getRegion() == clas12::FT) {
        // no filters specified
      }
      if (p->getRegion() == clas12::CD) {
        // no filters specified
      }

      // std::cout << std::format("Valid piplus.") << std::endl;
      return true;
    }

    bool _piminus(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
      // std::cout << std::format("Checking piminus...") << std::endl;
      // basic pion selection
      if (!cuts::generic::PID_cut(p, 221 && !cuts::generic::PID_cut(p, -221))) return false;
      if (!cuts::generic::charge_cut(p, -1)) return false;
      // still missing: pim_ele_reject_cut (checks if electron passed)

      if (p->getRegion() == clas12::FD) {
        if (!cuts::FD::EC_outer_vs_EC_inner_cut(p, tightness)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region1(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region2(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region3(p, inbending)) return false;
        if (!cuts::vertex::delta_vz_cut(p, reference_vertex_z)) return false;
      }
      if (p->getRegion() == clas12::FT) {
        // no filters specified
      }
      if (p->getRegion() == clas12::CD) {
        // no filters specified
      }

      // std::cout << std::format("Valid piminus.") << std::endl;
      return true;
    }

    bool _Kplus(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
      // std::cout << std::format("Checking Kplus...") << std::endl;
      if (!cuts::generic::PID_cut(p, 321)) return false;
      if (!cuts::generic::charge_cut(p, +1)) return false;
      if (!cuts::vertex::delta_vz_cut(p, reference_vertex_z)) return false;

      if (p->getRegion() == clas12::FD) {
        if (!cuts::FD::DC_fiducial_cut_edge_region1(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region2(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region3(p, inbending)) return false;
      }
      if (p->getRegion() == clas12::FT) {
        // no filters specified
      }
      if (p->getRegion() == clas12::CD) {
        // no filters specified
      }

      // std::cout << std::format("Valid Kplus.") << std::endl;
      return true;
    }

    bool _Kminus(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
      // std::cout << std::format("Checking Kminus...") << std::endl;
      if (!cuts::generic::PID_cut(p, 321)) return false;
      if (!cuts::generic::charge_cut(p, +1)) return false;
      // still missing: pim_ele_reject_cut (checks if electron passed)

      if (p->getRegion() == clas12::FD) {
        if (!cuts::FD::EC_outer_vs_EC_inner_cut(p, tightness)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region1(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region2(p, inbending)) return false;
        if (!cuts::FD::DC_fiducial_cut_edge_region3(p, inbending)) return false;
        if (!cuts::vertex::delta_vz_cut(p, 0.)) return false;
      }
      if (p->getRegion() == clas12::FT) {
        // no filters specified
      }
      if (p->getRegion() == clas12::CD) {
        if (!cuts::vertex::delta_vz_cut(p, reference_vertex_z)) return false;
      }

      // std::cout << std::format("Valid Kminus.") << std::endl;
      return true;
    }

    bool _photon(clas12::region_particle* p, bool inbending, int tightness, double reference_vertex_z) {
      // std::cout << std::format("Checking photon...") << std::endl;
      if (!cuts::generic::PID_cut(p, 22)) return false;
      if (!cuts::generic::charge_cut(p, 0)) return false;
      // still missing: pim_ele_reject_cut (checks if electron passed)

      if (p->getRegion() == clas12::FD) {
        if (!cuts::phot_beta_cut(p, tightness)) return false;
        if (!cuts::FD::phot_EC_outer_vs_EC_inner_cut(p)) return false;
        // reference code calls phot_EC_hit_position_fiducial_cut, but they are identical, bar one override for sector 4
        if (!cuts::FD::EC_hit_position_fiducial_cut_homogeneous(p, tightness, inbending)) return false;
      }
      if (p->getRegion() == clas12::FT) {
        // no filters specified
        if (!cuts::FT::FT_photid_FTCAL_fiducial_cut(p)) return false;
        if (!cuts::FT::FT_photid_beta_cut(p)) return false;
      }
      if (p->getRegion() == clas12::CD) {
        if (!cuts::vertex::delta_vz_cut(p, reference_vertex_z)) return false;
      }

      // std::cout << std::format("Valid photon.") << std::endl;
      return true;
    }
  }  // namespace impl
  static DecoratedCut<bool (*)(clas12::region_particle*, bool, int)> electron("selector::electron",
                                                                              selectors::impl::_electron);
  static DecoratedCut<bool (*)(clas12::region_particle*, bool, int, double)> proton("selector::proton",
                                                                                    selectors::impl::_proton);
  static DecoratedCut<bool (*)(clas12::region_particle*, bool, int, double)> neutron("selector::neutron",
                                                                                     selectors::impl::_neutron);
  static DecoratedCut<bool (*)(clas12::region_particle*, bool, int, double)> piplus("selector::piplus",
                                                                                    selectors::impl::_piplus);
  static DecoratedCut<bool (*)(clas12::region_particle*, bool, int, double)> piminus("selector::piminus",
                                                                                     selectors::impl::_piminus);
  static DecoratedCut<bool (*)(clas12::region_particle*, bool, int, double)> Kplus("selector::Kplus",
                                                                                   selectors::impl::_Kplus);
  static DecoratedCut<bool (*)(clas12::region_particle*, bool, int, double)> Kminus("selector::Kminus",
                                                                                    selectors::impl::_Kminus);
  static DecoratedCut<bool (*)(clas12::region_particle*, bool, int, double)> photon("selector::photon",
                                                                                    selectors::impl::_photon);
}  // namespace selectors