#pragma once
#include "cuts.h"

#include <Math/Vector3D.h>
#include <clas12defs.h>
#include <region_particle.h>

#include <numbers>

#include "cut_parameters.h"
#include "helpers.h"

using bounds = cuts::parameters::bounds;

namespace cuts {
  namespace generic::impl {
    bool _forward_detector_cut(clas12::region_particle* p) {
      int det = abs(p->par()->getStatus());
      // clas12::FD<= det && det < clas12::CD;  // enum missmatch?
      return (2000 <= det && det < 4000);
    }

    bool _central_detector_cut(clas12::region_particle* p) {
      int det = abs(p->par()->getStatus());
      return (4000 <= det && det < 5000);
    }

    bool _forward_tagger_cut(clas12::region_particle* p) {
      int det = abs(p->par()->getStatus());
      return (1000 <= det && det < 2000);
    }

    bool _PID_cut(clas12::region_particle* p, int pid) { return p->par()->getPid() == pid; }

    bool _charge_cut(clas12::region_particle* p, int charge) { return p->par()->getCharge() == charge; }

    // bool vertex_cut(clas12::region_particle* p, int runnum) {
    //   // from tbhayward/clas12_analysis_software
    //   int charge = p->par()->getCharge();
    //   float vz = p->par()->getVz();
    //
    //   //                          pos     neg   charge
    //   using bounds_pair = std::pair<bounds, bounds>;
    //   bounds_pair bounds_spring18_inb = {{-7.879, 1.515}, {-6.0606, 1.8182}};
    //   bounds_pair bounds_spring18_outb = {{-6.6667, 2.7273}, {-7.273, 0.9091}};
    //   bounds_pair bounds_fall18_inb = {{-8.485, 0.606}, {-6.364, 1.515}};
    //   bounds_pair bounds_fall18_outb = {{-6.970, 1.818}, {-7.879, 0.303}};
    //   bounds_pair bounds_spring19_inb = bounds_fall18_inb;
    //   bounds_pair bounds_summer22 = {{-9.394, -0.606}, {-7.576, 0.303}};
    //   bounds_pair bounds_fall22_spring23 = {{-8.788, 0.303}, {-5.758, 1.515}};
    //   bounds fallback = {-9, 2};
    //
    //   bounds_pair vz_bounds_pair;
    //   if (runnum == 11) {
    //     vz_bounds_pair = {{-10, 1.5}, {-9, 2}};
    //   } else if ((runnum >= 3173 && runnum >= 3293) || (runnum >= 3863 && runnum <= 3987))  // spring18 outb
    //   {
    //     vz_bounds_pair = bounds_spring18_outb;
    //   } else if (runnum >= 4003 && runnum <= 4325)  // spring18 inb
    //   {
    //     vz_bounds_pair = bounds_spring18_inb;
    //   } else if (runnum >= 5032 && runnum <= 5419)  // fall18 inbending
    //   {
    //     vz_bounds_pair = bounds_fall18_inb;
    //   } else if (runnum >= 5422 && runnum <= 5666)  // fall18 outbending
    //   {
    //     vz_bounds_pair = bounds_fall18_outb;
    //   } else if (runnum >= 6616 && runnum <= 6783)  // spring19 inbending
    //   {
    //     vz_bounds_pair = bounds_spring19_inb;
    //   } else if (runnum >= 16043 && runnum <= 16772)  // su22
    //   {
    //     vz_bounds_pair = bounds_summer22;
    //   } else if (runnum >= 16843 && runnum <= 17811)  // fa22 & sp32
    //   {
    //     vz_bounds_pair = bounds_fall22_spring23;
    //   }
    //
    //   bounds vz_bounds = charge > 0   ? vz_bounds_pair.first
    //                      : charge < 0 ? vz_bounds_pair.second
    //                                   : fallback;  // choose from charge or fallback for neutral
    //
    //   return vz_bounds.lower < vz && vz < vz_bounds.upper;
    // }
  }  // namespace generic::impl

  namespace FD::impl {
    bool _HTCC_nphe_cut(clas12::region_particle* p) {
      // original CC_nphe_cut
      double nphe_min = 2;
      return p->che(clas12::HTCC)->getNphe() > nphe_min;
    }

    // TODO: move boolean arguments to infering from region_particle
    bool _EC_sampling_fraction_cut(clas12::region_particle* p, bool inbending = false, bool simulation = false,
                                   bool spring19 = false) {
      // original EC_sampling_fraction_cut
      //////////////////////////////////////////////////////////////////////////////
      // band cut on SF PCAL (array entries are sectors)
      //////////////////////////////////////////////////////////////////////////////

      using namespace cuts::parameters::EC_sampling_fraction;

      // select the right LUT for the current conditions
      const band::BandParameterLUT* band_LUT;
      if (spring19)
        band_LUT = &band::spring2019;
      else if (simulation) {
        band_LUT = inbending ? &band::simulation_inb : &band::simulation_outb;
      } else {
        band_LUT = inbending ? &band::fall2018_inb : &band::fall2018_outb;
      }

      // helper variables to make the formula more readable
      int sector = p->cal(clas12::PCAL)->getSector() - 1;
      double P = p->par()->getP();
      double total_ECAL_energy =
          p->cal(clas12::PCAL)->getEnergy() + p->cal(clas12::ECIN)->getEnergy() + p->cal(clas12::ECOUT)->getEnergy();
      double sampling_fraction = total_ECAL_energy / P;

      [[maybe_unused]] double sigma_range = 3.5;

      // calculate band
      double mean = band_LUT->p0.mean[sector] * (1 + P / std::sqrt(P * P + band_LUT->p1.mean[sector])) +
                    band_LUT->p2.mean[sector] * P + band_LUT->p3.mean[sector] * P * P;
      double sigma = band_LUT->p0.sigma[sector] + band_LUT->p1.sigma[sector] / std::sqrt(P) +
                     band_LUT->p2.sigma[sector] * P + band_LUT->p3.sigma[sector] * P * P;
      // calulate cut limits
      bounds limits{.lower = mean - sigma_range * sigma, .upper = mean + sigma_range * sigma};

      bool pass_band = limits.lower <= sampling_fraction && sampling_fraction <= limits.upper;

      // std::cout << "SF band cut: " << limits.lower << " < " << sampling_fraction << " < " << limits.upper
      //           << " == " << pass_band << std::endl;
      //////////////////////////////////////////////////////////////////////////////
      // triangle cut on SF PCAL vs SF ECin (array entries are momentum bins)
      //////////////////////////////////////////////////////////////////////////////

      const triangle::TriangleParameterLUT* tri_LUT;
      if (spring19) {
        tri_LUT = &triangle::spring19;
      } else if (simulation) {
        tri_LUT = inbending ? &triangle::simulation_inb : &triangle::simulation_outb;
      } else {
        tri_LUT = inbending ? &triangle::fall18_inb : &triangle::fall18_outb;
      }

      int energy_bin = P <= 3 ? 0 : P <= 4 ? 1 : P <= 5 ? 2 : P <= 6 ? 3 : P <= 7 ? 4 : P <= 8 ? 5 : P <= 9 ? 6 : 7;
      int sec = p->cal(clas12::PCAL)->getSector() - 1;

      double total_EC_energy =
          p->cal(clas12::PCAL)->getEnergy() + p->cal(clas12::ECIN)->getEnergy() + p->cal(clas12::ECOUT)->getEnergy();

      double p0 = tri_LUT->sectors[sec]->p0[energy_bin];
      double p1 = tri_LUT->sectors[sec]->p1[energy_bin];

      double triangle_cut_min = (p1 - p0 * p->cal(clas12::ECIN)->getEnergy() / P);
      bool pass_triangle = triangle_cut_min < sampling_fraction;
      // std::cout << "SF triangle cut: " << triangle_cut_min << " < " << sampling_fraction << " == " << pass_triangle
      //           << std::endl;

      //////////////////////////////////////////////////////////////////////////////
      // threshold cut on SF PCAL
      //////////////////////////////////////////////////////////////////////////////
      bool pass_threshold = (p->cal(clas12::PCAL)->getEnergy() / P) > threshold::value;

      // // final decision

      if (pass_band) ++EC_SF_band;
      if (pass_triangle) ++EC_SF_triangle;
      if (pass_threshold) ++EC_SF_threshold;

      return pass_band && pass_triangle && pass_threshold;
    }

    bool _EC_hit_position_fiducial_cut_homogeneous(clas12::region_particle* p, tightness tightness, bool inbending) {
      // original EC_hit_position_fiducial_cut_homogeneous

      if (tightness != cuts::loose and tightness != cuts::medium and tightness != cuts::tight)
        throw std::runtime_error(
            "[EC_hit_position_fiducial_cut_homogeneous] tightness must be cuts::tightness::loose, "
            "cuts::tightness::medium, or cuts::tightness::tight.");
      // Cut using the natural directions of the scintillator bars/ fibers:
      double v = p->cal(clas12::PCAL)->getLv();
      double w = p->cal(clas12::PCAL)->getLw();
      [[maybe_unused]] double u = p->cal(clas12::PCAL)->getLu();

      /// v + w is going from the side to the back end of the PCAL, u is going from
      /// side to side 1 scintillator bar is 4.5 cm wide. In the outer regions
      /// (back) double bars are used. a cut is only applied on v and w

      struct vw_bounds {
        bounds v;
        bounds w;
        // double min_w;
        // double max_w;
      };

      using vw_bounds_LUT = std::vector<std::vector<vw_bounds>>;

      const vw_bounds_LUT LIMITS_LUT = {// inbending
                                        {
                                            // Loose
                                            {{9.0, 400.0}, {9.0, 400.0}},
                                            // Medium
                                            {{14.0, 400.0}, {14.0, 400.0}},
                                            // Tight
                                            {{19.0, 400.0}, {19.0, 400.0}},
                                        },
                                        // outbending
                                        {
                                            // Loose
                                            {{9.0, 400.0}, {9.0, 400.0}},
                                            // Medium
                                            {{14.0, 400.0}, {14.0, 400.0}},
                                            // Tight
                                            {{19.0, 400.0}, {19.0, 400.0}},
                                        }};

      // overrides for sector #4 (idx=3)
      const vw_bounds LOOSE_SECTOR4_OVERRIDES = {{13.5, 400.0}, {9.0, 400.0}};

      // select the right bounds from the LUT
      vw_bounds cut_bounds = LIMITS_LUT[inbending ? 0 : 1][tightness - 1];

      // std::cout << "v bounds: " << cut_bounds.v.lower << " < " << v << " < " << cut_bounds.v.upper
      //           << " == " << (cut_bounds.v.lower < v && v <= cut_bounds.v.upper) << std::endl;
      // std::cout << "w bounds: " << cut_bounds.w.lower << " < " << w << " < " << cut_bounds.w.upper
      //           << " == " << (cut_bounds.w.lower < w && w <= cut_bounds.w.upper) << std::endl;
      // Override for sector 4 with loose tightness (TODO: reason unkown) (not for photons)
      if (tightness == cuts::loose && p->cal(clas12::PCAL)->getSector() == 4 && p->par()->getPid() != 22) {
        return (LOOSE_SECTOR4_OVERRIDES.v.lower < v && v < LOOSE_SECTOR4_OVERRIDES.v.upper &&
                LOOSE_SECTOR4_OVERRIDES.w.lower < w && w < LOOSE_SECTOR4_OVERRIDES.w.upper);
      }

      return (cut_bounds.v.lower < v && v < cut_bounds.v.upper && cut_bounds.w.lower < w && w < cut_bounds.w.upper);
    }

    bool _EC_outer_vs_EC_inner_cut(clas12::region_particle* p, tightness tightness) {
      // original *_EC_outer_vs_EC_inner_cut
      if (tightness != cuts::loose and tightness != cuts::medium and tightness != cuts::tight)
        throw std::runtime_error(
            "[EC_outer_vs_EC_inner_cut] tightness must be cuts::tightness::loose, cuts::tightness::medium, or "
            "cuts::tightness::tight.");

      std::array<double, 3> edep_min = {0.06, 0.07, 0.09};
      // std::cout << "edep inner: " << p->cal(clas12::PCAL)->getEnergy() << std::endl;
      // std::cout << "edep_min for tightness " << tightness << ": " << edep_min[tightness - 1] << std::endl;

      return p->cal(clas12::PCAL)->getEnergy() > edep_min[tightness - 1];
    }

    bool _DC_fiducial_cut_edge(clas12::region_particle* p, int region, bool inbending) {
      // original DC_fiducial_cut_edge
      using DCEdgeCuts = struct {
        std::array<double, 3> inbending;
        std::array<double, 3> outbending;
      };
      const std::map<int, DCEdgeCuts> DCedge_LUT = {
          {11, {{5.0, 5.0, 10.0}, {3.0, 3.0, 10.0}}}, {2212, {{2.5, 2.5, 9.0}, {3.5, 3.0, 7.0}}},
          {211, {{2.5, 2.5, 9.0}, {3.5, 2.5, 6.5}}},  {-211, {{3.5, 3.0, 7.0}, {2.5, 2.5, 10.0}}},
          {321, {{2.5, 2.0, 9.0}, {3.5, 2.5, 6.5}}},  {-321, {{3.5, 2.5, 5.0}, {2.5, 2.5, 10.0}}}};

      if (!cuts::generic::impl::_PID_cut(p, 11) && !cuts::generic::impl::_PID_cut(p, 2212) &&
          !cuts::generic::impl::_PID_cut(p, 211) && !cuts::generic::impl::_PID_cut(p, -211) &&
          !cuts::generic::impl::_PID_cut(p, 321) && !cuts::generic::impl::_PID_cut(p, -321)) {
        throw std::domain_error(
            std::format("[DC_fiducial_cut_edge] Attempting cut on invalid PID '{}' for this cut.", p->getPid()));
      }
      // get iterator to arrays in the DCedge_LUT
      auto it = DCedge_LUT.find(p->getPid());
      if (it == DCedge_LUT.end())
        throw std::runtime_error("[DC_fiducial_cut_edge] AHhhhh, Panic! This should not be reachable.");
      // select edge cut based on inbending or outbending
      double edge_cut = inbending ? it->second.inbending[region - 1] : it->second.outbending[region - 1];

      double edge_val = 0;
      switch (region) {
        case 1:
          edge_val = p->traj(clas12::DC, clas12::DC1)->getEdge();
          break;
        case 2:
          edge_val = p->traj(clas12::DC, clas12::DC3)->getEdge();
          break;
        case 3:
          edge_val = p->traj(clas12::DC, clas12::DC6)->getEdge();
          break;
        default:
          throw std::domain_error("This particle's region \'" + std::to_string(region) +
                                  " is not 1,2 or 3. Check what getRegion really returns in code.");
          break;
      }

      return edge_val > edge_cut;
    }

    bool _DC_z_vertex_cut(clas12::region_particle* p, bool inbending) {
      int sector = p->cal(clas12::PCAL)->getSector();

      bounds vz_bounds_inb = {.lower = -8, .upper = 2};
      bounds vz_bounds_outb = {.lower = -11, .upper = 1};

      if (inbending) {
        return vz_bounds_inb.lower < p->par()->getVz() && p->par()->getVz() < vz_bounds_inb.upper;
      } else {
        return vz_bounds_outb.lower < p->par()->getVz() && p->par()->getVz() < vz_bounds_outb.upper;
      }
    }

    bool _phot_EC_sampling_fraction_cut(clas12::region_particle* p) {
      bounds limits = {0., 1.};

      double pcal = p->cal(clas12::PCAL)->getEnergy();
      double ecin = p->cal(clas12::ECIN)->getEnergy();
      double ecout = p->cal(clas12::ECOUT)->getEnergy();
      double total = pcal + ecin + ecout;

      double cutvalue = total / p->par()->getP();
      return (limits.lower < cutvalue && cutvalue < limits.upper);
    }

    bool _phot_EC_outer_vs_EC_inner_cut(clas12::region_particle* p) {
      double edep_min = 0.01;
      return (p->cal(clas12::ECIN)->getEnergy() + p->cal(clas12::ECOUT)->getEnergy()) > edep_min;
      // throw not_implemented_error("[phot_EC_outer_vs_EC_inner_cut] Not implemented yet.");
    }
  }  // namespace FD::impl

  namespace FT::impl {

    bool _FT_eid_FTCAL_fiducial_cut(clas12::region_particle* p) {
      double theta = std::acos(p->par()->getPz() / p->par()->getP()) * 180 / std::numbers::pi;

      return theta > 2.5 && theta < 4.5;
    }

    bool _FT_eid_FTTRK_fiducial_cut(clas12::region_particle* p) {
      // same as FT_eid_FTCAL_fiducial_cut
      return _FT_eid_FTCAL_fiducial_cut(p);
    }

    bool _FT_eid_FTHODO_fiducial_cut(clas12::region_particle* p) {
      // same as FT_eid_FTCAL_fiducial_cut
      return _FT_eid_FTCAL_fiducial_cut(p);
    }

    bool _FT_eid_energy_vs_radius_cut([[maybe_unused]] clas12::region_particle* _) { return true; }

    bool _FT_photid_FTCAL_fiducial_cut(clas12::region_particle* p) {
      double clusX = p->ft(clas12::FTCAL)->getX();
      double clusY = p->ft(clas12::FTCAL)->getY();
      ROOT::Math::XYZVector V3ECalPos(clusX, clusY, 0);

      bool res = true && V3ECalPos.R() > 8 && V3ECalPos.R() < 15 &&
                 std::pow(clusX + 8.5, 2) + std::pow(clusY - 10, 2) > 1.5 * 1.5 &&
                 std::pow(clusX + 10, 2) + std::pow(clusY + 5, 2) > 1.5 * 1.5 &&
                 std::pow(clusX + 6, 2) + std::pow(clusY + 13.5, 2) > 2 * 2 &&
                 std::pow(clusX - 4, 2) + std::pow(clusY + 6.7, 2) > 1.5 * 1.5 &&
                 std::pow(clusX - 6, 2) + std::pow(clusY + 6, 2) > 1;

      return res;
    }

    bool _FT_photid_beta_cut([[maybe_unused]] clas12::region_particle* _) { return true; }

  }  // namespace FT::impl

  namespace CD::impl {
    bool _CD_neutr_beta_cut(clas12::region_particle* p) { return cuts::impl::_neutr_beta_cut(p); }
  }  // namespace CD::impl

  namespace vertex::impl {

    // We need a good reference vertex z for the *_delta_vz_cuts.
    // The old way is searching for the highest momentum electron in the event
    // each time this cut is called. How can we make this smarter?
    // For now: we require a reference vertex and pass the problem of finding that to the future self.

    // <particle>_delta_vz_cuts and CD_<particle>_delta_vz_cuts have different parameters but in the end
    // they are stale code and have no effect. It is just checked if the reference vertex z is between -20 and 20.

    bool _delta_vz_cut(clas12::region_particle* p, double reference_vertex_z) {
      bounds tolerances{-20, 20};
      double delta_vz = reference_vertex_z - p->par()->getVz();
      return tolerances.lower < delta_vz && delta_vz < tolerances.upper;
    }

    // implementation for these are all shadowed by hardcoded bounds
    bool _CD_delta_vz_cut(clas12::region_particle* p, double reference_vertex_z) {
      throw not_implemented_error("[CD_delta_vz_cut] This cut is deprecated. Use idenical 'delta_vz_cut'.");
      // // p m=0.7486, s=3.237, min=-20, max=20
      // // n m=2.254, s=2.693, min=-20, max=20
      // // Pp m=-0.6183, s=3.684, min=-20, max=20
      // // Pm m=-0.5485, s=3.677, min=-20, max=20
      // // Kp m=1.658, s=2.52, min=-20, max=20
      // // Km m=-1.161, s=2.691, min=-20, max=20

      // struct ParamMeanStd {
      //   double mean;
      //   double sigma;
      // };
      // [[maybe_unused]] std::map<int, ParamMeanStd> params_per_pid = {
      //     {2212, {.mean = 0.7486, .sigma = 3.237}}, {2112, {.mean = 2.254, .sigma = 2.693}},
      //     {211, {.mean = -0.6183, .sigma = 3.684}}, {-211, {.mean = -0.5485, .sigma = 3.677}},
      //     {321, {.mean = 1.658, .sigma = 2.52}},    {-321, {.mean = -1.161, .sigma = 2.691}},
      // };

      // [[maybe_unsued]] int pid = p->par()->getPid();
      // [[maybe_unsued]] double mean = params_per_pid.at(pid).mean;
      // [[maybe_unsued]] double sigma = params_per_pid.at(pid).sigma;
      // [[maybe_unsued]] double dvz_min = mean - 3 * sigma;
      // [[maybe_unsued]] double dvz_max = mean + 3 * sigma;
      // [[maybe_unsued]] bounds dvz_bounds = {dvz_min, dvz_max};

      // // double delta_vz = p->par()->getVz(); // Is this the proper reference vertex?

      // bounds tolerances{-20, 20};
      // return tolerances.lower < reference_vertex_z && reference_vertex_z < tolerances.upper;
    }
  }  // namespace vertex::impl

  namespace impl {
    bool _phot_beta_cut(clas12::region_particle* p, tightness tightness) {
      if (tightness != cuts::loose and tightness != cuts::medium and tightness != cuts::tight)
        throw std::runtime_error(
            "[phot_beta_cut] tightness must be cuts::tightness::loose, cuts::tightness::medium, or "
            "cuts::tightness::tight.");

      std::array<bounds, 3> beta_cut_LUT = {{
          {0.9, 2.0},   // loose
          {0.9, 1.1},   // medium
          {0.95, 1.05}  // tight
      }};

      float beta = p->par()->getBeta();
      return beta_cut_LUT[tightness - 1].lower < beta && beta < beta_cut_LUT[tightness - 1].upper &&
             p->par()->getP() > 0.10;
    }

    bool _neutr_beta_cut(clas12::region_particle* p) {
      // very similar to phot_beta_cut
      bounds limits = {0., 0.95};
      return limits.lower < p->par()->getBeta() && p->par()->getBeta() < limits.upper;
    }

    bool _basic_FTOF_cut(clas12::region_particle* p) {
      return p->sci(clas12::FTOF1A)->getSector() != 0 || p->sci(clas12::FTOF1B)->getSector() != 0 ||
             p->sci(clas12::FTOF2)->getSector() != 0;  // equal to "is somewhere in FTOF?", I guess...
    }
  }  // namespace impl
}  // namespace cuts