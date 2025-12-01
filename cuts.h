#include <clas12defs.h>
#include <region_particle.h>

#include <numbers>

class not_implemented_error : public std::logic_error {
 public:
  explicit not_implemented_error(const std::string& what_arg) : std::logic_error(what_arg) {}
};

struct bounds {
  double upper;
  double lower;
};

namespace cuts {
bool PID_cut(clas12::region_particle* p, int pid) { return p->par()->getPid() == pid; }

bool charge_cut(clas12::region_particle* p, int charge) { return p->par()->getCharge() == charge; }

// We need a good reference vertex z for the proton delta vz cut.
// The old way is searching for the highest momentum electron in the event
// each time this cut is called. How can we make this smarter?
// For now: we require a reference vertex and pass the problem of finding that to the future self.
bool delta_vz_cut(clas12::region_particle* p, double reference_vertex_z) {
  bounds tolerances{-20, 20};
  double delta_vz = p->par()->getVz();
  return tolerances.lower < delta_vz && delta_vz < tolerances.upper;
}

////////////////////////////////////////////////////////////////////////////////
//                      ELECTRON CUTS (FORWARD DETECTOR)                      //
////////////////////////////////////////////////////////////////////////////////

bool CC_nphe_cut(clas12::region_particle* p) {
  double nphe_min = 2;
  return p->che(clas12::HTCC)->getNphe() > nphe_min;
}

bool EC_outer_vs_EC_inner_cut(clas12::region_particle* p, int tightness) {
  if (tightness > 1 or tightness < 3)
    throw std::runtime_error("[EC_outer_vs_EC_inner_cut] tightness must be 1, 2, or 3.");

  std::array<double, 3> edep_min = {0.06, 0.07, .09};

  return p->cal(clas12::PCAL)->getEnergy() > edep_min[tightness - 1];
}

// TODO: move boolean arguments to infering from region_particle
// TODO: Maybe remove LUT in favour to RCDB or CCDB, if possible ?
bool EC_sampling_fraction_cut(clas12::region_particle* p, bool inbending = false, bool simulation = false,
                              bool spring2019 = false) {

  //////////////////////////////////////////////////////////////////////////////
  // band cut on SF PCAL (array entries are sectors)
  //////////////////////////////////////////////////////////////////////////////

  using ParamterAllSectors = std::array<float, 6>;
  struct ParamMeanStdAllSectors {
    ParamterAllSectors mean;
    ParamterAllSectors sigma;
  };
  using BandParameterLUT = std::array<ParamMeanStdAllSectors, 4>;
  ParamMeanStdAllSectors band_p0, band_p1, band_p2, band_p3;

  // fall2018 inb:
  band_p0 = {.mean = {0.111767, 0.116619, 0.114606, 0.116586, 0.118251, 0.117391},
             .sigma = {-0.00497609, 0.0259435, 0.0296159, 0.0161445, 0.0239166, 0.0244309}};
  band_p1 = {.mean = {-0.0281943, 0.0662751, -0.0896597, 0.181465, 0.085993, 0.0186504},
             .sigma = {0.0275006, -0.000805156, -0.00449379, 0.0099462, 0.00192551, 0.00258059}};
  band_p2 = {.mean = {0.00711137, 0.00633334, 0.00912098, 0.00652068, 0.00416682, 0.00622289},
             .sigma = {0.00253641, -0.00386759, -0.00469883, -0.00182968, -0.00355973, -0.00398967}};
  band_p3 = {.mean = {-0.000878776, -0.000780257, -0.00108891, -0.0006957, -0.000485189, -0.000829729},
             .sigma = {-0.000173549, 0.00030325, 0.000380195, 0.00012328, 0.000302528, 0.000340911}};
  BandParameterLUT band_fall2018_inb = {band_p0, band_p1, band_p2, band_p3};

  // fall2018 outb:
  band_p0 = {.mean = {0.111919, 0.11244, 0.11457, 0.124517, 0.109132, 0.115026},
             .sigma = {-0.000828514, 0.019356, 0.023144, -0.000468566, 0.00500942, -0.00167471}};
  band_p1 = {.mean = {-0.00764253, 0.156704, 0.246338, 0.880436, -0.181137, 0.335205},
             .sigma = {0.023998, 0.00249064, -0.00118722, 0.0223872, 0.0167932, 0.0243664}};
  band_p2 = {.mean = {0.00937217, 0.00924749, 0.00931085, 0.000420224, 0.00935877, 0.00756135},
             .sigma = {0.000966279, -0.0022336, -0.00284786, 0.00143044, 0.000355477, 0.00124825}};
  band_p3 = {.mean = {-0.000948603, -0.00095255, -0.000993028, -0.000123388, -0.000825096, -0.000822597},
             .sigma = {-6.99914e-05, 0.000152666, 0.000164229, -0.000133009, -3.68797e-05, -0.000107538}};
  BandParameterLUT band_fall2018_outb = {band_p0, band_p1, band_p2, band_p3};

  //   // spring2019 inb:
  band_p0 = {.mean = {0.11253, 0.113735, 0.112401, 0.115128, 0.113048, 0.1147},
             .sigma = {0.0193473, 0.0351352, 0.0234448, 0.0238342, 0.0382829, 0.0125166}};
  band_p1 = {.mean = {-0.0689836, -0.044216, -0.160555, 0.108512, -0.153003, -0.0997027},
             .sigma = {0.00436399, -0.0100767, 0.00133137, 0.00193097, -0.0127937, 0.0139726}};
  band_p2 = {.mean = {0.00793526, 0.00835112, 0.010781, 0.00695328, 0.00814596, 0.00766677},
             .sigma = {-0.00222507, -0.00549204, -0.00353928, -0.00347335, -0.00623149, -0.00186593}};
  band_p3 = {.mean = {-0.000920149, -0.000977549, -0.00128293, -0.000728472, -0.000957354, -0.00106035},
             .sigma = {0.000165603, 0.000440174, 0.000274543, 0.000260158, 0.000508396, 0.000206116}};
  BandParameterLUT band_spring2019_inb = {band_p0, band_p1, band_p2, band_p3};

  //   // MC inb:
  band_p0 = {.mean = {0.118444, 0.118383, 0.118318, 0.118531, 0.117475, 0.119179},
             .sigma = {0.0204537, 0.0242836, 0.0320663, 0.0171258, 0.0236728, 0.0157762}};
  band_p1 = {.mean = {-0.0445042, -0.0326496, 0.00402908, 0.0384926, -0.0768068, 0.0045002},
             .sigma = {0.00334198, -0.000192939, -0.00718171, 0.00698743, 0.000511942, 0.00801015}};
  band_p2 = {.mean = {0.00522281, 0.00528503, 0.00542581, 0.00480286, 0.00601034, 0.00482612},
             .sigma = {-0.00299872, -0.00369136, -0.00537603, -0.00246959, -0.00371171, -0.00215882}};
  band_p3 = {.mean = {-0.000511299, -0.000504437, -0.000529496, -0.000449033, -0.000578208, -0.00046328},
             .sigma = {0.000232288, 0.000277151, 0.000425544, 0.00018372, 0.000290388, 0.000160679}};
  BandParameterLUT band_simulation_inb = {band_p0, band_p1, band_p2, band_p3};

  //   // MC outb:
  band_p0 = {.mean = {0.124075, 0.124086, 0.124071, 0.125947, 0.120091, 0.124457},
             .sigma = {0.00231891, 0.000686535, 0.000327404, -0.000165373, 0.00376051, 0.000856615}};
  band_p1 = {.mean = {0.259169, 0.257774, 0.286948, 0.449492, -0.0180208, 0.274625},
             .sigma = {0.0205562, 0.0222621, 0.0226224, 0.0233338, 0.0190044, 0.0220325}};
  band_p2 = {.mean = {0.00203391, 0.00204442, 0.00199255, 0.00074766, 0.00414372, 0.00178793},
             .sigma = {0.000189274, 0.00044114, 0.000487635, 0.000531129, -9.69477e-06, 0.000416313}};
  band_p3 = {.mean = {-0.000135587, -0.000134399, -0.000132802, -3.89331e-05, -0.000283801, -0.000114508},
             .sigma = {-3.11787e-05, -4.52271e-05, -4.89526e-05, -5.22239e-05, -2.02358e-05, -4.52287e-05}};
  BandParameterLUT band_simulation_outb = {band_p0, band_p1, band_p2, band_p3};
  // // Example access:
  // int parameter=1;
  // band_fall2018_inb[parameter].mean[sector];

  BandParameterLUT* band_LUT;
  if (spring2019)
    band_LUT = &band_spring2019_inb;
  else if (simulation) {
    band_LUT = inbending ? &band_simulation_inb : &band_simulation_outb;
  } else {
    band_LUT = inbending ? &band_fall2018_inb : &band_fall2018_outb;
  }

  // helper variables to make the formula more readable
  const auto& band_p0 = band_LUT->at(0);
  const auto& band_p1 = band_LUT->at(1);
  const auto& band_p2 = band_LUT->at(2);
  const auto& band_p3 = band_LUT->at(3);
  int sector = p->cal(clas12::PCAL)->getSector() - 1;
  double P = p->par()->getP();

  double sigma_range = 3.5;

  // calculate band
  double mean = band_p0.mean[sector] * (1 + P / std::sqrt(P * P + band_p1.mean[sector]));
  double sigma = band_p0.sigma[sector] + band_p1.sigma[sector] / std::sqrt(P) + band_p2.sigma[sector] * P +
                 band_p3.sigma[sector] * P * P;
  // calulate cut limits
  double upper_lim_total = mean + 3.5 * sigma;
  double lower_lim_total = mean - 3.5 * sigma;

  double sampling_fraction = p->cal(clas12::PCAL)->getEnergy() / P;
  bool pass_band = lower_lim_total <= sampling_fraction && sampling_fraction <= upper_lim_total;

  //////////////////////////////////////////////////////////////////////////////
  // triangle cut on SF PCAL vs SF ECin (array entries are momentum bins)
  //////////////////////////////////////////////////////////////////////////////
  using ParamAllEnergies = std::array<float, 8>;
  struct TwoParamsAllEnergies {
    ParamAllEnergies p0;
    ParamAllEnergies p1;
  };
  using TriangleParameterLUT = std::array<TwoParamsAllEnergies, 6>;

  TwoParamsAllEnergies tri_sec1, tri_sec2, tri_sec3, tri_sec4, tri_sec5, tri_sec6;
  tri_sec1 = {.p0 = {1.41582, 1.39934, 1.41204, 1.46385, 1.55892, 1.55892, 1.55892, 1.55892},
              .p1 = {0.212225, 0.215542, 0.217, 0.218279, 0.219881, 0.219881, 0.219881, 0.219881}};
  tri_sec2 = {.p0 = {1.44726, 1.44245, 1.47269, 1.53225, 1.61465, 1.61465, 1.61465, 1.61465},
              .p1 = {0.221991, 0.225772, 0.227888, 0.229099, 0.228898, 0.228898, 0.228898, 0.228898}};
  tri_sec3 = {.p0 = {1.38589, 1.3908, 1.42501, 1.48177, 1.57636, 1.57636, 1.57636, 1.57636},
              .p1 = {0.221492, 0.225738, 0.227955, 0.228604, 0.22836, 0.22836, 0.22836, 0.22836}};
  tri_sec4 = {.p0 = {1.38631, 1.38107, 1.39757, 1.44579, 1.54154, 1.54154, 1.54154, 1.54154},
              .p1 = {0.215784, 0.221511, 0.224982, 0.227812, 0.231076, 0.231076, 0.231076, 0.231076}};
  tri_sec5 = {.p0 = {1.50251, 1.52408, 1.52996, 1.49583, 1.39339, 1.39339, 1.39339, 1.39339},
              .p1 = {0.22202, 0.227163, 0.228794, 0.226487, 0.218168, 0.218168, 0.218168, 0.218168}};
  tri_sec6 = {.p0 = {1.51312, 1.52784, 1.57519, 1.67332, 1.85128, 1.85128, 1.85128, 1.85128},
              .p1 = {0.223651, 0.228082, 0.2305, 0.23241, 0.234238, 0.234238, 0.234238, 0.234238}};
  const TriangleParameterLUT TRI_FALL2018_INB = {tri_sec1, tri_sec2, tri_sec3, tri_sec4, tri_sec5, tri_sec6};

  //   // fall2018 outb:
  tri_sec1 = {.p0 = {1.35967, 1.33697, 1.34111, 1.39563, 1.49066, 1.49066, 1.49066, 1.49066},
              .p1 = {0.21934, 0.222755, 0.224377, 0.227803, 0.23137, 0.23137, 0.23137, 0.23137}};
  tri_sec2 = {.p0 = {1.36974, 1.36895, 1.39344, 1.46945, 1.61251, 1.61251, 1.61251, 1.61251},
              .p1 = {0.218443, 0.223847, 0.227189, 0.231625, 0.238017, 0.238017, 0.238017, 0.238017}};
  tri_sec3 = {.p0 = {1.31891, 1.30602, 1.33372, 1.41351, 1.54453, 1.54453, 1.54453, 1.54453},
              .p1 = {0.219256, 0.223559, 0.22568, 0.230189, 0.235643, 0.235643, 0.235643, 0.235643}};
  tri_sec4 = {.p0 = {1.34425, 1.31016, 1.28753, 1.29671, 1.33831, 1.33831, 1.33831, 1.33831},
              .p1 = {0.217914, 0.221843, 0.221103, 0.221309, 0.222932, 0.222932, 0.222932, 0.222932}};
  tri_sec5 = {.p0 = {1.42433, 1.42395, 1.40932, 1.40816, 1.39868, 1.39868, 1.39868, 1.39868},
              .p1 = {0.218131, 0.222749, 0.225099, 0.225572, 0.224091, 0.224091, 0.224091, 0.224091}};
  tri_sec6 = {.p0 = {1.43741, 1.41924, 1.43218, 1.51807, 1.64554, 1.64554, 1.64554, 1.64554},
              .p1 = {0.220976, 0.225786, 0.228382, 0.232594, 0.238174, 0.238174, 0.238174, 0.238174}};
  const TriangleParameterLUT TRI_FALL2018_OUTB = {tri_sec1, tri_sec2, tri_sec3, tri_sec4, tri_sec5, tri_sec6};

  //   // spring 2019:
  tri_sec1 = {.p0 = {1.39979, 1.38131, 1.40144, 1.45945, 1.57144, 1.57144, 1.57144, 1.57144},
              .p1 = {0.21633, 0.219878, 0.222027, 0.223849, 0.22592, 0.22592, 0.22592, 0.22592}};
  tri_sec2 = {.p0 = {1.55243, 1.55871, 1.61819, 1.72756, 1.89424, 1.89424, 1.89424, 1.89424},
              .p1 = {0.223711, 0.228475, 0.231997, 0.235231, 0.236813, 0.236813, 0.236813, 0.236813}};
  tri_sec3 = {.p0 = {1.36528, 1.37519, 1.42453, 1.50339, 1.60863, 1.60863, 1.60863, 1.60863},
              .p1 = {0.22, 0.224385, 0.227207, 0.228813, 0.228142, 0.228142, 0.228142, 0.228142}};
  tri_sec4 = {.p0 = {1.38535, 1.3697, 1.39661, 1.4662, 1.63342, 1.63342, 1.63342, 1.63342},
              .p1 = {0.214364, 0.219152, 0.222319, 0.225052, 0.229045, 0.229045, 0.229045, 0.229045}};
  tri_sec5 = {.p0 = {1.47796, 1.47884, 1.48836, 1.47034, 1.38339, 1.38339, 1.38339, 1.38339},
              .p1 = {0.219635, 0.223343, 0.224231, 0.221401, 0.210005, 0.210005, 0.210005, 0.210005}};
  tri_sec6 = {.p0 = {1.51755, 1.53504, 1.59927, 1.73397, 2.00518, 2.00518, 2.00518, 2.00518},
              .p1 = {0.222143, 0.226195, 0.228638, 0.23086, 0.233709, 0.233709, 0.233709, 0.233709}};
  const TriangleParameterLUT TRI_SPRING2019 = {tri_sec1, tri_sec2, tri_sec3, tri_sec4, tri_sec5, tri_sec6};

  //   // mc inb:
  tri_sec1 = {.p0 = {1.30263, 1.30977, 1.31412, 1.31338, 1.31648, 1.31648, 1.31648, 1.31648},
              .p1 = {0.224593, 0.227935, 0.229518, 0.22978, 0.22924, 0.22924, 0.22924, 0.22924}};
  tri_sec2 = {.p0 = {1.30425, 1.30269, 1.30435, 1.30472, 1.29839, 1.29839, 1.29839, 1.29839},
              .p1 = {0.224099, 0.226717, 0.22775, 0.227332, 0.224101, 0.224101, 0.224101, 0.224101}};
  tri_sec3 = {.p0 = {1.28919, 1.28002, 1.24232, 1.18734, 1.12886, 1.12886, 1.12886, 1.12886},
              .p1 = {0.223284, 0.225601, 0.224208, 0.220147, 0.214082, 0.214082, 0.214082, 0.214082}};
  tri_sec4 = {.p0 = {1.29025, 1.29296, 1.29721, 1.30699, 1.32036, 1.32036, 1.32036, 1.32036},
              .p1 = {0.221678, 0.224968, 0.226676, 0.227668, 0.227833, 0.227833, 0.227833, 0.227833}};
  tri_sec5 = {.p0 = {1.2883, 1.27655, 1.24203, 1.19403, 1.15325, 1.15325, 1.15325, 1.15325},
              .p1 = {0.224075, 0.226291, 0.225091, 0.22173, 0.217695, 0.217695, 0.217695, 0.217695}};
  tri_sec6 = {.p0 = {1.30049, 1.30766, 1.31197, 1.31372, 1.31596, 1.31596, 1.31596, 1.31596},
              .p1 = {0.224591, 0.227901, 0.229503, 0.23, 0.229422, 0.229422, 0.229422, 0.229422}};
  const TriangleParameterLUT TRI_SIMULATION_INB = {tri_sec1, tri_sec2, tri_sec3, tri_sec4, tri_sec5, tri_sec6};

  //   // mc outb:

  tri_sec1 = {.p0 = {1.34088, 1.34685, 1.34914, 1.3544, 1.35628, 1.35628, 1.35628, 1.35628},
              .p1 = {0.229433, 0.232865, 0.234838, 0.236458, 0.237532, 0.237532, 0.237532, 0.237532}};
  tri_sec2 = {.p0 = {1.33865, 1.34439, 1.34938, 1.35247, 1.35423, 1.35423, 1.35423, 1.35423},
              .p1 = {0.229232, 0.232698, 0.234938, 0.236421, 0.237479, 0.237479, 0.237479, 0.237479}};
  tri_sec3 = {.p0 = {1.33981, 1.3419, 1.3459, 1.3506, 1.35108, 1.35108, 1.35108, 1.35108},
              .p1 = {0.228797, 0.232055, 0.234322, 0.235882, 0.236784, 0.236784, 0.236784, 0.236784}};
  tri_sec4 = {.p0 = {1.34856, 1.34544, 1.34355, 1.34317, 1.34383, 1.34383, 1.34383, 1.34383},
              .p1 = {0.22849, 0.231908, 0.233365, 0.234226, 0.235005, 0.235005, 0.235005, 0.235005}};
  tri_sec5 = {.p0 = {1.33882, 1.34416, 1.34414, 1.34681, 1.34737, 1.34737, 1.34737, 1.34737},
              .p1 = {0.228946, 0.232535, 0.234382, 0.235812, 0.236718, 0.236718, 0.236718, 0.236718}};
  tri_sec6 = {.p0 = {1.34048, 1.3457, 1.34985, 1.35331, 1.35303, 1.35303, 1.35303, 1.35303},
              .p1 = {0.229432, 0.232767, 0.234884, 0.236384, 0.237264, 0.237264, 0.237264, 0.237264}};
  const TriangleParameterLUT TRI_SIMULATION_OUTB = {tri_sec1, tri_sec2, tri_sec3, tri_sec4, tri_sec5, tri_sec6};

  // calculate energy bin
  double P = p->par()->getP();
  int energy_bin = P <= 3 ? 0 : P <= 4 ? 1 : P <= 5 ? 2 : P <= 6 ? 3 : P <= 7 ? 4 : P <= 8 ? 5 : P <= 9 ? 6 : 7;

  const TriangleParameterLUT* tri_LUT;
  if (spring2019)
    tri_LUT = &TRI_SPRING2019;
  else if (simulation) {
    tri_LUT = inbending ? &TRI_SIMULATION_INB : &TRI_SIMULATION_OUTB;
  } else {
    tri_LUT = inbending ? &TRI_FALL2018_INB : &TRI_FALL2018_OUTB;
  }
  const TwoParamsAllEnergies& tri_params = tri_LUT->at(p->cal(clas12::PCAL)->getSector() - 1);
  const double& tri_p0 = tri_params.p0[energy_bin];
  const double& tri_p1 = tri_params.p1[energy_bin];

  bool pass_triangle =
      p->cal(clas12::PCAL)->getEnergy() / P > (tri_p1 - tri_p0 * p->cal(clas12::ECIN)->getEnergy() / P);


  //////////////////////////////////////////////////////////////////////////////
  // threshold cut on SF PCAL
  //////////////////////////////////////////////////////////////////////////////
  bool pass_threshold = (p->cal(clas12::PCAL)->getEnergy() / P) > 0.05;

  // final decision
  return pass_band && pass_triangle && pass_threshold;
}

bool EC_hit_position_fiducial_cut_homogeneous(clas12::region_particle* p, int tightness, bool inbending) {
  if (tightness > 1 or tightness < 3)
    throw std::runtime_error("[EC_hit_position_fiducial_cut_homogeneous] tightness must be 1, 2, or 3.");

  // Cut using the natural directions of the scintillator bars/ fibers:
  double v = p->cal(clas12::EC)->getLv();
  double w = p->cal(clas12::EC)->getLw();
  [[maybe_unused]] double u = p->cal(clas12::EC)->getLu();

  /// v + w is going from the side to the back end of the PCAL, u is going from
  /// side to side 1 scintillator bar is 4.5 cm wide. In the outer regions
  /// (back) double bars are used. a cut is only applied on v and w

  struct bounds_vw {
    bounds v{0, 0};
    bounds w{0, 0};
    // double min_w;
    // double max_w;
  };

  using EC_fiducial_bounds_vw_LUT = std::array<std::array<bounds_vw, 3>, 2>;

  const EC_fiducial_bounds_vw_LUT EC_FIDUCIAL_VW_LIMITS_LUT = {{// inbending
                                                         {{
                                                             // Loose
                                                             {{9.0, 400.0}, {9.0, 400.0}},
                                                             // Medium
                                                             {{14.0, 400.0}, {14.0, 400.0}},
                                                             // Tight
                                                             {{19.0, 400.0}, {19.0, 400.0}},
                                                         }},
                                                         // outbending
                                                         {{
                                                             // Loose
                                                             {{9.0, 400.0}, {9.0, 400.0}},
                                                             // Medium
                                                             {{14.0, 400.0}, {14.0, 400.0}},
                                                             // Tight
                                                             {{19.0, 400.0}, {19.0, 400.0}},
                                                         }}}};
  // overrides for sector #4 (idx=3)
  const bounds_vw LOOSE_SECTOR4_OVERRIDES = {13.5, 400.0, 9.0, 400.0};

  // select the right bounds from the LUT
  bounds_vw bounds = EC_FIDUCIAL_VW_LIMITS_LUT[static_cast<int>(inbending)][tightness - 1];

  // Override for sector 4 with loose tightness (TODO: reason unkown)
  if (tightness == 1 && p->cal(clas12::EC)->getSector() == 4) {
    bounds = LOOSE_SECTOR4_OVERRIDES;
  }

  return (v >= bounds.v.lower && v <= bounds.v.upper && w >= bounds.w.lower && w <= bounds.w.upper);
}

bool DC_fiducial_cut_edge(clas12::region_particle* p, int region, bool inbending) {
  // throw not_implemented_error("DC_fiducial_cut_edge not implemented");
  if (region < 1 or region > 3) throw std::runtime_error("[DC_fiducial_cut_edge] `region` must be 1, 2, or 3");

  using DCEdgeCuts = struct {
    std::array<double, 3> inbending;
    std::array<double, 3> outbending;
  };
  const std::map<int, DCEdgeCuts> DCedge_LUT = {
      {11, {{5.0, 5.0, 10.0}, {3.0, 3.0, 10.0}}}, {2122, {{2.5, 2.5, 9.0}, {3.5, 3.0, 7.0}}},
      {211, {{2.5, 2.5, 9.0}, {3.5, 2.5, 6.5}}},  {-211, {{3.5, 3.0, 7.0}, {2.5, 2.5, 10.0}}},
      {321, {{2.5, 2.0, 9.0}, {3.5, 2.5, 6.5}}},  {-321, {{3.5, 2.5, 5.0}, {2.5, 2.5, 10.0}}}};

  if (!PID_cut(p, 11) && !PID_cut(p, 2122) && !PID_cut(p, 211) && !PID_cut(p, -211) && !PID_cut(p, 321) &&
      !PID_cut(p, -321)) {
    throw std::runtime_error("[DC_fiducial_cut_edge] Attempting cut on unkown PID. Cut on PID first.");
  }
  // get iterator to arrays in the DCedge_LUT
  auto it = DCedge_LUT.find(p->getPid());
  if (it == DCedge_LUT.end())
    throw std::runtime_error("[DC_fiducial_cut_edge] AHhhhh, Panic! This should not be reachable.");
  // select edge cut based on inbending or outbending
  double edge_cut = inbending ? it->second.inbending[region - 1] : it->second.outbending[region - 1];

  // Is switching on region correct? And is clas12::DC[1,2,3] correct layers?
  // select_electron calls this cut with region = 1,2,3, does that translate to DC1,DC3, DC6 correctly?
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
      break;
  }

  return edge_val > edge_cut;
}

bool DC_z_vertex_cut(clas12::region_particle* p) {
  // isn't it already implemented in clas12::zVertexFilter ?
  throw not_implemented_error("[DC_z_vertex_cut] Not implemented yet.");

  // // pass 2 (adjusted to cross sections)

  // double vz_min_sect_inb[] = {-8, -8, -8, -8, -8, -8};
  // double vz_max_sect_inb[] = {2, 2, 2, 2, 2, 2};

  // double vz_min_sect_outb[] = {-11, -11, -11, -11, -11, -11};
  // double vz_max_sect_outb[] = {1, 1, 1, 1, 1, 1};

  // double vz_min_sect[6];
  // double vz_max_sect[6];

  // for (Int_t i = 0; i < 6; i++) {
  //   if (inbending == true) {
  //     vz_min_sect[i] = vz_min_sect_inb[i];
  //     vz_max_sect[i] = vz_max_sect_inb[i];
  //   }
  //   if (outbending == true) {
  //     vz_min_sect[i] = vz_min_sect_outb[i];
  //     vz_max_sect[i] = vz_max_sect_outb[i];
  //   }
  // }

  // double vz_min = 0;
  // double vz_max = 0;

  // for (Int_t k = 0; k < 6; k++) {
  //   if (part_Cal_PCAL_sector[j] - 1 == k) {
  //     vz_min = vz_min_sect[k];
  //     vz_max = vz_max_sect[k];
  //   }
  // }

  // if (part_vz[j] > vz_min && part_vz[j] < vz_max)
  //   return true;
  // else
  //   return false;
}

bool FT_eid_FTCAL_fiducial_cut(clas12::region_particle* p) {
  double theta = std::acos(p->par()->getPz() / p->par()->getP()) * 180 / std::numbers::pi;

  return theta > 2.5 && theta < 4.5;
}

bool FT_eid_FTTRK_fiducial_cut(clas12::region_particle* p) {
  // same as FT_eid_FTCAL_fiducial_cut
  return FT_eid_FTCAL_fiducial_cut(p);
}

bool FT_eid_FTHODO_fiducial_cut(clas12::region_particle* p) {
  // same as FT_eid_FTCAL_fiducial_cut
  return FT_eid_FTCAL_fiducial_cut(p);
}

bool FT_eid_energy_vs_radius_cut([[maybe_unused]] clas12::region_particle* _) { return true; }

bool phot_beta_cut(clas12::region_particle* p, int tightness) {
  if (tightness > 1 or tightness < 3) throw std::runtime_error("[phot_beta_cut] tightness must be 1, 2, or 3.");

  std::array<bounds, 3> beta_cut_LUT = {{
      {0.9, 2.0},   // loose
      {0.9, 1.1},   // medium
      {0.95, 1.05}  // tight
  }};

  float beta = p->par()->getBeta();
  return beta_cut_LUT[tightness - 1].lower < beta && beta < beta_cut_LUT[tightness - 1].upper &&
         p->par()->getP() > 0.10;
}

bool phot_EC_sampling_fraction_cut(clas12::region_particle* p) {
  double ECfrac_min = 0;
  double ECfrac_max = 1;

  double pcal = p->cal(clas12::PCAL)->getEnergy();
  double ecin = p->cal(clas12::ECIN)->getEnergy();
  double ecout = p->cal(clas12::ECOUT)->getEnergy();
  double total = pcal + ecin + ecout;

  double cutvalue = total / p->par()->getP();
  return (ECfrac_min < cutvalue && cutvalue < ECfrac_max);
}

bool phot_EC_outer_vs_EC_inner_cut(clas12::region_particle* p) {
  throw not_implemented_error("[phot_EC_outer_vs_EC_inner_cut] Not implemented yet.");
  double edep_min = 0.01;
  return (p->cal(clas12::ECIN)->getEnergy() + p->cal(clas12::ECOUT)->getEnergy()) > edep_min:
}

bool neutr_beta_cut(clas12::region_particle* p, int run) {
  // very similar to phot_beta_cut
  bounds limits = {0., 0.95};
  return limits.lower < p->par()->getBeta() && p->par()->getBeta() < limits.upper;
}

bool basic_FTOF_cut(clas12::region_particle* p) {
  return p->sci(clas12::FTOF1A)->getSector() != 0 || p->sci(clas12::FTOF1B)->getSector() != 0 ||
         p->sci(clas12::FTOF2)->getSector() != 0;  // equal to "is somewhere in FTOF?", I guess...
}

bool FT_photid_FTCAL_fiducial_cut(clas12::region_particle* p) {
  double clusX = p->ft(clas12::FTCAL)->getX();
  double clusY = p->ft(clas12::FTCAL)->getY();
  TVector3 V3ECalPos(clusX, clusY, 0);

  bool res = true && V3ECalPos.Mag() > 8 && V3ECalPos.Mag() < 15 &&
             TMath::Power(clusX + 8.5, 2) + TMath::Power(clusY - 10, 2) > 1.5 * 1.5 &&
             TMath::Power(clusX + 10, 2) + TMath::Power(clusY + 5, 2) > 1.5 * 1.5 &&
             TMath::Power(clusX + 6, 2) + TMath::Power(clusY + 13.5, 2) > 2 * 2 &&
             TMath::Power(clusX - 4, 2) + TMath::Power(clusY + 6.7, 2) > 1.5 * 1.5 &&
             TMath::Power(clusX - 6, 2) + TMath::Power(clusY + 6, 2) > 1;

  return res;
}

bool FT_photid_beta_cut([[maybe_unused]] clas12::region_particle* _) { return true; }

}  // namespace cuts