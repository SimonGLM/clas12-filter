/*
        new_filter.cxx

        A simple filter program that reads HIPO4 files using clas12reader,
        applies some filters and corrections using Iguana, and writes out
        selected variables into a ROOT RNTuple.

        Author: Simon Glennemeier-Marke (Justus-Liebig-University Giessen, 2025)

        Note:
        Code formatted using clang-format with {BasedOnStyle: Google, ColumnLimit: 120}
*/

// System
#include <fmt/format.h>

#include <chrono>
#include <cmath>
#include <string>

// ROOT
#include <Math/LorentzVector.h>
#include <Rtypes.h>
#include <TDatabasePDG.h>
#include <TTree.h>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RMiniFile.hxx>
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RSnapshotOptions.hxx>

// CLAS12
#include <clas12defs.h>
#include <clas12reader.h>
#include <region_particle.h>
using region_part_ptr = clas12::region_particle*;  // needed for compilation

#include <Iguana.h>                                // needs to be included after clas12reader.h
#include <iguana/algorithms/AlgorithmSequence.h>
#include <iguana/algorithms/clas12/EventBuilderFilter/Algorithm.h>

// Own
// #include "cuts.h" // future

using namespace std::chrono;
using FourVector = ROOT::Math::PxPyPzMVector;

std::string pdg_name(const int& pdg) {
  return pdg == 11     ? "ele"
         : pdg == 22   ? "phot"
         : pdg == 2212 ? "prot"
         : pdg == 2112 ? "neutr"
         : pdg == 211  ? "pip"
         : pdg == -211 ? "pim"
         : pdg == 321  ? "Kp"
         : pdg == -321 ? "Km"
                       : "unknown";
};

// A helper class that wraps an RNTupleModel and allows building it up dynamically
// by adding fields of various types identified by their names as strings.
// Also provides type-safe setting of the values of the fields via templated SetValue<T>() method.
// DynamicVarStore::MoveModel() must be used to move the ownership of the built model into
// the RNTupleWriter.
class DynamicVarStore {
 public:
  // Constructor that takes ownership of the model
  DynamicVarStore(std::unique_ptr<ROOT::RNTupleModel> model) : fModel(std::move(model)) {}

  // Build the model by adding fields of type T with the given name
  template <typename T>
  void AddField(const std::string& name) {
    // Create field in the model and store the shared_ptr<T> in the map
    fMap[name] = fModel->MakeField<T>(name);
  }

  template <typename T>
  void SetValue(const std::string& name, const T& value) {
    if (!IsField(name)) throw std::invalid_argument(std::format("Field '{}' does not exist in DynamicVarStore", name));
    // Get the shared_ptr<T> from the map, dereference it and set the value
    // no `std::get<...>(...).reset(value)` here as that would change the address that is pointed to
    *std::get<std::shared_ptr<T>>(fMap.at(name)) = value;
  }

  void ResetVectorFields() {
    // Get the shared_ptr<std::vector<T>> from the map and clear the vector
    for (const auto& [name, ptr] : fMap) {
      std::visit(
          [&](auto&& arg) {
            using U = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<U, std::shared_ptr<std::vector<int>>> ||
                          std::is_same_v<U, std::shared_ptr<std::vector<double>>>) {
              arg->clear();
            }
          },
          ptr);
    }
  }

  template <typename T>
  void AppendValue(const std::string& name, const T& value) {
    if (!IsField(name)) throw std::invalid_argument(std::format("Field '{}' does not exist in DynamicVarStore", name));
    // Get the shared_ptr<T> from the map, check type, then push_back
    std::visit(
        [&](auto&& arg) {
          using U = std::decay_t<decltype(arg)>;
          if constexpr (!std::is_same_v<U, std::shared_ptr<std::vector<T>>>) {
#if __GNUG__
            throw std::invalid_argument(std::format("Type missmatch between field '{}' and value of type '{}'", name,
                                                    abi::__cxa_demangle(typeid(arg).name(), NULL, NULL, nullptr)));
#else
            throw std::invalid_argument(std::format("Type missmatch between field '{}' and value type", name));
#endif
          } else {
            arg->push_back(value);
          }
        },
        fMap.at(name));
  }

  // Returns the reference to the Model pointer, so that ownership can be moved
  std::unique_ptr<ROOT::RNTupleModel>& Model() { return fModel; }
  // std::unique_ptr<ROOT::RNTupleModel>&& MoveModel() { return std::move(fModel); }

  void AddMomentaFieldsPerParticle(const std::vector<int>& particles) {
    for (int part : particles) {
      for (std::string&& var : {"px", "py", "pz", "E"}) {
        AddField<std::vector<double>>("p4_" + pdg_name(part) + "_" + var);
      }
      fParticlesOfInterestForMomentumFields.insert(part);
    }
  }

  void AddDetectorFieldsPerParticle(const std::vector<int>& particles) {
    for (int part : particles) {
      AddField<std::vector<int>>(pdg_name(part) + "_det");
      fParticlesOfInterestForDetectorFields.insert(part);
    }
  }

  void AddSectorFieldsPerParticle(const std::vector<int>& particles) {
    for (int part : particles) {
      AddField<std::vector<int>>(pdg_name(part) + "_sec");
      fParticlesOfInterestForSectorFields.insert(part);
    }
  }

  void AddChi2PidFieldsPerParticle(const std::vector<int>& particles) {
    for (int part : particles) {
      AddField<std::vector<double>>("p4_" + pdg_name(part) + "_chi2pid");
      fParticlesOfInterestForChi2PidFields.insert(part);
    }
  }

  void AddDCFields() {
    for (const char& chord : {'x', 'y', 'z'}) AddField<std::vector<double>>(std::format("p4_prot_dc{}1", chord));
    for (const int& i : {1, 2, 3}) AddField<std::vector<double>>(std::format("p4_ele_dcedge{}", i));
    fParticlesOfInterestForDCFields.insert(2212);
    fParticlesOfInterestForDCFields.insert(11);
  }

  std::set<int> fParticlesOfInterestForMomentumFields{};
  std::set<int> fParticlesOfInterestForDetectorFields{};
  std::set<int> fParticlesOfInterestForSectorFields{};
  std::set<int> fParticlesOfInterestForChi2PidFields{};
  std::set<int> fParticlesOfInterestForDCFields{};

  bool IsField(const std::string& name) { return fMap.contains(name); }

 private:
  // Each entry stores a shared_ptr<T> of one of the supported types
  using variant_type = std::variant<std::shared_ptr<int>, std::shared_ptr<float>, std::shared_ptr<double>,
                                    std::shared_ptr<std::vector<int>>, std::shared_ptr<std::vector<float>>,
                                    std::shared_ptr<std::vector<double>>>;
  using variant_map = std::unordered_map<std::string, variant_type>;

  variant_map fMap;
  std::unique_ptr<ROOT::RNTupleModel> fModel;
};

// No argument overload if used without arguments
void new_filter() { std::cout << "Called without arguments." << std::endl; }

void new_filter(std::string inFile, std::string outputfile = "/dev/null", uint numEvents = 0) {
  ROOT::EnableImplicitMT();
  auto c12 = std::make_unique<clas12::clas12reader>(inFile);
  c12->setVerbose();
  int events = numEvents != 0 ? numEvents : c12->getReader().getEntries();

  // hipo::dictionary dict = c12->getDictionary();
  // dict.show(); // Print all bank names

  // QADB
  //////////////////////////////////////////////////////////////////////////////
  clas12::clas12databases db;
  c12->connectDataBases(&db);
  c12->applyQA("pass2");
  // c12->db()->qadb_requireOkForAsymmetry(true);  // what is this?
  // Is this needed in general or specific for every ana task?

  // Clas12 Filters
  // c12->AddExactPid(11, 1);
  // c12->AddAtLeastPid(11, 1);
  // c12->AddZeroOfRestPid();


  // // Prepare Clas12 Iguana Filters
  //////////////////////////////////////////////////////////////////////////////
  clas12root::Iguana ig{};
  ig.SetClas12(&c12);
  // ------ Filters ------
  ig.GetFilters().Use("clas12::zVertexFilter");
  // ig.GetFilters().Use("clas12::FiducialFilter");   // Purged from source tree
  // ig.GetFilters().Use("clas12::PhotonGBTFilter");  // Purged from source tree
  // ------ Corrections ------
  // ig.GetTransformers().Use("clas12::MomentumCorrection");  // renamed in Iguana but not in clas12root::Iguana
  // ig.GetTransformers().Use("clas12::FTEnergyCorrection");  // therefore wholy inoperable
  // ------ Creators ------
  // ig.GetCreators().Use("physics::InclusiveKinematics");
  ig.SetOptionAll("log", "debug");
  ig.Start();

  // // Prepare Iguana Algorithm 
  //////////////////////////////////////////////////////////////////////////////
  iguana::AlgorithmSequence seq;
  seq.Add<iguana::clas12::ZVertexFilter>("clas12::ZVertexFilter"); // Filter for PIDs from EventBuilder
  seq.SetOption("pid_filter", "log", "info");
  seq.Start();

  // Preparations to have vec4's ready for iguana transforms
  auto pdg_db = TDatabasePDG::Instance();
  FourVector p4_el(0, 0, 0, pdg_db->GetParticle(11)->Mass());
  FourVector p4_prot(0, 0, 0, pdg_db->GetParticle(2212)->Mass());
  FourVector p4_neutr(0, 0, 0, pdg_db->GetParticle(2112)->Mass());
  FourVector p4_pip(0, 0, 0, pdg_db->GetParticle(211)->Mass());
  FourVector p4_pim(0, 0, 0, pdg_db->GetParticle(-211)->Mass());
  FourVector p4_Kp(0, 0, 0, pdg_db->GetParticle(321)->Mass());
  FourVector p4_Km(0, 0, 0, pdg_db->GetParticle(-321)->Mass());

  // // PREPARE OUTPUT
  //////////////////////////////////////////////////////////////////////////////
  std::unique_ptr<ROOT::RNTupleModel> model = ROOT::RNTupleModel::Create();
  DynamicVarStore vars(std::move(model));

  vars.AddField<int>("eventnumber");
  vars.AddField<int>("helicity");
  vars.AddField<float>("beam_charge");

  // The following lists determine the particles of interest for each group of fields
  // later they will come from a config file
  vars.AddMomentaFieldsPerParticle({11, 2212, 2112, 211, -211, 321, -321});
  vars.AddDetectorFieldsPerParticle({11, 22, 2212, 2112, 211, -211, 321, -321});
  vars.AddChi2PidFieldsPerParticle({11, 2212, 2112, 211, -211, 321, -321});
  vars.AddSectorFieldsPerParticle({11, 22, 2212, 2112, 211, -211, 321, -321});
  vars.AddDCFields();

  // Create Writer that takes ownership of the model
  std::string tempfile = "/tmp/rntuple.root";
  auto file = ROOT::RNTupleWriter::Recreate(std::move(vars.Model()), "nTupleName", tempfile);

  // Event loop
  //////////////////////////////////////////////////////////////////////////////
  int count = 0;
  uint progressInterval = 1;
  auto t0 = steady_clock::now();
  auto last_update = t0;

  fmt::print("Starting event loop for {} events...\n", events);
  while (c12->next() && ((events != -1 && count < events) || (events == -1)))  // 15.4s in Loop (c12.next() 8.8s)
  {
    // access variant stored in map with std::get<> and dereference the shared_ptr to set the that it holds.
    vars.SetValue("eventnumber", c12->runconfig()->getEvent());
    vars.SetValue("helicity", c12->event()->getHelicity());
    vars.SetValue("beam_charge", c12->event()->getBeamCharge());
    
    // Clear Vector Fields
    vars.ResetVectorFields();                 // clear and reserve all vector fields

    // ====================== EVENT CUTS ======================
    // do cuts on event level here if needed
    // something like EventBuilderFilter for example
    
    // maybe Iguana?
    // Filter here
    // ig.GetFilters().doAllFilters();

    // Correct here
    // ig.GetTransformers().doAllCorrections();
    
    auto particles = c12->getDetParticles();  // All detected particles in the event
    for (auto&& p : particles)                // 3.4s in loop
    {
      std::string name = pdg_name(p->getPid());
      if (name == "unknown") {
        // std::cout << "Unknown particle with PDG code " << p->getPid() << " found, skipping..." << std::endl;
        continue;
      }
      
      // Todo: helicity==inbending/outbending?
      // c12->runconfig()->getSolenoid(); // is this inbending/outbending?
      // c12->runconfig()->getTorus();    // or this?
      // c12->event()->getHelicity();     // or this?
      
      // ====================== PARTICLE CUTS ======================
      // if (!cuts::CC_nphe_cut(p)) continue; // future
      

      // ===========================================================
      //                        Collect data
      // ===========================================================
      TParticlePDG* pdg = pdg_db->GetParticle(p->getPid());
      if (pdg == nullptr) continue;

      // ====================== MOMENTUM FIELDS ======================
      // now we know particle_enum is a known particle, we can get away without std::optionals
      if (vars.fParticlesOfInterestForMomentumFields.contains(p->getPid())) {
        vars.AppendValue("p4_" + name + "_px", static_cast<double>(p->par()->getPx()));
        vars.AppendValue("p4_" + name + "_py", static_cast<double>(p->par()->getPy()));
        vars.AppendValue("p4_" + name + "_pz", static_cast<double>(p->par()->getPz()));
        vars.AppendValue("p4_" + name + "_E",
                         static_cast<double>(std::sqrt(std::pow(p->par()->getPx(), 2) + std::pow(p->par()->getPy(), 2) +
                                                       std::pow(p->par()->getPz(), 2) + std::pow(pdg->Mass(), 2))));
      }
      // ====================== DETECTOR FIELDS ======================
      if (vars.fParticlesOfInterestForDetectorFields.contains(p->getPid())) {
        // Old way used  1000<=abs(part_status)<2000 => FT: ele_det=1
        //               2000<=abs(part_Status)<4000 => FD: ele_det=2
        //               4000<=abs(part_Status)      => CD: ele_det=3
        // no idea if this is correct
        int status = std::abs(p->getStatus());
        int val = status >= 1000 && status < 2000 ? 1 : status >= 2000 && status < 4000 ? 2 : status >= 4000 ? 3 : -1;
        vars.AppendValue(name + "_det", val);
      }

      // ====================== SECTOR FIELDS ======================
      if (vars.fParticlesOfInterestForSectorFields.contains(p->getPid())) {
        vars.AppendValue(name + "_sec", static_cast<int>(p->getSector()));
      }

      // ====================== CHI2PID FIELDS ======================
      if (vars.fParticlesOfInterestForChi2PidFields.contains(p->getPid())) {
        vars.AppendValue("p4_" + name + "_chi2pid", static_cast<double>(p->getChi2Pid()));
      }

      // ====================== DC FIELDS ======================
      if (vars.fParticlesOfInterestForDCFields.contains(p->getPid())) {
        if (name == pdg_name(11)) {  // alternatively: if (p->getPid()==11)
          vars.AppendValue("p4_" + name + "_dcedge1", static_cast<double>(p->traj(clas12::DC, 6)->getEdge()));
          vars.AppendValue("p4_" + name + "_dcedge2", static_cast<double>(p->traj(clas12::DC, 18)->getEdge()));
          vars.AppendValue("p4_" + name + "_dcedge3", static_cast<double>(p->traj(clas12::DC, 36)->getEdge()));
        } else if (name == pdg_name(2122)) {
          vars.AppendValue("p4_" + name + "_dcx1", static_cast<double>(p->traj(clas12::DC, 6)->getX()));
          vars.AppendValue("p4_" + name + "_dcy1", static_cast<double>(p->traj(clas12::DC, 6)->getY()));
          vars.AppendValue("p4_" + name + "_dcz1", static_cast<double>(p->traj(clas12::DC, 6)->getZ()));
        }
      }
    }

    // tell the RNTupleWriter that the values of the shared_ptrs in
    // the fields have changed
    file->Fill();  // 2.7 s

    // progress update
    auto ti = steady_clock::now();
    if (duration_cast<seconds>(ti - last_update).count() >= progressInterval) {
      float time_remaining = duration_cast<milliseconds>(ti - t0).count() * (events * 1. / count - 1) / 1000.;
      float percent = 100. * count / events;
      std::string progress_bar(int(percent * 30) / 100, '=');
      progress_bar.append(">");
      progress_bar.resize(30, ' ');
      fmt::print("{:>9d}/{:d} ({:>3.0f}%) [{:s}] {:.1f}s remaining\n", count, events, percent, progress_bar,
                 time_remaining);
      last_update = steady_clock::now();
    }
    count++;
  }
  fmt::print("Processed a total of {} events in {:.1f} seconds.\n", count,
             duration_cast<milliseconds>(steady_clock::now() - t0).count() / 1000.);

  auto t1 = steady_clock::now();
  std::cout << "Writing to RNTuple to temporary file..." << std::flush;
  file.reset();  // close file, ~RNTupleWriter() writes to disk on destruction
  std::cout << std::format(" {:.1f}s", duration_cast<milliseconds>(steady_clock::now() - t1).count() / 1000.)
            << std::endl;
  // RNTupleWriter does not write the Tree structure to disk.
  // Work around this using tempfile and RDataFrame.Snapshot, which does.

  auto t2 = steady_clock::now();
  std::cout << "Reading temporary file into RDataFrame..." << std::flush;
  ROOT::RDF::RNode df = ROOT::RDF::FromRNTuple("nTupleName", tempfile.c_str());
  std::cout << std::format(" {:.1f}s", duration_cast<milliseconds>(steady_clock::now() - t2).count() / 1000.)
            << std::endl;

  auto t3 = steady_clock::now();
  std::cout << std::format("Writing RDataFrame into '{}'...", outputfile) << std::flush;
  ROOT::RDF::RSnapshotOptions opts;
  opts.fOverwriteIfExists = true;
  opts.fVector2RVec = false;
  opts.fOutputFormat = ROOT::RDF::ESnapshotOutputFormat::kTTree;
  df.Snapshot("out_tree", outputfile, df.GetColumnNames(), opts);  // 6.5s
  std::remove(tempfile.c_str());                                   // Delete the temporary file
  std::cout << std::format(" {:.1f}s", duration_cast<milliseconds>(steady_clock::now() - t3).count() / 1000.)
            << std::endl;

  std::cout << "done." << std::endl;

  //////////////////////////////////////////////////////////////////////////////
  // // Stale code
  // // Keep for now as reference

  // iguana::AlgorithmSequence seq;
  // seq.Add<iguana::clas12::EventBuilderFilter>("pid_filter"); // Filter for PIDs from EventBuilder
  // seq.SetOption("pid_filter", "log", "info");

  // hipo::bank particles(dict.getSchema("REC::Particle"));
  // for(int i =0; i<particles.getSize();i++){
  //   float px = particles.getFloat("px",i);
  //   float py = particles.getFloat("py",i);
  //   float pz = particles.getFloat("pz",i);
  //   float E = particles.getFloat("E",i);
  //   fmt::print("{:<20} {:<20} {:<20} {:<20}\n",px,py,pz,E);
  // }

  // seq.Start(banks);

  // auto b_config   = hipo::getBanklistIndex(banks, "RUN::config");
  // auto b_particle = hipo::getBanklistIndex(banks, "REC::Particle");

  // int iEvent = 0;
  // // hipo::event evt;
  // // while (reader.next()){
  // while(reader.next(banks) && (numEvents == 0 || iEvent++ < numEvents)) {

  //   auto& bank_config   = banks.at(b_config);
  //   auto& bank_particle = banks.at(b_particle);

  //   // print the event number
  //   fmt::print("===== EVENT {} =====\n", bank_config.getInt("event", 0));

  //   // print the particle bank before Iguana algorithms
  //   fmt::print("----- BEFORE IGUANA -----\n");
  //   bank_particle.show(); // the original particle bank

  //   // run the sequence of Iguana algorithms
  //   seq.Run(banks);

  //   // print the banks after Iguana algorithms
  //   fmt::print("----- AFTER IGUANA -----\n");
  //   bank_particle.show(); // the filtered particle bank, with corrected momenta

  //   // print a table; first the header
  //   fmt::print("----- Analysis Particles -----\n");
  //   fmt::print("  {:<20} {:<20} {:<20} {:<20}\n", "row == pindex", "PDG", "|p|", "sector");
  //   // then print a row for each particle
  //   // - use the `hipo::bank::getRowList()` method to loop over the bank rows that PASS the filter
  //   // - if you'd rather loop over ALL bank rows, iterate from `i=0` up to `i < hipo::bank::getRows()` instead
  //   for(auto const& row : bank_particle.getRowList()) {
  //     auto p = std::hypot(
  //         bank_particle.getFloat("px", row),
  //         bank_particle.getFloat("py", row),
  //         bank_particle.getFloat("pz", row));
  //     auto pdg = bank_particle.getInt("pid", row);
  //     fmt::print("  {:<20} {:<20} {:<20.3f}\n", row, pdg, p);
  //   }
  //   fmt::print("\n");
  // }
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: new_filter <input_file> [output_file] [num_events]" << std::endl;
    return 1;
  }

  std::string input_file = argv[1];
  std::string output_file = (argc >= 3) ? argv[2] : "output.root";
  uint num_events = (argc >= 4) ? std::stoi(argv[3]) : 0;

  new_filter(input_file, output_file, num_events);
  return 0;
}