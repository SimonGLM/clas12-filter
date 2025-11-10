// System
#include <cmath>
#include <chrono>
#include <fmt/format.h>

// ROOT
#include <Rtypes.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RMiniFile.hxx>

// Import classes from experimental namespace for the time being
using RNTupleModel = ROOT::RNTupleModel;
using RNTupleReader = ROOT::RNTupleReader;
using RNTupleWriter = ROOT::RNTupleWriter;

// CLAS12
#include <hipo4/dictionary.h>
#include <clas12reader.h>
#include <Iguana.h>
namespace pdg = iguana::particle; // for lightweight and short PDG code enum, alt.: TDatabasePDG::Instance()->GetParticle("name")->PdgCode()

// Own stuff

// A helper class that wraps an RNTupleModel and allows building it up dynamically
// by adding fields of various types identified by their names as strings.
// Also provides type-safe setting of the values of the fields via templated SetValue<T>() method.
// DynamicVarStore::MoveModel() must be used to move the ownership of the built model into
// the RNTupleWriter.
class DynamicVarStore
{
public:
  // Constructor that takes ownership of the model
  DynamicVarStore(std::unique_ptr<ROOT::RNTupleModel> model)
      : fModel(std::move(model)) {}

  // Build the model by adding fields of type T with the given name
  template <typename T>
  void AddField(const std::string &name)
  {
    // Create field in the model and store the shared_ptr<T> in the map
    fMap[name] = fModel->MakeField<T>(name);
  }

  template <typename T>
  void SetValue(const std::string &name, const T &value)
  {
    assert (IsField(name) && "Field does not exist in DynamicVarStore");
    // Get the shared_ptr<T> from the map, dereference it and set the value
    *std::get<std::shared_ptr<T>>(fMap.at(name)) = value;
  }

  void ResetVectorFields()
  {
    // Get the shared_ptr<std::vector<T>> from the map and clear the vector
    for (const auto &[name, ptr] : fMap)
    {
      std::visit([&](auto &&arg)
                 {
              using U = std::decay_t<decltype(arg)>;
              if constexpr (std::is_same_v<U, std::shared_ptr<std::vector<int>>> ||
                            std::is_same_v<U, std::shared_ptr<std::vector<float>>> ||
                            std::is_same_v<U, std::shared_ptr<std::vector<double>>>)
              {
                  arg->clear();
              } }, ptr);
    }
  }

  template <typename T>
 void AppendValue(const std::string &name, const T & value)
  {
    assert (IsField(name) && "Field does not exist in DynamicVarStore");
    // Get the shared_ptr<T> from the map, dereference it and return the value
    std::visit([&](auto &&arg)
               {
            using U = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<U, std::shared_ptr<std::vector<T>>>)
            {
                arg->push_back(value);
            } else {
                throw std::runtime_error("AppendValue called on non-vector field: " + name);
            } }, fMap.at(name));
    // std::get<std::shared_ptr<T>>(fMap.at(name)).get().push_back(value);
  }

  // template <typename T>
  // T* GetValue(const std::string &name)
  // {
  //   // Get the shared_ptr<T> from the map, dereference it and return the value
  //   return std::get<std::shared_ptr<T>>(fMap.at(name)).get();
  // }

  // Returns ownership of the model so the caller can move it into the writer
  std::unique_ptr<ROOT::RNTupleModel> &&MoveModel()
  {
    return std::move(fModel);
  }

  void addMomentaFieldsPerParticle(const std::vector<std::string> &particles)
  {
    for (std::string part : particles)
    {
      for (std::string &&var : {"px", "py", "pz"})
        AddField<std::vector<float>>("p4_" + part + "_" + var);
      AddField<std::vector<double>>("p4_"+part+"_E");
      fParticlesOfInterestForMomentumFields.insert(part);
    }
  }
  void addDetectorFieldsPerParticle(const std::vector<std::string> &particles)
  {
    for (std::string part : particles)
    {
      AddField<std::vector<int>>(part + "_det");
      fParticlesOfInterestForDetectorFields.insert(part);
    }
  }
  
  void addSectorFieldsPerParticle(const std::vector<std::string> &particles)
  {
    for (std::string part : particles)
    {
      AddField<std::vector<double>>("p4_"+part + "_sec");
      fParticlesOfInterestForSectorFields.insert(part);
    }
  }
  
  void addChi2PidFieldsPerParticle(const std::vector<std::string> &particles)
  {
    for (std::string part : particles)
    {
      AddField<std::vector<double>>("p4_"+part + "_chi2pid");
      fParticlesOfInterestForChi2PidFields.insert(part);
    }
  }

  std::set<std::string> fParticlesOfInterestForMomentumFields{};
  std::set<std::string> fParticlesOfInterestForDetectorFields{};
  std::set<std::string> fParticlesOfInterestForSectorFields{};
  std::set<std::string> fParticlesOfInterestForChi2PidFields{};

  bool IsField(const std::string &name){
    return fMap.find(name)!=fMap.end();
  }

private:
  // Each entry stores a shared_ptr<T> of one of the supported types
  typedef std::unordered_map<std::string, std::variant<
                                              std::shared_ptr<int>,
                                              std::shared_ptr<float>,
                                              std::shared_ptr<double>,
                                              // more future types can be added here e.g. TLorentzVector or std::vector
                                              std::shared_ptr<std::vector<int>>,
                                              std::shared_ptr<std::vector<float>>,
                                              std::shared_ptr<std::vector<double>>,
                                              std::shared_ptr<long long>>>
      variant_map;
  variant_map fMap;

  std::unique_ptr<ROOT::RNTupleModel> fModel;
};

// No argument overload if used without arguments
void new_filter()
{
  std::cout << "Called without arguments." << std::endl;
}

int new_filter(std::string inFile, std::string outputfile = "/dev/null", uint numEvents = -1)
{
  ROOT::EnableImplicitMT();
  clas12::clas12reader c12_reader(inFile.c_str(), {0});

  hipo::dictionary dict = c12_reader.getDictionary();
  int events = numEvents != -1 ? numEvents : c12_reader.getReader().getEntries();
  // dict.show(); // Print all bank names

  // // QADB
  //////////////////////////////////////////////////////////////////////////////
  // clas12::clas12databases db;
  // c12_reader.connectDataBases(&db);
  // c12_reader.applyQA("pass2");
  // c12_reader.db()->qadb_requireOkForAsymmetry(true); // what is this? Is this needed in general or specific for every ana task?

  // // Prepare Iguana Filters
  //////////////////////////////////////////////////////////////////////////////
  // clas12root::Iguana ig{};
  // ig.GetTransformers().Use("clas12::MomentumCorrection");
  // ig.GetFilters().Use("clas12::zVertexFilter");
  // ig.GetCreators().Use("physics::InclusiveKinematics");
  // ig.SetOptionAll("log", "debug");
  // // ig.Start();

  // // PREPARE OUTPUT
  //////////////////////////////////////////////////////////////////////////////
  std::unique_ptr<ROOT::RNTupleModel> model = ROOT::RNTupleModel::Create();
  DynamicVarStore vars(std::move(model));

  vars.AddField<int>("eventNumber");
  vars.AddField<int>("helicity");
  vars.AddField<float>("beam_charge");
  // The following lists determine the particles of interest for each group of fields
  vars.addMomentaFieldsPerParticle({ "ele", "prot", "pip", "pim", "Kp", "Km", "neutr"});
  vars.addDetectorFieldsPerParticle({"ele", "prot", "pip", "pim", "Kp", "Km", "neutr", "phot"});
  vars.addSectorFieldsPerParticle({  "ele", "prot", "pip", "pim", "Kp", "Km", "neutr", "phot"});
  vars.addChi2PidFieldsPerParticle({ "ele", "prot", "pip", "pim", "Kp", "Km", "neutr"});
  std::map<int, std::string> pdg_to_name{{11, "ele"},
                                         {2212, "prot"},
                                         {211, "pip"},
                                         {-211, "pim"},
                                         {321, "Kp"},
                                         {-321, "Km"},
                                         {22, "gamma"},
                                         {2112, "neutr"}};

  // Create Writer that takes ownership of the model
  std::string tempfile = "/tmp/rntuple.root";
  auto file = ROOT::RNTupleWriter::Recreate(vars.MoveModel(), "nTupleName", tempfile);

  // Event loop
  //////////////////////////////////////////////////////////////////////////////
  int count = 0;
  uint progressInterval = 1;
  auto t0 = std::chrono::steady_clock::now();
  auto last_update = t0;

  fmt::print("Starting event loop for {} events...\n", events);
  while (c12_reader.next() && ((events != -1 && count < events) || (events == -1))) {
    // access variant stored in map with std::get<> and dereference the shared_ptr to set the that it holds.
    vars.SetValue("eventNumber", c12_reader.runconfig()->getEvent());
    vars.SetValue("helicity", c12_reader.event()->getHelicity());
    vars.SetValue("beam_charge", c12_reader.event()->getBeamCharge());

    // Filter & Cuts here

    // Correct here

    vars.ResetVectorFields();                       // clear and reserve all vector fields
    auto particles = c12_reader.getDetParticles();  // All detected particles in the event
    for (auto &&p : particles){

      // Get PDG enum from namespace pdg=iguana::particle
      pdg::PDG particle_enum = pdg::PDG(p->getPid());
      if (particle_enum){
        std::cout << "Unkown particle with PDG Code '" << p->getPid() << "' is not in particle set."<<std::endl;
        continue;
      }

      // ====================== MOMENTUM FIELDS ======================
      // now we know particle_enum is a known particle, we can get away without std::optionals
      std::string name = pdg::name.at(particle_enum);
      if (vars.fParticlesOfInterestForMomentumFields.contains(name)){
        double m0 = pdg::mass.at(particle_enum); //in GeV
        vars.AppendValue("p4_" + name + "_px",p->par()->getPx());
        vars.AppendValue("p4_" + name + "_py",p->par()->getPy());
        vars.AppendValue("p4_" + name + "_pz",p->par()->getPz());
        vars.AppendValue("p4_" + name + "_E",std::sqrt(std::pow(p->par()->getPx(), 2) +
                                                       std::pow(p->par()->getPy(), 2) +
                                                       std::pow(p->par()->getPz(), 2) +
                                                       std::pow(m0, 2)));
      }
      // ====================== DETECTOR FIELDS ======================
      if (vars.fParticlesOfInterestForDetectorFields.contains(name)){
        // Determine which ... detector part the particle was detected in?
        // CD or FD?
        // Old way used  1000<=abs(part_status)<2000 => ele_det=1
        //               2000<=abs(part_Status)<4000 => ele_det=2
        // no idea if this is correct
        int val = p->getStatus()>=1000 && p->getStatus()<2000 ? 1 :
                  p->getStatus()>=2000 && p->getStatus()<4000 ? 2 : 0;
        vars.AppendValue(name + "_det", val);
      }

      // ====================== SECTOR FIELDS ======================
      if (vars.fParticlesOfInterestForSectorFields.contains(name)){
        vars.AppendValue(name + "_sec", p->getSector());
      }

      // ====================== CHI2PID FIELDS ======================
      if (vars.fParticlesOfInterestForChi2PidFields.contains(name)){
        vars.AppendValue("p4_"+name + "_chi2pid", p->getChi2Pid());
      }
    }

    // tell the RNTupleWriter that the values of the shared_ptrs in
    // the fields have changed
    file->Fill();

    // progress update
    auto ti = std::chrono::steady_clock::now();
    if (std::chrono::duration_cast<std::chrono::seconds>(ti - last_update).count() >= progressInterval)
    {
      float time_remaining = std::chrono::duration_cast<std::chrono::milliseconds>(ti - t0).count() * (events * 1. / count - 1) / 1000.;
      float percent = 100. * count / events;
      std::string progress_bar(int(percent*30)/100, '=');
      progress_bar.append(">");
      progress_bar.resize(30, ' ');
      fmt::print("{:>9d}/{:d} ({:>3.0f}%) [{:s}] {:.1f}s remaining\n", count, events, percent, progress_bar, time_remaining);
      last_update = std::chrono::steady_clock::now();
    }
    count++;
  }
  fmt::print("Processed a total of {} events in {:.1f} seconds.\n", count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t0).count() / 1000.);

  file.reset(); // close file, ~RNTupleWriter() writes to disk on destruction

  // RNTupleWriter does not write the Tree structure to disk.
  // Work around this using tempfile and RDataFrame.Snapshot, which does.
  ROOT::RDF::RNode df = ROOT::RDF::FromRNTuple("nTupleName", tempfile.c_str());
  ROOT::RDF::RSnapshotOptions opts;
  opts.fOverwriteIfExists = true;
  opts.fVector2RVec = false;
  opts.fOutputFormat = ROOT::RDF::ESnapshotOutputFormat::kTTree;
  df.Snapshot("out_tree", outputfile, df.GetColumnNames(), opts);
  std::remove(tempfile.c_str()); // Delete the temporary file

  std::cout << std::endl
            << "done." << std::endl;

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

  return 0;
}
