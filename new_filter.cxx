// System
#include <cmath>
#include <chrono>
#include <fmt/format.h>

// ROOT
#include <Rtypes.h>
#include <TTree.h>
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
// #include <ROOT/RNTupleReader.hxx>
#include <ROOT/RMiniFile.hxx>

// Import classes from experimental namespace for the time being
// using RNTupleModel = ROOT::RNTupleModel;
// using RNTupleReader = ROOT::Experimental::RNTupleReader;
// using RNTupleWriter = ROOT::RNTupleWriter;

// #include <ROOT/>

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "RDataFrame.hxx"
#endif
// #include <TDatabasePDG.h>

// CLAS12
#include <hipo4/dictionary.h>
#include <clas12reader.h>
#include <Iguana.h>
using namespace iguana::particle; // for lightweight and short PDG code enum, alt.: TDatabasePDG::Instance()->GetParticle("name")->PdgCode()

// Own
// #include "OutputColumnFactory.h"

using variant_pointer=std::variant<int *, float *, double *, std::vector<int> *, std::vector<double> *>;
using variant_map = std::unordered_map<std::string, variant_pointer>;

struct VarInfo
{
  std::string name;
  variant_pointer value;
};
// class BookKeeping
// {
// private:
//   /* data */
// public:
//   BookKeeping(/* args */);
//   ~BookKeeping();
//   template <typename T> 
//   addBranch(TTree *, variant_map &, const std::string &, T &);
// };

// Helper: create and register branches dynamically
template <typename T>
void addBranch(TTree *tree, variant_map &vars, const std::string &name, T &var)
{
  tree->Branch(name.c_str(), &var);
  vars[name] = &var;
}

// Helper: set elements inside the variant map
template <typename T>
void setVar(variant_map& vars, const std::string& name, const T& value) {
    if (auto p = std::get_if<T*>(&vars.at(name))) {
        **p = value;
    } else {
        throw std::runtime_error("Type mismatch for variable: " + name);
    }
}

void new_filter()
{
  std::cout << "Called without arguments." << std::endl;
}

int new_filter(std::string inFile, std::string outputfile = "/dev/null", uint numEvents = -1)
{
  clas12::clas12reader c12_reader(inFile.c_str(), {0});

  hipo::dictionary dict = c12_reader.getDictionary();
  int events = numEvents != -1 ? numEvents : c12_reader.getReader().getEntries();
  // dict.show(); // Print all bank names

  clas12::clas12databases db;
  c12_reader.connectDataBases(&db);
  c12_reader.applyQA("pass2");
  // c12_reader.db()->qadb_requireOkForAsymmetry(true); // what is this? Is this needed in general or specific for every ana task?

  // // Prepare Filters
  // clas12root::Iguana ig{};
  // ig.GetTransformers().Use("clas12::MomentumCorrection");
  // ig.GetFilters().Use("clas12::zVertexFilter");
  // ig.GetCreators().Use("physics::InclusiveKinematics");
  // ig.SetOptionAll("log", "debug");
  // // ig.Start();

  // std::vector<clas12::region_part_ptr> electron = c12_reader.getByID(PDG::electron);
  // std::vector<clas12::region_part_ptr> pi_m = c12_reader.getByID(PDG::pi_minus);
  // std::vector<clas12::region_part_ptr> pi_p = c12_reader.getByID(PDG::pi_plus);

  //////////////////////////////////////////////////////////////////////////////
  // // PREPARE OUTPUT
  //////////////////////////////////////////////////////////////////////////////
  // // OPTION 1: TreeFactory
  // // maybe do this depending on config file
  // TreeFactory *tf = new TreeFactory(); //maybe later add more metadata, like dtype (see filter_clas12_hipo.cxx:354)
  // // tf.add_momenta();  // p4_ele_px, ...
  // // tf.add_detected(); // ele_det, ...
  // // tf.add_sector();   // ele_sec, ...
  // // tf.add_chi2pid();  // p4_ele_chi2pid, ...
  // // // tf.add_<others>()  // whatever else is needed
  // // Then in evt loop:
  // // tf->get_buffer()->at("helicity") = c12_reader.event().getHelicity() // set writebuffer
  // // tf->get_tree()->Fill()                                              // fill buffer into tree
  // // Problems: What about different types?
  // //           In evt loop, writebuffer is a static typed `map<str,double>`

  // // // OPTION 2: RNTuple
  // std::unique_ptr<ROOT::RNTupleModel> model = ROOT::RNTupleModel::Create();
  // // create fields
  // auto fldEvent = model->MakeField<int>("eventNumber");
  // auto fldHelicity = model->MakeField<double>("helicity");
  // auto fldBeamCharge = model->MakeField<double>("beam_charge");

  // // tell the RNTupleWriter about the structure/model and define names
  // auto file = ROOT::RNTupleWriter::Recreate(std::move(model), "nTupleName","file.root");

  // // OPTION 3: std::variant based tree filling
  TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
  TTree *tree = new TTree("IamGroot", "Example Event Tree");
  variant_map vars;
  int fldEvent = 0;         // Some place in memory to store values
  int fldHelicity = 0;      // Some place in memory to store values
  float fldBeamCharge = 0.; // Some place in memory to store values
  addBranch(tree, vars, "eventnumber", fldEvent);       // Register variable and memory to tree and map book-keeping
  addBranch(tree, vars, "helicity", fldHelicity);       // Register variable and memory to tree and map book-keeping
  addBranch(tree, vars, "beam_charge", fldBeamCharge);  // Register variable and memory to tree and map book-keeping

  //////////////////////////////////////////////////////////////////////////////
  // Event loop
  int count = 0;
  uint progressInterval = 1;
  auto t0 = std::chrono::steady_clock::now();
  auto last_update = t0;

  fmt::print("Starting event loop for {} events...\n", events);
  while (c12_reader.next() && ((events != -1 && count < events) || (events == -1)))
  {
    // body
    
    setVar(vars, "eventnumber", c12_reader.runconfig()->getEvent());
    setVar(vars, "helicity", c12_reader.event()->getHelicity());
    setVar(vars, "beam_charge", c12_reader.event()->getBeamCharge()); 

    tree->Fill();
    // progress update
    if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - last_update).count() >= progressInterval)
    {
      fmt::print("Processed {:>7d}/{} ({:>3.0f}%)\n", count, events, 100. * count / events);
      last_update = std::chrono::steady_clock::now();
    }
    count++;
  }

  std::cout << std::endl
            << "done." << std::endl;
  tree->Write();
  outfile->Close();

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