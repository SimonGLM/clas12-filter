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
// legacy libformat consider std::format in future
#include <fmt/format.h>

#include <chrono>
#include <cmath>
#include <numeric>
#include <string>

// ROOT
#include <Math/LorentzVector.h>
#include <Math/Vector4D.h>
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

// Own
#include "cuts.h"
#include "dynamicvarstore.h"
#include "helpers.h"
#include "selector.h"
#include "selector_context.h"
#include "statistics_collector.h"

// Include cuts implementation for ACLiC
#include "cuts.cpp"

using namespace std::chrono;
using FourVector = ROOT::Math::PxPyPzMVector;

//////////////////////////////////////////////////////////////////////////////
// No argument overload if used without arguments
void new_filter() {
  std::cout << "Called without arguments. " << std::endl << std::endl;
  std::cout << "Usage: clas12root [options] 'new_filter.cxx+(<input_file>, <output_file>, <run_number>, <inbending>, "
               "<num_events>)'"
            << std::endl;
}
//////////////////////////////////////////////////////////////////////////////

void new_filter(std::string inFile, std::string outputfile = "/dev/null", bool isInbending = true, uint numEvents = 0) {
  bool verbose = false;

  ROOT::EnableImplicitMT();
  auto c12 = std::make_unique<clas12::clas12reader>(inFile);
  c12->setVerbose();
  int events = numEvents != 0 ? numEvents : c12->getReader().getEntries();

  // General conditions
  if (verbose) std::cout << "[C12] Setting up general run conditions..." << std::endl;
  // for some reasons the readQuickRunConfig reads the config and only gets the runnum
  int runnum = clas12::clas12reader::readQuickRunConfig(inFile);
  // but we also need the torus polarity, but c12->runconfig()->getTorus() segfaults
  bool inbending = isInbending;

  // arbitrary sequence to process particles in
  // this eases comparison of print statistics and ensures that we always
  // process electrons first (important for reference vertex)
  const std::vector<int> particle_sequence = {11, 2212, 2112, 211, -211, 321, -321, 22};

  // default tightness, overridden in selectors for some cuts <= this may change
  cuts::tightness tightness = cuts::tightness::loose;

  // specify evaluation mode for selectors
  // EarlyReturn: Immediatly return on rejected cut,
  // CompleteTrace: check ALL cuts, then return
  selectors::evaluation_mode = EvaluationMode::CompleteTrace;

  // flags to in/exclude detectors during selectors
  selectors::detector_flags[clas12::FT] = false;  // ignored for photons
  selectors::detector_flags[clas12::FD] = true;
  selectors::detector_flags[clas12::CD] = true;

  //////////////////////////////////////////////////////////////////////////////
  // CCDB + RCDB + QADB
  //////////////////////////////////////////////////////////////////////////////
  clas12::clas12databases db;
  c12->connectDataBases(&db);
  char* CCDB_CONNECTION = getenv("CCDB_CONNECTION");
  char* RCDB_CONNECTION = getenv("RCDB_CONNECTION");
  if (CCDB_CONNECTION == NULL) {
    std::cout << "[CCDB] Environment variable CCDB_CONNECTION not set! Trying remote connection..." << std::endl;
    c12->db()->SetCCDBRemoteConnection();
  } else {
    c12->db()->SetCCDBLocalConnection(CCDB_CONNECTION);
  }
  if (RCDB_CONNECTION == NULL) {
    std::cout << "[RCDB] Environment variable RCDB_CONNECTION not set! Trying remote connection..." << std::endl;
    c12->db()->SetRCDBRemoteConnection();
  } else {
    c12->db()->SetRCDBLocalConnection(RCDB_CONNECTION);
  }
  std::cout << "[QADB] Setting up QA requirements..." << std::endl;
  // c12->applyQA("latest");
  // c12->db()->qadb_requireOkForAsymmetry(true);  // Future: from config file // Deprecated
  // c12->db()->qadb_requireGolden(true);         // Future: from config file

  //////////////////////////////////////////////////////////////////////////////
  // PREPARE OUTPUT
  //////////////////////////////////////////////////////////////////////////////
  std::cout << "[DynamicVarStore] Initializing dynamic data structure..." << std::endl;
  std::unique_ptr<ROOT::RNTupleModel> model = ROOT::RNTupleModel::Create();
  DynamicVarStore vars(std::move(model));

  vars.AddField<int>("eventnumber");
  vars.AddField<int>("helicity");
  vars.AddField<float>("beam_charge");

  // The following lists in the arguments determine the particles of interest for each group of fields
  // later they will come from a config file
  vars.AddMomentaFieldsPerParticle({11, 2212, 2112, 211, -211, 321, -321});
  vars.AddDetectorFieldsPerParticle({11, 22, 2212, 2112, 211, -211, 321, -321});
  vars.AddChi2PidFieldsPerParticle({11, 2212, 2112, 211, -211, 321, -321});
  vars.AddSectorFieldsPerParticle({11, 22, 2212, 2112, 211, -211, 321, -321});
  vars.AddDCFields();

  // Create Writer that takes ownership of the model
  std::string tempfile = "/tmp/clas12_filter_intermediate_rntuple.root";
  if (verbose) std::cout << "[RNTuple] Creating RNTuple writer from model..." << std::endl;
  auto file = ROOT::RNTupleWriter::Recreate(std::move(vars.Model()), "nTupleName", tempfile);

  //////////////////////////////////////////////////////////////////////////////
  // Event loop
  //////////////////////////////////////////////////////////////////////////////
  int count = 0;
  size_t particle_count = 0;
  uint progressInterval = 1;
  auto t0 = steady_clock::now();
  uint last_count = 0;
  auto last_update = t0;
  auto pdg_db = TDatabasePDG::Instance();
  std::cout << std::format("Starting event loop for {} events...", events) << std::endl;
  while (c12->next() && ((events != -1 && count < events) || (events == -1))) {
    // progress bar and stats
    auto ti = steady_clock::now();
    if (duration_cast<seconds>(ti - last_update).count() >= progressInterval) {
      float time_remaining = duration_cast<milliseconds>(ti - t0).count() * (events * 1. / count - 1) / 1000.;
      float percent = 100. * count / events;
      std::string progress_bar(int(percent * 30) / 100, '=');
      progress_bar.append(">");
      progress_bar.resize(30, ' ');
      std::cout << std::format("{:>9d}/{:d} ({:>3.0f}%) [{:s}] {:.1f}s remaining ({:>6.0f} evts/s)\n", count, events,
                               percent, progress_bar, time_remaining, (count - last_count) * 1. / progressInterval)
                << std::flush;
      last_count = count;
      last_update = steady_clock::now();
    }
    count++;

    // Event properties.
    vars.SetValue("eventnumber", c12->runconfig()->getEvent());
    vars.SetValue("helicity", c12->event()->getHelicity());
    vars.SetValue("beam_charge", c12->event()->getBeamCharge());

    // Clear Vector Fields in wrapped RNTuple
    vars.ResetVectorFields();  // clear and reserve all vector fields
    double reference_vertex = std::numeric_limits<double>::quiet_NaN();

    // ====================== EVENT CUTS ======================
    // do cuts on event level here if needed
    // maybe Iguana in the future?
    // Filter: ig.GetFilters().doAllFilters();
    // Corrections: ig.GetTransformers().doAllCorrections();
    // ========================================================

    // ======================= WORK FLOW ========================
    //  1. Store particles in map to process in required order later
    //     1.1 If no electrons, skip event
    //  2. Loop over this map in particular particle order (e-, p, n, pi+, pi-, K+, K-, photons)
    //     2.1 Break loop and reject event if no valid electrons remain
    //     2.2 Sort particles of this type by momentum
    //     2.3 If reference vertex not set, set it to highest momentum electron
    //         This is only done once we are processing something other than electrons and since we have
    //         processed e- first, we can just take the momentum of the first (accepted) electron
    //  3. Loop over particles of this type and apply particle level cuts
    //     3.1 If particle fails cuts, remove it from map and skip it
    //     3.2 Otherwise, fill variables for this particle into RNTuple
    //  4. Fill RNTuple after all particles have been processed
    // ==========================================================

    // =========================== 1. ===========================
    if (verbose) std::cout << std::string(80, '=') << std::endl;
    if (verbose)
      std::cout << std::format("Event #{:<{}} has {} particles:", count, int(std::floor(std::log10(events))) + 1,
                               c12->getDetParticles().size())
                << std::endl;
    std::unordered_map<int, std::vector<clas12::region_part_ptr>> particlesByPDG;
    for (int pdg : particle_sequence) {
      particlesByPDG[pdg] = {};
    }
    for (auto& p : c12->getDetParticles()) {
      particlesByPDG[p->getPid()].push_back(p);
    }
    if (verbose)
      for (int pdg : particle_sequence)
        std::cout << std::format("{:<8s}: {:>2d}", pdg_name(pdg), particlesByPDG[pdg].size()) << std::endl;

    // --------------------------- 1.1 ---------------------------
    if (particlesByPDG.at(11).empty()) {
      if (verbose) std::cout << "No electrons in event. Skipping event..." << std::endl;
      continue;
    }

    // =========================== 2. ===========================
    // process all e- first, then nucleons, then mesons
    if (verbose) std::cout << "Processing particles..." << std::endl;
    for (int pdg : particle_sequence) {
      // skip early if there are no particles of this pdg in the event
      if (particlesByPDG[pdg].empty()) {
        // if (verbose) std::cout << "No " << pdg_name(pdg) << " in event. Skipping ..." << std::endl;
        continue;
      }

      // --------------------------- 2.1 ---------------------------
      // sort particles by momentum for each PDG code (if there are more than one particles to sort)
      if (particlesByPDG[pdg].size() > 1) {
        if (verbose)
          std::cout << "Sorting particles: " << pdg_name(pdg) << " " << particlesByPDG[pdg].size() << std::endl;
        std::sort(std::begin(particlesByPDG[pdg]), std::end(particlesByPDG[pdg]),
                  [](clas12::region_part_ptr a, clas12::region_part_ptr b) { return a->getP() > b->getP(); });
      }

      // --------------------------- 2.2 ---------------------------
      // after we processed electrons (pdg!=11), do we have electrons left? No? Reject event.
      if (pdg != 11 && particlesByPDG[11].empty()) {
        if (verbose) std::cout << "No valid electrons in event. Rejecting event..." << std::endl;
        break;
      }

      // --------------------------- 2.3 ---------------------------
      if (pdg != 11 && std::isnan(reference_vertex)) {
        // set reference_vertex to highest momentum electron
        // this is only done once, after all electrons have been processed,
        // and thus are only valid electrons
        if (verbose) std::cout << "Reference Vertex: " << std::flush;
        reference_vertex = particlesByPDG[11].at(0)->par()->getVz();
        if (verbose) std::cout << reference_vertex << std::endl;
      }

      // --------------------------- 3.1 ---------------------------
      // // ====================== PARTICLE CUTS ======================
      // check cuts with selector functions, for all particles first
      if (pdg == 11) {
        particlesByPDG[pdg].erase(std::remove_if(particlesByPDG[pdg].begin(), particlesByPDG[pdg].end(),
                                                 [&](auto&& p) { return !selectors::electron(p, inbending); }),
                                  particlesByPDG[pdg].end());
      }
      if (pdg == 2212) {
        particlesByPDG[pdg].erase(
            std::remove_if(particlesByPDG[pdg].begin(), particlesByPDG[pdg].end(),
                           [&](auto&& p) { return !selectors::proton(p, inbending, reference_vertex); }),
            particlesByPDG[pdg].end());
      }
      if (pdg == 2112) {
        particlesByPDG[pdg].erase(
            std::remove_if(particlesByPDG[pdg].begin(), particlesByPDG[pdg].end(),
                           [&](auto&& p) { return !selectors::neutron(p, inbending, reference_vertex); }),
            particlesByPDG[pdg].end());
      }
      if (pdg == 211) {
        particlesByPDG[pdg].erase(
            std::remove_if(particlesByPDG[pdg].begin(), particlesByPDG[pdg].end(),
                           [&](auto&& p) { return !selectors::piplus(p, inbending, reference_vertex); }),
            particlesByPDG[pdg].end());
      }
      if (pdg == -211) {
        particlesByPDG[pdg].erase(
            std::remove_if(particlesByPDG[pdg].begin(), particlesByPDG[pdg].end(),
                           [&](auto&& p) { return !selectors::piminus(p, inbending, reference_vertex); }),
            particlesByPDG[pdg].end());
      }
      if (pdg == 321) {
        particlesByPDG[pdg].erase(
            std::remove_if(particlesByPDG[pdg].begin(), particlesByPDG[pdg].end(),
                           [&](auto&& p) { return !selectors::Kplus(p, inbending, reference_vertex); }),
            particlesByPDG[pdg].end());
      }
      if (pdg == -321) {
        particlesByPDG[pdg].erase(
            std::remove_if(particlesByPDG[pdg].begin(), particlesByPDG[pdg].end(),
                           [&](auto&& p) { return !selectors::Kminus(p, inbending, reference_vertex); }),
            particlesByPDG[pdg].end());
      }
      if (pdg == 22) {
        particlesByPDG[pdg].erase(std::remove_if(particlesByPDG[pdg].begin(), particlesByPDG[pdg].end(),
                                                 [&](auto&& p) { return !selectors::photon(p, inbending); }),
                                  particlesByPDG[pdg].end());
      }

      // ========================== 3. ==========================
      // loop particles of this type
      for (auto&& p : particlesByPDG[pdg]) {
        std::string name = pdg_name(p->getPid());
        if (verbose) std::cout << std::format("Valid {} found!", name) << std::endl;
        particle_count++;

        // --------------------------- 3.2 ---------------------------
        // ===========================================================
        //                        Collect data
        // ===========================================================
        TParticlePDG* particle_pdg = pdg_db->GetParticle(p->getPid());
        if (particle_pdg == nullptr) continue;

        // ====================== MOMENTUM FIELDS ======================
        if (vars.fParticlesOfInterestForMomentumFields.contains(p->getPid())) {
          vars.AppendValue("p4_" + name + "_px", static_cast<double>(p->par()->getPx()));
          vars.AppendValue("p4_" + name + "_py", static_cast<double>(p->par()->getPy()));
          vars.AppendValue("p4_" + name + "_pz", static_cast<double>(p->par()->getPz()));
          vars.AppendValue(
              "p4_" + name + "_E",
              static_cast<double>(std::sqrt(std::pow(p->par()->getPx(), 2) + std::pow(p->par()->getPy(), 2) +
                                            std::pow(p->par()->getPz(), 2) + std::pow(particle_pdg->Mass(), 2))));
        }
        // ====================== DETECTOR FIELDS ======================
        if (vars.fParticlesOfInterestForDetectorFields.contains(p->getPid())) {
          // Old way used  1000<=abs(part_status)<2000 => FT: ele_det=1
          //               2000<=abs(part_Status)<4000 => FD: ele_det=2
          //               4000<=abs(part_Status)      => CD: ele_det=3
          // There might be a mismatch to clas12::FD (et.al.) constants, where clas12::CD=3000...
          // Keep the status check (instead of checking p->getRegion() against clas12 constants), and keep this mismatch
          // in mind.
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
    }

    // =========================== 4. ===========================
    file->Fill();
  }
  fmt::print("Processed a total of {} events with {} particles in {:.1f} seconds.\n", count, particle_count,
             duration_cast<milliseconds>(steady_clock::now() - t0).count() / 1000.);

  StatisticsCollector::print_hierarchical_statistics();

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
  std::cout << std::format(" {:>2.1f}s", duration_cast<milliseconds>(steady_clock::now() - t3).count() / 1000.)
            << std::endl;

  std::cout << "done." << std::endl;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: new_filter <input_file> [output_file] [isInbending] [num_events]" << std::endl;
    return 1;
  }

  std::string input_file = argv[1];
  std::string output_file = (argc >= 3) ? argv[2] : "output.root";
  bool isInbending = (argc >= 4) ? bool(argv[3]) : true;
  uint num_events = (argc >= 5) ? std::stoi(argv[4]) : 0;

  new_filter(input_file, output_file, isInbending, num_events);
  return 0;
}