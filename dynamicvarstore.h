#pragma once
#include <ROOT/RNTupleModel.hxx>

#include "helpers.h"

// A helper class that wraps an RNTupleModel and allows building it up dynamically
// by adding fields of various types identified by their names as strings.
// Also provides type-safe setting of the values of the fields via templated SetValue<T>() method.
// DynamicVarStore::MoveModel() must be used to move the ownership of the built model into
// the RNTupleWriter.

class DynamicVarStore {
 public:
  DynamicVarStore(std::unique_ptr<ROOT::RNTupleModel> model, bool verbose = false);
  ~DynamicVarStore() = default;
  std::unique_ptr<ROOT::RNTupleModel>& Model();
  template <typename T>
  void AddField(const std::string& name);
  template <typename T>
  void SetValue(const std::string& name, const T& value);
  template <typename T>
  void AppendValue(const std::string& name, const T& value);
  void ResetVectorFields();

  void AddMomentaFieldsPerParticle(const std::vector<int>& particles);
  void AddDetectorFieldsPerParticle(const std::vector<int>& particles);
  void AddSectorFieldsPerParticle(const std::vector<int>& particles);
  void AddChi2PidFieldsPerParticle(const std::vector<int>& particles);
  void AddDCFields();
  std::set<int> fParticlesOfInterestForMomentumFields{};
  std::set<int> fParticlesOfInterestForDetectorFields{};
  std::set<int> fParticlesOfInterestForSectorFields{};
  std::set<int> fParticlesOfInterestForChi2PidFields{};
  std::set<int> fParticlesOfInterestForDCFields{};

 private:
  bool IsField(const std::string&);
  bool verbose = false;
  // Each entry stores a shared_ptr<T> of one of the supported types
  using variant_type = std::variant<std::shared_ptr<int>, std::shared_ptr<float>, std::shared_ptr<double>,
                                    std::shared_ptr<std::vector<int>>, std::shared_ptr<std::vector<float>>,
                                    std::shared_ptr<std::vector<double>>>;  // add more vector types here if needed
  using variant_map = std::unordered_map<std::string, variant_type>;

  variant_map fMap;
  std::unique_ptr<ROOT::RNTupleModel> fModel;
};

// Constructor that takes ownership of the model
DynamicVarStore::DynamicVarStore(std::unique_ptr<ROOT::RNTupleModel> model, bool verbose)
    : fModel(std::move(model)), verbose(verbose) {
  if (this->verbose) std::cout << "[DynamicVarStore] Created with model." << std::endl;
}

// Check if name is a field registered by DVS
bool DynamicVarStore::IsField(const std::string& name) { return fMap.contains(name); }

// Returns the reference to the Model pointer, so that ownership can be moved
std::unique_ptr<ROOT::RNTupleModel>& DynamicVarStore::Model() { return fModel; }

// Build the model by adding fields of type T with the given name
template <typename T>
void DynamicVarStore::AddField(const std::string& name) {
  if (this->verbose) std::cout << std::format("[DynamicVarStore] Adding field '{}'.", name) << std::endl;
  // Create field in the model and store the shared_ptr<T> in the map
  fMap[name] = fModel->MakeField<T>(name);
}

// Set the content of the memory managed by the RNTupleModel
template <typename T>
void DynamicVarStore::SetValue(const std::string& name, const T& value) {
  if (this->verbose) std::cout << std::format("[DynamicVarStore] Setting field value '{}' .", name) << std::endl;
  if (!IsField(name)) throw std::invalid_argument(std::format("Field '{}' does not exist in DynamicVarStore", name));
  // Get the shared_ptr<T> from the map, dereference it and set the value
  // no `std::get<...>(...).reset(value)` here as that would change the address that is pointed to
  *std::get<std::shared_ptr<T>>(fMap.at(name)) = value;
}

// Append a value to the vector in the RNTupleModel managed memory
template <typename T>
void DynamicVarStore::AppendValue(const std::string& name, const T& value) {
  if (!IsField(name)) throw std::invalid_argument(std::format("Field '{}' does not exist in DynamicVarStore", name));
  // Get the shared_ptr<T> from the map, check type, then push_back
  std::visit(
      [&](auto&& arg) {
        using U = std::decay_t<decltype(arg)>;
        if constexpr (!std::is_same_v<U, std::shared_ptr<std::vector<T>>>) {
#if __GNUG__  // try to get sensible type names (only GLib)
          throw std::invalid_argument(std::format("Type missmatch between field '{}' and value of type '{}'", name,
                                                  abi::__cxa_demangle(typeid(arg).name(), NULL, NULL, nullptr)));
#else  // fall back to not providing a demangled type name
          throw std::invalid_argument(std::format("Type missmatch between field '{}' and value type", name));
#endif
        } else {
          arg->push_back(value);
        }
      },
      fMap.at(name));
}

// Reset the vectors
void DynamicVarStore::ResetVectorFields() {
  if (this->verbose) std::cout << "[DynamicVarStore] Resetting vector fields." << std::endl;
  // Get the shared_ptr<std::vector<T>> from the map and clear the vector
  for (const auto& [name, ptr] : fMap) {
    std::visit(
        [&](auto&& arg) {
          using U = std::decay_t<decltype(arg)>;
          if constexpr (std::is_same_v<U, std::shared_ptr<std::vector<int>>> ||
                        std::is_same_v<U, std::shared_ptr<std::vector<double>>>)  // also add more vector types here if
                                                                                  // needed
          {
            arg->clear();
          }
        },
        ptr);
  }
}

// Add various vector fields to the RNTupleModel for particles
void DynamicVarStore::AddMomentaFieldsPerParticle(const std::vector<int>& particles) {
  for (int part : particles) {
    for (std::string&& var : {"px", "py", "pz", "E"}) {
      AddField<std::vector<double>>("p4_" + pdg_name(part) + "_" + var);
    }
    fParticlesOfInterestForMomentumFields.insert(part);
  }
}

void DynamicVarStore::AddDetectorFieldsPerParticle(const std::vector<int>& particles) {
  for (int part : particles) {
    AddField<std::vector<int>>(pdg_name(part) + "_det");
    fParticlesOfInterestForDetectorFields.insert(part);
  }
}

void DynamicVarStore::AddSectorFieldsPerParticle(const std::vector<int>& particles) {
  for (int part : particles) {
    AddField<std::vector<int>>(pdg_name(part) + "_sec");
    fParticlesOfInterestForSectorFields.insert(part);
  }
}

void DynamicVarStore::AddChi2PidFieldsPerParticle(const std::vector<int>& particles) {
  for (int part : particles) {
    AddField<std::vector<double>>("p4_" + pdg_name(part) + "_chi2pid");
    fParticlesOfInterestForChi2PidFields.insert(part);
  }
}

void DynamicVarStore::AddDCFields() {
  for (const char& chord : {'x', 'y', 'z'}) AddField<std::vector<double>>(std::format("p4_prot_dc{}1", chord));
  for (const int& i : {1, 2, 3}) AddField<std::vector<double>>(std::format("p4_ele_dcedge{}", i));
  fParticlesOfInterestForDCFields.insert(2212);
  fParticlesOfInterestForDCFields.insert(11);
}
