#include <clas12reader.h>

#include <chrono>
#include <map>

#include "helpers.h"
#include "particle_selector.h"

using namespace std::chrono;
void count(const char* inFile) {
  const std::vector<int> particle_sequence = {11, 2212, 2112, 211, -211, 321, -321, 22};

  auto c12 = std::make_unique<clas12::clas12reader>(inFile);

  std::map<int, int> count;
  for (int pdg : particle_sequence) count[pdg] = 0;
  while (c12->next()) {
    std::unordered_map<int, std::vector<clas12::region_part_ptr>> particlesByPDG;
    for (auto& p : c12->getDetParticles()) {
      count[p->getPid()]++;
    }
  }
  for (auto&& [pdg, cnt] : count) {
    std::cout << std::format("PDG code {} ({}): {}\n", pdg, pdg_name(pdg), cnt);
  }
}