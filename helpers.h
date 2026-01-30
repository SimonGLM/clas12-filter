#pragma once
#include <string>

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