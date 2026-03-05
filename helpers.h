#pragma once
#include <string>

const std::string pdg_name(const int& pdg) {
  return pdg == 11     ? "ele"
         : pdg == 22   ? "phot"
         : pdg == 2212 ? "prot"
         : pdg == 2112 ? "neutr"
         : pdg == 211  ? "pip"
         : pdg == -211 ? "pim"
         : pdg == 321  ? "Kp"
         : pdg == -321 ? "Km"
                       : std::format("pdg: {}", pdg);
};

class not_implemented_error : public std::logic_error {
 public:
  explicit not_implemented_error(const std::string& what_arg) : std::logic_error(what_arg) {}
};
