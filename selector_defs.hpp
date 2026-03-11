#pragma once
#include "clas12defs.h"
#include "cut_defs.hpp"
#include "region_particle.h"
#include "selector.hpp"

namespace selectors {
  inline Selector<clas12::region_particle*> make_electron_selector(EvalMode mode = EvalMode::EarlyStopping) {
    Selector<clas12::region_particle*> selector("electron", mode);

    selector.add({"PID_cut", [&](clas12::region_particle* p) { return cuts::generic::PID_cut(p, 22); }});
    selector.add({"charge_cut", [&](clas12::region_particle* p) { return cuts::generic::charge_cut(p, -1); }});

    selector.add({"Forward Detector", [&](clas12::region_particle* p) {

                  }})
  }
}  // namespace selectors