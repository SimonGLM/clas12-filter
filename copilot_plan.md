# Plan: Context-Based Cut Tracking System

**Goal:** Make the cut statistics system maintainable for handoff while adding per-particle debugging capability.

## Problem Statement

Current decorator-based system has:

- Hidden global state (`call_stack()`) that's hard to understand
- No per-particle tracking (can't see which cuts each particle failed)
- No mode switching (can't choose early-return vs complete-evaluation)
- Preprocessor hack considered for mode switching (`#define custom_return`)

## Solution: Explicit Context Object

Replace decorators with `ParticleContext` that:

1. Tracks all cut results for each particle
2. Supports two evaluation modes (early-return vs complete-trace)
3. Makes data flow explicit (passed as parameter)
4. Separates particle-level tracking from global statistics

## Architecture

```
ParticleContext (per particle)
    ├─ EvaluationMode: EarlyReturn | CompleteTrace
    ├─ cut_history: vector<CutResult>
    └─ apply_cut(name, func) → bool (continue?)

StatisticsCollector (global, optional)
    ├─ Global counters for all cuts
    └─ print_summary() → current format output

Cut Functions (unchanged logic)
    └─ Pure functions: take particle, return bool
```

## Implementation Steps

### 1. Create `particle_context.h`

**New file** with:

- `enum class EvaluationMode { EarlyReturn, CompleteTrace }`
- `class ParticleContext` with:
  - `apply_cut(name, lambda)` or `check(name, func, args...)`
  - `passed()` returns overall result
  - `print_trace()` shows cut history
  - `get_failed_cuts()` for programmatic access

**Key behavior:**

```cpp
bool apply_cut(const char* name, auto&& cut_func) {
  bool result = cut_func();
  cut_history.push_back({name, result});
  
  if (!result) {
    overall_passed = false;
    return mode == EvaluationMode::CompleteTrace; // continue if tracing
  }
  return true; // continue
}
```

### 2. Simplify `cuts.h`

**Delete:**

- `StatisticsDecorator` class (entire ~130 lines)
- `DecoratedCut` template class
- `child_stats` tracking mechanism
- `call_stack()` static function

**Keep:**

- Namespace structure (`cuts::generic::`, `cuts::FD::`, etc.)
- Function declarations pointing to implementations
- Clean, simple function pointers

### 3. Keep `cuts.cpp` unchanged

- All cut implementations stay exactly as they are
- Functions remain pure: `bool _some_cut(particle*, args...)`
- No logic changes → no risk

### 4. Create `statistics_collector.h`

**New file** with:

- Singleton or static class for aggregate statistics
- `record_cut(name, result)` → increment counters
- `print_summary()` → current hierarchical output format
- Optional: `ParticleContext` auto-reports if enabled

### 5. Refactor `particle_selector.h`

**For each selector (`_electron`, `_proton`, etc.):**

**Current pattern:**

```cpp
bool _electron(clas12::region_particle* p, bool inbending, int tightness) {
  if (!cuts::generic::PID_cut(p, 11)) return false;
  if (!cuts::generic::charge_cut(p, -1)) return false;
  // ...
  return true;
}
```

**New pattern with helper macro:**

```cpp
#define CHECK_CUT(ctx, name, ...) \
  if (!ctx.check(name, __VA_ARGS__)) return ctx.passed()

bool _electron(ParticleContext& ctx, clas12::region_particle* p, 
               bool inbending, int tightness) {
  CHECK_CUT(ctx, "PID_cut(11)", cuts::generic::impl::_PID_cut, p, 11);
  CHECK_CUT(ctx, "charge_cut(-1)", cuts::generic::impl::_charge_cut, p, -1);
  // ... more CHECK_CUT calls
  return ctx.passed();
}
```

**Why this works:**

- EarlyReturn mode: `ctx.passed()` returns `false` immediately → early exit
- CompleteTrace mode: `ctx.passed()` returns `true` → continues to next cut
- Final `return ctx.passed()` gives correct overall result

### 6. Update call sites (`filter_clas12_hipo.cxx` or similar)

**Before:**

```cpp
if (selectors::impl::electron(p, inbending, tightness)) {
  // accepted
}
```

**After:**

```cpp
ParticleContext ctx(particle, EvaluationMode::EarlyReturn);
if (selectors::impl::electron(ctx, p, inbending, tightness)) {
  // accepted
} else {
  // Debug: ctx.print_trace() to see all failures
}

// At end of run:
StatisticsCollector::print_summary();
```

### 7. Fix include structure

**In `new_filter.cxx` (or similar):**

- Remove `#include "cuts.cpp"` (line 45)
- Add proper compilation: compile `cuts.cpp` separately and link
- Keep normal `#include "cuts.h"`

### 8. Cleanup

- Delete commented code blocks in `cuts.h` (lines 65-88, 103-123)
- Remove unused `bool passed = true;` variables in selectors
- Add brief comments explaining context usage

## Configuration Options

**Runtime mode switch:**

```cpp
auto mode = config.debug_trace ? EvaluationMode::CompleteTrace 
                                : EvaluationMode::EarlyReturn;
```

**Compile-time mode (zero overhead):**

```cpp
#ifdef DEBUG_TRACE_CUTS
  constexpr auto mode = EvaluationMode::CompleteTrace;
#else
  constexpr auto mode = EvaluationMode::EarlyReturn;
#endif
```

## Verification Plan

1. **Compile:** Ensure clean build with new structure
2. **Run:** Process test dataset, verify same acceptance counts
3. **Statistics:** Compare `StatisticsCollector::print_summary()` output to current
4. **Per-particle debug:** Inject failure, verify `ctx.print_trace()` shows all cuts
5. **Performance:** Profile both modes, ensure EarlyReturn has no overhead

## Benefits

### Maintainability (Primary Goal)

- **Explicit data flow:** Context passed visibly, no hidden globals
- **Self-documenting:** `CHECK_CUT(ctx, "cut_name", ...)` is clear
- **Easier to explain:** "Pass context, it tracks results" vs decorator magic
- **Testable:** Mock context for unit tests

### Per-Particle Debugging (Requested Feature)

```
Particle #1234 rejected:
  ✓ PID_cut(11)
  ✓ charge_cut(-1)
  ✗ HTCC_nphe_cut          <- First failure
  ✗ EC_sampling_fraction    <- Also failed
  ✓ DC_fiducial_region1
  ... (complete trace)
```

### Mode Switching (Clean Solution)

- No preprocessor hacks redefining `return`
- Type-safe, debuggable control flow
- Can vary per-particle if needed
- IDE-friendly code

### Preserved Features

- Global statistics (moved to StatisticsCollector)
- Hierarchical output format (selector → cuts)
- All cut logic unchanged
- Same performance in EarlyReturn mode

## Design Decisions

| Decision | Rationale |
|----------|-----------|
| Context over Decorator | Explicit > implicit for maintainability |
| Mode as enum parameter | Clear, type-safe, runtime or compile-time |
| Separate StatisticsCollector | Single responsibility, optional aggregation |
| Keep cut implementations | Risk mitigation, proven logic untouched |
| Macro for CHECK_CUT | Reduces boilerplate, standard pattern |

## Migration Notes

- **Breaking change:** Selector signatures gain `ParticleContext&` parameter
- **Call sites:** Need to create context before calling selectors
- **Statistics:** Access via `StatisticsCollector` instead of decorator registry
- **Backward compat:** Could wrap old API temporarily if needed

## Future Extensions

Once implemented, easy to add:

- Cut timing/profiling per particle
- Conditional cuts (skip if another failed)
- Cut result caching within event
- Custom logging/output formats
- Machine learning feature extraction from cut results
