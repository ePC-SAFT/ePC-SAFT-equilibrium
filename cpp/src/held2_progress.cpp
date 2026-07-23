#include "held2_progress.hpp"

#include <cmath>
#include <iomanip>
#include <ostream>

namespace epcsaft_equilibrium {
namespace {

void write_stage_heading(
    std::ostream& output,
    const Held2ProgressEvent& event,
    const char* suffix = nullptr
) {
    output << '\n' << event.stage;
    if (suffix != nullptr) {
        output << " - " << suffix;
    }
    output << '\n';
}

}  // namespace

Held2TerminalProgress::Held2TerminalProgress(std::ostream& output) noexcept
    : output_(output) {}

void Held2TerminalProgress::observe(const Held2ProgressEvent& event) {
    if (!healthy_) {
        return;
    }
    try {
        output_ << std::scientific << std::setprecision(6);
        switch (event.kind) {
        case Held2ProgressKind::CaseStart:
            output_ << "HELD2.0  case=" << event.case_id
                    << "  T=" << event.temperature_k << " K"
                    << "  P=" << event.pressure_pa << " Pa\n";
            break;
        case Held2ProgressKind::ReferenceStart:
            output_ << "\nREFERENCE PRESSURE ROOTS\n";
            break;
        case Held2ProgressKind::ReferenceRoot:
            output_ << "  root " << std::setw(3) << event.count
                    << "  V=" << std::setw(13) << event.volume
                    << "  P rel.=" << std::setw(13) << event.pressure_residual
                    << "  G/RT=" << std::setw(13) << event.objective
                    << "  " << event.mechanical_class << '\n';
            break;
        case Held2ProgressKind::ReferenceSelected:
            output_ << "  selected root " << event.count
                    << "  status=" << event.status << '\n';
            break;
        case Held2ProgressKind::StageStart:
            write_stage_heading(output_, event);
            break;
        case Held2ProgressKind::StageIEvaluation:
            output_ << "  eval " << std::setw(5) << event.iteration;
            if (event.status == "failed") {
                output_ << "  failed  reason=" << event.reason << '\n';
            } else {
                output_ << "  TPD=" << std::setw(13) << event.objective
                        << "  P rel.=" << std::setw(13)
                        << event.pressure_residual
                        << "  V=" << std::setw(13) << event.volume;
                if (std::isfinite(event.physical_total_ion_mole_fraction)
                    && std::isfinite(event.total_ion_mole_fraction_max)) {
                    output_ << "  x_ion=" << std::setw(13)
                            << event.physical_total_ion_mole_fraction
                            << "  limit=" << std::setw(13)
                            << event.total_ion_mole_fraction_max;
                }
                output_ << "  " << event.status << '\n';
            }
            break;
        case Held2ProgressKind::StageIIUpper:
            output_ << "  major " << std::setw(3) << event.major_iteration
                    << "  upper=" << std::setw(13) << event.upper_bound;
            if (std::isfinite(event.lower_bound)) {
                output_ << "  lower=" << std::setw(13) << event.lower_bound;
            } else {
                output_ << "  lower=unavailable";
            }
            if (std::isfinite(event.gap)) {
                output_ << "  gap=" << std::setw(13) << event.gap;
            } else {
                output_ << "  gap=unavailable";
            }
            output_
                    << "  primal=" << std::setw(13) << event.primal_residual
                    << "  dual=" << std::setw(13) << event.dual_residual
                    << "  active_cuts=" << event.count
                    << "  " << event.status << '\n';
            break;
        case Held2ProgressKind::LocalIteration:
            output_ << "  " << event.stage;
            if (event.major_iteration >= 0) {
                output_ << "  major=" << event.major_iteration;
            }
            if (event.attempt >= 0) {
                output_ << "  attempt=" << event.attempt;
            }
            output_ << "  iter=" << std::setw(5) << event.iteration
                    << "  objective=" << std::setw(13) << event.objective
                    << "  primal=" << std::setw(13) << event.primal_residual
                    << "  dual=" << std::setw(13) << event.dual_residual
                    << "  comp=" << std::setw(13) << event.complementarity
                    << "  step=" << std::setw(13) << event.step_norm
                    << "  alpha=" << std::setw(13) << event.primal_step
                    << "  ls=" << event.line_search_steps
                    << '\n';
            break;
        case Held2ProgressKind::Certificate:
            output_ << "  certificate";
            if (!event.stage.empty()) {
                output_ << "  stage=" << event.stage;
            }
            if (event.major_iteration >= 0) {
                output_ << "  major=" << event.major_iteration;
            }
            if (event.attempt >= 0) {
                output_ << "  attempt=" << event.attempt;
            }
            output_ << "  " << event.status;
            if (!event.reason.empty()) {
                output_ << "  reason=" << event.reason;
            }
            if (std::isfinite(event.primal_residual)) {
                output_ << "  primal=" << std::setw(13)
                        << event.primal_residual;
            }
            if (std::isfinite(event.dual_residual)) {
                output_ << "  dual=" << std::setw(13)
                        << event.dual_residual;
            }
            if (std::isfinite(event.complementarity)) {
                output_ << "  comp=" << std::setw(13)
                        << event.complementarity;
            }
            if (std::isfinite(event.gap)) {
                output_ << "  objective=" << std::setw(13)
                        << event.objective
                        << "  upper=" << std::setw(13)
                        << event.upper_bound
                        << "  gap=" << std::setw(13)
                        << event.gap;
            }
            output_ << '\n';
            break;
        case Held2ProgressKind::StageSkipped:
            write_stage_heading(output_, event, "SKIPPED");
            output_ << "  reason=" << event.reason << '\n';
            break;
        case Held2ProgressKind::Failure:
            output_ << "  FAILURE  stage=" << event.stage
                    << "  reason=" << event.reason << '\n';
            break;
        case Held2ProgressKind::Final:
            output_ << "\nFINAL  status=" << event.status;
            if (!event.reason.empty()) {
                output_ << "  reason=" << event.reason;
            }
            output_ << '\n';
            break;
        }
        output_.flush();
        if (!output_) {
            healthy_ = false;
        }
    } catch (...) {
        healthy_ = false;
    }
}

bool Held2TerminalProgress::healthy() const noexcept {
    return healthy_;
}

}  // namespace epcsaft_equilibrium
