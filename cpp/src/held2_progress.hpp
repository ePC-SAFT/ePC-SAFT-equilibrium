#pragma once

#include <iosfwd>
#include <limits>
#include <string>

namespace epcsaft_equilibrium {

enum class Held2ProgressKind {
    CaseStart,
    ReferenceStart,
    ReferenceRoot,
    ReferenceSelected,
    StageStart,
    StageIEvaluation,
    StageIIUpper,
    LocalIteration,
    Certificate,
    StageSkipped,
    Failure,
    Final,
};

struct Held2ProgressEvent {
    Held2ProgressKind kind = Held2ProgressKind::Failure;
    std::string case_id;
    std::string stage;
    std::string status;
    std::string reason;
    std::string mechanical_class;
    int iteration = -1;
    int major_iteration = -1;
    int attempt = -1;
    int count = -1;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    double objective = 0.0;
    double lower_bound = std::numeric_limits<double>::quiet_NaN();
    double upper_bound = std::numeric_limits<double>::quiet_NaN();
    double gap = std::numeric_limits<double>::quiet_NaN();
    double volume = 0.0;
    double pressure_residual = 0.0;
    double physical_total_ion_mole_fraction =
        std::numeric_limits<double>::quiet_NaN();
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    double primal_residual = std::numeric_limits<double>::quiet_NaN();
    double dual_residual = std::numeric_limits<double>::quiet_NaN();
    double complementarity = std::numeric_limits<double>::quiet_NaN();
    double step_norm = std::numeric_limits<double>::quiet_NaN();
    double dual_step = std::numeric_limits<double>::quiet_NaN();
    double primal_step = std::numeric_limits<double>::quiet_NaN();
    int line_search_steps = -1;
};

class Held2ProgressObserver {
public:
    virtual ~Held2ProgressObserver() = default;
    virtual void observe(const Held2ProgressEvent& event) = 0;
};

inline void observe_held2(
    Held2ProgressObserver* observer,
    const Held2ProgressEvent& event
) noexcept {
    if (observer != nullptr) {
        try {
            observer->observe(event);
        } catch (...) {
            // Progress reporting is observational and cannot alter solver flow.
        }
    }
}

class Held2TerminalProgress final : public Held2ProgressObserver {
public:
    explicit Held2TerminalProgress(std::ostream& output) noexcept;

    void observe(const Held2ProgressEvent& event) override;
    [[nodiscard]] bool healthy() const noexcept;

private:
    std::ostream& output_;
    bool healthy_ = true;
};

}  // namespace epcsaft_equilibrium
