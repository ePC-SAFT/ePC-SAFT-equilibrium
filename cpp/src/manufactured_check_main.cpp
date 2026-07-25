#include "flash.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_equilibrium::test {
void run_held2_step1_checks();
std::string run_held2_step2_checks(Held2ProgressObserver*);
}

namespace {

void run_workflow_check() {
    const epcsaft_equilibrium::Held2ManufacturedWorkflowResult result =
        epcsaft_equilibrium::solve_held2_manufactured_workflow(
            {0.0, 1.0, -1.0},
            {0.5, 0.25, 0.25}
        );
    if (
        result.workflow.outcome != "physical_equilibrium_accepted"
        || result.workflow.completed_step != 10
        || result.workflow.transitions.size() != 4
        || result.stage_i.outcome != "negative_witness_found"
        || result.stage_ii.outcome != "candidate_set"
        || result.stage_iii.physical_status != "accepted"
        || !std::isfinite(result.stage_iii.objective)
        || result.stage_iii.modified_balance_inf_norm >= 1.0e-9
        || result.stage_iii.ordinary_balance_inf_norm >= 1.0e-9
        || result.stage_iii.pressure_stationarity_inf_norm >= 1.0e-9
        || result.stage_iii.modified_potential_mixed_gap >= 1.0e-9
        || result.globality_certificate != "not_guaranteed"
    ) {
        throw std::runtime_error(
            "manufactured HELD2 workflow failed its scientific checks"
        );
    }
    std::cout << "workflow: pass\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 1
            && (argc < 3 || argc > 4
                || std::string(argv[1]) != "--case"
                || (argc == 4 && std::string(argv[3]) != "--trace"))) {
            throw std::invalid_argument(
                "usage: --case {step1|step2|workflow} [--trace]"
            );
        }
        const std::string check = argc == 1 ? "workflow" : argv[2];
        const bool trace = argc == 4;
        if (check == "step1") {
            epcsaft_equilibrium::test::run_held2_step1_checks();
            std::cout << "step1: pass\n";
        } else if (check == "step2") {
            epcsaft_equilibrium::Held2TerminalProgress progress(std::cout);
            const std::string hash =
                epcsaft_equilibrium::test::run_held2_step2_checks(
                    trace ? &progress : nullptr
                );
            std::cout << "step2: pass result_hash=" << hash << '\n';
        } else if (check == "workflow") {
            run_workflow_check();
        } else {
            throw std::invalid_argument("unknown manufactured check case");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "epcsaft-equilibrium-manufactured-check: "
                  << error.what() << '\n';
        return 1;
    }
}
