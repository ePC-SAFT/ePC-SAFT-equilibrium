#include "flash.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    try {
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
        std::cout << "manufactured HELD2 Steps 1-10 accepted\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "epcsaft-equilibrium-manufactured-check: "
                  << error.what() << '\n';
        return 1;
    }
}
