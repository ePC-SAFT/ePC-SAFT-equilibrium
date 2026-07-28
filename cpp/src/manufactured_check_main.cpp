#include "flash.hpp"
#include "held2_algorithm.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_equilibrium::test {
Held2StateEvaluation evaluate_manufactured_state(
    const Held2Coordinates&,
    const std::vector<double>&,
    double
);
void run_held2_step1_checks();
std::string run_held2_step2_checks(Held2ProgressObserver*);
void run_held2_step3_checks();
void run_held2_step4_checks();
void run_held2_step5_checks();
void run_held2_step6_checks();
void run_held2_step7_checks();
void run_held2_step8_checks();
void run_held2_step9_checks();
void run_held2_step10_checks();
}

namespace {

void run_workflow_check(bool trace) {
    using namespace epcsaft_equilibrium;
    const std::vector<double> charges{0.0, 1.0, -1.0};
    const Held2Coordinates coordinates = make_held2_coordinates(charges);
    Held2TerminalProgress progress(std::cout);
    const Held2AlgorithmResult result = run_held2_algorithm(
        {
            {"neutral", "cation", "anion"},
            charges,
            [coordinates](const auto& composition, double log_volume) {
                Held2StateEvaluation state = test::evaluate_manufactured_state(
                    coordinates, composition, log_volume
                );
                state.pressure_stationarity_relative *= -1.0;
                state.pressure_stationarity_derivative_log_volume *= -1.0;
                return state;
            },
            [](const auto&) {
                return std::array<double, 2>{0.5, 1.5};
            },
            [](const auto& composition, double) {
                return composition.front();
            },
            std::numeric_limits<double>::quiet_NaN(),
            {},
            {},
        },
        {298.15, 100000.0, {0.5, 0.25, 0.25}},
        {200, 20, 10},
        trace ? &progress : nullptr
    );
    if (
        result.outcome != "physical_equilibrium_accepted"
        || !result.step10
        || result.step10->status != "complete"
        || result.phases.size() != 2
        || !result.step10->final_certificate->accepted
        || result.globality_certificate != "not_guaranteed"
    ) {
        throw std::runtime_error(
            "manufactured HELD2 Steps 1-10 failed: "
            + result.failure_stage + ": " + result.failure_reason
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
                "usage: --case {step1|step2|step3|step4|step5|step6|step7|step8|step9|step10|workflow} [--trace]"
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
        } else if (check == "step3") {
            epcsaft_equilibrium::test::run_held2_step3_checks();
            std::cout << "step3: pass\n";
        } else if (check == "step4") {
            epcsaft_equilibrium::test::run_held2_step4_checks();
            std::cout << "step4: pass\n";
        } else if (check == "step5") {
            epcsaft_equilibrium::test::run_held2_step5_checks();
            std::cout << "step5: pass\n";
        } else if (check == "step6") {
            epcsaft_equilibrium::test::run_held2_step6_checks();
            std::cout << "step6: pass\n";
        } else if (check == "step7") {
            epcsaft_equilibrium::test::run_held2_step7_checks();
            std::cout << "step7: pass\n";
        } else if (check == "step8") {
            epcsaft_equilibrium::test::run_held2_step8_checks();
            std::cout << "step8: pass\n";
        } else if (check == "step9") {
            epcsaft_equilibrium::test::run_held2_step9_checks();
            std::cout << "step9: pass\n";
        } else if (check == "step10") {
            epcsaft_equilibrium::test::run_held2_step10_checks();
            std::cout << "step10: pass\n";
        } else if (check == "workflow") {
            run_workflow_check(trace);
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
