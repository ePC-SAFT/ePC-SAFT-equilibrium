#include "flash.hpp"
#include "held2_progress.hpp"
#include "provider.hpp"
#include "result_json.hpp"

#include <epcsaft/native_model_v1.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct Options {
    std::string model_config;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> feed;
    std::string output_path;
    bool trace = false;
};

[[noreturn]] void usage_error(const std::string& message) {
    throw std::invalid_argument(
        message
        + "\nusage: epcsaft-equilibrium-diagnostic"
          " --model-config PATH --temperature K --pressure PA"
          " --feed X1,X2[,X3...] [--trace] [--output PATH]"
    );
}

double parse_double(std::string_view text, const char* name) {
    std::size_t consumed = 0;
    double value = 0.0;
    try {
        value = std::stod(std::string(text), &consumed);
    } catch (const std::exception&) {
        usage_error(std::string(name) + " must be a number");
    }
    if (consumed != text.size() || !std::isfinite(value)) {
        usage_error(std::string(name) + " must be a finite number");
    }
    return value;
}

std::vector<double> parse_feed(std::string_view text) {
    std::vector<double> feed;
    std::size_t begin = 0;
    while (begin <= text.size()) {
        const std::size_t comma = text.find(',', begin);
        const std::size_t end =
            comma == std::string_view::npos ? text.size() : comma;
        if (end == begin) {
            usage_error("--feed contains an empty component");
        }
        feed.push_back(parse_double(
            text.substr(begin, end - begin), "--feed component"
        ));
        if (comma == std::string_view::npos) {
            break;
        }
        begin = comma + 1;
    }
    if (feed.empty()) {
        usage_error("--feed requires at least one component");
    }
    return feed;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument(argv[index]);
        if (argument == "--trace") {
            options.trace = true;
            continue;
        }
        if (argument == "--help") {
            std::cout
                << "usage: epcsaft-equilibrium-diagnostic"
                   " --model-config PATH --temperature K --pressure PA"
                   " --feed X1,X2[,X3...] [--trace] [--output PATH]\n";
            std::exit(0);
        }
        if (index + 1 >= argc) {
            usage_error(std::string(argument) + " requires a value");
        }
        const std::string_view value(argv[++index]);
        if (argument == "--model-config") {
            options.model_config = value;
        } else if (argument == "--temperature") {
            options.temperature_k = parse_double(value, "--temperature");
        } else if (argument == "--pressure") {
            options.pressure_pa = parse_double(value, "--pressure");
        } else if (argument == "--feed") {
            options.feed = parse_feed(value);
        } else if (argument == "--output") {
            options.output_path = value;
        } else {
            usage_error("unknown option: " + std::string(argument));
        }
    }
    if (options.model_config.empty()) {
        usage_error("--model-config is required");
    }
    if (options.temperature_k <= 0.0) {
        usage_error("--temperature must be positive");
    }
    if (options.pressure_pa <= 0.0) {
        usage_error("--pressure must be positive");
    }
    if (options.feed.empty()) {
        usage_error("--feed is required");
    }
    return options;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        std::array<char, 1024> error{};
        using NativeModel = std::unique_ptr<
            epcsaft_native_model_handle_v1,
            decltype(&epcsaft_native_model_destroy_v1)
        >;
        NativeModel model(
            epcsaft_native_model_load_v1(
                options.model_config.c_str(), error.data(), error.size()
            ),
            &epcsaft_native_model_destroy_v1
        );
        if (!model) {
            throw std::runtime_error(
                error[0] == '\0'
                    ? "Provider native model load failed"
                    : std::string(error.data())
            );
        }
        const epcsaft_native_sdk_v1* sdk =
            epcsaft_native_model_sdk_v1(model.get());
        const char* fingerprint =
            epcsaft_native_model_fingerprint_v1(model.get());
        if (sdk == nullptr) {
            throw std::runtime_error(
                "Provider native model returned a null SDK table"
            );
        }
        if (fingerprint == nullptr || fingerprint[0] == '\0') {
            throw std::runtime_error(
                "Provider native model returned an empty fingerprint"
            );
        }

        const epcsaft_equilibrium::ProviderContext provider(
            *sdk, fingerprint
        );
        epcsaft_equilibrium::Held2TerminalProgress progress(std::cerr);
        const epcsaft_equilibrium::FlashResult result =
            epcsaft_equilibrium::solve_tp_flash(
                provider,
                {
                    options.temperature_k,
                    options.pressure_pa,
                    options.feed,
                },
                {},
                options.trace ? &progress : nullptr
            );
        const std::string json =
            epcsaft_equilibrium::flash_result_to_json(result);
        if (options.output_path.empty()) {
            std::cout << json;
        } else {
            std::ofstream output(
                options.output_path,
                std::ios::binary | std::ios::trunc
            );
            if (!output) {
                throw std::runtime_error(
                    "cannot open output path: " + options.output_path
                );
            }
            output << json;
            if (!output) {
                throw std::runtime_error(
                    "cannot write output path: " + options.output_path
                );
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "epcsaft-equilibrium-diagnostic: "
                  << error.what() << '\n';
        return 2;
    }
}
