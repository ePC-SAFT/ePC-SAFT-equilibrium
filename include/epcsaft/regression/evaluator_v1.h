#ifndef EPCSAFT_REGRESSION_EVALUATOR_V1_H
#define EPCSAFT_REGRESSION_EVALUATOR_V1_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EPCSAFT_REGRESSION_EVALUATOR_V1_CAPSULE_NAME "epcsaft.regression.evaluator.v1"
#define EPCSAFT_REGRESSION_EVALUATOR_V1_ABI_VERSION 1u
#define EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY 256u

enum {
    EPCSAFT_REGRESSION_EVALUATOR_STATUS_OK_V1 = 0,
    EPCSAFT_REGRESSION_EVALUATOR_STATUS_UNAVAILABLE_V1 = 1,
    EPCSAFT_REGRESSION_EVALUATOR_STATUS_INVALID_INPUT_V1 = 2,
};
enum {
    EPCSAFT_REGRESSION_EVALUATOR_REQUEST_VALUES_ONLY_V1 = 0,
    EPCSAFT_REGRESSION_EVALUATOR_REQUEST_VALUES_AND_JACOBIAN_V1 = 1,
};

typedef struct epcsaft_regression_evaluator_row_result_v1 {
    uint32_t struct_size;
    int32_t status;
    char reason[EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY];
    char solver_status[EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY];
    char numerical_status[EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY];
    char physical_status[EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY];
    char derivative_status[EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY];
    char chart_topology[EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY];
    char provider_topology_fingerprint[
        EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY
    ];
    size_t kkt_dimension;
    size_t kkt_rank;
    double kkt_condition_number_inf;
} epcsaft_regression_evaluator_row_result_v1;

typedef struct epcsaft_regression_evaluator_result_v1 {
    uint32_t struct_size;
    int32_t status;
    size_t row_count;
    size_t parameter_count;
    size_t value_capacity;
    size_t jacobian_capacity;
    size_t row_result_capacity;
    int32_t request_mode;
    double* values;
    double* jacobian;
    epcsaft_regression_evaluator_row_result_v1* row_results;
    char provider_parameter_fingerprint[EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY];
    char artifact_identity[EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY];
    char error[EPCSAFT_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY];
} epcsaft_regression_evaluator_result_v1;

typedef int (*epcsaft_regression_evaluator_evaluate_v1)(
    void* context,
    const double* parameter_values,
    size_t parameter_count,
    int32_t request_mode,
    epcsaft_regression_evaluator_result_v1* result
);

typedef struct epcsaft_regression_evaluator_sdk_v1 {
    uint32_t abi_version;
    size_t table_size;
    void* model_context;
    size_t row_count;
    size_t parameter_count;
    const char* const* row_ids;
    const char* const* state_ids;
    const char* const* state_schema_ids;
    const char* const* observation_source_ids;
    const char* const* primitive_ids;
    const char* const* primitive_units;
    const char* const* transform_ids;
    const char* const* reference_ids;
    const char* const* reference_fingerprints;
    const char* const* parameter_ids;
    const char* const* parameter_units;
    const char* evaluator_identity;
    const char* capability_id;
    const char* capability_fingerprint;
    const char* provider_artifact_identity;
    const char* owner_artifact_identity;
    const char* contract_fingerprint;
    const char* model_fingerprint;
    const char* provider_parameter_fingerprint;
    const char* expected_provider_topology_fingerprint;
    const char* provider_sdk_capsule_name;
    uint32_t provider_sdk_abi_version;
    size_t provider_sdk_table_size;
    size_t provider_sdk_result_size;
    size_t provider_sdk_mixture_result_size;
    size_t provider_sdk_neutral_reference_result_size;
    size_t provider_sdk_neutral_reference_derivative_result_size;
    size_t provider_sdk_reacting_phase_parameter_result_size;
    const char* artifact_identity;
    int32_t single_thread_non_reentrant;
    int32_t value_only_avoids_derivative_work;
    size_t result_size;
    size_t row_result_size;
    epcsaft_regression_evaluator_evaluate_v1 evaluate;
} epcsaft_regression_evaluator_sdk_v1;

/*
 * The capsule owner retains the SDK table, context, and upstream model for the
 * complete fit. Calls are synchronous; the caller owns all result buffers.
 * The callback is native and does not require the Python GIL. A v1 consumer
 * must provide exact result and row-result struct sizes. The table's declared
 * single-thread/non-reentrant mode is authoritative.
 */

#ifdef __cplusplus
}
#endif

#endif
