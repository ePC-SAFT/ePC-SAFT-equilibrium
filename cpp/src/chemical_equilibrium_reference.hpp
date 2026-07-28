#pragma once

#include "chemical_equilibrium.hpp"

#include <cstddef>
#include <vector>

namespace epcsaft_equilibrium {

struct ReactionReferenceReconstruction {
    std::vector<double> reference;
    std::size_t rank = 0;
    double qr_diagonal_ratio = 0.0;
    double residual_inf_norm = 0.0;
};

[[nodiscard]] ReactionReferenceReconstruction reconstruct_reaction_reference(
    const DenseMatrix& reaction_matrix,
    const std::vector<double>& ln_k
);

}  // namespace epcsaft_equilibrium
