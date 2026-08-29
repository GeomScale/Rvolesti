// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Maria Zaza

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DIAGNOSTICS_MULTICHAIN_DIAGNOSTICS_HPP
#define DIAGNOSTICS_MULTICHAIN_DIAGNOSTICS_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include <Eigen/Eigen>
#include <unsupported/Eigen/FFT>

/// Potential scale reduction factors for multiple independent chains.
template <typename NT, typename VT>
struct MultichainPSRF
{
    VT marginal;
    NT multivariate;
};

/// Per-chain and summed effective sample sizes for multiple chains.
template <typename VT>
struct MultichainESS
{
    std::vector<VT> per_chain;
    VT total;
};

/// Complete PSRF/ESS diagnostic summary for multiple independent chains.
template <typename NT, typename VT>
struct MultichainDiagnostics
{
    VT marginal_psrf;
    NT multivariate_psrf;
    std::vector<VT> per_chain_ess;
    VT effective_sample_size;
    bool converged;
};

namespace multichain_diagnostics_detail
{

template <typename MT>
void validate_chains(
    std::vector<MT> const& chains,
    Eigen::Index minimum_draws)
{
    if (chains.size() < 2)
    {
        throw std::invalid_argument(
            "multichain diagnostics require at least two chains");
    }

    Eigen::Index const dimension = chains.front().rows();
    Eigen::Index const draws = chains.front().cols();

    if (dimension <= 0)
    {
        throw std::invalid_argument(
            "multichain diagnostics require a positive dimension");
    }

    if (draws < minimum_draws)
    {
        throw std::invalid_argument(
            "multichain diagnostics received too few draws per chain");
    }

    for (MT const& chain : chains)
    {
        if (chain.rows() != dimension || chain.cols() != draws)
        {
            throw std::invalid_argument(
                "multichain diagnostics require equal chain dimensions");
        }

        if (!chain.allFinite())
        {
            throw std::invalid_argument(
                "multichain diagnostics require finite samples");
        }
    }
}

template <typename NT, typename VT, typename MT>
void covariance_components(
    std::vector<MT> const& chains,
    MT& within,
    MT& between)
{
    Eigen::Index const dimension = chains.front().rows();
    Eigen::Index const draws = chains.front().cols();
    std::size_t const number_of_chains = chains.size();

    within = MT::Zero(dimension, dimension);
    between = MT::Zero(dimension, dimension);

    std::vector<VT> chain_means;
    chain_means.reserve(number_of_chains);

    VT pooled_mean = VT::Zero(dimension);

    for (MT const& chain : chains)
    {
        VT const mean = chain.rowwise().mean();
        chain_means.push_back(mean);
        pooled_mean += mean;

        MT const centered = chain.colwise() - mean;
        within.noalias() +=
            centered * centered.transpose() /
            static_cast<NT>(draws - 1);
    }

    within /= static_cast<NT>(number_of_chains);
    pooled_mean /= static_cast<NT>(number_of_chains);

    for (VT const& mean : chain_means)
    {
        VT const difference = mean - pooled_mean;
        between.noalias() += difference * difference.transpose();
    }

    between *=
        static_cast<NT>(draws) /
        static_cast<NT>(number_of_chains - 1);
}

template <typename NT, typename MT>
NT multivariate_psrf_from_covariances(
    MT const& within,
    MT const& between,
    std::size_t const number_of_chains,
    Eigen::Index const draws)
{
    Eigen::SelfAdjointEigenSolver<MT> within_solver(within);

    if (within_solver.info() != Eigen::Success)
    {
        throw std::runtime_error(
            "multichain PSRF failed to decompose within-chain covariance");
    }

    auto const& eigenvalues = within_solver.eigenvalues();
    MT const& eigenvectors = within_solver.eigenvectors();

    NT const largest_within =
        std::max(NT(0), eigenvalues.maxCoeff());

    NT const tolerance =
        NT(128) * std::numeric_limits<NT>::epsilon() *
        largest_within;

    if (eigenvalues.minCoeff() < -tolerance)
    {
        throw std::runtime_error(
            "multichain PSRF found an invalid within-chain covariance");
    }

    std::vector<Eigen::Index> positive_indices;
    std::vector<Eigen::Index> null_indices;

    for (Eigen::Index index = 0; index < eigenvalues.rows(); ++index)
    {
        if (eigenvalues(index) > tolerance)
        {
            positive_indices.push_back(index);
        }
        else
        {
            null_indices.push_back(index);
        }
    }

    if (!null_indices.empty())
    {
        // W^{-1} is undefined for a singular within-chain covariance.
        // Reporting infinity prevents a stuck or under-dispersed collection
        // of chains from being accepted as converged.
        return std::numeric_limits<NT>::infinity();
    }

    MT positive_basis(within.rows(), positive_indices.size());
    MT inverse_sqrt =
        MT::Zero(positive_indices.size(), positive_indices.size());

    for (std::size_t column = 0;
         column < positive_indices.size();
         ++column)
    {
        Eigen::Index const index = positive_indices[column];
        positive_basis.col(static_cast<Eigen::Index>(column)) =
            eigenvectors.col(index);

        inverse_sqrt(
            static_cast<Eigen::Index>(column),
            static_cast<Eigen::Index>(column)) =
                NT(1) / std::sqrt(eigenvalues(index));
    }

    MT whitened_between =
        inverse_sqrt * positive_basis.transpose() *
        between * positive_basis * inverse_sqrt;

    // Remove round-off asymmetry before using a self-adjoint solver.
    whitened_between =
        (whitened_between + whitened_between.transpose()) / NT(2);

    Eigen::SelfAdjointEigenSolver<MT> between_solver(
        whitened_between);

    if (between_solver.info() != Eigen::Success)
    {
        throw std::runtime_error(
            "multichain PSRF failed to solve the generalized eigenproblem");
    }

    NT const largest_generalized_eigenvalue =
        std::max(NT(0), between_solver.eigenvalues().maxCoeff());

    NT const variance_ratio =
        static_cast<NT>(draws - 1) / static_cast<NT>(draws) +
        (static_cast<NT>(number_of_chains + 1) /
         (static_cast<NT>(number_of_chains) *
          static_cast<NT>(draws))) *
            largest_generalized_eigenvalue;

    return std::sqrt(std::max(NT(0), variance_ratio));
}

/// FFT/Geyer effective sample size for one scalar chain.
///
/// The autocorrelations use the biased (divide by n) FFT estimator recommended
/// by Geyer and used by Stan. Adjacent autocorrelations are truncated with the
/// initial-positive sequence and converted to a monotone sequence. A zero
/// return value marks a constant chain, for which ESS is not defined.
template <typename NT, typename Derived>
NT single_chain_effective_sample_size(
    Eigen::MatrixBase<Derived> const& samples)
{
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<std::complex<NT>, Eigen::Dynamic, 1> CVT;

    Eigen::Index const draws = samples.size();

    VT centered = samples;
    centered.array() -= samples.mean();

    NT const centered_scale = centered.cwiseAbs().maxCoeff();

    if (!std::isfinite(centered_scale) || centered_scale <= NT(0))
    {
        return NT(0);
    }

    // Normalizing before the FFT avoids over/underflow without changing the
    // autocorrelation, so ESS remains invariant under nonzero rescaling.
    centered /= centered_scale;

    Eigen::Index const fft_size = 2 * draws;
    VT padded = VT::Zero(fft_size);
    padded.head(draws) = centered;

    Eigen::FFT<NT> fft;
    CVT frequency(fft_size);
    CVT autocovariance(fft_size);

    fft.fwd(frequency, padded);
    frequency = frequency.cwiseAbs2();
    fft.inv(autocovariance, frequency);

    NT const lag_zero = autocovariance(0).real();

    if (!std::isfinite(lag_zero) || lag_zero <= NT(0))
    {
        return NT(0);
    }

    std::vector<NT> pair_sums;
    pair_sums.reserve(static_cast<std::size_t>(draws / 2));

    for (Eigen::Index lag = 0; lag + 1 < draws; lag += 2)
    {
        NT even = autocovariance(lag).real() / lag_zero;
        NT odd = autocovariance(lag + 1).real() / lag_zero;

        // Suppress harmless FFT round-off outside the correlation range.
        even = std::max(NT(-1), std::min(NT(1), even));
        odd = std::max(NT(-1), std::min(NT(1), odd));

        NT pair_sum = even + odd;

        if (!std::isfinite(pair_sum) || pair_sum <= NT(0))
        {
            break;
        }

        if (!pair_sums.empty())
        {
            pair_sum = std::min(pair_sum, pair_sums.back());
        }

        pair_sums.push_back(pair_sum);
    }

    NT integrated_autocorrelation = NT(-1);

    for (NT const pair_sum : pair_sums)
    {
        integrated_autocorrelation += NT(2) * pair_sum;
    }

    // Stan permits ESS above n for antithetic chains, while bounding extreme
    // finite-sample estimates by n*log10(n).
    NT const log_draws = std::log10(static_cast<NT>(draws));
    NT const minimum_autocorrelation = NT(1) / log_draws;

    integrated_autocorrelation =
        std::max(integrated_autocorrelation, minimum_autocorrelation);

    return static_cast<NT>(draws) / integrated_autocorrelation;
}

} // namespace multichain_diagnostics_detail

/// Compute true multichain marginal and multivariate Gelman--Rubin PSRF.
///
/// Each input matrix contains one chain column-wise (dimension x draws).
/// All chains must be independent, target the same distribution, and have
/// equal post-burn-in lengths. The multivariate value is the standard
/// Brooks--Gelman square root of the largest generalized eigenvalue.
template <typename NT, typename VT, typename MT>
MultichainPSRF<NT, VT> compute_multichain_psrf(
    std::vector<MT> const& chains)
{
    multichain_diagnostics_detail::validate_chains(chains, 2);

    Eigen::Index const dimension = chains.front().rows();
    Eigen::Index const draws = chains.front().cols();
    std::size_t const number_of_chains = chains.size();

    MT within;
    MT between;

    multichain_diagnostics_detail::covariance_components<NT, VT>(
        chains,
        within,
        between);

    VT marginal(dimension);

    NT const within_weight =
        static_cast<NT>(draws - 1) / static_cast<NT>(draws);

    NT const between_weight =
        static_cast<NT>(number_of_chains + 1) /
        (static_cast<NT>(number_of_chains) *
         static_cast<NT>(draws));

    for (Eigen::Index coordinate = 0;
         coordinate < dimension;
         ++coordinate)
    {
        NT const within_variance = within(coordinate, coordinate);
        NT const between_variance = between(coordinate, coordinate);

        NT const scale = std::max(
            std::abs(within_variance),
            std::abs(between_variance));

        NT const tolerance =
            NT(128) * std::numeric_limits<NT>::epsilon() * scale;

        if (within_variance <= tolerance)
        {
            marginal(coordinate) =
                std::numeric_limits<NT>::infinity();

            continue;
        }

        NT const variance_ratio =
            within_weight +
            between_weight * between_variance / within_variance;

        marginal(coordinate) =
            std::sqrt(std::max(NT(0), variance_ratio));
    }

    NT const multivariate =
        multichain_diagnostics_detail::
            multivariate_psrf_from_covariances<NT>(
                within,
                between,
                number_of_chains,
                draws);

    return MultichainPSRF<NT, VT>{marginal, multivariate};
}

/// Compute Stan-style FFT/Geyer ESS for every chain and sum corresponding
/// coordinates. A constant coordinate receives ESS zero, making the complete
/// diagnostic non-converged rather than crediting a stuck chain with n draws.
template <typename NT, typename VT, typename MT>
MultichainESS<VT> compute_multichain_effective_sample_size(
    std::vector<MT> const& chains)
{
    multichain_diagnostics_detail::validate_chains(chains, 4);

    Eigen::Index const dimension = chains.front().rows();
    MultichainESS<VT> result;
    result.total = VT::Zero(dimension);
    result.per_chain.reserve(chains.size());

    for (MT const& chain : chains)
    {
        VT chain_ess(dimension);

        for (Eigen::Index coordinate = 0;
             coordinate < dimension;
             ++coordinate)
        {
            chain_ess(coordinate) =
                multichain_diagnostics_detail::
                    single_chain_effective_sample_size<NT>(
                        chain.row(coordinate).transpose());

            if (!std::isfinite(chain_ess(coordinate)) ||
                chain_ess(coordinate) < NT(0))
            {
                throw std::runtime_error(
                    "multichain ESS produced an invalid estimate");
            }
        }

        result.total += chain_ess;
        result.per_chain.push_back(chain_ess);
    }

    return result;
}

/// Evaluate the Task 4 stopping rule for multiple independent chains.
template <typename NT, typename VT, typename MT>
MultichainDiagnostics<NT, VT> evaluate_multichain_diagnostics(
    std::vector<MT> const& chains,
    NT const psrf_target = NT(1.1),
    NT const ess_target = NT(100))
{
    if (!std::isfinite(psrf_target) || psrf_target <= NT(1))
    {
        throw std::invalid_argument(
            "multichain diagnostics require a finite PSRF target above one");
    }

    if (!std::isfinite(ess_target) || ess_target <= NT(0))
    {
        throw std::invalid_argument(
            "multichain diagnostics require a positive finite ESS target");
    }

    MultichainPSRF<NT, VT> const psrf =
        compute_multichain_psrf<NT, VT>(chains);

    MultichainESS<VT> const ess =
        compute_multichain_effective_sample_size<NT, VT>(chains);

    bool every_chain_has_variation = true;

    for (VT const& chain_ess : ess.per_chain)
    {
        every_chain_has_variation =
            every_chain_has_variation &&
            chain_ess.minCoeff() > NT(0);
    }

    bool const converged =
        every_chain_has_variation &&
        psrf.marginal.maxCoeff() < psrf_target &&
        std::isfinite(psrf.multivariate) &&
        psrf.multivariate < psrf_target &&
        ess.total.minCoeff() >= ess_target;

    return MultichainDiagnostics<NT, VT>{
        psrf.marginal,
        psrf.multivariate,
        ess.per_chain,
        ess.total,
        converged};
}

#endif // DIAGNOSTICS_MULTICHAIN_DIAGNOSTICS_HPP
