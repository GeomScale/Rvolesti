// VolEsti (volume computation and sampling library)

// Contributed and/or modified by Apostolos Chalkis,
// as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_BILLIARD_GCW_WALK_HPP
#define RANDOM_WALKS_UNIFORM_BILLIARD_GCW_WALK_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include <Eigen/Eigen>

#include "random_walks/great_cycle_walk_helpers.hpp"


// Reflective Great Cycle Walk for the uniform distribution.
struct BilliardGCWalk
{
    BilliardGCWalk(double tau, unsigned int max_reflections)
        : param(tau, true, max_reflections, true)
    {
    }

    BilliardGCWalk(double tau)
        : param(tau, true, 0, false)
    {
    }

    BilliardGCWalk()
        : param(1.0, false, 0, false)
    {
    }

    struct parameters
    {
        parameters(double tau,
                   bool set_tau_,
                   unsigned int max_reflections_,
                   bool set_max_reflections_)
            : m_tau(tau),
              set_tau(set_tau_),
              m_max_reflections(max_reflections_),
              set_max_reflections(set_max_reflections_)
        {
        }

        double m_tau;
        bool set_tau;
        unsigned int m_max_reflections;
        bool set_max_reflections;
    };

    parameters param;


    template <
        typename Polytope,
        typename RandomNumberGenerator>
    struct Walk
    {
        typedef typename Polytope::PointType Point;
        typedef typename Point::FT NT;
        typedef typename Point::Coeff VT;
        typedef typename Polytope::MT MT;

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P,
             Point& p,
             NT const& tau,
             RandomNumberGenerator& rng)
            : _tau(tau),
              _max_reflections(500000u * P.dimension())
        {
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P,
             Point& p,
             NT const& tau,
             unsigned int max_reflections,
             RandomNumberGenerator& rng)
            : _tau(tau),
              _max_reflections(max_reflections)
        {
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P,
             Point& p,
             RandomNumberGenerator& rng,
             parameters const& params)
            : _tau(params.set_tau ? NT(params.m_tau) : NT(1)),
              _max_reflections(
                  params.set_max_reflections
                      ? params.m_max_reflections
                      : 500000u * P.dimension())
        {
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        inline void apply(GenericPolytope const& P,
                          Point& p,
                          unsigned int const& walk_length,
                          RandomNumberGenerator& rng)
        {
            for (unsigned int i = 0; i < walk_length; ++i)
            {
                advance(P, p, rng);
            }
        }

        inline void update_delta(NT tau)
        {
            _tau = tau;
        }

        inline void update_max_reflections(
            unsigned int max_reflections)
        {
            _max_reflections = max_reflections;
        }

    private:

        inline static bool normalize_on_sphere(
            VT& coefficients)
        {
            NT const norm =
                coefficients.norm();

            if (!std::isfinite(norm) ||
                norm <=
                    std::numeric_limits<NT>::epsilon())
            {
                return false;
            }

            coefficients /= norm;
            return true;
        }

        template <typename GenericPolytope>
        inline void initialize(GenericPolytope const& P,
                               Point& p,
                               RandomNumberGenerator& rng)
        {
            unsigned int const dimension = P.dimension();

            _Ar.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());

            _identity =
                MT::Identity(dimension, dimension);

            _projector.setZero(dimension, dimension);

            // Preserve the constructor behaviour of the branch walks:
            // constructing the walk performs the first trajectory.
            advance(P, p, rng);
        }

        template <typename GenericPolytope>
        inline void advance(GenericPolytope const& P,
                            Point& p,
                            RandomNumberGenerator& rng)
        {
            Point const initial_point = p;

            // Exponentially distributed trajectory length:
            //     L = -tau * log(eta), eta ~ Uniform(0, 1).
            NT eta = rng.sample_urdist();

            eta = std::max(
                eta,
                std::numeric_limits<NT>::min());

            NT remaining_length =
                -_tau * std::log(eta);

            _v =
                GetDirectionTangentPlane<Point>::apply(p, rng);

            unsigned int reflections = 0;

            while (true)
            {
                std::pair<NT, int> const boundary =
                    P.gc_intersect_positive(
                        p,
                        _v,
                        _Ar,
                        _Av);

                // The trajectory finishes before reaching a facet,
                // or no valid facet lies in the positive direction.
                if (remaining_length <= boundary.first ||
                    boundary.second < 0)
                {
                    VT const p_coeffs =
                        p.getCoefficients();

                    VT const v_coeffs =
                        _v.getCoefficients();

                    VT next =
                        std::cos(remaining_length) * p_coeffs +
                        std::sin(remaining_length) * v_coeffs;

                    if (!normalize_on_sphere(next))
                    {
                        p = initial_point;
                        return;
                    }

                    p = Point(next);
                    return;
                }

                NT const angle_tolerance =
                    NT(128) *
                    std::numeric_limits<NT>::epsilon();

                if (!std::isfinite(boundary.first) ||
                    boundary.first < -angle_tolerance)
                {
                    p = initial_point;
                    return;
                }

                // The requested trajectory requires more reflections
                // than permitted. Return the starting point.
                if (reflections >= _max_reflections)
                {
                    p = initial_point;
                    return;
                }

                // Move exactly to the facet before computing the reflection.
                NT const hit_angle =
                    std::max(NT(0), boundary.first);

                VT const p_coeffs =
                    p.getCoefficients();

                VT const v_coeffs =
                    _v.getCoefficients();

                // Tangent direction of the same great circle at the facet,
                // immediately before reflection.
                VT reflected_direction =
                    -std::sin(hit_angle) * p_coeffs +
                    std::cos(hit_angle) * v_coeffs;

                VT next =
                    std::cos(hit_angle) * p_coeffs +
                    std::sin(hit_angle) * v_coeffs;

                if (!normalize_on_sphere(next))
                {
                    p = initial_point;
                    return;
                }

                p = Point(next);
                remaining_length -= hit_angle;

                _projector =
                    _identity -
                    next * next.transpose();

                P.compute_reflection(
                    reflected_direction,
                    next,
                    _projector,
                    boundary.second);

                // Remove numerical radial drift and normalise.
                reflected_direction =
                    (_projector * reflected_direction).eval();

                NT const direction_norm =
                    reflected_direction.norm();

                if (!std::isfinite(direction_norm) ||
                    direction_norm <=
                        std::numeric_limits<NT>::epsilon())
                {
                    p = initial_point;
                    return;
                }

                reflected_direction /= direction_norm;
                _v = Point(reflected_direction);

                ++reflections;
            }
        }

        NT _tau;
        unsigned int _max_reflections;

        Point _v;

        MT _identity;
        MT _projector;

        VT _Ar;
        VT _Av;
    };
};

#endif
