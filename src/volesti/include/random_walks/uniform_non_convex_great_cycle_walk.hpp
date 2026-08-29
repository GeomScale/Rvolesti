// VolEsti (volume computation and sampling library)

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_NONCONVEX_GREAT_CYCLE_WALK_HPP
#define RANDOM_WALKS_UNIFORM_NONCONVEX_GREAT_CYCLE_WALK_HPP

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "random_walks/great_cycle_walk_helpers.hpp"

struct GCWalkOpt
{
    struct parameters
    {
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

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P,
             Point& p,
             RandomNumberGenerator& rng)
        {
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P,
             Point& p,
             RandomNumberGenerator& rng,
             parameters const&)
        {
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        inline void apply(GenericPolytope const& P,
                          Point& p,
                          unsigned int const& walk_length,
                          RandomNumberGenerator& rng)
        {
            for (unsigned int step = 0; step < walk_length; ++step)
            {
                _v = GetDirectionTangentPlane<Point>::apply(p, rng);

                _roots =
                    P.gc_intersect_all_roots(
                        p,
                        _v,
                        _Ar,
                        _Av,
                        _lambda);

                _lambda =
                    pick_next_angle(
                        P,
                        p,
                        _v,
                        _roots.first,
                        _roots.second,
                        rng);

                if (!UpdateGreatCirclePoint<Point>::apply(
                        p,
                        _v,
                        _lambda,
                        _Ar,
                        _Av))
                {
                    // Reject a numerically invalid move while keeping the
                    // cached state consistent with the current point.
                    _lambda = NT(0);
                }
            }
        }

    private:
        NT _lambda;
        VT _Ar;
        VT _Av;
        Point _v;
        std::pair<VT, VT> _roots;

        template <typename GenericPolytope>
        NT pick_next_angle(GenericPolytope const& P,
                           Point const& p,
                           Point const& v,
                           VT const& negative_roots,
                           VT const& positive_roots,
                           RandomNumberGenerator& rng) const
        {
            std::vector<NT> roots;
            roots.reserve(
                static_cast<std::size_t>(
                    negative_roots.rows() + positive_roots.rows()));

            for (int i = 0; i < negative_roots.rows(); ++i)
            {
                roots.push_back(negative_roots(i));
            }

            for (int i = 0; i < positive_roots.rows(); ++i)
            {
                roots.push_back(positive_roots(i));
            }

            std::sort(roots.begin(), roots.end());

            std::vector<std::pair<NT, NT>> feasible_intervals;
            std::vector<NT> interval_lengths;

            NT const pi = std::acos(NT(-1));
            NT left_endpoint = -pi;
            NT total_length = NT(0);

            for (std::size_t i = 0; i <= roots.size(); ++i)
            {
                NT right_endpoint =
                    (i < roots.size()) ? roots[i] : pi;

                if (right_endpoint > left_endpoint)
                {
                    NT midpoint =
                        (left_endpoint + right_endpoint) / NT(2);

                    Point candidate =
                        std::cos(midpoint) * p
                        + std::sin(midpoint) * v;

                    if (P.is_in(candidate, NT(1e-10)) == -1)
                    {
                        NT length = right_endpoint - left_endpoint;

                        feasible_intervals.push_back(
                            std::make_pair(left_endpoint, right_endpoint));

                        interval_lengths.push_back(length);
                        total_length += length;
                    }
                }

                left_endpoint = right_endpoint;
            }

            if (feasible_intervals.empty())
            {
                return NT(0);
            }

            NT target = rng.sample_urdist() * total_length;
            NT accumulated_length = NT(0);

            for (std::size_t i = 0;
                 i < feasible_intervals.size();
                 ++i)
            {
                accumulated_length += interval_lengths[i];

                if (target <= accumulated_length)
                {
                    NT left = feasible_intervals[i].first;
                    NT right = feasible_intervals[i].second;

                    return left
                           + rng.sample_urdist() * (right - left);
                }
            }

            return feasible_intervals.back().second;
        }

        template <typename GenericPolytope>
        inline void initialize(GenericPolytope const& P,
                               Point& p,
                               RandomNumberGenerator& rng)
        {
            _Ar.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());

            _v = GetDirectionTangentPlane<Point>::apply(p, rng);

            _roots =
                P.gc_intersect_all_roots(
                    p,
                    _v,
                    _Ar,
                    _Av);

            _lambda =
                pick_next_angle(
                    P,
                    p,
                    _v,
                    _roots.first,
                    _roots.second,
                    rng);

            if (!UpdateGreatCirclePoint<Point>::apply(
                    p,
                    _v,
                    _lambda,
                    _Ar,
                    _Av))
            {
                _lambda = NT(0);
            }
        }
    };
};

#endif
