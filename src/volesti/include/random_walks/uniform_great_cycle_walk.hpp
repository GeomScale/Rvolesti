// VolEsti (volume computation and sampling library)

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_GREAT_CYCLE_WALK_HPP
#define RANDOM_WALKS_UNIFORM_GREAT_CYCLE_WALK_HPP

#include <cmath>
#include <utility>

#include "random_walks/great_cycle_walk_helpers.hpp"

struct GCWalk
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

                std::pair<NT, NT> intersections =
                    P.gc_intersect(p, _v, _Ar, _Av, _lambda);

                _lambda =
                    rng.sample_urdist()
                        * (intersections.first - intersections.second)
                    + intersections.second;

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

        template <typename GenericPolytope>
        inline void initialize(GenericPolytope const& P,
                               Point& p,
                               RandomNumberGenerator& rng)
        {
            _Ar.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());

            _v = GetDirectionTangentPlane<Point>::apply(p, rng);

            std::pair<NT, NT> intersections =
                P.gc_intersect(p, _v, _Ar, _Av);

            _lambda =
                rng.sample_urdist()
                    * (intersections.first - intersections.second)
                + intersections.second;

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
