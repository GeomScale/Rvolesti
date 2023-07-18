// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// copula
Rcpp::NumericMatrix copula(Rcpp::Nullable<Rcpp::NumericVector> r1, Rcpp::Nullable<Rcpp::NumericVector> r2, Rcpp::Nullable<Rcpp::NumericMatrix> sigma, Rcpp::Nullable<unsigned int> m, Rcpp::Nullable<unsigned int> n, Rcpp::Nullable<double> seed);
RcppExport SEXP _volesti_copula(SEXP r1SEXP, SEXP r2SEXP, SEXP sigmaSEXP, SEXP mSEXP, SEXP nSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type r1(r1SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<unsigned int> >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<unsigned int> >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(copula(r1, r2, sigma, m, n, seed));
    return rcpp_result_gen;
END_RCPP
}
// direct_sampling
Rcpp::NumericMatrix direct_sampling(Rcpp::Nullable<Rcpp::List> body, Rcpp::Nullable<unsigned int> n, Rcpp::Nullable<double> seed);
RcppExport SEXP _volesti_direct_sampling(SEXP bodySEXP, SEXP nSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type body(bodySEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<unsigned int> >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(direct_sampling(body, n, seed));
    return rcpp_result_gen;
END_RCPP
}
// ess
Rcpp::NumericVector ess(Rcpp::NumericMatrix samples);
RcppExport SEXP _volesti_ess(SEXP samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type samples(samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(ess(samples));
    return rcpp_result_gen;
END_RCPP
}
// exact_vol
double exact_vol(Rcpp::Nullable<Rcpp::Reference> P);
RcppExport SEXP _volesti_exact_vol(SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::Reference> >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(exact_vol(P));
    return rcpp_result_gen;
END_RCPP
}
// frustum_of_simplex
double frustum_of_simplex(Rcpp::NumericVector a, double z0);
RcppExport SEXP _volesti_frustum_of_simplex(SEXP aSEXP, SEXP z0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type z0(z0SEXP);
    rcpp_result_gen = Rcpp::wrap(frustum_of_simplex(a, z0));
    return rcpp_result_gen;
END_RCPP
}
// geweke
bool geweke(Rcpp::NumericMatrix samples, Rcpp::Nullable<double> frac_first, Rcpp::Nullable<double> frac_last);
RcppExport SEXP _volesti_geweke(SEXP samplesSEXP, SEXP frac_firstSEXP, SEXP frac_lastSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type frac_first(frac_firstSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type frac_last(frac_lastSEXP);
    rcpp_result_gen = Rcpp::wrap(geweke(samples, frac_first, frac_last));
    return rcpp_result_gen;
END_RCPP
}
// inner_ball
Rcpp::NumericVector inner_ball(Rcpp::Reference P, Rcpp::Nullable<bool> lpsolve);
RcppExport SEXP _volesti_inner_ball(SEXP PSEXP, SEXP lpsolveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Reference >::type P(PSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type lpsolve(lpsolveSEXP);
    rcpp_result_gen = Rcpp::wrap(inner_ball(P, lpsolve));
    return rcpp_result_gen;
END_RCPP
}
// ode_solve
Rcpp::List ode_solve(Rcpp::Nullable<unsigned int> n, Rcpp::Nullable<double> step_size, Rcpp::Nullable<unsigned int> order, Rcpp::Nullable<unsigned int> dimension, Rcpp::Nullable<double> initial_time, Rcpp::Function F, Rcpp::Nullable<Rcpp::CharacterVector> method, Rcpp::Nullable<Rcpp::List> domains, Rcpp::Nullable<Rcpp::List> initial_conditions);
RcppExport SEXP _volesti_ode_solve(SEXP nSEXP, SEXP step_sizeSEXP, SEXP orderSEXP, SEXP dimensionSEXP, SEXP initial_timeSEXP, SEXP FSEXP, SEXP methodSEXP, SEXP domainsSEXP, SEXP initial_conditionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<unsigned int> >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<unsigned int> >::type order(orderSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<unsigned int> >::type dimension(dimensionSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type initial_time(initial_timeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type F(FSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::CharacterVector> >::type method(methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type domains(domainsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initial_conditions(initial_conditionsSEXP);
    rcpp_result_gen = Rcpp::wrap(ode_solve(n, step_size, order, dimension, initial_time, F, method, domains, initial_conditions));
    return rcpp_result_gen;
END_RCPP
}
// poly_gen
Rcpp::NumericMatrix poly_gen(int kind_gen, bool Vpoly_gen, bool Zono_gen, int dim_gen, int m_gen, Rcpp::Nullable<double> seed);
RcppExport SEXP _volesti_poly_gen(SEXP kind_genSEXP, SEXP Vpoly_genSEXP, SEXP Zono_genSEXP, SEXP dim_genSEXP, SEXP m_genSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type kind_gen(kind_genSEXP);
    Rcpp::traits::input_parameter< bool >::type Vpoly_gen(Vpoly_genSEXP);
    Rcpp::traits::input_parameter< bool >::type Zono_gen(Zono_genSEXP);
    Rcpp::traits::input_parameter< int >::type dim_gen(dim_genSEXP);
    Rcpp::traits::input_parameter< int >::type m_gen(m_genSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(poly_gen(kind_gen, Vpoly_gen, Zono_gen, dim_gen, m_gen, seed));
    return rcpp_result_gen;
END_RCPP
}
// psrf_multivariate
double psrf_multivariate(Rcpp::NumericMatrix samples);
RcppExport SEXP _volesti_psrf_multivariate(SEXP samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type samples(samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(psrf_multivariate(samples));
    return rcpp_result_gen;
END_RCPP
}
// psrf_univariate
Rcpp::NumericVector psrf_univariate(Rcpp::NumericMatrix samples, Rcpp::Nullable<std::string> method);
RcppExport SEXP _volesti_psrf_univariate(SEXP samplesSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(psrf_univariate(samples, method));
    return rcpp_result_gen;
END_RCPP
}
// raftery
Rcpp::NumericMatrix raftery(Rcpp::NumericMatrix samples, Rcpp::Nullable<double> q, Rcpp::Nullable<double> r, Rcpp::Nullable<double> s);
RcppExport SEXP _volesti_raftery(SEXP samplesSEXP, SEXP qSEXP, SEXP rSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type q(qSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type r(rSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(raftery(samples, q, r, s));
    return rcpp_result_gen;
END_RCPP
}
// rotating
Rcpp::NumericMatrix rotating(Rcpp::Reference P, Rcpp::Nullable<Rcpp::NumericMatrix> T, Rcpp::Nullable<int> seed);
RcppExport SEXP _volesti_rotating(SEXP PSEXP, SEXP TSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Reference >::type P(PSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type T(TSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(rotating(P, T, seed));
    return rcpp_result_gen;
END_RCPP
}
// rounding
Rcpp::List rounding(Rcpp::Reference P, Rcpp::Nullable<std::string> method, Rcpp::Nullable<double> seed);
RcppExport SEXP _volesti_rounding(SEXP PSEXP, SEXP methodSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Reference >::type P(PSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type method(methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(rounding(P, method, seed));
    return rcpp_result_gen;
END_RCPP
}
// sample_points
Rcpp::NumericMatrix sample_points(Rcpp::Nullable<Rcpp::Reference> P, Rcpp::Nullable<unsigned int> n, Rcpp::Nullable<Rcpp::List> random_walk, Rcpp::Nullable<Rcpp::List> distribution, Rcpp::Nullable<double> seed);
RcppExport SEXP _volesti_sample_points(SEXP PSEXP, SEXP nSEXP, SEXP random_walkSEXP, SEXP distributionSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::Reference> >::type P(PSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<unsigned int> >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type random_walk(random_walkSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_points(P, n, random_walk, distribution, seed));
    return rcpp_result_gen;
END_RCPP
}
// writeSdpaFormatFile
void writeSdpaFormatFile(Rcpp::Nullable<Rcpp::Reference> spectrahedron, Rcpp::Nullable<Rcpp::NumericVector> objectiveFunction, Rcpp::Nullable<std::string> outputFile);
RcppExport SEXP _volesti_writeSdpaFormatFile(SEXP spectrahedronSEXP, SEXP objectiveFunctionSEXP, SEXP outputFileSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::Reference> >::type spectrahedron(spectrahedronSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type objectiveFunction(objectiveFunctionSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type outputFile(outputFileSEXP);
    writeSdpaFormatFile(spectrahedron, objectiveFunction, outputFile);
    return R_NilValue;
END_RCPP
}
// loadSdpaFormatFile
Rcpp::List loadSdpaFormatFile(Rcpp::Nullable<std::string> inputFile);
RcppExport SEXP _volesti_loadSdpaFormatFile(SEXP inputFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type inputFile(inputFileSEXP);
    rcpp_result_gen = Rcpp::wrap(loadSdpaFormatFile(inputFile));
    return rcpp_result_gen;
END_RCPP
}
// volume
Rcpp::List volume(Rcpp::Reference P, Rcpp::Nullable<Rcpp::List> settings, Rcpp::Nullable<std::string> rounding, Rcpp::Nullable<double> seed);
RcppExport SEXP _volesti_volume(SEXP PSEXP, SEXP settingsSEXP, SEXP roundingSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Reference >::type P(PSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type settings(settingsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type rounding(roundingSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(volume(P, settings, rounding, seed));
    return rcpp_result_gen;
END_RCPP
}
// zono_approx
Rcpp::List zono_approx(Rcpp::Reference Z, Rcpp::Nullable<bool> fit_ratio, Rcpp::Nullable<Rcpp::List> settings, Rcpp::Nullable<double> seed);
RcppExport SEXP _volesti_zono_approx(SEXP ZSEXP, SEXP fit_ratioSEXP, SEXP settingsSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Reference >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type fit_ratio(fit_ratioSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type settings(settingsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(zono_approx(Z, fit_ratio, settings, seed));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_polytopes();
RcppExport SEXP _rcpp_module_boot_spectrahedron();

static const R_CallMethodDef CallEntries[] = {
    {"_volesti_copula", (DL_FUNC) &_volesti_copula, 6},
    {"_volesti_direct_sampling", (DL_FUNC) &_volesti_direct_sampling, 3},
    {"_volesti_ess", (DL_FUNC) &_volesti_ess, 1},
    {"_volesti_exact_vol", (DL_FUNC) &_volesti_exact_vol, 1},
    {"_volesti_frustum_of_simplex", (DL_FUNC) &_volesti_frustum_of_simplex, 2},
    {"_volesti_geweke", (DL_FUNC) &_volesti_geweke, 3},
    {"_volesti_inner_ball", (DL_FUNC) &_volesti_inner_ball, 2},
    {"_volesti_ode_solve", (DL_FUNC) &_volesti_ode_solve, 9},
    {"_volesti_poly_gen", (DL_FUNC) &_volesti_poly_gen, 6},
    {"_volesti_psrf_multivariate", (DL_FUNC) &_volesti_psrf_multivariate, 1},
    {"_volesti_psrf_univariate", (DL_FUNC) &_volesti_psrf_univariate, 2},
    {"_volesti_raftery", (DL_FUNC) &_volesti_raftery, 4},
    {"_volesti_rotating", (DL_FUNC) &_volesti_rotating, 3},
    {"_volesti_rounding", (DL_FUNC) &_volesti_rounding, 3},
    {"_volesti_sample_points", (DL_FUNC) &_volesti_sample_points, 5},
    {"_volesti_writeSdpaFormatFile", (DL_FUNC) &_volesti_writeSdpaFormatFile, 3},
    {"_volesti_loadSdpaFormatFile", (DL_FUNC) &_volesti_loadSdpaFormatFile, 1},
    {"_volesti_volume", (DL_FUNC) &_volesti_volume, 4},
    {"_volesti_zono_approx", (DL_FUNC) &_volesti_zono_approx, 4},
    {"_rcpp_module_boot_polytopes", (DL_FUNC) &_rcpp_module_boot_polytopes, 0},
    {"_rcpp_module_boot_spectrahedron", (DL_FUNC) &_rcpp_module_boot_spectrahedron, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_volesti(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}