function samples = build_correlated_load_samples(mpc, count, seed, options)
%BUILD_CORRELATED_LOAD_SAMPLES Correlated bounded load P/Q fluctuations.
% Every bus load is scaled by 1+CV*sqrt(3)*(2*F_t(h)-1), so each multiplier
% is bounded with unit-variance standardized deviation. The latent h shares a
% system factor, an area factor and one chi-square scale per sample; reactive
% load has its own local factor correlated with the active one. Zero base
% loads stay zero. Same law as data/protocol/input_model.json of the cloud run.
arguments
    mpc (1,1) struct
    count (1,1) double {mustBeInteger, mustBePositive}
    seed (1,1) double {mustBeInteger, mustBeNonnegative}
    options.CoefficientOfVariation (1,1) double {mustBePositive} = 0.05
    options.LatentDf (1,1) double {mustBePositive} = 12
    options.SystemFactorWeight (1,1) double = 0.4
    options.AreaFactorWeight (1,1) double = 0.3
    options.ReactiveLocalCorrelation (1,1) double = 0.6
end
define_constants;
n = size(mpc.bus, 1);
localWeight = sqrt(1 - options.SystemFactorWeight^2 - options.AreaFactorWeight^2);
assert(isreal(localWeight) && localWeight > 0, "TMM_Transmission:scenarios:FactorWeights", ...
    "System and area factor weights must leave a positive local variance.");
stream = RandStream("mt19937ar", "Seed", seed);
systemFactor = randn(stream, count, 1);
areaFactor = randn(stream, count, 1);
localP = randn(stream, count, n);
localQ = randn(stream, count, n);
scale = sqrt(options.LatentDf ./ chi2rnd_stream(stream, options.LatentDf, count));
rho = options.ReactiveLocalCorrelation;
shared = options.SystemFactorWeight*systemFactor + options.AreaFactorWeight*areaFactor;
latentP = (shared + localWeight*localP).*scale;
latentQ = (shared + localWeight*(rho*localP + sqrt(1 - rho^2)*localQ)).*scale;
multiplierP = 1 + options.CoefficientOfVariation*sqrt(3)*(2*tcdf(latentP, options.LatentDf) - 1);
multiplierQ = 1 + options.CoefficientOfVariation*sqrt(3)*(2*tcdf(latentQ, options.LatentDf) - 1);
samples = struct("loadP", mpc.bus(:, PD).'.*multiplierP, "loadQ", mpc.bus(:, QD).'.*multiplierQ, ...
    "count", count, "seed", seed, "unit", "MW / MVAr", ...
    "coefficientOfVariation", options.CoefficientOfVariation, "latentDf", options.LatentDf);
end

function values = chi2rnd_stream(stream, df, count)
assert(mod(df, 1) == 0, "TMM_Transmission:scenarios:LatentDf", ...
    "The latent chi-square scale is drawn as a sum of squares and needs an integer df.");
values = sum(randn(stream, count, df).^2, 2);
end
