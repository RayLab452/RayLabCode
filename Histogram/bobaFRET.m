function results = bobaFRET(perMolHist, centers, startGuess, nBoot, seed, nStartsFull, nStartsBoot)
%BOBAFRET  Bootstrap-based analysis (BOBA) of a molecule-weighted FRET histogram.
%
%   results = bobaFRET(perMolHist, centers, startGuess)
%   results = bobaFRET(perMolHist, centers, startGuess, nBoot, seed)
%   results = bobaFRET(perMolHist, centers, startGuess, nBoot, seed, nStartsFull, nStartsBoot)
%
%   Inputs
%     perMolHist  nMolecules x nBins. Row i is molecule i's histogram, normalized
%                 exactly as in SM_MakeHistogram_new (counts ./ number of frames).
%     centers     1 x nBins bin centers (histogramData(1,:) in the GUI).
%     startGuess  [mu1 sigma1 mu2 sigma2 ...]  starting values for the Gaussian fit.
%                 sigma is the Gaussian's standard deviation, the only width unit
%                 used anywhere in this file -- there is no FWHM conversion in the
%                 fitting math at all (FWHM = 2*sqrt(2*ln2)*sigma, if you ever need
%                 to report it, but nothing here computes or relies on it).
%     nBoot       number of bootstrap replicates (default 1000).
%     seed        optional integer for rng(), for reproducible results.
%     nStartsFull number of randomly-jittered restarts tried for the FULL-DATA fit
%                 (default 8). fminsearch only walks downhill from wherever it
%                 starts, so one bad startGuess can lock the whole analysis into
%                 the wrong solution (e.g. two states collapsing together). Trying
%                 several jittered starting points and keeping the best one makes
%                 the full fit much less sensitive to exactly where startGuess sits.
%     nStartsBoot number of jittered restarts tried for EACH bootstrap replicate
%                 (default 3). Kept deliberately small ("light" restarts): with
%                 nBoot=1000 replicates this already means up to 3000 fits, and the
%                 point here is only to keep each replicate from being rigidly
%                 anchored to the full-data solution, not to re-run a full search
%                 1000 times over.
%
%   Procedure
%     1. Fit the sum of K Gaussians to the mean histogram of all molecules, using
%        nStartsFull jittered restarts and keeping the best-converged one.
%     2. Repeat nBoot times: draw nMolecules molecules WITH replacement, average
%        their histograms (molecular weighting), refit starting from step 1's
%        solution AND nStartsBoot-1 light jitters of it (again keeping the best),
%        and record each state's center, sigma and population.
%     3. Report the mean, standard deviation and 95% percentile interval over
%        the replicates (sigma_BOBA).
%
%   Population of state k = area of Gaussian k / total area of all Gaussians, with
%   each area counted only across the histogram's own range (first to last bin edge),
%   so a curve's tail beyond the window does not count as population.
%   Amplitudes are solved by non-negative linear least squares; only centers and
%   widths are searched by fminsearch (as in the GUI, but with lsqnonneg instead
%   of a penalty for negative heights). States are sorted by center in every fit,
%   so state k means the same state in every replicate.
%
%   The per-bin bootstrap standard error of the histogram (results.hist.se) is what
%   the GUI draws as error bars on the histogram bars.
%
%   Output: struct with fields  full, boot, summary, hist, nMolecules  (see code).
%   results.full.nBasins reports how many of the nStartsFull restarts converged to
%   a distinct solution (by residual), as a rough sign of how multi-modal the fit is.

if nargin < 4 || isempty(nBoot),       nBoot       = 1000; end
if nargin < 6 || isempty(nStartsFull), nStartsFull = 8;    end
if nargin < 7 || isempty(nStartsBoot), nStartsBoot = 3;    end
if nargin >= 5 && ~isempty(seed), rng(seed); end

centers    = centers(:)';
startGuess = startGuess(:)';
if mod(numel(startGuess), 2) ~= 0
    error('startGuess must be [mu1 sigma1 mu2 sigma2 ...] (an even number of values).');
end
if size(perMolHist, 2) ~= numel(centers)
    error('perMolHist has %d columns but centers has %d elements.', size(perMolHist, 2), numel(centers));
end
nMol = size(perMolHist, 1);
K    = numel(startGuess) / 2;
if nMol < 2
    error('Need at least 2 molecules to bootstrap.');
end

%% 1. Fit to the full data set (thorough multi-start: don't trust a single startGuess)
yAll = mean(perMolHist, 1);
[full, okFull, nBasins] = fitGaussiansMultiStart(centers, yAll, startGuess, nStartsFull);
full.nBasins = nBasins;   % how many of the nStartsFull restarts landed on a distinct solution
if ~okFull
    warning('bobaFRET:fullFit', 'Fit to the full data set did not converge from any of the %d starting points; check startGuess.', nStartsFull);
elseif nBasins > 1
    warning('bobaFRET:multipleBasins', ...
        '%d of %d restarts on the full data converged to different solutions; the best one was kept, but this suggests the fit is sensitive to the starting guess.', nBasins, nStartsFull);
end
SS_res  = sum((yAll - full.yhat).^2);
SS_tot  = sum((yAll - mean(yAll)).^2);
full.R2 = 1 - SS_res / SS_tot;

%% 2. Bootstrap over molecules
pop = nan(nBoot, K);
mu  = nan(nBoot, K);
sig = nan(nBoot, K);
yBoot = zeros(nBoot, numel(centers));         % bootstrap histograms, for per-bin error bars
for b = 1:nBoot
    idx = randi(nMol, nMol, 1);               % resample molecules with replacement
    y   = mean(perMolHist(idx, :), 1);        % bootstrap cumulated histogram
    yBoot(b, :) = y;
    % Light jittered restarts: always try full.lambda itself, plus a couple of small
    % random perturbations of it, so this replicate isn't rigidly locked to the exact
    % full-data solution's basin (see nStartsBoot note above).
    [fit, ok] = fitGaussiansMultiStart(centers, y, full.lambda, nStartsBoot);
    if ok
        pop(b, :) = fit.pop;
        mu(b, :)  = fit.mu;
        sig(b, :) = fit.sigma;
    end
end
valid   = all(~isnan(pop), 2);
nValid  = sum(valid);
if nValid < 0.9 * nBoot
    warning('bobaFRET:failures', '%d of %d replicates did not converge; results may be unreliable.', nBoot - nValid, nBoot);
end
if nValid < 2
    error('Too few converged replicates to summarize.');
end

%% 3. Summarize
summary.popMean  = mean(pop(valid, :), 1);
summary.popStd   = std(pop(valid, :), 0, 1);          % sigma_BOBA of the population
summary.popCI95  = [pctl(pop(valid, :), 2.5); pctl(pop(valid, :), 97.5)];
summary.muMean    = mean(mu(valid, :), 1);
summary.muStd     = std(mu(valid, :), 0, 1);
summary.muCI95    = [pctl(mu(valid, :), 2.5); pctl(mu(valid, :), 97.5)];     % 95% percentile interval on each center
summary.sigmaMean = mean(sig(valid, :), 1);
summary.sigmaStd  = std(sig(valid, :), 0, 1);
summary.sigmaCI95 = [pctl(sig(valid, :), 2.5); pctl(sig(valid, :), 97.5)];   % 95% percentile interval on each sigma

results.full        = full;
results.boot.pop    = pop;
results.boot.mu     = mu;
results.boot.sigma  = sig;
results.boot.nBoot  = nBoot;
results.boot.nValid = nValid;
results.summary     = summary;
results.nMolecules  = nMol;

% Per-bin uncertainty of the histogram itself (all replicates; independent of the fits)
results.hist.centers = centers;
results.hist.mean    = yAll;
results.hist.se      = std(yBoot, 0, 1);                     % bootstrap standard error of each bin
results.hist.ci95    = [pctl(yBoot, 2.5); pctl(yBoot, 97.5)]; % 95% percentile interval of each bin

%% Print a short report
fprintf('\nBOBA FRET: %d molecules, %d/%d replicates converged. Full-data fit R^2 = %.4f\n', nMol, nValid, nBoot, full.R2);
fprintf('%-6s %-18s %-18s %-24s %-s\n', 'State', 'Center (mean+/-sd)', 'Sigma (mean+/-sd)', 'Population (mean+/-sd)', '95% interval');
for k = 1:K
    fprintf('%-6d %6.3f +/- %-9.3f %6.3f +/- %-9.3f %6.3f +/- %-15.3f [%.3f, %.3f]\n', k, ...
        summary.muMean(k), summary.muStd(k), summary.sigmaMean(k), summary.sigmaStd(k), ...
        summary.popMean(k), summary.popStd(k), summary.popCI95(1, k), summary.popCI95(2, k));
end
fprintf('\n');
end


%% ======================== local functions ========================

function [fit, ok, nBasins] = fitGaussiansMultiStart(t, y, lambda0, nStarts)
% Calls fitGaussians up to nStarts times: once from lambda0 exactly as given, and
% (nStarts-1) more times from randomly jittered versions of lambda0 (jitterLambda,
% below). Keeps whichever converged attempt has the lowest residual (fit.err).
%
% This exists because fminsearch is a local optimizer -- it only ever walks
% downhill from wherever it starts, so a single call can converge to a bad
% solution (e.g. two Gaussians collapsing onto each other) without any warning
% that a better fit existed nearby. Trying several nearby starting points and
% keeping the best is the standard, cheap defense against that.
%
% nBasins counts how many of the converged attempts landed on a solution more
% than 1% (relative) worse, in residual, than the best one -- a rough count of
% how many distinct "basins" the restarts found. nBasins==1 means every restart
% that converged agreed; higher means the fit is sensitive to the starting point.
bestErr = Inf;
fit     = [];
ok      = false;
errs    = [];
for s = 1:nStarts
    if s == 1
        lam0 = lambda0;              % always try the given starting point unperturbed
    else
        lam0 = jitterLambda(lambda0, t);
    end
    [cand, candOk] = fitGaussians(t, y, lam0);
    if candOk
        errs(end+1) = cand.err; %#ok<AGROW>
        if cand.err < bestErr
            bestErr = cand.err;
            fit     = cand;
            ok      = true;
        end
    end
end
if ~ok
    % nothing converged acceptably from any start; fall back to the plain,
    % unperturbed-start result so the caller still gets a consistent struct
    % (its own ok=false is what signals the failure).
    [fit, ok] = fitGaussians(t, y, lambda0);
    nBasins = 0;
    return
end
nBasins = sum(errs > 1.01 * bestErr) + 1;   % +1 counts the best solution's own basin
end

function lam = jitterLambda(lam0, t)
% One randomly-perturbed copy of lam0 = [mu1 sigma1 mu2 sigma2 ...], used to give
% fitGaussiansMultiStart a nearby-but-different place to start from. Centers are
% nudged by up to +/-7.5% of the data span; sigmas are scaled by a random factor
% in [0.7, 1.3]. Both are clamped back into the range fitError treats as valid.
K      = numel(lam0) / 2;
mu0    = lam0(1:2:end);
sigma0 = abs(lam0(2:2:end));
span   = t(end) - t(1);
dt     = t(2) - t(1);
mu     = mu0 + (rand(1, K) - 0.5) * 0.15 * span;
sigma  = sigma0 .* (0.7 + 0.6 * rand(1, K));
mu     = min(max(mu, t(1)), t(end));
sigma  = max(sigma, dt);
lam    = reshape([mu; sigma], 1, []);
end

function [fit, ok] = fitGaussians(t, y, lambda0)
% Fit a sum of Gaussians to histogram y (centers t). lambda = [mu1 sigma1 mu2 sigma2 ...].
% sigma is the standard deviation -- the only width unit anywhere in this function;
% there is no FWHM conversion in this fit at all. Returns states sorted by center.
opts = optimset('TolX', 1e-4, 'TolFun', 1e-4, 'MaxIter', 5000, ...
                'MaxFunEvals', 5000 * numel(lambda0), 'Display', 'off');
[lambda, ~, flag] = fminsearch(@(l) fitError(l, t, y), lambda0, opts);
[err, c] = fitError(lambda, t, y);
ok = (flag > 0) && (err < 1e2);               % err = 1e3 marks an out-of-bounds solution

mu    = lambda(1:2:end);
sigma = abs(lambda(2:2:end));
[mu, order] = sort(mu);
sigma = sigma(order);
c     = c(order)';
% Area under each Gaussian, counted only across the histogram's own range. The full
% closed-form area is c*sigma*sqrt(2*pi) (x from -Inf to +Inf); frac is the share of
% each curve that actually falls inside the first-to-last bin edges, so tails beyond
% the window (e.g. a broad state near FRET 1.2) are not counted as population.
dt    = t(2) - t(1);
lo    = t(1)   - dt/2;                % outer edge of the first bin
hi    = t(end) + dt/2;                % outer edge of the last bin
frac  = 0.5 * (erf((hi - mu) ./ (sigma * sqrt(2))) - erf((lo - mu) ./ (sigma * sqrt(2))));
area  = c .* sigma .* sqrt(2*pi) .* frac;

fit.lambda = reshape([mu; sigma], 1, []);      % sorted, ready to use as next start
fit.mu     = mu(:)';
fit.sigma  = sigma(:)';                        % native parameter -- not derived from anything
fit.amp    = c;
fit.pop    = area / sum(area);
fit.err    = err;   % residual norm of this fit, used by fitGaussiansMultiStart to pick the best restart
fit.yhat   = zeros(size(t));
for k = 1:numel(mu)
    fit.yhat = fit.yhat + c(k) * gaussianSigma(t, mu(k), sigma(k));
end
if sum(area) <= 0 || any(isnan(fit.pop))
    ok = false;
end
end

function [err, c] = fitError(lambda, t, y)
% Residual norm with amplitudes solved by non-negative least squares.
K     = numel(lambda) / 2;
mu    = lambda(1:2:end);
sigma = abs(lambda(2:2:end));
dt    = t(2) - t(1);
% Keep centers inside the histogram and sigma between one bin and the full span.
% (A sigma narrower than one bin can't be resolved by this histogram anyway; one
% wider than the whole span isn't really a localized state any more.)
if any(mu < t(1)) || any(mu > t(end)) || any(sigma < dt) || any(sigma > (t(end) - t(1)))
    err = 1e3;
    c   = zeros(K, 1);
    return
end
A = zeros(numel(t), K);
for k = 1:K
    A(:, k) = gaussianSigma(t, mu(k), sigma(k))';
end
c   = lsqnonneg(A, y(:));
err = norm(A * c - y(:));
end

function g = gaussianSigma(x, peakPosition, sigma)
% Gaussian parameterized directly by its standard deviation -- the native,
% only width unit used throughout this file. No FWHM conversion anywhere.
g = exp(-((x - peakPosition) .^ 2) ./ (2 .* sigma .^ 2));
end

function q = pctl(X, p)
% Percentile p (0-100) of each column, linear interpolation (no toolbox needed).
X   = sort(X, 1);
n   = size(X, 1);
pos = 1 + (n - 1) * p / 100;
lo  = floor(pos);
hi  = ceil(pos);
q   = X(lo, :) + (pos - lo) .* (X(hi, :) - X(lo, :));
end
