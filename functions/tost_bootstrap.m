function results = tost_bootstrap(x, mu0, Delta, alpha, statfun, B)
% Nonparametric one-sample TOST using percentile bootstrap
%
% Tests equivalence of statfun(x) to mu0 within ±Delta
%
% H0: stat <= mu0 - Delta  OR  stat >= mu0 + Delta
% H1: mu0 - Delta < stat < mu0 + Delta
%
% Inputs:
%   x        : data vector
%   mu0      : reference value (e.g. 1)
%   Delta    : equivalence margin (e.g. 0.05)
%   alpha    : significance level (default 0.05)
%   statfun  : @mean or @median (default @mean)
%   B        : number of bootstrap samples (default 10000)
%
% Output:
%   results  : structure with test summary

if nargin < 6 || isempty(B),       B = 10000;      end
if nargin < 5 || isempty(statfun), statfun = @mean; end
if nargin < 4 || isempty(alpha),   alpha = 0.05;   end

x = x(~isnan(x));
n = numel(x);

theta_hat = statfun(x);

% Bootstrap resampling
bootstat = zeros(B,1);
for b = 1:B
    xb = x(randi(n, n, 1));   % resample with replacement
    bootstat(b) = statfun(xb);
end

% Sort bootstrap distribution
bootstat = sort(bootstat);

% Percentile CI for equivalence: (1 - 2*alpha)
lo_idx = floor((alpha)*B);
hi_idx = ceil((1 - alpha)*B);

CI = [bootstat(lo_idx), bootstat(hi_idx)];

% Equivalence bounds
low = mu0 - Delta;
up  = mu0 + Delta;

% Equivalence decision
equivalent = (CI(1) > low) && (CI(2) < up);

% Approximate one-sided bootstrap p-values
p_lower = mean(bootstat <= low);
p_upper = mean(bootstat >= up);

% Store results
results.theta_hat = theta_hat;
results.n = n;
results.CI = CI;
results.Delta = Delta;
results.low = low;
results.up = up;
results.p_lower = p_lower;
results.p_upper = p_upper;
results.equivalent = equivalent;
results.statfun = func2str(statfun);
results.B = B;

% Display
fprintf('\nNonparametric Bootstrap TOST\n');
fprintf('Statistic (%s) = %.6f, n = %d\n', results.statfun, theta_hat, n);
fprintf('Equivalence bounds: [%.6f, %.6f]\n', low, up);
fprintf('%.1f%% CI: [%.6f, %.6f]\n', 100*(1-2*alpha), CI(1), CI(2));
fprintf('Bootstrap p_lower = %.4g\n', p_lower);
fprintf('Bootstrap p_upper = %.4g\n', p_upper);

if equivalent
    fprintf('=> Result: STATISTICALLY EQUIVALENT to %.3f\n', mu0);
else
    fprintf('=> Result: NOT EQUIVALENT to %.3f\n', mu0);
end
end