function [permBM, stats] = computeBethePermanent2_streaming(A, n, K_samples, varargin)
% computeBethePermanent2_streaming  Bethe permanent for M=2 via streaming covers.
%
%   [permBM, stats] = computeBethePermanent2_streaming(A, n, K_samples, ...)
%
% Computes the 2-cover Bethe permanent of A by averaging permanents over
% all (or sampled) expanded covers, then taking the 1/2 power.
% - Exact full enumeration if feasible, otherwise unbiased Monte Carlo.
% - Uses a memory-light Ryser–Gray permanent (requires perm_ryser_gray.m on path).
%
% INPUTS
%   A          : (n x n) nonnegative matrix.
%   n          : dimension of A (integer, consistent with size(A)).
%   K_samples  : if 0 or empty -> try exact enumeration when reasonable;
%                if >0         -> sample this many covers uniformly at random.
%
%   Name-Value options:
%     'MaxEnumBits' (default 22) :
%           Enumerate exhaustively only when Nexp=(n-1)^2 <= MaxEnumBits.
%     'ReturnCoverMean' (default false):
%           If true, stats.meanPerm stores the average permanent over covers.
%
% OUTPUTS
%   permBM     : Bethe-2 permanent estimate, i.e., sqrt( mean_c perm(A_c) ).
%   stats      : struct with diagnostic info:
%                  .mode   : 'enum' or 'sample'
%                  .Nexp   : number of free 2x2 blocks = (n-1)^2
%                  .Nsamp  : actual number of evaluated covers
%                  .meanPerm (if requested): mean of perm(A_c)
%                  .stderr (sampling only): std(err) of sqrt(mean) via delta
%
% EXAMPLE
%   % A = rand(4); n=4; exact enumeration (since (n-1)^2=9 <= 22)
%   [pb, st] = computeBethePermanent2_streaming(A, 4, 0);
%
%   % Large n -> sampling 1000 covers
%   [pb2, st2] = computeBethePermanent2_streaming(A, 12, 1000);
%
% Author: Binghong Wu (original idea), refactored by ChatGPT helper
% Date  : 2025-09-03

% -------------------- sanity checks --------------------
if nargin < 3 || isempty(K_samples), K_samples = 0; end
p = inputParser;
p.addParameter('MaxEnumBits', 22, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('ReturnCoverMean', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
MaxEnumBits     = p.Results.MaxEnumBits;
ReturnCoverMean = logical(p.Results.ReturnCoverMean);

if size(A,1) ~= n || size(A,2) ~= n
    error('A must be n-by-n with the given n.');
end
if exist('perm_ryser_gray','file') ~= 2
    error('perm_ryser_gray.m not found on path. Please add it or place alongside this file.');
end

% -------------------- problem dimensions --------------------
Nexp = (n-1)*(n-1);    % number of free (non-first-row/col) 2x2 blocks
nM   = 2*n;            % expanded dimension for M=2

% -------------------- reusable buffers --------------------
Aexp = zeros(nM, nM);           % expanded matrix buffer
P0 = eye(2); P1 = [0 1; 1 0];   % 2x2 permutation blocks
r0 = (0:n-1)*2 + 1; r1 = r0+1;  % row-block indices
c0 = r0;           c1 = r1;     % col-block indices

% -------------------- builders --------------------
    function build_from_bits_scalar(bits)
        % bits: nonnegative scalar (double), read kth bit by floor/div/mod
        for ii = 1:n
            for jj = 1:n
                R = r0(ii):r1(ii); C = c0(jj):c1(jj);
                if ii==1 || jj==1
                    Aexp(R,C) = A(ii,jj) * P0;
                else
                    k = (ii-2)*(n-1) + (jj-2); % 0-based index
                    useSwap = mod(floor(bits / 2^k), 2) ~= 0;
                    if useSwap
                        Aexp(R,C) = A(ii,jj) * P1;
                    else
                        Aexp(R,C) = A(ii,jj) * P0;
                    end
                end
            end
        end
    end

    function build_from_bitvec(bv)
        % bv: logical/0-1 row vector of length Nexp
        for ii = 1:n
            for jj = 1:n
                R = r0(ii):r1(ii); C = c0(jj):c1(jj);
                if ii==1 || jj==1
                    Aexp(R,C) = A(ii,jj) * P0;
                else
                    k = (ii-2)*(n-1) + (jj-2); % 0-based
                    if bv(k+1)
                        Aexp(R,C) = A(ii,jj) * P1;
                    else
                        Aexp(R,C) = A(ii,jj) * P0;
                    end
                end
            end
        end
    end

% -------------------- choose enum vs sample --------------------
use_enum = (K_samples==0) && (Nexp <= MaxEnumBits);

sum_perm = 0.0;
ns       = 0;      % actual number of evaluated covers

if use_enum
    total = 2^Nexp;  % safe since Nexp <= MaxEnumBits
    for bits = 0:total-1
        build_from_bits_scalar(bits);
        sum_perm = sum_perm + perm_ryser_gray(Aexp);
    end
    ns = total;
    avg_perm = sum_perm / ns;           % exact mean over all covers
    mode = 'enum';
else
    if K_samples <= 0
        % Auto-choose a sample size if user asked for full but infeasible
        K_samples = min(2000, max(200, 50*Nexp));
    end
    % One-pass sampling, no giant integers: sample Nexp Bernoulli(1/2) bits
    sum_perm_sq = 0.0;  % for variance/SE estimate
    for t = 1:K_samples
        bv = rand(1, Nexp) < 0.5;    % random cover bits
        build_from_bitvec(bv);
        val = perm_ryser_gray(Aexp);
        sum_perm    = sum_perm    + val;
        sum_perm_sq = sum_perm_sq + val*val;
    end
    ns = K_samples;
    avg_perm = sum_perm / ns;        % unbiased estimate of mean over covers
    mode = 'sample';
end

% -------------------- Bethe-2 aggregation --------------------
permBM = sqrt(avg_perm);             % (mean)^(1/2) for M=2

% -------------------- stats --------------------
if nargout >= 2
    stats = struct('mode',mode,'Nexp',Nexp,'Nsamp',ns);
    if ReturnCoverMean
        stats.meanPerm = avg_perm;
    end
    if strcmp(mode,'sample')
        % SE for sqrt(E[X]) via delta method around sample mean
        mu_hat = avg_perm;
        var_hat = (sum_perm_sq/ns - mu_hat^2) / max(1, ns-1);
        % d/dmu sqrt(mu) = 1/(2 sqrt(mu))
        if mu_hat > 0
            stats.stderr = 0.5 / sqrt(mu_hat) * sqrt(var_hat/ns);
        else
            stats.stderr = NaN;
        end
    end
end
end
