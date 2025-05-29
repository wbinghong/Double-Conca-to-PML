function [A, q, S, mu_row] = generate_random_mu_random_q_matrix(n, alpha, mu1, mu2, lambda)
%GENERATE_ROWTYPE_MATRIX Generate a random matrix under the row-type model
%
%   This function generates an n x n matrix A where:
%     - Each column is based on a common i.i.d. random base vector q_j ~ Exp(lambda)
%     - Each row is assigned a type S_i ∈ {1, 2} with probability alpha and 1 - alpha
%     - The i-th row of A is [q_1^{mu_{S_i}}, ..., q_n^{mu_{S_i}}]
%
%   INPUTS:
%       n      - size of the square matrix (matrix is n x n)
%       alpha  - probability Pr(S_i = 1) for choosing exponent mu1
%       mu1    - exponent for row type S_i = 1
%       mu2    - exponent for row type S_i = 2
%       lambda - rate parameter of the exponential distribution Exp(lambda)
%
%   OUTPUTS:
%       A       - generated matrix of size n x n, with A(i,j) = q_j^{mu_{S_i}}
%       q       - base vector q_j ~ Exp(lambda), sorted in descending order
%       S       - row type vector S_i ∈ {1, 2}, sorted to match descending mu
%       mu_row  - vector of per-row exponents [mu1 or mu2], descending order

    format long g;

    % --- Step 1: Sample base vector q_j ~ Exp(lambda), sort descending ---
    q = -log(rand(n, 1)) / lambda;    % inverse transform sampling
    q = sort(q, 'descend');           % for reproducibility and structure
    % q = q / sum(q);                 % normalization disabled intentionally

    % --- Step 2: Sample row types S_i ∈ {1, 2} with Pr(S_i = 1) = alpha ---
    S = datasample([1, 2], n, 'Weights', [alpha, 1 - alpha]);

    % --- Step 3: Map row types to exponent vector mu_row, sort by mu descending ---
    mu_row = mu1 * (S == 1) + mu2 * (S == 2);      % mu_row(i) = mu1 or mu2
    [mu_row, sort_idx] = sort(mu_row, 'descend'); % sort: larger exponents first
    S = S(sort_idx);                              % ensure S and mu_row aligned

    % --- Step 4: Construct A(i,j) = q_j^{mu_row(i)} ---
    Qmat   = repmat(q.', n, 1);          % replicate q_j across rows
    Mumat  = repmat(mu_row.', 1, n);       % replicate mu_row across columns
    epsilon = 1e-12;                     % safeguard against numerical zero
    A = Qmat .^ Mumat + epsilon;         % element-wise exponentiation
end
