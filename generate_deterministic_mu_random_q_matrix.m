function [A, q, S, mu_row] = generate_deterministic_mu_random_q_matrix(n, alpha, mu1, mu2, lambda)
%GENERATE_DETERMINISTIC_ROWTYPE_MATRIX Generate a matrix with deterministic row-type split
%
%   This function generates an n x n matrix A where:
%     - Each q_j ~ Exp(lambda), sorted in descending order
%     - First floor(alpha * n) rows are assigned mu1
%     - Remaining rows are assigned mu2
%     - A(i,j) = q_j^{mu_row(i)} where mu_row is fixed
%
%   INPUTS:
%       n      - size of square matrix
%       alpha  - fraction of rows to be assigned mu1 (remaining gets mu2)
%       mu1    - exponent for first block
%       mu2    - exponent for second block
%       lambda - exponential rate parameter (mean = 1/lambda)
%
%   OUTPUTS:
%       A       - generated matrix of size n x n
%       q       - q_j vector, sampled ~ Exp(lambda), sorted descending
%       S       - row type vector, 1 for mu1 rows, 2 for mu2 rows
%       mu_row  - row-wise exponent vector [mu1 or mu2], fixed deterministic

    format long g;

    % Step 1: Sample q_j ~ Exp(lambda), sorted descending
    q = -log(rand(n, 1)) / lambda;
    q = sort(q, 'descend');

    % Step 2: Construct deterministic mu_row vector
    num_mu1 = floor(alpha * n);
    num_mu2 = n - num_mu1;
    mu_row = [repmat(mu1, num_mu1, 1); repmat(mu2, num_mu2, 1)];
    S = [repmat(1, num_mu1, 1); repmat(2, num_mu2, 1)];  % row type indicator

    % Step 3: Construct A(i,j) = q_j^{mu_row(i)}
    Qmat   = repmat(q.', n, 1);         % replicate q across rows
    Mumat  = repmat(mu_row, 1, n);      % replicate mu_row across columns
    epsilon = 1e-12;
    A = Qmat .^ Mumat + epsilon;        % element-wise exponentiation
end
