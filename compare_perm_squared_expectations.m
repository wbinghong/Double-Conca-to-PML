clear; clc;
format long g;

% --- Global Parameters
alpha       = 0.5;
mu1         = 2;
mu2         = 1;
lambda      = 1;
num_trials  = 100;
n_range     = 3:8;

% Some fixed quantities
C11 = alpha * factorial(2*mu1) / lambda^(2*mu1);
C22 = (1 - alpha) * factorial(2*mu2) / lambda^(2*mu2);
C12 = sqrt(alpha * (1 - alpha)) * factorial(mu1 + mu2) / lambda^(mu1 + mu2);
S = [C11, C12;
     C12, C22];
lambda_max = max(eig(S));

% Storage to a table
summary_table = zeros(length(n_range), 5);  % columns: n, E[perm(M)^2], E[perm(A)^2], tr(S^k), lambda^k

for idx = 1:length(n_range)
    n = n_range(idx);
    fprintf('\n===== n = %d =====\n', n);

    % Step 1: \expval{\perm(\matr{M})^2}
    perm_M_list = zeros(num_trials, 1);
    for t = 1:num_trials
        [M, ~, ~, ~] = generate_deterministic_mu_random_q_matrix(n, alpha, mu1, mu2, lambda);
        perm_M_list(t) = perm(M, n);
    end
    mean_perm_M2 = mean(perm_M_list.^2);

    % Step 2: \expval{\perm(\matr{A})^2}
    perm_A_list = zeros(num_trials, 1);
    for t = 1:num_trials
        [A, ~, ~, ~] = generate_random_mu_random_q_matrix(n, alpha, mu1, mu2, lambda);
        perm_A_list(t) = perm(A, n);
    end
    mean_perm_A2 = mean(perm_A_list.^2);

    % Step 3: trace(S^k) and lambda_max^k
    trace_S = zeros(1, n);
    trace_S_pf = lambda_max.^(1:n);
    S_k = eye(2);
    for k = 1:n
        S_k = S_k * S;
        trace_S(k) = trace(S_k);
    end

    % Step 4: Cycle index substitution
    Z = cycle_index_Sn(n);
    Z_syms = sym('z', [1 n]);

    Z_substituted_trace = subs(Z, Z_syms, trace_S);
    Z_substituted_pf    = subs(Z, Z_syms, trace_S_pf);

    perm_from_S_trace = factorial(n)^2 * double(Z_substituted_trace);
    perm_from_S_pf    = factorial(n)^2 * double(Z_substituted_pf);

    % Store values
    summary_table(idx, :) = [n, mean_perm_M2, mean_perm_A2, perm_from_S_trace, perm_from_S_pf];
end

% --- Plotting -------------------------------------------------------------
figure;
semilogy(summary_table(:,1), summary_table(:,2), 'o-', 'LineWidth', 2, 'DisplayName', 'E[perm(M)^2]');
hold on;
semilogy(summary_table(:,1), summary_table(:,3), 's-', 'LineWidth', 2, 'DisplayName', 'E[perm(A)^2]');
semilogy(summary_table(:,1), summary_table(:,4), '^-', 'LineWidth', 2, 'DisplayName', 'Symbolic via tr(S^k)');
semilogy(summary_table(:,1), summary_table(:,5), 'd--', 'LineWidth', 2, 'DisplayName', 'Symbolic via \lambda_{max}^k');
grid on;
xlabel('n');
ylabel('Estimate of E[perm(.)^2] (log scale)');
title('Comparison of perm^2 Approximations vs n');
legend('Location', 'northwest');
set(gca, 'FontSize', 12);

% --- Save figure ----------------------------------------------------------
output_dir = 'results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

set(gcf, 'Units', 'inches', 'Position', [1, 1, 6, 4]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6, 4]);
set(gcf, 'PaperPosition', [0, 0, 6, 4]);
set(gcf, 'PaperPositionMode', 'manual');

filename = sprintf('squared_perm_comparison_alpha_%.2f_mu1_%d_mu2_%d_lambda_%.2f_n_%dto%d_trials_%d.pdf', ...
                   alpha, mu1, mu2, lambda, n_range(1), n_range(end), num_trials);
save_path = fullfile(output_dir, filename);
print(gcf, save_path, '-dpdf', '-r300');
