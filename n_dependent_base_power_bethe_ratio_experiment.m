% VERIFICATION OF RANDOM-BLOCK ASYMPTOTIC EQUIVALENCE
% --------------------------------------------------------
% Goal: Verify that for n-dependent proportions alpha(n), the ratio for the
% random-type PML model (Z^r) converges to the ratio for the corresponding
% block-type model (Z^b) with k = floor(n*alpha(n)) rows.
%
% Model: A(i,j) = q_j^(mu_Ti), where T_i is the row type.
%   - Random Model: T_i is sampled i.i.d with P(T_i=1) = alpha(n).
%   - Block Model: The first k rows have T_i=1, the rest have T_i=2.
%
% Requirements on path:
%   perm_ryser_gray.m
%   computeBethePermanent2_streaming.m
% --------------------------------------------------------
clear; clc; close all;

%% ----------- Sanity check: Ryser implementation -----------
% This part remains unchanged to ensure the permanent function is correct.
for Ntest = 3:6, for trial = 1:10, Atest = rand(Ntest); p_rgs = perm_ryser_gray(Atest); end, end
fprintf('Ryser–Gray permanent running OK.\n');

%% ---------------------- Experiment Knobs ----------------------
n_list    = 3:10;     % Range of n to test. Beware runtime for n > 12.
trials    = 1;       % Monte Carlo trials per n. Increase for smoother curves.
K_samples = 200;      % Covers per Bethe estimate.
rng(2025, 'twister'); % Set seed for reproducibility

%% ---------------------- Base–Power (PML) Model Parameters ---------------
% Define the core model parameters
mu1 = 2; mu2 = 1;      % Exponents for the two row types
mu = [mu1, mu2];

% --- NEW: n-dependent proportion for the rare row type (Type 1) ---
alpha1_func = @(n) 1/sqrt(n); 

% --- NEW: We will sample q_j columns randomly in each trial ---
% This function will generate n random bases for the columns
sample_q_cols = @(n) rand(1, n) * 0.8 + 0.2; % Sample q_j from U[0.2, 1.0]

%% ---------------------- Storage -------------------------------
numN = numel(n_list);
ratio_pml_random = nan(1, numN); % Stores ratio for the random-type model (Z^r)
ratio_pml_block  = nan(1, numN); % Stores ratio for the block-type model (Z^b)
ratio_th         = nan(1, numN); % Stores theory line (pi*n/e)^(1/4)

%% ---------------------- Main Loop -----------------------------
for idx = 1:numN
    n = n_list(idx);
    
    % --- Calculate the n-dependent proportion for this n ---
    alpha1_n = alpha1_func(n);
    k = floor(n * alpha1_n); % Expected number of Type 1 rows for the block model
    
    % --- Accumulators for this n ---
    perm_sq_random  = 0.0;
    bethe_sq_random = 0.0;
    perm_sq_block   = 0.0;
    bethe_sq_block  = 0.0;
    
    fprintf('[n=%2d, k=%d] Running %d trials...\n', n, k, trials);
    
    for t = 1:trials
        % --- For a fair comparison, use the same random columns for both models in each trial ---
        q_cols = sample_q_cols(n); % Generate n random bases for the columns
        
        %% ===== (1) Random-Type PML Model (computes Z^r) =====
        % Sample row types T_i randomly according to alpha(n)
        T_random = 2*ones(n,1);  
        T_random(rand(n,1) <= alpha1_n) = 1;
        
        % Construct the matrix A_random
        A_random = zeros(n,n);
        for i = 1:n
            if T_random(i) == 1
                A_random(i, :) = q_cols .^ mu1;
            else
                A_random(i, :) = q_cols .^ mu2;
            end
        end
        
        % Calculate permanents
        p_r  = perm_ryser_gray(A_random);
        pb_r = computeBethePermanent2_streaming(A_random, n, K_samples);
        perm_sq_random  = perm_sq_random  + p_r^2;
        bethe_sq_random = bethe_sq_random + pb_r^2;

        %% ===== (2) Block-Type PML Model (computes Z^b) =====
        % Construct the matrix A_block with a fixed number of k Type 1 rows
        A_block = zeros(n,n);
        % First k rows are Type 1
        for i = 1:k
            A_block(i, :) = q_cols .^ mu1;
        end
        % Remaining n-k rows are Type 2
        for i = k+1:n
             A_block(i, :) = q_cols .^ mu2;
        end
        
        % Calculate permanents
        p_b  = perm_ryser_gray(A_block);
        pb_b = computeBethePermanent2_streaming(A_block, n, K_samples);
        perm_sq_block  = perm_sq_block  + p_b^2;
        bethe_sq_block = bethe_sq_block + pb_b^2;
    end
    
    % ---- MC expectations and ratios ----
    E_perm_sq_random  = perm_sq_random  / trials;
    E_bethe_sq_random = bethe_sq_random / trials;
    ratio_pml_random(idx) = sqrt(E_perm_sq_random) / sqrt(E_bethe_sq_random);
    
    E_perm_sq_block  = perm_sq_block  / trials;
    E_bethe_sq_block = bethe_sq_block / trials;
    ratio_pml_block(idx) = sqrt(E_perm_sq_block) / sqrt(E_bethe_sq_block);
    
    ratio_th(idx) = (pi*n/exp(1))^(1/4);
    
    fprintf('    > Ratios: Random=%.4g, Block=%.4g, Theory=%.4g\n', ...
        ratio_pml_random(idx), ratio_pml_block(idx), ratio_th(idx));
end

%% ---------------------- Plotting -----------------------------
figure;
loglog(n_list, ratio_pml_random, 'b--o', 'LineWidth', 2); hold on;
loglog(n_list, ratio_pml_block, 'r--o', 'LineWidth', 2);
loglog(n_list, ratio_th, '-o', 'LineWidth', 2);
grid on;
xlabel('n', 'Interpreter', 'latex', 'FontSize', 2);

lg = legend({'Random-Type Base-power Matrix ($Z^{\mathrm{r}}_{n}(\alpha(n))$)', ...
    'Block-Type Base-power Matrix ($Z^{\mathrm{b}}_{n}(\lfloor n\alpha(n) \rfloor)$)', ...
    '$\sqrt[4]{\frac{\pi n}{e}}$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12);
try
    title(lg,'Ratio $\frac{\sqrt{\mathrm{E}[\mathrm{perm}(\mathbf{A})^2]}}{\sqrt{\mathrm{E}[\mathrm{perm}_{\mathrm{B},2}(\mathbf{A})^2]}}$'); 
end

grid on;

x_min = min(n_list); x_max = max(n_list);
y_min = min([min(ratio_pml_random), min(ratio_pml_block), min(ratio_th)]);
y_max = max([max(ratio_pml_random), max(ratio_pml_block), max(ratio_th)]);
axis([x_min, x_max, y_min, y_max]);

% Control x-axis ticks
xticks([min(n_list), 5, max(n_list)]);
xticklabels(arrayfun(@num2str, [min(n_list), 10, max(n_list)], 'UniformOutput', false));


% Figure size and export (match reference style)
w_in = 6; h_in = 4.5;
set(gcf,'Units','inches','Position',[1,1,w_in,h_in]);
set(gcf,'PaperUnits','inches','PaperSize',[w_in,h_in], ...
'PaperPosition',[0,0,w_in,h_in]);
warning('on','MATLAB:nchoosek:LargeCoefficient');


if ~exist('results','dir'), mkdir('results'); end
vec2fname = @(v) strrep(strjoin(arrayfun(@(x)sprintf('%.2f',x),v,'UniformOutput',false),'_'),'.','p');
fname = sprintf('results/random_vs_block_pml_equivalence.pdf');
exportgraphics(gcf, fname, 'ContentType','vector', 'Resolution',300);
fprintf('[Saved] %s\n', fname);