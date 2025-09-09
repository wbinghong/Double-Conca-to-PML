% BLOCK vs BASE–POWER (PML) — BETHE RATIO VS n (log–log)
% --------------------------------------------------------
% Goal: On the same figure, plot the empirical curves for two matrix models:
%   (1) Two-block i.i.d. row model (Uniform blocks)
%   (2) Two-type Base–Power (Random) PML matrix:  A(i,j) = q_{T_i}^{mu_{S_j}}
% and compare them against the theory line (pi*n/e)^(1/4) using yellow circles.
%
% Requirements on path:
%   perm_ryser_gray.m
%   computeBethePermanent2_streaming.m   % (M=2, streaming covers)

clear; clc; close all;

%% ----------- Sanity check: Ryser == recursive perm (small N) -----------
for Ntest = 3:6
    for trial = 1:10
        Atest = rand(Ntest);
        % If you have your recursive perm(A,n), you can uncomment to verify equality
        % p_rec = perm(Atest, Ntest);
        p_rgs = perm_ryser_gray(Atest);
        % if abs(p_rec - p_rgs) > 1e-8 * max(1, p_rec)
        %     error('Mismatch at N=%d: rec=%g vs ryser=%g', Ntest, p_rec, p_rgs);
        % end
    end
end
fprintf('Ryser–Gray permanent running OK.\n');

%% ---------------------- Experiment knobs ----------------------
n_list    = 3:12;     % beware runtime for large n (perm on 2n x 2n)
trials    = 1;       % MC trials per n (increase for smoother)
K_samples = 200;     % covers per Bethe estimate (set 0 for auto/enum when tiny)

rng(7, 'twister');

%% ---------------------- Two-block (Uniform) params ------------
% A_top ~ U[L_A,U_A], A_bottom ~ U[L_B,U_B]
L_A = 0.6; U_A = 1.4;    % block-A (top)  mean 1.0, width 0.8
L_B = 0.3; U_B = 0.9;    % block-B (bottom) mean 0.6, width 0.6
theta = 0.6;             % n1 = floor(theta*n)

sample_DA = @(r,c) (U_A - L_A) .* rand(r,c) + L_A;
sample_DB = @(r,c) (U_B - L_B) .* rand(r,c) + L_B;

%% ---------------------- Base–Power (PML) params ---------------
%   A(i,j) = q_{T_i}^{mu_{S_j}},  S_j ~ Cat(alpha),  T_i ~ Cat(beta)
mu1 = 2; mu2 = 1;            % row-type exponents
q1  = 0.7; q2  = 0.3;        % column-type bases
alpha1 = 0.3;                 % P(S=1)=alpha1, P(S=2)=1-alpha1
beta1  = 0.4;                 % P(T=1)=beta1,  P(T=2)=1-beta1
mu = [mu1, mu2];
q  = [q1,  q2];

%% ---------------------- Storage -------------------------------
numN        = numel(n_list);
ratio_block = nan(1, numN);   % empirical ratio for block model
ratio_pml   = nan(1, numN);   % empirical ratio for base–power model
ratio_th    = nan(1, numN);   % theory (pi n / e)^(1/4)

%% ---------------------- Main loop -----------------------------
for idx = 1:numN
    n = n_list(idx);

    % ---- Two-block accumulators ----
    perm_sq_block  = 0.0;
    bethe_sq_block = 0.0;

    % ---- PML accumulators ----
    perm_sq_pml  = 0.0;
    bethe_sq_pml = 0.0;

    for t = 1:trials
        %% ===== (1) TWO-BLOCK UNIFORM =====
        n1 = floor(theta * n); n2 = n - n1;
        A_block = zeros(n,n);
        A_block(1:n1,     :) = sample_DA(n1, n);
        A_block(n1+1:end, :) = sample_DB(n2, n);

        p_b  = perm_ryser_gray(A_block);
        pb_b = computeBethePermanent2_streaming(A_block, n, K_samples);
        perm_sq_block  = perm_sq_block  + p_b^2;
        bethe_sq_block = bethe_sq_block + pb_b^2;

        %% ===== (2) BASE–POWER (PML) =====
        % Sample types (ensure shapes): T is n×1 (rows), S is 1×n (cols)
        T = 2*ones(n,1);  T(rand(n,1) <= beta1) = 1;   % P(T_i=1)=beta1
        S = 2*ones(1,n);  S(rand(1,n) <= alpha1) = 1;  % P(S_j=1)=alpha1
        qT  = q(T);                                    % n×1
        muS = mu(S);                                   % 1×n
        for i = 1:n
            for j = 1:n
                A_pml(i, j) = qT(i) ^ muS(j);
            end
        end
        p_p  = perm_ryser_gray(A_pml);
        pb_p = computeBethePermanent2_streaming(A_pml, n, K_samples);
        perm_sq_pml  = perm_sq_pml  + p_p^2;
        bethe_sq_pml = bethe_sq_pml + pb_p^2;
    end

    % ---- MC expectations and ratios ----
    E_perm_sq_block  = perm_sq_block  / trials;
    E_bethe_sq_block = bethe_sq_block / trials;
    ratio_block(idx) = sqrt(E_perm_sq_block) / sqrt(E_bethe_sq_block);

    E_perm_sq_pml  = perm_sq_pml  / trials;
    E_bethe_sq_pml = bethe_sq_pml / trials;
    ratio_pml(idx) = sqrt(E_perm_sq_pml) / sqrt(E_bethe_sq_pml);

    ratio_th(idx) = (pi*n/exp(1))^(1/4);

    fprintf('[n=%2d]  ratio_block=%.4g,  ratio_pml=%.4g,  theory=%.4g\n', ...
        n, ratio_block(idx), ratio_pml(idx), ratio_th(idx));
end

%% ---------- Plot ----------
figure;
loglog(n_list, ratio_pml, 'b--o','LineWidth',2); hold on;
loglog(n_list, ratio_block,'r--o','LineWidth',2);
loglog(n_list, ratio_th, '-o', 'LineWidth', 2);

lg = legend({'Two-block i.i.d. Matrix', ...
'Two-type Base-power Random Matrix', ...
'$\sqrt[4]{\frac{\pi n}{e}}$'}, ...
'Interpreter','latex','Location','northwest');
% Legend title per request
try title(lg,'Interested Ratio $\frac{\sqrt{\mathrm{E}[\mathrm{perm}(\mathbf{A})^2]}}{\sqrt{\mathrm{E}[\mathrm{perm}_{\mathrm{B},2}(\mathbf{A})^2]}}$ for'); end

grid on;

% Set axes to data range precisely
x_min = min(n_list); x_max = max(n_list);
y_min = min([min(ratio_pml), min(ratio_block), min(ratio_th)]);
y_max = max([max(ratio_pml), max(ratio_block), max(ratio_th)]);
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
fname = sprintf('results/bethe_ratio_block_vs_pml.pdf');
exportgraphics(gcf, fname, 'ContentType','vector', 'Resolution',300);
fprintf('[Saved] %s\n', fname);


