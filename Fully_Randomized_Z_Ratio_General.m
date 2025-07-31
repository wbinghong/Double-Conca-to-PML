% ============================================================
%  Multitype (u × v)   指定 (α⃗ , β⃗ ) 组合的 Z_B / Z 曲线
%  • 采用  Tα·Z·Tβᵀ = Zᵣ  直接左右求解（无 kron/lsqminnorm，显著提速）
%  • μ、q 皆降序；结果保存为 PDF
% ============================================================
clear; clc; warning('off','MATLAB:nchoosek:LargeCoefficient');

% ---------- 用户参数 ---------- -----------------------------------------
u = 3;  v = 3;                       % 行 / 列 类型数
alpha_target = [0.6 0.2 0.2];        % ⃗α  (长度=u, 和=1)
beta_target  = [0.4 0.4 0.2];        % ⃗β  (长度=v, 和=1)
n_max = 30;                          % 最大矩阵阶
mu_vec = u:-1:1;                     % μ_i 降序：3,2,1
q_vec  = linspace(0.8,0.2,v);        % q_j 降序：0.8→0.2
% -------------------------------------------------------------------------

assert(length(alpha_target)==u && abs(sum(alpha_target)-1)<1e-12,...
       'alpha_target 长度 / 和 数错误');
assert(length(beta_target)==v  && abs(sum(beta_target) -1)<1e-12,...
       'beta_target  长度 / 和 数错误');

n_list = [];  ratio_r = [];  ratio_b = [];

for n = 2:n_max
    % ===== 1. Weak compositions =====
    K = compositions(n,u);     Nu = size(K,1);      % Nu × u
    L = compositions(n,v);     Nv = size(L,1);      % Nv × v
    
    % ===== 2. Bernstein 基方阵 =====
    T_alpha = buildT(K,n);                       % Nu × Nu  sparse
    T_beta  = buildT(L,n);                       % Nv × Nv  sparse
    
    % ===== 3. 计算 Z^r / Z_B^r 整表 =====
    Zr  = zeros(Nu,Nv);   ZBr = zeros(Nu,Nv);
    for i = 1:Nu
        alpha_vec = K(i,:)./n;
        for j = 1:Nv
            beta_vec = L(j,:)./n;
            Eq = @(p) sum(beta_vec .* (q_vec.^p));
            % ---- S (u×u) ----
            S = zeros(u,u);
            for a = 1:u
                S(a,a) = alpha_vec(a)*Eq(2*mu_vec(a));
                for b = a+1:u
                    val = sqrt(alpha_vec(a)*alpha_vec(b))*Eq(mu_vec(a)+mu_vec(b));
                    S(a,b)=val; S(b,a)=val;
                end
            end
            eigvals = eig(S);
            lambda_max = max(eigvals);  lambda_min = min(eigvals);
            r = lambda_min / lambda_max;
            % ---- Z^r ----
            if abs(1-r)<1e-12
                Zr(i,j) = lambda_max^n * (n+1);
            else
                Zr(i,j) = lambda_max^n * (1-r^(n+1))/(1-r);
            end
            % ---- Z_B^r ----
            zb = 0;
            for m2 = 0:n
                term_m = (0.5*(1+r))^m2 / factorial(m2);
                inner  = 0;
                for k2 = 0:n-m2
                    inner = inner + nchoosek(2*k2,k2)* ...
                                   nchoosek(2*(n-m2-k2),n-m2-k2)/4^(n-m2) * r^k2;
                end
                zb = zb + term_m*inner;
            end
            ZBr(i,j) = lambda_max^n * zb;
        end
    end
    
    % ===== 4. 直接左右求解  Tα·Z·Tβᵀ = Zᵣ  =====
    Zb  = (T_alpha \ Zr ) / T_beta.';    % Z^b
    ZBb = (T_alpha \ ZBr) / T_beta.';    % Z_B^b
    
    % ===== 5. 找到与目标 α⃗,β⃗ 对应的格子 =====
    k_target = round(alpha_target * n);
    k_target(end) = n - sum(k_target(1:end-1));
    l_target = round(beta_target  * n);
    l_target(end) = n - sum(l_target(1:end-1));
    
    idxK = find(all(K==k_target,2),1);
    idxL = find(all(L==l_target,2),1);
    if isempty(idxK) || isempty(idxL)
        error('n=%d 时找不到对应组合，可增大 n 或调整 target。',n);
    end
    
    % ===== 6. 记录比值 =====
    ratio_r(end+1) = abs( ZBr(idxK,idxL) / Zr(idxK,idxL) );
    ratio_b(end+1) = abs( ZBb(idxK,idxL) / Zb(idxK,idxL) );
    n_list(end+1)  = n;
    
    fprintf('[n=%2d]  k=%s  l=%s  |ZB^r/Z^r|=%.3e  |ZB^b/Z^b|=%.3e\n',...
            n, mat2str(k_target), mat2str(l_target), ...
            ratio_r(end), ratio_b(end));
end

% ===== 7. 画图 =====
figure;
loglog(n_list, ratio_r,'b-','LineWidth',2); hold on;
loglog(n_list, ratio_b,'r--','LineWidth',2);
loglog(n_list, 1./sqrt(n_list),'k-','LineWidth',2);
xlabel('$n$','Interpreter','latex');
ylabel(sprintf('Ratios at $\\alpha=%s$, $\\beta=%s$', ...
       mat2str(alpha_target), mat2str(beta_target)), 'Interpreter','latex');
legend({'$|Z_{\mathrm B}^r/Z^r|$', '$|Z_{\mathrm B}^b/Z^b|$', '$1/\\sqrt{n}$'}, ...
       'Interpreter','latex','Location','southwest'); grid on;

% ---------- save as PDF ----------
w_in = 6; h_in = 4.5;
set(gcf,'Units','inches','Position',[1 1 w_in h_in]);
set(gcf,'PaperUnits','inches','PaperSize',[w_in h_in],'PaperPosition',[0 0 w_in h_in]);

if ~exist('results','dir'), mkdir('results'); end
vec2fname = @(v) strrep(strjoin(arrayfun(@(x)sprintf('%.2f',x),v,'UniformOutput',false),'_'),'.','p');
fname = sprintf('results/Z_ratio_alpha_%s_beta_%s.pdf', vec2fname(alpha_target), vec2fname(beta_target));
exportgraphics(gcf, fname, 'ContentType','vector', 'Resolution',300);
fprintf('[Saved] %s\n', fname);

% ============================================================
%                       helper functions
% ============================================================
function K = compositions(n,u)
% All weak compositions of n into u parts
    if u==1, K = n; return; end
    C = nchoosek(1:n+u-1,u-1);  % bar positions
    m = size(C,1);   K = zeros(m,u);
    for idx = 1:m
        cuts     = [0, C(idx,:), n+u];
        K(idx,:) = diff(cuts)-1;
    end
end

function T = buildT(K,n)
% Sparse multi-index Bernstein-Vandermonde matrix
    N   = size(K,1);           mf = factorial(n);   F = factorial(0:n);
    rows = []; cols = []; vals = [];
    for r = 1:N
        a = K(r,:)./n;                         % sample point
        for c = 1:N
            k = K(c,:);
            val = mf/prod(F(k+1))*prod(a.^k);  % n!/k! * a^k
            if val
                rows(end+1)=r; cols(end+1)=c; vals(end+1)=val;
            end
        end
    end
    T = sparse(rows,cols,vals,N,N);
end
