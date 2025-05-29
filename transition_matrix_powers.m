clear; clc;

syms X11 X22 X12 real

alpha = 0.5;

% Entries
C11 = alpha * X11;
C22 = (1 - alpha) * X22;
C12 = sqrt(alpha * (1 - alpha)) * X12;

% S
S = [C11, C12;
     C12, C22];
disp('S =');
disp(S);

% interested k (which is not have to be even in our step calculation)
k_list = [2, 3, 4];

% S^k and trace(S^k)
for k = k_list
    Sk = simplify(S^k);
    tr_Sk = simplify(trace(Sk));
    
    % save to workzone on right
    assignin('base', sprintf('S%d', k), Sk);
    assignin('base', sprintf('tr_S%d', k), tr_Sk);
    
    % fprintf('S^%d =\n', k);
    % disp(Sk);
    % 
    % fprintf('trace(S^%d) =\n', k);
    % disp(tr_Sk);
end