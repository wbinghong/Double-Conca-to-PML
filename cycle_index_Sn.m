function Z = cycle_index_Sn(n)
    % Use symbolic to calculate the cycle index of symmetry group n
    syms t
    z = sym('z', [1 n]);   % 定义 z1, z2, ..., zn
    assume(t, 'real');
    
    % egf
    expr = 0;
    for k = 1:n
        expr = expr + (z(k) * t^k) / k;
    end
    C = exp(expr);

    % n th deri at t = 0
    dC = diff(C, t, n);
    Z_sym = subs(dC, t, 0) / factorial(n);

    % output
    Z = expand(Z_sym);
end
