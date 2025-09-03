function p = perm_ryser_gray(A)
% Permanent via Ryser's formula with Gray code (exact)
% perm(A) = sum_{S != empty} (-1)^(n - |S|) prod_i sum_{j in S} a_{ij}
    [N, ~] = size(A);
    cols = cell(1,N);
    for j = 1:N, cols{j} = A(:,j); end

    rowSums = zeros(N,1);
    p = 0.0;
    prevGray = uint64(0);
    subsetSize = 0;                   % |S|
    totalStates = bitshift(uint64(1), N);  % 2^N

    for g = uint64(1) : totalStates-1  % skip empty set (product=0)
        gray = bitxor(g, bitshift(g, -1));
        diff = bitxor(gray, prevGray);
        % find flipped bit index j (1..N)
        j = find(bitget(diff, 1:N), 1, 'first');
        if bitget(gray, j)
            rowSums = rowSums + cols{j};  subsetSize = subsetSize + 1;
        else
            rowSums = rowSums - cols{j};  subsetSize = subsetSize - 1;
        end
        prevGray = gray;

        % sign = (-1)^(N - |S|)
        sgn = 1.0; if mod(N - subsetSize, 2) ~= 0, sgn = -1.0; end

        % product_i rowSums(i)
        prodTerm = 1.0;
        for i = 1:N
            si = rowSums(i);
            if si == 0, prodTerm = 0.0; break; end
            prodTerm = prodTerm * si;
        end
        p = p + sgn * prodTerm;
    end
end
