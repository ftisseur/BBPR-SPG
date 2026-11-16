function product = apply_lyapunov_correction(Lam,L,R,V)
%APPLY_LYAPUNOV_CORRECTION Multiply a Lyapunov correction by V.
%   product = H*V, where
%   diag(Lam)*H + H*diag(conj(Lam)) = L*R' + R*L'.
%   Row blocks avoid forming the full dense H for large problems.

n = size(L,1);
block_size = 256;
if n <= block_size
    H = L*R';
    H = H + H';
    H = H./(Lam+Lam');
    product = H*V;
    return
end

product = zeros(n,size(V,2),'like',V);
for first = 1:block_size:n
    rows = first:min(first+block_size-1,n);
    Hblock = L(rows,:)*R';
    Hblock = Hblock + (L*R(rows,:)')';
    Hblock = Hblock./(Lam(rows)+Lam');
    product(rows,:) = Hblock*V;
end
end
