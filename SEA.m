function [z,averageit] = SEA(R,M1,D1,K1,z,locked)
%function [z,averageit] = SEA(R,M1,D1,K1,z,locked)
% SEA: Ehrlich-Aberth to compute the eigenvalues of
% diag(M1)*lambda^2 +(D1+RR^T)*lambda + diag(K1)
% R: n x r matrix  describing external low-rank damping
% M1: n x 1 vector describing diagonal mass matrix
% D1: n x 1 vector describing diagonal internal damping
% K1: n x 1 vector describing diagonal stiffness matrix
% z      : 2n x 1 vector (starting points)
% locked : 2n x 1 binary vector


[n,r] = size(R);
N = 2*n;
I = eye(r);

warning('off','MATLAB:nearlySingularMatrix') % We expect nearly singular matrices
warning('off','MATLAB:singularMatrix')       % or even singular

maxIt = 100*N;
counter = 0;
nlocked = sum(locked);
epsi = eps;
it = 0;
while 1
    it = it+1;
    if mod(it,50)==0
        epsi = epsi*10; %To relax stoppoing criterion. Recommended by L. Taslaman
    end
    for k=1:N
        if locked(k), continue, end
        counter = counter+1;
        zk = z(k);
        %pdp = 1/lowrankupd(zk)

        A = M1*(zk*zk) + D1*zk + K1;   % n-by-1 diagonal vector
        si = find(~abs(A)); A(si,1)=eps;
        AR = R ./ A; % Columnwise division
        RHAR = R.' * AR; % Nonconjugate transpose
        idplus = eye(r) + zk * RHAR;
        M1AinvR = M1 .* AR;
        D1AinvR = D1 .* AR;
        B1 = AR.' * M1AinvR;
        B3 = AR.' * D1AinvR;
 %To remove
        if any(~isfinite(A))
           fprintf('Nonfinite A: k=%d, it=%d, zk=%g%+gi\n', ...
            k,it,real(zk),imag(zk));
           %keyboard
        end

        if any(~isfinite(idplus(:)))
           fprintf('Nonfinite idplus: k=%d, it=%d, zk=%g%+gi\n', ...
           k,it,real(zk),imag(zk));
           %keyboard
        end
 %up to here
        temp1 = (idplus.' \ B1.').';
        temp2 = RHAR * (idplus \ RHAR);
        temp3 = (idplus.' \ B3.').';
        stemp3 = warning('query','MATLAB:nearlySingularMatrix');
        ppdp = 2*zk*sum(M1 ./ A) + trace(RHAR) - 2*zk^2*trace(temp1) ...
            - zk*trace(temp2) - zk*trace(temp3) + sum(D1 ./ A);
        pdp = 1/ppdp;

        if isinf(pdp)
            display('entered isinf(pdp)')
            pdp = abs(zk)*(rand+1i*rand)*1e-3;
        end
        summation = sum(1./(zk-z([1:k-1 k+1:N]))); % O(2n)
        upd = pdp/(1-pdp*summation);
        if abs(upd)<epsi*abs(zk) || isnan(upd)
            locked(k) = 1;
            nlocked = nlocked+1;
            if isnan(upd), continue, end
        end
        z(k) = zk - upd;
    end
    if nlocked==N || counter == maxIt
        break
    end
end
if (nlocked<N)
    warning('Convergence not obtained after %d iters for %d evals\n',maxIt,N-nlocked)
end
averageit = counter/N;
warning('on')

%---------

% function dpp = lowrankupd(x) % total flop count: O(4n+2rn+4r^2n)
%
% A = M1*(x*x) + D1*x + K1;   % n-by-1 diagonal vector
% AR = R ./ A; % Columnwise division
% RHAR = R.' * AR; % Nonconjugate transpose
% idplus = eye(r) + x * RHAR;
% M1AinvR = M1 .* AR;
% D1AinvR = D1 .* AR;
% B1 = AR.' * M1AinvR;
% B3 = AR.' * D1AinvR;
% temp1 = (idplus.' \ B1.').';
% temp2 = RHAR * (idplus \ RHAR);
% temp3 = (idplus.' \ B3.').';
% dpp = 2*x*sum(M1 ./ A) + trace(RHAR) - 2*x^2*trace(temp1) ...
%     - x*trace(temp2) - x*trace(temp3) + sum(D1 ./ A);
% end

end
