function [T,Lam,Tinv] = EA_eig(Dint,R1,Om) 
%EA_EIG Consolidation of calleig_EA_qn, eig_EA_qn and loreqes.
% Requires SEA.m. Preserves the original arithmetic and random-number calls.
 [VR,Lam,VL,DVR] = compute_eigenpairs(Dint,R1,Om);
       T = [Om.*VR; VR.*Lam.'];% same as  X = [diag(Om)*VR;VR*diag(Lam)]
         %     Tinv = [(diag(Lam)*VL'+VL'*(diag(Dint)+R1*R1'))*(diag(1./Om)) VL']; 
       s = Lam.*VL';      
       ss = VL'.*Dint.';
       sss = VL'*R1;
       sss4 = sss*R1';
       mm = s+ss+sss4;
       tmp = 1./Om;
       tmp0 = mm.*tmp.';
       Tinv = [tmp0 VL'];
end

function [VR,Lam,VL,DVR] = compute_eigenpairs(Dint,R,Om)
% Computes right eigenvectors VR, eigenvalues Lam, and left eigenvectors VL of
% quadratic matrix polynomial l^2 I_n + l(Dint+RR^T)+Om^2 with Dint and Om
% diagonal, and R is n x kd with kd<<n.

% This is a modification of Leo's code to include internal damping.
% Must be applied to diagonal quadratics + low rank.

[n,r] = size(R);
locked = zeros(2*n,1);
warning('off','MATLAB:nearlySingularMatrix') % We expect nearly singular matrices
warning('off','MATLAB:singularMatrix')       % or even singular

M1 = ones(n,1);
K1 = Om.^2;
[Lam,VR] = low_rank_quadratic(M1,Dint,R,K1);
%berr = bwerr(diag(M1),D,diag(K1),V,Lam)
VL = conj(VR); %left e'vecs of Q: nx2n
%berrVL = bwerr(diag(M1),diag(Dint)+R*R.',diag(K1),VL,conj(z0)) %this gives small be
%need some form of normalization
%D = diag(Dint)+R*R';
D1 = Dint.*VR;
D2 = R'*VR;
D2 = R*D2;
DVR = D1+D2;
VLDVR = sum(conj(VL).*DVR,1);
VLDVR = transpose(VLDVR);
VLT=conj(VL);
 stmp = sum(VLT.*VR,1);
 s=stmp.'; 
 s = 2*(s.*Lam);
 ctest = s+ VLDVR;

for i=1:2*n %This can be made more efficiently  
    c=ctest(i); %  c = VL(:,i)'*(2*Lam(i)*diag(M1)+D)*VR(:,i); improved by Qingna
    VL(:,i) = (1/conj(c))*VL(:,i);   
%    VL(:,i)'*(2*E(i)*diag(M1)+D)*VR(:,i)  %should be 1
end
end

function [z,V] = low_rank_quadratic(M1,D1,R,K1)
% Modified from Leo's code to include internal damping. Exclusively related to
% our problem

n = length(M1);
r = size(R,2);
warning('off','MATLAB:nearlySingularMatrix') % We expect nearly singular matrices
warning('off','MATLAB:singularMatrix')       % or even singular

normM = max(abs(M1));
normK = max(abs(K1));
D = R*R.';
if n>100
   normD = normest(diag(D1)+D,1e-3);
else
   normD = norm(diag(D1)+D);
end

% Construct initial eigenvalues
z = zeros(2*n,1)+1i*zeros(2*n,1);
signD1 = sign(D1); %this Matlab function does'nt work well, next line needed.
if norm(signD1) == 0, signD1 =ones(n,1);end
z(1:n,1) = .5*(-D1+signD1.*sqrt(D1.^2-4*K1));
z(n+1:end,1) = K1./z(1:n,1);

%Initial locking of eigenvalues
locked = zeros(2*n,1);
ae = abs(z(1:n,1))'; %only pick up the first n eigenvalues
%f = sqrt(ae) ./ (normM*ae + normD*sqrt(ae) + normK);
f = ae ./ (normM*(ae.^2) + normD*ae + normK);
J = find(sqrt(sum(D.*D)) .* f < n*eps);
locked([J J+n])=1;
%[(sqrt(sum(D.*D)) .* f)' locked(1:n)], pause
%sumlocked = sum(locked)

%% Perturb unlocked starting points
nlocked = find(~locked);
z(nlocked) = z(nlocked) - 1e-8*abs(z(nlocked)).*(rand(length(nlocked),1)+1i*rand(length(nlocked),1));

%% Find eigenvalues using special Ehrlich-Aberth iteration
[z,averageit] = SEA(R,M1,D1,K1,z,locked);
%fprintf('On avarage the eigenvalues required %0.1f iterations to converge.\n',averageit);

%% Retrieve eigenvectors
warning('off')
V = [eye(n) eye(n)];
I = eye(r);
for i = 1:2*n
    if locked(i)
        continue
    end
    A = M1*(z(i)^2) + D1*z(i) + K1;
    amax = max(abs(A));
    smallA = find(abs(A)<amax*eps);
    A(smallA) = amax*(rand+1i*rand)*eps; 
    AinvR = A.\R;% AinvR = bsxfun(@ldivide,A,R);
    cap = (R.'*AinvR)*z(i) + I;
    [~,s,v] = svd(cap);
%     if isnan(v)
%        ss = s
%        vv = v
%     end
    V(:,i) = AinvR*v(:,end);
    if s(1,1)/s(r,r)<1/eps && s(1,1)>realmin % do one iteration of inverse iteration
        x = conj(V(:,i));
        % Sherman-Morrison-Woodbury O(nr^2)      
        x = A.\x; %x = bsxfun(@ldivide,A,x);   
        V(:,i) = x - z(i)*(AinvR*(cap\(R.'*x)));
    end
end
V = normalize(V,'norm');
warning('on')

end
