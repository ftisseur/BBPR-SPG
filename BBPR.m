function [nu,info] = BBPR(opts)
%BBRMA  BB-step projected residual algorithm for optimal damping.
%  nu = BBRMA(opts) solves an optimal damping problem with nonnegative constraints
%  It finds the vector nu that minimizes the total energy of the system with
%  corresponding quadratic l^2*M+l*(Dint+D*diag(nu)*D')+K with D of size nxk.
%  opts is a structure with fields
%  Phi: nxn modal basis of K-l*M such that Phi'*M*Phi = I_n.
%  Om: nx1 vector containing the square root of the eigenvalues of K-l*M
%  Gamma: nx1 such that Gamma = diag(Phi'*Dint*Phi)
%  PhiD: nxk matrix Phi'*D
%  s: number of eigenmodes to dampen. Default is n
%  nu0: initial viscosities
%  tolres: tolerence for residual function h(\nu).
%  tolnu: tolerance for accumulation point nu.
%  f: set to true to evaluate the objective function.



%  Based on the BBPR algorithm described in
%  "Fast Algorithms for Optimal Damping in Mechanical Systems
%  by Q Li, F Tisseur, arXiv.2601.05404
%  last updated by FT on 04/08/2026

%Set up the parameters.
eta = 2;
iter_max = 1000;
eta_min = 1.0e-5;
eta_max = 1.0e10;
rng(1) % add setting random seed, added by Qingna Li, 5 August 2025

[tmprY,tmprW,sum1,W0hat,LeftY0] = getYW0(opts);
nu_prev = opts.nu0;
nu = nu_prev;
j = 0;
res = zeros(1,iter_max);
while j<iter_max
   [lY1,rY1,lW1,rW1,TinvU,TadjU,Lam,Tinv]  = getDWYEA_v3(nu,opts,tmprY,tmprW);
   [h] = getHDelta0(nu,lY1,rY1,lW1,rW1,Lam,sum1,TadjU,TinvU);
   nh = norm(h);
   res(j+1)=nh;
   if j>2 && (nh < opts.tolres && res_nu < opts.tolnu)
      break;
   end
   if j>0
      deltah = h - h0;
      denominator = deltanu'*deltah;
      if ~isfinite(denominator) || denominator <= 0
         etaBB = eta_max;
      else
         etaBB = (deltanu'*deltanu)/denominator;
         if ~isfinite(etaBB) || etaBB <= 0; etaBB = eta_max; end
      end
      eta = max(eta_min,min(eta_max,etaBB));
   end
   nu = max(0,real(nu_prev-eta*h));
   deltanu = nu - nu_prev;
   norm_delta = norm(deltanu);
   norm_prev = norm(nu_prev);
   if norm_prev > 0
      res_nu = norm_delta/norm_prev;
   elseif norm_delta == 0
      res_nu = 0;
   else
      res_nu = Inf;
   end
   h0 = h;
   nu_prev = nu;
   j = j+1;
end

if j~=iter_max && isfield(opts,'f') && opts.f %Evaluate the objective function
   [n,k] = size(opts.PhiD);
   s = opts.s;
   tmptr = W0hat*LeftY0;
   trace1w =real( sum(sum((conj(LeftY0)).*tmptr)));
   DWhat = lW1*rW1';
   DWhat = DWhat+ DWhat';
   DWhat = DWhat./(conj(Lam)+Lam.');%
   %% calculate the second part of objective function: trace(ZW)
   LeftY = zeros(2*n,2*s);
   LeftY(:,1:s) = Tinv(:,1:s);
   LeftY(:,s+1:2*s)= Tinv(:,n+1:n+s);
   tmptr = DWhat*LeftY;
   trace2w =real( sum(sum((conj(LeftY)).*tmptr)));
   trace2w = trace2w/(2*s);
   trss = trace1w + trace2w;
   info.f = trss;
else
   info.f = [];
end

res = res(1:min(j+1,iter_max));
info.reshist = res;
info.res = res(end);
info.iter = j;
info.it_eig=numel(res);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tmprY,tmprW,sum1,W0hat,LeftY0] = getYW0(opts)
Gamma = opts.gamma;
Omega = opts.omega;
[n,k] = size(opts.PhiD);
s = opts.s;
D = opts.PhiD;
vargamma = sqrt(Gamma.*Gamma-4.0.*Omega.*Omega);

sum1 = Gamma + vargamma;
diff1 = Gamma - vargamma;
lambda0 = -0.5 *[sum1; diff1];

t11 = -0.5*(diff1./Omega);
t12 = -0.5 .*(sum1 ./Omega);
% T0 = [diag(t11) diag(t12); diag(t21) diag(t21)];

a11 = Omega./vargamma;
a12 = 0.5 *(sum1 ./vargamma);
a22 = -0.5 *(diff1./vargamma);
T0inv = [diag(a11) diag(a12);
         diag(-a11) diag(a22) ];
LeftY0(:,1:s) = T0inv(:,1:s);
LeftY0(:,s+1:2*s)= T0inv(:,n+1:n+s);

LeftY0 = 1/sqrt(2*s).* LeftY0;
Y0hat = -LeftY0*LeftY0';
Y0hat = Y0hat./(lambda0+lambda0');

b11 = conj(t11).*t11 + 1;
b12 = conj(t11).*t12 + 1;
b21 = conj(t12).*t11 + 1;
b22 = conj(t12).*t12 + 1;
W0hat = -[diag(b11) diag(b12); diag(b21) diag(b22)];
W0hat = W0hat./(conj(lambda0)+lambda0.');
%W0hat = W0hat./(repmat(conj(lambda0),1,2*n)+repmat(lambda0.',2*n,1));
Y0hat = (Y0hat+Y0hat')*0.5;
W0hat = (W0hat+W0hat')*0.5;
tmp1 = [D;D]; % T0'*[zeros(size(D));D]
tmp2 =[ a12.*D;a22.*D]; % T0inv*[zeros(size(D));D]
left1 = Y0hat*tmp1;
right1 = W0hat'*tmp2;

sum1 = sum(conj(left1).*right1,1).';
tmprY = [t11.*left1(1:n,:)+t12.*left1(n+1:2*n,:);
         left1(1:n,:)+left1(n+1:2*n,:)];
tmprW = [conj(a11).*right1(1:n,:)-conj(a11).*right1(n+1:2*n,:);
        conj(a12).*right1(1:n,:)+conj(a22).*right1(n+1:2*n,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lY1,rY1,lW1,rW1,TinvU,TadjU,Lam,Tinv] = getDWYEA_v3(nu,opts,tmprY,tmprW)
idx = find(abs(nu)>=0.1*opts.tolres);
PhiD = opts.PhiD;
[n,k] = size(opts.PhiD);
R = PhiD(:,idx)*diag(sqrt(nu(idx)));
[T,Lam,Tinv] = EA_eig(opts.gamma,R,opts.omega);
TinvU = Tinv(:,n+1:2*n) * opts.PhiD;
TadjU = T(n+1:2*n,:)' * opts.PhiD;
lY1 = TinvU*diag(nu);
lW1 = TadjU*diag(nu);
rY1 = Tinv*tmprY;  % T^{-1}Y0^H(U^T)^H
rW1 = T'*tmprW;     % T^{H}W0^H(U^T)^H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = getHDelta0(nu,lY1,rY1,lW1,rW1,Lam,sum1,TadjU,TinvU)
k = size(lY1,2);
d = zeros(k,1);
py = apply_lyapunov_correction(Lam,lY1,rY1,TadjU);
pw = apply_lyapunov_correction(conj(Lam),lW1,rW1,TinvU);

for i = 1:k
   sum2 = py(:,i)'*pw(:,i);
sum3 = rY1(:,i)'*pw(:,i);
sum4 = py(:,i)'*rW1(:,i);
d(i) = real(sum1(i) + sum2 + sum3 + sum4);
d(i) = real(d(i) + d(i)');
end
h = nu - max(nu+d,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
