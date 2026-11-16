function [nu,info] = SPG(opts)
%SPG  Spectral Projection Gradient method for optimal damping.
%  nu = SPG(opts) solves an optimal damping with nonnegative constraints,
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

doprint = 0;
% initialization
eta_min = 1.0e-5;
eta_max = 1.0e10;
eta_fallback = eta_max; %1.0e5;
j = 0;
iter_max = 1000;
tol = opts.tolres;
tol_nu = opts.tolnu;
nu0 = opts.nu0;
eta = 2;
[n,k] = size(opts.PhiD);
dhat = zeros(k,1);
ss = zeros(k,1);
%% extra for line search
M_ls = 5; %10;
sigma = 1.0e-4;
ratiorho = 0.5;
res = zeros(1,iter_max);
f_hist = zeros(1,iter_max);
tracef = zeros(1,iter_max);
nuhist = zeros(k,iter_max);
iter_eig = 0; %% number of eigenvalue decompositions
iter_ls = 0; %% number of iterations that call line search

rng(1) % add setting random seed, added by Qingna Li, 5 August 2025
[tmprY,tmprW,sum1,trace1w] = getYW0_v5(opts);
nu = nu0;
[trss,T,Tinv,Lam,TadjU,DWhat,rW] = f_w(nu,opts,tmprW,trace1w);
[h,ep] = grad_w(nu,opts,T,Tinv,Lam,TadjU,DWhat,rW,tmprY,tmprW,sum1);
nh = norm(h);
ep1 = ep;
if doprint
  fprintf('-----------Start the Main Loop-------------------\n')
  fprintf('--iter:        res        eta      nu                          trace     t \n')
end
iter_eig = iter_eig + 1;
while j<iter_max
  res(j+1) = nh;
  f_hist(j+1) = trss;
  if j>2&& (nh < tol  &&res_nu< tol_nu)
     if doprint
       fprintf('--%d-th:  %3.2e --  %5.5f\n',j,res(j+1),trss)
       fprintf('successfully stopped\n')
     end
     break;
  end
  if j>0
     denominator = deltanu'*deltag;
     if ~isfinite(denominator) || denominator <= 0
        etaBBg = eta_fallback;
     else
        etaBBg = (deltanu'*deltanu)/denominator;
        if ~isfinite(etaBBg) || etaBBg <= 0; etaBBg = eta_fallback; end
     end
     eta =max(eta_min,min(eta_max,etaBBg));
  end
  for i = 1:k
      dhat(i) = max(0,real(nu(i)-eta*ep(i)))-nu(i) ; %added for modify BB on 05 December 2024
      ss(i) = nu(i)+ dhat(i);  %added for modify BB on 05 December 2024
  end
  descent_a=ep'*dhat;%,pause
  [trss,Tnew,Tinvnew,Lamnew,TadjUnew,DWhatnew,rWnew] = f_w(ss,opts,tmprW,trace1w);
  iter_eig = iter_eig + 1;
  if j<M_ls
     fmax=max(f_hist(1:j+1));
  else
     fmax = max((f_hist(j-M_ls+2:j+1)));
  end
  fdif = real(trss-fmax);
  rightdes=real(descent_a*sigma);
  difsearch = fdif-rightdes;%,pause
  gamma = 1.0;
  j_ls = 0;
  while difsearch>0.0 %&& linesearch_flag==1
     gamma = gamma*ratiorho;
     j_ls=j_ls+1;
     ss= nu+ gamma*dhat;
     [trss,Tnew,Tinvnew,Lamnew,TadjUnew,DWhatnew,rWnew] = f_w(ss,opts,tmprW,trace1w);
     rightdes=rightdes*ratiorho;
     fdif = trss-fmax;
     difsearch = fdif-rightdes;
     if j_ls >= 11 && difsearch > 0
        error('SPG:LineSearchFailure', ...
              'Line-search inequality still fails after 11 step reductions.');
     end
  end
  if j_ls>0; iter_ls = iter_ls+1;end
  [hnew,ep] = grad_w(ss,opts,Tnew,Tinvnew,Lamnew,TadjUnew,DWhatnew,rWnew,tmprY,tmprW,sum1);
  iter_eig= iter_eig + j_ls;
  h = hnew;
  nh = norm(h);
  nu0 = nu;
  nu = ss;
  nuhist(:,j+1) = nu(:);
  deltanu = nu-nu0;
  norm_delta = norm(deltanu);
  norm_prev = norm(nu0);
  if norm_prev > 0
     res_nu = norm_delta/norm_prev;
  elseif norm_delta == 0
     res_nu = 0;
  else
     res_nu = Inf;
  end
  deltag = ep-ep1;
  ep1 = ep;
  if doprint
     fprintf('--%d-th: %d    %3.2e %5.5f %5.5f\n ',j,j_ls,res(j+1),eta,trss);
  end
  tracef(j+1) = trss;
  j = j+1;
end
nres = min(j+1,iter_max);
res = res(1:nres);
f_hist = f_hist(1:nres);
tracef = tracef(1:j);
nuhist = nuhist(:,1:j);
info.reshist = res;
info.res = res(end);
info.iter = j;
info.it_eig=iter_eig;
info.trace = trss;
info.nu = nu;
info.nuhist = nuhist;
info.tracef = tracef;
info.fhist = f_hist;
info.iter_ls = iter_ls;%% number of iterations that call line search strategy;


%-------
function [tmprY,tmprW,sum1,trace1w] = getYW0_v5(opts)
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
%T0 = [diag(t11) diag(t12); diag(t21) diag(t21)];
a11 = Omega./vargamma;
a12 = 0.5 *(sum1 ./vargamma);
a22 = -0.5 *(diff1./vargamma);
T0inv = [diag(a11) diag(a12);diag(-a11) diag(a22) ];
LeftY0(:,1:s) = T0inv(:,1:s);
LeftY0(:,s+1:2*s)= T0inv(:,n+1:n+s);
LeftY0 = 1/sqrt(2*s).* LeftY0;
Y0hat = -LeftY0*LeftY0';
Y0hat = Y0hat./(lambda0+lambda0');
b11 = conj(t11).*t11 + 1;
b12 = conj(t11).*t12 + 1;
b21 = conj(t12).*t11 + 1;
b22 = conj(t12).*t12 + 1;
% W0hat = -T0adj*T0adj';
W0hat = -[diag(b11) diag(b12); diag(b21) diag(b22)];
%W0hat = W0hat./(repmat(conj(lambda0),1,2*n)+repmat(lambda0.',2*n,1));
W0hat = W0hat./(conj(lambda0)+lambda0.');
Y0hat = (Y0hat+Y0hat')*0.5;
W0hat = (W0hat+W0hat')*0.5;
tmp1 = [D;D];      % T0'*[zeros(size(D));D]
tmp2 =[ a12.*D;a22.*D];        % T0inv*[zeros(size(D));D]
left1 = Y0hat*tmp1;
right1 = W0hat'*tmp2;
sum1=zeros(k,1);
for i = 1 : k % # can be improved by Cauchylike*Matrix
   sum1(i)= left1(:,i)'*right1(:,i);
end
%tmprY = T0*left1;
tmprY = [t11.*left1(1:n,:)+t12.*left1(n+1:2*n,:);left1(1:n,:)+left1(n+1:2*n,:)];
% tmprW = T0inv'*right1;
tmprW = [conj(a11).*right1(1:n,:)-conj(a11).*right1(n+1:2*n,:); conj(a12).*right1(1:n,:)+conj(a22).*right1(n+1:2*n,:)];
%% add first part in function value, that is, trace(ZW0)
tmptr = W0hat*LeftY0;
trace1w =real( sum(sum((conj(LeftY0)).*tmptr)));
% trace1wtest = trace(LeftY0'*tmptr)
% trace1w - trace1wtest,pause

function  [trss,T,Tinv,Lam,TadjU,DWhat,rW1] = f_w(nu,opts,tmprW,trace1w)
%what does it do?
[n,k] = size(opts.PhiD);
idx = find(abs(nu)>=0.1*opts.tolres);
PhiD = opts.PhiD;
R = PhiD(:,idx)*diag(sqrt(nu(idx)));
[T,Lam,Tinv] = EA_eig(opts.gamma,R,opts.omega);
TadjU = T(n+1:2*n,:)' * opts.PhiD;
lW1 = TadjU*diag(nu);
rW1 = T'*tmprW;     % T^{H}W0^H(U^T)^H
DWhat = lW1*rW1';
DWhat = DWhat+ DWhat';
DWhat = DWhat./(conj(Lam)+Lam.');%
s = opts.s;

% calculate the second part of objective function: trace(ZW)
LeftY = zeros(2*n,2*s);
LeftY(:,1:s) = Tinv(:,1:s);
LeftY(:,s+1:2*s)= Tinv(:,n+1:n+s);
tmptr = DWhat*LeftY;
trace2w =real( sum(sum((conj(LeftY)).*tmptr)));
trace2w = trace2w/(2*s);
trss = trace1w + trace2w;



function  [h,d_ep] = grad_w(nu,opts,T,Tinv,Lam,TadjU,DWhat,rW1,tmprY,tmprW,sum1);
[n,k] = size(opts.PhiD);
d = zeros(k,1);
TinvU = Tinv(:,n+1:2*n) * opts.PhiD;
lY1 = TinvU*diag(nu);
rY1 = Tinv*tmprY;  % T^{-1}Y0^H(U^T)^H
pw = DWhat*TinvU;
py = apply_lyapunov_correction(Lam,lY1,rY1,TadjU);
for i = 1:k
    sum2 = py(:,i)'*pw(:,i);
    sum3 = rY1(:,i)'*pw(:,i);
    sum4 = py(:,i)'*rW1(:,i);
    d(i) = real(sum1(i) + sum2 + sum3 + sum4);
    d(i) = real(d(i) + d(i)');
end
d_ep = -d;
h = nu - max(nu+d,0);
