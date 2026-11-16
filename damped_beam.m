function [M,D,K] = damped_beam(n,opts)
%DAMPED_BEAM  QEP from simply supported beam damped in the middle.
%  [M,D,K] = damped_beam(N) constructs an N-by-N
%  quadratic matrix polynomial lambda^2*M + lambda*D + K from a finite
%  element model of a beam clamped at both ends with a damper in the middle.
%  N/2 is the number of finite elements.
%  N should be even. Otherwise N-1 is used. Default: N=200.
%  Half of the eigenvalues of the problem are pure imaginary
%  and are eigenvalues of the undamped problem (D = 0).

%  Reference:
%  N. J. Higham, D. S. Mackey, F. Tisseur and S. D. Garvey.
%  Scaling, sensitivity and stability in the numerical solution of
%  quadratic eigenvalue problems. Internat. J. Numer. Methods Eng.,
%  73(3):344-360, 2008.


% Geometric and material properties of the beam.
if nargin < 2
   width = 0.05;
   height = 0.005;
   glength = 1; %set back to 1
else
   width = opts.width;
   height = opts.height;
   glength = opts.length;
end
nele = floor(n/2);
n = 2*nele;
dlength = glength/nele;
E = 7e10;
I = width*height^3/12;
area = width*height;
rho =  0.674/(area*glength);

damp = 5;

K1 = [12 6*dlength; 6*dlength 4*dlength^2];
K2 = [-12 6*dlength; -6*dlength 2*dlength^2];
K3 = [12 -6*dlength; -6*dlength 4*dlength^2];

M1 = [156 22*dlength; 22*dlength 4*dlength^2];
M2 = [54 -13*dlength ;13*dlength -3*dlength^2];
M3 = [156 -22*dlength; -22*dlength 4*dlength^2];

% u = ones(nele,1);
v1 = [ones(nele,1); 0];
v2 = [0; ones(nele,1)];
Du = spdiags(ones(nele+1,1),1,nele+1,nele+1);
Dv1 = spdiags(v1,0,nele+1,nele+1);
Dv2 = spdiags(v2,0,nele+1,nele+1);

% K = kron(diag(v1),K1)+kron(diag(v2),K3)+kron(diag(u,1),K2)+kron(diag(u,-1),K2');
K = kron(Dv1,K1) + kron(Dv2,K3) + kron(Du,K2) + kron(Du',K2');
ind = [(2:n) n+2];  % Takes care of boundary conditions.
K = (E*I/dlength^3)*K(ind,ind);

% M = kron(diag(v1),M1)+kron(diag(v2),M3)+kron(diag(u,1),M2)+kron(diag(u,-1),M2');
M = kron(Dv1,M1) + kron(Dv2,M3) + kron(Du,M2) + kron(Du',M2');
M = (rho*area*dlength/420)*M(ind,ind);

% D = zeros(n); D(nele,nele) = damp;
D=sparse(nele,nele,damp,n,n);

end
