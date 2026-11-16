% Perform the tests for Table 2 in the paper
% "Fast Algorithms for Optimal Damping in Mechanical Systems
%  by Q Li, F Tisseur, arXiv.2601.05404
%  last updated by FT on 04/08/2026


opts.tolres = 1e-8;
opts.tolnu = 1e-5;
opts.f = true;

fprintf('Problem    n  k   s   nu_0  Alg. #iter   #ls  #eig   f(nu_*) \n');
fprintf('------------------------------------------------------------\n');


for  pb = (1:9)
  switch pb
    case 1
      %fprintf('damp1 with n=4, s=4, k=1\n')
      pb_name = 'damp1-a';
      n = 4;
      opts.s = 4;
      M = diag(1:n);
      K = diag(10*ones(n,1))+diag(-5*ones(n-1,1),-1)+diag(-5*ones(n-1,1),1);
      k = 1;
      R = zeros(n,1);  %one damper
      R(2,1) = 1;
      damping = 1;
      a = 0.01; %critical damping
      opts.nu0 = 1.0*ones(k,1); %initial viscosities
      nu_0 = 'e';

    case 2
      %fprintf('damp1 with n=20, s=20, k=1\n')
      pb_name = 'damp1-b';
      n = 20;
      opts.s = 20;
      M = diag(1:n);
      K = diag(50*ones(n,1))-diag(25*ones(n-1,1),-1)-diag(25*ones(n-1,1),1);
      k = 1;
      R = zeros(n,1);  %one damper
      R(2,1) = 1;
      damping = 1;
      a = 0.01; %critical damping
      opts.nu0 = 1.0*ones(k,1); %initial viscosities
      nu_0 = 'e';

    case 3
      %fprintf('damp1 with n=20, s=20, k=2\n')
      pb_name = 'damp1-c';
      n = 20;
      opts.s = 20;
      M = diag(1:n);
      K = diag(50*ones(n,1))-diag(25*ones(n-1,1),-1)-diag(25*ones(n-1,1),1);
      k = 2;
      R = zeros(n,2);  %two dampers
      R(2,1) = 1;
      R(19,2) = 1;
      damping = 1;
      a = 0.01;  %critical damping
      opts.nu0 = 1.0*ones(k,1); %initial viscosities
      nu_0 = 'e';

   case 4
      %fprintf('damp1 with n=20, s=20, k=2\n')
      pb_name = 'damp1-c';
      n = 20;
      opts.s = 20;
      M = diag(1:n);
      K = diag(50*ones(n,1))-diag(25*ones(n-1,1),-1)-diag(25*ones(n-1,1),1);
      k = 2;
      R = zeros(n,2);  %two dampers
      R(2,1) = 1;
      R(19,2) = 1;
      damping = 1;
      a = 0.01;   %critical damping
      opts.nu0 = 10.0*ones(k,1); %initial viscosities
      nu_0 = '10e';

   case 5
      %fprintf('damp2 with n=801, s=27, k=3\n')
      pb_name = 'damp2-a';
      t  = 400;
      n = 2*t+1;
      m1 = 5*t/2;  m2 = 3*t/2+600;
      k1 = 100; k2 = 150; k3 = 200;
      M = diag([[m1-4*(1:t/2)] [3*(t/2+1:t)-t] [500+(1:t)] m2]);
      K = zeros(n,n);
      Kt = diag(-ones(t-1,1),-1)+ 2*diag(ones(t,1)) - diag(ones(t-1,1),1);
      K(1:t,1:t) = k1*Kt;
      K(t+1:2*t,t+1:2*t) = k2*Kt;
      K(n,t) = -k1; K(n,2*t) = -k2;
      K(t,n) = -k1; K(2*t,n) = -k2;
      K(n,n) = k1+k2+k3;
      k = 3; %three dampers
      R = zeros(n,3);
      l = [50;550;220];
      R(l(1),1) = 1;
      R(l(2),2) = 1; R(l(3)+t,2) = -1;
      R(l(3),3) = 1;
      damping = 1;
      a = 0.02; %critical damping
      opts.s = 27;
      opts.nu0 = ones(k,1); %initial viscosities
      nu_0 = 'e';

   case 6
      %fprintf('damp2 with n=1601, s=27, k=3\n')
      pb_name = 'damp2-b';
      t  = 800;
      n = 2*t+1;
      m1 =5*t/2; m2 = 3*t/2+600;
      k1 = 100; k2 = 150; k3 = 200;
      M = diag([[m1-4*(1:t/2)] [3*(t/2+1:t)-t] [500+(1:t)] m2]);
      K = zeros(n,n);
      Kt = diag(-ones(t-1,1),-1)+ 2*diag(ones(t,1)) - diag(ones(t-1,1),1);
      K(1:t,1:t) = k1*Kt;
      K(t+1:2*t,t+1:2*t) = k2*Kt;
      K(n,t) = -k1; K(n,2*t) = -k2;
      K(t,n) = -k1; K(2*t,n) = -k2;
      K(n,n) = k1+k2+k3;
      k = 3; %three dampers
      R = zeros(n,3);
      l = [50;950;120];
      R(l(1),1) = 1;
      R(l(2),2) = 1; R(l(3)+t,2) = -1;
      R(l(3),3) = 1;
      damping = 1;
      a = 0.02; %critical damping
      opts.s = 27;
      opts.nu0 = ones(k,1); %initial viscosities
      nu_0 = 'e';

   case 7
      %fprintf('damp2 with n=2001, s=20, k=3\n')
      pb_name = 'damp2-c';
      t  = 1000;
      n = 2*t+1;
      m2 = 2000;
      k1 = 100; k2 = 150; k3 = 200;
      M = diag([1000*ones(1,t) 1500*ones(1,t)  m2]);
      K = zeros(n,n);
      Kt = diag(-ones(t-1,1),-1)+ 2*diag(ones(t,1)) - diag(ones(t-1,1),1);
      K(1:t,1:t) = k1*Kt;
      K(t+1:2*t,t+1:2*t) = k2*Kt;
      K(n,t) = -k1; K(n,2*t) = -k2;
      K(t,n) = -k1; K(2*t,n) = -k2;
      K(n,n) = k1+k2+k3;
      k = 3; %three dampers
      R = zeros(n,3);
      l = [850;1950;120];
      R(l(1),1) = 1;
      R(l(2),2) = 1; R(l(3)+t,2) = -1;
      R(l(3),3) = 1;
      damping = 1;
      a = 0.02; %critical damping
      opts.s = 20;
      opts.nu0 = ones(k,1); %initial viscosities
      nu_0 = 'e';

   case 8
      %fprintf('beam-a with n=200, s=40, k=3\n')
      pb_name = 'beam-a';
      n = 200 ;
      opts.s = 40; %number of evals/vecs to dampen
      l = [50 100 150]; % Dampers' positions
      damping = 1;
      a = 0.2; %critical damping
      [M,D,K]= damped_beam(n);
      M = full(M); D = full(D); K = full(K);
      k = 3;
      R = zeros(n,3);
      R(l(1),1) = 1;
      R(l(2),2) = 1;
      R(l(3),3) = 1;
      opts.nu0 = ones(3,1); %initial viscosities
      nu_0 = 'e';


   case 9
      %fprintf('beam-b with n=1000, s=40, k=5\n')
      pb_name = 'beam-b';
      n = 1000;
      opts.s = 80; %number of evals/vecs to dampen
      l = [150 300 500 700 850]; % Dampers' positions
      damping = 1;
      a = 0.2; %critical damping
      [M,D,K]= damped_beam(n);
      M = full(M); D = full(D); K = full(K);
      k = 5;
      R = zeros(n,k);
      for i = 1:length(l)
        R(l(i),i) = 1;
      end
      opts.nu0 = ones(k,1); %initial viscosities
      nu_0 = 'e';

  end

  %% Solve (Mx^2+K)v=0
  [Phi,K1] = eig(K,M,'chol');
  opts.Phi = Phi;
  Om2 = diag(K1); %vector of eigenvalues
  opts.omega = sqrt(Om2);
  if damping == 1
     opts.gamma = a*opts.omega; %Critical damping
  elseif damping == 2
     opts.gamma = a*ones(n,1)+ b*omega; %Rayleigh
  end
  opts.PhiD = Phi.'*R; %This is Dbar for Qingna
  II = eye(2*n);
  G = [II(:,1:opts.s) II(:,n+1:n+opts.s)];
  opts.G = G./sqrt(2*opts.s);
  opts.f = true;
  clear info_BBPR
  clear info_SPG

  [nu_BBPR,info_BBPR] = BBPR(opts); info_BBPR.iter_ls =[];

  [nu_SPG,info_SPG] = SPG(opts);

  fprintf('%-7s  %2.0f %3.0f  %3.0f %4s %-6s  %5.0f %4.0f  %4.0f   %5.1e \n',...
  pb_name,n,k,opts.s,nu_0,'  BBPR',info_BBPR.iter,info_BBPR.iter_ls,info_BBPR.it_eig,info_BBPR.f);
  fprintf('%-7s  %2.0f %3.0f  %3.0f %4s  %-6s %5.0f %4.0f  %4.0f   %5.1e \n',...
  pb_name,n,k,opts.s,nu_0,'  SPG ',info_SPG.iter,info_SPG.iter_ls,info_SPG.it_eig,info_SPG.trace);

end

