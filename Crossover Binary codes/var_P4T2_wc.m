function y = var_P4T2_wc(ps,mu,alpha,tau,gam,R)
p = 4;

% AABB BBAA ABBA BAAB ABAB BABA / ABBB BABB BBAB BBBA BAAA ABAA AABA AAAB
% AAAA BBBB
%ps = [25 25 25 25];
tau1 = [1 1 -1 -1]';
tau2 = [-1 -1 1 1]';
tau3 = [1 -1 -1 1]';
tau4 = [-1 1 1 -1]';
tau5 = [1 -1 1 -1]';
tau6 = [-1 1 -1 1]';

tau7 = [1 -1 -1 -1]';
tau8 = [-1 1 -1 -1]';
tau9=[-1 -1 1 -1]';
tau10 = [-1 -1 -1 1]';
tau11 = [-1 1 1 1]';
tau12 = [1 -1 1 1]';
tau13 = [1 1 -1 1]';
tau14 = [1 1 1 -1]';
tau15 = [1 1 1 1]';
tau16 = [-1 -1 -1 -1]';


GM = [0 zeros(1,p-1); eye(p-1) zeros(p-1,1)];
gam1 = GM*tau1;
gam2 = GM*tau2;
gam3 = GM*tau3;
gam4 = GM*tau4;
gam5 = GM*tau5;
gam6 = GM*tau6;
gam7 = GM*tau7;
gam8 = GM*tau8;
gam9 = GM*tau9;
gam10 = GM*tau10;
gam11 = GM*tau11;
gam12 = GM*tau12;
gam13 = GM*tau13;
gam14 = GM*tau14;
gam15 = GM*tau15;
gam16 = GM*tau16;





X1 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau1 gam1];
X2 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau2 gam2];
X3 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau3 gam3];
X4 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau4 gam4];
X5 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau5 gam5];
X6 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau6 gam6 ];

X7 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau7 gam7];
X8 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau8 gam8];
X9 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau9 gam9];
X10 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau10 gam10];
X11 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau11 gam11];
X12 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau12 gam12];
X13 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau13 gam13];
X14 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau14 gam14];
X15 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau15 gam15];
X16 = [ones(p,1) [0,0,0;1,0,0;0,1,0;0,0,1] tau16 gam16];




beta = [mu alpha tau gam]';
Z1 = X1*beta;
eta11 = exp(Z1)./(1+exp(Z1));
eta12 = 1-eta11;
etad1 = eta12.*eta11;

Z2 = X2*beta;
eta21 = exp(Z2)./(1+exp(Z2));
eta22 = 1-eta21;
etad2 = eta21.*eta22;

Z3 = X3*beta;
eta31 = exp(Z3)./(1+exp(Z3));
eta32 = 1-eta31;
etad3 = eta31.*eta32;

Z4 = X4*beta;
eta41 = exp(Z4)./(1+exp(Z4));
eta42 = 1-eta41;
etad4 = eta41.*eta42;

Z5 = X5*beta;
eta51 = exp(Z5)./(1+exp(Z5));
eta52 = 1-eta51;
etad5 = eta51.*eta52;

Z6 = X6*beta;
eta61 = exp(Z6)./(1+exp(Z6));
eta62 = 1-eta61;
etad6 = eta61.*eta62;


Z7 = X7*beta;
eta71 = exp(Z7)./(1+exp(Z7));
eta72 = 1-eta71;
etad7 = eta72.*eta71;

Z8 = X8*beta;
eta81 = exp(Z8)./(1+exp(Z8));
eta82 = 1-eta81;
etad8 = eta81.*eta82;

Z9 = X9*beta;
eta91 = exp(Z9)./(1+exp(Z9));
eta92 = 1-eta91;
etad9 = eta91.*eta92;

Z10 = X10*beta;
eta101 = exp(Z10)./(1+exp(Z10));
eta102 = 1-eta101;
etad10 = eta101.*eta102;

Z11 = X11*beta;
eta111 = exp(Z11)./(1+exp(Z11));
eta112 = 1-eta111;
etad11 = eta111.*eta112;

Z12 = X12*beta;
eta121 = exp(Z12)./(1+exp(Z12));
eta122 = 1-eta121;
etad12 = eta121.*eta122;



Z13 = X13*beta;
eta131 = exp(Z13)./(1+exp(Z13));
eta132 = 1-eta131;
etad13 = eta131.*eta132;

Z14 = X14*beta;
eta141 = exp(Z14)./(1+exp(Z14));
eta142 = 1-eta141;
etad14 = eta141.*eta142;

Z15 = X15*beta;
eta151 = exp(Z15)./(1+exp(Z15));
eta152 = 1-eta151;
etad15 = eta151.*eta152;

Z16 = X16*beta;
eta161 = exp(Z16)./(1+exp(Z16));
eta162 = 1-eta161;
etad16 = eta161.*eta162;




W1 = diag(sqrt(etad1));
W2 = diag(sqrt(etad2));
W3 = diag(sqrt(etad3));
W4 = diag(sqrt(etad4));
W5 = diag(sqrt(etad5));
W6 = diag(sqrt(etad6));

W7 = diag(sqrt(etad7));
W8 = diag(sqrt(etad8));
W9 = diag(sqrt(etad9));
W10 = diag(sqrt(etad10));
W11 = diag(sqrt(etad11));
W12 = diag(sqrt(etad12));
W13 = diag(sqrt(etad13));
W14 = diag(sqrt(etad14));
W15 = diag(sqrt(etad15));
W16 = diag(sqrt(etad16));


% V1 = sigb*W1*W1;
% V2 = sigb*W2*W2;
% V3 = sigb*W3*W3;


%R = rho*ones(p)+(1-rho)*eye(p);
R1 = eye(p)/(R);
A1 = W1*R1*W1;
A2 = W2*R1*W2;
A3 = W3*R1*W3;
A4 = W4*R1*W4;
A5 = W5*R1*W5;
A6 = W6*R1*W6;

A7 = W7*R1*W7;
A8 = W8*R1*W8;
A9 = W9*R1*W9;
A10 = W10*R1*W10;
A11 = W11*R1*W11;
A12 = W12*R1*W12;
A13 = W13*R1*W13;
A14 = W14*R1*W14;
A15 = W15*R1*W15;
A16 = W16*R1*W16;


%V = sigmb*eye(p);

A = ps(1)*X1'*A1*X1 + ps(2)*X2'*A2*X2 + ps(3)*X3'*A3*X3 + ps(4)*X4'*A4*X4 + ps(5)*X5'*A5*X5 + ps(6)*X6'*A6*X6 + ps(7)*X7'*A7*X7 + ps(8)*X8'*A8*X8 + ps(9)*X9'*A9*X9 + ps(10)*X10'*A10*X10 + ps(11)*X11'*A11*X11 + ps(12)*X12'*A12*X12 + ps(13)*X13'*A13*X13 + ps(14)*X14'*A14*X14 + ps(15)*X15'*A15*X15 + ps(16)*X16'*A16*X16;

B = inv(A);
%y=eig(A)
%y=rank(A)
y = B(6,6);
