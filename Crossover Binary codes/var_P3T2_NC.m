function y = var_P3T2_NC(ps,mu,alpha,tau,R)
p = 3;
%D2 is ABB ABA AAB BAA BAB BBA AAA BBB
%ps = [25 25 25 25];
tau1 = [1 -1 -1]';
tau2 = [1 -1 1]';
tau3 = [1 1 -1]';
tau4 = [-1 1 1]';
tau5 = [-1 1 -1]';
tau6 = [-1 -1 1]';
tau7 = [1,1,1]';
tau8 = [-1,-1,-1]';




X1 = [ones(p,1) [0,0;1,0;0,1] tau1];
X2 = [ones(p,1) [0,0;1,0;0,1] tau2];
X3 = [ones(p,1) [0,0;1,0;0,1] tau3];
X4 = [ones(p,1) [0,0;1,0;0,1] tau4];
X5 = [ones(p,1) [0,0;1,0;0,1] tau5];
X6 = [ones(p,1) [0,0;1,0;0,1] tau6];
X7 = [ones(p,1) [0,0;1,0;0,1] tau7];
X8 = [ones(p,1) [0,0;1,0;0,1] tau8];

beta = [mu alpha tau]';
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
etad7 = eta71.*eta72;

Z8 = X8*beta;
eta81 = exp(Z8)./(1+exp(Z8));
eta82 = 1-eta81;
etad8 = eta81.*eta82;





W1 = diag(sqrt(etad1));
W2 = diag(sqrt(etad2));
W3 = diag(sqrt(etad3));
W4 = diag(sqrt(etad4));
W5 = diag(sqrt(etad5));
W6 = diag(sqrt(etad6));
W7 = diag(sqrt(etad7));
W8 = diag(sqrt(etad8));
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


%V = sigmb*eye(p);

A = ps(1)*X1'*A1*X1 + ps(2)*X2'*A2*X2 + ps(3)*X3'*A3*X3 + ps(4)*X4'*A4*X4 + ps(5)*X5'*A5*X5 + ps(6)*X6'*A6*X6 + ps(7)*X7'*A7*X7 + ps(8)*X8'*A8*X8; 


B = inv(A);
%y=rank(A)
y = log(B(4,4));
