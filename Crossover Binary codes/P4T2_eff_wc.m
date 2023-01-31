clear all;
% correlation
p=4;
% compound Symmetric
al = 0.4;

% R=al*ones(p,p);
% for i=1:p
%     R(i,i)=1;
% end

% AR(1) correlation

for i=1:p
    for j=1:p
        R(i,j) = al^(abs(i-j));
    end
end

% point estimate for CS
% mu = -0.253;
% alpha(1) = 0.404;
% alpha(2) = -0.526;
% alpha(3) = -0.265;
% tau = -0.143;
% gam = -0.484;

% point estimates for AR1

mu = -0.263;
alpha(1) = 0.411;
alpha(2) = -0.518;
alpha(3) = -0.251;
tau = -0.16;
gam = -0.469;



fun_wc = @(ps)var_P4T2_wc(ps,mu,alpha,tau,gam,R)
 


Aeq = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
beq = 1;
lb = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
ub = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0];
A = [];
b = [];


ps0=[1/16,1/16,1/16,1/16,1/16,1/16,1/16,1/16,1/16,1/16,1/16,1/16,1/16,1/16,1/16,1/16];

% AABB BBAA ABBA BAAB ABAB BABA

ps_opt = fmincon(fun_wc,ps0,A,b,Aeq,beq,lb,ub)



for i=1:5000
    % intervals for CS correlation
    
% mu = -1.129+rand(1)*(0.623+1.129);
% alpha(1) = -0.684+rand(1)*(1.492+0.684);
% alpha(2) = -1.518+rand(1)*(0.466+1.518);
% alpha(3) = -1.133+rand(1)*(0.603+1.133);
% tau = -0.566+rand(1)*(0.280+0.566);
% gam = -1.0132+rand(1)*(0.0452+1.0132);


  % intervals for AR correlation
  
mu = -1.141+rand(1)*(0.615+1.141);
alpha(1) = -0.671+rand(1)*(1.493+0.671);
alpha(2) = -1.498+rand(1)*(0.462+1.498);
alpha(3) = -1.117+rand(1)*(0.615+1.117);
tau = -0.638+rand(1)*(0.318+0.638);
gam = -1.047+rand(1)*(0.109+1.047);









var_opt_wc= var_P4T2_wc(ps_opt,mu,alpha,tau,gam,R);

ps_1 = [1/4,1/4,1/4,1/4,0,0,0,0,0,0,0,0,0,0,0,0];
ps_2 = [1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
ps_3 = [0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0];
ps_4 = [0,0,0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0];
var_ex_1(i)= var_P4T2_wc(ps_1,mu,alpha,tau,gam,R);
var_ex_2(i)= var_P4T2_wc(ps_2,mu,alpha,tau,gam,R);
var_ex_3(i)= var_P4T2_wc(ps_3,mu,alpha,tau,gam,R);
var_ex_4(i)= var_P4T2_wc(ps_4,mu,alpha,tau,gam,R);


eff_1(i) = (exp(var_opt_wc)./exp(var_ex_1(i)))^(1/6);
eff_2(i) = (exp(var_opt_wc)./exp(var_ex_2(i)))^(1/6);
eff_3(i) =(exp(var_opt_wc)./exp(var_ex_3(i)))^(1/6);
eff_4(i) = (exp(var_opt_wc)./exp(var_ex_4(i)))^(1/6);


end
%med = median(var_opt_nc);
figure(1)
boxplot([eff_1',eff_2',eff_3',eff_4'],'labels',{'$$\Gamma_1$$','$$\Gamma_2$$','$$\Gamma_3$$','$$\Gamma_4$$'})
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
ylabel('Efficiency') 
%hold on
%hold on
% plot(med)