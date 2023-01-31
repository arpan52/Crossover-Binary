clear all;
% correlation
p=3;
% compound Symmetric
al = 0.4;
% 
R=al*ones(p,p)
for i=1:p
    R(i,i)=1;
end

% AR(1) correlation

% R=al*ones(p,p);
% R(1,3)=al*al;
% R(3,1)=al*al;
% for i=1:p
%     R(i,i)=1;
% end

% point estimates for CS correlation
mu = -0.2116;
alpha(1) = 0.2121;
alpha(2) = -0.4377;
tau = 0.4833;
gam = 0.527;

% point estimates for AR correlation

% mu = -0.212;
% alpha(1) = 0.213;
% alpha(2) = -0.437;
% tau = 0.484;
% gam = 0.532;

fun_wc = @(ps)var_P3T2_wc(ps,mu,alpha,tau,gam,R)
 


Aeq = [1,1,1,1,1,1,1,1];
beq = 1;
lb = [0,0,0,0,0,0,0,0];
ub = [1,1,1,1,1,1,1,1];
A = [];
b = [];


ps0=[1/8,1/8,1/8,1/8,1/8,1/8,1/8,1/8];

% ABB ABA AAB BAA BAB BBA AAA BBB

ps_opt = fmincon(fun_wc,ps0,A,b,Aeq,beq,lb,ub)


for i=1:5000
    % intervals for CS correlation
mu = -0.8210+rand(1)*(0.3978+0.8210);
alpha(1) = -0.5991+rand(1)*(1.0233+0.5991);
alpha(2) = -1.3311+rand(1)*(0.4557+1.3311);
tau = 0.1178+rand(1)*(0.8488-0.1178);
gam = 0.0976+rand(1)*(0.9564-0.0976);


  % intervals for AR correlation
  
% mu = -0.8216+rand(1)*(0.3976+0.8216);
% alpha(1) = -0.5984+rand(1)*(1.0244+0.5984);
% alpha(2) = -1.3308+rand(1)*(0.4568+1.3308);
% tau = 0.1175+rand(1)*(0.8505-0.1175);
% gam = 0.0988+rand(1)*(0.9652-0.0988);





var_opt_wc= var_P3T2_wc(ps_opt,mu,alpha,tau,gam,R);

ps_1 = [1/2,0,0,1/2,0,0,0,0];
ps_2 = [1/4,0,1/4,1/4,0,1/4,0,0];
ps_3 = [1/4,1/4,0,1/4,1/4,0,0,0];
var_ex_1(i)= var_P3T2_wc(ps_1,mu,alpha,tau,gam,R);
var_ex_2(i)= var_P3T2_wc(ps_2,mu,alpha,tau,gam,R);
var_ex_3(i)= var_P3T2_wc(ps_3,mu,alpha,tau,gam,R);


eff_1(i) = (exp(var_opt_wc)./exp(var_ex_1(i)))^(1/5);
eff_2(i) = (exp(var_opt_wc)./exp(var_ex_2(i)))^(1/5);
eff_3(i) = (exp(var_opt_wc)./exp(var_ex_3(i)))^(1/5);


end
%med = median(var_opt_nc);
figure(1)
boxplot([eff_1',eff_2',eff_3'],'labels',{'$$\Gamma_1$$','$$\Gamma_2$$','$$\Gamma_3$$'})
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';

%hold on
% plot(med)