clear all;
% correlation
p=2;
% compound Symmetric
al = 0.1;
% compound Symmetric
R=al*ones(p,p)
for i=1:p
    R(i,i)=1;
end

% point estimates of parameters
mu =0.717;
alpha = -0.185;
tau = 0.227;
gam = 0.185

fun_wc = @(ps)var_P2T2_wc(ps,mu,alpha,tau,gam,R)
 


Aeq = [1,1,1,1];
beq = 1;
lb = [0,0,0,0];
ub = [1,1,1,1];
A = [];
b = [];


ps0=[1/4,1/4,1/4,1/4];

% AB AA BA BB

ps_opt = fmincon(fun_wc,ps0,A,b,Aeq,beq,lb,ub)


for i=1:5000
mu = 0.2976+rand(1)*(1.1364-0.2976);
alpha = -0.5652+rand(1)*(0.1352+0.5652);
tau = -0.1924+rand(1)*(0.6464+0.1924);
gam = -0.5441+rand(1)*(0.9141+0.5441);
%mu=0.5;
%alpha=[1];
%tau=-2;


var_opt_wc= var_P2T2_wc(ps_opt,mu,alpha,tau,gam,R);

ps_1 = [1/2,0,1/2,0];
ps_2 = [1/4,1/4,1/4,1/4];
var_ex_1(i)= var_P2T2_wc(ps_1,mu,alpha,tau,gam,R);
var_ex_2(i)= var_P2T2_wc(ps_2,mu,alpha,tau,gam,R);

eff_1(i) = (exp(var_opt_wc)./exp(var_ex_1(i)))^(1/4);
eff_2(i) = (exp(var_opt_wc)./exp(var_ex_2(i)))^(1/4);


end
%med = median(var_opt_nc);
figure(1)
boxplot([eff_1',eff_2'],'labels',{'$$\Gamma_1$$','$$\Gamma_2$$'});

bp = gca;

bp.XAxis.TickLabelInterpreter = 'latex';
ylabel('Efficiency') 
%hold on
% plot(med)