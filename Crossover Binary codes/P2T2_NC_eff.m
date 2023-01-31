clear all;
% correlation
p=2;

al = 0.4;
% compound Symmetric
R=al*ones(p,p)
for i=1:p
    R(i,i)=1;
end


% point estimates of parameters
mu =0.7125;
alpha = -0.1794;
tau = 0.1333;

fun_nc = @(ps)var_P2T2_NC(ps,mu,alpha,tau,R)
 %ps=[1/4,1/4,1/4,1/4];
%fun = P2T2(ps,mu,alpha,tau,gam,R)
Aeq = [1,1,1,1];
beq = 1;
lb = [0,0,0,0];
ub = [1,1,1,1];
A = [];
b = [];


ps0=[1/4,1/4,1/4,1/4];

% AB AA BA BB

ps_opt = fmincon(fun_nc,ps0,A,b,Aeq,beq,lb,ub)


for i=1:5000
mu = 0.2997+rand(1)*(1.1253-0.2997);
alpha = -0.5600+rand(1)*(0.2012+0.5600);
tau = -0.0572+rand(1)*(0.3238+0.0572);

%mu=0.5;
%alpha=[1];
%tau=-2;



var_opt_nc= var_P2T2_NC(ps_opt,mu,alpha,tau,R);

ps_1 = [1/2,0,1/2,0];
ps_2 = [1/4,1/4,1/4,1/4];
var_ex_1(i)= var_P2T2_NC(ps_1,mu,alpha,tau,R);
var_ex_2(i)= var_P2T2_NC(ps_2,mu,alpha,tau,R);

eff_1(i) = (exp(var_opt_nc)./exp(var_ex_1(i)))^(1/3);
eff_2(i) = (exp(var_opt_nc)./exp(var_ex_2(i)))^(1/3);



end

ps_opt
figure(1)
boxplot([eff_1',eff_2'],'labels',{'$$\Gamma_1$$','$$\Gamma_2$$'});

bp = gca;

bp.XAxis.TickLabelInterpreter = 'latex';
ylabel('Efficiency') 
%hold on 
%hold on
% plot(med)