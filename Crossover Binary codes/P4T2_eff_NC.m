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


% point estimates for CS correlation
%     
% mu = -0.295;
% alpha(1) = 0.385;
% alpha(2) = -0.383;
% alpha(3) = -0.205;
% tau = 0.11;

% point estimates for AR correlation
    
mu = -0.299;
alpha(1) = 0.385;
alpha(2) = -0.378;
alpha(3) = -0.208;
tau = 0.128;





fun_wc = @(ps)var_P4T2_NC(ps,mu,alpha,tau,R)
 


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
%     
% mu = -1.16+rand(1)*(0.57+1.16);
% alpha(1) = -0.683+rand(1)*(1.453+0.683);
% alpha(2) = -1.381+rand(1)*(0.615+1.381);
% alpha(3) = -1.05+rand(1)*(0.64+1.05);
% tau = -0.227+rand(1)*(0.447+0.227);


  % intervals for AR correlation
  
mu = -1.161+rand(1)*(0.563+1.161);
alpha(1) = -0.685+rand(1)*(1.455+0.685);
alpha(2) = -1.374+rand(1)*(0.618+1.374);
alpha(3) = -1.057+rand(1)*(0.641+1.057);
tau = -0.242+rand(1)*(0.498+0.242);







var_opt_NC= var_P4T2_NC(ps_opt,mu,alpha,tau,R);

ps_1 = [1/4,1/4,1/4,1/4,0,0,0,0,0,0,0,0,0,0,0,0];
ps_2 = [1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
ps_3 = [0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0];
ps_4 = [0,0,0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0];
var_ex_1(i)= var_P4T2_NC(ps_1,mu,alpha,tau,R);
var_ex_2(i)= var_P4T2_NC(ps_2,mu,alpha,tau,R);
var_ex_3(i)= var_P4T2_NC(ps_3,mu,alpha,tau,R);
var_ex_4(i)= var_P4T2_NC(ps_4,mu,alpha,tau,R);


eff_1(i) = (exp(var_opt_NC)./exp(var_ex_1(i)))^(1/5);
eff_2(i) = (exp(var_opt_NC)./exp(var_ex_2(i)))^(1/5);
eff_3(i) = (exp(var_opt_NC)./exp(var_ex_3(i)))^(1/5);
eff_4(i) = (exp(var_opt_NC)./exp(var_ex_4(i)))^(1/5);



end
%med = median(var_opt_nc);
figure(1)
boxplot([eff_1',eff_2',eff_3',eff_4'],'labels',{'$$\Gamma_1$$','$$\Gamma_2$$','$$\Gamma_3$$','$$\Gamma_4$$'})
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
ylabel('Efficiency') 
%hold on
% plot(med)