clear all;
% correlation
p=3;
% compound Symmetric
al = 0;
 
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
mu = -.207;
alpha(1) = 0.211;
alpha(2) = -0.326;
tau = 0.449;


% point estimates for AR correlation

% mu = -0.2106;
% alpha(1) = 0.2143;
% alpha(2) = -0.3199;
% tau = 0.4268;


fun_wc = @(ps)var_P3T2_NC(ps,mu,alpha,tau,R)
 



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
    
mu = -0.8185+rand(1)*(0.4045+0.8185);
alpha(1) = -0.6396+rand(1)*(1.0616+0.6396);
alpha(2) = -1.1237+rand(1)*(0.4717+1.1237);
tau = 0.1021+rand(1)*(0.7959-0.1021);


  % intervals for AR correlation
  
% mu = -0.8231+rand(1)*(0.4019+0.8231);
% alpha(1) = -0.6334+rand(1)*(1.0620+0.6334);
% alpha(2) = -1.1162+rand(1)*(0.4764+1.1162);
% tau = 0.0722+rand(1)*(0.7814-0.0722);




var_opt_NC= var_P3T2_NC(ps_opt,mu,alpha,tau,R);

ps_1 = [1/2,0,0,1/2,0,0,0,0];
ps_2 = [1/4,0,1/4,1/4,0,1/4,0,0];
ps_3 = [1/4,1/4,0,1/4,1/4,0,0,0];
var_ex_1(i)= var_P3T2_NC(ps_1,mu,alpha,tau,R);
var_ex_2(i)= var_P3T2_NC(ps_2,mu,alpha,tau,R);
var_ex_3(i)= var_P3T2_NC(ps_3,mu,alpha,tau,R);


eff_1(i) = (exp(var_opt_NC)./exp(var_ex_1(i)))^(1/4);
eff_2(i) = (exp(var_opt_NC)./exp(var_ex_2(i)))^(1/4);
eff_3(i) = (exp(var_opt_NC)./exp(var_ex_3(i)))^(1/4);



end
%med = median(var_opt_nc);
figure(1)
boxplot([eff_1',eff_2',eff_3'],'labels',{'$$\Gamma_1$$','$$\Gamma_2$$','$$\Gamma_3$$'})
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
ylabel('Efficiency') 
%hold on
% plot(med)