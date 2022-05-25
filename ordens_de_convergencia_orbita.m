tic
y_exato = [-1.316567703298366e+01     1.172592500454820e+01];
erros_RK = [];
erros_BF = [];
erros_BS = [];
b = 0.1:-0.001:0.001;
for h = b
    [y_final_RK] = planetasRK([0.8,0.8],[-2.55,0.5],h,30,false);
    [y_final_BF] = planetasBF([0.8,0.8],[-2.55,0.5],h,30,false);
    [y_final_BS] = planetasBS([0.8,0.8],[-2.55,0.5],h,30,1e-8,h,0.1,false);
    erro1 = norm(y_exato - y_final_RK);
    erro2 = norm(y_exato - y_final_BF);
    erro3 = norm(y_exato - y_final_BS);
    erros_RK = [erros_RK erro1];
    erros_BF = [erros_BF erro2];
    erros_BS = [erros_BS erro3];
end
c = log(b);

log_RK = log(erros_RK);
log_BF = log(erros_BF);
log_BS = log(erros_BS);

Y1 = log_RK';
Y2 = log_BF';
Y3 = log_BS';
X(:,1) = ones(1,length(Y1));
X(:,2) = c;

r1 = linsolve(X.' * X,X.' * Y1); 
r2 = linsolve(X.' * X,X.' * Y2); 
r3 = linsolve(X.' * X,X.' * Y3); 


for i = 1:length(c)
   reta1(i) = r1(1) + r1(2) * c(i);
   reta2(i) = r2(1) + r2(2) * c(i);
   reta3(i) = r3(1) + r3(2) * c(i);
end
 ordem_de_convergencia_de_runge_kutta = r1(2)
 ordem_de_convergencia_de_Adams_Bashforth = r2(2)
 ordem_de_convergencia_de_Bogacki_shampine = r3(2)
 
subplot(1,3,1)
plot(c, log_RK,'.', c, reta1, 'red')
legend('logaritmo do erro em função do logaritmo de h', 'aproximação dos mínimos quadrados');


subplot(1,3,2)
plot(c, log_BF,'.', c, reta2, 'red')
legend('logaritmo do erro em função do logaritmo de h', 'aproximação dos mínimos quadrados');

subplot(1,3,3)
plot(c, log_BS,'.', c, reta3, 'red')
legend('logaritmo do erro em função do logaritmo de h', 'aproximação dos mínimos quadrados');
toc

