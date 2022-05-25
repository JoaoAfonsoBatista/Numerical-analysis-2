function [output1,output2,output3,output4] = cossenoRK()
erros = [];
erros_log = [];
f = @(t,y,v) v;
g = @(t,y,v) -y;
a = cos(10);
b = 0.1:-0.00005:0.00005;
for h = b
    malha = 0:h:10;
    passo = h;
    y(1) = 1;
    v(1) = 0;

    for i = 1:length(malha)-1
       k_1 = passo * f(malha(i),y(i),v(i));

       k_2 = passo * f(malha(i)+passo/2, y(i) + k_1/2,v(i)+k_1/2);

       k_3 = passo * f(malha(i)+passo/2, y(i) + k_2/2,v(i)+k_2/2);

       k_4 = passo * f(malha(i)+passo, y(i) + k_3,v(i)+k_3);

       y(i+1) = y(i) + k_1/6 + k_2/3 + k_3/3 + k_4/6;



       k_1 = passo * g(malha(i),y(i),v(i));

       k_2 = passo * g(malha(i)+passo/2, y(i) + k_1/2,v(i)+k_1/2);

       k_3 = passo * g(malha(i)+passo/2, y(i) + k_2/2,v(i)+k_2/2);

       k_4 = passo * g(malha(i)+passo, y(i) + k_3,v(i)+k_3);

       v(i+1) = v(i) + k_1/6 + k_2/3 + k_3/3 + k_4/6;
    end
    erro_y = abs(y(end) - a);
    erros = [erros erro_y];
    erros_log = [erros_log log(erro_y)];
end
 Y = erros_log';
 X(:,1) = ones(1,length(Y));
 X(:,2) = log(b);
 
 r = linsolve(X.' * X,X.' * Y);
 

 
 for i = 1:length(log(b))
    reta(i) = r(1) + r(2) * log(b(i));
 end
 output1 = log(b);
 output2 = erros_log;
 output3 = reta;
 output4 = r(2);
end
