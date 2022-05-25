function [output1,output2,output3,output4] = cossenoBS()
erros = [];
erros_log = [];
f = @(t,y,v) v;
g = @(t,y,v) -y;
a = cos(10);
b = 0.1:-0.0001:0.0001;
limite = 10;
epsilon = 1e-12;
passo_max = 0.1;

for passo_min = b
    h(1) = passo_min;
    passo = passo_min;
    y(1) = 1;
    v(1) = 0;
    i = 1;
    while h(i) < limite
       ky_1 = f(h(i),y(i),v(i));

       ky_2 = f(h(i)+passo/2, y(i) + passo * ky_1/2,v(i)+ passo * ky_1/2);

       ky_3 = f(h(i)+3*passo/4, y(i) + 3*passo * ky_2/4,v(i)+ 3*passo * ky_2/4);

       y_aux = y(i) + (2/9) * passo * ky_1 + (1/3) * passo * ky_2 + (4/9) * passo * ky_3;


       kv_1 = g(h(i),y(i),v(i));

       kv_2 = g(h(i)+passo/2, y(i) + passo * kv_1/2,v(i)+ passo * kv_1/2);

       kv_3 = g(h(i)+3*passo/4, y(i) + 3*passo * kv_2/4,v(i)+ 3*passo * kv_2/4);

       v_aux = v(i) + (2/9) * passo * kv_1 + (1/3) * passo * kv_2 + (4/9) * passo * kv_3;



       ky_4 = f(h(i) + passo, y_aux, v_aux);

       zy = y(i) + passo*(7/24)*ky_1 + (passo/4)*ky_2 + (passo/3)*ky_3 + (passo/8)*ky_4;

       kv_4 = g(h(i) + passo, y_aux, v_aux);

       zv = v(i) + passo*(7/24)*kv_1 + (passo/4)*kv_2 + (passo/3)*kv_3 + (passo/8)*kv_4;


       erro_y = norm(y_aux - zy);

       erro_v = norm(v_aux - zv);

       passo1 = passo * sqrt(epsilon/erro_y);

       passo2 = passo * sqrt(epsilon/erro_v);

       %entao recalculamos o y desta iteraçao com o nosso
       %passo________________________________________________________________________________________________________________________________________

       passo_aux1 = min([passo1,passo2,passo_max]);
       passo = max([passo_min,passo_aux1]);

       ky_2 = f(h(i)+passo/2, y(i) + passo * ky_1/2,v(i)+ passo * ky_1/2);

       ky_3 = f(h(i)+3*passo/4, y(:) + 3*passo * ky_2/4,v(i)+ 3*passo * ky_2/4);

       y_aux = y(i) + (2/9) * passo * ky_1 + (1/3) * passo * ky_2 + (4/9) * passo * ky_3;


       kv_2 = g(h(i)+passo/2, y(i) + passo * kv_1/2,v(i)+ passo * kv_1/2);

       kv_3 = g(h(i)+3*passo/4, y(i) + 3*passo * kv_2/4,v(i)+ 3*passo * kv_2/4);

       v_aux = v(i) + (2/9) * passo * kv_1 + (1/3) * passo * kv_2 + (4/9) * passo * kv_3;


       y(i+1) = y_aux;
       v(i+1) = v_aux;
       h(i+1) = h(i) + passo;
       i = i+1;

    end
    erro_y = abs(y(end) - cos(h(end)));
    erros = [erros erro_y];
    erros_log = [erros_log log(erro_y)];
end
 Y = erros_log';
 X(:,1) = ones(1,length(Y));
 X(:,2) = log(b);
 
 r = linsolve(X.' * X, X.' * Y);
 

 
 for i = 1:length(log(b))
    reta(i) = r(1) + r(2) * log(b(i));
 end
 output1 = log(b);
 output2 = erros_log;
 output3 = reta;
 output4 = r(2);
end
