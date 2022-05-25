function [output1,output2,output3,output4] = cossenoBF()
erros = [];
erros_log = [];
f = @(t,y,v) v;
g = @(t,y,v) -y;
a = cos(10);
b = 0.1:-0.00005:0.00005;
for passo = b
    h = 0:passo:10;
    y(1) = 1;
    v(1) = 0;

    for i = 1:3
       f1(i) = f(h(i),y(i),v(i));
       k_1 = passo * f1(i);

       k_2 = passo * f(h(i)+passo/2, y(i) + k_1/2,v(i)+k_1/2);

       k_3 = passo * f(h(i)+passo/2, y(i) + k_2/2,v(i)+k_2/2);

       k_4 = passo * f(h(i)+passo, y(i) + k_3,v(i)+k_3);


       y_aux = y(i) + k_1/6 + k_2/3 + k_3/3 + k_4/6;
       y(i+1) = y_aux;

       g1(i) = g(h(i),y(i),v(i));

       k_1 = passo * g1(i);

       k_2 = passo * g(h(i)+passo/2, y(i) + k_1/2,v(i)+k_1/2);

       k_3 = passo * g(h(i)+passo/2, y(i) + k_2/2,v(i)+k_2/2);

       k_4 = passo * g(h(i)+passo, y(i) + k_3,v(i)+k_3);

       v_aux = v(i) + k_1/6 + k_2/3 + k_3/3 + k_4/6;

       v(i+1) = v_aux;


    end

    %algortimo de adams-bashforth de ordem 4

    for i = 4:length(h)-1
       f1_aux = f(h(i), y(i),v(i));
       y_aux = y(i) + (passo/24) * ( 55 * f1_aux - 59 * f1(3) + 37 * f1(2) -9 * f1(1)  );
       y(i+1) = y_aux;
       f1 = [f1(2:end) f1_aux];

       g1_aux = g(h(i),y(i),v(i));
       v_aux = v(i) + (passo/24) * ( 55 * g1_aux - 59 * g1(3) + 37 * g1(2) -9 * g1(1)  );

       v(i+1) = v_aux;
       g1 = [g1(2:end)  g1_aux];

    end
    erro_y = abs(y(end) - a);
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
