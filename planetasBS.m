function [y_final] = planetasBS(d,v_0,passo,limite,epsilon,passo_min,passo_max,plot1)
%malha temporal:
h = [0];
passos = [passo];

planeta1_massa = 2;
planeta1_centro = [10,0];
planeta1_r = 2;
planeta1_p = 0.1;
planeta1_tau = -10;

planeta2_massa = 5;
planeta2_centro = [-10,-10];
planeta2_r = 2;
planeta2_p = 0.1;
planeta2_tau = 5;

planeta3_massa = 10;
planeta3_centro = [0,0];
planeta3_r = 25;
planeta3_p = 0.1;
planeta3_tau = 2.9;

planeta_pos = @(t, centro,r,p,tau)  centro + r * [cos(p*(t - tau)) , sin(p*(t - tau))];

planeta1_pos(1,:) = planeta_pos(0,planeta1_centro,planeta1_r,planeta1_p,planeta1_tau);
planeta2_pos(1,:) = planeta_pos(0,planeta2_centro,planeta2_r,planeta2_p,planeta2_tau);
planeta3_pos(1,:) = planeta_pos(0,planeta3_centro,planeta3_r,planeta3_p,planeta3_tau);

f = @(t,y,v) v;

g = @(t,y,v,planeta1_pos,planeta2_pos,planeta3_pos) planeta3_massa * ((planeta3_pos-y) / (abs(norm(planeta3_pos-y))^3)) + planeta2_massa * ((planeta2_pos-y) / (abs(norm(planeta2_pos-y))^3)) + planeta1_massa * ((planeta1_pos-y) / (abs(norm(planeta1_pos-y))^3));

y(1,:) = planeta1_pos + d;
v(1,:) =  v_0;
i = 1;

distancia1(1) = norm(y(1,:) - planeta1_pos(1,:));
distancia2(1) = norm(y(1,:) - planeta2_pos(1,:));
distancia3(1) = norm(y(1,:) - planeta3_pos(1,:));

distancia4(1) = norm(planeta1_pos(1,:) - planeta2_pos(1,:));
distancia5(1) = norm(planeta1_pos(1,:) - planeta3_pos(1,:));
distancia6(1) = norm(planeta2_pos(1,:) - planeta3_pos(1,:));
while h(i) < limite
   ky_1 = f(h(i),y(i,:),v(i,:));
   
   ky_2 = f(h(i)+passo/2, y(i,:) + passo * ky_1/2,v(i,:)+ passo * ky_1/2);
   
   ky_3 = f(h(i)+3*passo/4, y(i,:) + 3*passo * ky_2/4,v(i,:)+ 3*passo * ky_2/4);
    
   y_aux = y(i,:) + (2/9) * passo * ky_1 + (1/3) * passo * ky_2 + (4/9) * passo * ky_3;
   

   kv_1 = g(h(i),y(i,:),v(i,:),planeta1_pos(i,:),planeta2_pos(i,:),planeta3_pos(i,:));
   
   kv_2 = g(h(i)+passo/2, y(i,:) + passo * kv_1/2,v(i,:)+ passo * kv_1/2,planeta1_pos(i,:),planeta2_pos(i,:),planeta3_pos(i,:));
   
   kv_3 = g(h(i)+3*passo/4, y(i,:) + 3*passo * kv_2/4,v(i,:)+ 3*passo * kv_2/4,planeta1_pos(i,:),planeta2_pos(i,:),planeta3_pos(i,:));
    
   v_aux = v(i,:) + (2/9) * passo * kv_1 + (1/3) * passo * kv_2 + (4/9) * passo * kv_3;
   

   
   ky_4 = f(h(i) + passo, y_aux, v_aux);
   
   zy = y(i,:) + passo*(7/24)*ky_1 + (passo/4)*ky_2 + (passo/3)*ky_3 + (passo/8)*ky_4;
   
   kv_4 = g(h(i) + passo, y_aux, v_aux,planeta1_pos(i,:),planeta2_pos(i,:),planeta3_pos(i,:));
   zv = v(i,:) + passo*(7/24)*kv_1 + (passo/4)*kv_2 + (passo/3)*kv_3 + (passo/8)*kv_4;
   
   
   erro_y = norm(y_aux - zy);
   
   erro_v = norm(v_aux - zv);
   
   passo1 = passo * sqrt(epsilon/erro_y);

   passo2 = passo * sqrt(epsilon/erro_v);

   %entao recalculamos o y desta iteraçao com o nosso
   %passo________________________________________________________________________________________________________________________________________

   passo_aux1 = min([passo1,passo2,passo_max]);
   passo = max([passo_min,passo_aux1]);

   ky_2 = f(h(i)+passo/2, y(i,:) + passo * ky_1/2,v(i,:)+ passo * ky_1/2);

   ky_3 = f(h(i)+3*passo/4, y(i,:) + 3*passo * ky_2/4,v(i,:)+ 3*passo * ky_2/4);

   y_aux = y(i,:) + (2/9) * passo * ky_1 + (1/3) * passo * ky_2 + (4/9) * passo * ky_3;


   kv_2 = g(h(i)+passo/2, y(i,:) + passo * kv_1/2,v(i,:)+ passo * kv_1/2,planeta1_pos(i,:),planeta2_pos(i,:),planeta3_pos(i,:));

   kv_3 = g(h(i)+3*passo/4, y(i,:) + 3*passo * kv_2/4,v(i,:)+ 3*passo * kv_2/4,planeta1_pos(i,:),planeta2_pos(i,:),planeta3_pos(i,:));

   v_aux = v(i,:) + (2/9) * passo * kv_1 + (1/3) * passo * kv_2 + (4/9) * passo * kv_3;
   
   
   y = [y ; y_aux];
   v = [v ; v_aux];
   h(i+1) = h(i) + passo;
   passos(i+1) = passo;
   planeta1_pos(i+1,:) = planeta_pos(h(i+1),planeta1_centro,planeta1_r,planeta1_p,planeta1_tau);
   planeta2_pos(i+1,:) = planeta_pos(h(i+1),planeta2_centro,planeta2_r,planeta2_p,planeta2_tau);
   planeta3_pos(i+1,:) = planeta_pos(h(i+1),planeta3_centro,planeta3_r,planeta3_p,planeta3_tau);
   
   distancia1(i+1) = norm(y(i+1,:) - planeta1_pos(i+1,:));
   distancia2(i+1) = norm(y(i+1,:) - planeta2_pos(i+1,:));
   distancia3(i+1) = norm(y(i+1,:) - planeta3_pos(i+1,:));
   
   distancia4(i+1) = norm(planeta1_pos(i+1,:) - planeta2_pos(i+1,:));
   distancia5(i+1) = norm(planeta1_pos(i+1,:) - planeta3_pos(i+1,:));
   distancia6(i+1) = norm(planeta2_pos(i+1,:) - planeta3_pos(i+1,:));
   
   
   i = i+1;
   
end
y_final = y(end,:);
if plot1
subplot(2,2,1)
plot(planeta1_centro(1), planeta1_centro(2),'*', planeta2_centro(1), planeta2_centro(2),'*', planeta3_centro(1), planeta3_centro(2),'*', planeta1_pos(:,1), planeta1_pos(:,2),'red', planeta2_pos(:,1), planeta2_pos(:,2),'blue', planeta3_pos(:,1), planeta3_pos(:,2),'magenta', y(:,1),y(:,2))
legend('centro planeta 1', 'centro planeta 2', 'centro planeta 3', 'trajetória planeta 1', 'trajetória planeta 2', 'trajetória planeta 3', 'trajetória astronauta');

subplot(2,2,2)
plot(h,distancia1,'red', h, distancia2,'blue', h, distancia3,'magenta')
legend('distância entre astronauta e planeta 1', 'distância entre astronauta e planeta 2', 'distância entre astronauta e planeta 3');

subplot(2,2,3)
plot(h, distancia4,'red', h, distancia5,'blue', h, distancia6,'magenta')
legend('distância entre planetas 1 e 2', 'distância entre planetas 1 e 3', 'distância entre planetas 2 e 3');

subplot(2,2,4)
plot3(h, y(:,1), y(:,2), h, planeta1_pos(:,1), planeta1_pos(:,2),'red', h, planeta2_pos(:,1), planeta2_pos(:,2),'blue', h, planeta3_pos(:,1), planeta3_pos(:,2),'magenta')
legend('trajetória astronauta', 'trajetória planeta 1', 'trajetória planeta 2', 'trajetória planeta 3');
grid on
end
end