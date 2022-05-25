function [y_final] = planetasRK(d,v_0,passo,limite,plot1)
%malha temporal:
h = 0:passo:limite;

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

distancia1(1) = norm(y(1,:) - planeta1_pos(1,:));
distancia2(1) = norm(y(1,:) - planeta2_pos(1,:));
distancia3(1) = norm(y(1,:) - planeta3_pos(1,:));

distancia4(1) = norm(planeta1_pos(1,:) - planeta2_pos(1,:));
distancia5(1) = norm(planeta1_pos(1,:) - planeta3_pos(1,:));
distancia6(1) = norm(planeta2_pos(1,:) - planeta3_pos(1,:));

for i = 1:length(h)-1
   k_1 = passo * f(h(i),y(i,:),v(i,:));
   
   k_2 = passo * f(h(i)+passo/2, y(i,:) + k_1/2,v(i,:)+k_1/2);
   
   k_3 = passo * f(h(i)+passo/2, y(i,:) + k_2/2,v(i,:)+k_2/2);
    
   k_4 = passo * f(h(i)+passo, y(i,:) + k_3,v(i,:)+k_3);
   
   
   y_aux = y(i,:) + k_1/6 + k_2/3 + k_3/3 + k_4/6;
   y = [y ; y_aux];
   
   
   k_1 = passo * g(h(i),y(i,:),v(i,:),planeta1_pos(i,:),planeta2_pos(i,:),planeta3_pos(i,:));
   
   k_2 = passo * g(h(i)+passo/2, y(i,:) + k_1/2,v(i,:)+k_1/2,planeta1_pos(i,:),planeta2_pos(i,:),planeta3_pos(i,:));
   
   k_3 = passo * g(h(i)+passo/2, y(i,:) + k_2/2,v(i,:)+k_2/2,planeta1_pos(i,:),planeta2_pos(i,:),planeta3_pos(i,:));
    
   k_4 = passo * g(h(i)+passo, y(i,:) + k_3,v(i,:)+k_3,planeta1_pos(i,:),planeta2_pos(i,:),planeta3_pos(i,:));
   
   v_aux = v(i,:) + k_1/6 + k_2/3 + k_3/3 + k_4/6;
   
   v = [v ; v_aux];
   
   planeta1_pos(i+1,:) = planeta_pos(h(i+1),planeta1_centro,planeta1_r,planeta1_p,planeta1_tau);
   planeta2_pos(i+1,:) = planeta_pos(h(i+1),planeta2_centro,planeta2_r,planeta2_p,planeta2_tau);
   planeta3_pos(i+1,:) = planeta_pos(h(i+1),planeta3_centro,planeta3_r,planeta3_p,planeta3_tau);
   
   distancia1(i+1) = norm(y(i+1,:) - planeta1_pos(i+1,:));
   distancia2(i+1) = norm(y(i+1,:) - planeta2_pos(i+1,:));
   distancia3(i+1) = norm(y(i+1,:) - planeta3_pos(i+1,:));
   
   distancia4(i+1) = norm(planeta1_pos(i+1,:) - planeta2_pos(i+1,:));
   distancia5(i+1) = norm(planeta1_pos(i+1,:) - planeta3_pos(i+1,:));
   distancia6(i+1) = norm(planeta2_pos(i+1,:) - planeta3_pos(i+1,:));
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