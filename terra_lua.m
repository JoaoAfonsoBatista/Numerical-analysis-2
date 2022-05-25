function [] = terra_lua(v_0,passo,limite,epsilon,passo_min,passo_max,velocidades)
%malha temporal:
d = [7,0];
h = [0];
passos = [passo];

terra_massa = 81 * ((9+2)/8) * 10^5;
terra_centro = [7,0];
terra_r = 0;
terra_p = 0;
terra_tau = 0;

lua_massa = ((9+2)/8) * 10^5;
lua_centro = [7,0];
lua_r = 384;
lua_p = pi/14;
lua_tau = 0;



planeta_pos = @(t, centro,r,p,tau)  centro + r * [cos(p*(t - tau)) , sin(p*(t - tau))];

terra_pos(1,:) = planeta_pos(0,terra_centro,terra_r,terra_p,terra_tau);
lua_pos(1,:) = planeta_pos(0,lua_centro,lua_r,lua_p,lua_tau);

f = @(t,y,v) v;

g = @(t,y,v,terra_pos,lua_pos) lua_massa * ((lua_pos-y) / (norm(lua_pos-y))^3) + terra_massa * ((terra_pos-y) / (norm(terra_pos-y))^3);

y(1,:) = terra_pos + d;
v(1,:) =  v_0;
i = 1;

distancia1(1) = norm(y(1,:) - terra_pos(1,:));
distancia2(1) = norm(y(1,:) - lua_pos(1,:));

while h(i) < limite
    %calculo de RK-4
   ky_1 = passo * f(h(i),y(i,:),v(i,:));
   
   ky_2 = passo * f(h(i)+passo/2, y(i,:) + ky_1/2,v(i,:)+ky_1/2);

   ky_3 = passo * f(h(i)+passo/2, y(i,:) + ky_2/2,v(i,:)+ky_2/2);
    
   ky_4 = passo * f(h(i)+passo, y(i,:) + ky_3,v(i,:)+ky_3);
   
   
   y_aux = y(i,:) + ky_1/6 + ky_2/3 + ky_3/3 + ky_4/6;
   

   kv_1 = passo * g(h(i),y(i,:),v(i,:),terra_pos(i,:),lua_pos(i,:));
   
   kv_2 = passo * g(h(i)+passo/2, y(i,:) + kv_1/2,v(i,:)+kv_1/2,terra_pos(i,:),lua_pos(i,:));
   
   kv_3 = passo * g(h(i)+passo/2, y(i,:) + kv_2/2,v(i,:)+kv_2/2,terra_pos(i,:),lua_pos(i,:));
    
   kv_4 = passo * g(h(i)+passo, y(i,:) + kv_3,v(i,:)+kv_3,terra_pos(i,:),lua_pos(i,:));
   
   v_aux = v(i,:) + kv_1/6 + kv_2/3 + kv_3/3 + kv_4/6;
   
%_________________________________________________________________________________________________________________________
%calculo de RK-3
   ky1_1 = ky_1 / passo;
   
   ky1_2 = ky_2 / passo;
   
   ky1_3 = f(h(i)+passo,y(i,:) - ky1_1*passo + 2*ky1_2*passo, v(i,:) - ky1_1*passo + 2*ky1_2*passo);
         
   y1_aux = y(i,:) + (passo/6)*(ky1_1 + 4*ky1_2 + ky1_3);
   
   
   kv1_1 = kv_1 / passo;
   
   kv1_2 = kv_2 / passo;
   
   kv1_3 = g(h(i)+passo,y(i,:) - kv1_1*passo + 2*kv1_2*passo, v(i,:) - kv1_1*passo + 2*kv1_2*passo,terra_pos(i,:),lua_pos(i,:));
         
   v1_aux = v(i,:) + (passo/6)*(kv1_1 + 4*kv1_2 + kv1_3);
   
   
   
   erro_y = norm(y_aux - y1_aux);
   
   erro_v = norm(v_aux - v1_aux);
   
   passo1 = passo * sqrt(epsilon/erro_y);

   passo2 = passo * sqrt(epsilon/erro_v);

   %entao recalculamos o y desta iteraçao com o nosso
   %passo________________________________________________________________________________________________________________________________________

   passo_aux1 = min([passo1,passo2,passo_max]);
   passo = max([passo_min,passo_aux1]);

   ky_1 = passo * f(h(i),y(i,:),v(i,:));
   
   ky_2 = passo * f(h(i)+passo/2, y(i,:) + ky_1/2,v(i,:)+ky_1/2);
   
   ky_3 = passo * f(h(i)+passo/2, y(i,:) + ky_2/2,v(i,:)+ky_2/2);
    
   ky_4 = passo * f(h(i)+passo, y(i,:) + ky_3,v(i,:)+ky_3);
   
   
   y_aux = y(i,:) + ky_1/6 + ky_2/3 + ky_3/3 + ky_4/6;
   

   kv_1 = passo * g(h(i),y(i,:),v(i,:),terra_pos(i,:),lua_pos(i,:));
   
   kv_2 = passo * g(h(i)+passo/2, y(i,:) + kv_1/2,v(i,:)+kv_1/2,terra_pos(i,:),lua_pos(i,:));
   
   kv_3 = passo * g(h(i)+passo/2, y(i,:) + kv_2/2,v(i,:)+kv_2/2,terra_pos(i,:),lua_pos(i,:));
    
   kv_4 = passo * g(h(i)+passo, y(i,:) + kv_3,v(i,:)+kv_3,terra_pos(i,:),lua_pos(i,:));
   
   v_aux = v(i,:) + kv_1/6 + kv_2/3 + kv_3/3 + kv_4/6;
   
   if not(isempty(velocidades))
       if round(velocidades(1,1),3) == round(h(i),3)
           v_aux = v_aux + [velocidades(1,2),velocidades(1,3)];
           velocidades = velocidades(2:end,:);
       end
   end
   
   y = [y ; y_aux];
   v = [v ; v_aux];
   h(i+1) = h(i) + passo;
   passos(i+1) = passo;
   terra_pos(i+1,:) = planeta_pos(h(i+1),terra_centro,terra_r,terra_p,terra_tau);
   lua_pos(i+1,:) = planeta_pos(h(i+1),lua_centro,lua_r,lua_p,lua_tau);
   
   distancia1(i+1) = norm(y(i+1,:) - terra_pos(i+1,:));
   distancia2(i+1) = norm(y(i+1,:) - lua_pos(i+1,:));
   i = i+1;
   
end


subplot(1,3,1)
plot(terra_centro(1), terra_centro(2),'*', lua_pos(:,1), lua_pos(:,2),'blue', y(:,1), y(:,2))
legend('centro Terra', 'trajetória Lua', 'trajetória astronauta');

subplot(1,3,2)
plot(h, distancia1,'red', h, distancia2,'blue')
legend('distância entre astronauta e Terra', 'distância entre astronauta e Lua');

subplot(1,3,3)
plot3(h,y(:,1),y(:,2), h,terra_pos(:,1),terra_pos(:,2),'red', h,lua_pos(:,1),lua_pos(:,2),'blue')
grid on
legend('trajetória astronauta', 'trajetória Terra', 'trajetória Lua');

end