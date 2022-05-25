tic
[a1,b1,c1,r1] = cossenoRK();

[a2,b2,c2,r2] = cossenoBF();

[a3,b3,c3,r3] = cossenoBS();
ordem_de_convergencia_de_runge_kutta = r1
ordem_de_convergencia_de_Adams_Bashforth = r2
ordem_de_convergencia_de_Bogacki_shampine = r3

subplot(1,3,1)
plot(a1, b1,'.', a1, c1, 'red')
legend('logaritmo do erro em função do logaritmo de h', 'aproximação dos mínimos quadrados');

subplot(1,3,2)
plot(a2, b2,'.', a2, c2, 'red')
legend('logaritmo do erro em função do logaritmo de h', 'aproximação dos mínimos quadrados');

subplot(1,3,3)
plot(a3, b3,'.', a3, c3, 'red')
legend('logaritmo do erro em função do logaritmo de h', 'aproximação dos mínimos quadrados');
toc