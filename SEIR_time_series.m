clear all
clc
fixed_param;

[t,y1] = ode15s(@SEIR_model,1:1:5000,y0,[]);
S=y1(:,1:4:(4*n)); E=y1(:,2:4:(4*n));
I=y1(:,3:4:(4*n)); R=y1(:,4:4:(4*n));

plot(I, 'LineWidth',2)
xlabel('Time')
legend('I_1(t)','I_2(t)','I_3(t)','I_4(t)','I_5(t)')
title('(A)')