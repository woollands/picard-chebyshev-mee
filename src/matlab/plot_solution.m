% Codeto plot Trajectory and Hamiltonian
% Robyn Woollands

clear
% close all
clc

cd ../
load output.txt
cd matlab
% load rk1210Out.txt
% output = rk1210Out;
Period = 4.306316113361824e+04;

figure(1)
plot(output(:,1)./Period,output(:,2:4),'o-','Linewidth',2)
title('Position')
xlabel('Time (s)')
ylabel('Position (km)')
legend('x','y','z')

figure(2)
plot(output(:,1)./Period,output(:,5:7),'o-','Linewidth',2)
title('Velocity')
xlabel('Time (s)')
ylabel('Velocity (km/s)')
legend('vx','vy','vz')

figure(3)
semilogy(output(:,1)./Period,output(:,8),'ro-')
grid on
title('Hamiltonian')
xlabel('Time (Orbits)')
ylabel('| E - E_o | / | E_o |')
set(gca, 'FontName', 'Helvetica','FontSize',16)
