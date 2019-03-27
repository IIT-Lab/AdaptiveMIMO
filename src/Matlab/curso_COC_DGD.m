
close all
clear all
clc

w0 = 0.1*pi; % Frecuencia de Oscilación del Pulso1 de Banda Angosta 1
M = 200;      % Cantidad de puntos de los pulsos de BA no modulados
n = 0:M;     % Indice pulso no modulado

fn = (1-1.0*cos(2*pi.*n./M)); % Pulso de Banda Angosta no modulado
fn = fn/max(fn);

phi=pi/2*1.0;
chi=pi/4*1.0;


s_x = [zeros(1,M) cos(w0*n)*cos(chi) zeros(1,M)]; 
s_y = [zeros(1,M) cos(w0*n+phi)*sin(chi) zeros(1,M)]; 

f_s_x = [zeros(1,M) fn.*cos(w0*n)*cos(chi) zeros(1,M)]; 
f_s_y = [zeros(1,M) fn.*cos(w0*n+phi)*sin(chi) zeros(1,M)]; 

den=1;
t_x=s_x;
f_t_x=f_s_x;
z=[1:length(f_s_x)]/M;
for DGD=0:1:M-40
    num=[zeros(1,round(DGD)) 1];
    f_t_y=filter(num,den,f_s_y);
    t_y=filter(num,den,s_y);

    subplot 121
    plot3(z,f_t_x,0*f_t_x,z,0*f_t_y,f_t_y,'r');hold off
    ylabel('E_x');zlabel('E_y');xlabel('z');
    grid
    axis([1 4 -1 1 -1 1])

    subplot 122
    plot(t_x(M+DGD+1:2*M),t_y(M+DGD+1:2*M),'r');hold off
    axis([-1 1 -1 1])
    pause(.1)
    if DGD==0
        pause(1)
    end
    if DGD==1
        pause(2)
    end
end

h= sinc(.25*w0*(-50:50));
h=h/sum(h);
r=filter(h, 1, f_s_x.*(2*s_x)+f_s_y.*(2*s_y));

