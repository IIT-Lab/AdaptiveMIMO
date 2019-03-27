%============================================
% Curso de Comunicaciones Opticas Coherente
% Dr. Mario R. Hueda
% Contenido: variacion del SOP versus angulos "phi" y "chi"
% Version: 21_10_17
%============================================

%close all
clear
figure(1)
hold off
figure(2)
hold off

wt=[-pi:.01:pi];
phi=pi/3*.0;
chi=pi/8*1.0;
flag_chi=1;

N=50;
for n=0:1:N
    if flag_chi==1
        chi=n*pi/N;
    else
        phi=n*pi/N*2;
    end

    % Angulos de la elipse
    fv = [cos(chi); sin(chi)*exp(j*phi)];
    [alfa epsilon ar rs] = polellip(fv);    
    alfa=alfa/180*pi;
    epsilon=epsilon/180*pi;

    % Componentes temporales del campo electrico
    Ex=cos(chi)*cos(wt);
    Ey=sin(chi)*cos(wt+phi);
    
    % Generacion de fig. con el SOP
    figure(1)
    hold off
    h=plot(Ex,Ey,'r');
    axis([-1 1 -1 1])
    set(h,'Linewidth',2);
    set(h,'Markersize',14);
    set(gca,'XScale','lin','FontWeight','bold','FontSize',18,'YScale','lin','FontWeight','bold','FontSize',18);
    set(gca,'Linewidth',2);
    h=title(['$\chi=$' num2str(chi,'%1.2f') '/' '$\phi=$' num2str(phi,'%1.2f') '/' '$\alpha=$' num2str(alfa,'%1.2f') '/' '$\epsilon=$' num2str(epsilon,'%1.2f')]);
    set(h,'Interpreter','latex'); 
    set(h, 'FontName', 'Helvetica', 'FontWeight','Bold', 'FontSize', 20);

    % Representacion en la esfera de Poincare
    figure(2)
    % Vector de Jones
    fv = [cos(chi); sin(chi)*exp(j*phi)];
    % Esfera de Poincare
    stokes(fv)
    hold on 
end