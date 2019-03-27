%============================================
% Curso de Comunicaciones Opticas Coherente
% Dr. Mario R. Hueda
% Material relacionado: Pags. 21-24 
% Contenido: variacion del SOP versus angulos "phi" y "chi"
% Version: 09_10_17
%============================================

%close all
clear
figure(1)
hold off
figure(2)
hold off

wt=[-pi:.01:pi];
phi=pi/2*1.0;
chi=pi/8*1.0;
flag_chi=1;

N=50;
for n=0:1:2*N
    if flag_chi==1
        chi=n*pi/N/2;
    else
        phi=n*pi/N;
    end
    
    % Vector de Stokes
    S=[cos(2*chi) sin(2*chi)*cos(phi) sin(2*chi)*sin(phi)];

    % Coordenadas esfericas de S
    [varphi,tita,r] = cart2sph(S(1),S(2),S(3));
    
    % Angulos de la elipse
    alfa=0.5*varphi;
    epsilon=0.5*tita;

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