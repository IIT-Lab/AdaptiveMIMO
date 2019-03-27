%============================================
% Curso de Comunicaciones Opticas Coherente
% Dr. Mario R. Hueda
% Material relacionado: Pags. 21-24 
% Contenido: variacion del SOP versus angulos "phi" y "chi"
%============================================

close all

wt=[-pi:.01:pi];
phi=pi/2*1.0;
chi=pi/4*1.0;
flag_chi=1;
N=50;
for n=0:1.*N
    if flag_chi==1
        chi=n*pi/N/2+1e-6;
    else
        phi=n*pi/N/2+1e-6;
    end
    alfa=0.5*atan(tan(2*chi)*cos(phi));
    epsilon=0.5*asin(sin(2*chi)*sin(phi));

    Ex=cos(chi)*cos(wt);
    Ey=sin(chi)*cos(wt+phi);
    h=plot(Ex,Ey,'r');
    axis([-1 1 -1 1])
    set(h,'Linewidth',2);
    set(h,'Markersize',14);
    set(gca,'XScale','lin','FontWeight','bold','FontSize',18,'YScale','lin','FontWeight','bold','FontSize',18);
    set(gca,'Linewidth',2);
    h=title(['$\chi=$' num2str(chi,'%1.2f') ' / ' '$\phi=$' num2str(phi,'%1.2f') ' / ' '$\alpha=$' num2str(alfa,'%1.2f') ' / ' '$\epsilon=$' num2str(epsilon,'%1.2f')]);
    set(h,'Interpreter','latex'); 
    set(h, 'FontName', 'Helvetica', 'FontWeight','Bold', 'FontSize', 20);
    pause(.1)
end