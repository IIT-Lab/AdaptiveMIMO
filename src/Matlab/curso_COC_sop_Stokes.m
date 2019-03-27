%============================================
% Curso de Comunicaciones Opticas Coherente
% Dr. Mario R. Hueda
% Contenido: variacion del SOP versus angulos "alfa" y "epsilon"
% Version: 21_10_17
%============================================

%close all
clear
figure(1)
hold off
figure(2)
hold off

wt=[-pi:.01:pi];
alfa=-pi/10*1.0;
epsilon=pi/3*1.0;
flag_alfa=0;

N=50;
for n=0:1:N
    if flag_alfa==1
        alfa=n*pi/N;
    else
        epsilon=n*pi/N;
    end
    
    % Vector de Stokes
    S=[cos(2*epsilon)*cos(2*alfa) cos(2*epsilon)*sin(2*alfa) sin(2*epsilon)];


    % Obtencion del vector de Jones a partir del vector de Stokes 
    s_dot_sigma=[S(1), S(2)-j*S(3);S(2)+j*S(3), -S(1)];
    [V v]=eig(s_dot_sigma);
    m=find(abs(diag(v)-1)<1e-5); % Busqueda del autovector corresp. al autovalor=1
    J=V(:,m);
    J=exp(-j*angle(J(1)))*J; % Vector Jones
    
    % Angulos del vector de Jones
    phi=angle(J(2));
    chi=atan(abs(J(2))/real(J(1)));
    
    % Redifinicion de los angulos de la elipse a partir del Vector de Jones
    fv = [cos(chi); sin(chi)*exp(j*phi)];
    [alfa_2 epsilon_2 ar rs] = polellip(fv);  
    alfa_2=alfa_2/180*pi;
    epsilon_2=epsilon_2/180*pi;
    
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
    h=title(['$\chi=$' num2str(chi,'%1.2f') '/' '$\phi=$' num2str(phi,'%1.2f') '/' '$\alpha=$' num2str(alfa_2,'%1.2f') '/' '$\epsilon=$' num2str(epsilon_2,'%1.2f')]);
    set(h,'Interpreter','latex'); 
    set(h, 'FontName', 'Helvetica', 'FontWeight','Bold', 'FontSize', 20);

    % Representacion en la esfera de Poincare
    figure(2)
    % Vector de Jones
    fv = J;
    % Esfera de Poincare
    stokes(fv)
    hold on 
end
