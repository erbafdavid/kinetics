%%% Program gasdynamics
 
%%% Dynamics of a 2D gas
%%% this program is mainly made to study gas kinetics (can be developed for other purposes: fluide mechanics ,thermodynamics) :
%%% the parameters are are in the microscopique range

 %% Included simulations :
 
  %%% 1.Maxwell-Boltzmann distribution (perfect gas kinetics)
  %%% 2.Validation of ergodicite hypothesis
  %%% 3.Diffusion from center
  %%% 4.Diffusion of two species 
  %%% 5.Mean free path
  %%% 6.Shock wave
  %%% 7.VAN DER WAALS molecular interactions
  
%%% USAGE :

%%% 1/ Launch program
%%% 2/ choose simulation
%%% 3/ Displace the figures on screen to avoid overlap
%%% 4/ press "enter" to run the simulation

close all ;clear ; clc ;
disp(' ');

string = [' Program gasdynamics.m \n  Real GAZ kinetics inside box  \n ' ...
' choose between following simulation  \n  [1] --> MaxwellBoltzmann Kinetic theory distribution \n ' ...
' [2] --> Ergodicity hypothesis validation \n  [3] --> Diffusion of two population \n ' ...
' [4] --> Diffusion from center \n  [5] --> shock wave \n ' ...
' [6] --> Mean free path \n  [7] --> VAN DER WAALS molecular interactions \n ' ] ;

fprintf (string)

choice =input('select choice number \n ');


%% CHARACTERISTICS

%%% universal constant

m =1.6737236*10^-27; % mass of hydrogen atom
diameter = 1.06*10^-11 ;% diameter of hydrogen atom
kb=1.38064*10^-23; % Boltzmann constant
Na=6.022140*10^23; % Avogadro constant
R=kb*Na; % Perfect Gaz constant
g= 9.80665 ; %gravity  constant

%% propreties

% Dimensions

Lx = 10^-9; %% box length
Ly = 10^-9; %% box width

A=(Lx*Ly) ;% section of the box 
sigma=pi*diameter;  % hydrogen effective section

%% Properties of particles 1

T1 = 298;
Nparticles_1=400;
Nparticles_2 =0 ;vq2=0 ;Vm2=0;Um2=0 ; % for a simple simulation
X1min = 0; X1max = Lx; %min and max position of domain initially occupied
 
T2=0 ; %
vq1 =sqrt(3*kb*T1/m);
Um1 = 0; Vm1 = 0; % mean horizontal and vertical velocities

Nparticles=Nparticles_1+Nparticles_2;

%% Properties of particles 2 (initially in left half domain or on the center )

X_fictive_surface= 0.5;

%% plotting parameters
Nbin = 10; % number of vertical 'bins' for statistics
fig1 = 1; % set to 1 for figure of particule positions in figure 1
fig2 = 1; % set to 1 to plot 'pressure' force on right boundary and temperature in figure 2
fig3 = 1; % set to 1 to plot statistics in vertical bins in figure 3;
fig4 = 0; % set to 1 to plot density of particles 1 at successive instants in figure 4
fig5 = 0; % set to 1 to plot vertical velocities at successive instants in figure 5
fig6 = 0;
fig7 = 0 ;
stepfig = 1; % number of time steps between plot refreshing (for fig. 2,3) (raise to accelerate)
stepfig4 = 20; % idem for fig. (4,5)

%% Numerical parameters
dt = 3*10^-15; % time step dt should be smaller than diam/vq to resolve properly collisions
Tmax=400*dt ; % end of simulation
Naveraging = 1500; % for time-averaging of pressure, etc... if large value average since beginning of simulation
pausemode = 0; % set to 0 fhor automatic time advance, and to 1 for manual (press enter to advance)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of parameter settings ; the following should not be modified
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch choice 

    case {1,2,3,4,5,6}
        

%% Initialisations : initial positions and velocities of particles



%% ergodicity hypothesis

if (choice==2)
T=input('entrez une Periode T suffisament grande(>300) T ='); % Averaging period 
Tmax=T*dt ; % end of simulation
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add gravity
gravity=input('Would you like to add gravity effect [Y/N] ','s') ; % add gravity effect


%% CHARACTERISTICS PARTICLES 1 

X(1:Nparticles_1) = rand(1,Nparticles_1)*Lx;
Y(1:Nparticles_1) = rand(1,Nparticles_1)*Ly;
U(1:Nparticles_1) = Um1+(vq1/sqrt(2))*randn(1,Nparticles_1);
V(1:Nparticles_1) = Vm1+(vq1/sqrt(2))*randn(1,Nparticles_1);


%% CHARACTERISTICS PARTICLES 2
%% case of linear diffusion

if (choice ==3)

Nparticles_2=input('choose number particles species 2 de diffuse max:2000 = ');

Nparticles=Nparticles_1+Nparticles_2;

X2min = 0.5*Lx; X2max = Lx; %min and max position of domain initially occupied
T2 = 400; % 'temperature of particles species 2'
vq2 =sqrt(3*kb*T2/m);
Um2 = 0; Vm2 = 0; % mean horizontal and vertical velocities  

X(Nparticles_1+1:Nparticles) =X2min+(X2max-X2min)*rand(1,Nparticles_2);
Y(Nparticles_1+1:Nparticles) = rand(1,Nparticles_2)*Ly;
U(Nparticles_1+1:Nparticles) = Um1+vq2/sqrt(2)*randn(1,Nparticles_2);
V(Nparticles_1+1:Nparticles) = Vm2+vq2/sqrt(2)*randn(1,Nparticles_2);

end 


%% DIFFUSION FROM CENTER + %% SHOCK WAVE
if choice ==4 || choice==5
switch choice
    case 4

    Nparticles_2=round(Nparticles_1/10); %% low concentration of particles 
    radius = Lx/20; 
  
    case 5

   Nparticles_2=Nparticles_1;
   vq1=3*vq1;
   radius =(Lx/200); %%% high concentration of particles 
end

Nparticles = Nparticles_1+Nparticles_2;
 
center = [Lx/2,Ly/2]; % center coordinates of the circle [x0,y0]
 
%%% in case of diffusion or explosion 
%%% initial particles position in middle %% radial distribution  

angle = 2*pi*rand(Nparticles_2,1);
r = radius*sqrt(rand(Nparticles_2,1));
X(Nparticles_1+1:Nparticles) = (r.*cos(angle)+ center(1))';
Y(Nparticles_1+1:Nparticles) = (r.*sin(angle)+ center(2))';
U(Nparticles_1+1:Nparticles) = Um2+(vq2/sqrt(2))*randn(1,Nparticles_2);
V(Nparticles_1+1:Nparticles) = Vm2+(vq2/sqrt(2))*randn(1,Nparticles_2);
end



%% MEAN FREE PATH INIALISATION

if (choice==6) 
Tmax=300*dt;

X1=X;Y1=Y ;U1=U ;V1=V;
kr=round(Nparticles*rand(1)); %% picking a random particle
Lpx=zeros(1,Nparticles) ;Lpy=zeros(1,Nparticles); s=ones(1,Nparticles); n=1 ;
fig1=0;

end

%%%%%

%% Various initialisations before simulation

%% construction of two vectors used for drawing 'bins' in figure 1

for i=1:2:Nbin+1
Xcage(2*i-1) = (i-1)*Lx/Nbin;
Xcage(2*i) = Xcage(2*i-1);
Xcage(2*i+1) = i*Lx/Nbin;
Xcage(2*i+2) = Xcage(2*i+1);
Ycage(2*i-1) = -1;
Ycage(2*i) = Ly+1;
Ycage(2*i+1) = Ly+1;
Ycage(2*i+2) = -1;
end

t_tab = [];
Fx_tab = [];Fy_tab = [];P_tab = [];Tauxy_tab = [];
flux_momentum_x_tab=[];flux_momentum_y_tab=[];
flux_momentum_x_av_tab=[];flux_momentum_y_av_tab=[];
T_tab = []; time=[];


for i=1:Nbin
x_bin(2*i-1)=(i-1)*Lx/10;
x_bin(2*i) = i*Lx/10;
x_bin_c(2*i-1)=(i-.5)*Lx/10;
x_bin_c(2*i)=(i-.5)*Lx/10;
end



if(fig4==1)
figure(4);
hold off;
end
if(fig5==1)
figure(5);
hold off;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beginning of time loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for it=1:Tmax/dt
t = (it-1)*dt;

if (choice ==1 || choice==2 )

Ustat(it,:)=U; % Updated velocities of particles will be used for statistics
Vstat(it,:)=V; % Updated velocities of particles will be used for statistics
Vmodulus(it,:)=sqrt(U.^2+V.^2); % Updated velocity modulus of particles will be used for statistics

end

%% balistic part of the trajectory
X = X+dt*U;
 if (gravity == 'Y' )
Y = Y+dt*V-g*dt^2/2 ; % gravitational effect
 else
Y = Y+dt*V ;
 end
 
%% collisions with boundaries Perfect case

Fx_right_boundary = 0;Fy_right_boundary = 0;


for i=1:Nparticles

    if (choice==1 || choice==2 || choice ==3 || choice==4 || choice==5)

  if(X(i)<0 && U(i)<0) %% collision with left boundary
U(i)=-U(i);
X(i)=-X(i);
  end

  if(X(i)>Lx && U(i)>0) %% collision with right boundary

Uians = U(i);
Vians = V(i);
X(i)=Lx-(X(i)-Lx);
U(i)=-U(i);
Fx_right_boundary = Fx_right_boundary-m*(U(i)-Uians)/dt;
Fy_right_boundary = Fy_right_boundary+m*(V(i)-Vians)/dt;
  end

  if(Y(i)<0 && V(i)<0) %% collision with down boundary
V(i)=-V(i);
Y(i)=-Y(i);
  end

  if (Y(i)>Ly && V(i)>0) %% collision with top boundary
Y(i)=Ly-(Y(i)-Ly);
V(i)=-V(i);
  end
    end
    
%% collision with wall  periodicity %%%

if choice == 6 

    if(X(i)< 0) 
X(i)=X(i)+Lx;
   end

   if(X(i)>Lx) 
X(i)=X(i)-Lx;
   end

   if(Y(i)<0 ) 
Y(i)=Y(i)+Ly;
   end
   
   if (Y(i)> Ly) 
Y(i)=Y(i)-Ly;
   end
end
end


%% elastic collisions between particles

for i = 1:Nparticles
 for j=i+1:Nparticles
distance = sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2);
  
   if( (distance<= diameter) && (X(j)-X(i))*(U(j)-U(i))+(Y(j)-Y(i))*(V(j)-V(i)) < 0)
tx = Y(j)-Y(i);
ty = X(i)-X(j);
U1a = U(j);
U2a = U(i);
V1a = V(j);
V2a = V(i);
U(j) = (U1a*tx*tx+U2a*ty*ty+V1a*tx*ty-V2a*tx*ty)/(tx*tx+ty*ty);
U(i) = (U1a*ty*ty+U2a*tx*tx-V1a*tx*ty+V2a*tx*ty)/(tx*tx+ty*ty);
V(j) = (U1a*tx*ty-U2a*tx*ty+V1a*ty*ty+V2a*tx*tx)/(tx*tx+ty*ty);
V(i) = -(U1a*tx*ty-U2a*tx*ty-V1a*tx*tx-V2a*ty*ty)/(tx*tx+ty*ty);

   end
  end
end


  if (choice ==1|| choice ==3|| choice ==4 ||choice==5 )

%% comptabilisation of flux across mid-surface
flux_momentum_x = 0;flux_momentum_y = 0;

for i = 1:Nparticles

  if((X(i)-dt*U(i)-X_fictive_surface)*(X(i)-X_fictive_surface)<0)
      
   if(U(i)>0) % particule crossing in positive direction
flux_momentum_x = flux_momentum_x+U(i)/dt;
flux_momentum_y = flux_momentum_y+V(i)/dt;

   else  % particule crossing in negative direction
flux_momentum_x = flux_momentum_x-U(i)/dt;
flux_momentum_y = flux_momentum_y-V(i)/dt;

   end
   
  end
  
end

%% statistics

t_tab = [t_tab t];

% force exerted on right boundary

Fx_tab = [Fx_tab Fx_right_boundary];
P = mean(Fx_tab(max(1,it-Naveraging):it))/Lx;
P_tab = [P_tab P];

Fy_tab = [Fy_tab Fy_right_boundary];
Tauxy = mean(Fy_tab(max(1,it-Naveraging):it));
Tauxy_tab = [Tauxy_tab Tauxy];

% flux of momentum on mid plane

flux_momentum_x_tab = [flux_momentum_x_tab flux_momentum_x];
meanx = mean(flux_momentum_x_tab(max(1,it-Naveraging):it));
flux_momentum_x_av_tab = [flux_momentum_x_av_tab meanx];

flux_momentum_y_tab = [flux_momentum_y_tab flux_momentum_y];
meany = mean(flux_momentum_y_tab(max(1,it-Naveraging):it));
flux_momentum_y_av_tab = [flux_momentum_y_av_tab meany];


Temp = mean(U.^2+V.^2);
T_tab = [T_tab Temp];


%% plots

%% Figure 1 : instantaneous positions of particles
 %% Figure 1 : instantaneous positions of particles

 if(fig1==1)
        figure(1);
        hold off;
        plot([Lx Lx],[0 Ly]);
        hold on;
           
     %% Figure 3 : statistics in vertical bins   
 
     if(fig3==1)
           plot(Xcage,Ycage,'k--') % draws the bins
        end
        if (fig6+fig2>0)
            plot([Lx Lx],[0 Ly],'m');
            plot([X_fictive_surface X_fictive_surface],[0 Ly],'c--');
        end
        
h = plot(X(1:Nparticles_1),Y(1:Nparticles_1),'o','Color','red','MarkerFaceColor','blue','Markersize',5);
h2 = plot(X(Nparticles_1+1:Nparticles),Y(Nparticles_1+1:Nparticles),'o','Color','red','MarkerFaceColor','red','MarkerSize',4);

title(['Time : ',num2str(it),' (unit)'],'color','r')
set(figure(1),'position',[10 300 600 400])
        
        axis([0 Lx 0 Ly])
       
        set(figure(3),'position',[700 300 600 400])
    end

%% Figure 2 : Pressure on right boundary

if(mod(t/dt,stepfig)==0)

    if (fig2==1)
figure(2);

plot(t_tab,Fx_tab,'k:',t_tab,P_tab,'m');
title('Pressure Force on right boundary (instantaneous and averaged)');
xlabel('Time');ylabel('N/m');

    end
    
end


%% Figure 6 : STRESS on right boundary
if(mod(t/dt,stepfig)==0)

    if (fig6==1)
figure(6);
subplot(3,1,1);
plot(t_tab,Fy_tab,'k:',t_tab,Tauxy_tab,'m');
title('Tangential stress on right boundary (instantaneous and averaged)');

subplot(3,1,2);
plot(t_tab,flux_momentum_y_tab,'k:',t_tab,flux_momentum_y_av_tab,'c');
title('Flux of vertical momentum across mid-plane (instantaneous and averaged)');

subplot(3,1,1);
plot(t_tab,Fy_tab,'k:',t_tab,Tauxy_tab,'m');
title('Tangential stress on right boundary (instantaneous and averaged)');

subplot(3,1,3);
plot(t_tab,Tauxy_tab,'m',t_tab,flux_momentum_y_av_tab,'c');
title('averaged values for mid-plane flux and right boundary');
set(figure(6),'position',[500 40 600 250])

    end
end


%% Figure 3

for i=1:Nbin

Vcar_moy_bin(2*i-1) = sum(V(X>(i-1)*Lx/10&X<i*Lx/10).^2)/sum(X>(i-1)*Lx/10&X<i*Lx/10);
Vcar_moy_bin(2*i) = Vcar_moy_bin((2*i-1));
Ucar_moy_bin(2*i-1) = sum(U(X>(i-1)*Lx/10&X<i*Lx/10).^2)/sum(X>(i-1)*Lx/10&X<i*Lx/10);
Ucar_moy_bin(2*i) = Ucar_moy_bin(2*i-1);

Vmoy_bin(2*i-1) = sum(V(X>(i-1)*Lx/10&X<i*Lx/10))/sum(X>(i-1)*Lx/10&X<i*Lx/10);
Vmoy_bin(2*i) = Vmoy_bin(2*i-1);
Umoy_bin(2*i-1) = sum(U(X>(i-1)*Lx/10&X<i*Lx/10))/sum(X>(i-1)*Lx/10&X<i*Lx/10);
Umoy_bin(2*i) = Umoy_bin(2*i-1);

%% T_bin(2*i-1) = (sum( (V(X>(i-1)/10&X<i/10)-Vmoy_bin(2*i)).^2 ...
%                      + (U(X>(i-1)/10&X<i/10)-Umoy_bin(2*i)).^2 )/sum(X>(i-1)/10&X<i/10));
%%
T_bin(2*i-1)=((Vcar_moy_bin(2*i-1)+Ucar_moy_bin(2*i-1))*m)/(3*kb);
T_bin(2*i) =T_bin(2*i-1);
N1_bin(2*i-1) = sum(X(1:Nparticles_1)>(i-1)/10&X(1:Nparticles_1)<i/10);
N1_bin(2*i) = N1_bin(2*i-1);

end

%% Figure 3 :

if (fig3~=0)
figure(3);
subplot(2,2,1);


if (choice==3)  || (choice==4) || (choice==5 )
histogram(X(1:Nparticles_1),[0:Lx/10:Lx]);
hold on;
histogram(X(Nparticles_1+1:Nparticles),[0:Lx/10:Lx]);
histogram(X,[0:Lx/10:Lx],'Displaystyle','stairs');

else % one species

    histogram(X(1:Nparticles),[0:Lx/10:Lx]);
end

title('Histogram of particule positions')
axis([0 Lx 0 Nparticles/Nbin*1.2]);
hold off;

subplot(2,2,2);
plot(x_bin,T_bin,'g');
title(' "Temperature" (mean square velocity)')
axis([0 Lx 0 max(T1,T2)*1.5]);
hold off;

subplot(2,2,3);
plot(x_bin,Vmoy_bin,'c');
axis([0 Lx min([Vm1 Vm2 ])-0.2*max(vq1,vq2) max([Vm1 Vm2 ])+0.2*max(vq1,vq2) ]);
title('mean vertical velocity');
hold off;

subplot(2,2,4);
plot(x_bin,Umoy_bin,'m');
title('mean horizontal velocity')
axis([0 Lx min(Um1,Um2)-0.2*max(vq1,vq2) max(Um1,Um2)+0.2*max(vq1,vq2) ]);
hold off;

end

if(fig4==1&&mod(t/dt,stepfig4)==0)

    figure(4);
plot(x_bin_c,N1_bin,'b');
%    axes([0 1 min(Vm1,Vm2)-0.2*max(vq1,vq2) max(Vm1,Vm2)+0.2*max(vq1,vq2) ]);
title('repartition of blue particles at successive instants');
hold on;
end
%% figure 5

  if (fig5==1&&mod(t/dt,stepfig4)==0)
figure(5);
plot(x_bin_c,Vmoy_bin,'c');
%  axes([0 1 min(Vm1,Vm2)-0.2*max(vq1,vq2) max(Vm1,Vm2)+0.2*max(vq1,vq2) ]);
title('mean vertical velocity at successive instants');
hold on;

  end
  end
  
  %%% mean free path 
  
  if (choice==6)
      
  for i=1:Nparticles 
  
    if (U1(i)==U(i)) && (V1(i)==V(i) ) % if No collision

format longE
dX(i)=abs(X(i)-X1(i));  %% calculating X vector for each particle
dY(i)=abs(Y(i)-Y1(i)) ; %% ..........  Y .........................
Lpx(i)=Lpx(i)+dX(i); Lpy(i)=Lpy(i)+dY(i); % ****sum of small displacement*****
Lsimt(i,s(i))=sqrt((Lpx(i))^2+(Lpy(i))^2); %calculating distance
X1(i)=X(i);Y1(i)=Y(i) ;% updating positions
format short

    else %% if there's a "collison" then we select another path

format longE
s(i)=s(i)+1;
dX(i)=0 ;dY(i)=0;Lpx(i)=0;Lpy(i)=0;% reinitialisation for another path
U1(i)=U(i);V1(i)=V(i); % updating velocities
X1(i)=X(i);Y1(i)=Y(i) ;% updating positions
format short
    end
  end
  clc
  disp(it)
  end
  
  
  
  
  %% statistics to be used  for diffusion 
 if (choice==3) || choice==4 || (choice==5)
Xsta(it,:)=X(Nparticles_1+1:Nparticles)-Lx/2;
Ysta(it,:)=Y(Nparticles_1+1:Nparticles)-Ly/2; 
 end

%% end of time step
if choice ==2
    clc
disp(it)
end

if(it~=1 && pausemode==0)
    pause(0.001); % laisse 1 miliseconde pour remettre a jour les figures.
else
pause()
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%% POST TREATEMENT %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (choice == 1)
%% comparing with Boltzmann Kinetic Theory Probability distribution
%% histograms of vertical & horizental Velocities

figure (10)
hold on
edges = -(3/2)*max(max(abs(Ustat))):max(max(abs(Ustat)))/20:(3/2)*max(max(abs(Ustat)));

hu = histogram(Ustat,edges,'Normalization','pdf');
hv = histogram(Vstat,edges,'Normalization','pdf');
f1 =@(y) sqrt(m/(2*pi*kb*T1))*exp(-m*(y^2)/(2*kb*T1));
fplot(f1,[-(3/2)*max(max(abs(Ustat))) (3/2)*max(max(abs(Ustat)))],'LineWidth',1.5)

xlabel('Speed (units/time)')
ylabel('Probability Distribution')
legend('Simulation for x-velocity', 'Simulation for y-velocity',' Maxwell-Boltzmann  Distribution')
hold off
TO=298;
%% histograms of velocity modulus

figure(11)
hold on
edges = min(min(Vmodulus)):max(max(Vmodulus))/40: max(max(Vmodulus));
hvm = histogram(Vmodulus ,edges,'Normalization','pdf');


f2 = @(y) (4)*pi*((m/(2*pi*kb*T1))^(3/2))*(y^2)*exp(-m*(y^2)/(2*kb*T1));

fplot(f2,[0 3*max(max(Vmodulus))],'LineWidth',1.5)
xlabel('Speed (units/time)')
ylabel('Probability Distribution')
legend('Simulation', 'Maxwell-Boltzmann  Distribution')
hold off

end

 if (choice ==2)

%% this section is for testing the validity of ergodicity hypothesis
k=round(Nparticles*rand(1,5)); % on tir  10 une particule au hasard

for i=1:Tmax/dt
 for j=1:5
Vmtp(i,j)=mean(Vmodulus(1:i,k(j)));

 end
end

T=1:1:Tmax/dt;

figure(4)
plot(T,Vmtp,'LineWidth',1.1)
xlabel('Averaging Period')
ylabel('\fontsize{12}\color{blue} Statistical & temporal variation of velocity modulus')
title('\fontsize{14}\color{blue} Validation of ergodicity hypothesis' )
legend('particle A','particle B','particle C','particle D','particle E')
end


%% DIFFUSION FROM CENTER 

if (choice==4)  
figure(5)
T=1:1:(Tmax/dt)-1;

 %%%%% for diffusion 

    for j=1:Tmax/dt-1
Xmeanrootsq(j)=(rms(Xsta(j,:)));
Ymeanrootsq(j)=(rms(Ysta(j,:)));

    end

Xsqmean=Xmeanrootsq.^2;
Ysqmean=Xmeanrootsq.^2;

Rsq=(Xsqmean+Ysqmean) ; %%%% Mean Square Distance 
Rsqmax=(((Lx)^2)/6).*ones(1,length(T)); %% maximum mean square distance  % theorical %

%%%% Plots

p1=plot(T,Rsq,'r','LineWidth',1.5); xlabel('Time '),ylabel('\fontsize{12}\color{red} Mean square distance R²'),title('diffusion law');

hold on
p2=plot(T,Rsqmax,'k');
legend('\fontsize{12} simulation ' ,'\fontsize{12} maximum mean square distance')
grid on ; grid minor ;

end
 
if choice==5

    figure(5)
T=1:1:(Tmax/dt)-1;

%%% for shock wave visualisation 
for j=1:Tmax/dt-1
Xmean(j)=(mean(Xsta(j,:)));
Ymean(j)=(mean(Ysta(j,:)));

end

R=sqrt(Xmean.^2+Ymean.^2);% Mean Distance (Shock wave evolution )
P=loglog(T,R,'b --','LineWidth',1.5); xlabel('Log Time'),ylabel('\fontsize{14}\color{red} Log(Mean particle distance R)'),title('\fontsize{14}\color{blue} Detonation');

grid on ; grid minor ;

end

%% MEAN FREE PATH 
if (choice ==6)

for i =1:Nparticles
m(i)=sum(Lsimt(i,1:s(i)))/s(i) ; %  between 2 (the first collision until the last one 
end
Lsim=mean(m);

%%% free mean path theorical
Ltheo=A/(sqrt(2)*Nparticles*diameter) ;

fprintf(' MEAN FREE PATH theorique \n ')
disp(Ltheo)
fprintf(' MEAN FREE PATH simulation \n')
disp(Lsim)

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VAN DER WAALS GAZ SIMULATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 7
    
%%% propreties
%%% Dimensions

Xmin = 0; %% 
Xmax= Lx; %% 
Ymin=0;
Ymax=Ly;

%%%%%%%%% Properties of particles %%%%%%%%%%%
Nparticles=50;
mh=m;
sigma=diameter;    
output= '';
epsilon = 500*kb;      %material parameter
time = 6*10^-14;           %total time in seconds
dt =8*10^-17;             %time step in seconds
method = 'Verlet';       %integration method
borders = [Xmin Xmax Ymin Ymax];


%% Initialisations : initial positions and velocities of particles

X=zeros(Nparticles,time/dt+1);
Y=zeros(Nparticles,time/dt+1);
U=zeros(Nparticles,time/dt+1);
V=zeros(Nparticles,time/dt+1);

m=(ones(Nparticles,time/dt+1))*mh;

X(1:Nparticles,1) = Xmin+(Xmax-Xmin)*rand(1,Nparticles);
Y(1:Nparticles,1)= Ymin+(Ymax-Ymin)*rand(1,Nparticles);

U(1:Nparticles,1) =Um1+(vq1/sqrt(2))*randn(1,Nparticles);
V(1:Nparticles,1) =Vm1+(vq1/sqrt(2))*randn(1,Nparticles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beginning of time loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:time/dt+1
    for i=1:Nparticles
        for j=i+1:Nparticles %set force/energy interactions for each particle
            
               
                %compute the distance
                rij = sqrt((X(i,t)-X(j,t))^2+ (Y(i,t)-Y(j,t))^2);
               
                %compute unit normal vector, n(1)=n_x, n(2)=n_y
                n = [ X(i,t)-X(j,t) Y(i,t)-Y(j,t) ]/norm([ X(i,t)-X(j,t) Y(i,t)-Y(j,t) ]);
               
                %for x direction
 
                f_x(j) = n(1)*(-4*epsilon)*((6*sigma^6)/rij^7 - (12*sigma^12)/rij^13);
               
                %for y direction
                f_y(j) = n(2)*(-4*epsilon)*((6*sigma^6)/rij^7 - (12*sigma^12)/rij^13);
               
                %potential energy
                V_j(j) = 4*epsilon*((sigma/rij)^12-(sigma/rij)^6);
           
            end
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
   %%        Compute the position and velocity using Euler or Verlet
        if (strcmp('Euler',method) || t==1) %first timestap is always Euler  
            % Forward Euler x
            X(i,t+1)= X(i,t) + V(i,t)*dt;
            V(i,t+1) = V(i,t) + (sum(f_x)/m(i))*dt;
          
        
           
            % Forward Euler y
            U(i,t+1) = U(i,t) + (sum(f_y)/m(i))*dt;
            Y(i,t+1) = Y(i,t) + U(i,t)*dt;
           
       
            
           
        elseif(strcmp('Verlet',method))
            %Verlet algorithm x
            X(i,t+1) = -X(i,t-1) + 2*X(i,t) + sum(f_x)/m(i)*dt^2;
            V(i,t) = (X(i,t-1) - X(i,t+1))/(2*dt);
            
           
            %Verlet algorithm y
            Y(i,t+1) = -Y(i,t-1) + 2*Y(i,t) + sum(f_y)/m(i)*dt^2;
            U(i,t) = (Y(i,t-1) - Y(i,t+1))/(2*dt);
        end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
        %%% borders condition 
            if (X(i,t+1) <= borders(1)+diameter/2) || (X(i,t+1) >= borders(2)-diameter/2) %x
                X(i,t+1) = 2*X(i,t) - X(i,t+1);
                if (strcmp('Euler',method))
                    V(i,t+1) = -V(i,t+1);
                elseif(strcmp('Verlet',method))
                    V(i,t) = -V(i,t);
                end    
            end
            if (Y(i,t+1) <= borders(3)+diameter/2) || (Y(i,t+1) >= borders(4)-diameter/2) %y
                Y(i,t+1) = 2*Y(i,t) - Y(i,t+1);
                if (strcmp('Euler',method))
                    U(i,t+1) = -U(i,t+1);
                elseif(strcmp('Verlet',method))
                    U(i,t) = -U(i,t);
                end  
            end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%   Energy calculations   %%%%%%%%%%%%%%%   
        Vj(i,t) = sum(V_j); %Potential energy: sum of influence of other particles
        T(i,t) = 1/2*m(1)*(V(i,t)^2+U(i,t)^2); %Kinetic energy: 1/2mv^2 
    end
end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(1)
for it=1:time/dt+1 
  
    plot([Xmin Xmax],[Ymin Ymax]);
    axis([Xmin Xmax Ymin Ymax])
h = plot( X(:,it),Y(:,it),'o','Color','red','MarkerFaceColor','blue','Markersize',8);
pause(0.001)
set(gcf, 'Position', get(0, 'Screensize'));

end



%%%%%%%%%% Calculate and plot energy %%%%%%%%%%%%%%%%%%%%%%%%
for t=1:time/dt+1
    Vj_sum(t)=sum(Vj(:,t)) ;%sum the potential energy of all particles 
    T_sum(t) = sum(T(:,t)) ;%1/2.*m'.*(sum(u.^2)+sum(v.^2)); %sum from all particles: 0.5*m*v^2
end 



figure(2);

plot(0:dt:time,T_sum,'r')
hold on
plot(0:dt:time,Vj_sum,'b')
plot(0:dt:time,(T_sum(:)+Vj_sum(:)),'k')
title('The tota energy of the system','FontSize',7);
xlabel('Time [s]','FontSize',12);
ylabel('Energy [J]','FontSize',12);
legend('Kinetic energy','potential energy','total energy')
hold off

%%% Create file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp('', output) %check if output file is provided
    createfile(X, Y, V, U, output);
end
end