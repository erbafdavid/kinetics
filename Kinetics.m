%%% Program Kinetics.m
%%% 
%%% Dynamics of a 2D gas of hard spheres 
%%% Pedagogical program to illustrate the concepts of thermodynamics 
%%% and fluid mechanics.
%%%
%%% Version 1.1, 23/01/2017
%%% D. Fabre
%%%
%%% USAGE :
%%% 1/ Set parameters (lines 25-48 of the program)
%%% 2/ Launch program
%%% The program will open the following figures :
%%%      Figure 1 : Instantaneous position of particules
%%%      Figure 2 : Statistics of particules in vertical 'bins'
%%%      Figure 3 : Force exerted on right boundary due to collision and its averaged value ("Pressure") 
%%% 3/ Displace the figures on screen to avoid overlap
%%% 4/ press "enter" to run the simulation
clear all;

%%% Version beta : thermostat (non validé)
%%% 
%%% PARAMETERS :
%%%

disp(' '); 
disp(' Program Kinetics.m ');
disp(' '); 
close all;

disp(' Which kind of simulation do you want to do ?')

disp(' [1] -> state equation of the gas')
disp(' [2] -> diffusion process with two populations of particules') 
disp(' [3] -> diffusion of vertical momentum')
disp(' [4] -> diffusion of temperature')
disp(' [5] -> shock wave');
disp(' [many other possibilities, under work]');

choice = input(' your choice ?') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choice 1 
if(choice==1)

% Dimensions
Lx = 1;
Ly = 1;
particulediameter = 0.01; % 0.01 is a good value if N is of order 1000
g= 0.;
% Properties of particules 1 (initially in left half domain)
Nparticules_1 = 1000; 
X1min = 0; X1max = 1; %min and max position of domain initially occupied
T1 = 2; % 'temperature' 
vq1 = sqrt(T1);
Um1 = 0; Vm1 = 0; % mean horizontal and vertical velocities
color_1 = 'b'; % color used for plots
% Properties of particules 2 (initially in left half domain)
Nparticules_2 = 000; 
X2min = 0.5; X2max = 1; %min and max position of domain initially occupied
T2 = 2; % 'temperature' 
vq2 = sqrt(T2); 
Um2 = 0; Vm2 = 0; % mean horizontal and vertical velocities
color_2 = 'r'; % color used for plots 
% Properties of the boundaries
BC_left = 'wall'; % allowed values are 'wall' and 'thermostat'
BC_updown = 'periodic' % allowed values are 'wall' and 'periodic'
BC_right = 'wall';
% in case left boundary is a thermostat, one can impose the mean velocity
% and temperature in an oscillating way
Twall_m = 0;
Twall_a = 0;
Uwall_m = 0;
Uwall_a = 0;
Vwall_m = 1;
Vwall_a = 0;
omega = 0.1;
X_fictive_surface= 0.5;
% plotting parameters
Nbin = 10; % number of vertical 'bins' for statistics
fig1 = 1; % set to 1 for figure of particule positions in figure 1
fig2 = 1; % set to 1 to plot 'pressure' force on right boundary and temperature in figure 2 
fig3 = 1; % set to 1 to plot statistics in vertical bins in figure 3;  
fig4 = 0; % set to 1 to plot density of particules 1 at successive instants in figure 4 
fig5 = 0; % set to 1 to plot vertical velocities at successive instants in figure 5
stepfig = 2; % number of time steps between plot refreshing (for fig. 2,3) (raise to accelerate)
stepfig4 = 20; % idem for fig. (4,5)
% Numerical parameters
Tmax = 1000; % end of simulation 
dt = 0.005; % time step ; dt should be smaller than diam/vq to resolve properly collisions
diameter_for_plots = 345*particulediameter; % diameter or particules in pixels for plots
Naveraging = 100000; % for time-averaging of pressure, etc... if large value average since beginning of simulation
pausemode = 0; % set to 0 for automatic time advance, and to 1 for manual (press enter to advance) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choice 1 
elseif(choice==2)

% Dimensions
Lx = 1;
Ly = 1;
particulediameter = 0.01; % 0.01 is a good value if N is of order 1000
g= 0.;
% Properties of particules 1 (initially in left half domain)
Nparticules_1 = 500; 
X1min = 0; X1max = 0.5; %min and max position of domain initially occupied
T1 = 2; % 'temperature' 
vq1 = sqrt(T1);
Um1 = 0; Vm1 = 0; % mean horizontal and vertical velocities
color_1 = 'b'; % color used for plots
% Properties of particules 2 (initially in left half domain)
Nparticules_2 = 500; 
X2min = 0.5; X2max = 1; %min and max position of domain initially occupied
T2 = 2; % 'temperature' 
vq2 = sqrt(T2); 
Um2 = 0; Vm2 = 0; % mean horizontal and vertical velocities
color_2 = 'r'; % color used for plots 
% Properties of the boundaries
BC_left = 'wall'; % allowed values are 'wall' and 'thermostat'
BC_updown = 'periodic' % allowed values are 'wall' and 'periodic'
BC_right = 'wall';
% in case left boundary is a thermostat, one can impose the mean velocity
% and temperature in an oscillating way
Twall_m = 0;
Twall_a = 0;
Uwall_m = 0;
Uwall_a = 0;
Vwall_m = 1;
Vwall_a = 0;
omega = 0.1;
X_fictive_surface= 0.5;
% plotting parameters
Nbin = 10; % number of vertical 'bins' for statistics
fig1 = 1; % set to 1 for figure of particule positions in figure 1
fig2 = 0; % set to 1 to plot 'pressure' force on right boundary and temperature in figure 2 
fig3 = 1; % set to 1 to plot statistics in vertical bins in figure 3;  
fig4 = 1; % set to 1 to plot density of particules 1 at successive instants in figure 4 
fig5 = 0; % set to 1 to plot vertical velocities at successive instants in figure 5
stepfig = 2; % number of time steps between plot refreshing (for fig. 2,3) (raise to accelerate)
stepfig4 = 20; % idem for fig. (4,5)
% Numerical parameters
Tmax = 1000; % end of simulation 
dt = 0.005; % time step ; dt should be smaller than diam/vq to resolve properly collisions
diameter_for_plots = 345*particulediameter; % diameter or particules in pixels for plots
Naveraging = 100000; % for time-averaging of pressure, etc... if large value average since beginning of simulation
pausemode = 0; % set to 0 for automatic time advance, and to 1 for manual (press enter to advance) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choice 3 
elseif(choice==3)

% Dimensions
Lx = 1;
Ly = 1;
particulediameter = 0.01; % 0.01 is a good value if N is of order 1000
g= 0.;
% Properties of particules 1 (initially in left half domain)
Nparticules_1 = 500; 
X1min = 0; X1max = 0.5; %min and max position of domain initially occupied
T1 = 2; % 'temperature' 
vq1 = sqrt(T1);
Um1 = 0; Vm1 = 1; % mean horizontal and vertical velocities
color_1 = 'b'; % color used for plots
% Properties of particules 2 (initially in left half domain)
Nparticules_2 = 500; 
X2min = 0.5; X2max = 1; %min and max position of domain initially occupied
T2 = 2; % 'temperature' 
vq2 = sqrt(T2); 
Um2 = 0; Vm2 = 0; % mean horizontal and vertical velocities
color_2 = 'b'; % color used for plots 
% Properties of the boundaries
BC_left = 'wall'; % allowed values are 'wall' and 'thermostat'
BC_updown = 'periodic' % allowed values are 'wall' and 'periodic'
BC_right = 'wall';
% in case left boundary is a thermostat, one can impose the mean velocity
% and temperature in an oscillating way
Twall_m = 0;
Twall_a = 0;
Uwall_m = 0;
Uwall_a = 0;
Vwall_m = 1;
Vwall_a = 0;
omega = 0.1;
X_fictive_surface= 0.5;
% plotting parameters
Nbin = 10; % number of vertical 'bins' for statistics
fig1 = 1; % set to 1 for figure of particule positions in figure 1
fig2 = 0; % set to 1 to plot 'pressure' force on right boundary and temperature in figure 2 
fig3 = 1; % set to 1 to plot statistics in vertical bins in figure 3;  
fig4 = 0; % set to 1 to plot density of particules 1 at successive instants in figure 4 
fig5 = 1; % set to 1 to plot vertical velocities at successive instants in figure 5
stepfig = 2; % number of time steps between plot refreshing (for fig. 2,3) (raise to accelerate)
stepfig4 = 20; % idem for fig. (4,5)
% Numerical parameters
Tmax = 1000; % end of simulation 
dt = 0.005; % time step ; dt should be smaller than diam/vq to resolve properly collisions
diameter_for_plots = 345*particulediameter; % diameter or particules in pixels for plots
Naveraging = 100000; % for time-averaging of pressure, etc... if large value average since beginning of simulation
pausemode = 0; % set to 0 for automatic time advance, and to 1 for manual (press enter to advance) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choice 4 
elseif(choice==4)

% Dimensions
Lx = 1;
Ly = 1;
particulediameter = 0.01; % 0.01 is a good value if N is of order 1000
g= 0.;
% Properties of particules 1 (initially in left half domain)
Nparticules_1 = 500; 
X1min = 0; X1max = 1; %min and max position of domain initially occupied
T1 = 1; % 'temperature' 
vq1 = sqrt(T1);
Um1 = 0; Vm1 = 0; % mean horizontal and vertical velocities
color_1 = 'b'; % color used for plots
% Properties of particules 2 (initially in left half domain)
Nparticules_2 = 500; 
X2min = 0.5; X2max = 1; %min and max position of domain initially occupied
T2 = 2; % 'temperature' 
vq2 = sqrt(T2); 
Um2 = 0; Vm2 = 0; % mean horizontal and vertical velocities
color_2 = 'b'; % color used for plots 
% Properties of the boundaries
BC_left = 'wall'; % allowed values are 'wall' and 'thermostat'
BC_updown = 'periodic' % allowed values are 'wall' and 'periodic'
BC_right = 'wall';
% in case left boundary is a thermostat, one can impose the mean velocity
% and temperature in an oscillating way
Twall_m = 0;
Twall_a = 0;
Uwall_m = 0;
Uwall_a = 0;
Vwall_m = 1;
Vwall_a = 0;
omega = 0.1;
X_fictive_surface= 0.5;
% plotting parameters
Nbin = 10; % number of vertical 'bins' for statistics
fig1 = 1; % set to 1 for figure of particule positions in figure 1
fig2 = 0; % set to 1 to plot 'pressure' force on right boundary and temperature in figure 2 
fig3 = 1; % set to 1 to plot statistics in vertical bins in figure 3;  
fig4 = 0; % set to 1 to plot density of particules 1 at successive instants in figure 4 
fig5 = 0; % set to 1 to plot vertical velocities at successive instants in figure 5
stepfig = 2; % number of time steps between plot refreshing (for fig. 2,3) (raise to accelerate)
stepfig4 = 20; % idem for fig. (4,5)
% Numerical parameters
Tmax = 1000; % end of simulation 
dt = 0.005; % time step ; dt should be smaller than diam/vq to resolve properly collisions
diameter_for_plots = 345*particulediameter; % diameter or particules in pixels for plots
Naveraging = 100000; % for time-averaging of pressure, etc... if large value average since beginning of simulation
pausemode = 0; % set to 0 for automatic time advance, and to 1 for manual (press enter to advance) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choice 5 
elseif(choice==5)

% Dimensions
Lx = 1;
Ly = 1;
particulediameter = 0.01; % 0.01 is a good value if N is of order 1000
g= 0.;
% Properties of particules 1 (initially in left half domain)
Nparticules_1 = 1000; 
X1min = 0; X1max = 0.5; %min and max position of domain initially occupied
T1 = 2; % 'temperature' 
vq1 = sqrt(T1);
Um1 = 0; Vm1 = 0; % mean horizontal and vertical velocities
color_1 = 'b'; % color used for plots
% Properties of particules 2 (initially in left half domain)
Nparticules_2 = 300; 
X2min = 0.5; X2max = 1; %min and max position of domain initially occupied
T2 = 2; % 'temperature' 
vq2 = sqrt(T2); 
Um2 = 0; Vm2 = 0; % mean horizontal and vertical velocities
color_2 = 'b'; % color used for plots 
% Properties of the boundaries
BC_left = 'wall'; % allowed values are 'wall' and 'thermostat'
BC_updown = 'periodic' % allowed values are 'wall' and 'periodic'
BC_right = 'wall';
% in case left boundary is a thermostat, one can impose the mean velocity
% and temperature in an oscillating way
Twall_m = 0;
Twall_a = 0;
Uwall_m = 0;
Uwall_a = 0;
Vwall_m = 0;
Vwall_a = 0;
omega = 0.1;
X_fictive_surface= 0.5;
% plotting parameters
Nbin = 10; % number of vertical 'bins' for statistics
fig1 = 1; % set to 1 for figure of particule positions in figure 1
fig2 = 0; % set to 1 to plot 'pressure' force on right boundary and temperature in figure 2 
fig3 = 1; % set to 1 to plot statistics in vertical bins in figure 3;  
fig4 = 0; % set to 1 to plot density of particules 1 at successive instants in figure 4 
fig5 = 0; % set to 1 to plot vertical velocities at successive instants in figure 5
stepfig = 2; % number of time steps between plot refreshing (for fig. 2,3) (raise to accelerate)
stepfig4 = 20; % idem for fig. (4,5)
% Numerical parameters
Tmax = 1000; % end of simulation 
dt = 0.005; % time step ; dt should be smaller than diam/vq to resolve properly collisions
diameter_for_plots = 345*particulediameter; % diameter or particules in pixels for plots
Naveraging = 100000; % for time-averaging of pressure, etc... if large value average since beginning of simulation
pausemode = 0; % set to 0 for automatic time advance, and to 1 for manual (press enter to advance) 

end
        
        
        

%%
%% End of parameter settings ; the following should not be modified
%%

%% Initialisations : initial positions and velocities of particules
Nparticules = Nparticules_1+Nparticules_2;

X(1:Nparticules_1) = X1min+(X1max-X1min)*rand(1,Nparticules_1);
Y(1:Nparticules_1) = rand(1,Nparticules_1);
U(1:Nparticules_1) = Um1+vq1/sqrt(2)*randn(1,Nparticules_1);
V(1:Nparticules_1) = Vm1+vq1/sqrt(2)*randn(1,Nparticules_1);

if(Nparticules_2>0) 
X(Nparticules_1+1:Nparticules) = X2min+(X2max-X2min)*rand(1,Nparticules_2);
Y(Nparticules_1+1:Nparticules) = rand(1,Nparticules_2);
U(Nparticules_1+1:Nparticules) = Um1+vq2/sqrt(2)*randn(1,Nparticules_2);
V(Nparticules_1+1:Nparticules) = Vm2+vq2/sqrt(2)*randn(1,Nparticules_2);
end

%% Various initialisations before simulation

% construction of two vectors used for drawing 'bins' in figure 1
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
    F_tab = [];P_tab = [];
    flux_momentum_x_tab=[];flux_momentum_y_tab=[]; 
    flux_momentum_x_av_tab=[];flux_momentum_y_av_tab=[]; 
    T_tab = [];
 
for i=1:Nbin
    x_bin(2*i-1)=(i-1)/10;
    x_bin(2*i) = i/10;
    x_bin_c(2*i-1)=(i-.5)/10;
    x_bin_c(2*i)=(i-.5)/10;
 end

    if(fig4==1)
    figure(4);
    hold off;
    end
    if(fig5==1)
    figure(5);
    hold off;
    end
    
    
%%
%% Beginning of time loop
%%

for it=1:Tmax/dt
    t = (it-1)*dt;
    
    % balistic part of the trajectory
    X = X+dt*U;%-g*dt^2/2;
    Y = Y+dt*V;
    
    % collisions with boundaries
    Force_right_boundary = 0;
    for i=1:Nparticules
        if(X(i)<0&U(i)<0) %% collision with left boundary
                if(strcmp(BC_left,'wall')==1) %left boundary treated as a wall 
                    U(i)=-U(i);
                    X(i)=-X(i);
                elseif(strcmp(BC_left,'thermostat')==1) % left boundary treated as a "thermostat"
                    X(i) = -X(i);
                    U(i) = (Uwall_m+Uwall_a*sin(omega*t))...
                    +sqrt(Twall_m+Twall_a*sin(omega*t))/sqrt(2)*abs(randn);
                    V(i) = (Vwall_m+Vwall_a*sin(omega*t))...
                    +sqrt(Twall_m+Twall_a*sin(omega*t))/sqrt(2)*(randn);
                end
        end
         if(X(i)>Lx&U(i)>0) %% collision with right boundary
            U(i)=-U(i);
            X(i)=Lx-(X(i)-Lx);
           Force_right_boundary = Force_right_boundary-2*U(i)/dt;
         end
         if(Y(i)<0)
            Y(i)=Y(i)+Ly;
         end
         if(Y(i)>1)
            Y(i)=Y(i)-Ly;
         end
    end
    
    % comptabilisation of flux across mid-surface
    flux_momentum_x = 0;flux_momentum_y = 0;
    for i = 1:Nparticules
    if((X(i)-dt*U(i)-X_fictive_surface)*(X(i)-X_fictive_surface)<0);
        if(U(i)>0) % particule crossing in positive direction
            flux_momentum_x = flux_momentum_x+U(i)/dt;
            flux_momentum_y = flux_momentum_y+(V(i)-0)/dt; 
        else  % particule crossing in negative direction
            flux_momentum_x = flux_momentum_x+(0-U(i))/dt;
            flux_momentum_y = flux_momentum_y+(0-V(i))/dt; 
        end
    end
    end
   
    % elastic collisions between particules    
    for i = 1:Nparticules
        for j=i+1:Nparticules
            distance = sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2);
            if( (distance<particulediameter)&(X(j)-X(i))*(U(j)-U(i))+(Y(j)-Y(i))*(V(j)-V(i)) < 0)
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

    %% statistics
    t_tab = [t_tab t]; 
    F_tab = [F_tab Force_right_boundary];
    P = mean(F_tab(max(1,it-Naveraging):it));
    P_tab = [P_tab P];
    
    flux_momentum_x_tab = [flux_momentum_x_tab flux_momentum_x];
    meanx = mean(flux_momentum_x_tab(max(1,it-Naveraging):it));
    flux_momentum_x_av_tab = [flux_momentum_x_av_tab meanx];
    
    flux_momentum_y_tab = [flux_momentum_y_tab flux_momentum_y];
    meany = mean(flux_momentum_y_tab(max(1,it-Naveraging):it));
    flux_momentum_y_av_tab = [flux_momentum_y_av_tab meany];
    
    
    Temp = mean(U.^2+V.^2);
    T_tab = [T_tab Temp];

    
    %% plots
    
   
     
    if(fig1==1)
        figure(1);
        hold off;
        plot([Lx Lx],[0 Ly]);
        hold on;
           %% Figure 1 : instantaneous positions of particules
        if(fig3==1)
           plot(Xcage,Ycage,'k--') % draws the bins
        end
        if (fig2==1)
            plot([Lx Lx],[0 Ly],'m');
            plot([X_fictive_surface X_fictive_surface],[0 Ly],'c--');
        end
        h = plot(X(1:Nparticules_1),Y(1:Nparticules_1),[color_1 'o'],'MarkerSize',diameter_for_plots);
        set(h,'MarkerFaceColor',color_1)
        h2 = plot(X(Nparticules_1+1:Nparticules),Y(Nparticules_1+1:Nparticules),[color_2 'o'],'MarkerSize',diameter_for_plots);
        set(h2,'MarkerFaceColor',color_2)
        axis square;
        axis([0 Lx 0 Ly])
        hold off;
    end
   
    if(mod(t/dt,stepfig)==0)     
    if (fig2==1)
    figure(2);
    subplot(6,1,1);
    plot(t_tab,F_tab,'m');
    title('Pressure Force on right boundary (instantaneous)');
%    axes([0 1 0 1]);
    
    subplot(6,1,2);
    plot(t_tab,P_tab,'m');
    title('Pressure Force on right boundary (averaged)');
    
    subplot(6,1,3);
    plot(t_tab,flux_momentum_x_tab,'c');
    title('Flux of horizontal momentum across mid-plane (instantaneous)');
    
    subplot(6,1,4);
    plot(t_tab,flux_momentum_x_av_tab,'c');
    title('Flux of horizontal momentum across mid-plane (averaged)');
    
    subplot(6,1,5);
    plot(t_tab,flux_momentum_y_tab,'g');
    title('Flux of vertical momentum across mid-plane (instantaneous)');
    
    subplot(6,1,6);
    plot(t_tab,flux_momentum_y_av_tab,'g');
    title('Flux of vertical momentum across mid-plane (averaged)');
    
%    axes([0 1 0 1]); 
    
%    subplot(3,1,3);
%    plot(t_tab,T_tab);
%    title('Temperature');

    
    end
    
    
    for i=1:Nbin
    
    Vmoy_bin(2*i-1) = sum(V(X>(i-1)/10&X<i/10))/sum(X>(i-1)/10&X<i/10);
    Vmoy_bin(2*i) = Vmoy_bin(2*i-1);
    Umoy_bin(2*i-1) = sum(U(X>(i-1)/10&X<i/10))/sum(X>(i-1)/10&X<i/10);
    Umoy_bin(2*i) = Umoy_bin(2*i-1);
    T_bin(2*i-1) = (sum( (V(X>(i-1)/10&X<i/10)-Vmoy_bin(2*i)).^2 ...
                           + (U(X>(i-1)/10&X<i/10)-Umoy_bin(2*i)).^2 )/sum(X>(i-1)/10&X<i/10));
    T_bin(2*i) =T_bin(2*i-1);
    N1_bin(2*i-1) = sum(X(1:Nparticules_1)>(i-1)/10&X(1:Nparticules_1)<i/10);
    N1_bin(2*i) = N1_bin(2*i-1);
    
    end
    
    if (fig3~=0)
    figure(3); 
    subplot(2,2,1);
    if(color_1==color_2) % one species
    histogram(X(1:Nparticules),[0:.1:1]);
    else % two species
    histogram(X(1:Nparticules_1),[0:.1:1]);
    hold on;
    histogram(X(Nparticules_1+1:Nparticules),[0:.1:1]);
    histogram(X,[0:.1:1],'Displaystyle','stairs');
    end
    title('Histogram of particule positions')
    axis([0 1 0 Nparticules/Nbin*1.2]);
    hold off;
    
    subplot(2,2,2);
    plot(x_bin,T_bin,'g');
    title(' "Temperature" (mean square velocity)')
    axis([0 1 0 max(T1,T2)*1.5]);
    hold off;

    subplot(2,2,3);
    plot(x_bin,Vmoy_bin,'c');
    axis([0 1 min(Vm1,Vm2)-0.2*max(vq1,vq2) max([Vm1 Vm2 Vwall_m+Vwall_a])+0.2*max(vq1,vq2) ]);
    title('mean vertical velocity');
    hold off;
    
    subplot(2,2,4);
    plot(x_bin,Umoy_bin,'m');
    title('mean horizontal velocity')
    axis([0 1 min(Um1,Um2)-0.2*max(vq1,vq2) max(Um1,Um2)+0.2*max(vq1,vq2) ]);
    hold off;
    
    end
  
    if(fig4==1&&mod(t/dt,stepfig4)==0)
    figure(4);
    plot(x_bin_c,N1_bin,'b');
%    axes([0 1 min(Vm1,Vm2)-0.2*max(vq1,vq2) max(Vm1,Vm2)+0.2*max(vq1,vq2) ]);
    title('repartition of blue particules at successive instants');
    hold on;
    end
    
    if(fig5==1&&mod(t/dt,stepfig4)==0)
    figure(5);
    plot(x_bin_c,Vmoy_bin,'c');
%    axes([0 1 min(Vm1,Vm2)-0.2*max(vq1,vq2) max(Vm1,Vm2)+0.2*max(vq1,vq2) ]);
    title('mean vertical velocity at successive instants');
    hold on;
    end
    
    
%% end of time step    
    if(it~=1&pausemode==0)
        pause(0.001); % laisse 1 miliseconde pour remettre à jour les figures.
    else
        pause % attend qu'on appuie sur entrée pour lancer l'animation
        disp(' Position the figures and press "enter" to launch the simulation : ');
    end
    end

end