clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')


q = 1.60217662e-19; 
electron_conc = 1e15*100^2; 
m0 = 9.10938356e-31; 
effective_m = 0.26*m0; 
Temperature = 300; 
Boltz_const = 1.38064852e-23; 
Thermal_v = sqrt(2*Boltz_const*Temperature/effective_m); 



widithx = 200e-9; % width
lengthy = 100e-9; % length
time_step = lengthy/Thermal_v/100;
Voltageinx = 0.1; % 

%question a
electricfield = Voltageinx/widithx;
fprintf('The electric field of the charge is %0.2f V/m.\n',electricfield);
%question b
force= electricfield*q;
fprintf('The force on each electron is %0.2E. \n ',force);
%question c  a=f/m
acc= force/effective_m;
fprintf('The acceleration on each electron is %f .\n',acc);
timestep=100;
t_mn = 0.2e-12; 
dt=1e-14;
BoundaryX = 200e-9;                    % X boundary
BoundaryY = 100e-9;                    % Y boundary
n=1000;
T=300;
C.kb = 1.3806504e-23;
C.m_0 = 9.10938215e-31;             % electron mass
Em = 0.26 * C.m_0;                    % Mass of the Electron
vel= sqrt(2*C.kb*T/Em);
Pscat = 1 - exp(-(dt/t_mn));
ppx=rand(n,1)*BoundaryX;
ppy=rand(n,1)*BoundaryY;
vx=randn(n,1)*vel/sqrt(2);
vy=randn(n,1)*vel/sqrt(2);
randomvalue= randi(n,[10,1]);

 J = zeros(1,timestep-1);
for i=2:timestep
  % pxold=ppx;
    %pyold=ppy;
    
    st=Pscat> rand(n,1);
    vx(st)= randn(sum(st),1)*vel/sqrt(2);
    vy(st)=randn(sum(st),1)*vel/sqrt(2);
    pxold=ppx;
    pyold=ppy;
    ppx=pxold+ vx*dt;
    ppy=pyold+vy*dt;
    m=ppx<0;
    m1=ppx>BoundaryX;
    ppx(m) = ppx(m)+ BoundaryX;
    pxold(m)=BoundaryX;

    ppx(m1)=ppx(m1)- BoundaryX;
     pxold(m1)=0;
    k=ppy<0;
    k1=ppy>BoundaryY;
    vy(k)=-vy(k);
    vy(k1)=-vy(k1);
    
    %set colour
    myColors = ['r' 'b' 'g' 'y' 'm' ];
     myColorTyp = j;
    for j=1:10
        subplot(3,1,1);
        %title('The 2-D plot of particle trajectories');
       plot([pxold(randomvalue(j)) ppx(randomvalue(j))], [pyold(randomvalue(j)) ppy(randomvalue(j))]);
       hold on
    end
    pause(0.1)

    title('The 2-D plot of particle trajectories');
   
   average=(mean(sqrt((vx.^2) + (vy.^2))));
 
    %Tavgold=TemperatureAvg;
   TemperatureAvg = ((average.^2) * Em)/(2 * C.kb);
   Tavgold=TemperatureAvg;
  
    %TAvgp1 = TemperatureAvg;
   subplot(3,1,2);
   plot([i-1 i],[Tavgold TemperatureAvg],'r');
   xlim([0 timestep]);
   %ylim([0 400]);
   pause(0.1)
   hold on
   title('The Average Temperature');
    % J(i) = numElectrons * mean((average)) * q;
   J(i) = n * mean((vx)) * q;

end
% 
%question d
%current formula is J = vnqNy

figure(2)

plot(linspace(2,timestep,timestep),J);
title('The Current Plot');
ylabel('J');
xlabel('TimeSteps');
%question e
newaverage = sqrt(vx.^2 + vy.^2);
TemperatureAvgnew = (Em * (newaverage.^2))./(2 * C.kb);
TheppX = linspace(min(ppx), max(ppx), 100);
TheppY = linspace(min(ppy), max(ppy), 50);
[X,Y] = meshgrid(TheppX, TheppY);
Tsurf = griddata(ppx,ppy,TemperatureAvgnew,X,Y);
figure(3)
surf(Tsurf);%Tsurf
title('Temperature Map');

