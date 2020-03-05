function [ A,V,rho,time ] = TrajectTimeBased( Cd,S,TempI,rhoI,Vi,Xi,m,dt )
% Rocketry at Virignia Tech
% Ishan Arora
%
% The following function calculates rocket trajectory using numerical
% integration over small time intervals. For each time interval
% acceleration and air density is assumed constant, then updated for the
% next time interval
%
% The function is used in the CdAvsCdRSolver.m program.

%*****Requires metric inputs*******

%% Inputs
  %TempI = intial temperature (K)
  %rhoI = intial rho (kg/m^3)
  %Vi = intial velocity (m/s)
  %Xi = intial altitude (m)
  %Cd = drag coeffcient
  %S = rocket reference area (m^2)
  %dt = time step (s)
  %m = rocket mass (kg)

%% Outputs
  %A = altiude (ft)
  %V = velocity vect (ft/s)
  %time = time vect (s)
  %rho = air density
  
% gravitational constant
g = 9.81; %m/s^2
%lapse rate for temperature change in the troposphere
lpRate = -0.0065; %k/m
%universal gas constant
R = 287.05; %N*m/kg*K
%loop conter
iter = 1;
%sets Velocity final equal to Vi
Vf = Vi; %m/s
%sets final position equal to intial position
Xf = Xi; %m
%imtialized rhoN and TempN
TempN = TempI; %K
rhoN = rhoI; %kg/m^3
%variable for time
time(iter) = 0; %sec

%while loop that runs while velocity final is greater than 0 (b4 burnout)
while Vf > 0
  %calculates acceleration
  a = -((g)+(((0.5*rhoN*(Vi^2)*Cd*S))/m)); %m/s^2

  %finds new velocity
  Vf = Vi+(a*dt); %m/s
  %finds new altiude
  Xf = ((0.5*a*(dt^2))+(Vi*dt))+Xi; %m
  %updates temp
  TempF = TempI+(lpRate*(Xf-Xi)); %K
  %updats rho
  rhoN = rhoI*((TempF/TempI)^(-((g/(lpRate*R))+1))); %units = kg/m^3
  %resets values
  Xi = Xf; %m
  Vi = Vf;%m/s
  TempI = TempF; %K
  rhoI = rhoN; %kg/m^3
  time(iter) = dt*iter; %sec
  %logs velocity and position values
  V(iter) = Vi; %m/s
  A(iter) = Xi; %m
  % logs rho values
  rho(iter) = rhoI;
  iter = iter+1;
end

%converts to U.S. units
V = V*3.28084; %ft/s
A = A.*3.28084; %ft

end
