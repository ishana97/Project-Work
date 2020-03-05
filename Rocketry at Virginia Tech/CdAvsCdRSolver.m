% Rocketry at Virginia Tech
%
% Ishan Arora Fall 2018
%
% Coefficient of Drag Required vs Coefficient of Drag Actual Solver
%
% The following program the following program plots required Cd vs actual rocket Cd over time.
%
% The following program is built for a rocket competing in the Spaceport America Cup Competiton. One of the 
% main features of the rocket is that it needs to be able to reach 10,000 ft. AGL with little margin. To do this the team
% implemented flaps (air brakes) which would deploy to slow the rocket down. The idea being the rocket would intially 
% overshoot the target altitude and then use the flaps to fine tune its final altitude. 
%
% The program is designed to simulate what would actually occur onboard. First the program calculated the initial 
% trajectory after the solid rocket motor burnout. Following this, is the main section of the program. A while loop 
% continuously runs while the predicted final velocity is greater at the end of the time step.  Next the program calculates 
% the required coefficient of drag (with initial conditions) to reach 10,000 ft. precisely. Next the program, will examine 
% the difference between the actual coefficient of drag and the required coefficient of drag to make a decision. The program 
% will continue to run and make commands while the final velocity is above 0. 

clc
clear

%% Inputs
% begins by intializing all the variables needed to find the Required Coefficient of Drag (CdR)

% intial altitude
Xi = 1190; %(m)
%intial velocity
Vi = 228; %m/s
%mass
m = 19.692; %kg
%air density
rho = 0.95; %kg/m^3
%rocket ref area
S = 0.0197; %m^2
%target apogee
apTar = 3048; %m
%lower Cd limit for root finding
CdL = 0.1;
%upper Cd limit for root finding
CdU = 1;
% error tolerance for finding CdS
errTol = 0.05; %percent
%max iterations fo finding CdS
maxIt = 50;
%value for actual Cd of rocket
CdA = 0.335;
%base drag coeffcient for when flaps are fully closed
CdB = 0.335;
%sets upper and lower CdA bound based on rocket max Cd and min Cd
CdAL = CdB; % Coefficient of Drag Actual Lower Bound
CdAU = 0.5;% Coefficient of Drag Actual Upper Bound
%value for the time step to be run in time based trajectory analysis
dt = 0.005; %sec

%variable to represent the amount the actual rocket drag coeffcient changes
%during a given flap adjustment
dCdA = 0.1;
%cource flap adjust
dCdAC = 0.2;
%fine flap adjust
dCdAF = 0.0088;

%variable to represent fin response time, ie time elapsed between command
%to adjust fins. does not include actual time elapsed while motor moves
%flaps into position, fine response time
finRespF = 0.1; %sec
%flap position response time course
finRespC = 0.5; %sec
finResp = finRespF; %sec
%percent error threshold in which to use fine fin adjustment
pFineErrThr = 10; %percent

%runs function to get CdR
[ CdR, pError, time ] = JJCdRSolver( Xi,Vi,m,rho,S,apTar,CdL,CdU,errTol,
maxIt );

% runs trajectory in in a single step at a time to allow for CdA to be adjusted

%Requires metric inputs
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
%%
TempI = 260; %K
rhoI = rho; %kg/m^3
% gravitational constant
g = 9.81; %m/s^2
%lapse rate for temperature change in the troposphere
lpRate = -0.0065; %k/m
%universal gas constant
R = 287.05; %N*m/kg*K
%loop counter
iter = 1;
%sets Velocity final equal to Vi
Vf = Vi; %m/s
%sets final position equal to intial position
Xf = Xi; %m
%initialized rhoN and TempN
TempN = TempI; %K

rhoN = rhoI; %kg/m^3
%variable for time
time(iter) = 0; %sec

%% Plots trajectory of intial CdA
% Tells us the scenario right after the rocket motor burns out.
[At,Vt,rhot,timet] = TrajectTimeBased( CdA,S,270,rho,Vi,Xi,m,dt);
figure (2)
plot(timet,At);
title('Trajectory with Intial CdA');
xlabel('Time (sec)');
ylabel('Altitude (ft)');
grid on;

%%
%Variable for velocity cut off point, where active drag closes
VcutOff = 25; %m/s
%variable to represent time elapsed since last fin response
timeElap = 0; %sec
%counts the number of chances for fin actuation there are
finActCountChance = 0;
%number of actual fin actuations
finActNum = 0;

%while loop that runs while velocity final is greater than 0 (before motor burnout)
while Vf > 0
  % Calculates CdR for given initial conditions
  [CdR(iter), pError, timeCDS ] = JJCdRSolver( Xi,Vi,m,rhoI,S,apTar,
  CdL,CdU,errTol,maxIt );
  % calculates acceleration
  a = -((g)+(((0.5*rhoN*(Vi^2)*CdA*S))/m)); %m/s^2
  %finds new velocity
  Vf = Vi+(a*dt); %m/s
  %finds new altiude
  Xf = ((0.5*a*(dt^2))+(Vi*dt))+Xi; %m
  %updates temperature
  TempF = TempI+(lpRate*(Xf-Xi)); %K
  %updats density
  rhoN = rhoI*((TempF/TempI)^(-((g/(lpRate*R))+1))); %units = kg/m^3
  %resets values - for next loop
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

  %INSERT ABILITY TO CHANGE CdA HERE
  if timeElap > finResp && Vi > VcutOff
    % Course fin response
    if ((abs(CdR-CdA)/CdR)*100) > pFineErrThr
      dCdA = dCdAC;
      finResp = finRespC;
    end
    % Fine fin response 
    if ((abs(CdR-CdA)/CdR)*100) < pFineErrThr
      dCdA = dCdAF;
      finResp = finRespF;
    end
    % Scenario if need higher Cd than what is physically possible. 
    if CdR(iter) > CdA && timeElap > 0 && CdA < CdAU
      CdA = CdA+dCdA;
      if CdA > CdAU
        CdA = CdAU;
      end
      timeElap = 0;
      finActNum = finActNum+1;
    end

    if CdR(iter) < CdA && timeElap > 0 && CdA > CdAL
      CdA = CdA - dCdA;
      if CdA < CdAL
        CdA = CdAL;
      end
      timeElap = 0;
      finActNum = finActNum+1;
    end

    timeElap = 0;
    finActCountChance = finActCountChance+1;
  end

  %commands active drag shut down, keeping same response rate
  if Vi < VcutOff && timeElap > finResp
    finResp = finRespC;
    if CdA > CdB %CdB is equivalent to CdAL
      CdA = CdA - dCdAC;
    else
      disp('Full Shut Down');
    end
    CdA = CdAL;
  end

  timeElap = 0;
  finActNum = finActNum+1;
  disp(time(iter));
  disp('Command Shut Down');
  end

  CdAU_Plot(iter) = CdAU;
  CdAL_Plot(iter) = CdAL;
  CdA_Plot(iter) = CdA;
  timeElap = timeElap+dt; %sec
  iter = iter+1;
end

%% Outputs - Figures & Plots -> ATTATCHED IN CdAvsCdROutput.pdf

% converts to U.S. units
V = V*3.28084; %ft/s
A = A.*3.28084; %ft

figure (1);
plot(time,CdR,'r');
hold on
plot(time,CdA_Plot,'b');
hold on
plot(time,CdAL_Plot,'k');
hold on
plot(time,CdAU_Plot,'c');
%hold on;
%plot(time,V,'c');
%hold on
%plot(time,A,'k');
legend('CdR','CdA','CdA Lower Limit','CdA Upper Limit');
xlabel('time (s)');
ylabel('Cd');
title('CdR vs CdA');
grid on;

figure (3)
plot(time,A);
title('Trajectory with CdA Actuations');
xlabel('Time (s)');
ylabel('Altitude (ft)');
grid on
display(A(end));
display(finActCountChance);
display(finActNum);
