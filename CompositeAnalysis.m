% Composite Analysis Code
% Ishan Arora
% Fall 2018
%
% The following program calculates composite properties of a flat
composite
% laminate

clear
clc
format long
%% Input Material Properties
E1=1.55*10^11;
E2=1.21*10^10;
E3=12.1*10^9;

Nu23=0.458;
Nu13=0.248;
Nu12=0.248;

G23=3.2*10^9;
G13=4.4*10^9;
G12=4.4*10^9;

Nu21=(E2*Nu12)/E1;

alpha_1 = -0.01800E-6;
alpha_2 = 24.3E-6;
alpha_3 = alpha_2;

beta_1 = 146E-6;
beta_2 = 4770E-6;
beta_3 = beta_2;
%% S - Reduced Compliance
% Calculates components of the reduced compliance matrix
S11 = 1/E1;
S12 = Nu12/E1;
S13 = Nu13/E1;
S21 = S12;
S22 = 1/E2;
S23 = -Nu23/E2;
S31 = S13;
S32 = S23;
S33 = 1/E3;
S44 = 1/G23;
S55 = 1/G13;
S66 = 1/G12;
 %% Q - Reduced Stiffness

% Calculates components of the reduced stiffness matrix
Q11 = E1/(1 -Nu12*Nu21);
Q12 = (E2*Nu12)/(1 - Nu12*Nu21);
Q22 = E2/(1 - Nu12*Nu21);
Q66 = G12;
%% Transformed S
% Calculates components of the transformed S matrix)
theta = input('What is the weave angle? (must input a vector and repeat
every number once)');
m = cosd(theta);
 n = sind(theta);
Sb11 = S11.*(m.^4) + (2.*S12 + S66).*(n.^2).*(m.^2) + S22.*(n.^4);
Sb12 = (S11 + S22 - S66).*(n.^2).*(m.^2) + S12.*((n.^4) + (m.^4));
Sb16 = (2.*S11 - 2.*S12 - S66).*(n.*m.^3) - (2.*S22 - 2.*S12 - S66).*(n.
^3).*m;
Sb22 = S11.*(n.^4) + (2.*S12 + S66).*(n.^2).*(m.^2) + S22.*(m.^4);
Sb26 = (2.*S11 - 2.*S12 - S66).*(n.^3).*m - (2.*S22 - 2.*S12 - S66).*n.*
(m.^3);
Sb66 = 2.*(2.*S11 + 2.*S22 - 4.*S12 - S66).*(n.^2).*(m.^2) + S66.*((n.^4)
+ (m.^4));
%% Transformed Q
% Calculates components of the transformed Q matrix)
Qb11 = Q11.*(m.^4) + 2.*(Q12 + 2.*Q66)*(n.^2).*(m.^2) + Q22.*(n.^4);
Qb12 = (Q11 + Q22 - 4.*Q66).*(n.^2).*(m.^2) + Q12.*((n.^4) + (m.^4));
Qb16 = (Q11 - Q12 - 2.*Q66).*n.*(m.^3) + (Q12 - Q22 + 2.*Q66).*(n.^3).*m;
Qb22 = Q11.*(n.^4) + 2.*(Q12 + 2.*Q66).*(n.^2).*(m.^2) + Q22.*(m.^4);
Qb26 = (Q11 - Q12 - 2.*Q66).*(n.^3).* m + (Q12 - Q22 + 2.*Q66).*n.*(m.
^3);
Qb66 = (Q11 + Q22 - 2.*Q12 - 2.*Q66).*(n.^2).*(m.^2) + Q66.*((n.^4) + (m.
^4));
ee=1;

% Builidn 3 dimensional Qb matrix per ply of composite material.
while ee<=length(Qb11)
Qb(:,:,ee) = [Qb11(ee) Qb12(ee) Qb16(ee);Qb12(ee) Qb22(ee) Qb26(ee); Qb16
(ee) Qb26(ee) Qb66(ee)];
ee=ee+1;
end

%% Coefficients of thermal & moisture deformation
alpha_x = (alpha_1)*(cosd(theta).^2)+(alpha_2)*(sind(theta).^2);
alpha_y = (alpha_1)*(sind(theta).^2)+(alpha_2)*(cosd(theta).^2);
alpha_xy = 2*(alpha_1-alpha_2)*cosd(theta).*sind(theta);
beta_x = (beta_1)*(cosd(theta).^2)+(beta_2)*(sind(theta).^2);
beta_y = (beta_1)*(sind(theta).^2)+(beta_2)*(cosd(theta).^2);
beta_xy = 2*(beta_1-beta_2)*cosd(theta).*sind(theta);

 %% Kirchhofff

thickness = input('What is the thickness?');
Plys = input('How many plys are there?');
deltaT = input('What is the change in temperature?');
deltaM = input('What is the change in moisture?');

delta = thickness/Plys;
for i=2:2:Plys*2 %This accounts for the bottom & top of every ply. Thus
the number of plys has effectively doubled to account for this.
z(1)=-(thickness/2);
z(i)=z(i-1)+delta;
z(i+1)=z(i);
end
z(end)=[];
%% Composite Lamination Theory
% Calcualte components of the ABD matrix
A=zeros(3,3);
B=zeros(3,3);
D=zeros(3,3);
for v=1:3
for j=1:3
for k=1:Plys
A(v,j)=A(v,j)+Qb(v,j,2*k)*(z(2*k)-z(2*k-1));
B(v,j)=B(v,j)+0.5*(Qb(v,j,2*k)*(z(2*k)^(2)-z(2*k-1)^(2)));
D(v,j)=D(v,j)+0.33*(Qb(v,j,2*k)*(z(2*k)^(3)-z(2*k-1)^(3)));
end
end
end

% Creates the ABD Matrix for the laminate:
ABD = [A(1,1) A(1,2) A(1,3) B(1,1) B(1,2) B(1,3);
A(1,2) A(2,2) A(2,3) B(1,2) B(2,2) B(2,3);
A(1,3) A(2,3) A(3,3) B(1,3) B(2,3) B(3,3);
B(1,1) B(1,2) B(1,3) D(1,1) D(1,2) D(1,3);
B(1,2) B(2,2) B(2,3) D(1,2) D(2,2) D(2,3);
B(1,3) B(2,3) B(3,3) D(1,3) D(2,3) D(3,3);];

% Creates the abd matrix for the laminate.
abd = inv(ABD);
a = [abd(1,1) abd(1,2) abd(1,3);
abd(2,1) abd(2,2) abd(2,3);
abd(3,1) abd(3,2) abd(3,3);];
b = [abd(1,4) abd(1,5) abd(1,6);
abd(2,4) abd(2,5) abd(2,6);
abd(3,4) abd(3,5) abd(3,6);];
d = [abd(4,4) abd(4,5) abd(4,6);
abd(5,4) abd(5,5) abd(5,6);
abd(6,4) abd(6,5) abd(6,6);];


%% Force and Moment Resultants
M=sym('M')
%Applied Loads <------------- Comment out if we do not know know applied
%loads
N_x = 0; N_y = 0; N_xy = 0; M_x = M; M_y = (2/5)*M; M_xy = (-1/8)*M;
%Given midsurface strain and curvature <------- Comment out below if
known
% rr=2;
% while rr<=length(sigma_x)
% N_x = N_x +((sigma_x(rr)+sigma_x(rr-1))/2)*(z(rr)-z(rr-1));
% N_y = N_y +((sigma_y(rr)+sigma_y(rr-1))/2)*(z(rr)-z(rr-1));
% N_xy = N_xy +((tau_xy(rr)+tau_xy(rr-1))/2)*(z(rr)-z(rr-1));
% fun1 = @(LL)((sigma_x(rr)+sigma_x(rr-1))/2)*LL;
% fun2 = @(LL)((sigma_y(rr)+sigma_y(rr-1))/2)*LL;
% fun3 = @(LL)((tau_xy(rr)+tau_xy(rr-1))/2)*LL;
% M_x = M_x + integral(fun1,z(rr-1),z(rr));
% M_y = M_y + integral(fun2,z(rr-1),z(rr));
% M_xy = M_xy + integral(fun3,z(rr-1),z(rr));
% rr = rr +2;
% end
% T7 = table([N_x],[N_y],[N_xy],[M_x],[M_y],[M_xy],'VariableNames',
{'N_x','N_y','N_xy','M_x','M_y','M_xy'});
% disp(T7)
%
%
%
% Initialize all Force and Moment Resulatants:
N_x_t =0; N_y_t =0; N_xy_t =0; M_x_t =0; M_y_t = 0;M_xy_t = 0;N_x_m = 0;
N_y_m =0; N_xy_m
= 0;M_x_m =0; M_y_m =0; M_xy_m =0;

for kk=2:2:2*Plys % Calculates the Force & Moment Resulatants at the top
of every ply.
N_x_t = N_x_t + (Qb(1,1,kk)*alpha_x(kk)+Qb(1,2,kk)*alpha_y(kk)+Qb(1,3,kk)
*alpha_xy
(kk)).*(z(kk)-z(kk-1));
N_y_t = N_y_t + (Qb(1,2,kk)*alpha_x(kk)+Qb(2,2,kk)*alpha_y(kk)+Qb(2,3,kk)
*alpha_xy
(kk)).*(z(kk)-z(kk-1));
N_xy_t = N_xy_t + (Qb(1,3,kk)*alpha_x(kk)+Qb(2,3,kk)*alpha_y(kk)+Qb(3,3,
kk)*alpha_xy
(kk)).*(z(kk)-z(kk-1));
M_x_t = M_x_t + 0.5*(Qb(1,1,kk)*alpha_x(kk)+Qb(1,2,kk)*alpha_y(kk)+Qb
(1,3,kk)
*alpha_xy(kk)).*(z(kk)^(2)-z(kk-1)^(2));


 M_y_t = M_y_t +0.5*(Qb(1,2,kk)*alpha_x(kk)+Qb(2,2,kk)*alpha_y(kk)+Qb(2,3,

kk)*alpha_xy
(kk)).*(z(kk)^(2)-z(kk-1)^(2));
M_xy_t = M_xy_t + 0.5*(Qb(1,3,kk)*alpha_x(kk)+Qb(2,3,kk)*alpha_y(kk)+Qb
(3,3,kk)
*alpha_xy(kk)).*(z(kk)^(2)-z(kk-1)^(2));
N_x_m = N_x_m + (Qb(1,1,kk)*beta_x(kk)+Qb(1,2,kk)*beta_y(kk)+Qb(1,3,kk)
*beta_xy(kk)).
*(z(kk)-z(kk-1));
N_y_m = N_y_m + (Qb(1,2,kk)*beta_x(kk)+Qb(2,2,kk)*beta_y(kk)+Qb(2,3,kk)
*beta_xy(kk)).
*(z(kk)-z(kk-1));
N_xy_m = N_xy_m + (Qb(1,3,kk)*beta_x(kk)+Qb(2,3,kk)*beta_y(kk)+Qb(3,3,kk)
*beta_xy
(kk)).*(z(kk)-z(kk-1));
M_x_m = M_x_m + 0.5*(Qb(1,1,kk)*beta_x(kk)+Qb(1,2,kk)*beta_y(kk)+Qb(1,3,
kk)*beta_xy
(kk)).*(z(kk)^(2)-z(kk-1)^(2));
M_y_m = M_y_m +0.5*(Qb(1,2,kk)*beta_x(kk)+Qb(2,2,kk)*beta_y(kk)+Qb(2,3,
kk)*beta_xy
(kk)).*(z(kk)^(2)-z(kk-1)^(2));
M_xy_m = M_xy_m + 0.5*(Qb(1,3,kk)*beta_x(kk)+Qb(2,3,kk)*beta_y(kk)+Qb
(3,3,kk)*beta_xy
(kk)).*(z(kk)^(2)-z(kk-1)^(2));
end

 % Final Values:
 NMTHat = [N_x_t,N_y_t,N_xy_t,M_x_t,M_y_t,M_xy_t]';
 NMMHat = [N_x_m,N_y_m,N_xy_m,M_x_m,M_y_m,M_xy_m]';
 NMT = NMTHat*deltaT;
 NMM = NMMHat*deltaM;

%% Midsurface strains and resultants
NM = [N_x+NMT(1)+NMM(1) N_y+NMT(2)+NMM(2) N_xy+NMT(3)+NMM(3) M_x+NMT(4)
+NMM(4) M_y+NMT(5)
+NMM(5) M_xy+NMT(6)+NMM(6)]';
mdsc = abd*NM;
disp('The thermal force and moment resultants')
T11 = table([N_x_t*deltaT],[N_y_t*deltaT],[N_xy_t*deltaT],[M_x_t*deltaT],
[M_y_t*deltaT],
[M_xy_t*deltaT],'VariableNames';{'Nxt','Nyt','Nxyt','Mxt','Myt','Mxyt'});
T12 = table([N_x_m*deltaM],[N_y_m*deltaM],[N_xy_m*deltaM],[M_x_m*deltaM],
[M_y_m*deltaM],
[M_xy_m*deltaM],'VariableNames';{'Nxm','Nym','Nxym','Mxm','Mym','Mxym'});

disp(T11) % Prints values
disp(T12)
eps_O_x =mdsc(1,1);
eps_O_y = mdsc(2,1);

gam_O_xy = mdsc(3,1);
kap_O_x = mdsc(4,1);
kap_O_y = mdsc(5,1);
kap_O_xy = mdsc(6,1);
%
% ****Comment out below if values are unknown and need to be calculated
% eps_O_x =0.00127495;
% eps_O_y = -0.00148638;
% gam_O_xy = 4.9e-20;
% kap_O_x = 1.59e-15;
% kap_O_y = -2.03e-15;
% kap_O_xy = -5.39e-16;

eps_x = eps_O_x + z*kap_O_x - alpha_x*deltaT - beta_x*deltaM;
eps_y = eps_O_y + z*kap_O_y - alpha_y*deltaT - beta_y*deltaM;
gam_xy = gam_O_xy + z*kap_O_xy - alpha_xy*deltaT - beta_xy*deltaM;

%Laminate Stresses
sigma_x = Qb11.*(eps_x) + Qb12.*(eps_y) + Qb16.*(gam_xy);
sigma_y = Qb12.*(eps_x) + Qb22.*(eps_y) + Qb26.*(gam_xy);
tau_xy = Qb16.*(eps_x) + Qb26.*(eps_y) + Qb66.*(gam_xy);

%Transformations
sigma_1 = (m.^2).*sigma_x + (n.^2).*sigma_y + (2.*m.*n).*tau_xy;
sigma_2 = (n.^2).*sigma_x + (m.^2).*sigma_y + (-2.*m.*n).*tau_xy;
tau_12 = (-m.*n).*sigma_x + (m.*n).*sigma_y + ((m.^2)-(n.^2)).*tau_xy;
eps_1 = (m.^2).*eps_x + (n.^2).*eps_y + (2.*m.*n).*0.5.*gam_xy;
eps_2 = (n.^2).*eps_x + (m.^2).*eps_y + (-2.*m.*n).*0.5.*gam_xy;
gam_12 = 2.*((-m.*n).*eps_x + (m.*n).*eps_y + ((m.^2) - (n.^2)).*0.5.
*gam_xy);

%% Effective Engineering Properties
H = thickness;
E_bar_x = (A(1,1)*A(2,2)-(A(1,2)^2))/(A(2,2)*H);
E_bar_y = (A(1,1)*A(2,2)-(A(1,2)^2))/(A(1,1)*H);
G_bar_xy = (A(3,3))/H;
v_bar_xy = (A(1,2))/(A(2,2));
v_bar_yx = (A(1,2))/(A(1,1));

%% Failure
fail_1_ten = 1500E6;
fail_1_comp = -1250E6;
fail_2_ten = 50E6;
fail_2_comp = -200E6;
fail_tau_12 = 100E6;

for i = 1:length(sigma_1) % Loop solves for critical values for Laminate
eq1(i) = sigma_1(i) == fail_1_ten;

eq2(i) = sigma_1(i) == fail_1_comp;
eq3(i) = sigma_2(i) == fail_2_ten;
eq4(i) = sigma_2(i) == fail_2_comp;
eq5(i) = tau_12(i) == fail_tau_12;
M1Tension(i) =(solve(eq1(i), M));
M1Compression(i) = solve(eq2(i),M);
M2Tension(i) = solve(eq3(i),M);
M2Compression(i) = solve(eq4(i),M);
MShear12(i) = solve(eq5(i),M);
end

%% Output - Graphs & Material Properties
ss =1;
zz =1;
disp('Q_bar Matrices:')
while ss<=Plys*2
disp('Layer')
disp(zz)
T6 = table([Qb(:,:,ss)],'VariableNames',{'Q_bar'});
disp(T6)
zz = zz+1;
ss=ss+2;
end
ss=1;
zz=1;
disp('alpha_x, alpha_y, and alpha_xy values:')

while ss<=Plys*2
disp('Layer')
disp(zz)
T7 = table(alpha_x(ss), alpha_y(ss), alpha_xy(ss),'VariableNames',
{'alpha_x','alpha_y','alpha_xy'});
disp(T7)
zz = zz + 1;
ss = ss + 2;
end

disp('ABD Matrix:')
disp(ABD)
disp('abd Matrix:')
disp(abd)
disp('Effective Laminate Properties')
T5 = table([E_bar_x],[E_bar_y],[G_bar_xy],[v_bar_xy],
[v_bar_yx],'VariableNames',
{'E_Bar_x','E_Bar_y','G_Bar_xy','V_bar_xy','V_bar_yx'});
disp(T5)
disp('Unit thermal force and moment resultants:')
T8 = table([N_x_t],[N_y_t],[N_xy_t],[M_x_t],[M_y_t],

[M_xy_t],'VariableNames',
{'Nxt_hat','Nyt_hat','Nxyt_hat','Mxt_hat','Myt_hat','Mxyt_hat'});
% T9 = table([N_x_m],[N_y_m],[N_xy_m],[M_x_m],[M_y_m],
[M_xy_m],'VariableNames',
{'Nxm_hat','Nym_hat','Nxym_hat','Mxm_hat','Mym_hat','Mxym_hat'});
disp(T8)
% disp(T9)
disp('The thermal force and moment resultants')
T11 = table([N_x_t*deltaT],[N_y_t*deltaT],[N_xy_t*deltaT],[M_x_t*deltaT],
[M_y_t*deltaT],
[M_xy_t*deltaT],'VariableNames',{'Nxt','Nyt','Nxyt','Mxt','Myt','Mxyt'});
T12 = table([N_x_m*deltaM],[N_y_m*deltaM],[N_xy_m*deltaM],[M_x_m*deltaM],
[M_y_m*deltaM],
[M_xy_m*deltaM],'VariableNames',{'Nxm','Nym','Nxym','Mxm','Mym','Mxym'});
disp(T11)
disp(T12)
disp('Resulting mid-surface strains and curvatures')
T10 = table([eps_O_x],[eps_O_y],[gam_O_xy],[kap_O_x],[kap_O_y],
[kap_O_xy],'VariableNames',
{'eps_O_x','eps_O_y','gam_O_xy','kap_O_x','kap_O_y','kap_O_xy'});
disp(T10)
Ply = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];% <----------Change based on how
many layers
disp('Stresses and Strains in XYZ Coordinate System:')
% T1 = table([Ply.'],[z.'],[sigma_x.'],[sigma_y.'],
[tau_xy.'],'VariableNames',
{'Ply','z','SigmaX','SigmaY','TauXY'});
% T2 = table([Ply.'],[z.'],[eps_x.'],[eps_y.'],
[gam_xy.'],'VariableNames',
{'Ply','z','EpsilonX','EpsilonY','GammaXY'});
% disp(T1);disp(T2)
% disp('Stresses and Strains in 123 Coordinate System:')
% T3 = table([Ply.'],[z.'],[sigma_1.'],[sigma_2.'],
[tau_12.'],'VariableNames',
{'Ply','z','Sigma1','Sigma2','Tau12'});
% T4 = table([Ply.'],[z.'],[eps_1.'],[eps_2.'],
[gam_12.'],'VariableNames',
{'Ply','z','Epslon1','Epsilon2','Gamma12'});
% disp(T3);disp(T4);
disp('M to cause failure in 1D tension')
ta1 = vpa(M1Tension,5);
[M1,I1] = min(abs((ta1)));
disp(ta1(I1))
disp('in layer:')
disp(I1)
disp('M to cause failure in 1D compression')

ta2 = vpa(M1Compression,5);
[M2,I2] = min(abs((ta2)));
disp(ta2(I2))
disp('in layer:')
disp(I2)
disp('M to cause failure in 2D tension')
ta3 = vpa(M2Tension,5);
[M3,I3] = min(abs((ta3)));
disp(ta3(I3))
disp('in layer:')
disp(I3)
disp('M to cause failure in 2D compression')
ta4 = vpa(M2Compression,5);
[M4,I4] = min(abs((ta4)));
disp(ta4(I4))
disp('in layer:')
disp(I4/2)
disp('M to cause failure 12D shear')

ta5 = vpa(MShear12,5);
[M5,I5] = min(abs((ta5)));
disp(ta5(I5))
disp('in layer:')
disp(I5/2)
