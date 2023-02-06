clc
clear
clear all

%% V-n Flight Envelope for retracted flap regime

W = 1633*9.81;   % Maximum take-off weight (MTOW - unit - Newton
S = 12.8;        % Surface area of wing, units-m^2
CL_max = 2.258;   % CLmax from tornado calculation-CL v alpha graph
nCL_max = -2.659; % -CLmax from tornado calculation-CL v alpha graph
Vc = 93.11;      % Maximum cruise velocity units-m/s
rho = 1.225;     % density at sea level (Kg/m^3)
V_d = 1.5.*Vc;   % Dive speed (m/s)

%% V-n coordinates for positive load factors
% n=L/w=((0.5*rho*S*CL)/W)*V^2
% Constant c = (0.5*rho*S*CL)/W
c = (0.5.*rho.*S.*CL_max)./W;

% V-n coordinates for stall limit - positive load factors
% max.load factor n = 4.4
n_a = 4.4;
V_a = sqrt(n_a./c); % Where V_a - design maneuvering speed

V1 = linspace(0,V_a);
n1 = c.*(V1.^2);

%% V-n coordinates for negative load factors
% Note: -n=-L/w=((0.5*rho*S*nCL_max)/W)*V^2
% Constant c2 = (0.5*rho*S*-nCL_max)/W
c2 = (0.5.*rho.*S.*nCL_max)./W;

% V-n coordinates for stall limit - negative load factors
% min.load factor n = -1.8
n_b = -1.8;
V_b = sqrt(n_b./c2); %Where V_b - negative maneuvering speed

V2 = linspace(0,V_b);
n2 = c2.*(V2.^2);

%% Plot1
plot(V1,n1, 'b'); hold on;
% V-n coordinates for structual damage limit - positive load factors
plot([V_a V_d],[n_a n_a],'k'); hold on;
plot([V_d V_d],[0 n_a],'k');hold on;
% V-n coordinates for stall limit - negative load factors
plot(V2,n2,'b'); hold on;
% V-n coordinates for structual damage limit - negative load factors
plot([V_b Vc],[n_b n_b],'k');hold on;
plot([Vc V_d],[n_b 0],'k'); grid on

%% V-n Flight Envelope for deployed flap regime

S = 12.8;          % Surface area of wing, units-m^2
CL_max1 = 3.5;    % CLmax from tornado calculation-CL v alpha graph
nCL_max1 = 0;      % -CLmax from tornado calculation-CL v alpha graph
Vf = 60.19;        % Maximum flap-down velocity units-m/s
V_D = 1.5.*Vf;     % Dive speed (m/s)

%% V-n coordinates for positive load factors (Flapped flight)
% n=L/w=((0.5*rho*S*CL)/W)*V^2
% Constant c = (0.5*rho*S*CL)/W
c3 = (0.5.*rho.*S.*CL_max1)./W;

% V-n coordinates for stall limit - positive load factors
% max. n = 2.0
n_c = 2.0;
V_A = sqrt(n_c./c3); % Where V_a - design maneuvering speed

V3 = linspace(0,V_A);
n3 = c3.*(V3.^2);

%% Plot2
plot(V3,n3, 'r'); hold on;
% V-n coordinates for structual damage limit - positive load factors
plot([V_A V_D],[n_c n_c],'k'); hold on;
plot([V_D V_D],[0 n_c],'k');hold on;
plot([0 V_D],[0 0], 'k'); grid on; hold on

%% Graph Annotation
plot([Vc Vc],[n_b n_a], 'g--'); hold on;  %Cruise velocity, Unflapped position
plot([Vf Vf],[0 n_c], 'g--');hold on    %Flap down 
Vn_1 = sqrt((2.*W)/(rho.*S.*CL_max1));
plot([Vn_1 V_d],[1 1], 'y--');
plot([Vn_1 Vn_1],[0 1], 'm--');
Vn_2 = sqrt((2.*W)/(rho.*S.*CL_max));
plot([Vn_2 Vn_2],[0 1], 'k--');
plot([V_a V_a],[0 n_a], 'k--');
plot([V_A V_A],[0 n_c], 'm--');
ylabel('Load Factor, n');
xlabel('Velocity (m/s)');
title('V-n Flight Envelope');
xlim([0,150]);
ylim([-3,5.5]);
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
