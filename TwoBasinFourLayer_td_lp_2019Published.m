% This is a time-dependent, two-basin model with four layers.  
% The code is based on the two-basin steady state code.
%
% AFT, SH, September 11, 2018
%
% DIAGNOSTICS TO SOLVE FOR:
% y1A, y2A, y3A, z1A, z2A, z3A, y1P, y2P, y3P, z1P, z2P, z3P
%
% EQUATIONS TO SOLVE:
% d(z1A)/dt = -T_diff_1A + T_NADW     - T_res_1A - T_zonal_1                        
% d(z2A)/dt = -T_diff_2A + phi*T_NADW - T_res_2A - T_zonal_1 - T_zonal_2            
% d(z3A)/dt = -T_diff_3A +            - T_res_3A - T_zonal_1 - T_zonal_2 - T_zonal_3
% d(z1P)/dt = -T_diff_1P              - T_res_1P + T_zonal_1                         
% d(z2P)/dt = -T_diff_2P              - T_res_2P + T_zonal_1 + T_zonal_2
% d(z3P)/dt = -T_diff_3P              - T_res_3P + T_zonal_1 + T_zonal_2 + T_zonal_3 
%
% d(y1A)/dt = T_wind_1A + T_eddy_1A + T_res_1A
% d(y2A)/dt = T_wind_2A + T_eddy_2A + T_res_2A
% d(y3A)/dt = T_wind_3A + T_eddy_3A + T_res_3A
% d(y1P)/dt = T_wind_1P + T_eddy_1P + T_res_1P
% d(y2P)/dt = T_wind_2P + T_eddy_2P + T_res_2P
% d(y3P)/dt = T_wind_3P + T_eddy_3P + T_res_3P
%
% Input:  tfin_yr = Length of simulation in years
%         dt_day  = Time step in days
%         save_freq = Frequency of saving data in years
%    FOR HYSTERESIS LOOP:
%         TNADWin = input NADW flux (in m^3 s^-1)
%         phi_in = input NADW density (between 0 and 1; if phi_in = 1, all
%           NADW goes to layer 3, if phi_in = 0, all NADW goes to layer 2)
%         y_si_in = input sea ice extent (between 0 and -ell)
%         AABWratio_in = input for relative strength of bottom water
%           formation in Atlantic versus Pacific (default = 1)
% 
% Output file has 14 rows: time in increments of save frequency (days)
%                          y1A (m)
%                          y2A (m)
%                          y3A (m)
%                          y1P (m)
%                          y2P (m)
%                          y3P (m)
%                          Z1A (m)
%                          Z2A (m)
%                          Z3A (m)
%                          Z1P (m)
%                          Z2P (m)
%                          Z3P (m)
%                          NADW flux (m^3 s^-1)
%
% SH 09/12/19

function [dens_ts_saver] = TwoBasinFourLayer_td_lp_2019Published(tfin_yr,dt_day,save_freq,TNADWin,phi_in,y_si_in,AABWratio_in)


% MODEL PARAMETERS
global kappa K Ly LA LP H ell U_0 tau0 F_0 TNADW y_si AABWratio hm rho0 f b1 b2 b3 b4 tauatm Deltabatm phi

kappa = 1e-4;       % diapycnal diff. [m^2 s^-1]
K = 1000;           % isopycnal diff. [m^2 s^-1]
Ly = 1e7;           % meridional basin extent [m]
LA = 3e6;           % Atlantic zonal extent [m]
LP = 10e6;          % Pacific zonal extent [m]
H = 4000;           % ocean depth [m]
ell = 2e6;          % ACC meridional extent [m]
U_0 = .05;          % ACC zonal velocity [m s^-1]
tau0 = .1;          % Reference wind stress [N m^-2]
F_0 = 2e-8;         % Reference surface buoyancy flux [m^2 s^-3]
TNADW = TNADWin;    % Transformation in NADW [m^3 s^-1]
y_si = y_si_in;     % Sea ice extent
AABWratio = AABWratio_in; % Ratio of bottom water formation in Atlabtic versus Pacific
hm = 100;           % Mixed layer depth [m]
rho0 = 1025;        % Reference density [kg m^-3]
f = -1e-4;          % Coriolis parameter [s^-1]
b1 = 5e-3;          % Buoyancy of the AAIW class [m s^{-2}]
b2 = 1e-3;          % Buoyancy of the UCDW/PDW class [m s^{-2}]
b3 = -1e-3;         % Buoyancy of the LCDW/NADW class [m s^{-2}]
b4 = -2.5e-3;       % Buoyancy of the AABW class [m s^{-2}]
tauatm = 2*2.6e6;   % Relaxation time, surface buoyancy forcing [s; 2.6e6 = 1 month]
Deltabatm = 1.3e-3; % Range of surface temperature/buoyancy; [m s^{-2}] *** Corresponds to 50 W/m^2 at N bdy ***
phi = phi_in;       % Partitioning of NADW between layers 2 and 3 (phi = 1 --> layer 3; phi = 0 --> layer 2)


% INITIAL CONDITIONS

load init.mat       % Initial conditions
y1A = init(1); y2A = init(2); y3A = init(3);
y1P = init(4); y2P = init(5); y3P = init(6);
z1A = init(7); z2A = init(8); z3A = init(9);
z1P = init(10);z2P = init(11);z3P = init(12);


% TIME STEPPING

tfin = tfin_yr*3.15e7;
dt = dt_day*8.64e4;
save_freq_s = save_freq*3.15e7;
saver = ceil(tfin/save_freq_s);

% FIRST TIME STEP

clear ttt
ttt = 0:dt:tfin; tl = length(ttt);
tl 

%clear dens_ts
dens_ts_saver = zeros(9,saver+1);
dens_ts_saver(1,1) = 0;
dens_ts_saver(2:13,1) = [y1A y2A y3A y1P y2P y3P z1A z2A z3A z1P z2P z3P];
dens_ts = [y1A y2A y3A y1P y2P y3P z1A z2A z3A z1P z2P z3P];
dens_ts = dens_ts + dt*rhs(0,dens_ts);

% TIME INTEGRATION (Using Runge Kutta 4 method)

for i=2:(tl-1)                              % calculation loop
    
    [kky1A kky2A kky3A kky1P kky2P kky3P kkz1A kkz2A kkz3A kkz1P kkz2P kkz3P] = rhs( ttt(i)     , dens_ts           );
    kk1 = [kky1A,kky2A,kky3A kky1P,kky2P,kky3P kkz1A,kkz2A,kkz3A kkz1P,kkz2P kkz3P];

    [kky1A kky2A kky3A kky1P kky2P kky3P kkz1A kkz2A kkz3A kkz1P kkz2P kkz3P] = rhs( ttt(i)+dt/2, dens_ts + dt*kk1/2);
    kk2 = [kky1A,kky2A,kky3A kky1P,kky2P,kky3P kkz1A,kkz2A,kkz3A kkz1P,kkz2P kkz3P];

    [kky1A kky2A kky3A kky1P kky2P kky3P kkz1A kkz2A kkz3A kkz1P kkz2P kkz3P] = rhs( ttt(i)+dt/2, dens_ts + dt*kk2/2);
    kk3 = [kky1A,kky2A,kky3A kky1P,kky2P,kky3P kkz1A,kkz2A,kkz3A kkz1P,kkz2P kkz3P];

    [kky1A kky2A kky3A kky1P kky2P kky3P kkz1A kkz2A kkz3A kkz1P kkz2P kkz3P] = rhs( ttt(i)+dt  , dens_ts + dt*kk3  );
    kk4 = [kky1A,kky2A,kky3A kky1P,kky2P,kky3P kkz1A,kkz2A,kkz3A kkz1P,kkz2P kkz3P];
    
    dens_ts = dens_ts + (1/6)*(kk1 + 2*kk2 + 2*kk3 + kk4)*dt;
    
    if mod(i,floor((tl-1)/saver))==0
        dens_ts_saver(1  ,i/(floor((tl-1)/saver))+1 ) = ttt(i)/3.1536e7;
        dens_ts_saver(2:13,i/(floor((tl-1)/saver))+1 ) = dens_ts;
        TNADW_t = TNADW;
        dens_ts_saver(14, i/(floor((tl-1)/saver))+1 ) = TNADW_t;
        percent_completed = round((i+1)/tl*100)
    end
    
%     kk1'
%     kk2'
%     kk3'
%     dens_ts(:,i+1)'
%     pause
    
end

end



% Calculate fluxes

function [rhs_y1A,rhs_y2A,rhs_y3A,rhs_y1P,rhs_y2P,rhs_y3P,rhs_z1A,rhs_z2A,rhs_z3A,rhs_z1P,rhs_z2P,rhs_z3P] = rhs(tt,yy);
global kappa K Ly LA LP H ell U_0 tau0 F_0 TNADW y_si AABWratio hm rho0 f b1 b2 b3 b4 tauatm Deltabatm phi


y1A = yy(1); y2A = yy(2); y3A = yy(3); y1P = yy(4);  y2P = yy(5);  y3P = yy(6);
z1A = yy(7); z2A = yy(8); z3A = yy(9); z1P = yy(10); z2P = yy(11); z3P = yy(12);
Deltab1 = b1-b2; Deltab2 = b2 - b3; Deltab3 = b3 - b4;

% Calculate the isopycnal slopes
sl0 = tau0/rho0/f/K;
sl1A = -z1A/y1A; 
sl2A = -z2A/y2A;
sl3A = -z3A/y3A;
sl1P = -z1P/y1P;
sl2P = -z2P/y2P;
sl3P = -z3P/y3P;


% Calculate the surface buoyancy fluxes
F1Aa=0; F1Af=0; F1As=0; F2Aa=0; F2Af=0; F2As=0; F3Aa=0; F3Af=0; F3As=0; 
F1Pa=0; F1Pf=0; F1Ps=0; F2Pa=0; F2Pf=0; F2Ps=0; F3Pa=0; F3Pf=0; F3Ps=0;

% Sea ice positions
y_si_A = y_si;
y_si_P = y_si;

% First add the atmospheric component (north of sea ice edge)
if y1A>y_si_A; bhat = Deltabatm*( (y1A-y_si_A)/ell ) + (b1+b2)/2; F1Aa = hm/tauatm*(bhat - (b1+b2)/2); end
if y2A>y_si_A; bhat = Deltabatm*( (y2A-y_si_A)/ell ) + (b1+b2)/2; F2Aa = hm/tauatm*(bhat - (b2+b3)/2); end
if y3A>y_si_A; bhat = Deltabatm*( (y3A-y_si_A)/ell ) + (b1+b2)/2; F3Aa = hm/tauatm*(bhat - (b3+b4)/2); end
if y1P>y_si_P; bhat = Deltabatm*( (y1P-y_si_P)/ell ) + (b1+b2)/2; F1Pa = hm/tauatm*(bhat - (b1+b2)/2); end
if y2P>y_si_P; bhat = Deltabatm*( (y2P-y_si_P)/ell ) + (b1+b2)/2; F2Pa = hm/tauatm*(bhat - (b2+b3)/2); end
if y3P>y_si_P; bhat = Deltabatm*( (y3P-y_si_P)/ell ) + (b1+b2)/2; F3Pa = hm/tauatm*(bhat - (b3+b4)/2); end

% Next add the freshwater flux
siwidth = 2e5;       % Width of the sea ice melt region
F1Af = .4*F_0*exp( -((y1A-y_si_A)/siwidth)^2 ); 
F2Af = .4*F_0*exp( -((y2A-y_si_A)/siwidth)^2 );
F3Af = .4*F_0*exp( -((y3A-y_si_A)/siwidth)^2 );
F1Pf = .4*F_0*exp( -((y1P-y_si_P)/siwidth)^2 );
F2Pf = .4*F_0*exp( -((y2P-y_si_P)/siwidth)^2 );
F3Pf = .4*F_0*exp( -((y3P-y_si_P)/siwidth)^2 );

% Then add brine rejection  
% AABWratio = 1;  % Ratio defined as Atlantic/Pacific
brwidth = 1.6e5;    % Width of the sea ice formation/brine rejection region
Fbmaxwidth = 1e4;   % Extra negative buoyancy forcing by southern boundary prevents model from running out of domain
if y1A<y_si_A; F1As = -AABWratio              *F_0*exp( -((y1A+ell)/brwidth)^2 ) - 4*F_0*exp( -((y1A+ell)/Fbmaxwidth)^2 ); end
if y2A<y_si_A; F2As = -AABWratio              *F_0*exp( -((y2A+ell)/brwidth)^2 ) - 4*F_0*exp( -((y2A+ell)/Fbmaxwidth)^2 ); end
if y3A<y_si_A; F3As = -AABWratio              *F_0*exp( -((y3A+ell)/brwidth)^2 ) - 4*F_0*exp( -((y3A+ell)/Fbmaxwidth)^2 ); end
if y1P<y_si_P; F1Ps = -(1+LA/LP*(1-AABWratio))*F_0*exp( -((y1P+ell)/brwidth)^2 ) - 4*F_0*exp( -((y1P+ell)/Fbmaxwidth)^2 ); end
if y2P<y_si_P; F2Ps = -(1+LA/LP*(1-AABWratio))*F_0*exp( -((y2P+ell)/brwidth)^2 ) - 4*F_0*exp( -((y2P+ell)/Fbmaxwidth)^2 ); end
if y3P<y_si_P; F3Ps = -(1+LA/LP*(1-AABWratio))*F_0*exp( -((y3P+ell)/brwidth)^2 ) - 4*F_0*exp( -((y3P+ell)/Fbmaxwidth)^2 ); end

% Total surface buoyancy forcing
F1A = F1Aa + F1Af + F1As; F2A = F2Aa + F2Af + F2As; F3A = F3Aa + F3Af + F3As;
F1P = F1Pa + F1Pf + F1Ps; F2P = F2Pa + F2Pf + F2Ps; F3P = F3Pa + F3Pf + F3Ps;

% Calculate residual transport
T_res_1A = F1A*( 0     - y2A/2 )/Deltab1*LA;
T_res_2A = F2A*( y1A/2 - y3A/2 )/Deltab2*LA;
T_res_3A = F3A*( y2A/2 + ell/2 )/Deltab3*LA;
T_res_1P = F1P*( 0     - y2P/2 )/Deltab1*LP;
T_res_2P = F2P*( y1P/2 - y3P/2 )/Deltab2*LP;
T_res_3P = F3P*( y2P/2 + ell/2 )/Deltab3*LP;


% Determine the RHS for the dy/dt equations
rhs_y1A = 1/hm*( -tau0/rho0/f + K*sl1A - T_res_1A/LA);
rhs_y2A = 1/hm*( -tau0/rho0/f + K*sl2A - T_res_2A/LA);
rhs_y3A = 1/hm*( -tau0/rho0/f + K*sl3A - T_res_3A/LA);
rhs_y1P = 1/hm*( -tau0/rho0/f + K*sl1P - T_res_1P/LP);
rhs_y2P = 1/hm*( -tau0/rho0/f + K*sl2P - T_res_2P/LP);
rhs_y3P = 1/hm*( -tau0/rho0/f + K*sl3P - T_res_3P/LP);


% Determine the diffusive transport in each basin
difftrans = -1600;      % Depth of kink in kappa profile
diffthick = 300;        % Vertical depth over which kappa changes from high in deep ocean to low in upper ocean
Deltakap = 8e-5;        % Change in kappa between deep and upper ocean
    
holdkap = kappa - Deltakap/2 * (tanh( (z1A-(difftrans))/diffthick ) + 1);
    T_diff_1A = holdkap*Ly*LA/(0   - z1A);
holdkap = kappa - Deltakap/2 * (tanh( (z2A-(difftrans))/diffthick ) + 1);
    T_diff_2A = holdkap*Ly*LA/(z1A - z2A);
holdkap = kappa*(1 - exp((-z3A-H)/100)) - Deltakap/2 * (tanh( (z3A-(difftrans))/diffthick ) + 1);  % kappa goes to zero exponentailly over bottom 100m to prevent model from running out of domain.
    T_diff_3A = holdkap*Ly*LA/(z2A - z3A);
holdkap = kappa - Deltakap/2 * (tanh( (z1P-(difftrans))/diffthick ) + 1);
    T_diff_1P = holdkap*Ly*LP/(0   - z1P);
holdkap = kappa - Deltakap/2 * (tanh( (z2P-(difftrans))/diffthick ) + 1);
    T_diff_2P = holdkap*Ly*LP/(z1P - z2P);
holdkap = kappa*(1 - exp((-z3P-H)/100)) - Deltakap/2 * (tanh( (z3P-(difftrans))/diffthick ) + 1);  % kappa goes to zero exponentailly over bottom 100m to prevent model from running out of domain.
    T_diff_3P = holdkap*Ly*LP/(z2P - z3P);

% For constant kappa profile uncomment one of the "holdkap" values and the code below and comment the code above under "Determine the diffusive transport in the basin"

 % holdkap = kappa - Deltakap;       % low
 % holdkap = kappa - Deltakap/2;     % medium
 % holdkap = kappa;                  % high

 % T_diff_1A = holdkap*Ly*LA/(0   - z1A);
 % T_diff_2A = holdkap*Ly*LA/(z1A - z2A);
 % T_diff_3A = holdkap*Ly*LA/(z2A - z3A);
 % T_diff_1P = holdkap*Ly*LP/(0   - z1P);
 % T_diff_2P = holdkap*Ly*LP/(z1P - z2P);
 % T_diff_3P = holdkap*Ly*LP/(z2P - z3P);

% For linear kappa profile uncomment the following code and comment the code above under "Determine the diffusive transport in the basin"

 % holdkap = -Deltakap/4000 * z1A + 2e-5;
 %     T_diff_1A = holdkap*Ly*LA/(0   - z1A);
 % holdkap = -Deltakap/4000 * z2A + 2e-5;
 %     T_diff_2A = holdkap*Ly*LA/(z1A - z2A);
 % holdkap = -Deltakap/4000 * z3A + 2e-5;
 %     T_diff_3A = holdkap*Ly*LA/(z2A - z3A);
 % holdkap = -Deltakap/4000 * z1P + 2e-5;
 %     T_diff_1P = holdkap*Ly*LP/(0   - z1P);
 % holdkap = -Deltakap/4000 * z2P + 2e-5;
 %     T_diff_2P = holdkap*Ly*LP/(z1P - z2P);
 % holdkap = -Deltakap/4000 * z3P + 2e-5;
 %     T_diff_3P = holdkap*Ly*LP/(z2P - z3P);

% Determine the NADW transport (constant).
TNADW_t = TNADW;

%%% Note for zonal transfers that positive indicates transfer from the
%%% Pacific into the Atlantic; negative is Atlantic --> Pacific!
T_zonal_1 = -U_0*( (z1A*y1A/2 - z1P*y1P/2) - 0 );
T_zonal_2 = -U_0*( (z2A*y2A/2 - z2P*y2P/2) - (z1A*y1A/2 - z1P*y1P/2) );
T_zonal_3 = -U_0*( (z3A*y3A/2 - z3P*y3P/2) - (z2A*y2A/2 - z2P*y2P/2) ); 
T_zonal_4 = -U_0*( 0 - (z3A*y3A/2 - z3P*y3P/2) );

% Calculate the RHS for the dz/dt equations
rhs_z1A = 1/Ly/LA*( -T_diff_1A + TNADW_t     - T_res_1A - T_zonal_1                         );
rhs_z2A = 1/Ly/LA*( -T_diff_2A + phi*TNADW_t - T_res_2A - T_zonal_1 - T_zonal_2             );
rhs_z3A = 1/Ly/LA*( -T_diff_3A               - T_res_3A - T_zonal_1 - T_zonal_2 - T_zonal_3 );
rhs_z1P = 1/Ly/LP*( -T_diff_1P               - T_res_1P + T_zonal_1 			            );
rhs_z2P = 1/Ly/LP*( -T_diff_2P               - T_res_2P + T_zonal_1 + T_zonal_2             );
rhs_z3P = 1/Ly/LP*( -T_diff_3P               - T_res_3P + T_zonal_1 + T_zonal_2 + T_zonal_3 );


end