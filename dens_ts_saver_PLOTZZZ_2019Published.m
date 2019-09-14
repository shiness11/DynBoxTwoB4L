% This script makes plots of all dens_ts_saver timeseries from constant phi runs
% 
% Run this script in the same folder with the output from running "NADWloop_2019Published.m"
% 
% It finds all variables with the prefix "dens_ts_saver" and makes a plot for each one.
%     Top panel of the three-panel figure is NADW flux versus time
%     Middle panel of three-panel figure has outcrop position for y1A, y2A, y3A, y1P, y2P, y3P versus time
%     Bottom panel of three-panel figure has basin depth for z1A, z2A, z3A, z1P, z2P, z3P versus time
% 
% Script also calculates various other parameters: residual transport (T_res) between each layer, diffusive transport (T_diff) between each 
%   layer, and zonal transport (T_zonal) within each layer. When running the model with a different vertical diffusivity profile, there are
%   commented sections that need to be uncommented (see below).
%
% SH 09/12/19


close all
%%
Timeseries = who('-regexp','dens_ts_saver');

rate_A = zeros(length(Timeseries),100);
rate_P = zeros(length(Timeseries),100);
ind = zeros(length(Timeseries),2);
time = zeros(length(Timeseries),2);
%%
for i = 1:length(Timeseries)
    
    figure(i)
    subplot(3,1,1)
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(14,:))'])
    subplot(3,1,2)
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(2,:)./2e6)'])
    hold on
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(3,:)./2e6)'])
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(4,:)./2e6)'])
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(5,:)./2e6)'])
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(6,:)./2e6)'])
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(7,:)./2e6)'])
    subplot(3,1,3)
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(8,:))'])
    hold on
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(9,:))'])
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(10,:))'])
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(11,:))'])
    eval(['plot(' Timeseries{i} '(1,:),' Timeseries{i} '(12,:))'])
    
    %% Equilibration time calculations
    
    for j = 1:100
        eval(['rate_A(i,j) = (' Timeseries{i} '(3,j+1) - ' Timeseries{i} '(3,j))/50;'])
        eval(['rate_P(i,j) = (' Timeseries{i} '(6,j+1) - ' Timeseries{i} '(6,j))/50;'])
    end
    
    [C1,i1] = max(abs(rate_A(i,:)));
    [C2,i2] = max(abs(rate_P(i,:)));
    
    ind(i,1) = i1;
    ind(i,2) = i2;
    
    eval(['time(i,1) = ' Timeseries{i} '(1,i1+1);'])
    eval(['time(i,2) = ' Timeseries{i} '(1,i2+1);'])
    
    %% Residual transport calculations
    
    eval(['y1A = ' Timeseries{i} '(2,:);'])
    eval(['y2A = ' Timeseries{i} '(3,:);'])
    eval(['y3A = ' Timeseries{i} '(4,:);'])
    eval(['y1P = ' Timeseries{i} '(5,:);'])
    eval(['y2P = ' Timeseries{i} '(6,:);'])
    eval(['y3P = ' Timeseries{i} '(7,:);'])

    % Definitions
    y_si_A = -1e6;    % Change this if sea ice positions is changed
    y_si_P = -1e6;
    F_0 = 2e-8;
    b1 = 5e-3;
    b2 = 1e-3;
    b3 = -1e-3;
    b4 = -2.5e-3;
    tauatm = 2*2.6e6;
    Deltabatm = 1.3e-3;
    ell = 2e6;
    Deltab1 = b1-b2; Deltab2 = b2 - b3; Deltab3 = b3 - b4;
    LA = 3e6;
    LP = 10e6;
    
    % Calculate the surface buoyancy fluxes
    F1Aa=0; F1Af=0; F1As=0; F2Aa=0; F2Af=0; F2As=0; F3Aa=0; F3Af=0; F3As=0; 
    F1Pa=0; F1Pf=0; F1Ps=0; F2Pa=0; F2Pf=0; F2Ps=0; F3Pa=0; F3Pf=0; F3Ps=0;
    
    % First add the atmospheric component 
    if y1A>y_si_A; bhat = Deltabatm.*( (y1A-y_si_A)./ell ) + (b1+b2)/2; F1Aa = hm./tauatm.*(bhat - (b1+b2)/2); end
    if y2A>y_si_A; bhat = Deltabatm.*( (y2A-y_si_A)./ell ) + (b1+b2)/2; F2Aa = hm./tauatm.*(bhat - (b2+b3)/2); end
    if y3A>y_si_A; bhat = Deltabatm.*( (y3A-y_si_A)./ell ) + (b1+b2)/2; F3Aa = hm./tauatm.*(bhat - (b3+b4)/2); end
    if y1P>y_si_P; bhat = Deltabatm.*( (y1P-y_si_P)./ell ) + (b1+b2)/2; F1Pa = hm./tauatm.*(bhat - (b1+b2)/2); end
    if y2P>y_si_P; bhat = Deltabatm.*( (y2P-y_si_P)./ell ) + (b1+b2)/2; F2Pa = hm./tauatm.*(bhat - (b2+b3)/2); end
    if y3P>y_si_P; bhat = Deltabatm.*( (y3P-y_si_P)./ell ) + (b1+b2)/2; F3Pa = hm./tauatm.*(bhat - (b3+b4)/2); end
    
    % Next add the freshwater flux
    siwidth = 2e5;
    F1Af = .4.*F_0.*exp( -((y1A-y_si_A)./siwidth).^2 ); 
    F2Af = .4.*F_0.*exp( -((y2A-y_si_A)./siwidth).^2 );
    F3Af = .4.*F_0.*exp( -((y3A-y_si_A)./siwidth).^2 );
    F1Pf = .4.*F_0.*exp( -((y1P-y_si_P)./siwidth).^2 );
    F2Pf = .4.*F_0.*exp( -((y2P-y_si_P)./siwidth).^2 );
    F3Pf = .4.*F_0.*exp( -((y3P-y_si_P)./siwidth).^2 );
    
    % Then add brine rejection  
    AABWratio = 1;  % Ratio defined as Atlantic/Pacific
    brwidth = 1.6e5;
    if y1A<y_si_A; F1As = -AABWratio              .*F_0.*exp( -((y1A+ell)./brwidth).^2 ); end
    if y2A<y_si_A; F2As = -AABWratio              .*F_0.*exp( -((y2A+ell)./brwidth).^2 ); end
    if y3A<y_si_A; F3As = -AABWratio              .*F_0.*exp( -((y3A+ell)./brwidth).^2 ); end
    if y1P<y_si_P; F1Ps = -(1+LA/LP*(1-AABWratio)).*F_0.*exp( -((y1P+ell)./brwidth).^2 ); end
    if y2P<y_si_P; F2Ps = -(1+LA/LP*(1-AABWratio)).*F_0.*exp( -((y2P+ell)./brwidth).^2 ); end
    if y3P<y_si_P; F3Ps = -(1+LA/LP*(1-AABWratio)).*F_0.*exp( -((y3P+ell)./brwidth).^2 ); end
    
    
    F1A = F1Aa + F1Af + F1As; F2A = F2Aa + F2Af + F2As; F3A = F3Aa + F3Af + F3As;
    F1P = F1Pa + F1Pf + F1Ps; F2P = F2Pa + F2Pf + F2Ps; F3P = F3Pa + F3Pf + F3Ps;
    
    T_res_1A(:,i) = F1A.*( 0     - y2A/2 )./Deltab1.*LA;
    T_res_2A(:,i) = F2A.*( y1A/2 - y3A/2 )./Deltab2.*LA;
    T_res_3A(:,i) = F3A.*( y2A/2 + ell/2 )./Deltab3.*LA;
    T_res_1P(:,i) = F1P.*( 0     - y2P/2 )./Deltab1.*LP;
    T_res_2P(:,i) = F2P.*( y1P/2 - y3P/2 )./Deltab2.*LP;
    T_res_3P(:,i) = F3P.*( y2P/2 + ell/2 )./Deltab3.*LP;
    
    %% Diffusive transport calculations
    
    eval(['z1A = ' Timeseries{i} '(8,:);'])
    eval(['z2A = ' Timeseries{i} '(9,:);'])
    eval(['z3A = ' Timeseries{i} '(10,:);'])
    eval(['z1P = ' Timeseries{i} '(11,:);'])
    eval(['z2P = ' Timeseries{i} '(12,:);'])
    eval(['z3P = ' Timeseries{i} '(13,:);'])
    
    Ly = 1e7;
    kappa = 1e-4;
    difftrans = -1600;       % Change this if kink depth changes
    diffthick = 300;         % Change this if kink width changes
    Deltakap = 8e-5;         % Change this if difference in kappa between upper and deep ocean changes
    
    
holdkap = kappa - Deltakap./2. * (tanh( (z1A-(difftrans))./diffthick ) + 1);
    T_diff_1A(:,i) = holdkap.*Ly.*LA./(0   - z1A);
holdkap = kappa - Deltakap./2. * (tanh( (z2A-(difftrans))./diffthick ) + 1);
    T_diff_2A(:,i) = holdkap.*Ly.*LA./(z1A - z2A);
holdkap = kappa - Deltakap./2. * (tanh( (z3A-(difftrans))/diffthick ) + 1);
    T_diff_3A(:,i) = holdkap.*Ly.*LA./(z2A - z3A);
holdkap = kappa - Deltakap./2. * (tanh( (z1P-(difftrans))./diffthick ) + 1);
    T_diff_1P(:,i) = holdkap.*Ly.*LP./(0   - z1P);
holdkap = kappa - Deltakap./2. * (tanh( (z2P-(difftrans))./diffthick ) + 1);
    T_diff_2P(:,i) = holdkap.*Ly.*LP./(z1P - z2P);
holdkap = kappa - Deltakap./2. * (tanh( (z3P-(difftrans))./diffthick ) + 1);
    T_diff_3P(:,i) = holdkap.*Ly.*LP./(z2P - z3P);
    
% !!! For runs with constant or linear kappa, the above code must be commented out and some of the code below should be uncommented (see main TwoBasinFourLayer_td_lp_2019Published script)
% For constant kappa profile uncomment one of the "holdkap" values and the code below and comment the code above under "Determine the diffusive transport in the basin"

 % holdkap = kappa - Deltakap;       % low
 % holdkap = kappa - Deltakap/2;     % medium
 % holdkap = kappa;                  % high

 % T_diff_1A(:,i) = holdkap*Ly*LA/(0   - z1A);
 % T_diff_2A(:,i) = holdkap*Ly*LA/(z1A - z2A);
 % T_diff_3A(:,i) = holdkap*Ly*LA/(z2A - z3A);
 % T_diff_1P(:,i) = holdkap*Ly*LP/(0   - z1P);
 % T_diff_2P(:,i) = holdkap*Ly*LP/(z1P - z2P);
 % T_diff_3P(:,i) = holdkap*Ly*LP/(z2P - z3P);

% For linear kappa profile uncomment the following code and comment the code above under "Determine the diffusive transport in the basin"

 % holdkap = -Deltakap/4000 * z1A + 2e-5;
 %     T_diff_1A(:,i) = holdkap*Ly*LA/(0   - z1A);
 % holdkap = -Deltakap/4000 * z2A + 2e-5;
 %     T_diff_2A(:,i) = holdkap*Ly*LA/(z1A - z2A);
 % holdkap = -Deltakap/4000 * z3A + 2e-5;
 %     T_diff_3A(:,i) = holdkap*Ly*LA/(z2A - z3A);
 % holdkap = -Deltakap/4000 * z1P + 2e-5;
 %     T_diff_1P(:,i) = holdkap*Ly*LP/(0   - z1P);
 % holdkap = -Deltakap/4000 * z2P + 2e-5;
 %     T_diff_2P(:,i) = holdkap*Ly*LP/(z1P - z2P);
 % holdkap = -Deltakap/4000 * z3P + 2e-5;
 %     T_diff_3P(:,i) = holdkap*Ly*LP/(z2P - z3P);

    %% Zonal transport
    U_0 = .05;
    T_zonal_1(:,i) = -U_0.*( (z1A.*y1A./2 - z1P.*y1P./2) - 0 );
    T_zonal_2(:,i) = -U_0.*( (z2A.*y2A./2 - z2P.*y2P./2) - (z1A.*y1A./2 - z1P.*y1P./2) );
    T_zonal_3(:,i) = -U_0.*( (z3A.*y3A./2 - z3P.*y3P./2) - (z2A.*y2A./2 - z2P.*y2P./2) ); 
    T_zonal_4(:,i) = -U_0.*( 0 - (z3A.*y3A./2 - z3P.*y3P./2) );

end