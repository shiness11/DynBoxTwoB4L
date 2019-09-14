% This script iteratively runs the function "TwoBasinFourLayer_td_lp_2019Published.m" with varying NADW flux values
% 
% This example loop runs with phi_in = 1, y_si_in = -ell/2 and will reproduce the curve in Figure 4a
%     It uses the input file "Example_phi1_NADW4.mat"
%     
%     Output from this script is a set of variables with timeseries of:
%       time, y1A, y2A, y3A, z1A, z2A, z3A, y1P, y2P, y3P, z1P, z2P, z3P, and NADW flux.
%       Files are labeled "dens_ts_saver_ii" where "ii" is the loop iteration number, and they are saved in the current directory.
%             
% *** Note: this will take a long time to run!! (~24 hours) ***
% 
% SH 09/12/19

clear all

load Example_phi1_NADW4.mat            % .mat file with initial conditions to start loop
init(1) = dens_ts_saver_01(2,end);     % this is phi = 1, TNADW = 4 Sv, y_si = -ell/2
init(2) = dens_ts_saver_01(3,end);
init(3) = dens_ts_saver_01(4,end);
init(4) = dens_ts_saver_01(5,end);
init(5) = dens_ts_saver_01(6,end);
init(6) = dens_ts_saver_01(7,end);
init(7) = dens_ts_saver_01(8,end);
init(8) = dens_ts_saver_01(9,end);
init(9) = dens_ts_saver_01(10,end);
init(10) = dens_ts_saver_01(11,end);
init(11) = dens_ts_saver_01(12,end);
init(12) = dens_ts_saver_01(13,end);
save('init.mat','init')

clear all

% TNADWvec is vector of NADW values to loop through in Sv
TNADWvec = [4 5 6 7 7.25 7.5 7.75 8 8.25 8.5 8.75 9 9.25 9.5 10 11 12 13 14 15 16 17 18 19 20 19 18 17 16 15 14 13 12 11 10 9.5 9.25 9 8.75 8.5 8.25 8 7.75 7.5 7.25 7 6 5 4];
outmat = zeros(15,101,length(TNADWvec));

for ii = 1:length(TNADWvec)
    TNADWin = TNADWvec(ii).*10^6;   % Convert Sv to m^3 s^-1
    AABWratio_in = 1;
    phi_in = 1;
    y_si_in = -1e6;    % Sea ice edge position = -ell/2
    
    % Initialize with high temporal resolution run for 10 years at each step
    [dens_ts_saver] = TwoBasinFourLayer_td_lp_2019Published(10,.1,.1,TNADWin,phi_in,y_si_in,AABWratio_in);
    init(1) = dens_ts_saver(2,end);
    init(2) = dens_ts_saver(3,end);
    init(3) = dens_ts_saver(4,end);
    init(4) = dens_ts_saver(5,end);
    init(5) = dens_ts_saver(6,end);
    init(6) = dens_ts_saver(7,end);
    init(7) = dens_ts_saver(8,end);
    init(8) = dens_ts_saver(9,end);
    init(9) = dens_ts_saver(10,end);
    init(10) = dens_ts_saver(11,end);
    init(11) = dens_ts_saver(12,end);
    init(12) = dens_ts_saver(13,end);
    save('init.mat','init')
    clear dens_ts_saver init

    % Full time run
    eval(['dens_ts_saver_' num2str(ii) ' = TwoBasinFourLayer_td_lp_2019Published(10000,.5,50,TNADWin,phi_in,y_si_in,AABWratio_in);'])
   
    ii
end

