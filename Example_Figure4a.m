% Run this script after running "NADWloop_2019Published.m" (in the same directory) to reproduce Figure 4a from the paper.
% 
% SH 09/12/19

NADWu_phi1 = TNADWvec(1:25);
NADWd_phi1 = TNADWvec(25:end);

up_phi1(1) = dens_ts_saver_1(3,end)./2e6;
up_phi1(2) = dens_ts_saver_2(3,end)./2e6;
up_phi1(3) = dens_ts_saver_3(3,end)./2e6;
up_phi1(4) = dens_ts_saver_4(3,end)./2e6;
up_phi1(5) = dens_ts_saver_5(3,end)./2e6;
up_phi1(6) = dens_ts_saver_6(3,end)./2e6;
up_phi1(7) = dens_ts_saver_7(3,end)./2e6;
up_phi1(8) = dens_ts_saver_8(3,end)./2e6;
up_phi1(9) = dens_ts_saver_9(3,end)./2e6;
up_phi1(10) = dens_ts_saver_10(3,end)./2e6;
up_phi1(11) = dens_ts_saver_11(3,end)./2e6;
up_phi1(12) = dens_ts_saver_12(3,end)./2e6;
up_phi1(13) = dens_ts_saver_13(3,end)./2e6;
up_phi1(14) = dens_ts_saver_14(3,end)./2e6;
up_phi1(15) = dens_ts_saver_15(3,end)./2e6;
up_phi1(16) = dens_ts_saver_16(3,end)./2e6;
up_phi1(17) = dens_ts_saver_17(3,end)./2e6;
up_phi1(18) = dens_ts_saver_18(3,end)./2e6;
up_phi1(19) = dens_ts_saver_19(3,end)./2e6;
up_phi1(20) = dens_ts_saver_20(3,end)./2e6;
up_phi1(21) = dens_ts_saver_21(3,end)./2e6;
up_phi1(22) = dens_ts_saver_22(3,end)./2e6;
up_phi1(23) = dens_ts_saver_23(3,end)./2e6;
up_phi1(24) = dens_ts_saver_24(3,end)./2e6;
up_phi1(25) = dens_ts_saver_25(3,end)./2e6;
down_phi1(1) = dens_ts_saver_25(3,end)./2e6;
down_phi1(2) = dens_ts_saver_26(3,end)./2e6;
down_phi1(3) = dens_ts_saver_27(3,end)./2e6;
down_phi1(4) = dens_ts_saver_28(3,end)./2e6;
down_phi1(5) = dens_ts_saver_29(3,end)./2e6;
down_phi1(6) = dens_ts_saver_30(3,end)./2e6;
down_phi1(7) = dens_ts_saver_31(3,end)./2e6;
down_phi1(8) = dens_ts_saver_32(3,end)./2e6;
down_phi1(9) = dens_ts_saver_33(3,end)./2e6;
down_phi1(10) = dens_ts_saver_34(3,end)./2e6;
down_phi1(11) = dens_ts_saver_35(3,end)./2e6;
down_phi1(12) = dens_ts_saver_36(3,end)./2e6;
down_phi1(13) = dens_ts_saver_37(3,end)./2e6;
down_phi1(14) = dens_ts_saver_38(3,end)./2e6;
down_phi1(15) = dens_ts_saver_39(3,end)./2e6;
down_phi1(16) = dens_ts_saver_40(3,end)./2e6;
down_phi1(17) = dens_ts_saver_41(3,end)./2e6;
down_phi1(18) = dens_ts_saver_42(3,end)./2e6;
down_phi1(19) = dens_ts_saver_43(3,end)./2e6;
down_phi1(20) = dens_ts_saver_44(3,end)./2e6;
down_phi1(21) = dens_ts_saver_45(3,end)./2e6;
down_phi1(22) = dens_ts_saver_46(3,end)./2e6;
down_phi1(23) = dens_ts_saver_47(3,end)./2e6;
down_phi1(24) = dens_ts_saver_48(3,end)./2e6;
down_phi1(25) = dens_ts_saver_49(3,end)./2e6;

figure(1);clf
plot(NADWu_phi1,up_phi1,'bo--')
hold on
plot(NADWd_phi1,down_phi1,'bo-')
xlabel('NADW flux (Sv)')
ylabel('y2A outcrop position')
title('Figure 4a')