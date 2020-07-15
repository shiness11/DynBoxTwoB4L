# Code and example initial conditions file for “The role of the Southern Ocean in abrupt transitions and hysteresis in glacial ocean circulation” by Sophia K.V. Hines, Andrew F. Thompson, and Jess F. Adkins, Paleoceanography and Paleoclimatology, 2019.

## The following scripts can be used to produce Figures 4, 5, 6, 7, 8, 9, 11

**List of contents:**
 - TwoBasinFourLayer_td_lp_2019Published.m
 - NADWloop_2019Published.m
 - Example_Figure4a.m
 - dens_ts_saver_PLOTZZZ_2019Published.m
 - Example_phi1_NADW4.mat

**TwoBasinFourLayer_td_lp_2019Published.m :**
 
This script solves the time-dependent two-basin four-layer model. The input parameters are: “tfin_yr” the length of the simulation in years, “dt_day” the time step in days, “save_freq” the saving frequency in years, “TNADW_in” the NADW flux in m^3 s^-1, “phi_in” the NADW density, “y_si_in” the sea ice extent, and “AABWratio_in” the relative strength of AABW formation in the Atlantic versus the Pacific (AABWratio). To make the hysteresis loop, the model is run iteratively using the script “NADWloop_2019Published.m”.


**NADWloop_2019Published.m :**

This script iteratively runs the script “TwoBasinFourLayer_td_lp_2019Published.m” at different NADW flux values. The output of each run is saved and used as initial conditions for the subsequent run in the loop to create the hysteresis loop. An example initial condition is provided (“Example_phi1_NADW4.mat”) for phi = 1, NADW = 4. ***NOTE: This script will take a long time to run (~24 hours) ***


**Example_Figure4a.m :**

Run this script after running “NADWloop_2019Published.m” from the same directory to reproduce Figure 4a.


**dens_ts_saver_PLOTZZZ_2019Published.m :**

This script plots each time-dependent run for the hysteresis loop (so that the run can be evaluated to determine whether or not it has come to steady state). It also calculates T_res, T_diff, and T_zonal for all runs.


**Example_phi1_NADW4.mat :**

Initial conditions example—see NADWloop_2019Published.m
