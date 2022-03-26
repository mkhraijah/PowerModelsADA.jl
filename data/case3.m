function mpc = case3
%CASE3BUS_P6_6  Case of 3 bus system.
%   From Problem 6.6 in book 'Computational
%   Methods for Electric Power Systems' by Mariesa Crow
%   created by Rui Bo on 2007/11/12

%   MATPOWER
 %% MATPOWER Case Format : Version 2
mpc.version = '2';

 %%-----  Power Flow Data  -----%%
 %% system MVA base
mpc.baseMVA = 100;

%% bus data
%    bus_i    type    Pd    Qd    Gs    Bs    area    Vm    Va    baseKV    zone    Vmax    Vmin
mpc.bus = [
    1   3   350     1.00    0    0    1    1.05    0    230    1    1.10    0.90;
    2   2   400     2.50    0    0    2    1    0    230    1    1.10    0.90;
    3   2   250     1.00    0    0    3    1    0    230    1    1.10    0.90;
];

%% generator data
% Note:
% 1)It's better of gen to be in number order, otherwise gen and genbid
% should be sorted to make the lp solution output clearly(in number order as well)
% 2)set Pmax to nonzero. set to 999 if no limit
% 3)If change the order of gen, then must change the order in genbid
% accordingly
%    bus    Pg    Qg    Qmax    Qmin    Vg    mBase    status    Pmax Pmin Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
   1    1.8218   0    6.00    -600    1.05       100    1    600    0	0	0	0	0	0	0	0	0	0	0	0;
   2    2.7277   0    4.00    -400    1.02       100    1    400    0	0	0	0	0	0	0	0	0	0	0	0;
   3    5.4505   0    1.00    -100    1.02       100    1    100    0	0	0	0	0	0	0	0	0	0	0	0;
];
%gen(:, 9) = 999; % inactive the Pmax constraints

%% branch data
%    fbus    tbus    r    x    b    rateA    rateB    rateC    ratio    angle    status angmin	angmax
mpc.branch = [
    1    2    0.01     0.1    0.050    999    100    100    0    0    1     -360    360;
    1    3    0.05     0.1    0.025    999    100    100    0    0    1     -360    360;
    2    3    0.05     0.1    0.025    999    100    100    0    0    1     -360    360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%    2    startup    shutdown    n    c(n-1)    ...    c0
mpc.gencost = [
   2    0    0    3    1.5     1    0;
   2    0    0    3    1       2    0;
   2    0    0    3    0.5     2.5  0;
];
