clc;
clear;
close all;
% Read result
    KQ_aco=csvread('MSUMR_aco.txt');
    KQ_ssa=csvread('MSUMR_ssa.txt');
% Plot SSA
plot(KQ_aco, '-x','LineWidth',2);
hold on
% Plot ACO
plot(KQ_ssa, '-o','LineWidth',2);
 xlabel('CMax');
 ylabel('Max Sum Reward');
legend ('ACO','SSA', 'Location', 'southeast');
 grid on;
hold all;
hold off