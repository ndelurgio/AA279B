%% Script settings
clear all; close all; clc;

% Latex formatting
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% CONFIGURATION
ti_utc = datetime(2028, 1, 1, 1, 1, 1); % time at start of Earth re-entry


%% Constants
mu_Sun = 1.3271244004193938*10^11; % [km^3/s^2]

