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
mu_Earth = 398600.4418; % [km^3/s^2]
mu_Mars = 4.282837e4; 

% Planetary radii
% Source: https://science.nasa.gov/resource/solar-system-sizes/
R_Earth = 6371; % [km]
R_Mars = 3390; 

% Visualization
earth_blue = [11,98,212]/255;
mars_red = [252,80,3]/255;