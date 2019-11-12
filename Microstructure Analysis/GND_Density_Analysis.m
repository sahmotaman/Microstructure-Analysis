%%=========================================================================
% MTEX (5.1.1) script for dislocation density density analysis
%
%
% Author: S. Amir H. Motaman
%         Steel Institute (IEHK), RWTH Aachen University
%
% Reference:
%    Motaman, S.A.H.; Roters, F.; Haase, C.; 2019. Acta Materialia.
%    Anisotropic polycrystal plasticity due to microstructural heterogeneity:
%    A multi-scale experimental and numerical study on additively
%    manufactured metallic materials.
%%=========================================================================

%% Input parameters

% Starting the parallel pool
% parpool(4)

% Specifying the EBSD object
ebsd = ebsd_12_oBD;

% Specifying the phase
phase = 'fcc';

% Minimum high-angle grain boundary misorientation (in degrees)
min_GB_misorientation = 10;

% Working directory
working_directory = pwd;

% Output file name
output_file_name = 'LPBF_X30MnAl19-1_12.5%Strained_OrthogonaltoBD';

%% Smoothing the EBSD data and resolving the grains

[grains,ebsd.grainId] = calcGrains(ebsd,'angle',min_GB_misorientation*degree)
ebsd(grains(grains.grainSize <= 5)) = [];
F = halfQuadraticFilter;
ebsd_smoothed = smooth(ebsd(phase),F,'fill')
grains = calcGrains(ebsd_smoothed(phase),'angle',min_GB_misorientation*degree);
grains = smooth(grains,4);

%% Calculation of curvature tensor and dislocation density

ebsd_gridified = ebsd_smoothed(phase).gridify;
kappa = ebsd_gridified.curvature;
ds = dislocationSystem.fcc(ebsd_smoothed(phase).CS);
ds_rot = ebsd_gridified.orientations * ds;
[rho,factor] = fitDislocationSystems(kappa,ds_rot);
rho_sum = factor * sum(abs(rho .* ds_rot.u),2);
rho_mean = nanmean(rho_sum)
log_rho_sum = log10(rho_sum);
rho_log_mean = 10^(nanmean(log_rho_sum))

%% Plotting and storing the dislocation density map

% Plotting the GND density map
figure
plot(ebsd_smoothed,rho_sum)
mtexColorMap('parula')
mtexColorbar
% set(gca,'ColorScale','log');
set(gca,'CLim',[0 1e15]);
hold on
plot(grains.boundary,'linewidth',1.5,'linecolor','k')
hold off

% Saving the plotted figure
% saveas(gcf,[output_file_name,'_GND.fig'])
% saveas(gcf,[output_file_name,'_GND.jpg'])

%% Storing the workspace variables

% clearvars ds_rot rho
% save([output_file_name,'_GND.mat'])
