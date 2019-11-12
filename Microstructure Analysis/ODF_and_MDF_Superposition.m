%%=========================================================================
% MTEX (5.1.1) script for sampling and superposition of ODFs and correlated
% MDFs measured on two orthogonal planes
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

% Working directory
working_directory = pwd;

% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

% Matrix/main phase
phase = 'fcc';

% Minimum grain size (euivalent diameter) to be considered (in micron)
min_grain_size = 2.8;

% Number of bins in grain analysis
n_bins = 30;

% Minimum high-angle grain boundary misorientation (in degrees)
min_GB_misorientation = 10;

% Assigning weights to the orthogonal planes
weight_oBD = 0.5;
% weight_pBD = 1 - weight_oBD;

% Number of total orientations in the superposed discretized ODF file
n_orientations = 5000;

% Number of total items in the superposed discretized MDF file
n_misorientations = 200;

% Name of the superposed discretized ODF file
odf_file_name = 'ODF.txt';

% Name of the superposed discretized MDF file
mdf_file_name = 'MDF.txt';

%% Analysis of orhthogonal-to-BD (oBD) EBSD data

% Grain anaylsis
grains_oBD = calcGrains(ebsd_asbuilt_oBD(phase),'angle',min_GB_misorientation*degree);
grains_oBD = smooth(grains_oBD,4);

% ODF anaylsis
odf_ebsd_oBD = calcODF(ebsd_asbuilt_oBD(phase).orientations);

% Rotation of the ODF if necessary
% rot = rotation('axis',xvector,'angle',90*degree);
% odf_ebsd_oBD = rotate(odf_ebsd_oBD,rot);

%% Analysis of parallel-to-BD (pBD) EBSD data

% Grain anaylsis
grains_pBD = calcGrains(ebsd_asbuilt_pBD(phase),'angle',min_GB_misorientation*degree);
grains_pBD = smooth(grains_pBD,4);

% ODF anaylsis
odf_ebsd_pBD = calcODF(ebsd_asbuilt_pBD(phase).orientations);

% Rotation of the ODF if necessary
rot = rotation('axis',xvector,'angle',90*degree);
odf_ebsd_pBD = rotate(odf_ebsd_pBD,rot);

%% ODF analysis of orhthogonal-to-BD (oBD) XRD data

% Rotation of data about an axis
rot = rotation('axis',xvector,'angle',90*degree);
pf_asbuilt_oBD = rotate(pf_asbuilt_oBD,rot);

% ODF anaylsis
odf_xrd_oBD = calcODF(pf_asbuilt_oBD);

%% ODF analysis of parallel-to-BD (pBD) XRD data

% Rotation of data about an axis
% rot = rotation('axis',xvector,'angle',90*degree);
% pf_asbuilt_pBD = rotate(pf_asbuilt_pBD,rot);

% ODF anaylsis
odf_xrd_pBD = calcODF(pf_asbuilt_pBD);

%% Generating a superposed discretized ODF

% Calculation of the number of orientations to be sampled from each ODF
n_orientations_oBD = round(weight_oBD * n_orientations);
n_orientations_pBD = n_orientations - n_orientations_oBD;

% Setting the ODF of orthogonal planes
odf_oBD = odf_xrd_oBD;
odf_pBD = odf_xrd_pBD;

% Sampling orientations from the orhthogonal-to-BD (oBD) plane
sampled_orientations_oBD = discreteSample(odf_oBD,n_orientations_oBD,'withoutReplacement');
euler_angles_oBD = zeros(n_orientations_oBD,3);
euler_angles_oBD(:,1) = sampled_orientations_oBD.phi1 / degree;
euler_angles_oBD(:,2) = sampled_orientations_oBD.Phi  / degree;
euler_angles_oBD(:,3) = sampled_orientations_oBD.phi2 / degree;

% Sampling orientations from the parallel-to-BD (pBD) plane
sampled_orientations_pBD = discreteSample(odf_pBD,n_orientations_pBD,'withoutReplacement');
euler_angles_pBD = zeros(n_orientations_pBD,3);
euler_angles_pBD(:,1) = sampled_orientations_pBD.phi1 / degree;
euler_angles_pBD(:,2) = sampled_orientations_pBD.Phi  / degree;
euler_angles_pBD(:,3) = sampled_orientations_pBD.phi2 / degree;

% Superposition of the sampled orientations associated with the orthogonal planes
euler_angles = [euler_angles_oBD;euler_angles_pBD];
euler_angles(:,4) = 1;
euler_angles(:,5) = 1;

% Writing the superposed discretized ODF file
dlmwrite([working_directory '\' odf_file_name],['Angle Count:' num2str(n_orientations)],'delimiter','');
dlmwrite([working_directory '\' odf_file_name],euler_angles,'-append','delimiter',' ');

%% Reading and plotting the superposed discretized ODF file

% Crystal symmetry
cs_superposed = crystalSymmetry('m-3m', [1 1 1]);

% Specimen symmetry
ss_superposed = specimenSymmetry('1');

% Specifying kernel
psi = deLaValeePoussinKernel('halfwidth',10*degree);

% Reading the superposed discretized ODF file
odf_superposed = loadODF([working_directory '\' odf_file_name],cs_superposed,ss_superposed,'density','kernel',psi,'resolution',5*degree,...
                 'interface','generic','ColumnNames',{'phi1' 'Phi' 'phi2' 'weight'},'Bunge','Degree');

% Plotting the superposed discretized ODF file
figure
plotIPDF(odf_superposed,[xvector,yvector,zvector],'antipodal','contours',20,'LineStyle','-','fill','on','linecolor','w')
mtexColorbar
CLim(gcm,'equal')
% CLim(gcm,[0.0,4.2])
mtexColorbar
mtexColorbar

%% Generating a superposed discretized MDF

% Calculation of the number of misorientations to be sampled from each MDF
n_misorientations_oBD = round(weight_oBD * n_misorientations);
n_misorientations_pBD = n_misorientations - n_misorientations_oBD;

% Sampling misorientations from the orhthogonal-to-BD (oBD) plane
condition = grains_oBD.equivalentRadius > (min_grain_size / 2.0);
correlated_mdf_oBD = grains_oBD(condition).boundary(phase,phase).misorientation;
mis_axis_angle_oBD = [correlated_mdf_oBD.angle,correlated_mdf_oBD.axis.hkl];
sampled_mis_axis_angle_oBD = datasample(mis_axis_angle_oBD,n_misorientations_oBD);

% Sampling misorientations from the parallel-to-BD (pBD) plane
condition = grains_pBD.equivalentRadius > (min_grain_size / 2.0);
correlated_mdf_pBD = grains_pBD(condition).boundary(phase,phase).misorientation;
mis_axis_angle_pBD = [correlated_mdf_pBD.angle,correlated_mdf_pBD.axis.hkl];
sampled_mis_axis_angle_pBD = datasample(mis_axis_angle_pBD,n_misorientations_pBD);

% Superposition of the misorientations associated with the orthogonal planes
superposed_mis_axis_angle = [mis_axis_angle_oBD;mis_axis_angle_pBD];
superposed_sampled_mis_axis_angle = [sampled_mis_axis_angle_oBD;sampled_mis_axis_angle_pBD];
superposed_sampled_mis_axis_angle(:,5) = 1;

% Writing the superposed discretized MDF file
dlmwrite([working_directory '\' mdf_file_name],n_misorientations);
dlmwrite([working_directory '\' mdf_file_name],superposed_sampled_mis_axis_angle,'-append','delimiter',' ');

% Plotting the misorientation angle distribution associated with the
% superposed and sampled-superposed discretized MDFs
figure
histogram(superposed_sampled_mis_axis_angle(1:end,1) / degree,n_bins,'Normalization','probability');
[counts,centers] = hist(superposed_mis_axis_angle(1:end,1) / degree,n_bins);
hold on
plot (centers,counts/sum(counts))
hold off
