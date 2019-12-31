%%=========================================================================
% MTEX (5.1.1) script for assignment of crystallographic texture to the
% grains in a morphological volume element
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

% SPPARKS file name
file_name = 'Experimental_RVE';

% The measured/experimental ODF to be used for sampling the orientations
odf_experimental = odf_xrd_asbuilt_oBD;

% Applying the necessary rotation on the ODF
% rot = rotation('axis',xvector,'angle',90*degree);
% odf_experimental = rotate(odf_experimental,rot);

%% Loading the morphological volume element

% Loading the geometry of the morphological volume element
sites = dlmread([file_name '.spparks'],'',10);

% Grid resolution of the morphological volume element
n_grid = round((size(sites,1))^(1/3));

% Writing the geometry matrix
geometry_matrix = zeros((n_grid * n_grid),n_grid);
for i = 1:n_grid^2
    geometry_matrix(i,:) = sites(((i-1)  * n_grid + 1):(i * n_grid),2);
end

% Finding the total number of grains in the volume element
n_grains = max(sites(:,2));
counts = histc(sites(:,2),1:n_grains);

% Calculating the volume fraction of each grain in the volume element
f_grains = counts / n_grid^3;

%% Assigning file names

damask_geometry_file_name = [working_directory '\geometry_' file_name '_NG' int2str(n_grains) '_GR' int2str(n_grid) '.geom'];
orientations_file_name = [working_directory '\sampled_orientations_' file_name '_NG' int2str(n_grains) '.txt'];
microstructure_file_name = [working_directory '\microstructure_NG' int2str(n_grains) '.config'];
texture_file_name   = [working_directory '\texture_' file_name '_NG' int2str(n_grains) '.config'];
workspace_file_name = [working_directory '\workspace_' file_name '_NG' int2str(n_grains) '.mat'];
IPFs_fig_file_name = [working_directory '\IPFs_' file_name '_NG' int2str(n_grains) '.fig'];
IPFs_automatic_colorbar_jpg_file_name = [working_directory '\IPFs_automatic_colorbar_' file_name '_NG' int2str(n_grains) '.jpg'];
IPFs_fixed_colorbar_jpg_file_name = [working_directory '\IPFs_fixed_colorbar_' file_name '_NG' int2str(n_grains) '.jpg'];
IPFs_individual_orientations_jpg_file_name = [working_directory '\IPFs_individual_orientations_' file_name '_NG' int2str(n_grains) '.jpg'];

%% Creating DAMASK's geometry file

% Wrting the headers
damask_geom_file = fopen(damask_geometry_file_name,'w');
fprintf(damask_geom_file,'5\theader\n');
fprintf(damask_geom_file,'grid\ta %u\tb %u\tc %u\n',[n_grid n_grid n_grid]);
fprintf(damask_geom_file,'size\tx 1.000000\ty 1.000000\tz 1.000000\n');
fprintf(damask_geom_file,'origin\tx 0.000000\ty 0.000000\tz 0.000000\n');
fprintf(damask_geom_file,'microstructures %u\n',n_grains);
fprintf(damask_geom_file,'homogenization\t1\n');
fclose(damask_geom_file);

% Writing the sites
dlmwrite(damask_geometry_file_name,geometry_matrix,'-append','delimiter','\t');

%% Sampling discrete orientations

% Sampling a discrete set of orientations from the available ODF
orientations = discreteSample(odf_experimental,n_grains,'withoutReplacement','resolution',3*degree);

% Assembling the matrix of Bunge-Euler angles
euler_angles = zeros(n_grains,3);
euler_angles(:,1) = orientations.phi1 / degree;
euler_angles(:,2) = orientations.Phi  / degree;
euler_angles(:,3) = orientations.phi2 / degree;

%% Creating DAMASK's microstrcture file

microstructure_file = fopen(microstructure_file_name,'w');
for i = 1 : n_grains
    fprintf(microstructure_file,'[Grain%u]\n', i);
    fprintf(microstructure_file,'crystallite 1\n');
    fprintf(microstructure_file,'(constituent)   phase 1   texture %2u   fraction 1.0\n', i);
end
fclose(microstructure_file);

%% Creating DAMASK's texture file

texture_file = fopen(texture_file_name,'w');
for i = 1 : n_grains
    fprintf(texture_file,'[Grain%u]\n', i);
    fprintf(texture_file,'(gauss)   phi1 %3.1f   Phi %3.1f   phi2 %3.1f   scatter 0.0   fraction 1.0\n', euler_angles(i,:));
end
fclose(texture_file);

% dlmwrite(orientations_file_name,euler_angles,' ');

%% Creating the ODF, plotting the IPFs and storing everything

odf_bulk = calcODF(orientations,f_grains,'noFourier');
% odf_bulk = calcODF(orientations);
% odf_bulk.components{1,1}.weights = f_grains;
ODF_deviation_error = calcError(odf_bulk,odf_experimental)

figure
plotIPDF(odf_bulk,[xvector,yvector,zvector],'antipodal','contours',20,'LineStyle','-','fill','on','linecolor','w')
mtexColorbar
CLim(gcm,'equal')
% CLim(gcm,[0.85,1.45])
mtexColorbar
mtexColorbar
% saveas(gcf,IPFs_fixed_colorbar_jpg_file_name)
% saveas(gcf,IPFs_fixed_colorbar_jpg_file_name)

% figure
% plotIPDF(odf_bulk,[xvector,yvector,zvector],'antipodal','contours',20,'LineStyle','-','fill','on','linecolor','w')
% mtexColorbar
% saveas(gcf,IPFs_automatic_colorbar_jpg_file_name)
% saveas(gcf,IPFs_automatic_colorbar_jpg_file_name)

% figure
% plotIPDF(orientations,[xvector,yvector,zvector])
% saveas(gcf,IPFs_individual_orientations_jpg_file_name)
% saveas(gcf,IPFs_individual_orientations_jpg_file_name)

% save(workspace_file_name)
