%%=========================================================================
% MTEX (5.1.1) script for twin fraction analysis
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

% Specifying the EBSD object
ebsd = ebsd_12_oBD;

% Specifying the phase
phase = 'fcc';

%% Smoothing the EBSD data

% Selecting the filter type
%F = medianFilter;
%F = splineFilter;
F = halfQuadraticFilter;

% Defining the size of the window to be used for finding the median
%F.numNeighbours = 3;

% smooth the data
ebsd = smooth(ebsd_original(phase),F,'fill');

% plot the smoothed data
ipfKey = ipfHSVKey(ebsd(phase));
ipfKey.inversePoleFigureDirection = zvector;
color = ipfKey.orientation2color(ebsd(phase).orientations);
figure;
plot(ebsd(phase),color,'micronbar','off')

%% Twin Boundary analysis

% Segmenting the grains
[grains,ebsd(phase).grainId,ebsd(phase).mis2mean] = calcGrains(ebsd(phase),'angle',10*degree);

% Removing the two pixel grains
ebsd(grains(grains.grainSize<=2)) = [];

% smoothing the grain boundaries
grains = grains.smooth(5);

% visualize the grains
figure;
plot(grains,grains.meanOrientation)

% Storing the crystal symmetry
cs = grains.CS;

% Restricting the grain boundaries to a specific phase transistion
grain_boundary = grains.boundary(phase,phase);

% Detecting the twin (CSL-3) boundaries based on an orientation relationsship
twin_boundary_definition = CSL(3,ebsd(phase).CS);
% twinning = orientation('map',Miller(1,1,-1,cs),Miller(-5,-1,1,cs),Miller(1,-1,-1,cs),Miller(-1,1,5,cs));

% Printing the rotational angle-axis of the corresponding twin boundaries
round(twin_boundary_definition.axis)
twin_boundary_definition.angle / degree

% Identifying the twin boundaries with an angular tolerance of 5 degrees
is_twin_boundary = angle(grain_boundary.misorientation,twin_boundary_definition) < 5 * degree;
twin_boundary = grain_boundary(is_twin_boundary);

% Plotting the twinning boundaries
ipfKey = ipfColorKey(grains);
ipfKey.inversePoleFigureDirection = zvector;
color = ipfKey.orientation2color(grains.meanOrientation);
figure;
plot(grains,color,'micronbar','off')
hold on
%plot(gB,angle(gB.misorientation,twinning),'linewidth',4)
plot(twin_boundary,'linecolor','b','linewidth',1,'displayName','twin boundary')
hold off

% Merging the twins along twin boundaries
% [mergedGrains,parentId] = merge(grains,twin_boundary,'calcMeanOrientation');
% hold on
% plot(mergedGrains.boundary,'linecolor','k','linewidth',2.0,'linestyle','-','displayName','merged grains','micronbar','off')
% hold off

%% Automatic selection of twinned regiones from matrix

unique_twin_boundary = unique(twin_boundary.grainId,'rows');

n = size(unique_twin_boundary,1);
twin_ID = zeros(n,1);
for i = 1:n
    if grains(unique_twin_boundary(i,1)).grainSize < grains(unique_twin_boundary(i,2)).grainSize
        twin_ID(i,1) = unique_twin_boundary(i,1);
    else
        twin_ID(i,1) = unique_twin_boundary(i,2);
    end
end
twin_ID = unique(twin_ID,'rows');
new_twin_id = [];

%% Manual removal of parent grains

twin_ID = [twin_ID; new_twin_id'];
twin_ID = unique(twin_ID,'rows');

ipfKey = ipfColorKey(grains(twin_ID));
ipfKey.inversePoleFigureDirection = zvector;
color = ipfKey.orientation2color(grains(twin_ID).meanOrientation);
figure;
plot(grains(twin_ID),color,'micronbar','off')
hold on
plot(twin_boundary,'linecolor','b','linewidth',1)
hold off

points = [];
new_matrix_id = [];
for i=1:100
	temp = impoint;
	points(i,:) = temp.getPosition;
	new_matrix_id(i) = grains(points(i,1),points(i,2)).id;
end

%% Manual addition of of twins

for i = 1:size(new_matrix_id,2)
    twin_ID = twin_ID(twin_ID ~= new_matrix_id(i));
end

grains_ID = grains.id;

for i = 1:size(twin_ID,1)
    grains_ID = grains_ID(grains_ID ~= twin_ID(i));
end

ipfKey = ipfColorKey(grains(grains_ID));
ipfKey.inversePoleFigureDirection = zvector;
color = ipfKey.orientation2color(grains(grains_ID).meanOrientation);
figure;
plot(grains(grains_ID),color,'micronbar','off')
hold on
plot(twin_boundary,'linecolor','b','linewidth',1)
hold off

points = [];
new_twin_id = [];
for i=1:100
	temp = impoint;
	points(i,:) = temp.getPosition;
	new_twin_id(i) = grains(points(i,1),points(i,2)).id;
end

%% Calculation of twin area fraction

twin_ID = [twin_ID; new_twin_id'];
twin_ID = unique(twin_ID,'rows');

ipfKey = ipfColorKey(grains(twin_ID));
ipfKey.inversePoleFigureDirection = zvector;
color = ipfKey.orientation2color(grains(twin_ID).meanOrientation);
figure;
plot(grains(twin_ID),color,'micronbar','off')
hold on
plot(twin_boundary,'linecolor','b','linewidth',1,'displayName','twin boundary')
hold off

twins = grains(twin_ID);
twin_area_fraction = sum(area(twins)) / sum(area(grains))

