%%=========================================================================
% MTEX (5.1.1) script for grain size, morphology and texture anaylsis
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

% Matrix/main phase
phase = 'fcc';

% Minimum grain size (euivalent diameter) to be considered (in micron)
min_grain_size = 2.8;

% Number of bins in grain analysis
n_bins = 30;

% Minimum high-angle grain boundary misorientation (in degrees)
min_GB_misorientation = 10;

% Specifying the EBSD data object
ebsd = ebsd_asbuilt_oBD;

% Specifying the XRD pole figure object
pf = pf_asbuilt_oBD;

%% Evaluating the ODF from the EBSD data

% psi = calcKernel(grains(phase).meanOrientation);
% odf = calcODF(ebsd(phase).orientations,'kernel',psi);
odf_ebsd = calcODF(ebsd(phase).orientations);

%% Evaluating the ODF from the XRD data

% Rotation of data about an axis
% rot = rotation('axis',xvector,'angle',90*degree);
% pf = rotate(pf,rot);

% Calculation of ODF
odf_xrd = calcODF(pf);

%% Plotting the IPFs

odf = odf_xrd;

figure
plotIPDF(odf,[xvector,yvector,zvector],'antipodal','contours',20,'LineStyle','-','fill','on','linecolor','w')
mtexColorbar
CLim(gcm,'equal')
%CLim(gcm,[0,4])
mtexColorbar
mtexColorbar

% annotate([Miller(1,-1,-1,odf.CS)],'label',{'111'},'BackgroundColor','w','Tag','setAboveMarker','FontSize',30)
% annotate([Miller(0, 0, 1,odf.CS)],'label',{'001'},'BackgroundColor','w','Tag','setBelowMarker','FontSize',30)
% annotate([Miller(0, 1, 1,odf.CS)],'label',{'011'},'BackgroundColor','w','Tag','setBelowMarker','FontSize',30)

%% Smoothing EBSD data and grains

[grains,ebsd.grainId] = calcGrains(ebsd,'angle',min_GB_misorientation*degree);
ebsd(grains(grains.grainSize < 10)) = [];

F = halfQuadraticFilter;
ebsd = smooth(ebsd,F,'fill');

%% Grain analysis

grains = calcGrains(ebsd(phase),'angle',min_GB_misorientation*degree);
grains = smooth(grains,4);

figure
plot(grains)

%% Grains filteration (optional removing the incomplete grains)

outerBoundary_id = any(grains.boundary.grainId == 0,2);
grain_id = grains.boundary(outerBoundary_id).grainId;
grain_id(grain_id == 0) = [];
grains(grains.id == grain_id) = [];
grains(grains.equivalentRadius <= (min_grain_size / 2.0)) = [];

figure
plot(grains)

%% Plotting EBSD IPF orientation map

ipfKey = ipfHSVKey(ebsd(phase));
% figure
% plot(ipfKey)

ipfKey.inversePoleFigureDirection = xvector;
% ipfKey.inversePoleFigureDirection = yvector;
% ipfKey.inversePoleFigureDirection = zvector;
color = ipfKey.orientation2color(ebsd(phase).orientations);

figure
plot(ebsd(phase),color)

hold on
plot(grains.boundary,'linewidth',1.0)
hold off

%% Grain size analysis

condition = grains.equivalentRadius > (min_grain_size / 2.0);
grain_area = grains(condition).area;
total_area = sum(grain_area);
grain_size = 2.0 * grains(condition).equivalentRadius;
n_grains = length(grain_size);

[bin_counts,bin_centers] = hist(grain_size,n_bins);
bin_counts = bin_counts';
bin_centers = bin_centers';
bin_width = bin_centers(3) - bin_centers(2);

bin_low = zeros(n_bins,1);
bin_high = zeros(n_bins,1);
for i = 1 : n_bins
    bin_low(i) = bin_centers(i) - bin_width / 2.0;
    bin_high(i) = bin_centers(i) + bin_width / 2.0;
end
max_grain_size = max(grain_size);
min_grain_size = min(grain_size);
bin_low(1) = min_grain_size;
bin_high(n_bins) = max_grain_size;

area_per_bin = zeros(n_bins,1);
area_counts_per_bin = zeros(n_bins,1);
for i = 1 : n_grains
      for j = 1 : n_bins
           if (grain_size(i) >= bin_low(j)) && (grain_size(i) <= bin_high(j))
              area_per_bin(j) = area_per_bin(j) + grain_area(i);
              area_counts_per_bin(j) = area_counts_per_bin(j) + 1;
           end
      end
end

bin_area_fraction = zeros(n_bins,1);
bin_area_fraction = 100.0 * area_per_bin / total_area;
bin_number_fraction = zeros(n_bins,1);
bin_number_fraction = 100.0 * bin_counts / n_grains;

mean_grain_size_area   = sum(bin_centers .* bin_area_fraction / 100.0)
mean_grain_size_number = sum(bin_centers .* bin_number_fraction / 100.0)
effective_grain_size = (mean_grain_size_area + mean_grain_size_number) / 2.0

figure
h = bar(bin_centers,[bin_area_fraction,bin_number_fraction],'histc','FaceColor','flat');
h(1,1).FaceColor = 'red';
h(1,2).FaceColor = 'blue';

l = cell(1,2);
l{1} = 'Area fraction';
l{2} = 'Number fraction';

legend(h,l);
set(gca,'FontSize',15)
% xticks(0:20:(bin_centers(n_bins)+ xticks_resolution));
% ylim([0 30])
% yticks(0:5:30);

xlabel('Grain size (\it{d}\rm) [\mum]')
ylabel('Fraction [%]')

% Fitting a log-normal function to the histogram
logn_fit = lognfit(grain_size);
mu = logn_fit(1)
sigma = logn_fit(2)
x = 0:0.1:round(max_grain_size);
hold on
p = plot(x,100.0 * bin_width * lognpdf(x,mu,sigma))
hold off
set(p,'DisplayName','Lognormal fit','LineStyle','--','LineWidth',2.5)

%% Grain aspect ratio analysis

condition = grains.equivalentRadius > (min_grain_size / 2.0);
grain_area = grains(condition).area;
total_area = sum(grain_area);
grain_size = 2.0 * grains(condition).equivalentRadius;
n_grains = length(grain_size);
aspect_ratio = 1.0 ./ grains(condition).aspectRatio;

[bin_counts,bin_centers] = hist(aspect_ratio,n_bins);
bin_counts = bin_counts';
bin_centers = bin_centers';
bin_width = bin_centers(3) - bin_centers(2);

bin_low  = zeros(n_bins,1);
bin_high = zeros(n_bins,1);
for i = 1 : n_bins
    bin_low(i) = bin_centers(i) - bin_width / 2.0;
    bin_high(i) = bin_centers(i) + bin_width / 2.0;
end
bin_low(1) = 0.0;
bin_high(n_bins) = 1.0;

area_per_bin = zeros(n_bins,1);
area_counts_per_bin = zeros(n_bins,1);
for i = 1 : n_grains
      for j = 1 : n_bins
           if (aspect_ratio(i) >= bin_low(j)) && (aspect_ratio(i) <= bin_high(j))
              area_per_bin(j) = area_per_bin(j) + grain_area(i);
              area_counts_per_bin(j) = area_counts_per_bin(j) + 1;
           end
      end
end

bin_area_fraction = zeros(n_bins,1);
bin_area_fraction = 100.0 * area_per_bin / total_area;
bin_number_fraction = zeros(n_bins,1);
bin_number_fraction = 100.0 * bin_counts / n_grains;

mean_aspect_ratio_area   = sum(bin_centers .* bin_area_fraction / 100.0)
mean_aspect_ratio_number = sum(bin_centers .* bin_number_fraction / 100.0)
effective_aspect_ratio = (mean_aspect_ratio_area + mean_aspect_ratio_number) / 2.0

figure
h = bar(bin_centers,[bin_area_fraction,bin_number_fraction],'histc','FaceColor','flat');
h(1,1).FaceColor = 'red';
h(1,2).FaceColor = 'blue';

l = cell(1,2);
l{1} = 'Area fraction';
l{2} = 'Number fraction';

legend(h,l);
set(gca,'FontSize',15)
xticks(0:0.1:1);
% ylim([0 35])
% yticks(0:5:35);

xlabel('Grain aspect ratio (\itm\rm) [-]')
ylabel('Fraction [%]')

% Fitting a normal function to the histogram
[mu,sigma] = normfit(aspect_ratio)
x = 0:0.01:1;
hold on
p = plot(x,100.0 * bin_width * normpdf(x,mu,sigma))
hold off
set(p,'DisplayName','Normal fit','LineStyle','--','LineWidth',2.5)

%% Grain shape (ellipse axis) orientation analysis

condition = grains.equivalentRadius > (min_grain_size / 2.0);
grain_area = grains(condition).area;
total_area = sum(grain_area);
grain_size = 2.0 * grains(condition).equivalentRadius;
n_grains = length(grain_size);

[omega,a,b] = fitEllipse(grains(condition));
omega = omega / degree;
theta = zeros(n_grains,1);
for i = 1:n_grains
    if (omega(i) >= 0.0) && (omega(i) <=90.0)
        theta(i) = 90.0 - omega(i);
    else
        theta(i) = 270.0 - omega(i);
    end
end


[bin_counts,bin_centers] = hist(theta,n_bins);
bin_counts = bin_counts';
bin_centers = bin_centers';
bin_width = bin_centers(3) - bin_centers(2);

bin_low  = zeros(n_bins,1);
bin_high = zeros(n_bins,1);
for i = 1 : n_bins
    bin_low(i) = bin_centers(i) - bin_width / 2.0;
    bin_high(i) = bin_centers(i) + bin_width / 2.0;
end
bin_low(1) = 0.0;
bin_high(n_bins) = 180.0;

area_per_bin = zeros(n_bins,1);
area_counts_per_bin = zeros(n_bins,1);
for i = 1 : n_grains
      for j = 1 : n_bins
           if (theta(i) >= bin_low(j)) && (theta(i) <= bin_high(j))
              area_per_bin(j) = area_per_bin(j) + grain_area(i);
              area_counts_per_bin(j) = area_counts_per_bin(j) + 1;
           end
      end
end

bin_area_fraction = zeros(n_bins,1);
bin_area_fraction = 100.0 * area_per_bin / total_area;
bin_number_fraction = zeros(n_bins,1);
bin_number_fraction = 100.0 * bin_counts / n_grains;

mean_angle_area   = sum(bin_centers .* bin_area_fraction / 100.0)
mean_angle_number = sum(bin_centers .* bin_number_fraction / 100.0)
effective_angle = (mean_angle_area + mean_angle_number) / 2.0

figure
h = bar(bin_centers,[bin_area_fraction,bin_number_fraction],'histc','FaceColor','flat');
h(1,1).FaceColor = 'red';
h(1,2).FaceColor = 'blue';

l = cell(1,2);
l{1} = 'Area fraction';
l{2} = 'Number fraction';

legend(h,l);
set(gca,'FontSize',15)
xticks(0:20:180);
% ylim([0 35])
% yticks(0:5:35);

xlabel('Grain major axis angle (\it\omega\rm) [ {\circ} ]')
ylabel('Fraction [%]')

% Fitting a normal function to the histogram
[mu,sigma] = normfit(theta)
x = 0:0.1:180;
hold on
p = plot(x,100.0 * bin_width * normpdf(x,mu,sigma))
hold off
set(p,'DisplayName','Normal fit','LineStyle','--','LineWidth',2.5)


%% Misorientation analysis

condition = grains.equivalentRadius > (min_grain_size / 2.0);
correlated_mdf = grains(condition).boundary(phase,phase).misorientation;
uncorrelated_mdf = calcMDF(odf_xrd,odf_xrd);
% uncorrelated_mdf = calcMisorientation(ebsd(phase));

figure
plotAngleDistribution(correlated_mdf,'DisplayName','Correlated (boundary misorientation derived from EBSD data)','resolution',2*degree)
hold on
plotAngleDistribution(uncorrelated_mdf,'DisplayName','Uncorrelated (derived from XRD macro-texture data)','Color','r','resolution',0.1*degree)
hold off
legend('-dynamicLegend','Location','northwest')

xlabel('Misorientation angle [ {\circ} ]')
ylabel('Frequency [%]')
set(gca,'FontSize',15)
% ylim([0 10])
% yticks(0:2:10)

figure
plotAxisDistribution(correlated_mdf,'smooth')
mtexColorbar;
% CLim(gcm,[0.0,1.8])
% annotate([Miller(-1,1,1,correlated_mdf.CS)],'label',{'111'},'BackgroundColor','w','Tag','setAboveMarker','FontSize',30)
% annotate([Miller(0, 0,1,correlated_mdf.CS)],'label',{'001'},'BackgroundColor','w','Tag','setBelowMarker','FontSize',30)
% annotate([Miller(0, 1,1,correlated_mdf.CS)],'label',{'011'},'BackgroundColor','w','Tag','setBelowMarker','FontSize',30)
mtexTitle('Correlated','FontSize',30)

figure
plotAxisDistribution(uncorrelated_mdf,'smooth')
mtexColorbar;
% CLim(gcm,[0.0,3.7])
% annotate([Miller(-1,1,1,uncorrelated_mdf.CS)],'label',{'111'},'BackgroundColor','w','Tag','setAboveMarker','FontSize',30)
% annotate([Miller(0, 0,1,uncorrelated_mdf.CS)],'label',{'001'},'BackgroundColor','w','Tag','setBelowMarker','FontSize',30)
% annotate([Miller(0, 1,1,uncorrelated_mdf.CS)],'label',{'011'},'BackgroundColor','w','Tag','setBelowMarker','FontSize',30)
mtexTitle('Uncorrelated','FontSize',30)

%% Taylor factor analysis

% Define a family of slip systems
sS = slipSystem.fcc(grains.CS);
sS = sS.symmetrise;

% Strain tensor
epsilon = strainTensor(diag([1 -0.5 -0.5]));    % uniaxial tension along x direction
%epsilon = strainTensor(diag([-0.5 1 -0.5]));    % uniaxial tension along y direction

% Apply Taylor model
condition = grains.equivalentRadius > (min_grain_size / 2.0);
[taylor_factor,b,uncorrelated_mdf] = calcTaylor(grains(condition).meanOrientation\epsilon,sS);

% Equivalent grain size
grain_size = 2.0 * grains(condition).equivalentRadius;

% Weighted average Taylor factor
mean_taylor_factor = sum((grain_size.^2).* taylor_factor) / sum(d.^2)
