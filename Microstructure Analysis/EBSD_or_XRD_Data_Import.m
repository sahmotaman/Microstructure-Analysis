%%=========================================================================
% MTEX (5.1.1) script for importing the EBSD and XRD data
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

% Working directory
working_directory = [pwd '\Microstructure Characterization Data\'];

% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Importing the as-built orhthogonal-to-BD (oBD) XRD pole figure data

% Crystal symmetry
cs = crystalSymmetry('m-3m', [1 1 1]);

% Specimen symmetry
ss = specimenSymmetry('1');

% Spcifying the files to be imported
file_names = {...
  [working_directory 'LPBF_X30MnAl19-1_AsBuilt_OrthogonaltoBD_01_002.gpol'],...
  [working_directory 'LPBF_X30MnAl19-1_AsBuilt_OrthogonaltoBD_01_022.gpol'],...
  [working_directory 'LPBF_X30MnAl19-1_AsBuilt_OrthogonaltoBD_01_111.gpol'],...
  };

% Specify Miller indice
h = { ...
  Miller(0,0,2,cs),...
  Miller(0,2,2,cs),...
  Miller(1,1,1,cs),...
  };

% Creating a pole figure variable containing the data
pf_asbuilt_oBD = loadPoleFigure(file_names,h,cs,ss,'interface','gpol');

% Rotation of data about an axis
rot = rotation('axis',xvector,'angle',90*degree);
pf_asbuilt_oBD = rotate(pf_asbuilt_oBD,rot);

% ODF anaylsis
odf_xrd_asbuilt_oBD = calcODF(pf_asbuilt_oBD);

%% Importing the as-built parallel-to-BD (pBD) XRD pole figure data

% Crystal symmetry
cs = crystalSymmetry('m-3m', [1 1 1]);

% Specimen symmetry
ss = specimenSymmetry('1');

% Spcifying the files to be imported
file_names = {...
  [working_directory 'LPBF_X30MnAl19-1_AsBuilt_ParalleltoBD_01_002.gpol'],...
  [working_directory 'LPBF_X30MnAl19-1_AsBuilt_ParalleltoBD_01_022.gpol'],...
  [working_directory 'LPBF_X30MnAl19-1_AsBuilt_ParalleltoBD_01_111.gpol'],...
  };

% Specify Miller indice
h = { ...
  Miller(0,0,2,cs),...
  Miller(0,2,2,cs),...
  Miller(1,1,1,cs),...
  };

% Creating a pole figure variable containing the data
pf_asbuilt_pBD = loadPoleFigure(file_names,h,cs,ss,'interface','gpol');

% Rotation of data about an axis
% rot = rotation('axis',xvector,'angle',90*degree);
% pf_asbuilt_pBD = rotate(pf_asbuilt_pBD,rot);

% ODF anaylsis
odf_xrd_asbuilt_pBD = calcODF(pf_asbuilt_pBD);

%% Importing the as-built orhthogonal-to-BD (oBD) EBSD data

% Crystal symmetry
cs = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [2.87 2.87 2.87], 'mineral', 'bcc', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [3.65 3.65 3.65], 'mineral', 'fcc', 'color', 'light green'),...
  crystalSymmetry('6/mmm', [2.53 2.53 4.079], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'hcp', 'color', 'light red')};

% Specifying the file to be imported
file_name = 'LPBF_X30MnAl19-1_AsBuilt_OrthogonaltoBD';
file_extension = '.ang';

% Creating an EBSD variable containing the data
ebsd_asbuilt_oBD = loadEBSD([working_directory,file_name,file_extension],cs,...
                   'interface',strip(file_extension,'left','.'),'convertEuler2SpatialReferenceFrame');

%% Importing the as-built parallel-to-BD (pBD) EBSD data

% Crystal symmetry
cs = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [2.87 2.87 2.87], 'mineral', 'bcc', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [3.65 3.65 3.65], 'mineral', 'fcc', 'color', 'light green'),...
  crystalSymmetry('6/mmm', [2.53 2.53 4.079], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'hcp', 'color', 'light red')};

% Specifying the file to be imported
file_name = 'LPBF_X30MnAl19-1_AsBuilt_ParalleltoBD';
file_extension = '.ang';

% Creating an EBSD variable containing the data
ebsd_asbuilt_pBD = loadEBSD([working_directory,file_name,file_extension],cs,...
                   'interface',strip(file_extension,'left','.'),'convertEuler2SpatialReferenceFrame');

%% Importing the 12.5% strained orhthogonal-to-BD (oBD) EBSD data

% Crystal symmetry
cs = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [3.66 3.66 3.66], 'mineral', 'fcc', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', 'bcc', 'color', 'light green'),...
  crystalSymmetry('6/mmm', [2.53 2.53 4.079], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'hcp', 'color', 'light red'),...
  crystalSymmetry('mmm', [5.108 6.777 4.54], 'mineral', 'Fe3C', 'color', 'cyan')};

% Specifying the file to be imported
file_name = 'LPBF_X30MnAl19-1_12.5%Strained_OrthogonaltoBD';
file_extension = '.ctf';

% Creating an EBSD variable containing the data
ebsd_12_oBD = loadEBSD([working_directory,file_name,file_extension],cs,...
              'interface',strip(file_extension,'left','.'),'convertEuler2SpatialReferenceFrame');

%% Importing the 12.5% strained parallel-to-BD (pBD) EBSD data

% Crystal symmetry
cs = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [3.66 3.66 3.66], 'mineral', 'fcc', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', 'bcc', 'color', 'light green'),...
  crystalSymmetry('6/mmm', [2.53 2.53 4.079], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'hcp', 'color', 'light red')};

% Specifying the file to be imported
file_name = 'LPBF_X30MnAl19-1_12.5%Strained_ParalleltoBD';
file_extension = '.ctf';

% Creating an EBSD variable containing the data
ebsd_12_pBD = loadEBSD([working_directory,file_name,file_extension],cs,...
              'interface',strip(file_extension,'left','.'),'convertEuler2SpatialReferenceFrame');

%% Importing the 25% strained orhthogonal-to-BD (oBD) EBSD data

% Crystal symmetry
cs = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', 'bcc', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [3.66 3.66 3.66], 'mineral', 'fcc', 'color', 'light green'),...
  crystalSymmetry('6/mmm', [2.53 2.53 4.079], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'hcp', 'color', 'light red')};

% Specifying the file to be imported
file_name = 'LPBF_X30MnAl19-1_25%Strained_OrthogonaltoBD';
file_extension = '.ctf';

% Creating an EBSD variable containing the data
ebsd_25_oBD = loadEBSD([working_directory,file_name,file_extension],cs,...
              'interface',strip(file_extension,'left','.'),'convertEuler2SpatialReferenceFrame');

%% Importing the 25% strained parallel-to-BD (pBD) EBSD data

% Crystal symmetry
cs = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', 'bcc', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [3.66 3.66 3.66], 'mineral', 'fcc', 'color', 'light green'),...
  crystalSymmetry('6/mmm', [2.53 2.53 4.079], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'hcp', 'color', 'light red')};

% Specifying the file to be imported
file_name = 'LPBF_X30MnAl19-1_25%Strained_ParalleltoBD';
file_extension = '.ctf';

% Creating an EBSD variable containing the data
ebsd_25_pBD = loadEBSD([working_directory,file_name,file_extension],cs,...
              'interface',strip(file_extension,'left','.'),'convertEuler2SpatialReferenceFrame');

%% Importing the 37.5% strained orhthogonal-to-BD (oBD) EBSD data

% Crystal symmetry
cs = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', 'bcc', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [3.66 3.66 3.66], 'mineral', 'fcc', 'color', 'light green'),...
  crystalSymmetry('6/mmm', [2.53 2.53 4.079], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'hcp', 'color', 'light red')};

% Specifying the file to be imported
file_name = 'LPBF_X30MnAl19-1_37.5%Strained_OrthogonaltoBD';
file_extension = '.ctf';

% Creating an EBSD variable containing the data
ebsd_37_oBD = loadEBSD([working_directory,file_name,file_extension],cs,...
              'interface',strip(file_extension,'left','.'),'convertEuler2SpatialReferenceFrame');

%% Importing the 37.5% strained parallel-to-BD (oBD) EBSD data

% Crystal symmetry
cs = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [3.66 3.66 3.66], 'mineral', 'fcc', 'color', 'light blue'),...
  crystalSymmetry('m-3m', [2.866 2.866 2.866], 'mineral', 'bcc', 'color', 'light green'),...
  crystalSymmetry('6/mmm', [2.53 2.53 4.079], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'hcp', 'color', 'light red')};

% Specifying the file to be imported
file_name = 'LPBF_X30MnAl19-1_37.5%Strained_ParalleltoBD';
file_extension = '.ctf';

% Creating an EBSD variable containing the data
ebsd_37_pBD = loadEBSD([working_directory,file_name,file_extension],cs,...
              'interface',strip(file_extension,'left','.'),'convertEuler2SpatialReferenceFrame');
 