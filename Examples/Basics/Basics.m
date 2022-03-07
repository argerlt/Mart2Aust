% ======================================================================= %
%       Basics.m
%       Shows the simplest method for reconstructing the Austenite phase
%       from an EBSD scan of Martensite
% ======================================================================= %

% Create "AF96.ang" file from provided sample data if not already made
if size(dir('Af96.ang'),1) == 0
    data = h5read("../../Resources/EBSD/EBSD.hdf5","/AF96_small/AF001");
    header = h5readatt("../../Resources/EBSD/EBSD.hdf5","/AF96_small",'header');
    f = fopen('AF96.txt','w');
    fprintf(f,header);
    fclose(f);
    writematrix(data','AF96.txt','delimiter',' ','WriteMode','append')
    movefile AF96.txt AF96.ang
end

% load some metadata about folder locations and do a flight check
Mart2Aust_Parent_folder = "C:\Users\agerlt\workspace\Mart2Aust";
meta.Data_folder = Mart2Aust_Parent_folder + filesep +'Resources'+...
    filesep + 'EBSD'+ filesep +'AF96_small';
meta.MTEX_folder = Mart2Aust_Parent_folder + filesep + 'MTEX';
meta.Functions_folder = Mart2Aust_Parent_folder + filesep + 'Mart2Aust';
meta.current_folder = string(pwd);
meta.MTEX_Version = "mtex-5.7.0";
addpath(genpath(meta.Functions_folder));

check_AusRecon_loaded;
try
    disp(['Compatable MTEX version ', check_MTEX_version(meta.MTEX_Version), ' detected'])
catch
    disp('Incompatable or missing version of MTEX.')
    disp('Attempting to load local copy ... ')
    addpath(meta.MTEX_folder + filesep + meta.MTEX_Version);
    startup_mtex
end
clear Aus_Recon_Parent_folders split_loc f data header Mart2Aust_Parent_folder


%TODO: everything above this should 3 three functions:
%   create_default_data
%   whereami
%   flight_check


% Load default reconstruction options and file location 
options = load_options("default");
fname = dir('AF96.ang');
% Use these two choices to generate some other location information
name = split(string(fname.name),'.');
name = name(end-1);
location = [fname.folder, filesep, fname.name];

% Load a Martensite/Austenite Orientation Relationship and MDF halfwidth
% that was previously measured on this scan (remove this line to recalculate
% the values. See some of the other examples for more details on defining
% and calculating the orientation relationship)
options.OR_ksi =  [2.9057   7.7265   8.2255];
options.OR_noise =  0.0281;

% load the ebsd
original_ebsd =  EBSD.load(location, ...
    'convertEuler2SpatialReferenceFrame','setting 2');
% Not all scans list phases in the same order: this command fixes that
ebsd = prep_for_Recon(original_ebsd,options);
clear original_ebsd
% at this point, we have identical scans with the following phase IDs:
% 0 : Unindexed
% 1 : Untransformed Parent (High temperature, or HT)
% 2 : Transformed Child (Low Temperature, or LR)
% 3 : Reconstructed Parent (starts empty) (Reconstructed, or R)

% Determine the Orientation Relationship, if not already given
CS_HT =ebsd.CSList{1};
CS_LT =ebsd.CSList{2};
[ebsd.opt.OR,ebsd.opt.HW,~] = AutoOR_estimation(ebsd,options);
% Determine the Misorientation Distribution function for the LT phase
[ebsd.opt.LT_MDF,ebsd.opt.psi] = calc_LT_MDF(CS_HT, CS_LT, ebsd.opt.HW, ebsd.opt.OR);

% Perform the actual Reconstruction
Recon_ebsd = Call_Reconstruction(ebsd,options);

% Segment the variants using the low temp and high temp maps
variant_int_map = variants(ebsd, Recon_ebsd,options);
