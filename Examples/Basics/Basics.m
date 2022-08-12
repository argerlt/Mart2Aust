% ======================================================================= %
%       Basics.m
%       Shows the simplest method for reconstructing the Austenite phase
%       from an EBSD scan of Martensite
% ======================================================================= %
tic

% load some metadata about folder locations and do a flight check
Mart2Aust_Parent_folder = "C:\Users\ALEBRUS\Desktop\Software\Matlab\Mart2Aust-main";
meta.Data_folder = Mart2Aust_Parent_folder + filesep +'Resources'+...
    filesep + 'EBSD'+ filesep +'AF96_small';
meta.MTEX_folder = Mart2Aust_Parent_folder + filesep + 'Mart2Aust' + filesep + 'MTEX';
meta.Functions_folder = Mart2Aust_Parent_folder + filesep + 'Mart2Aust';
meta.current_folder = string(pwd);
meta.MTEX_Version = "mtex-5.7.0";
addpath(genpath(meta.Functions_folder));

% check_AusRecon_loaded;
try
    disp(['Compatable MTEX version ', check_MTEX_version(meta.MTEX_Version), ' detected'])
catch
    disp('Incompatable or missing version of MTEX.')
    disp('Attempting to load local copy ... ')
    addpath(meta.MTEX_folder + filesep + meta.MTEX_Version);
    startup_mtex
end
setMTEXpref('xAxisDirection','east')
setMTEXpref('zAxisDirection','outOfPlane')
clear Aus_Recon_Parent_folders split_loc f data header Mart2Aust_Parent_folder


%TODO: everything above this should 3 three functions:
%   create_default_data
%   whereami
%   flight_check


% Load default reconstruction options and file location 
options = load_options("default");
pname   = 'C:\Users\ALEBRUS\Documents\Microstructure_Work\Exxon_2022\Data';
fname   = '55-HR_EBSD1.ang';
tname   = fullfile(pname,fname);

% Load a Martensite/Austenite Orientation Relationship and MDF halfwidth
% that was previously measured on this scan (remove this line to recalculate
% the values. See some of the other examples for more details on defining
% and calculating the orientation relationship)
%% Alter Default Options
% options.OR_ksi =  [2.7 9.2 9.4];
% % options.OR_ksi   = [2.4284 8.6749 8.8738];
% % options.OR_ksi   = [3.89,8.84,9.17];
% options.OR_noise =  0.029;
% options.RGC_in_plane_m = 6;
% options.RGC_in_plane_b = 12;
% options.RGC_post_pre_m = 1.75e-1;
% options.RGC_post_pre_b = 0.6;
% options.OR_sampling_size = 2000;

options.OR_plot_PAG_Mart_Guess = 1;
options.OR_plot_ODF_of_PAG_OR_guess = 1;

%% Load EBSD and Format for Reconstruction
% Load original ebsd
original_ebsd =  EBSD.load(tname, ...
    'convertEuler2SpatialReferenceFrame','setting 2');

% Truncate ebsd if flagged
truncate = 1;
if truncate
    ebsd = truncate_ebsd(original_ebsd,10,125,7,107);
else
    ebsd = original_ebsd;
end

% Not all scans list phases in the same order: this command fixes that
ebsd = prep_for_Recon(ebsd,options);

% clear original_ebsd
% at this point, we have identical scans with the following phase IDs:
% 0 : Unindexed
% 1 : Untransformed Parent (High temperature, or HT)
% 2 : Transformed Child (Low Temperature, or LR)
% 3 : Reconstructed Parent (starts empty) (Reconstructed, or R)

%% Auto OR Computation

% Determine the Orientation Relationship, if not already given
CS_HT = ebsd.CSList{1};
CS_LT = ebsd.CSList{2};
[ebsd.opt.OR,ebsd.opt.HW,~] = AutoOR_estimation(ebsd,options);

%% Calculate Misorientation Distribution Function

% Determine the Misorientation Distribution function for the LT phase
[ebsd.opt.LT_MDF,ebsd.opt.psi] = calc_LT_MDF(CS_HT, CS_LT, ebsd.opt.HW, ebsd.opt.OR);

%% Perform Reconstruction
tic

[Recon_ebsd, Recon_Likelihood] = Call_Reconstruction_Austin(ebsd,options);

[a,AusGrns,E2627,E112_13] = CalcGrainSize(Recon_ebsd,3,1);

Recon_Time = toc/60;
%%
tic
options.RGC_in_plane_m = 5;
options.RGC_in_plane_b = 3;
options.RGC_post_pre_m = 1.75e-1;
options.RGC_post_pre_b = 2;

% Perform the actual Reconstruction
[Recon_ebsdOrig, Recon_LikelihoodOrig] = Call_Reconstruction(ebsd,options);
Recon_TimeOrig = toc/60;
%% Segment Variants

% Segment the variants using the low temp and high temp maps
variant_int_map = variants(ebsd, Recon_ebsd,options);
