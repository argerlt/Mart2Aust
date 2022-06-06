% ======================================================================= %
%       Basics.m
%       Shows the simplest method for reconstructing the Austenite phase
%       from an EBSD scan of Martensite
% ======================================================================= %

%% Create "AF96.ang" file from provided sample data if not already made
%NOTE: skip this section if you already have an ang file in this example
%folder you would like to investigate
if size(dir('Af96.ang'),1) == 0
    data = h5read("../../Resources/EBSD/EBSD.hdf5","/AF96_small/AF001");
    header = h5readatt("../../Resources/EBSD/EBSD.hdf5","/AF96_small",'header');
    f = fopen('AF96.txt','w');
    fprintf(f,header);
    fclose(f);
    writematrix(data','AF96.txt','delimiter',' ','WriteMode','append')
    movefile AF96.txt AF96.ang
end

%% load some metadata about folder locations and do a flight check
% Skip this if MTEX 5.7 is already loaded, Mart2Aust has been added to
% your Matlab Path, and you don't need access to the AF96 examples
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

% ======================================================================= %
%                        Reconstruction starts here                       %
% ======================================================================= %
% Load default reconstruction options 
options = load_options("default");
% define the name of the file we are reconstructing
fname = dir('AF96.ang');
% use fname to generate a test name and file location
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
[Recon_ebsd, likelihoods] = Call_Reconstruction(ebsd,options);

% Segment out the subgrains using the low temp and high temp maps
subgrain_IDs = subgrains(ebsd, Recon_ebsd,options);

% Save the results. Its easiest to make EBSD objects from saved .ang files,
% But MTEX's auto-reader breaks when reading it's own ang files back in.
% So, as a fix, we use our own file writer.
M2A_export_ang(Recon_ebsd,'Reconstruction.ang')
% we will also do the same with the original where the phase names have
% been changed
M2A_export_ang(ebsd,'Original.ang')
% Finally, we can save out all the data (including the variant map) as its
% own text file for loading into other programs like HEXRD,Dream3D, or Orix
GenText(name, ebsd, Recon_ebsd, likelihoods, subgrain_IDs)