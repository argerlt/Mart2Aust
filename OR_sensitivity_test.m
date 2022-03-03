% Example for how to perfomr a reconstruction on a stakc of EBSD scans

%============================================
% Flight Check
%============================================
% Change this line or everything breaks
Aus_Recon_Parent_folder = "C:\Users\agerlt\workspace\Aus_Recon";
% make struct of where things are
meta.Data_folder = Aus_Recon_Parent_folder + filesep +...
    'EBSD'+ filesep +'AF96_321x';
meta.MTEX_folder = Aus_Recon_Parent_folder + filesep + 'MTEX';
meta.Functions_folder = string(pwd) + filesep + 'Functions';
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
clear Aus_Recon_Parent_folder

% ================================================================ %
%   Load default reconstruction options and ebsd file locations    %
% ================================================================ %
%%%%%%%   Load options  %%%%%%%%
options = load_options("default");
options.OR_sampling_size = 500;
options.OR_fit_TolFun = 1e-4;
options.OR_fit_TolX = 1e-4;
%%%%%%%   Create task list of EBSD scans  %%%%%%%%
% Make a list of the EBSD text files you want to run through AusRecon
fnames = dir(meta.Data_folder + '\'+'*.ang');
fnames = fnames(1:10);
%%%%%%%% Build Tasks structure for storing inputs and oututs %%%%%%%
% NOTE: this isn't required, "Tasks" iteslf never gets passed into
% functions, it's just convenient for storing data, assuming Tasks doesnt
% get so big it causes OOM errors. (For big stuff, just load and save one
% at a time)
for i = length(fnames):-1:1
    name = split(string(fnames(i).name),'.');
    Tasks(i).name  = name(end-1);
    Tasks(i).location = [fnames(i).folder, filesep, fnames(i).name];
    Tasks(i).options = options;
    Tasks(i).stage = 0;
end
% remove persistant variables
clear name i options fnames data_foldername;

tic
for i = 1:length(Tasks)

    % not all scans are created equal. attempt to homogenize phase names
    original_ebsd =  EBSD.load(Tasks(i).location,...
        'convertEuler2SpatialReferenceFrame'...
        ,'setting 2');
    reformatted_ebsd = prep_for_Recon(original_ebsd,Tasks(i).options);
    Tasks(i).ebsd = reformatted_ebsd;
    clear original_ebsd reformatted_ebsd
    % at this point, we have identical scans with the following phase IDs:
    %-1 : Unindexed
    % 1 : Untransformed Parent (High temperature, or HT)
    % 2 : Transformed Child (Low Temperature, or LR)
    % 3 : Reconstructed Parent (starts empty) (Reconstructed, or R)
    % 4 : Variant Child(starts empty) (Variant, or R)

    % Determine the Orientation Relationship, if not already given
    CS_HT =Tasks(i).ebsd.CSList{1};
    CS_LT =Tasks(i).ebsd.CSList{2};
    [OR,HW,metadata] = AutoOR_estimation(Tasks(i).ebsd,Tasks(i).options);
    % save those values
    Tasks(i).ebsd.opt.OR = OR;
    Tasks(i).ebsd.opt.HW = HW;
end
time = toc;