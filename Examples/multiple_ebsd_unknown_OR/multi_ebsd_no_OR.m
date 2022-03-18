% ======================================================================= %
%       Multi_ebsd_no_OR.m
%       Example for reconstructing multiple EBSD scans where the 
%       orientation relationship is not known a priori
% ======================================================================= %

% Create Several "AF96.ang" files from provided sample data if not already 
% available.
% NOTE: users wishing to use their own data should skip this section, add
% their own ebsd files to this directory, and skip to the "flight check"
% section of the code.
if size(dir('AF96_001.ang'),1) == 0
    header = h5readatt("../../Resources/EBSD/EBSD.hdf5","/AF96_small",'header');
    for i = 1:9
        dname = "/AF96_small/AF" + num2str(i,'%.3d');
        fname = "AF96_" + num2str(i,'%.3d')+".txt";
        angname = "AF96_" + num2str(i,'%.3d')+".ang";
        data = h5read("../../Resources/EBSD/EBSD.hdf5",dname);
        f = fopen(fname,'w');
        fprintf(f,header);
        writematrix(data',fname,'delimiter',' ','WriteMode','append')
        fclose(f);
        movefile(fname, angname,"f")
    end

end



%============================================
% Flight Check
%============================================
% Change this line or everything breaks
Mart2Aust_Parent_folder = "C:\Users\agerlt\workspace\Mart2Aust";
% make struct of where things are
% for the sake of explicit clarity, either place a copy of MTEX 5.7.0 
% inside the Mart2Aust "MTEX" folder, or change the following line to the
% correct absolute path for MTEX 5.7.0
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
clear Aus_Recon_Parent_folder

% ================================================================ %
%   Load default reconstruction options and ebsd file locations    %
% ================================================================ %
%%%%%%%   Load options  %%%%%%%%
options = load_options("default");
%%%%%%%   Create task list of EBSD scans  %%%%%%%%
% Make a list of the EBSD text files you want to run through AusRecon
fnames = dir('*.ang');
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

% ==================================================================== %
% OR values from previous calculation (represents most of the runtime) %
% ==================================================================== %
% Tasks(1).options.OR_ksi =  [2.9057   7.7265   8.2255];
% Tasks(2).options.OR_ksi =  [2.9057   7.7714   8.1489];
% Tasks(3).options.OR_ksi =  [3.0124   8.2746   8.6362];
% Tasks(4).options.OR_ksi =  [2.8199   8.3107   8.6233];
% Tasks(5).options.OR_ksi =  [ 4.8555  9.6552   9.9376];
% Tasks(6).options.OR_ksi =  [3.0825   7.9865   8.4129];
% Tasks(7).options.OR_ksi =  [ 3.0624  7.8404   8.2899];
% Tasks(8).options.OR_ksi =  [ 2.7885  9.7678   9.8262];
% Tasks(9).options.OR_ksi =  [ 3.2479  8.3699   8.7926];
% Tasks(10).options.OR_ksi = [ 2.9300  9.0332   9.2412];
% Tasks(1).options.OR_noise =  0.0281;
% Tasks(2).options.OR_noise =  0.0275;
% Tasks(3).options.OR_noise =  0.0333;
% Tasks(4).options.OR_noise =  0.0332;
% Tasks(5).options.OR_noise =  0.0514;
% Tasks(6).options.OR_noise =  0.0310;
% Tasks(7).options.OR_noise =  0.0293;
% Tasks(8).options.OR_noise =  0.0402;
% Tasks(9).options.OR_noise =  0.0304;
% Tasks(10).options.OR_noise = 0.0340;

% ================================================================ %
%                  prep, reconstruct, and segment                  %
% ================================================================ %
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
    % Also the Misorientation Distribution function for the LT phase
    [LT_MDF,psi] = calc_LT_MDF(CS_HT, CS_LT, HW, OR);
    % save those values
    Tasks(i).ebsd.opt.OR = OR;
    Tasks(i).ebsd.opt.HW = HW;
    Tasks(i).OR_metadata = metadata;
    Tasks(i).ebsd.opt.LT_MDF = LT_MDF;
    Tasks(i).ebsd.opt.psi = psi;
    Tasks(i).stage = 2;
    clear psi OR metadata M LT_MDF HW

    % Perform the actual Reconstruction
    Tasks(i).Recon_ebsd = Call_Reconstruction(Tasks(i).ebsd,Tasks(i).options);

    % Segment the variants using the low temp and high temp maps
    Tasks(i).variant_int_map = variants(Tasks(i).ebsd, ...
        Tasks(i).Recon_ebsd,Tasks(i).options);
end