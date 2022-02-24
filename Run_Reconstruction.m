% Script for running the Austenite reconstruction algorithm
% NOTE: all the comments here are NOT for final release, they are just
% either personal reminders, or Notes for/by Austin/Steve.

%============================================
% Flight Check
%============================================
clear all
close all
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

%============================================
% Setting up Recon Jobs
%============================================
%%%%%%%   Load options  %%%%%%%%
% when we run AusRecon, we pass around a LOT of options that have had a
% habit of getting hard-coded into places they shouldn't be. instead,
% Lets just make an options structure, where people can load a default,
% alter what they want, then save their own defaults
options = load_options("default");
% Can load other options like this:
options = load_options("debug");
% then change default values once loaded like this:
options.OR_ksi = [3.09,8.10,8.48];
options.OR_noise = 1.7*degree;
% or create multiple option structures for different expected jobs, or edit
% them on the fly in a for loop, or whatever else. For now, delete and load
% defaults (no given OR, no auto_OR plots, no txt_out, no segmentation)
clear options
options = load_options();

% Options should NOT be edited by any following default functions. Also,
% metadata should NOT be written to it; that belongs in the EBSD struct,
% similar to how grain size is done in MTEX

%%%%%%%   Create task list of EBSD scans  %%%%%%%%
% Make a list of the EBSD text files you want to run through AusRecon
fnames = dir(meta.Data_folder + '\'+'*.ang');
fnames = fnames(1:1);
% delete this last line later, for now just grabs the twinned grain to test
% twinning edge cases
%fnames(6) = dir('../EBSD/AF96_Large/4D-XIII-A_cleaned.ang');
% Use that list to make a non-scalar struct object for where files are, what
% they are named, which options file to use, and recording what stage they
% made it to before exiting (change as needed)
for i = 1:length(fnames)
    name = split(string(fnames(i).name),'.');
    Tasks(i).name  = name(end-1);
    Tasks(i).location = [fnames(i).folder, filesep, fnames(i).name];
    Tasks(i).options = options;
    Tasks(i).stage = 0;
end

% remove persistant variables
clear name i options fnames data_foldername;
%============================================
% Pre-Recon Setup
%============================================
%%%%%%%   Load scans into MTEX EBSD structures  %%%%%%%%
% NOTE: For now, I am storing the EBSD objects in the Tasks struct. THIS IS
% A REALLY BAD PRACTICE FOR LARGE TESTS!!! DO NOT COPY THIS!!!!!!! This is
% for DEBUGGING. Normally, you save these ebsd objects as .mat or .txt
for i = 1:length(Tasks)
    % NOTE: at this point, there are just too many ebsd files to generalize
    % a loading step. this code assumes you find a method for loading your
    % files using MTEX. here is a link to a good starting point:
    % https://mtex-toolbox.github.io/EBSDImport.html
    % let users choose how they want to import the ebsd with the MTEX
    % loader. Trying to predict every possible fringe scenario is a suckers
    % game (both Alex and I wasted literal weeks on this)
    original_ebsd =  EBSD.load(Tasks(i).location,...
        'convertEuler2SpatialReferenceFrame'...
        ,'setting 2');
    % NOTE: if the command above gives results that you cannot easily align
    % with the old examples, here is the older (but maybe less correct? way
    % to load them):
    %original_ebsd =  EBSD.load(Tasks(i).location,'wizard')

    %%%%%%%   SECOND MAJOR PROBLEM FUNCTION   %%%%%%%%
    % EBSD files come in all shapes and sizes. We need to make a single style
    % where we say which phase is parent, which is child, and fix the ordering
    % of them.
    reformatted_ebsd = prep_for_Recon(original_ebsd,Tasks(i).options);
    Tasks(i).ebsd = reformatted_ebsd;
    Tasks(i).stage = 1;
end
clear original_ebsd reformatted_ebsd i

% at this point, we have identical scans with the following phase IDs:
%-1 : Unindexed
% 1 : Untransformed Parent (High temperature, or HT)
% 2 : Transformed Child (Low Temperature, or LR)
% 3 : Reconstructed Parent (starts empty) (Reconstructed, or R)

%============================================
% DETERMINING ORIENTATION RELATIONSHIP AND CALCULATING THE LT_MDF
%============================================
% Since I've already ran this for the five scans of interest, load the
% already calculated values here so the next step auto-skips the OR part

% ------
%Tasks(1).options.OR_ksi = [2.9606    7.8468    8.2958];
%Tasks(2).options.OR_ksi = [2.9311    7.9114    8.2708];
%Tasks(3).options.OR_ksi = [3.0933    11.9621   11.9633];
%Tasks(4).options.OR_ksi = [ 2.8688   8.3493    8.6775];
%Tasks(5).options.OR_ksi = [  4.0538  7.9958    8.8122];
%Tasks(1).options.OR_noise = 0.0278;
%Tasks(2).options.OR_noise = 0.0281;
%Tasks(3).options.OR_noise = 0.0311;
%Tasks(4).options.OR_noise =0.0336;
%Tasks(5).options.OR_noise = 0.0299;
%Tasks(6).options.OR_ksi = [  2.7651  9.5816    9.7156];
%Tasks(7).options.OR_ksi = [  3.4433  8.2166    8.6507];
%Tasks(6).options.OR_noise = 0.0361;
%Tasks(7).options.OR_noise = 0.0375;
% ------
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

for i = 1:length(Tasks)
    CS_HT =Tasks(1).ebsd.CSList{1};
    CS_LT =Tasks(1).ebsd.CSList{2};
%    try
        % First calculate the correct OR and HW
        % NOTE: ask steve for why the HW is calculated during this step
        [OR,HW,metadata] = AutoOR_estimation(Tasks(i).ebsd,Tasks(i).options);
        Tasks(i).ebsd.opt.OR = OR;
        Tasks(i).ebsd.opt.HW = HW;
        Tasks(i).OR_metadata = metadata;
        %    Tasks(i).ebsd.opt.OR_metadata = metadata;
        % Use those values to find the MDF which will be used for populating the
        % out of plane weights
        [LT_MDF,psi] = calc_LT_MDF(CS_HT, CS_LT, ...
            Tasks(i).ebsd.opt.HW,...
            Tasks(i).ebsd.opt.OR);
        Tasks(i).ebsd.opt.LT_MDF = LT_MDF;
        Tasks(i).ebsd.opt.psi = psi;
        Tasks(i).stage = 2;
 %   catch
        disp('beans!!!')
    end
%end
clear psi OR metadata M LT_MDF i HW

% %============================================
% % PRIOR AUSTENITE RECONSTRUCTION
% %============================================
for i = 1:length(Tasks)
    Tasks(i).Recon_ebsd = Call_Reconstruction(Tasks(i).ebsd,Tasks(i).options);
    [Aus_gb,S3,S9] = find_twin_gbs(Tasks(i).Recon_ebsd,3);

    in = Tasks(i).ebsd;
    out = Tasks(i).Recon_ebsd;
    figure()
    plot(out(out.phaseId == 3),out(out.phaseId == 3).orientations)
    hold on
    plot([S3 S9], 'linecolor','w','linewidth',2)
    plot(S3,'linecolor','r','linewidth',0.5,'displayName','S3 Twin')
    plot(S9,'linecolor','b','linewidth',0.5,'displayName','S9 Twin')
    plot(Aus_gb.boundary,'linecolor','k','linewidth',1,'displayName','Prior')
    saveas(gcf,[Tasks(i).location(1:end-4) '_twins_out.png'])
    close all
    %    [grains,Tasks(i).Recon_ebsd.grainId] = ...
    %        calcGrains(Tasks(i).Recon_ebsd('indexed'),'angle',1*degree);

    % 2/18/22 Note: I timed 10 of these, took 1879 seconds if you skipp the
    % autoOR, . approx. twice as long otherwise (accidentally deleted exact
    % timing). 743 of those seconds were spend on nfsoftmex calls, so either
    % need to make them more efficient, or intelligently plan the calls so
    % fewer are made but more data is sent per call. (or just ignore this all
    % together and move on with our lives)

end
% 
% %%
% % ===== Restart from here for Segmentation Troubleshooting ===== %
% clear all
% close all
% load misc/Post_OR_10_recon_pass.mat
% % reset options in case they change
% for i = 1:length(Tasks)
%     Tasks(i).options = load_options;
% end
% 
% %============================================
% % VARIANT SEGMENTATION
% %============================================
% for i = length(Fnames)
%     % THIS ALSO NEEDS A MAJOR REWRITE. I have spent easily 100 hours trying
%     % to understand this function, we need to just state out loud what it
%     % does and how it does it, and write it out on the board as well
%     Tasks(i).segmented = Segmentation(Tasks(i).ebsd, ...
%         Tasks(i).Recon_ebsd,Tasks(i).options);
% end
% 
% %============================================
% % Done
% %============================================
% % Time to go get lunch
% 



