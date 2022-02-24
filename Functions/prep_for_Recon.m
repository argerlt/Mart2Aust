function ebsd = prep_for_Recon(ebsd,options,HT_Id,LT_Id)
%NOTE: CHANGE NAMES FROM PRE AND POST TO HT, LT, AND R
% UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if string(options.material) == "Steel"
    if exist('HT_Id','var') && exist('LT_Id','var')
        ebsd = steel_prep(ebsd,options,HT_Id,LT_Id);
    else
        ebsd = steel_prep(ebsd,options);
    end
elseif string(options.material) == "Titanium"
    error("Error: Titanium not yet implemented")
else
    error("Error: The only valid options.material values are 'Titanium' and 'Steel'")
end
% add specimen symmetry
ebsd.opt.SS = specimenSymmetry(options.specimen_symmetry);
% The EBSD was sucessfully loaded and prepped. update the progress marker
ebsd.opt.step = "Loaded";
end

function CSList = make_Steel_CSList(options)
% Build the CSList from values in options
% just seperated out to avoid workspace clutter
assert(options.Low_Temp_phase_symm == "m-3m","Steel Reconstruction has only been tested using m-3m symmetry (comment out this line to suppress this warning)")
assert(options.High_Temp_phase_symm == "m-3m","Steel Reconstruction has only been tested using m-3m symmetry (comment out this line to suppress this warning)")
%Someday, we should code in other cases (IE, ferrite instead of Mart,
%different symmetries, etc) but for now, everything becomes 432


% need to add switch here later to allow for grabbing unit cell dimensions
% from ebsd scan or from options input (annoying switch/catch loop, do
% later)
M = options.Low_Temp_lattice_parameters;
A = options.High_Temp_lattice_parameters;
CS_HT = crystalSymmetry(options.High_Temp_phase_symm);
CS_HT.mineral = options.High_Temp_phase_name;
CS_HT.color = options.High_Temp_phase_color;
CS_HT.axes = vector3d([A(1),0,0],[0,A(2),0],[0,0,A(3)]);

CS_LT = crystalSymmetry(options.Low_Temp_phase_symm);
CS_LT.mineral = options.Low_Temp_phase_name;
CS_LT.color = options.Low_Temp_phase_color;
CS_LT.axes = vector3d([M(1),0,0],[0,M(2),0],[0,0,M(3)]);

CS_R = crystalSymmetry(options.High_Temp_phase_symm);
CS_R.mineral = options.Reconstructed_phase_name;
CS_R.color = options.Reconstructed_phase_color;
CS_R.axes = vector3d([A(1),0,0],[0,A(2),0],[0,0,A(3)]);

CS_V = crystalSymmetry(options.Low_Temp_phase_symm);
CS_V.mineral = options.Variant_phase_name;
CS_V.color = options.Variant_phase_color;
CS_V.axes = vector3d([M(1),0,0],[0,M(2),0],[0,0,M(3)]);

CSList = {CS_HT,CS_LT,CS_R,CS_V,'notIndexed'};
end

function [HT_Id,LT_Id] = determine_HT_LT_Steel(old_CSList,options)
%scan the CSList for mineral names or cell dimensions that make sense
% NOTE: At some point, some edge case WILL break this searcher. at that
% point, people are going to have to start choosing the HT/LT phases
% themselves, or write their edge case into this loop.

% grab the searchable data
n = length(old_CSList);
names = strings(1,n); % Used to store phase names
cell_dims = zeros(1,n); % Used to store cell dimensions
phase_linspace = 1:n; 
for i = 1:n
    try
        [a,b,c] =old_CSList{i}.axes.double;
        names(1,i) = string(old_CSList{i}.mineral);
        cell_dims(1,i) = mean(a+b+c);
    catch
    end
end

% First try searching by phase names. Any phase that matches, remove from
% remaining searches
% Try Martensite names
LT_Id = min(phase_linspace(contains(names,'art') == 1));
if ~isempty(LT_Id)
    names(LT_Id)= "";
    cell_dims(LT_Id) = 0;
end
%Try Austenite names
HT_Id = min(phase_linspace(contains(names,'ust') == 1));
if ~isempty(HT_Id)
    names(HT_Id)= "";
    cell_dims(HT_Id) = 0;
end

% If either of those failed, time to instead try matching phase options
% with the expected cell dimensions taken from the options object

% Martensite first
if isempty(LT_Id)
    LT_abc = options.Low_Temp_lattice_parameters; % target LT dimensions
    [~,LT_Id] = min((cell_dims-LT_abc(1)).^2);
    names(LT_Id)= "";
    cell_dims(LT_Id) = 0;
end
% Then Austenite
if isempty(HT_Id)
    HT_abc = options.High_Temp_lattice_parameters; % target HT dimensions
    [delta,HT_Id] = min((cell_dims-HT_abc(1)).^2);
    if delta >0.05
        HT_Id = 60;
    end
end

% because of how the last two dimension searches are written to find the
% minimum delta, these will ALWAYS assign one phase to Martensite and one
% to Austenite. Users are still responsible for verifying this was done
% correctly however.
end


function ebsd = steel_prep(ebsd,options,HT_Id,LT_Id)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
phaseIds = ebsd.phaseId+100;
phaseMap = ebsd.phaseMap;
CSList = make_Steel_CSList(options);
% If HT_Id and LT_Id are explicitly given, things are easy:
if exist('HT_Id','var') && exist('LT_Id','var')
    assert(isinteger(LT_Id),"LT_Id must be an integer value")
    assert(isinteger(HT_Id),"HT_Id must be an integer value")
    assert(ismember(LT_Id,phaseMap),'LT_Id must be an existing phaseId')
    if ismember(LT_Id,phaseMap)
    else
        HT_Id = 60;
    end
else
% Otherwise, scan the Phase data for mineral names or cell dimensions that 
% make sense.
[HT_Id,LT_Id] = determine_HT_LT_Steel(ebsd.CSList,options);
end

% using the HT_Id and LT_Id, fix the phaseIds array
phaseIds(phaseIds == (100+HT_Id)) = 1;
phaseIds(phaseIds == (100+LT_Id)) = 2;
phaseIds(phaseIds > 80) = 0;

% At this point, we want to rearrange the ebsd phase data to match the
% standardized format. first delete old data in such a way that MTEX wont
% complain
ebsd.phaseId = phaseIds*0;
ebsd.phaseMap = [0];
%Then repopulate
ebsd.CSList = CSList;
ebsd.phaseMap = [1,2,3,4,0];
ebsd.phaseId = phaseIds;

% at this point, we have identical scans with the following phase IDs:
%-1 : Unindexed
% 1 : Untransformed Parent
% 2 : Transformed Child
% 3 : Reconstructed Parent (starts empty)
% 4 : Idealized Variants (starts empty)
% 0 : Not Indexed (will also include phases like pearlite or cemetite that aren't part of the resonstruction)
end






