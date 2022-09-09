function out = recondata(name,varargin)
% Load data provided with Mart2Aust and used in the examples
%
% Syntax
%   recondata        % displays a list of available loading routines
%   recondata name   % loads specified data set
%
% NOTE: this is a near carbon copy of MTEX 5.7.0 mtexdata.m, just pointing
% at a different database

% copy the list of available Recon_dataset files into the mtex data folder
list_name = fullfile(mtexDataPath,"Recon_summary.txt");
url = 'https://raw.githubusercontent.com/mesoOSU/Reconstruction_Data/main/Recon_summary.txt';
disp(' Checking for available files at')
disp(['  <a href="' url '">' url '</a>'])
websave(list_name,url);


% Make list of available ebsd files
list = readtable(fullfile(mtexDataPath,'Recon_summary.txt'),'ReadRowNames',true);
type2var = containers.Map({'PoleFigure', 'EBSD', 'grain2d'}, {'pf','ebsd','grains'});

% if/then for if no input provided (ie, list the inputs available)
if nargin < 1
    if nargout == 0
        disp('available loading routines for mtex sample data');
        disp(list)
    else
        out = list.Properties.RowNames;
    end
    return
end

% Catch for if nonexistant data is requested
if isempty(strmatch(name,list.Properties.RowNames))
    warning('mtex:missingData','data not found, please choose one listed below')
    disp(list)
    return
end

type = char(list(name,:).type);
% change warning to error to make it catchable
w = warning('error','MATLAB:load:cannotInstantiateLoadedVariable');

% try to load as mat file
try

    matFile = fullfile(mtexDataPath,[ lower(name) '.mat']);
    load(matFile,'out');

catch
    % If no .mat exists, check to see if .ang or .ctf exists
    fName = fullfile(mtexDataPath,type,char(list(name,:).files));
    if isempty(dir(fName))
        % If not, download the original and save it in the mtex data folder
        url = ['https://media.githubusercontent.com/media/mesoOSU/Reconstruction_Data/main/' type '/' char(list(name,:).files)];
        disp('  downloading data from ')
        disp(' ');
        disp(['   <a href="' url '">' url '</a>'])
        disp(' ');
        disp(['  and saving it to ',fName]);
        websave(fName,url);
    end
    % raw data is downloaded, but must now be loaded into a .mat instance.
    % This means accessing data the first time is slow, but much faster from
    % then on.

    % Each file is different, so do switch cases for each file.
    switch lower(name)

        %Switch case for all the AF96 datasets
        case {'steellarge','steel_1','steel_2','steel_3','steel_4','steel_5',...
                'steel_6','steel_7','steel_8','steel_9','steel_10'}

            CS = {...
                'notIndexed',...
                crystalSymmetry('m-3m', [3.65 3.65 3.65], 'mineral', 'Austenite', 'color', [0.5647 0.5647 0.5647]),...
                crystalSymmetry('I4/mmm', [2.87 2.87 2.87], 'mineral', 'Martensite', 'color', [0.5451 0 0])...
                };
            out = EBSD.load(fName,...
                'CS',CS,...
                'ColumnNames', {'Euler 1' 'Euler 2' 'Euler 3' 'X' 'Y' 'IQ' 'CI' 'Phase' 'SEM_signal' 'Fit'}, ...
                'convertEuler2SpatialReferenceFrame'...
                );

        
        case  {'tilarge', 'ti_1'}

            CS = {...
                'notIndexed',...
                crystalSymmetry('m-3m', [3.192 3.192 3.192], 'mineral', 'HT Beta', 'color', [0.5294    0.8078    0.9804]),...
                crystalSymmetry('6/mmm', [2.954 2.954 4.729], 'mineral', 'LT Alpha', 'color', [0.5608    0.7373    0.5608])...
                };
            out = EBSD.load(fName,...
                'CS',CS,...
                'ColumnNames',{'Phase' 'X' 'Y' 'Bands' 'Error' 'Euler 1' 'Euler 2' 'Euler 3' 'MAD' 'BC' 'BS'}, ...
                'convertEuler2SpatialReferenceFrame'...
                );
    end
    disp([' saving data to ' matFile])
    save(matFile,'out');
end

if nargout == 0
    assignin('base',type2var(type),out);
    if ~check_option(varargin,'silent')
        evalin('base',type2var(type));
    end
    clear out;
end

% restore warning style
warning(w);

end
