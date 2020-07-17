%function resultsAsGroups=placeCellAnalysis(runMode)
%
% Plot and compute place cells mean and sem (by environment) for
% parameters.
%
% Arguments:
% none
%
% Returns:
% none
%
% File Writes:
% Plots of mean+/-sem by environment for each parameter for all pyr cells
% and for place cells.
% CSV file of mean+/-sem by environment for each parameter for place cells.
% Databases: filenames and data analysis results

tic
addpath(genpath('/Users/scott/Dropbox/Documents/MATLAB/jsonlab'));
nargin = 0
if nargin == 0
    test_mode = 0;
else
    if strcmp(runMode, 'test')
        test_mode = 1;
    end
end

read_from = 'json';

firstpass_results_path = '/Users/scott/Dropbox/MER_Data/FirstPassAnalysisResults/';

% Check SSH connection to web server (and active VPN if off-campus)
fprintf('Checking connection to web server\n')
[s,r]=system('ping -c 1  htep.bumc.bu.edu');
if s > 0
    dlgAnswer = questdlg('No ssh connection to web server, continue?');
    if strcmp(dlgAnswer, 'No'); return; elseif strcmp(dlgAnswer, 'Cancel'); return; end
end

% Load database of filenames
fnameDB = load('/Users/scott/Dropbox/MER_Data/db/AnalysisResultsFilenamesDB.mat');

if strcmp(read_from, 'excel')
%% Read and parse xlsx file
[~,~,c] = xlsread('/Users/scott/Dropbox/MER_Data/db/PyrUnits2.xlsx');
postSortData = struct();
% Build struct of file contents
for i=2:size(c,2)  % row names in first column
    postSortData(i-1).date = c{1,i};
    postSortData(i-1).postsortReviewFile = c{2,i};
    for j=3:size(c(:,i),1)
        if ~isnan(c{j,i})
            postSortData(i-1).units{j-2} = c{j,i};
        end
    end
end
% postSortData
% stop
fprintf('Hand-selected pyr cell database has %i experiments\n', length(postSortData))


% Construct dialog to choose file to analyze
% Based on experiments with hend-selected pyramidal cells.
% Match postsortreview filename in postsortReviewFnameDB with selected pyr
% cells since postsortReviewFnameDB is newer and has only subset.
% firstpass filename is dependent on autosort filename--note, could be
% multiple firstpass filenames for a particular autosort (not done yet).
selector = [];
selector_dlg = [];
selectorIndx = 0;
for i=1:length(postSortData)
    postsortIndx = find(strcmp({fnameDB.postsortReviewFnameDB.postsortReviewFilename},postSortData(i).postsortReviewFile));
    %postSortData(i).postsortReviewFile
    %stop
    if ~isempty(postsortIndx)  % match between database from running postsortreview and database from hand-selected pyr cells
%         selectorIndx = selectorIndx + 1;
%         selector(selectorIndx).postsortReviewFnameDB = postSortData(i).postsortReviewFile;
        autosortIndx = find(strcmp({fnameDB.autosortFnameDB.autosortFilename},...
            fnameDB.postsortReviewFnameDB(postsortIndx).autosortFilename));
        
        if ~isempty(autosortIndx)  % some old stuff might not be in autosort filename database; this will resolve over time
            firstpassIndx = find(strcmp({fnameDB.firstpassFnameDB.autosortFilename},...
                fnameDB.postsortReviewFnameDB(postsortIndx).autosortFilename));
            if ~isempty(firstpassIndx) % ensure that firstpass analysis has been run for this experiment
            selectorIndx = selectorIndx + 1;
        selector(selectorIndx).postsortReviewFilename = postSortData(i).postsortReviewFile;
        selector(selectorIndx).autosortFilename = fnameDB.autosortFnameDB(autosortIndx).autosortFilename;
        selector(selectorIndx).firstpassFilename = fnameDB.firstpassFnameDB(firstpassIndx).firstpassFilename;
        selector(selectorIndx).recordingDate = fnameDB.autosortFnameDB(autosortIndx).recordingDate;
        selector(selectorIndx).datafileNames = fnameDB.autosortFnameDB(autosortIndx).datafileNames;
        selector(selectorIndx).experiment = fnameDB.autosortFnameDB(autosortIndx).experiment;
        selector(selectorIndx).animal = fnameDB.autosortFnameDB(autosortIndx).animal;
        selector(selectorIndx).experimenter = fnameDB.autosortFnameDB(autosortIndx).experimenter;
        %firstpassIndx = find(strcmp({fnameDB.firstpassFnameDB.autosortFilename},fnameDB.autosortFnameDB(autosortIndx).autosortFilename));
        %selector(selectorIndx).firstpassFilename = fnameDB.firstpassFnameDB(firstpassIndx).firstpassFilename;
        selector(selectorIndx).dialog = [selector(selectorIndx).recordingDate '  |  ' selector(selectorIndx).experiment '  |  '...
            num2str(selector(selectorIndx).animal) '  |  ' selector(selectorIndx).experimenter '  |  '...
            selector(selectorIndx).datafileNames '  |  ' selector(selectorIndx).postsortReviewFilename];
        selector_dlg{selectorIndx} = selector(selectorIndx).dialog;
            end
        end
    end 
end
%selector
%selector
%selector_dlg'
selector_dlg = sort(selector_dlg);  % sort by recording date
%selector_dlg'


%%%but, only show if firstpass datafile exists!!!

% popupstr = [];
% for i=1:length({fnameDB.postsortReviewFnameDB.postsortReviewFilename})
%     autosortIndx = find(strcmp({fnameDB.autosortFnameDB.autosortFilename},fnameDB.postsortReviewFnameDB(i).autosortFilename));
%     popupstr(i).recordingDate = fnameDB.autosortFnameDB(autosortIndx).recordingDate;
%     popupstr(i).postsortReviewFilename = fnameDB.postsortReviewFnameDB(i).postsortReviewFilename;
%     popupstr(i).firstpassFilename = fnameDB.postsortReviewFnameDB(i).firstpassFilename;
%     popupstr(i).presentation = [popupstr(i).recordingDate '  |  ' popupstr(i).postsortReviewFilename];
%     %fnameDB.autosortFnameDB.datafileNames
%     popupstr(i).datafileNames = fnameDB.autosortFnameDB(autosortIndx).datafileNames;
%     popupstr(i).presentation2 = [popupstr(i).recordingDate '  |  ' popupstr(i).datafileNames];
% end


%% Popup list of files to select from
set(0, 'DefaultUICOntrolFontSize', 20);  % Kludge: change default font size

% [selectedIndx,v] = listdlg('PromptString', 'Select experiment (recording date | raw datafiles)', 'SelectionMode', 'single',...
%     'ListSize', [1200 500], 'ListString', {selector.dialog});
[selectedIndx,v] = listdlg('PromptString',...
    'Select experiment (recording date | experiment | animal | who | raw datafiles | postsortReview name)',...
    'SelectionMode', 'single',...
    'ListSize', [1200 500], 'ListString', selector_dlg);

if v == 0  % canceled by user
    return
end
%popupstr(selectedIndx)

selected_str = selector_dlg(selectedIndx);
selected_str_split = strsplit(char(selected_str), ' | ');
fname_postsortReviewFilename = char(selected_str_split(end));  % autosort name follows |
fname_postsortReviewFilename = fname_postsortReviewFilename(2:end);  % remove the space following |
psrIndx = find(strcmp({fnameDB.postsortReviewFnameDB.postsortReviewFilename}, fname_postsortReviewFilename));
autosortIndx = find(strcmp({fnameDB.autosortFnameDB.autosortFilename}, fnameDB.postsortReviewFnameDB(psrIndx).autosortFilename));
fname_autosort = fnameDB.autosortFnameDB(autosortIndx).autosortFilename;
firstpassIndx = find(strcmp({fnameDB.firstpassFnameDB.autosortFilename}, fname_autosort));
fname_firstpassFilename = fnameDB.firstpassFnameDB(firstpassIndx).firstpassFilename;
            
elseif strcmp(read_from, 'json')
    
    % hack in json database
    fprintf('Loading hand-selected pyr cells from json database\n')
    
    if ~test_mode
    % Download database from htep
    local_path = '/Users/scott/Dropbox/MER_Data/db/';
    remote_file = 'sdowning@htep.bumc.bu.edu:/var/www/html/mer/selectedPyrUnitsDB.json';
    %scp sdowning@htep.bumc.bu.edu:/var/www/html/mer/selectedPyrUnitsDB.json /Users/scott/Dropbox/MER_Data/db/
    
    [scp_status,scp_cmdout] = system(['scp' ' ' remote_file ' ' local_path]);
    if scp_status ~= 0
        fprintf('***file %s not downloaded from server***\n', remote_file)
        beep
        %return
    end
    end
    
    json_file = loadjson('/Users/scott/dropbox/MER_Data/db/selectedPyrUnitsDB.json');
    lastIndx = length(json_file);
    j55 = loadjson('/Users/scott/dropbox/MER_Data/db/selectedPyrUnitsDBhack.json');
    %json_file(lastIndx+1) = j55
    json_file = [json_file j55];

    
    if ~test_mode
    % Make dialog string
    dlgSelect = [];
    for i=1:length(json_file)
        dlgSelect{i} = [json_file{i}.recordingDate ' | ' json_file{i}.experiment ' | ' json_file{i}.animal ' | '...
            json_file{i}.users   ' | ' json_file{i}.postsortReviewName];
    end
    
    [dlgSelect,map_sort_unsorted] = sort(dlgSelect);  % sort by recording date
    
    % Popup list of files to select from
    set(0, 'DefaultUICOntrolFontSize', 20);  % Kludge: change default font size
    [selectedIndx,v] = listdlg('PromptString', 'Select: recording date | experiment | animal | who |', 'SelectionMode', 'single',...
        'ListSize', [1200 500], 'ListString', dlgSelect);
    
    if v == 0  % canceled by user
        return
    end
    
    
    %json_file{map_sort_unsorted(selectedIndx)}
    json_file_indx = map_sort_unsorted(selectedIndx)
    %stop
    
    % Sanity check
    selected_str = dlgSelect(selectedIndx);
    selected_str_split = strsplit(char(selected_str), ' | ');
    the_date = char(selected_str_split(1));
    assert(strcmp(json_file{map_sort_unsorted(selectedIndx)}.recordingDate, the_date));
    else
        json_file_indx = 20;
    end
    
    selectedPyrUnits = strsplit(json_file{json_file_indx}.units, ',');
    
    fname_postsortReviewFilename = json_file{json_file_indx}.postsortReviewName;
    fname_autosort = json_file{json_file_indx}.autosortName;
    firstpassIndx = find(strcmp({fnameDB.firstpassFnameDB.autosortFilename}, fname_autosort));
    if isempty(firstpassIndx)
        beep
        fprintf('*** Need to run firstpass analysis before running placecell analysis for %s\n', dlgSelect{selectedIndx})
        stop
        return
        %TEST the following code, then integrate above the return
        %statement, then remove return statement if answer yes
        
        dlg_answer = questdlg('Run firstpass analysis now?','','No')
        if strcmp(dlg_answer, 'Yes')
            % Start firstpass analysis
            options.useEnvAligned = 1;  % align environments before determining parameters
            options.useMorphing = 1;  % morph non-square environments to square
            notes = '';  % firstpass notes
            merAnalysisFirstPass4(fname_autosort,notes,options);
            %***check that have correct names now, may have to
            %reloadfnameDB and stuff.
        else
            return
        end
        
    end
    fname_firstpassFilename = fnameDB.firstpassFnameDB(firstpassIndx).firstpassFilename;
    
end   
            
  if ~test_mode          
% fname_postsortReviewFilename = selector(selectedIndx).postsortReviewFilename;
fprintf('Recording date: %s\n', json_file{map_sort_unsorted(selectedIndx)}.recordingDate)
end

fprintf('Postsort review filename: %s\n', fname_postsortReviewFilename);
fprintf('Autosort filename: %s\n', fname_autosort);
% fname_firstpassFilename = selector(selectedIndx).firstpassFilename;
fname_firstpassFilename_absolute = [firstpass_results_path fname_firstpassFilename '.mat'];
fprintf('Loading FirstPass Data from: %s\n', fname_firstpassFilename_absolute)
data = load(fname_firstpassFilename_absolute);
data = data.data;
% data
% stop
% data.autoSortUnits(2)
% stop
assert(strcmp(fname_firstpassFilename, data.firstpassAnalysisName))
%   data
%   data.autoSortUnits(1)
%   stop
  
% Determine index into data for pyr cells
[C,ia,ib]=intersect({data.autoSortUnits.name},selectedPyrUnits);
ia = sort(ia);  % intersect output not in numerical order
pyrUnits = {data.autoSortUnits(ia).name};

% Sanity check, compare this list with those from postSortData(postSortData_indx)
assert(isempty(setdiff(selectedPyrUnits,pyrUnits)))

% Display some things to check of loaded data struct
fprintf('Recording date: %s\n', data.recordingDate);
fprintf('Experiment: %s\n', data.experiment);
fprintf('Animal: %i\n', data.animal);
fprintf('User: %s\n', data.experimenter);
fprintf('Raw datafile names: %s\n', data.rawDatafileNames);

current_date = datestr(now);
%fileName2SaveBase = [char(data.project) '_' char(data.experiment) '_' num2str(data.animal) '_' char(data.drug) '_'...
fileName2SaveBase = [data.project '_' data.experiment '_' num2str(data.animal) '_' data.drug '_'...
    datestr(current_date, 'YYYY-mmm-dd--HH-MM-SS')];

% Determine index into postSortData for specific postsortanalysis file

%postSortData_indx = find(strcmp({postSortData.postsortReviewFile}, fname_postsortReviewFilename))


%% Runtime constants
titles = containers.Map;  % python-like dict
titles('MeanFRbyenvAllbinsSmoothed') = 'Overall Mean Firing Rate';
titles('FRmaxenvsmoothed') = 'Peak Firing Rate';
titles('InfoperSpike') = 'Info Per Spike';
titles('largestSizePFinBins') = 'Place Field Size';
titles('infieldMeans') = 'Mean Infield Firing Rate';
titles('infieldPeaks') = 'Peak Infield Firing Rate';
titles('correlations') = 'Inter-env Correlation';
titles('inField') = 'In-Field Mean Firing Rate';
parameter_tables = {'MeanFRbyenvAllbinsSmoothed', 'FRmaxenvsmoothed', 'InfoperSpike', 'largestSizePFinBins', 'correlations'};
parameter_tables = {'MeanFRbyenvAllbinsSmoothed', 'FRmaxenvsmoothed', 'InfoperSpike', 'largestSizePFinBins', 'correlations', 'inField'};
%labels_xticks = {'F','F','F','F','F','F','F','F','F','F'};
labels_xticks_corr = {'14' '12' '23' '34'};
labels_xticks_corr = {'1-4' '1-2' '2-3' '3-4'};
ylabels = {'spikes/s','spikes/s','bits/spike','cm^2','p','spikes/s','spikes/s'};
if strcmpi(char(data.experiment),'remapping')
    labels_xticks = {'F','N','N','F','F','N','N','F','F','F'};
else
    labels_xticks = {'F','F','F','F','F','F','F','F','F','F'};
end


%%  Loop over parameters to compute stats and plot
arrayIndx = 0;
figure(1);clf
figure(2);clf
stats = struct();
correlations_linear_pyr = [];
correlations_linear_pc = [];

if 1==2
% Determine index into data for pyr cells
%[C,ia,ib]=intersect({data.autoSortUnits.name},postSortData(psr_indxs).units);
[C,ia,ib]=intersect({data.autoSortUnits.name},postSortData(postSortData_indx).units);
ia = sort(ia);  % intersect output not in numerical order
pyrUnits = {data.autoSortUnits(ia).name};

% Sanity check, compare this list with those from postSortData(postSortData_indx)
assert(isempty(setdiff(postSortData(postSortData_indx).units,pyrUnits)))
end

% Determine index into data for place cells (for in-field parameters)
ia2 = [];
for iaa = 1:length(ia)
    if any(data.autoSortUnits(ia(iaa)).numberPFfound) >= 1
    %if data.autoSortUnits(ia(iaa)).numberPFfound >= 1
        ia2 = [ia2 ia(iaa)];
    end
    %tmp2 = cell2mat({data.autoSortUnits(ia2).(char(parameter_table))}');
end
placecellUnits = {data.autoSortUnits(ia2).name};

% Determine index for stable PCs (based on F1F1 corr > 0.5)
ia_stablePCs = [];
ia_UnstablePCs = [];
minCorr = 0.5;
for iaa2 = 1:length(ia2)  % ia2 are place cell indicies
    rho = corr([reshape(data.autoSortUnits(ia2(iaa2)).allPFs2{1},[],1) reshape(data.autoSortUnits(ia2(iaa2)).allPFs2{2},[],1)]);
    if rho(1,2) > minCorr  % was 0.5
        ia_stablePCs = [ia_stablePCs ia2(iaa2)];
    else
        ia_UnstablePCs = [ia_UnstablePCs ia2(iaa2)];
    end
end
PCcondition = ['Select PCs by ' num2str(minCorr) ' min F1F2 corr'];

ia_stablePCs
ia_UnstablePCs

if 1==2
% Determine index for stable PCs (based on aveFR in F1 > 0.27)
ia_stablePCs = [];
ia_UnstablePCs = [];
for iaa2 = 1:length(ia2)  % ia2 are place cell indicies
    if data.autoSortUnits(ia2(iaa2)).MeanFRbyenvAllbinsSmoothed(1) > 1  % was 0.5
        ia_stablePCs = [ia_stablePCs ia2(iaa2)];
        else
        ia_UnstablePCs = [ia_UnstablePCs ia2(iaa2)];
    end
end
    ia_stablePCs
ia_UnstablePCs
end

% ia_stablePCs
% ia_UnstablePCs

% Determine index for pyr cells that are not any type of placecell
ia_pyrNotPC = setdiff(ia,ia2);

% Compute
% tmp = cell2mat({data.autoSortUnits(ia).(char(parameter_table))}');
% bar_means = mean(tmp,1, 'omitnan');
% bar_sems = std(tmp,0,1, 'omitnan')/sqrt(size(tmp,1));

resultsAsGroups = [];
numEnv = length(data.env_xCenterCoordinate);
for parameter_table = parameter_tables
    fprintf('..parameter: %s\n', char(parameter_table))
    tmp = [];  % array to hold pyr data table (unit x environment) for single parameter
    tmp2 = [];  % array to hold pc data table (unit x environment) for single parameter
    tmp3 = [];
    tmp4 = [];
    tmp5 = [];
    arrayIndx = arrayIndx + 1;
    
    if strcmp(parameter_table,'MeanFRbyenvAllbinsSmoothed') | ...
            strcmp(parameter_table,'FRmaxenvsmoothed') | ...
            strcmp(parameter_table,'InfoperSpike')
        
        tmp = cell2mat({data.autoSortUnits(ia).(char(parameter_table))}');
        tmp2 = cell2mat({data.autoSortUnits(ia2).(char(parameter_table))}');
        tmp3 = cell2mat({data.autoSortUnits(ia_stablePCs).(char(parameter_table))}');
        tmp4 = cell2mat({data.autoSortUnits(ia_UnstablePCs).(char(parameter_table))}');
        tmp5 = cell2mat({data.autoSortUnits(ia_pyrNotPC).(char(parameter_table))}');
        
%         resultsAsGroups(arrayIndx).parameter = (char(parameter_table))
%         resultsAsGroups(arrayIndx).allPyrMean = mean(tmp,1, 'omitnan')
%         resultsAsGroups(arrayIndx).allPyrSem = std(tmp,0,1, 'omitnan')/sqrt(size(tmp,1))
%         stop
        
        numPyrCells = size(tmp, 1);
        numPCs = size(tmp2, 1);
        
    elseif strcmp(parameter_table,'correlations')
        % grab correlation data for pyr cells; convert correlation matrix (for each unit) to linear
        correlations_linear_pyr = [];
        for iaa = 1:length(ia)
            number_elements = numel(data.autoSortUnits(ia(iaa)).rho);
            correlations_linear_pyr(iaa,:) = reshape(data.autoSortUnits(ia(iaa)).rho, [1, number_elements]);
        end
        if number_elements == 16
            tmp = correlations_linear_pyr(:,[13,5,10,15]);  % E1E4, E1E2, E2E3, E3E4
        elseif number_elements == 4
            tmp = correlations_linear_pyr(:,3);  % E1E2
        end
        
        % grab correlation data for place cells
        correlations_linear_pyr = [];
        for iaa2 = 1:length(ia2)  % ia2 are place cell indicies
            number_elements = numel(data.autoSortUnits(ia2(iaa2)).rho);
            correlations_linear_pc(iaa2,:) = reshape(data.autoSortUnits(ia2(iaa2)).rho, [1, number_elements]);
        end
        if number_elements == 16
            tmp2 = correlations_linear_pc(:,[13,5,10,15]);  % E1E4, E1E2, E2E3, E3E4
        elseif number_elements == 4
            tmp2 = correlations_linear_pc(:,3);  % E1E2
            %reshape(posFRfiltered{i},[],1)]
            for iaa2 = 1:length(ia2)  % ia2 are place cell indicies
            rho = corr([reshape(data.autoSortUnits(ia2(iaa2)).allPFs2{1},[],1) reshape(data.autoSortUnits(ia2(iaa2)).allPFs2{2},[],1)]);
            rho2(iaa2,1) = rho(1,2);
            aveFR2(iaa2) = data.autoSortUnits(ia2(iaa2)).MeanFRbyenvAllbinsSmoothed(1);
            peakFR2(iaa2) = data.autoSortUnits(ia2(iaa2)).FRmaxenvsmoothed(1);
            infoPerSpike2(iaa2) = data.autoSortUnits(ia2(iaa2)).InfoperSpike(1);
            largestPFsizeBins2(iaa2) = data.autoSortUnits(ia2(iaa2)).largestSizePFinBins(1);
            end
            rho2;
aveFR2;
figure(200);clf
plot(aveFR2,rho2, 'o','MarkerFaceColor', [0.365 0.647 0.855], 'MarkerEdgeColor', 'none')
xlabel('aveFR in F1')
ylabel('F1F2 correlation')
title(['Compare correlation with average FR for ' data.recordingDate ' recording']) 
set(gca, 'FontSize', 16)
grid;box off

figure(201);clf
plot(peakFR2,rho2, 'o','MarkerFaceColor', [0.365 0.647 0.855], 'MarkerEdgeColor', 'none')
xlabel('peakFR in F1')
ylabel('F1F2 correlation')
title(['Compare correlation with peak FR for ' data.recordingDate ' recording']) 
set(gca, 'FontSize', 16)
grid;box off

figure(202);clf
plot(infoPerSpike2,rho2, 'o','MarkerFaceColor', [0.365 0.647 0.855], 'MarkerEdgeColor', 'none')
xlabel('SIC in F1')
ylabel('F1F2 correlation')
title(['Compare correlation with SIC for ' data.recordingDate ' recording']) 
set(gca, 'FontSize', 16)
grid;box off

figure(203);clf
plot(largestPFsizeBins2,rho2, 'o','MarkerFaceColor', [0.365 0.647 0.855], 'MarkerEdgeColor', 'none')
xlabel('F1 placefield size')
ylabel('F1F2 correlation')
title(['Compare correlation with placefield size for ' data.recordingDate ' recording']) 
set(gca, 'FontSize', 16)
grid;box off

%             tmp2
            %stop
        end
        
        %**
        % grab correlation data for stable PCs; convert correlation matrix (for each unit) to linear
        if ~isempty(ia_stablePCs)
        correlations_linear_pyr = [];
        for iaa = 1:length(ia_stablePCs)
            number_elements = numel(data.autoSortUnits(ia_stablePCs(iaa)).rho);
            correlations_linear_pyr(iaa,:) = reshape(data.autoSortUnits(ia_stablePCs(iaa)).rho, [1, number_elements]);
        end
        if number_elements == 16
            tmp3 = correlations_linear_pyr(:,[13,5,10,15]);  % E1E4, E1E2, E2E3, E3E4
        elseif number_elements == 4
            tmp3 = correlations_linear_pyr(:,3);  % E1E2
        end
        else
            tmp3 = 0;
        end
        
        % grab correlation data for unstable PCs; convert correlation matrix (for each unit) to linear
        correlations_linear_pyr = [];
        for iaa = 1:length(ia_UnstablePCs)
            number_elements = numel(data.autoSortUnits(ia_UnstablePCs(iaa)).rho);
            correlations_linear_pyr(iaa,:) = reshape(data.autoSortUnits(ia_UnstablePCs(iaa)).rho, [1, number_elements]);
        end
        if number_elements == 16
            tmp4 = correlations_linear_pyr(:,[13,5,10,15]);  % E1E4, E1E2, E2E3, E3E4
        elseif number_elements == 4
            if ~isempty(correlations_linear_pyr)
            tmp4 = correlations_linear_pyr(:,3);  % E1E2
            else
                tmp4 = 0;
            end
        end
        
        % grab correlation data for pyr not PCs; convert correlation matrix (for each unit) to linear
        correlations_linear_pyr = [];
        %fprintf('length(ia_pyrNotPC): %i\n', length(ia_pyrNotPC))
        %ia_pyrNotPC
        if ~isempty(ia_pyrNotPC)
        for iaa = 1:length(ia_pyrNotPC)
            number_elements = numel(data.autoSortUnits(ia_pyrNotPC(iaa)).rho);
            correlations_linear_pyr(iaa,:) = reshape(data.autoSortUnits(ia_pyrNotPC(iaa)).rho, [1, number_elements]);
        end
        if number_elements == 16
            tmp5 = correlations_linear_pyr(:,[13,5,10,15]);  % E1E4, E1E2, E2E3, E3E4
        elseif number_elements == 4
            tmp5 = correlations_linear_pyr(:,3);  % E1E2
        end
        else
            tmp5 = 0;
        end
        
    elseif strcmp(parameter_table,'inField')
        infieldMeans = [];
        infieldPeaks = [];
        % compute infield mean and peak for place fields
        for loopIndx = 1:length(ia2)  % ia2 are place cell indicies
            for envIndx = 1:numEnv
                infieldFR = data.autoSortUnits(ia2(loopIndx)).allPFs2{envIndx};
                % build table of data (unit x env) to do stats with
                infieldMeans(loopIndx,envIndx) = mean(infieldFR(find(infieldFR>0)));
                
                if isempty(infieldFR(find(infieldFR>0)))
                    infieldPeaks(loopIndx,envIndx) = NaN;  % max function doesn't deal well with empty
                else
                    infieldPeaks(loopIndx,envIndx) = max(infieldFR(find(infieldFR>0)));
                end
            end
        end
        tmp2 = infieldMeans;
        mean(infieldPeaks);
        
        infieldMeans = [];
        for loopIndx = 1:length(ia_stablePCs)  
            for envIndx = 1:numEnv
                infieldFR = data.autoSortUnits(ia_stablePCs(loopIndx)).allPFs2{envIndx};
                infieldMeans(loopIndx,envIndx) = mean(infieldFR(find(infieldFR>0)));
            end
        end
        tmp3 = infieldMeans;    
        
        infieldMeans = [];
        for loopIndx = 1:length(ia_UnstablePCs)  
            for envIndx = 1:numEnv
                infieldFR = data.autoSortUnits(ia_UnstablePCs(loopIndx)).allPFs2{envIndx};
                infieldMeans(loopIndx,envIndx) = mean(infieldFR(find(infieldFR>0)));
            end
        end
        tmp4 = infieldMeans; 
        
    elseif strcmp(parameter_table,'largestSizePFinBins')
        tmp = cell2mat({data.autoSortUnits(ia).(char(parameter_table))}');
        tmp = tmp * 3.5 * 3.5;  % scale by bin area
        tmp2 = cell2mat({data.autoSortUnits(ia2).(char(parameter_table))}');
        tmp2 = tmp2 * 3.5 * 3.5;  % scale by bin area
        tmp3 = cell2mat({data.autoSortUnits(ia_stablePCs).(char(parameter_table))}');
        tmp3 = tmp3 * 3.5 * 3.5;  % scale by bin area
        tmp4 = cell2mat({data.autoSortUnits(ia_UnstablePCs).(char(parameter_table))}');
        tmp4 = tmp4 * 3.5 * 3.5;  % scale by bin area
        tmp5 = cell2mat({data.autoSortUnits(ia_pyrNotPC).(char(parameter_table))}');
        tmp5 = tmp5 * 3.5 * 3.5;  % scale by bin area
    end
    
    % arrange results in struct array (parameter ~ type of cell)
        resultsAsGroups(arrayIndx).parameter = (char(parameter_table));
        resultsAsGroups(arrayIndx).allPyrMean = mean(tmp,1, 'omitnan');
        resultsAsGroups(arrayIndx).allPyrSem = std(tmp,0,1, 'omitnan')/sqrt(size(tmp,1));
        resultsAsGroups(arrayIndx).allPyrn = size(tmp,1);
        
        resultsAsGroups(arrayIndx).allPCsMean = mean(tmp2,1, 'omitnan');
        resultsAsGroups(arrayIndx).allPCsSem = std(tmp2,0,1, 'omitnan')/sqrt(size(tmp2,1));
        resultsAsGroups(arrayIndx).allPCsn = size(tmp2,1);
        
        resultsAsGroups(arrayIndx).stablePCsMean = mean(tmp3,1, 'omitnan');
        resultsAsGroups(arrayIndx).stablePCsSem = std(tmp3,0,1, 'omitnan')/sqrt(size(tmp3,1));
        resultsAsGroups(arrayIndx).stablePCsn = size(tmp3,1);
        
        resultsAsGroups(arrayIndx).unStablePCsMean = mean(tmp4,1, 'omitnan');
        resultsAsGroups(arrayIndx).unStablePCsSem = std(tmp4,0,1, 'omitnan')/sqrt(size(tmp4,1));
        resultsAsGroups(arrayIndx).unStablePCsn = size(tmp4,1);
        
        resultsAsGroups(arrayIndx).pyrNotPCsMean = mean(tmp5,1, 'omitnan');
        resultsAsGroups(arrayIndx).pyrNotPCsSem = std(tmp5,0,1, 'omitnan')/sqrt(size(tmp5,1));
        resultsAsGroups(arrayIndx).pyrNotPCsn = size(tmp5,1);
        %stop
    
    for pltNum = 1:2
        if pltNum == 1
            figure(1)
            bar_means = mean(tmp,1, 'omitnan');
            bar_sems = std(tmp,0,1, 'omitnan')/sqrt(size(tmp,1));
            n = size(tmp,1);
            page_title = 'Pyr cells';
        else
            figure(2)
            bar_means = mean(tmp2,1, 'omitnan');
            bar_sems = std(tmp2,0,1, 'omitnan')/sqrt(size(tmp2,1));
            n = size(tmp2,1);
            page_title = 'Place cells';
            stats(arrayIndx).parameter = char(parameter_table);
            stats(arrayIndx).meanByEnv = bar_means;
            stats(arrayIndx).semByEnv = bar_sems;
            stats(arrayIndx).n = size(tmp2,1); 
            stats(arrayIndx).table = tmp2;
        end
        
        if ~isempty(bar_means)  % prevent graph display if no data
            subplot(2,3,arrayIndx)
            hold on
            plt=bar(bar_means,'FaceColor',[0.365 0.647 0.855], 'EdgeColor', 'none');
            errorbar(bar_means,bar_sems, 'Color', [0.365 0.647 0.855],'linestyle', 'none')
            %             if strcmp(parameter_table,'inField')
            %                 yyaxis right
            %                 plot(mean(infieldPeaks),'o-r')  % add max infield FR to mean infield FR plot
            %                 ylabel('peak (spikes/s)', 'Color', 'r')
            %                 yyaxis left
            %             end
            hold off
            
            grid
            box off
            xticks(1:size(bar_means,2));
            if ~strcmp(parameter_table,'correlations')
                xticklabels(labels_xticks);
            else
                if numEnv == 2
                    xticklabels('1-2');
                else
                    xticklabels(labels_xticks_corr);
                end
            end
            ylabel(ylabels{arrayIndx})
            set(gca, 'FontSize', 12)
            title(titles(char(parameter_table)))
            v=axis;
            axis([.25 size(bar_means,2)+0.75 v(3) v(4)]);
            
            if arrayIndx == 1
                text(0.05,0.95,['n=' num2str(n)],'Units','Normalized')
            end
        end
        if arrayIndx == length(parameter_tables)
            annot=annotation('textbox',[0,0,1,1],'String',page_title,'LineStyle','none','VerticalAlignment','top',...
                'HorizontalAlignment','Center','FontSize',12,'Color','Blue','Interpreter','none');
            annot=annotation('textbox',[0,0,1,1],'String',mfilename,'LineStyle','none','VerticalAlignment','top',...
                'HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
% %             annot=annotation('textbox',[0,0,1,1],'String',fname_postsortReviewFilename,'LineStyle','none',...
% %                 'VerticalAlignment','Bottom','HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
id_str = [num2str(data.animal) '/' data.experiment '/' data.project '/' data.recordingDate '/' data.experimenter];
annot=annotation('textbox',[0,0,1,1],'String',id_str,'LineStyle','none',...
                'VerticalAlignment','Bottom','HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
            annot=annotation('textbox',[0,0,1,1],'String',datestr(now, 'YYYY-mmm-dd'),'LineStyle','none',...
                'VerticalAlignment','Bottom','HorizontalAlignment','Right','FontSize',8,'Interpreter','none');
            
            fileName2save = ['/Users/scott/Dropbox/MER_Data/Results2/SingleFiles/' fileName2SaveBase '_' num2str(pltNum) '.png'];
            picNameStub = ['/Users/scott/Dropbox/MER_Data/Results2/SingleFiles/' fileName2SaveBase '_']; [picNameStub '1.png'];
            if ~test_mode;print('-dpng', fileName2save);end
        end
        %end
    end
end

% numUnits = length(data.autoSortUnits);
% length(autoSortResult.autoSortUnits)-length(trodes);  % first cluster of each trode is noise

%numTetrodes = length(unique(str2num(char(extractBetween({data.autoSortUnits.name}, 'T', '_')))));
numTetrodes = length(unique(char(extractBetween({data.autoSortUnits.name}, 'T', '_'))));  % this change needed for apr 3 recording
numUnits = length(data.autoSortUnits) - numTetrodes;
fprintf('Number of units: %i\n', numUnits)
fprintf('Number of pyr cells: %i\n', numPyrCells)
fprintf('Number of place cells: %i\n', numPCs)
fprintf('Plots and Stats table located at %s\n', '/Users/scott/Dropbox/MER_Data/Results2/SingleFiles/')

%resultsAsGroups

figure(4);clf
for i=1:4
    subplot(2,2,i)
% bar([resultsAsGroups(1).allPyrMean;resultsAsGroups(1).allPCsMean;resultsAsGroups(1).stablePCsMean]')
% title(resultsAsGroups(1).parameter)
%bar([resultsAsGroups(i).allPyrMean;resultsAsGroups(i).allPCsMean;resultsAsGroups(i).stablePCsMean;resultsAsGroups(i).
%unStablePCsMean;resultsAsGroups(i).pyrNotPCsMean]')
hB=bar([resultsAsGroups(i).allPyrMean;resultsAsGroups(i).allPCsMean;resultsAsGroups(i).pyrNotPCsMean;...
    resultsAsGroups(i).stablePCsMean;resultsAsGroups(i).unStablePCsMean;]', 'EdgeColor', 'none')
title(resultsAsGroups(i).parameter)
            %xticks(1:size(bar_means,2));
            xticks(1:size(resultsAsGroups(1).allPyrMean,2));
            if ~strcmp(parameter_table,'correlations')
                xticklabels(labels_xticks);
            else
                if numEnv == 2
                    xticklabels('1-2');
                else
                    xticklabels(labels_xticks_corr);
                end
            end
ylabel(ylabels{i})
grid
            set(gca, 'FontSize', 12)
alpha(.5)  % set face transparency of all bars to permit viewing legend when directly over a bar
box off

hAx=gca;            % get a variable for the current axes handle
% str = {'E1', 'E2'};
% hAx.XTickLabel=str; % label the ticks
hT=[];              % placeholder for text object handles
for j=1:length(hB)  % iterate over number of bar objects
%   hT=[hT text(hB(j).XData+hB(j).XOffset, hB(j).YData, num2str(hB(j).YData.','%.3f'), ...
%                           'VerticalAlignment','bottom','horizontalalign','center')];  % label each bar with magnitude                
                          end

if i == 1
    ns = [resultsAsGroups(1).allPyrn resultsAsGroups(1).allPCsn resultsAsGroups(1).stablePCsn...
        resultsAsGroups(1).unStablePCsn resultsAsGroups(1).pyrNotPCsn];
    % add number of units to each bar (first group only)
    for j=1:length(hB)  % iterate over number of bar objects
    text(hB(j).XData(1)+hB(j).XOffset, hB(j).YData(1), num2str(ns(j)),'VerticalAlignment','bottom','horizontalalign','center')
    end
end
if i == 2
    %text(0.05,0.95,['n=' num2str(n)]
    legend('allPyr','allPCs','pyrNotPCs','stablePCs','unStablePCs','Orientation','horizontal','FontWeight','bold')
    legend('boxoff')
end
end

pageTitle = ['Stable PCs = ' PCcondition];
            annot=annotation('textbox',[0,0,1,1],'String',pageTitle,'LineStyle','none','VerticalAlignment','top',...
                'HorizontalAlignment','Center','FontSize',12,'Color','Blue','Interpreter','none');
            annot=annotation('textbox',[0,0,1,1],'String',mfilename,'LineStyle','none','VerticalAlignment','top',...
                'HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
% %             annot=annotation('textbox',[0,0,1,1],'String',fname_postsortReviewFilename,'LineStyle','none',...
% %                 'VerticalAlignment','Bottom','HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
id_str = [num2str(data.animal) '/' data.experiment '/' data.project '/' data.recordingDate '/' data.experimenter];
annot=annotation('textbox',[0,0,1,1],'String',id_str,'LineStyle','none',...
                'VerticalAlignment','Bottom','HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
            annot=annotation('textbox',[0,0,1,1],'String',datestr(now, 'YYYY-mmm-dd'),'LineStyle','none',...
                'VerticalAlignment','Bottom','HorizontalAlignment','Right','FontSize',8,'Interpreter','none');


figure(5);clf
if 1==2
for i=5:6
    subplot(2,1,i-4)
    hB=bar([resultsAsGroups(i).allPyrMean;resultsAsGroups(i).allPCsMean;resultsAsGroups(i).pyrNotPCsMean;...
    resultsAsGroups(i).stablePCsMean;resultsAsGroups(i).unStablePCsMean;]', 'EdgeColor', 'none')
% pause
% hB=bar([resultsAsGroups(i).allPyrMean;resultsAsGroups(i).allPCsMean;resultsAsGroups(i).pyrNotPCsMean;...
%     resultsAsGroups(i).stablePCsMean;resultsAsGroups(i).unStablePCsMean;], 'EdgeColor', 'none')
% stop
if i == 5
   sameColors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
   hB.FaceColor = 'flat';
   hB.CData(1,:) = sameColors(1,:)
   hB.CData(2,:) = sameColors(2,:)
   hB.CData(3,:) = sameColors(3,:)
   hB.CData(4,:) = sameColors(4,:)
   hB.CData(5,:) = sameColors(5,:)

end
if i ==6
%        legend('allPyr','allPCs','pyrNotPCs','stablePCs','unStablePCs','Orientation','horizontal','FontWeight','bold')
%     legend('boxoff')    
        for j=1:length(hB)  % iterate over number of bar objects
    %text(hB(j).XData(1)+hB(j).XOffset, 0.85*hB(j).YData(1), num2str(ns(j)),'VerticalAlignment','bottom','horizontalalign','center')
    end
end
title(resultsAsGroups(i).parameter)
grid
end
end
if 1==2  % problem in this section with feb 9 older exp
subplot(2,1,1)
hold on  % must be better way to get each bar different color and legend
b1=bar(1,resultsAsGroups(5).allPyrMean, 'EdgeColor', 'none')
errorbar(1,resultsAsGroups(5).allPyrMean,resultsAsGroups(5).allPyrSem, 'Color', [0 0 0],'linestyle', 'none')
b2=bar(2,resultsAsGroups(5).allPCsMean, 'EdgeColor', 'none')
errorbar(2,resultsAsGroups(5).allPCsMean,resultsAsGroups(5).allPCsSem, 'Color', [0 0 0],'linestyle', 'none')
b3=bar(3,resultsAsGroups(5).pyrNotPCsMean, 'EdgeColor', 'none')
errorbar(3,resultsAsGroups(5).pyrNotPCsMean,resultsAsGroups(5).pyrNotPCsSem, 'Color', [0 0 0],'linestyle', 'none')
b4=bar(4,resultsAsGroups(5).stablePCsMean, 'EdgeColor', 'none')
errorbar(4,resultsAsGroups(5).stablePCsMean,resultsAsGroups(5).stablePCsSem, 'Color', [0 0 0],'linestyle', 'none')
b5=bar(5,resultsAsGroups(5).unStablePCsMean, 'EdgeColor', 'none')
errorbar(5,resultsAsGroups(5).unStablePCsMean,resultsAsGroups(5).unStablePCsSem, 'Color', [0 0 0],'linestyle', 'none')
hold off
ylabel('r')
title(resultsAsGroups(5).parameter)
xticks('')
grid
    legend([b1,b2,b3,b4,b5],'allPyr','allPCs','pyrNotPCs','stablePCs','unStablePCs','Orientation','horizontal','FontWeight','bold')
    legend('boxoff')
    box off
    set(gca, 'FontSize', 12)
alpha(.5)  % set face transparency of all bars to permit viewing legend when directly over a bar

subplot(2,1,2)
i=6;
% hB=bar([resultsAsGroups(i).allPyrMean;resultsAsGroups(i).allPCsMean;resultsAsGroups(i).pyrNotPCsMean;...
%     resultsAsGroups(i).stablePCsMean;resultsAsGroups(i).unStablePCsMean;]', 'EdgeColor', 'none')
% hB=bar([[0 0];resultsAsGroups(i).allPCsMean;[0 0];...
%     resultsAsGroups(i).stablePCsMean;resultsAsGroups(i).unStablePCsMean;]', 'EdgeColor', 'none')  % this works but no errorbars

bars = [[0 0];resultsAsGroups(i).allPCsMean;[0 0];...
     resultsAsGroups(i).stablePCsMean;resultsAsGroups(i).unStablePCsMean]
 hB=bar(bars', 'EdgeColor', 'none')
end

 if 1==2
%hold on  % must be better way to get each bar different color and legend
b1=bar(1,resultsAsGroups(i).allPCsMean, 'EdgeColor', 'none')
%errorbar(1,resultsAsGroups(i).allPyrMean,resultsAsGroups(5).allPyrSem, 'Color', [0 0 0],'linestyle', 'none')
b2=bar(2,resultsAsGroups(i).allPCsMean, 'EdgeColor', 'none')
errorbar(2,resultsAsGroups(i).allPCsMean,resultsAsGroups(5).allPCsSem, 'Color', [0 0 0],'linestyle', 'none')
b3=bar(3,resultsAsGroups(i).pyrNotPCsMean, 'EdgeColor', 'none')
errorbar(3,resultsAsGroups(i).pyrNotPCsMean,resultsAsGroups(5).pyrNotPCsSem, 'Color', [0 0 0],'linestyle', 'none')
b4=bar(4,resultsAsGroups(i).stablePCsMean, 'EdgeColor', 'none')
errorbar(4,resultsAsGroups(i).stablePCsMean,resultsAsGroups(5).stablePCsSem, 'Color', [0 0 0],'linestyle', 'none')
b5=bar(5,resultsAsGroups(i).unStablePCsMean, 'EdgeColor', 'none')
errorbar(5,resultsAsGroups(i).unStablePCsMean,resultsAsGroups(5).unStablePCsSem, 'Color', [0 0 0],'linestyle', 'none')
hold off
 end

        for j=1:length(hB)  % iterate over number of bar objects
    %text(hB(j).XData(1)+hB(j).XOffset, 0.85*hB(j).YData(1), num2str(ns(j)),'VerticalAlignment','bottom','horizontalalign','center')
        end
        ylabel('spikes/s')
    title(resultsAsGroups(i).parameter)
    xticklabels({'F','F'})
grid
box off
set(gca, 'FontSize', 12)
alpha(.5)  % set face transparency of all bars to permit viewing legend when directly over a bar

pageTitle = ['Stable PCs = ' PCcondition];
            annot=annotation('textbox',[0,0,1,1],'String',pageTitle,'LineStyle','none','VerticalAlignment','top',...
                'HorizontalAlignment','Center','FontSize',12,'Color','Blue','Interpreter','none');
            annot=annotation('textbox',[0,0,1,1],'String',mfilename,'LineStyle','none','VerticalAlignment','top',...
                'HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
% %             annot=annotation('textbox',[0,0,1,1],'String',fname_postsortReviewFilename,'LineStyle','none',...
% %                 'VerticalAlignment','Bottom','HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
id_str = [num2str(data.animal) '/' data.experiment '/' data.project '/' data.recordingDate '/' data.experimenter];
annot=annotation('textbox',[0,0,1,1],'String',id_str,'LineStyle','none',...
                'VerticalAlignment','Bottom','HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
            annot=annotation('textbox',[0,0,1,1],'String',datestr(now, 'YYYY-mmm-dd'),'LineStyle','none',...
                'VerticalAlignment','Bottom','HorizontalAlignment','Right','FontSize',8,'Interpreter','none');
            

% return
% stop

%% How many place cells if inter-environment correlation for hand-selected pyr cells must be 0.8?
% Find first subsequent environment that matches the first environment
envs = strsplit(data.environments, ',');
envs_match = 0;
for i=2:length(envs)
    if strcmp(envs{i},envs{1})
        envs_match = i;
        break  % end after match found
    end
end

if envs_match
    fprintf('Environment %i matches the first environment\n',envs_match)
    
    pc_corr8 = [];
    pc_corrs = [];
    diff_sum = [];
    diff_sum_limit = [];  % only keep those with FR > 0.27
    for iaa = 1:length(ia)  %ia = index into struct for pyr cells
        
        pf20s = [];
            for envIndx = 1:numEnv
                peakFR = max(max(data.autoSortUnits(ia(iaa)).posFRfiltered{envIndx}));
                pfMinFR = 0.2*peakFR;
                pf20 = data.autoSortUnits(ia(iaa)).posFRfiltered{envIndx};
                pf20(pf20 < pfMinFR) = 0;
                pf20s{envIndx} = pf20;
            end
            fr1 = mean(mean(pf20s{1}));
            fr2 = mean(mean(pf20s{2}));
            diff_sum(iaa) = (abs(fr1-fr2)) / (fr1 + fr2);
            if data.autoSortUnits(ia(iaa)).MeanFRbyenvAllbinsSmoothed(1) > 0.27
            %if fr1 > 0.27 & fr2 > 0.27
                diff_sum_limit = [diff_sum_limit, (abs(fr1-fr2)) / (fr1 + fr2)];
            end
        
        pc_corr = data.autoSortUnits(ia(iaa)).rho;
        if isequal(size(pc_corr),[2,2])
            if pc_corr(1,2) > 0.8
                pc_corr8 = [pc_corr8 ia(iaa)];
            end
            pc_corrs(iaa) = pc_corr(1,2);
            
%             pf20s = [];
%             for envIndx = 1:numEnv
%                 peakFR = max(max(data.autoSortUnits(ia(iaa)).posFRfiltered{envIndx}));
%                 pfMinFR = 0.2*peakFR;
%                 pf20 = data.autoSortUnits(ia(iaa)).posFRfiltered{envIndx};
%                 pf20(pf20 < pfMinFR) = 0;
%                 pf20s{envIndx} = pf20;
%             end
%             fr1 = mean(mean(pf20s{1}));
%             fr2 = mean(mean(pf20s{2}));
%             diff_sum(iaa) = (abs(fr1-fr2)) / (fr1 + fr2);
%             if data.autoSortUnits(ia(iaa)).MeanFRbyenvAllbinsSmoothed(1) > 0.27
%             %if fr1 > 0.27 & fr2 > 0.27
%                 diff_sum_limit = [diff_sum_limit, (abs(fr1-fr2)) / (fr1 + fr2)];
%             end
            
            
        elseif isequal(size(pc_corr),[4,4])
            if pc_corr(1,envs_match) > 0.8
                pc_corr8 = [pc_corr8 ia(iaa)];
            end
            pc_corrs(iaa) = pc_corr(1,envs_match);
            
        else
            if pc_corr(1,2) > 0.8
                pc_corr8 = [pc_corr8 ia(iaa)];
            end
            pc_corrs(iaa) = pc_corr(1,2);
        end
    end
    fprintf('Number of place cells with inter-environment correlation > 0.8: %i\n', length(pc_corr8))
    diff_sum
    diff_sum_limit
    
    %% Plot the histogram of E1E2 correlation for each place cell
    figure(3);clf
    histogram(pc_corrs,20)
    %xlabel('Place cell E1E2 correlation')
    xlabel(['E1E' num2str(envs_match) ' correlation'])
    ylabel('Frequency')
    title('Pyramidal cell stability')
    ys = get(gca, 'YLim');
    axis([-0.02 1 0 ys(2)+0.5])
    box off
    grid
    annot=annotation('textbox',[0,0,1,1],'String',mfilename,'LineStyle','none','VerticalAlignment','top',...
                'HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
            id_str = [num2str(data.animal) '/' data.experiment '/' data.project '/' data.recordingDate '/' data.experimenter];
annot=annotation('textbox',[0,0,1,1],'String',id_str,'LineStyle','none',...
                'VerticalAlignment','Bottom','HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
            annot=annotation('textbox',[0,0,1,1],'String',datestr(now, 'YYYY-mmm-dd'),'LineStyle','none',...
                'VerticalAlignment','Bottom','HorizontalAlignment','Right','FontSize',8,'Interpreter','none');
    
    
    fileName2save = ['/Users/scott/Dropbox/MER_Data/Results2/SingleFiles/' fileName2SaveBase '_' num2str(3) '.png'];
    if ~test_mode;print('-dpng', fileName2save);end
    
    
    %% Plot the histogram of E1E2 firing rate change for each place cell
    figure(6);clf
    hold on
    histogram(diff_sum,20)
    histogram(diff_sum_limit,20)
    hold off
    legend('all pyr cells','pyr cells FR>0.27')
    legend boxoff
    %xlabel(['Place cell E1E' num2str(envs_match) ' firing rate change'])
    xlabel(['E1E2 firing rate change'])
    ylabel('Frequency')
    title('Pyramidal cell rate remapping')
    ys = get(gca, 'YLim');
    axis([-0.02 1 0 ys(2)+1])
    box off
    grid
    text(0.6,0.7,'Mean FR of PF (20% peak FR)','Units','Normalized')
    text(0.6,0.67,'Method: difference / sum','Units','Normalized')
            
    annot=annotation('textbox',[0,0,1,1],'String',mfilename,'LineStyle','none','VerticalAlignment','top',...
                'HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
            id_str = [num2str(data.animal) '/' data.experiment '/' data.project '/' data.recordingDate '/' data.experimenter];
annot=annotation('textbox',[0,0,1,1],'String',id_str,'LineStyle','none',...
                'VerticalAlignment','Bottom','HorizontalAlignment','Left','FontSize',8,'Interpreter','none');
            annot=annotation('textbox',[0,0,1,1],'String',datestr(now, 'YYYY-mmm-dd'),'LineStyle','none',...
                'VerticalAlignment','Bottom','HorizontalAlignment','Right','FontSize',8,'Interpreter','none');
    
    fileName2save = ['/Users/scott/Dropbox/MER_Data/Results2/SingleFiles/' fileName2SaveBase '_' num2str(4) '.png'];
    if ~test_mode;print('-dpng', fileName2save);end
    %if test_mode;print('-dpng', fileName2save);end  %%*** temporary
    
end


%% Dump stats to CSV file
%StatsTable_filename = ['/Users/scott/Dropbox/MER_Data/Results2/SingleFiles/'  postSortData(psr_indxs).postsortReviewFile '.csv'];
StatsTable_filename = ['/Users/scott/Dropbox/MER_Data/Results2/SingleFiles/'  fileName2SaveBase '.csv'];
if ~test_mode
    fid = fopen(StatsTable_filename, 'w');
else
    fid = 1;
end

fprintf(fid, 'Stats for place cells for %s\n',fname_postsortReviewFilename);
fprintf(fid,'Parameter,');
for i=1:length(stats(1).meanByEnv)
    fprintf(fid,'E%s,',num2str(i));
end
fprintf(fid,'\n');

for i=1:length(stats)
    fprintf(fid,stats(i).parameter);
    fprintf(fid,'\nMeans,');
    fprintf(fid,'%6.4f,',stats(i).meanByEnv);
    fprintf(fid,'\nSems,');
    fprintf(fid,'%6.4f,',stats(i).semByEnv);
    fprintf(fid,'\nn,');
    fprintf(fid,'%i',stats(i).n);
    fprintf(fid,'\n');
end
if fid ~= 1;fclose(fid);end

if 1==2
%% Dump PC table to CSV
% todo
% (1) rename parameters, (2) change to bin size, (3) add column headings
% (correlation different than the rest), (4) maybe add CA1/3, (5) change
% Nan to 0, (6) add column means--at least for a check
fid = 1;
for i=1:length(stats)
    %     size(stats(i).table)
    %     size(placecellUnits)
    %     size(ia2)
    sums = zeros(1,size(stats(i).table,2))
    for j=1:size(stats(i).table,1)
        %fprintf(fid, '%s,', stats(i).parameter);
        %titles
        fprintf(fid, '%s,', titles(stats(i).parameter));
        % put unit number
        fprintf(fid, '%s,', placecellUnits{j});
        % print table row
        for k=1:size(stats(i).table,2)
            %fprintf(fid, '%s\n', mat2str(stats(i).table(j,:)));
            %%%if ~isnan(stats(i).table(j,k))
            fprintf(fid, '%5.2f', stats(i).table(j,k));
            sums(k) = sums(k) + stats(i).table(j,k);
            %%%else
               %%% fprintf(fid, '%5.2f', 0);
                %%%sums(k) = sums(k) + 0;
            %%%end
            %sums(k) = sums(k) + stats(i).table(j,k);
        end
        fprintf(fid, '\n')
        
        %stop
    end
    
    fprintf(fid, '\n')
    %sums
    test_mode99 = 1;
%     if test_mode99
%         fprintf(fid, '%s,%s,', 'Mean', ' ')
%         for k=1:size(stats(i).table,2)
%             fprintf(fid, '%5.2f', sums(k) / size(stats(i).table,1))
%         end
%     end
%     sums / size(stats(i).table,1)
%     mean(stats(i).table)
    mean(stats(i).table,1, 'omitnan')
    %stop
end
stop
end


%%
if test_mode  % don't update database if test mode
    return
end

dlg_answer = questdlg('Update database?', '', 'yes',...
    'no', 'yes');  % note that last string is default selection, thus a repeat
if strcmp(dlg_answer, 'no')
    return
end


%% Update analysis results database

fname = '/Users/scott/Dropbox/MER_Data/db/AnalysisResults.mat';

% Define current record
dbase_record = struct();

dbase_record.recordingDate = data.recordingDate;
dbase_record.postsortReviewFilename = fname_postsortReviewFilename;
%dbase_record.firstpassFilename = fname_firstpassFilename;
dbase_record.firstpassFilename = data.firstpassAnalysisName;


dbase_record.placecellAnalysisFilename = fileName2SaveBase;
dbase_record.placecellAnalysisDate = datestr(current_date, 'YYYY-mmm-dd, HH-MM');
dbase_record.placecellAnalysisNotes = '';
dbase_record.numUnits = numUnits;
dbase_record.numPyramidalCells = numPyrCells;
dbase_record.numPlaceCells = numPCs;
dbase_record.placecellParameterStats = stats;

dbase_record.placecellCriteria = 'pyr cell with PF in at least 1 env';
dbase_record.correlationCriteria = 'pearson correlation of binned firing rate maps (between envs) of placecells';
dbase_record.remappingCriteria = '';

dbase_record.selectedPyramidalUnits = strjoin(pyrUnits,',');
dbase_record.computedPlacecellUnits = strjoin(placecellUnits,',');


% Update and save database
db_tables = whos('-file', fname);
if find(strcmp({db_tables.name}, 'placecellAnalysisResults'))
    %if exist(fname)
    placecellAnalysisResults = load(fname, 'placecellAnalysisResults');
    placecellAnalysisResults = placecellAnalysisResults.placecellAnalysisResults;
    nextIncr = length(placecellAnalysisResults) + 1;
    placecellAnalysisResults(nextIncr) = dbase_record;
    %save(fname, '-append', 'placecellAnalysis')
else
    placecellAnalysisResults = dbase_record;
    %save(fname, 'placecellAnalysis')
end
save(fname, '-append', 'placecellAnalysisResults')


%% Update analysis results filename database
% Database already exists (established in autosort), add to table if it
% exists, otherwise create it.
fname = '/Users/scott/Dropbox/MER_Data/db/AnalysisResultsFilenamesDB.mat';

% Define current record
dbase_record = struct();

dbase_record.postsortReviewFilename = fname_postsortReviewFilename;  % 1 level back
dbase_record.firstpassFilename = data.firstpassAnalysisName;  % 2 levels back

dbase_record.placeCellAnalysisFilename = fileName2SaveBase;
dbase_record.placeCellAnalysisDate = datestr(current_date, 'YYYY-mmm-dd, HH-MM');
dbase_record.placeCellAnalysisNotes = '';
dbase_record.numUnits = numUnits;

dbase_record.numPyramidalCells = numPyrCells;
dbase_record.numPlaceCells = numPCs;
dbase_record.pyrcellAnalysisResults_name = [picNameStub '1.png'];
dbase_record.placecellAnalysisResults_name = [picNameStub '2.png'];
dbase_record.placecellStabilityHistogram_name = [picNameStub '3.png'];
dbase_record.placecellAnalysisStatsResults_name = StatsTable_filename;
dbase_record.placecellRateRemappingHistogram = [picNameStub '4.png'];

dbase_record.placeCellAnalysisFile = mfilename;

[~, dbase_record.placeCellAnalysisComputer] = system('hostname');
dbase_record.placeCellAnalysisRuntimes = toc;

% Update and save database
db_tables = whos('-file', fname);
if find(strcmp({db_tables.name}, 'placeCellAnalysisFnameDB'))
    placeCellAnalysisFnameDB = load(fname, 'placeCellAnalysisFnameDB');
    placeCellAnalysisFnameDB = placeCellAnalysisFnameDB.placeCellAnalysisFnameDB;
    nextIncr = length(placeCellAnalysisFnameDB) + 1;
    placeCellAnalysisFnameDB(nextIncr) = dbase_record;
else
    placeCellAnalysisFnameDB = dbase_record;
end
save(fname, '-append', 'placeCellAnalysisFnameDB')


%% Upload plots to server
floc = 'sdowning@htep.bumc.bu.edu:/var/www/html/mer/paths';
local_path = '/Users/scott/Dropbox/MER_Data/Results2/SingleFiles/';
remote_path = 'sdowning@htep.bumc.bu.edu:/var/www/html/mer/analysisResults/';
for i=1:4
[scp_status,scp_cmdout] = system(['scp' ' ' picNameStub num2str(i) '.png' ' ' remote_path]);
if scp_status ~= 0
    fprintf('***file %s not transfered to server***\n', [picNameStub num2str(i) '.png'])
end
end
%StatsTable_filename
[scp_status,scp_cmdout] = system(['scp' ' ' StatsTable_filename ' ' remote_path]);
if scp_status ~= 0
    fprintf('***file %s not transfered to server***\n', StatsTable_filename)
end


%% Update analysis results webpage
displayAnalysisResultsTable2('n');


%% Update Pipelinelog
task = 'placeCellAnalysis';
description = [char(data.project) '/' char(data.experiment) '/' num2str(data.animal) '/' char(data.experimenter)];
task_output = fname_postsortReviewFilename;
task_inputs = [fname_firstpassFilename '.mat'];
updatePipelineLog(task,description,task_output,task_inputs)


