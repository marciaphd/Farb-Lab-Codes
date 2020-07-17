%function remapAnalysisV2()
%remapAnalysisV2
% TODO
% put animals names on parameter plots

addpath(genpath('/Users/scott/Dropbox/Documents/MATLAB/jsonlab'));
addpath(genpath('/Users/scott/Dropbox/Documents/MATLAB/FileExchange/comp_struct'));

nargin=0
if nargin == 0
    project_name = 'PPG';
    experiment_name = 'Remapping';
    drug_name = 'TB';
    cogStats = {'Y', 'AU', 'AI'};
    %project_name = 'BuiltIn_Test';
    %experiment_name = 'Test-Remapping';
    %drug_name = 'DrugA';
end

tic
debugMode = 0;
verbose = 0;
PCMINE1E2CORRELATION = 0.8;
PCMINE1E2CORRELATION = 0.6;
PYRCELLMINMEANFRALLENV = 0.1; % 0.1 Hz minimum mean FR in each environment
PYRCELLMAXMEANFR = 5;  % maximum FR in any environment
env1 = 1; env2 = 4;  % significant correlation between 2 environments establishes placefield


% Read file locations from configuration file
config=loadjson('/Users/scott/dropbox/mer_data/configuration/config.json');

% Load json database file of experimental details
experimentDB = loadjson(config.experimentsDB,'SimplifyCell',1);

% Load file that shows progress of analysis by file
epAnalysisStatus = loadjson(config.epAnalysisStatusDB,'SimplifyCell',1);
%config.merAnalysisFirstPassResultsPath = '/Users/scott/Dropbox/MER_Data/Test2/';  %***change this*** and use config struct

% Load CA1/3 info
ca1_3 = loadjson('/Users/scott/Dropbox/MER_Data/db/ppgCA1_3.json','SimplifyCell',1);
for i=1:length(ca1_3)
    messy_name = ca1_3(i).Datafile(1:end-4);
    messy_name(messy_name==' ')='';
    ca1_3(i).filingName = messy_name;
end


%% Determine experimentDB indices for project and experiment
indx_project = find(strcmp(project_name, {experimentDB.Project}));
indx_experiment = find(strcmp(experiment_name,{experimentDB.Experiment}));
indx_analyze = find(strcmp('y', {experimentDB.Analyze}));
indx_project_experiment_analyze = intersect(intersect(indx_project,indx_experiment),indx_analyze);


%% Main loop
% each dose in separate file
aggData = struct([]);
parameters = {'aveFRtableAllCells','aveFRtable','peakFRtable','sictable','pfSizetable'};
parameter_tables = {'MeanFRbyenvAllbinsSmoothed', 'FRmaxenvsmoothed', 'InfoperSpike', 'pfSize'};

parameter_tables = {'MeanFRbyenvAllbinsSmoothed', 'FRmaxenvsmoothed', 'InfoperSpike', 'pfSize', 'aveFRallPyrCells'};
%parameter_tables = parameter_tables(1:2);
for cogStat = cogStats
    aggDataIncr = length(aggData)+1;
    aggData(aggDataIncr).cogStat = cogStat;
    animalsAgg = {};
    fprintf('\nCogStat: %s\n', char(cogStat))
    %avePyrFRIncr = length(avePyrFR)+1;
    indx_cogStat_all = find(strcmp({experimentDB.Status}, char(cogStat)));
    indx_cogStat = intersect(indx_cogStat_all,indx_project_experiment_analyze);
    animals = unique({experimentDB(indx_cogStat).Animal});
    
    % ensure merAnalysisFirstPass results exist
    animals = merAnalysisFirstPassResultsExist(animals,experimentDB,epAnalysisStatus,indx_project_experiment_analyze,project_name);
    fprintf(' Animals: %s\n',mat2str(char(animals)))
    aggData(aggDataIncr).animals = mat2str(char(animals));
    
    for dr_order = 0:1  % 0=Veh, 1=TB
        fprintf(' dr_order: %i\n',dr_order);
        
        for ca_loc = [1 3]
            cnt_ca = 0;
            fprintf('  CA type: %i\n',ca_loc)
            for parameter_table = parameter_tables
                fprintf('  Parameter: %s\n', char(parameter_table))
            %for parameter = 1
                
                tmp_agg = [];  % array to hold aggregate table
                for animalIndx = 1:length(animals)
                    animal = animals{animalIndx};
                    fprintf('   Animal: %s\n',char(animal));
                    animalsAgg = [animalsAgg char(animal)];
                    tmp = [];  % array to hold growing table for specific datafile
                    
                    indx_expDB = findIndex(animal,dr_order, drug_name, project_name, experiment_name, experimentDB);
                    assert(length(indx_expDB)==1)  % single datafile for animal
                    if debugMode;fprintf('   epDB index: %i\n', indx_expDB);end
                    %if isempty(indx_expDB);stop;continue;end
                    fprintf('   Datafile: %s\n', experimentDB(indx_expDB).Datafile)
                    [data,filingName] = loadMerAnalysisFirstPassResults(experimentDB,indx_expDB);
                    assert(strcmp(num2str(dr_order), data.d_rOrderIndx));
                    
                    ca13ByTrode = ca13PositionByTrode(filingName,ca1_3);
                    
                    fprintf('   %i units in merAnalysisFirstPass results\n', length({data.autoSortUnits.MeanFRbyenvAllbinsSmoothed}))
                    indx_lowpass = lowPassFilter(data,PYRCELLMAXMEANFR,debugMode);
                    indx_nonNoiseClusters = removeNoiseClusters(data,debugMode);
                    
                    %env1 = 1; env2 = 4;  % significant correlation between 2 environments establishes placefield
                    indx2 = allExclusions(data,PYRCELLMAXMEANFR,PCMINE1E2CORRELATION,env1,env2,debugMode);
                    
                    if strcmp(char(parameter_table),'aveFRallPyrCells')
                        indx_keep = intersect(indx_lowpass,indx_nonNoiseClusters);
                    else
                    indx_keep = indx2.include;
                    end
                    fprintf('   %i units survived after exclusions\n',length(indx_keep));
                    %indx_keep =
                    %intersect(indx_lowpass,indx_nonNoiseClusters);  % this
                    %is for all pyr cells
                    %fprintf('   %i units survived after removing interneurons and noise clusters\n',length(indx_keep));
                    
                    for i = 1:length(indx_keep)
                        unit_name = data.autoSortUnits(indx_keep(i)).name;
                        trode = str2double(extractBetween(unit_name,'T','_'));
                        ca_type = ca13ByTrode(trode);
                        if ca_type == ca_loc
                            cnt_ca = cnt_ca + 1;
                            % the (1:4) is a kludge, because some files erroneously
                            % have 5 environments, probably should merge to 4.
                            
                            if strcmp(char(parameter_table),'aveFRallPyrCells')
                                tmp = [tmp; data.autoSortUnits(indx_keep(i)).MeanFRbyenvAllbinsSmoothed(1:4)];
                            else
                                tmp = [tmp; data.autoSortUnits(indx_keep(i)).(char(parameter_table))(1:4)];
                            end 
                            
                        end
                    end
                    
                    fprintf('   %i CA%i units found\n',cnt_ca, ca_loc)
                    tmp_agg = [tmp_agg; tmp];
                end  % animal
            %end  % parameter
            
            % save parameter mean, sem, n
            %par = 'aveFRallPyrCells';
            aggData(aggDataIncr).(['mean_' char(parameter_table) '_CA' num2str(ca_loc) '_D' num2str(dr_order)]) = mean(tmp_agg);
            aggData(aggDataIncr).(['sem_' char(parameter_table) '_CA' num2str(ca_loc) '_D' num2str(dr_order)]) = sem(tmp_agg);
            aggData(aggDataIncr).(['n_' char(parameter_table) '_CA' num2str(ca_loc) '_D' num2str(dr_order)]) = size(tmp_agg, 1);
            end  % parameter
        end  % ca
    end  % dr_order
end  % cogStat


%% Plot all 4 parameters on a page for each specific cogStat
% Each subplot: single parameter, single cogStat, all doses, all CA1/3, all
% environments.
yAxisLabels = {'Spikes/sec', 'Spikes/sec', 'Bits/spike', 'cm^2'};
inputs.xlabel = 'Environment';
inputs.xAxisTickLabels = {'F', 'N', 'N', 'F'};
inputs.legend = {'CA1-Veh','CA1-TB','CA3-Veh','CA3-TB'};
%inputs.pageFooterLeft = [num2str(PYRCELLMINMEANFRALLENV) 'Hz~' num2str(PYRCELLMAXMEANFR) 'Hz, r=' num2str(PCMINE1E2CORRELATION)];
%inputs.pageFooterLeft = [num2str(PYRCELLMINMEANFRALLENV) 'Hz~' num2str(PYRCELLMAXMEANFR) 'Hz, r_' num2str(env1) num2str(env2) '=' num2str(PCMINE1E2CORRELATION)];
inputs.pageFooterLeft = [num2str(PYRCELLMINMEANFRALLENV) 'Hz~' num2str(PYRCELLMAXMEANFR) 'Hz, r_{' num2str(env1) ',' num2str(env2) '}=' num2str(PCMINE1E2CORRELATION)];
inputs.pageFooterRight = ['~/Results2/agg/PPG'];

% Each row defines a group (=env), each column are elements in group
% first col: Veh-CA1, TB-CA1, Veh-CA3, TB-CA3 for env 1
for plt_page = 1:3
    inputs.figNum = 1 + plt_page;
    inputs.pageTitle = [project_name ':' experiment_name ':' drug_name ':' char(aggData(plt_page).cogStat)];
    inputs.pltName = ['/Users/scott/Dropbox/MER_Data/Results2/agg/PPG/TB_SingleConcentration_' char(aggData(plt_page).cogStat) '.png'];
    figure(inputs.figNum);clf;
    parameter_indx = 0;
    %for parameter_table = parameter_tables{1:4}  % particular subplot
    for parameter_table_indx = 1:length(parameter_tables)-1  % particular subplot
        parameter_table = parameter_tables(parameter_table_indx);
        parameter_indx = parameter_indx + 1;
        inputs.ylabel = char(yAxisLabels(parameter_indx));
        % Build subplot
        inputs.pltTitles = char(parameter_table);
        plt.means = zeros(4,4);plt.sems = zeros(4,4);plt.ns = zeros(1,4);
        ss_indx = 0;
        for ca = [1 3]
            for dr_order = 0:1
                ss_indx = ss_indx + 1;
                plt.means(:,ss_indx) = (aggData(plt_page).(['mean_' char(parameter_table) '_CA' num2str(ca) '_D' num2str(dr_order)]))';
                plt.sems(:,ss_indx) = (aggData(plt_page).(['sem_' char(parameter_table) '_CA' num2str(ca) '_D' num2str(dr_order)]))';
                plt.ns(ss_indx) = (aggData(plt_page).(['n_' char(parameter_table) '_CA' num2str(ca) '_D' num2str(dr_order)]))';
            end
        end
        barChartErrorbarsV2(plt, parameter_indx, inputs);
    end
    print('-dpng', inputs.pltName)
end



%% Plot average FR for all pyr cells on single page
% 1 page, subplot for each cog status
% Set plot parameters
inputs.figNum = 1;
inputs.xlabel = 'Environment';
inputs.ylabel = 'Firing Rate (Hz)';
inputs.xAxisTickLabels = {'F', 'N', 'N', 'F'};
inputs.legend = {'CA1-Veh','CA1-TB','CA3-Veh','CA3-TB'};
inputs.pltName = ['/Users/scott/Dropbox/MER_Data/Results2/agg/PPG/TBaveFRAllPyrCells' '.png'];
inputs.pageTitle = [project_name ':' experiment_name ':' drug_name '--Average Firing Rate, All Pyramidal Cells'];
inputs.pageFooterLeft = [num2str(PYRCELLMINMEANFRALLENV) 'Hz~' num2str(PYRCELLMAXMEANFR) 'Hz'];
inputs.pageFooterRight = ['~/Results2/agg/PPG/TBaveFRAllPyrCells.png'];
% Make plot titles
% for i = 1:length(avePyrFR)
%     plt(i).pltTitles = [avePyrFR(i).cogStat '-- ' mat2str(char(avePyrFR(i).animals))];
% end

%barChartErrorbars(plt, inputs);


% Generate data for plot
% each row defines a group (=env), each column are elements in group
% first col: Veh-CA1, TB-CA1, Veh-CA3, TB-CA3 for env 1
figure(inputs.figNum);clf;
for plt_indx = 1:3
    inputs.pltTitles = [char(aggData(plt_indx).cogStat) '-- ' aggData(plt_indx).animals];
    plt.means = zeros(4,4);plt.sems = zeros(4,4);plt.ns = zeros(1,4);
    ss_indx = 0;
    for ca = [1 3]
        for dr_order = 0:1
            ss_indx = ss_indx + 1;
            plt.means(:,ss_indx) = (aggData(plt_indx).(['mean_' 'aveFRallPyrCells' '_CA' num2str(ca) '_D' num2str(dr_order)]))';
            %plt.sems(:,ss_indx) = (aggData(1).(['sem_' par '_CA' num2str(ca) '_D' num2str(dr_order)]))';
            plt.sems(:,ss_indx) = (aggData(plt_indx).(['sem_' 'aveFRallPyrCells' '_CA' num2str(ca) '_D' num2str(dr_order)]))';
            plt.ns(ss_indx) = (aggData(plt_indx).(['n_' 'aveFRallPyrCells' '_CA' num2str(ca) '_D' num2str(dr_order)]))';
        end
    end
    barChartErrorbarsV2(plt, plt_indx, inputs);
end
print('-dpng', inputs.pltName)


%% Built in test
%if strcmp(project_name, 'BuiltIn_Test')
    aggData_test = aggData(1);  % avoid comparing Nan, []
    %save ('/Users/scott/Dropbox/MER_Data/Test2/remapAnalysis_test.mat', 'aggData_test')
    test_struct = load('/Users/scott/Dropbox/MER_Data/Test2/remapAnalysis_test.mat');
    [same, er1, er2] = comp_struct(test_struct.aggData_test,aggData_test);
    if isempty(er1) && isempty(er2)
        fprintf('***Built-in test passed***\n')
    else
        er1
        er2
        error('Built-in test failed')
    end
%end

fprintf('Wow, that took %5.2f seconds!\n', toc)

function indx = lowPassFilter(data,PYRCELLMAXMEANFR,debugMode)
% find index of units with peak firing rate less than PYRCELLMAXMEANFR in all environments
indx = [];
for i = 1:length(data.autoSortUnits)
    if isempty(find(data.autoSortUnits(i).MeanFRbyenvAllbinsSmoothed > PYRCELLMAXMEANFR))
        %if isempty(find(data.autoSortUnits(i).FRmaxenvsmoothed > PYRCELLMAXMEANFR))
        indx = [indx i];
    end
end
if debugMode;fprintf('%i units have average overall FR exceeding %i Hz\n',...
        length({data.autoSortUnits.MeanFRbyenvAllbinsSmoothed}) - length(indx),PYRCELLMAXMEANFR);end
end

%function indx = placeCellFilter(data,PCMINE1E2CORRELATION,debugMode)
function indx = placeCellFilter(data,PCMINE1E2CORRELATION,env1,env2,debugMode)
% Find index of units that are valid place cells
indx_sansNoiseClusters = removeNoiseClusters(data,debugMode);  % remove noise clusters first-rho undefined for empty cluster 
%env1 = 1; env2 = 4;  % PUT IN BEGINNING, AND DISPLAY***
indx = [];
%for i = 1:length(data.autoSortUnits)
    for i = 1:length(indx_sansNoiseClusters)
        indx_ = indx_sansNoiseClusters(i);
                        %unit_name = data.autoSortUnits(indx_keep(i)).name;
    %if data.autoSortUnits(i).rho(1,2) > PCMINE1E2CORRELATION
    if data.autoSortUnits(indx_).rho(env1,env2) > PCMINE1E2CORRELATION
        indx = [indx indx_];
        %if debugMode;fprintf('++Unit %s E1E2 corr: %5.2f\n',data.autoSortUnits(i).name,data.autoSortUnits(i).rho(1,2));end
        if debugMode;fprintf('++Unit %s E%iE%i corr: %5.2f\n',data.autoSortUnits(indx_).name,env1,env2,data.autoSortUnits(indx_).rho(env1,env2));end
        else
            %if debugMode;fprintf('--Unit %s E1E2 corr: %5.2f\n',data.autoSortUnits(i).name,data.autoSortUnits(i).rho(1,2));end
            if debugMode;fprintf('--Unit %s E%iE%i corr: %5.2f\n',data.autoSortUnits(indx_).name,env1,env2,data.autoSortUnits(indx_).rho(env1,env2));end
    end
end
end

function indx = removeNoiseClusters(data,debugMode)
indx = [];
for i = 1:length(data.autoSortUnits)
    c = strsplit(data.autoSortUnits(i).name,'_');
    assert(length(c) == 2)
    if ~strcmp('1',c{2})  % not noise cluster
        indx = [indx i];
    end
end
if debugMode;fprintf('%i units are noise clusters\n',...
        length({data.autoSortUnits.MeanFRbyenvAllbinsSmoothed}) - length(indx));end
end

function indx = highPassFilter(data)
% find index of units with mean overall FR > 0.1 Hz in at least one env.
indx = [];
for i = 1:length(data.autoSortUnits)
    if length(find(data.autoSortUnits(i).MeanFRbyenvAllbinsSmoothed < 0.1)) < 4
        indx = [indx i];
    end
end
end

function indx = allExclusions(data,PYRCELLMAXMEANFR,PCMINE1E2CORRELATION,env1,env2,debugMode)
indx.sansNoiseClusters = removeNoiseClusters(data,debugMode);
indx.lowPass = lowPassFilter(data,PYRCELLMAXMEANFR,debugMode);
indx.pcs = placeCellFilter(data,PCMINE1E2CORRELATION,env1,env2,debugMode);
indx.highPass = highPassFilter(data);
%indx.include = intersect(intersect(indx.lowPass,indx.pcs), indx.highPass);
indx.include = intersect(intersect(intersect(indx.lowPass,indx.pcs), indx.highPass),indx.sansNoiseClusters);
end

function indx = findIndex(animal,dr_order, drug, project, experiment, experimentDB)
a = find(strcmp({experimentDB.Animal},char(animal)));
dr_ = find(strcmp({experimentDB.dr_order},num2str(dr_order)));
drg = find(strcmp({experimentDB.Drug},drug));
pr = find(strcmp({experimentDB.Project},project));
ex = find(strcmp({experimentDB.Experiment},experiment));
indx = intersect(intersect(intersect(intersect(a, dr_), drg), pr), ex);
end

function animals = merAnalysisFirstPassResultsExist(animals,experimentDB,epAnalysisStatus,indx_project_experiment_analyze,project_name)
% ensure merAnalysisFirstPass results exist
animals_new = {};
for animal = animals
    indxs_animal = find(strcmp({experimentDB.Animal}, char(animal)));
    indxs_animal = intersect(indxs_animal, indx_project_experiment_analyze);
    datafilenames = {experimentDB(indxs_animal).Datafile};
    for d = datafilenames
        if ~isempty(find(strcmp(char(d),{epAnalysisStatus.Datafile})))
            %if exist(epAnalysisStatus(find(strcmp(char(d),{epAnalysisStatus.Datafile}))).Analyed)
            if strcmp(epAnalysisStatus(find(strcmp(char(d),{epAnalysisStatus.Datafile}))).Analyed, 'y')
                animals_new = [animals_new; animal];
            else
                fprintf('  No merAnalysisFirstPass results for %s\n', char(d))
            end
        elseif strcmp(project_name, 'BuiltIn_Test')
            animals_new = [animals_new; animal];
        end
    end
end
animals_new;
unique(animals_new);
animals = unique(animals_new);
end





function ca13ByTrode = ca13PositionByTrode(filingName,ca1_3)
% detemine CA1 or CA3 for each tetrode of given datafile
indx_CA13 = find(strcmp(filingName, {ca1_3.filingName}));
%fprintf('CA13DB index: %i\n', indx_CA13)
ca1=str2num(char(strsplit(ca1_3(indx_CA13).CA1)));
ca3=str2num(char(strsplit(ca1_3(indx_CA13).CA3)));
ca13ByTrode=zeros(24,1);
ca13ByTrode(ca1)=1;
ca13ByTrode(ca3)=3;
assert(isempty(intersect(ca1, ca3)))
end

function [data,filingName] = loadMerAnalysisFirstPassResults(experimentDB,indx_expDB)
%fprintf('experimentDB index: %i\n', indx_expDB)
[~,fileName] = fileparts(experimentDB(indx_expDB).Datafile);
filingName = fileName;
filingName(filingName==' ') = '';
data_in = load(['/Users/scott/Dropbox/MER_Data/FirstPassAnalysisResults/' filingName '.mat']);
data = data_in.data;
assert(strcmp(filingName, data.datafile));  % ensure read correct merAnalysisFirstPass results
end