%=================================================================================================================================%
%=================================================================================================================================%
% Benchmark model for the paper: "Modeling the health impact of gross wastewater contamination events in drinking water networks"

% Step 1. Generate demands based on the STREaM model
% Step 2. Assign new network demands (generated from STREaM) to L-Town.inp
% Step 3. Run contaminations
% Step 4. Calculate downstream population
% Step 5. Calculate infection risk

% The default is that steps 1,2, and 4 are commented so the user can generate the results as seen in the paper. To generate your own results, please uncomment the steps and follow the instructions in the code.

%{
 Copyright (c) 2023 [KWR Water Research Institute] (https://www.kwrwater.nl/en/) and [KIOS Research and Innovation Centre of Excellence, University of Cyprus] (www.kios.org.cy)
 
 Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence");
 You may not use this work except in compliance with the Licence.
 You may obtain a copy of the Licence at: https://joinup.ec.europa.eu/collection/eupl/eupl-text-11-12
 
 Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the Licence for the specific language governing permissions and limitations under the Licence.
 
 Author(s)     : Sotirios Paraskevopoulos [sotirios.paraskevopoulos@kwrwater.nl]
                 Stelios Vrachimis [vrachimis.stelios@ucy.ac.cy]
                 Marios Kyriakou [kiriakou.marios@ucy.ac.cy]
 
 Last revision : November 2023
%}

%% Step 1. Generate Stream demands
try 
d.unload
catch ERR
end 
fclose all;clear class;clear all;clc;close all;
% Start EPANET MATLAB TOOLKIT
addpath(genpath(pwd));
disp('EPANET-MATLAB Toolkit Paths Loaded.'); 
% inpname = 'L-TOWN.inp';
% 
% d = epanet(inpname);
% d.setTimeReportingStep(300)
% d.setTimeHydraulicStep(300)
% d.setTimePatternStep(300)
% hydraulics = d.getComputedTimeSeries;
% 
% %%% Calculate hydraulics
% Demand= hydraulics.Demand(:,d.getNodeJunctionIndex);
% Demand= Demand.*1000;  %CMH to L/h
% K= round(size(Demand, 1)/7);
% Dt= double(d.getTimeHydraulicStep)/3600;
% Vpd= 150;
% V_sc=cumsum(Demand(1:K,:),1)*Dt; %Calculate volume consumed per node in 24 hours
% People_per_node= V_sc(end,:)/Vpd; % People per node 
% 
% %%% ::: STREAM :::
% 
% % AVAILABLE FIXTURES
% % 1. toilet
% % 2. shower
% % 3. faucet
% % 4. clothes washer
% % 5. dishwasher
% % 6. bathtub
% %%% ::: INPUT SETTINGS :::
% %%% ::: LOADING COMPLETE DATABASE :::
% load database.mat
% 
% People_per_node_rnd= round(People_per_node);
% Population_unc=0.1;
% Node=d.getNodeNameID;
% Scenario={};
% % for 100 scenarios change the following for loop (1:100)
% for scen=1
%     scen
%     for i=1:length(People_per_node_rnd)
%         disp(['Simulating Node ',num2str(i),' of ',num2str(length(People_per_node_rnd))])
%         % Create +-10% uncertainty
%         Population=People_per_node_rnd(i);
%         Population_l=Population-Population_unc*Population;
%         Population_u=Population+Population_unc*Population;
%         Population=Population_l+rand(1,length(Population)).*(Population_u-Population_l);
%         Population = round(Population);
% 
%     %input population
%     home=0;
%     %initialization
%     StToilet=0;
%     StShower=0;
%     StFaucet=0;
%     StClothesWasher=0;
%     StDishwasher=0;
%     TOTAL=0;
% 
%     while Population>0
% 
%     % --- A. Household size setting
%     home=home+1;
%     param.HHsize = randi(5,1); % This parameter should be in the interval (1,6).
%     % From 1 to 5, it indicates the number of people living in the house. 6 means ">5".
% 
%     Population=Population-param.HHsize;
%     % --- B. Water consuming fixtures selection
%     % Legend:
%     % 0 = not present
%     % 1 = present
% 
%     param.appliances.StToilet = 1;
%     param.appliances.HEToilet = 0;
% 
%     param.appliances.StShower = 1;
%     param.appliances.HEShower = 0;
% 
%     param.appliances.StFaucet = 1;
%     param.appliances.HEFaucet = 0;
% 
%     param.appliances.StClothesWasher = 1;
%     param.appliances.HEClothesWasher = 0;
% 
%     param.appliances.StDishwasher = 1;
%     param.appliances.HEDishwasher = 0;
% 
%     param.appliances.StBathtub = 1;
%     param.appliances.HEBathtub = 0;
% 
%     % --- C. Time horizon length setting
%     param.H = 3; % It is measured in [days]
% 
%     % --- D. Time sampling resolution
%     param.ts = 30; % It is measured in [10 seconds] units. The maximum resolution allowed is 10 seconds (param.ts = 1).
% 
%     % Setting the seed
%     % rng(1);
% 
%     % Parameters structure settings and check
%     % Checking input consistency
%     temp=checkInput(param);
%     % clearvars -except param
% 
%     %%% ::: WATER END-USE TIME SERIES GENERATION :::
% 
%     % Initialization
%     outputTrajectory = initializeTrajectories(param);
%     % Include the first step
%     appNames = fieldnames(outputTrajectory);
%     for app=appNames'
%         outputTrajectory.(char(app))=zeros(1,length(outputTrajectory.TOTAL)+30);
%     end
% 
%     % End-use water use time series generation
%     outputTrajectory = generateConsumptionEvents(outputTrajectory, param, database);
%     % disp('End-use consumption trajectories created');
% 
%     % Total water use time series aggregation
%     outputTrajectory = sumToTotal(outputTrajectory);
%     % disp('Total consumption trajectory created');
% 
%     % Data scaling to desired sampling resolution
%     outputTrajectory = aggregateSamplingResolution(outputTrajectory, param);
%     % disp('Data scaled to desired sampling resolution');
% 
%     StToilet=outputTrajectory.StToilet+StToilet;
%     StShower=outputTrajectory.StShower+StShower;
%     StFaucet=outputTrajectory.StFaucet+StFaucet;
%     StClothesWasher=outputTrajectory.StClothesWasher+StClothesWasher;
%     StDishwasher=outputTrajectory.StDishwasher+StDishwasher;
%     TOTAL=outputTrajectory.TOTAL+TOTAL;
%     clear outputTrajectory;
%     end
% 
%     output.output(scen).(Node{i}).StToilet=StToilet;
%     output.output(scen).(Node{i}).StShower=StShower;
%     output.output(scen).(Node{i}).StFaucet=StFaucet;
%     output.output(scen).(Node{i}).StClothesWasher=StClothesWasher;
%     output.output(scen).(Node{i}).StDishwasher=StDishwasher;
%     output.output(scen).(Node{i}).TOTAL=TOTAL;
% 
%     end
% 
%     %%% Assign stream demand for total water use
%     for j=1:length(People_per_node_rnd)
%     Stream_demand_tot(:,j)= output.output(scen).(Node{j}).TOTAL;
%     end
%     Stream_demand_tot= (Stream_demand_tot.*12)/1000; %convert from L/5min to CMH
%     Stream.Scenario{scen}= Stream_demand_tot;
% 
%     % Assign stream demand for injestion of water
%     for j=1:length(People_per_node_rnd)
%     Stream_demand_Faucet(:,j)= output.output(scen).(Node{j}).StFaucet;
%     end
%     Stream_demand_Faucet= (Stream_demand_Faucet.*12)/1000; %convert from L/5min to CMH
%     Stream_faucet.Scenario{scen}= Stream_demand_Faucet;
% end
% 
% % Saving
% save ('./Stream_demands.mat','output','Stream','Stream_faucet')


%=================================================================================================================================%
%=================================================================================================================================%

%% Step 2. Assign new network demands
% try 
% d.unload
% catch ERR
% end 
% fclose all;clear class;
% close all;clear all;clc;
% % Start EPANET MATLAB TOOLKIT
% %%% Load Network:
% inpname = 'L-TOWN.inp';
% dispname = 'L-TOWN';
% d=epanet(inpname);
% nj=d.getNodeJunctionCount;
% nn=d.getNodeCount;
% 
% % Assign demands to network:
% %%% Get existing patterns:
% %%%     d.getNodeDemandPatternIndex{categoryIndex}(nodeIndex)
% demInd1 = double(d.getNodeDemandPatternIndex{1}(:)); 
% demInd1G = demInd1;
% demInd2 = double(d.getNodeDemandPatternIndex{2}(:));
% demInd3 = double(d.getNodeDemandPatternIndex{3}(:));
% baseDemInd1 = double(d.getNodeBaseDemands{1}(:));
% baseDemInd2 = double(d.getNodeBaseDemands{2}(:));
% baseDemInd3 = double(d.getNodeBaseDemands{3}(:));
% 
% %%% Zero base demands
% for i=1:nn
%     disp(['Zero base demand ',num2str(i)])
%     d.setNodeBaseDemands(i, 1, 0)
%     d.deleteNodeJunctionDemand(i,2)
%     d.deleteNodeJunctionDemand(i,2)
% end
% 
% %%% Delete three(3) old patterns:
% d.deletePattern(1)
% d.deletePattern(1) % pattern number 2 become 1 after deletion of 1
% d.deletePattern(1)
% 
% %%% Save new input file:
% emptyInpName = ['networks\',dispname,'_empty','.inp'];
% d.saveInputFile(emptyInpName);
% disp('EMPTY NETWORK READY!')
% 
% %%% Load new actual demands:
% load stream_demands_paper.mat % Stream demands as generated in the paper
% % load stream_demands.mat % Uncomment here to use the demands mat file you just created in Step 1
% 
% %%% Calculate and assign new base demands:
% % for 100 scenarios change the following for loop (1:100)
% for m=1
%     d=epanet(emptyInpName);
%     demInd1 = demInd1G;
%     base_Stream_demand = mean(Stream.Scenario{m});
%     for i=1:nj
%         disp(['Assign base demand ',num2str(i)])
%         d.setNodeBaseDemands(i, 1, base_Stream_demand(i))
%     end
% 
% %%% Calculate new patterns:
%     pattern_Stream_demand = Stream.Scenario{m}/base_Stream_demand;
%     pattern_Stream_demand(isnan(pattern_Stream_demand))=0;
% 
%     %%% Add new patterns:
%     for i = 1:nj
%         demands = pattern_Stream_demand(:);
%         resDemPatInd(i)=d.addPattern(['P-Res-',num2str(i)],demands);
%         disp(['Creating pattern Residential ',num2str(i)])
%     end
% 
%     for i=1:nj
%         disp(['Indexing pattern ',num2str(i),' category 1'])
%         if demInd1(i)==0
%             continue 
%         elseif demInd1(i)==1 % Residential
%             demInd1(i)=i;
%         else
%             error('unknown demand pattern')
%         end
%     end
% 
%     %%% Assign new patterns:
%     for i=1:nn
%         disp(['Assigning pattern ',num2str(i),' out of ',num2str(nn)])
%         d.setNodeDemandPatternIndex(i, 1, demInd1(i))
%     end
% 
%     % Correct times:
%     d.setTimeReportingStep(300)
%     d.setTimeHydraulicStep(300)
%     d.setTimePatternStep(300)
% 
% % Save new input file:
%     newInpname = ['networks\',dispname,'_stream',num2str(m),'.inp'];
%     d.saveInputFile(newInpname);
%     disp('NETWORK READY!')
% end

%=================================================================================================================================%
%=================================================================================================================================%

%% Step 3. Run contamination
try 
d.unload
catch ERR
end 
fclose all;clear class;clear all;clc;close all;

% Define msx file
msx={};
msx.FILENAME = './networks/L_town_Microbial_Cont_new.msx';

% Section Title
msx.TITLE = {'Modeling the health impact of gross wastewater contamination events in drinking water networks'};

% Section Options
msx.AREA_UNITS = 'FT2'; %AREA_UNITS FT2/M2/CM2
msx.RATE_UNITS = 'HR'; %TIME_UNITS SEC/MIN/HR/DAY
msx.SOLVER = 'RK5'; %SOLVER EUL/RK5/ROS2
msx.TIMESTEP = 300; %TIMESTEP in seconds
msx.RTOL = 0.001;  %RTOL value
msx.ATOL = 0.001;  %ATOL value

% Section Species
% <type> <specieID> <units> (<atol> <rtol>)
msx.SPECIES(1) = {'BULK CL2 MG 0.01 0.001'}; % Chlorine
msx.SPECIES(2) = {'BULK P CFU 0.01 0.001'}; % Pathogen
msx.SPECIES(3) = {'BULK C_FRA MG 0.01 0.001'}; %Fast reducing agent
msx.SPECIES(4) = {'BULK C_SRA MG 0.01 0.001'}; % Slow reducing agent

% Section Coefficients
% CONSTANT name value % PARAMETER name value
msx.COEFFICIENTS(1) = {'CONSTANT T 12'}; % Water temperature in degree C
msx.COEFFICIENTS(2) = {'CONSTANT Kp 265.8'}; % Inactivation rate
msx.COEFFICIENTS(3) = {'CONSTANT A 1'}; % Aplification factor (dm/h), decimeter Monteiro et al 2020
msx.COEFFICIENTS(4) = {'CONSTANT B 14 '}; % Rate coefficient (L/MG)

% Section Terms
% <termID> <expression>
msx.TERMS(1) = {'Km (1.5826e-04 * RE^0.88 )/ D'}; %Mass transport coeff (dm/h), taken from epanet-msx manual(ft/h)
msx.TERMS(2) = {'KWAL A*EXP(-B*CL2)'};
msx.TERMS(3) = {'K_FAST 0.2808'}; % Monteiro
msx.TERMS(4) = {'K_SLOW 0.0071'}; % Monteiro

% Section Pipes
% EQUIL <specieID> <expression> % RATE <specieID> <expression> % FORMULA <specieID> <expression>
msx.PIPES(1) = {'RATE CL2 -K_FAST*C_FRA*CL2-K_SLOW*C_SRA*CL2-((4/D)*(KWAL/(1+KWAL/Km)))*CL2'};
msx.PIPES(2) = {'RATE P -Kp*P*CL2'};
msx.PIPES(3) = {'RATE C_FRA -K_FAST*C_FRA*CL2'};
msx.PIPES(4) = {'RATE C_SRA -K_SLOW*C_SRA*CL2'};

% Section Tanks
% EQUIL <specieID> <expression> % RATE <specieID> <expression> % FORMULA <specieID> <expression>
msx.TANKS(1) = {'RATE CL2 -K_FAST*C_FRA*CL2-K_SLOW*C_SRA*CL2'};
msx.TANKS(2) = {'RATE P -Kp*P*CL2'};
msx.TANKS(3) = {'RATE C_FRA -K_FAST*C_FRA*CL2'};
msx.TANKS(4) = {'RATE C_SRA -K_SLOW*C_SRA*CL2'};

% Section Sources
% <type> <nodeID> <specieID> <strength> (<patternID>)
msx.SOURCES(1) = {'SETPOINT R1 CL2 0.5 CL2PAT'};
msx.SOURCES(2) = {'SETPOINT R1 C_SRA 1.85 C_SRA_REPAT'}; % Monteiro
msx.SOURCES(3) = {'SETPOINT R2 CL2 0.5 CL2PAT'};
msx.SOURCES(4) = {'SETPOINT R2 C_SRA 1.85 C_SRA_REPAT'}; % Monteiro

% Section Patterns
% <patternID> <multiplier> <multiplier> 
msx.PATTERNS(1) = {'PathPAT 1'};
msx.PATTERNS(2) = {'CL2PAT 1'};
msx.PATTERNS(3) = {'C_FRAPAT 1'};
msx.PATTERNS(4) = {'C_SRAPAT 1'};
msx.PATTERNS(5) = {'C_SRA_REPAT 1'};

% Section Report
msx.REPORT(1) = {'NODES ALL'}; % Write results for all species at all nodes
msx.REPORT(2) = {'LINKS ALL'}; % Write results for all species at all links
 
% Load a network
for i=1
    % i=1;
    tmpinp = 'L-TOWN_stream_paper'; % The network as generated in the paper
    % tmpinp = 'L-TOWN_stream'; % Uncomment here to use the network (.inp file) that was created in Step 2
    inpname = ['networks\',tmpinp, num2str(i),'.inp'];
    d = epanet(inpname, 'loadfile');
    t_d = 3;
    % Define scenario parameters
    d.setTimeHydraulicStep(300);
    d.setTimeReportingStep(300);
    Tts= t_d*24*60*60/300;
    duration_hrs = t_d*24;
    duration_sec = duration_hrs*60*60;
    d.setTimeSimulationDuration(duration_sec)   

    % Run hydraulics
    hydraulics = d.getComputedTimeSeries;
    Demand= hydraulics.Demand;

    % Write MSX File
    % Section Parameters
    % PIPE <pipeID> <paramID> <value> % TANK <tankID> <paramID> <value>
    % Calculation of wall chlorine decay parameters. User can skip this step and use default values for A
    nl=double(d.getLinkCount);
    B=14;
    for m =1:nl
        indp = d.getLinkNameID{m};
        if d.getLinkRoughnessCoeff(m)>=140
            A=0.01;
        else
            A=rand(1);
        end
        msx.PARAMETERS(m) = {['PIPE ',indp,' A ',num2str(A)]};
        msx.PARAMETERS(nl+m) = {['PIPE ',indp,' B ',num2str(B)]};
    end

    % Write MSX
    d.writeMSXFile(msx);
    
    % Printing
    fprintf('MSX File written successfully.\n');
    
    % Load MSX File
    d.loadMSXFile(msx.FILENAME);

    % Set contaminant concentration and a time profile (start and end time)
    Campylobacter = 9.02e+06; % Initial concentration of Campylobacter
    % Campylobacter_HIGH = 6.2e+07; % Higher initial concentration of Campylobacter
    % Enterovirus = 1.39e+06; % Initial concentration of enterovirus
    % Enterovirus_HIGH = 2.08e+07; % Higher initial concentration of enterovirus
    % Cryptosporidium = 3.54e+07; % Initial concentration of Cryptosporidium
    % Cryptosporidium_HIGH = 5.4e+08; % Higher initial concentration of Cryptosporidium
    Pathogen_concentration = Campylobacter; % CFU/L or PFU/L or oocysts/L
    TOC_concentration=140; % 140 mg/L for normal scenario, 250 mg/L for high scenario
    C_FRA_fraction=0.4; % fraction of fast reducing agent 
    C_SRA_fraction=0.6; % fraction of slow reducing agent 
    Injection_rate=100; % liter/h
    C_FRA=C_FRA_fraction*TOC_concentration*Injection_rate;
    C_SRA=C_SRA_fraction*TOC_concentration*Injection_rate;
    Pathogen_mass = Pathogen_concentration*Injection_rate; % Injection_concentration*injection_rate = mass/hour
    Injection_start_time = 384;
    Injection_end = 480 ; % 1 hour=396, 2 hours= 408, 8 hours=480, 24 hours= 672
    scenario = 1;
    % SOURCES
    Source_nodes = d.getNodeNameID([112,775,44]);
    injection_conc_P = Pathogen_mass; 
    injection_conc_C_FRA = C_FRA; 
    injection_conc_C_SRA = C_SRA;
    Injection_stop_time = Injection_end; % 1 hour=12, 8 hours=96, 24 hours=288    
    % Initialize the contamination vector 
    Path_pat = zeros(1, Tts); % initialize the vector
    C_FRAPAT= zeros(1, Tts);
    C_SRAPAT= zeros(1, Tts);
    Path_pat(Injection_start_time:Injection_stop_time) = 1; % create the injection pattern 
    C_FRAPAT(Injection_start_time:Injection_stop_time) = 1; % create the injection pattern
    C_SRAPAT(Injection_start_time:Injection_stop_time) = 1; % create the injection pattern
    % Define injection node and run all scenarios
    for nn = 1:3
        disp(['Simulating contamination scenario ',num2str(nn),' of 3...'])
        injection_node_P= Source_nodes(nn);
        injection_node_C_FRA=Source_nodes(nn);
        injection_node_C_SRA=Source_nodes(nn);
        % Simulate the pathogen contamination
        % Specify Pathogen injection source
        d.setMSXSources(injection_node_P, 'P', 'MASS', injection_conc_P, 'PathPAT');
        d.setMSXSources(injection_node_C_FRA, 'C_FRA', 'MASS', injection_conc_C_FRA, 'C_FRAPAT');
        d.setMSXSources(injection_node_C_SRA, 'C_SRA', 'MASS', injection_conc_C_SRA, 'C_SRAPAT');
        % Set pattern of injection
        d.setMSXPattern('PathPAT',Path_pat);
        d.setMSXPattern('C_FRAPAT',C_FRAPAT);
        d.setMSXPattern('C_SRAPAT',C_SRAPAT);
        % Solve MSX quality dynamics
        Qmsx_species_P = d.getMSXComputedQualitySpecie('P');
        % Qmsx_species_C_FRA = d.getMSXComputedQualitySpecie('C_FRA');
        % Qmsx_species_C_SRA = d.getMSXComputedQualitySpecie('C_SRA');
        % Qmsx_species_CL2 = d.getMSXComputedQualitySpecie('CL2');
        d.setMSXSources(injection_node_P, 'P', 'MASS', 0, 'PathPAT'); % Reset injection source
        d.setMSXSources(injection_node_C_FRA, 'C_FRA', 'MASS', 0, 'C_FRAPAT'); % Reset injection source
        d.setMSXSources(injection_node_C_SRA, 'C_SRA', 'MASS', 0, 'C_SRAPAT'); % Reset injection source
        NodeQuality_strm_P{scenario} = Qmsx_species_P.NodeQuality;
        % NodeQuality_strm_CL2{scenario} = Qmsx_species_CL2.NodeQuality;
        % NodeQuality_strm_C_FRA{scenario} = Qmsx_species_C_FRA.NodeQuality;
        % NodeQuality_strm_C_SRA{scenario} = Qmsx_species_C_SRA.NodeQuality;
        scenario=scenario+1;
    end
    save(['./Campylobacter_8h_',num2str(i)],'NodeQuality_strm_P', 'msx','hydraulics','d')
end

%=================================================================================================================================%
%=================================================================================================================================%


%% Step 4. Calculate downstream population
% Calculate the downstream population for the three contamination locations (Loc_L, Loc-M, Loc-S).
% try 
% d.unload
% catch ERR
% end 
% fclose all;clear class;clear all;clc;close all;
% 
% %%% Load
% load stream_demands_paper.mat % Stream demands as generated in the paper
% % load stream_demands.mat % Uncomment here to use the demands mat file you just created in Step 1
% 
% tmpinp = 'L-TOWN_stream_paper'; % The network as generated in the paper
% % tmpinp = 'L-TOWN_stream'; % Uncomment here to use the network (.inp file) that was created in Step 2
% i=1;
% inpname = ['networks\',tmpinp, num2str(i),'.inp'];
% d = epanet(inpname, 'loadfile');
% hydraulics = d.getComputedTimeSeries;
% Dt= double(d.getTimeHydraulicStep)/3600;
% 
% %%% Initialize
% NodeNum= d.getNodeCount;
% NodeIDs=d.getNodeJunctionNameID;
% nj = double(d.getNodeJunctionCount);
% A=zeros(NodeNum);
% distance=zeros(nj);
% 
% %%% Create digraph
% avg_flows = mean(hydraulics.Flow);
% Flowsign= sign(avg_flows);
% Nidx= d.getLinkNodesIndex;
% for i=1:size(Nidx,1)
%     if Flowsign(i)==1
%         A(Nidx(i,1),Nidx(i,2)) = 1;
%     else
%         A(Nidx(i,2),Nidx(i,1)) = 1;
%     end
% end
% g = digraph(A);
% 
% %%% Find downstream nodes for each node
% for i=d.getNodeJunctionIndex
%     Downstream.(NodeIDs{i})= nearest(g, i, Inf);
% end
% Exclude= [783 784 785];
% for i=d.getNodeJunctionIndex
%     Downstream.(NodeIDs{i}) = Downstream.(NodeIDs{i})(find(ismember(Downstream.(NodeIDs{i}),Exclude)==0));
% end
% 
% Downstream_high=Downstream.n112;
% Downstream_mid=Downstream.n775;
% Downstream_low=Downstream.n44;
% b=1;
% Stream.Scenario{b}=Stream.Scenario{b}.*1000;
% Vpd= 150;
% V_sc=cumsum(Stream.Scenario{b}(385:672,:),1)*Dt; %Calculate volume consumed per node in 24 hours
% People_per_node= (V_sc(288,:)./Vpd); 
% People_per_node=round(People_per_node);
% Pop_high=sum(People_per_node(Downstream_high)');
% Pop_mid=sum(People_per_node(Downstream_mid)');
% Pop_low=sum(People_per_node(Downstream_low)');
% save('./Downstream_population.mat', "Pop_low","Pop_mid","Pop_high","Downstream_low","Downstream_mid","Downstream_high","People_per_node")

%=================================================================================================================================%
%=================================================================================================================================%


%% Step 5. Calculate infection risk
% ASSUMPTION: 1 liter per person in 1 day and the liter is consumped as 0.25 liter or less in one timestep 
try 
d.unload
catch ERR
end 
fclose all;clear class;clear all;clc;close all 
load stream_demands_paper.mat % Stream demands as generated in the paper
% load stream_demands.mat % Uncomment here to use the demands mat file you just created in Step 1

% Enterovirus

% load Enterovirus_8h_1.mat % Load contamination results for enterovirus
% load Enterovirus_8h_noCL21.mat % Load contamination results for enterovirus without CL2

% Campylobacter

load Campylobacter_8h_1.mat % Load contamination results for Campylobacter
% load Campylobacter_8h_noCL21.mat % Load contamination results for Campylobacter without CL2

% Cryptosporidium

% load Cryptosporidium_8h_1.mat % Load contamination results for Ctyptosporidium

%% Sensitivity analysis
% load Enterovirus_2h_1.mat % Load enterovirus contamination results for 2h contamination scenario
% load Enterovirus_24h_1.mat % Load enterovirus contamination results for 24h contamination scenario
% load Enterovirus_8h_HIGH_1.mat % Load enterovirus contamination results for high initial concentration contamination scenario
% load Enterovirus_8h_lowinact_1.mat % Load enterovirus contamination results for low inactivation rate contamination scenario

% load Campylobacter_2h_1.mat % Load Campylobacter contamination results for 2h contamination scenario
% load Campylobacter_24h_1.mat % Load Campylobacter contamination results for 24h contamination scenario
% load Campylobacter_8h_HIGH_1.mat % Load Campylobacter contamination results for high initial concentration contamination scenario
% load Campylobacter_8h_lowinact_1.mat % Load Campylobacter contamination results for low inactivation rate contamination scenario

load Downstream_population_paper.mat %  Load the downstream population as generated in the paper
% load Downstream_population.mat %  Uncomment here to use the downstream population mat file that was created in Step 4

Downstream_pop=[Pop_high Pop_mid Pop_low];

sensenum=1;
for b = 1:sensenum
    disp(['Sensitivity scenario ',num2str(b)])
    if sensenum==1
        inpname = 'L-TOWN_stream_paper1.inp';
    else
        inpname = ['L-TOWN_stream' num2str(b) '.inp'];
    end

    % load(matFileName);
    d = epanet(inpname, 'loadfile-ph');
    for i=1:3
        NodeQuality_strm_P{i}(:, end-2:end) = [];
    end
    Demand= hydraulics.Demand(:, d.getNodeJunctionIndex);
    Demand= Demand.*1000;  %CMH to L/h
    Stream.Scenario{b}=Stream.Scenario{b}.*1000;
    Stream_faucet.Scenario{b}=Stream_faucet.Scenario{b}.*1000; %CMH to L/h
    Stream_tap=Stream_faucet.Scenario{b}(385:672,:);
    Stream_tap=round(Stream_tap);
    t_d = 3;
    K= t_d*24*60*60/d.getTimeHydraulicStep;
    Dt= double(d.getTimeHydraulicStep)/3600;
    
    %Create consumed stream matrix
    junctionIndex = d.getNodeJunctionIndex;
    numJunctions = numel(junctionIndex);
    people_tot=[];
    stream_tot=[];
    fraction=[];
    Consumed_Stream=[];
    
    for i=1:numJunctions
        stream_tot(i)=sum(Stream_tap(:,i));
        people_tot(i)=People_per_node(i);
        fraction(i)= people_tot(i)/stream_tot(i);
        Consumed_Stream(:,i)=Stream_tap(:,i).*fraction(i);
        % liter=sum(Consumed_Stream(:,1))/(People_per_node(1)*3);
    end
    % Initialize Volume matrix to calculate water consumption distribution per person per node
    Volume = cell(1, numJunctions);
    for m = 1:numJunctions
        Volume{m} = zeros(size(Stream_faucet.Scenario{b}(385:672,:),1), People_per_node(m));
    end
 
    for m=1:numJunctions
    disp(['Calculating volume consumption for node ',num2str(m),' out of 782...'])
        for n=1:People_per_node(m)
            for i = 1:size(Consumed_Stream,1)
                while Consumed_Stream(i,m) > 0
                    if Consumed_Stream(i,m) < 0.00001
                        break
                    end
                    % If the current individual hasn't drunk 1 liter yet
                    if sum(Volume{m}(:, n)) ~= 1
                        if Consumed_Stream(i,m) >= 0.25
                            if 1 - sum(Volume{m}(:, n)) >= 0.25
                                Volume{m}(i, n) = 0.25;
                                Consumed_Stream(i,m) = Consumed_Stream(i,m) - Volume{m}(i, n);
                            else
                                Volume{m}(i, n) = 1 - sum(Volume{m}(:, n));
                                Consumed_Stream(i,m) = Consumed_Stream(i,m) - Volume{m}(i, n);
                            end
                        else
                            if Consumed_Stream(i,m) + sum(Volume{m}(:, n))<1
                                Volume{m}(i, n) = Consumed_Stream(i,m);
                                Consumed_Stream(i,m) = Consumed_Stream(i,m) - Volume{m}(i, n);
                            else
                                Volume{m}(i, n)= 1-sum(Volume{m}(:, n));
                                Consumed_Stream(i,m) = Consumed_Stream(i,m) - Volume{m}(i, n);
                            end
                        end
                    end
    
                    % Move to the next individual
                    n = n + 1;
                    if n > People_per_node(m)
                        n = 1; % That way we ensure that all volume has been consumed in every timestep
                    end
                end
            end
        end
    end
    
    Total_water_pp=cell(1, numJunctions);
    for m = 1:numJunctions
        Total_water_pp{m} = zeros(length(Stream_faucet.Scenario{1}), People_per_node(m));
    end
    for m=1:length(Total_water_pp)
        Total_water_pp{m}=sum(Volume{m});
    end
    %%% Calculate Dose
    % Initialize size of dose matrix
    % Dose=cell(1, numJunctions);
    for m = 1:numJunctions
        Dose{m} = zeros(size(Stream_tap,1), People_per_node(m));
    end
    
    % Initialize contamination matrix for 1 day
    Cont_matrix=[];
    for m=1:length(NodeQuality_strm_P)
        Cont_matrix{m}=NodeQuality_strm_P{m}(385:672,:);
    end
    
    % Calculate dose for each contamination scenario
    Scenario_Dose = cell(1, 3);
    
    for i=1:3
        for m=1:length(Volume)
            Dose{m}=Volume{m}.*Cont_matrix{i}(:,m);
        end
        Scenario_Dose{i}=Dose;
    end
    
    %%% Risk of infection (dose-response)
    % Initialize variables
    Risk_025 = cell(1, 3);
    
    % Campylobacter parameters
    alpha_campy = 0.38;
    beta_campy=0.51;
    
    % Cryptosporidium parameters
    alpha_crypto= 0.106;
    beta_crypto= 0.295;
    
    % Enterovirus parameters
    r_entero=0.014472;
    
    % Loop over the 3 cells in Dose
    for k = 1:3 % Loc-L=1, Loc-M=2, Loc-S=3
        % Initialize cell array to store risk matrices for this Dose cell
        Risk_025{k} = cell(1, 782);
    
        % Loop over the 782 cells in this Scenario_Dose cell
        for m = 1:782
        disp(['Calculating dose for node ',num2str(m),' out of 782',' in scenario',num2str(k)])
            % Extract the 288xindividuals matrix from the current cell
            dose_matrix = Scenario_Dose{k}{1, m};
    
            % Calculate Risk_025 for this matrix
            % Campylobacter Beta-Poisson with hypergeometric function model
            risk_matrix_campy = 1 - hypergeom(alpha_campy,alpha_campy+beta_campy,-dose_matrix);
    
            % Enterovirus and Cryptosporidium exponential model
            % risk_matrix_entero= 1-exp(-r_entero*dose_matrix);
            % risk_matrix_crypto= 1 - hypergeom(alpha_crypto,alpha_crypto+beta_crypto,-dose_matrix);
    
            % Store the risk matrix in the current cell
            Risk_025{k}{1, m} = risk_matrix_campy;
        end
    end
    
    % Number of infections per timestep
    % Initialize variables
    Total_Infections_ts = cell(1, 3);
    
    % Loop over the 3 cells in Dose
    for k = 1:3 % Loc-L=1, Loc-M=2, Loc-S=3
    disp(['Calculating number of infections and infection risk for scenario ',num2str(k),' out of 3...'])
        % Initialize cell array to store total infections for this Dose cell
        Total_Infections_ts{k} = [];
    
        % Loop over the 782 nodes
        for m = 1:782
            % Extract the risk matrix from the current cell
            risk_matrix = Risk_025{k}{1, m};
            risk_matrix_total{k}{m}=risk_matrix;
            if isempty(risk_matrix)
                continue;  % Skip the computations for this node
            else
                risk_matrix = risk_matrix';
                % calculate total risk per person for all timesteps
                for p=1:size(risk_matrix,1)
                    Total_risk_per_person{k}{m}(p)=1-prod(1-risk_matrix(p,:));
                end
                % calculate total risk per person for each timestep
                Total_risk_per_person_for_each_ts{k}{m}=zeros(size(risk_matrix));
                Total_risk_per_person_for_each_ts{k}{m}(:,1)=risk_matrix(:,1);
                for n=2:288
                    for p=1:size(risk_matrix,1)
                        Total_risk_per_person_for_each_ts{k}{m}(p,n)=1 - (1-Total_risk_per_person_for_each_ts{k}{m}(p,n-1))*(1-risk_matrix(p,n));
                    end
                end
                % calculate total infections per timestep
                for n=1:288
                    Total_infections_per_timestep{k}(m,n)=sum(Total_risk_per_person_for_each_ts{k}{m}(:,n),1);
                end
            end
        end
        Total_Infections_day{k}=sum(cell2mat(Total_risk_per_person{k})); % Calculate the total number of infections in a day
        Total_risk_of_infection{k}=(Total_Infections_day{k}/sum(People_per_node))*100; % Calculate the risk of infection expressed as the percentage of people being infected
        Affected_Risk{k}=(Total_Infections_day{k}/Downstream_pop(k))*100; % Calculate the risk of infection for the downstream affected population. Results might deviate slightly (e.g. 100.13 %) because of population round at each node.
        Total_infections_per_timestep_aggregated{k} = sum(Total_infections_per_timestep{k}, 1);
    end
end
save('./Infection_risk_Campylobacter_8h','Volume','Total_water_pp','Consumed_Stream','Dose','Risk_025','Total_Infections_ts','Total_risk_of_infection','Total_Infections_day','Affected_Risk','People_per_node','Total_infections_per_timestep_aggregated','Total_infections_per_timestep','Total_risk_per_person_for_each_ts')

%% Plotting Risk of infection through time
% Here you can plot the paper or your own results. The default is the paper results.

load Infection_risk_Enterovirus_8h.mat
% load Infection_risk_Enterovirus_8h_noCL2.mat
load Infection_risk_Campylobacter_8h.mat
% load Infection_risk_Campylobacter_8h_noCL2.mat
% load Infection_risk_Cryptosporidium_8h.mat

for k=1:3
    Total_infections_per_timestep_aggregated{k}=Total_infections_per_timestep_aggregated{k}./sum(People_per_node);
end

% Plot here
figure;
    for k=1:3
        x_axis = (1:numel(Total_infections_per_timestep_aggregated{k}))'* 5 / 60;  % convert to hours
        plot(x_axis, 100*Total_infections_per_timestep_aggregated{k}, 'LineWidth', 4); hold on;  % convert to percentage
    end
xlabel('Time (hours)', 'FontSize', 30); % Increased font size for x-axis
ylabel('Risk of infection (%)', 'FontSize', 30);     % Increased font size for y-axis, also convert to percentage
xline(8, '--r', 'LineWidth', 4);
% xline(16, '--r', 'LineWidth', 4);
h_legend = legend('Loc-L', 'Loc-M', 'Loc-S','Location', 'northwest');
set(h_legend, 'FontSize', 30);

% Increase the size of the axes
ax = gca;  % get current axes
ax.FontSize = 30; % Increase size of axis values
ylim([0 12]);
hold off;

%% Uncomment to plot the results of the Sensitivity analysis. You need first to have generated the different scenarios (e.g. contamination of 2 hours, contamination of 24 hours, etc.).

% Load the mat file for each sensitivity analysis parameter (2h, 24h, High concentration, low inactivation) and assign it to each pathogens category
% (e.g. Campy2{k}=Total_infections_per_timestep_aggregated{k};). Then divide each category with the total population and plot the results.

% Campylobacter
k=1; % Change to 2 or 3 to select contamination locations. Default is 1= Loc-L
% load Infection_risk_Campylobacter_2h.mat
% load Infection_risk_Campylobacter_8h.mat
% load Infection_risk_Campylobacter_24h.mat
% load Infection_risk_Campylobacter_8h_HIGH.mat
% load Infection_risk_Campylobacter_8h_lowinact.mat

% Campy2{k}=Total_infections_per_timestep_aggregated{k};
% Campy8{k}=Total_infections_per_timestep_aggregated{k};
% Campy24{k}=Total_infections_per_timestep_aggregated{k};
% CampyHigh{k}=Total_infections_per_timestep_aggregated{k};
% Campylowinact{k}=Total_infections_per_timestep_aggregated{k};

% Campy2{k}=Campy2{k}./sum(People_per_node);
% Campy8{k}=Campy8{k}./sum(People_per_node);
% Campy24{k}=Campy24{k}./sum(People_per_node);
% CampyHigh{k}=CampyHigh{k}./sum(People_per_node);
% Campylowinact{k}=Campylowinact{k}./sum(People_per_node);

% Enterovirus

% load Infection_risk_Enterovirus_2h.mat
% load Infection_risk_Enterovirus_8h.mat
% load Infection_risk_Enterovirus_24h.mat
% load Infection_risk_Enterovirus_8h_HIGH.mat
% load Infection_risk_Enterovirus_8h_lowinact.mat
% 
% k=1; % Change to 2 or 3 to select contamination locations. Default is 1= Loc-L
% Entero2{k}=Total_infections_per_timestep_aggregated{k};
% Entero8{k}=Total_infections_per_timestep_aggregated{k};
% Entero24{k}=Total_infections_per_timestep_aggregated{k};
% EnteroHigh{k}=Total_infections_per_timestep_aggregated{k};
% Enterolowinact{k}=Total_infections_per_timestep_aggregated{k};

% Entero2{k}=Entero2{k}./sum(People_per_node);
% Entero8{k}=Entero8{k}./sum(People_per_node);
% Entero24{k}=Entero24{k}./sum(People_per_node);
% EnteroHigh{k}=EnteroHigh{k}./sum(People_per_node);
% Enterolowinact{k}=Enterolowinact{k}./sum(People_per_node);

% Plot here
% figure;
% x_axis = (1:numel(Total_infections_per_timestep_aggregated{k}))'* 5 / 60;  % convert to hours
% plot(x_axis, 100*Entero2{k}, ':m', 'LineWidth', 4); hold on;  % convert to percentage, red dashed line
% plot(x_axis, 100*Entero24{k}, '--m', 'LineWidth', 4); hold on;  % convert to percentage, blue dotted line
% plot(x_axis, 100*EnteroHigh{k}, '-.g', 'LineWidth', 4); hold on;  % convert to percentage, green dash-dot line
% plot(x_axis, 100*Enterolowinact{k}, '-c', 'LineWidth', 4); hold on;  % convert to percentage, cyan solid line
% plot(x_axis, 100*Entero8{k}, '-m', 'LineWidth', 4); hold on;  % convert to percentage, magenta solid line

% xlabel('Time (hours)', 'FontSize', 30); % Increased font size for x-axis
% ylabel('Risk of infection (%)', 'FontSize', 30);     % Increased font size for y-axis, also convert to percentage
% xline(8, '-r', 'LineWidth', 4);
% xline(24, '--r', 'LineWidth', 4);

% define labels for the legend
% legend_labels = {'2 hours', '24 hours', 'High concentration', 'Low inactivation','Original scenario'};
% h_legend = legend(legend_labels, 'Location', 'northwest');
% set(h_legend, 'FontSize', 26);

% Increase the size of the axes
% ax = gca;  % get current axes
% ax.FontSize = 30; % Increase size of axis values
% ylim([0 30]);
% hold off;







