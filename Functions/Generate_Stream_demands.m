%% Step 1. Generate Stream demands
try 
d.unload
catch ERR
end 
fclose all;clear class;clear all;clc;close all;
% Start EPANET MATLAB TOOLKIT
start_toolkit;
% Load a network.
inpname = 'L-TOWN.inp';

d = epanet(inpname);
d.setTimeReportingStep(300)
d.setTimeHydraulicStep(300)
d.setTimePatternStep(300)
hydraulics = d.getComputedTimeSeries;

%%% Calculate hydraulics
Demand= hydraulics.Demand(:,d.getNodeJunctionIndex);
Demand= Demand.*1000;  %CMH to L/h
K= round(size(Demand, 1)/7);
Dt= double(d.getTimeHydraulicStep)/3600;
Vpd= 150;
V_sc=cumsum(Demand(1:K,:),1)*Dt; %Calculate volume consumed per node in 24 hours
People_per_node= V_sc(end,:)/Vpd; % People per node 

%%% ::: STREAM :::

% AVAILABLE FIXTURES
% 1. toilet
% 2. shower
% 3. faucet
% 4. clothes washer
% 5. dishwasher
% 6. bathtub
%%% ::: INPUT SETTINGS :::
%%% ::: LOADING COMPLETE DATABASE :::
homeFolder = pwd;
addpath([homeFolder '/_DATA']); % Path to the folder where the database.mat file is stored
load database.mat

People_per_node_rnd= round(People_per_node);
Population_unc=0.1;
Node=d.getNodeNameID;
Scenario={};
% for 100 scenarios change the following for loop (1:100)
for scen=1
    scen
    for i=1:length(People_per_node_rnd)
        disp(['Simulating Node ',num2str(i),' of ',num2str(length(People_per_node_rnd))])
        % Create +-10% uncertainty
        Population=People_per_node_rnd(i);
        Population_l=Population-Population_unc*Population;
        Population_u=Population+Population_unc*Population;
        Population=Population_l+rand(1,length(Population)).*(Population_u-Population_l);
        Population = round(Population);
    
    %input population
    home=0;
    %initialization
    StToilet=0;
    StShower=0;
    StFaucet=0;
    StClothesWasher=0;
    StDishwasher=0;
    TOTAL=0;
    
    while Population>0
    
    % --- A. Household size setting
    home=home+1;
    param.HHsize = randi(5,1); % This parameter should be in the interval (1,6).
    % From 1 to 5, it indicates the number of people living in the house. 6 means ">5".
    
    Population=Population-param.HHsize;
    % --- B. Water consuming fixtures selection
    % Legend:
    % 0 = not present
    % 1 = present
    
    param.appliances.StToilet = 1;
    param.appliances.HEToilet = 0;
    
    param.appliances.StShower = 1;
    param.appliances.HEShower = 0;
    
    param.appliances.StFaucet = 1;
    param.appliances.HEFaucet = 0;
    
    param.appliances.StClothesWasher = 1;
    param.appliances.HEClothesWasher = 0;
    
    param.appliances.StDishwasher = 1;
    param.appliances.HEDishwasher = 0;
    
    param.appliances.StBathtub = 1;
    param.appliances.HEBathtub = 0;
    
    % --- C. Time horizon length setting
    param.H = 3; % It is measured in [days]
    
    % --- D. Time sampling resolution
    param.ts = 30; % It is measured in [10 seconds] units. The maximum resolution allowed is 10 seconds (param.ts = 1).
    
    % Setting the seed
    % rng(1);
    
    % Parameters structure settings and check
    % Checking input consistency
    temp=checkInput(param);
    % clearvars -except param
    
    %%% ::: WATER END-USE TIME SERIES GENERATION :::
    
    % Initialization
    outputTrajectory = initializeTrajectories(param);
    % Include the first step
    appNames = fieldnames(outputTrajectory);
    for app=appNames'
        outputTrajectory.(char(app))=zeros(1,length(outputTrajectory.TOTAL)+30);
    end
    
    % End-use water use time series generation
    outputTrajectory = generateConsumptionEvents(outputTrajectory, param, database);
    % disp('End-use consumption trajectories created');
    
    % Total water use time series aggregation
    outputTrajectory = sumToTotal(outputTrajectory);
    % disp('Total consumption trajectory created');
    
    % Data scaling to desired sampling resolution
    outputTrajectory = aggregateSamplingResolution(outputTrajectory, param);
    % disp('Data scaled to desired sampling resolution');
    
    StToilet=outputTrajectory.StToilet+StToilet;
    StShower=outputTrajectory.StShower+StShower;
    StFaucet=outputTrajectory.StFaucet+StFaucet;
    StClothesWasher=outputTrajectory.StClothesWasher+StClothesWasher;
    StDishwasher=outputTrajectory.StDishwasher+StDishwasher;
    TOTAL=outputTrajectory.TOTAL+TOTAL;
    clear outputTrajectory;
    end
    
    output.output(scen).(Node{i}).StToilet=StToilet;
    output.output(scen).(Node{i}).StShower=StShower;
    output.output(scen).(Node{i}).StFaucet=StFaucet;
    output.output(scen).(Node{i}).StClothesWasher=StClothesWasher;
    output.output(scen).(Node{i}).StDishwasher=StDishwasher;
    output.output(scen).(Node{i}).TOTAL=TOTAL;
    
    end

    %%% Assign stream demand for total water use
    for j=1:length(People_per_node_rnd)
    Stream_demand_tot(:,j)= output.output(scen).(Node{j}).TOTAL;
    end
    Stream_demand_tot= (Stream_demand_tot.*12)/1000; %convert from L/5min to CMH
    Stream.Scenario{scen}= Stream_demand_tot;

    % Assign stream demand for injestion of water
    for j=1:length(People_per_node_rnd)
    Stream_demand_Faucet(:,j)= output.output(scen).(Node{j}).StFaucet;
    end
    Stream_demand_Faucet= (Stream_demand_Faucet.*12)/1000; %convert from L/5min to CMH
    Stream_faucet.Scenario{scen}= Stream_demand_Faucet;
end

% Saving
save ('./Stream_demands.mat','output','Stream','Stream_faucet')