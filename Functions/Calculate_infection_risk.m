%% Step 5. Calculate infection risk
% ASSUMPTION: 1 liter per person in 1 day and the liter is consumped as 0.25 liter or less in one timestep 
try 
d.unload
catch ERR
end 
fclose all;clear class;clear all;clc;close all
% Start EPANET MATLAB TOOLKIT
start_toolkit;
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

for b = 1
    b
    inpname = ['L-TOWN_stream' num2str(b) '.inp'];

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
    % Initialize Volume
    Volume = cell(1, numJunctions);
    for m = 1:numJunctions
        Volume{m} = zeros(size(Stream_faucet.Scenario{b}(385:672,:),1), People_per_node(m));
    end
 
    for m=1:numJunctions
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
    for k = 1:3
        % Initialize cell array to store risk matrices for this Dose cell
        Risk_025{k} = cell(1, 782);
    
        % Loop over the 782 cells in this Scenario_Dose cell
        for m = 1:782
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
    for k = 1:3 % Loc-L, Loc-M, Loc-S
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
