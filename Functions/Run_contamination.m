%% Step 3. Run contamination
try 
d.unload
catch ERR
end 
fclose all;clear class;clear all;clc;close all;
% Start EPANET MATLAB TOOLKIT
start_toolkit; 

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
    fprintf('The MSX File written successfully.\n');
    
    % Load MSX File
    d.loadMSXFile(msx.FILENAME);

    % Set contaminant concentration and a time profile (start and end time)
    Campylobacter = 9.02e+06; % Initial concentration of Campylobacter
    Campylobacter_HIGH = 6.2e+07; % Higher initial concentration of Campylobacter
    Enterovirus = 1.39e+06; % Initial concentration of enterovirus
    Enterovirus_HIGH = 2.08e+07; % Higher initial concentration of enterovirus
    Cryptosporidium = 3.54e+07; % Initial concentration of Cryptosporidium
    Cryptosporidium_HIGH = 5.4e+08; % Higher initial concentration of Cryptosporidium
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
        scenario
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
        Qmsx_species_C_FRA = d.getMSXComputedQualitySpecie('C_FRA');
        Qmsx_species_C_SRA = d.getMSXComputedQualitySpecie('C_SRA');
        Qmsx_species_CL2 = d.getMSXComputedQualitySpecie('CL2');
        d.setMSXSources(injection_node_P, 'P', 'MASS', 0, 'PathPAT'); % Reset injection source
        d.setMSXSources(injection_node_C_FRA, 'C_FRA', 'MASS', 0, 'C_FRAPAT'); % Reset injection source
        d.setMSXSources(injection_node_C_SRA, 'C_SRA', 'MASS', 0, 'C_SRAPAT'); % Reset injection source
        NodeQuality_strm_P{scenario} = Qmsx_species_P.NodeQuality;
        NodeQuality_strm_CL2{scenario} = Qmsx_species_CL2.NodeQuality;
        NodeQuality_strm_C_FRA{scenario} = Qmsx_species_C_FRA.NodeQuality;
        NodeQuality_strm_C_SRA{scenario} = Qmsx_species_C_SRA.NodeQuality;
        scenario=scenario+1;
    end
    save(['./Campylobacter_8h_',num2str(i)],'NodeQuality_strm_P','NodeQuality_strm_CL2','NodeQuality_strm_C_FRA','NodeQuality_strm_C_SRA','msx','hydraulics','d')
end
