%% Step 2. Assign new network demands
try 
d.unload
catch ERR
end 
fclose all;clear class;
close all;clear all;clc;
% Start EPANET MATLAB TOOLKIT
start_toolkit;
addpath(genpath(pwd))

%%% Load Network:
inpname = 'L-TOWN.inp';
dispname = 'L-TOWN';
d=epanet(inpname);
nj=d.getNodeJunctionCount;
nn=d.getNodeCount;

% Assign demands to network:
%%% Get existing patterns:
%%%     d.getNodeDemandPatternIndex{categoryIndex}(nodeIndex)
demInd1 = double(d.getNodeDemandPatternIndex{1}(:)); 
demInd1G = demInd1;
demInd2 = double(d.getNodeDemandPatternIndex{2}(:));
demInd3 = double(d.getNodeDemandPatternIndex{3}(:));
baseDemInd1 = double(d.getNodeBaseDemands{1}(:));
baseDemInd2 = double(d.getNodeBaseDemands{2}(:));
baseDemInd3 = double(d.getNodeBaseDemands{3}(:));

%%% Zero base demands
for i=1:nn
    disp(['Zero base demand ',num2str(i)])
    d.setNodeBaseDemands(i, 1, 0)
    d.deleteNodeJunctionDemand(i,2)
    d.deleteNodeJunctionDemand(i,2)
end

%%% Delete three(3) old patterns:
d.deletePattern(1)
d.deletePattern(1) % pattern number 2 become 1 after deletion of 1
d.deletePattern(1)

%%% Save new input file:
emptyInpName = ['networks\',dispname,'_empty','.inp'];
d.saveInputFile(emptyInpName);
disp('EMPTY NETWORK READY!')

%%% Load new actual demands:
homeFolder = pwd;
addpath([homeFolder '/Paper_results']);
load stream_demands_paper.mat % Stream demands as generated in the paper
% load stream_demands.mat % Uncomment here to use the demands mat file you just created in Step 1

%%% Calculate and assign new base demands:
% for 100 scenarios change the following for loop (1:100)
for m=1
    d=epanet(emptyInpName);
    demInd1 = demInd1G;
    base_Stream_demand = mean(Stream.Scenario{m});
    for i=1:nj
        disp(['Assign base demand ',num2str(i)])
        d.setNodeBaseDemands(i, 1, base_Stream_demand(i))
    end

%%% Calculate new patterns:
    pattern_Stream_demand = Stream.Scenario{m}/base_Stream_demand;
    pattern_Stream_demand(isnan(pattern_Stream_demand))=0;
    
    %%% Add new patterns:
    for i = 1:nj
        demands = pattern_Stream_demand(:);
        resDemPatInd(i)=d.addPattern(['P-Res-',num2str(i)],demands);
        disp(['Creating pattern Residential ',num2str(i)])
    end
    
    for i=1:nj
        disp(['Indexing pattern ',num2str(i),' category 1'])
        if demInd1(i)==0
            continue 
        elseif demInd1(i)==1 % Residential
            demInd1(i)=i;
        else
            error('unknown demand pattern')
        end
    end
    
    %%% Assign new patterns:
    for i=1:nn
        disp(['Assigning pattern ',num2str(i),' out of ',num2str(nn)])
        d.setNodeDemandPatternIndex(i, 1, demInd1(i))
    end
    
    % Correct times:
    d.setTimeReportingStep(300)
    d.setTimeHydraulicStep(300)
    d.setTimePatternStep(300)

% Save new input file:
    newInpname = ['networks\',dispname,'_stream',num2str(m),'.inp'];
    d.saveInputFile(newInpname);
    disp('NETWORK READY!')
end
