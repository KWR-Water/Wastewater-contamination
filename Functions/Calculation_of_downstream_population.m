%% Step 4. Calculate downstream population
% Calculate the downstream population for the three contamination locations (Loc_L, Loc-M, Loc-S).
try 
d.unload
catch ERR
end 
fclose all;clear class;clear all;clc;close all;
start_toolkit;

%%% Load
homeFolder = pwd;
addpath([homeFolder '/Paper_results']);
load stream_demands_paper.mat % Stream demands as generated in the paper
% load stream_demands.mat % Uncomment here to use the demands mat file you just created in Step 1

tmpinp = 'L-TOWN_stream_paper'; % The network as generated in the paper
% tmpinp = 'L-TOWN_stream'; % Uncomment here to use the network (.inp file) that was created in Step 2
i=1;
inpname = ['networks\',tmpinp, num2str(i),'.inp'];
d = epanet(inpname, 'loadfile');
hydraulics = d.getComputedTimeSeries;
Dt= double(d.getTimeHydraulicStep)/3600;

%%% Initialize
NodeNum= d.getNodeCount;
NodeIDs=d.getNodeJunctionNameID;
nj = double(d.getNodeJunctionCount);
A=zeros(NodeNum);
distance=zeros(nj);

%%% Create digraph
avg_flows = mean(hydraulics.Flow);
Flowsign= sign(avg_flows);
Nidx= d.getLinkNodesIndex;
for i=1:size(Nidx,1)
    if Flowsign(i)==1
        A(Nidx(i,1),Nidx(i,2)) = 1;
    else
        A(Nidx(i,2),Nidx(i,1)) = 1;
    end
end
g = digraph(A);

%%% Find downstream nodes for each node
for i=d.getNodeJunctionIndex
    Downstream.(NodeIDs{i})= nearest(g, i, Inf);
end
Exclude= [783 784 785];
for i=d.getNodeJunctionIndex
    Downstream.(NodeIDs{i}) = Downstream.(NodeIDs{i})(find(ismember(Downstream.(NodeIDs{i}),Exclude)==0));
end

Downstream_high=Downstream.n112;
Downstream_mid=Downstream.n775;
Downstream_low=Downstream.n44;
b=1;
Stream.Scenario{b}=Stream.Scenario{b}.*1000;
Vpd= 150;
V_sc=cumsum(Stream.Scenario{b}(385:672,:),1)*Dt; %Calculate volume consumed per node in 24 hours
People_per_node= (V_sc(288,:)./Vpd); 
People_per_node=round(People_per_node);
Pop_high=sum(People_per_node(Downstream_high)');
Pop_mid=sum(People_per_node(Downstream_mid)');
Pop_low=sum(People_per_node(Downstream_low)');
save('./Downstream_population.mat', "Pop_low","Pop_mid","Pop_high","Downstream_low","Downstream_mid","Downstream_high","People_per_node")