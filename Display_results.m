
%% Plotting Risk of infection through time
% Here you can plot the paper or your own results. The default are the paper results that first need to be generated.

load Infection_risk_Enterovirus_8h.mat
% load Infection_risk_Enterovirus_8h_noCL2.mat
% load Infection_risk_Campylobacter_8h.mat
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




%% Sensitivity analysis

% Load the mat file for each sensitivity analysis parameter (2h, 24h, High
% concentration, low inactivation) and assign it in each pathogens category
% (e.g. Campy2{k}=Total_infections_per_timestep_aggregated{k};). Then
% devide each category with the total population and plot the results.

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

Campy2{k}=Campy2{k}./sum(People_per_node);
Campy8{k}=Campy8{k}./sum(People_per_node);
Campy24{k}=Campy24{k}./sum(People_per_node);
CampyHigh{k}=CampyHigh{k}./sum(People_per_node);
Campylowinact{k}=Campylowinact{k}./sum(People_per_node);

% Enterovirus

load Infection_risk_Enterovirus_2h.mat
load Infection_risk_Enterovirus_8h.mat
load Infection_risk_Enterovirus_24h.mat
load Infection_risk_Enterovirus_8h_HIGH.mat
load Infection_risk_Enterovirus_8h_lowinact.mat
% 
k=1; % Change to 2 or 3 to select contamination locations. Default is 1= Loc-L
Entero2{k}=Total_infections_per_timestep_aggregated{k};
Entero8{k}=Total_infections_per_timestep_aggregated{k};
Entero24{k}=Total_infections_per_timestep_aggregated{k};
EnteroHigh{k}=Total_infections_per_timestep_aggregated{k};
Enterolowinact{k}=Total_infections_per_timestep_aggregated{k};

Entero2{k}=Entero2{k}./sum(People_per_node);
Entero8{k}=Entero8{k}./sum(People_per_node);
Entero24{k}=Entero24{k}./sum(People_per_node);
EnteroHigh{k}=EnteroHigh{k}./sum(People_per_node);
Enterolowinact{k}=Enterolowinact{k}./sum(People_per_node);

% Plot here
figure;
x_axis = (1:numel(Total_infections_per_timestep_aggregated{k}))'* 5 / 60;  % convert to hours
plot(x_axis, 100*Entero2{k}, ':m', 'LineWidth', 4); hold on;  % convert to percentage, red dashed line
plot(x_axis, 100*Entero24{k}, '--m', 'LineWidth', 4); hold on;  % convert to percentage, blue dotted line
plot(x_axis, 100*EnteroHigh{k}, '-.g', 'LineWidth', 4); hold on;  % convert to percentage, green dash-dot line
plot(x_axis, 100*Enterolowinact{k}, '-c', 'LineWidth', 4); hold on;  % convert to percentage, cyan solid line
plot(x_axis, 100*Entero8{k}, '-m', 'LineWidth', 4); hold on;  % convert to percentage, magenta solid line

xlabel('Time (hours)', 'FontSize', 30); % Increased font size for x-axis
ylabel('Risk of infection (%)', 'FontSize', 30);     % Increased font size for y-axis, also convert to percentage
xline(8, '-r', 'LineWidth', 4);
xline(24, '--r', 'LineWidth', 4);

% define labels for the legend
legend_labels = {'2 hours', '24 hours', 'High concentration', 'Low inactivation','Original scenario'};
h_legend = legend(legend_labels, 'Location', 'northwest');
set(h_legend, 'FontSize', 26);

% Increase the size of the axes
ax = gca;  % get current axes
ax.FontSize = 30; % Increase size of axis values
ylim([0 30]);
hold off;
