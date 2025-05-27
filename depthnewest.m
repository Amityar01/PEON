%% depthnewest.m
% Analyzes and plots the laminar and areal distribution of identified PEONs
% (Probability Encoding Omission Neurons) and all recorded neurons.
% This script generates figures corresponding to Figure 2E (laminar distribution)
% and Figure 2F (areal distribution) of the manuscript (Yaron et al., 2025).
%
% Workflow:
% 1. Defines cortical layer boundaries (supragranular, granular, infragranular)
%    based on depth measurements.
% 2. Takes neuron depth data (e.g., from 'T.depth') and PEON indices
%    (e.g., 'PEONs_training' from 'FindPEONS.m') as input.
% 3. Filters neurons by depth criteria (e.g., 0-1800 µm).
% 4. Counts the total number of selected neurons and the number of PEONs
%    within each defined cortical layer and within predefined depth intervals (e.g., 100µm bins).
% 5. Generates a horizontal bar plot (histogram) showing the laminar distribution
%    of all selected neurons and PEONs. Annotates with PEON percentages per layer.
% 6. Takes neuron areal location data (e.g., from 'T.loc') to identify neurons
%    in Primary Auditory Cortex (A1), Ventral Auditory Field (VAF), and
%    Anterior Auditory Field (AAF).
% 7. Counts total selected neurons and PEONs in each auditory area.
% 8. Generates pie charts illustrating the proportion of PEONs within A1, VAF, and AAF.
%
% Inputs (expected to be in the MATLAB workspace):
%   - T: A table or structure containing neuron metadata, specifically:
%        - T.depth: A vector with the cortical depth for each neuron/channel.
%        - T.loc: A cell array or categorical vector with the auditory area
%                 location (e.g., 'A1', 'VAF', 'AAF') for each neuron/channel.
%   - omind: Vector of indices for all neurons initially selected for analysis
%            (e.g., satisfying basic depth and activity criteria). This is used
%            to filter 'T.depth' and 'T.loc'.
%   - PEONs_training: Vector of indices (relative to the original neuron list)
%                     for neurons classified as PEONs by 'FindPEONS.m'.
%
% Outputs:
%   - Figure 1: Horizontal bar plot showing laminar distribution of all neurons
%               and PEONs, saved as PDF.
%   - Figure 2: A figure with three subplots, each a pie chart showing the
%               proportion of PEONs in A1, VAF, and AAF, saved as PDF.
%
% Dependencies:
%   - Standard MATLAB plotting functions.
%
close all;

% Define layer boundaries to match figure legend
supragranular_upper = 600;  % Supragranular: 0-600μm
granular_lower = 600;       % Granular: 600-900μm
granular_upper = 900;
infragranular_lower = 900;  % Infragranular: 900-1400μm
infragranular_upper = 1400; % Upper boundary of infragranular layer
lo=T.loc(1:8:end)';
deep=T.depth(1:8:end)';
dp=deep(omind);
% Prepare datasets
dataset1 = dp;

dataset3 = dp(PEONs_training);


% Create consistent 100μm bins for a clean look
depth_intervals = 0:100:1800; 

% Count neurons within each interval
counts1 = histcounts(dataset1, depth_intervals);
counts3 = histcounts(dataset3, depth_intervals);

% Calculate midpoints for plotting
depth_midpoints = depth_intervals(1:end-1) + diff(depth_intervals) / 2;

% Create figure with appropriate size ratio
figure('Position', [100, 100, 500, 600], 'Color', 'w');

% Plot histograms with the clean style
barh(depth_midpoints, counts1, 'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', 'none');
hold on;
barh(depth_midpoints, counts3, 'FaceColor', [0.3, 0.45, 0.7], 'EdgeColor', 'none');

% Customize axes
set(gca, 'YDir', 'reverse'); % Ensures 0 is at the top of vertical axis
set(gca, 'TickDir', 'out');  % Places markers outside axis as requested
set(gca, 'FontSize', 11);
set(gca, 'Box', 'off');
set(gca, 'Layer', 'top');

% Set axis limits
ax_lim = 250;  % Fixed at 250 like reference image
ylim([0 1800]);
xlim([0 ax_lim]);

% Set yticks to match the correct layer boundaries
set(gca, 'YTick', [0, 600, 900, 1400, 1800]);  % Include all key boundaries

% Add horizontal dotted lines at all layer boundaries
plot([0 ax_lim], [600, 600], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);  % Supragranular/Granular
plot([0 ax_lim], [900, 900], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);  % Granular/Infragranular
plot([0 ax_lim], [1400, 1400], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);  % Infragranular upper boundary

% Add layer labels at the middle of each layer
text(ax_lim*0.85, 300, 'I-III', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', 11);
text(ax_lim*0.85, 750, 'IV', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', 11);
text(ax_lim*0.85, 1150, 'V/VI', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', 11);

% Clean labels
ylabel('Depth (\mum)', 'FontSize', 11);
xlabel('Number of neurons', 'FontSize', 11);

% Calculate counts by layer using the EXACT layer boundaries from the legend
count_supra1 = sum(dataset1 < supragranular_upper);
count_gran1 = sum(dataset1 >= granular_lower & dataset1 <= granular_upper);
count_infra1 = sum(dataset1 > infragranular_lower & dataset1 <= infragranular_upper);

count_supra3 = sum(dataset3 < supragranular_upper);
count_gran3 = sum(dataset3 >= granular_lower & dataset3 <= granular_upper);
count_infra3 = sum(dataset3 > infragranular_lower & dataset3 <= infragranular_upper);

% Calculate percentages
percent_supra3 = count_supra3 / count_supra1 * 100;
percent_gran3 = count_gran3 / count_gran1 * 100;
percent_infra3 = count_infra3 / count_infra1 * 100;

% Create a simple legend
legend('All neurons', 'PEON', 'Location', 'northeast', 'FontSize', 10, 'Box', 'off');

% Add percentage text matching reference style
text(ax_lim*0.82, 1080, {'PEON:', ...
                      ['Supragranular: ' num2str(percent_supra3,'%.1f') '%'], ...
                      ['Granular: ' num2str(percent_gran3,'%.1f') '%'], ...
                      ['Infragranular: ' num2str(percent_infra3,'%.1f') '%']}, ...
     'FontSize', 10, 'VerticalAlignment', 'top');

hold off;

% Use painters renderer and export the figure as a PDF with vector content
set(gcf, 'Renderer', 'painters');
exportgraphics(gcf, 'neuron_distribution_by_depthn.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%%
ALL1=sum(lo(omind)=='1');
PEON1=sum(lo(PEONs_training)=='1');
ALL2=sum(lo(omind)=='2');
PEON2=sum(lo(PEONs_training)=='2');
ALL3=sum(lo(omind)=='3');
PEON3=sum(lo(PEONs_training)=='3');
% %%
% Actual counts of omission-sensitive neurons (ORNs and PEONs)

actual_counts_peons = [PEON1, PEON2, PEON3];

% Empirical distribution of all neurons across layers
empirical_distribution = [ALL1, ALL2, ALL3];
empirical_probabilities = empirical_distribution / sum(empirical_distribution);
% Define neuron types
neuronType = {'All Neurons'; 'PEON'};

% Counts for each area
Area1_Count = [ALL1; PEON1];
Area2_Count = [ALL2;  PEON2];
Area3_Count = [ALL3;  PEON3];

% Percentages for each area (normalized relative to "All Neurons")
Area1_Percent = [100;  (PEON1 / ALL1) * 100];
Area2_Percent = [100;  (PEON2 / ALL2) * 100];
Area3_Percent = [100;  (PEON3 / ALL3) * 100];

% Create the table
Tb = table(neuronType, Area1_Count, Area1_Percent, ...
            Area2_Count, Area2_Percent, ...
            Area3_Count, Area3_Percent, ...
            'VariableNames', {'NeuronType', 'Area1_Count', 'Area1_Percent', ...
            'Area2_Count', 'Area2_Percent', 'Area3_Count', 'Area3_Percent'});

% Display table
disp(Tb);

% Number of bootstrap samples
n_samples = 300000;

% Function to perform bootstrap sampling
bootstrap_samples_peons = mnrnd(sum(actual_counts_peons), empirical_probabilities, n_samples);

% Function to calculate p-values
calculate_p_values = @(actual_counts, bootstrap_samples) ...
    arrayfun(@(i) mean(bootstrap_samples(:, i) >= actual_counts(i)), 1:length(actual_counts));

% Calculate p-values for ORNs and PEONs
p_values_peons = calculate_p_values(actual_counts_peons, bootstrap_samples_peons);

% Output the results
disp('P-values for PEONs:');
disp(p_values_peons);

% Optional: Plot the bootstrap distributions for visualization


for i = 1:length(actual_counts_peons)
    subplot(3, 2, i*2);
    histogram(bootstrap_samples_peons(:, i), 50, 'Normalization', 'probability');
    hold on;
    xline(actual_counts_peons(i), 'r--', 'LineWidth', 2);
    title(sprintf('Bootstrap Distribution for PEONs - Layer %d', i));
    xlabel('Count');
    ylabel('Frequency');
end
%%
% --- Data Definition ---
% Define the names of the auditory fields
areaNames = {'A1', 'VAF', 'AAF'};
% Define the total number of neurons recorded in each area
totalNeurons = [632, 241, 117];
% Define the number of PEONs identified in each area
peonCounts = [93, 23, 7];
% Define the p-value for significance (for annotation)
p_value_A1 = 0.0036; % From your bootstrap analysis

% --- Calculations ---
% Calculate the percentage of PEONs in each area
peonPercentages = (peonCounts ./ totalNeurons) * 100;

% --- Plotting ---
% Create a new figure window
figure;

% Create the bar chart
hBar = bar(peonPercentages);

% Get the handle to the axes
ax = gca;

% Set the x-axis tick labels to the area names
ax.XTickLabel = areaNames;

% Set the y-axis label
ylabel('Percentage of PEONs (%)');

% Set the title for the plot (adjust 'F)' as needed for figure paneling)
title('F) Distribution Across Auditory Fields');

% Adjust y-axis limits for better visualization (optional)
ylim([0, max(peonPercentages) + 5]); % Add some padding above the highest bar

% --- Annotation for Significance ---
% Add an asterisk above the A1 bar if its p-value is significant
if p_value_A1 < 0.05
    % Get the x-coordinate for the center of the first bar (A1)
    x_coord_A1 = hBar.XData(1);
    % Get the y-coordinate for the top of the first bar (A1)
    y_coord_A1 = hBar.YData(1);
    % Place the asterisk slightly above the bar
    text(x_coord_A1, y_coord_A1 + 1, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
    % You could also add the p-value text if desired:
    % text(x_coord_A1, y_coord_A1 + 2, sprintf('p=%.4f', p_value_A1), 'HorizontalAlignment', 'center');
end

% Optional: Improve appearance
box off; % Turn off the plot box outline
set(ax, 'TickDir', 'out'); % Set tick direction outwards

%%
% --- Data Definition ---
% Assume the following variables have been previously calculated and exist in the workspace:
% ALL1, ALL2, ALL3: Total neuron counts for Area 1, 2, 3
% PEON1, PEON2, PEON3: PEON counts for Area 1, 2, 3
% p_values_peons = [p_val_area1, p_val_area2, p_val_area3]; % Calculated p-values

% Define the names of the auditory fields
areaNamesShort = {'A1', 'VAF', 'AAF'};

% Define the total number of neurons recorded in each area using existing variables
totalNeurons = [ALL1, ALL2, ALL3];

% Define the number of PEONs identified in each area using existing variable
peonCounts = [PEON1, PEON2, PEON3];

% Define the p-value for significance for each area
p_values = p_values_peons;

% --- Prettifying Options ---
% Define colors matching your figure
peonColor = [0.2 0.5 0.8]; % More similar to the blue in your figure
otherColor = [0.9 0.9 0.9]; % Very light gray for "Other"

% Define explosion factor (offset the first slice - PEONs)
explode = [0.1, 0]; % Subtle explosion for the PEON slice

% --- Create Figure ---
figure('Position', [100, 100, 900, 300]); % Width to fit 3 pie charts side by side

% --- Panel Label ---
annotation('textbox', [0.01, 0.9, 0.05, 0.05], 'String', 'F', ...
    'EdgeColor', 'none', 'FontSize', 14, 'FontWeight', 'bold');

% --- Area 1 (A1) Pie Chart ---
subplot(1, 3, 1);
dataA1 = [peonCounts(1), totalNeurons(1) - peonCounts(1)]; % Data: [PEONs, Others]
perc_peon_A1 = (peonCounts(1) / totalNeurons(1)) * 100;
perc_other_A1 = 100 - perc_peon_A1;

% Create pie chart with clean labels
hPieA1 = pie(dataA1, explode);
title(sprintf('%s (n=%d)', areaNamesShort{1}, totalNeurons(1)), 'FontSize', 11);

% Apply custom colors with thin black outline
hPieA1(1).FaceColor = peonColor;
hPieA1(1).EdgeColor = 'k';
hPieA1(1).LineWidth = 0.5;
hPieA1(3).FaceColor = otherColor;
hPieA1(3).EdgeColor = 'k';
hPieA1(3).LineWidth = 0.5;

% Set text labels
hPieA1(2).String = sprintf('PEONs\n(%d, %.1f%%)', peonCounts(1), perc_peon_A1);
hPieA1(4).String = sprintf('Other\n(%d, %.1f%%)', totalNeurons(1) - peonCounts(1), perc_other_A1);
hPieA1(2).FontSize = 9;
hPieA1(4).FontSize = 9;

% Add p-value annotation
if p_values(1) < 0.05
    text(0, -1.3, sprintf('p = %.4f*', p_values(1)), 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'Color', 'r');
else
    text(0, -1.3, sprintf('p = %.4f', p_values(1)), 'HorizontalAlignment', 'center', ...
        'FontSize', 9);
end

% --- Area 2 (VAF) Pie Chart ---
subplot(1, 3, 2);
dataA2 = [peonCounts(2), totalNeurons(2) - peonCounts(2)];
perc_peon_A2 = (peonCounts(2) / totalNeurons(2)) * 100;
perc_other_A2 = 100 - perc_peon_A2;

% Create pie chart
hPieA2 = pie(dataA2, explode);
title(sprintf('%s (n=%d)', areaNamesShort{2}, totalNeurons(2)), 'FontSize', 11);

% Apply custom colors
hPieA2(1).FaceColor = peonColor;
hPieA2(1).EdgeColor = 'k';
hPieA2(1).LineWidth = 0.5;
hPieA2(3).FaceColor = otherColor;
hPieA2(3).EdgeColor = 'k';
hPieA2(3).LineWidth = 0.5;

% Set text labels
hPieA2(2).String = sprintf('PEONs\n(%d, %.1f%%)', peonCounts(2), perc_peon_A2);
hPieA2(4).String = sprintf('Other\n(%d, %.1f%%)', totalNeurons(2) - peonCounts(2), perc_other_A2);
hPieA2(2).FontSize = 9;
hPieA2(4).FontSize = 9;

% Add p-value annotation
if p_values(2) < 0.05
    text(0, -1.3, sprintf('p = %.4f*', p_values(2)), 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'Color', 'r');
else
    text(0, -1.3, sprintf('p = %.4f', p_values(2)), 'HorizontalAlignment', 'center', ...
        'FontSize', 9);
end

% --- Area 3 (AAF) Pie Chart ---
subplot(1, 3, 3);
dataA3 = [peonCounts(3), totalNeurons(3) - peonCounts(3)];
perc_peon_A3 = (peonCounts(3) / totalNeurons(3)) * 100;
perc_other_A3 = 100 - perc_peon_A3;

% Create pie chart
hPieA3 = pie(dataA3, explode);
title(sprintf('%s (n=%d)', areaNamesShort{3}, totalNeurons(3)), 'FontSize', 11);

% Apply custom colors
hPieA3(1).FaceColor = peonColor;
hPieA3(1).EdgeColor = 'k';
hPieA3(1).LineWidth = 0.5;
hPieA3(3).FaceColor = otherColor;
hPieA3(3).EdgeColor = 'k';
hPieA3(3).LineWidth = 0.5;

% Set text labels
hPieA3(2).String = sprintf('PEONs\n(%d, %.1f%%)', peonCounts(3), perc_peon_A3);
hPieA3(4).String = sprintf('Other\n(%d, %.1f%%)', totalNeurons(3) - peonCounts(3), perc_other_A3);
hPieA3(2).FontSize = 9;
hPieA3(4).FontSize = 9;

% Add p-value annotation
if p_values(3) < 0.05
    text(0, -1.3, sprintf('p = %.4f*', p_values(3)), 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'Color', 'r');
else
    text(0, -1.3, sprintf('p = %.4f', p_values(3)), 'HorizontalAlignment', 'center', ...
        'FontSize', 9);
end

% --- Overall Figure Title ---
sgtitle('Distribution of PEONs Across Auditory Fields', 'FontSize', 12);

% --- Set consistent axis properties for all subplots ---
for i = 1:3
    subplot(1, 3, i);
    axis equal tight;
end
set(gcf, 'Renderer', 'painters');
exportgraphics(gcf, 'PEON_distribution_arean.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

