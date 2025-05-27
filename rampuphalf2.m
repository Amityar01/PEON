%% rampuphalf2.m
% Analyzes and plots the trial-by-trial buildup dynamics of omission responses
% for identified PEONs (Probability Encoding Omission Neurons).
% This script corresponds to the analysis shown in Figure 3A of the manuscript
% (Yaron et al., 2025). It bins omission responses by their sequence position
% (derived from 'allind') and fits a logistic growth curve to the mean responses.
%
% Inputs (expected to be in the MATLAB workspace, typically from 'FindPEONS.m'):
%   - PEONs_training: Vector of indices for neurons classified as PEONs.
%   - preferred_tone_direction: Vector indicating the preferred omission
%                               direction (-1 or 1) for each PEON. Used to
%                               align 'allind' data.
%   - allind: Cell array (neurons x probability_conditions) where each cell
%             contains a vector of sequence positions (timing information, e.g.,
%             trial number within a block) for each stimulus/omission. This is
%             used for binning responses. It's aligned based on preferred_tone_direction.
%   - PEON_responses_ordered: 3D matrix (identified_PEONs x num_testing_trials x num_probabilities)
%                             containing aligned trial-by-trial omission responses for
%                             PEONs on the testing trial set.
%   - testing_indices: Vector of trial indices to be used from
%                      PEON_responses_ordered (part of split-half).

%
% Outputs:
%   - Generates one or more figures displaying:
%     1. Mean omission responses binned by sequence position, with error bars.
%     2. A logistic growth curve fitted to these binned mean responses.
%     3. Annotations including the logistic function equation and R-squared value.
%   - Workspace variables created during execution (e.g., beta_hat, rsq,
%     binnedResponseMap).
%
% Dependencies:
%   - MATLAB's Statistics and Machine Learning Toolbox (for 'nlinfit').
%
% Example Usage:
%   % Ensure all required input variables from FindPEONS.m are in the workspace.
%   rampuphalf2;

% Author: Amit Yaron
%% --- PART 1: Logistic Fit to Binned Data ---
clear tim

% Process timing data for each condition and neuron
for i = 1:8
    for j = 1:length(PEONs_training)
        if preferred_tone_direction(j) == 1
            tim(i,j,:) = mod(allind{PEONs_training(j),i}, 1000);
        else
            tim(9-i,j,:) = mod(allind{PEONs_training(j),i}, 1000);
        end
    end
end
tim = tim(:,:,testing_indices);

% Data dimensions check
[numConditions, numNeurons, numTimePoints] = size(tim);
[numNeuronsData, numTimePointsData, numConditionsData] = size(PEON_responses_ordered);
if numNeurons ~= numNeuronsData || numTimePoints ~= numTimePointsData || numConditions ~= numConditionsData
    error('Dimensions of tim and PEON_responses_ordered do not match.');
end

close all

% Define dynamic bin edges (20 bins) and calculate bin centers
minTime = min(tim(:));
maxTime = max(tim(:));
binEdges = linspace(minTime, maxTime, 21); % 20 bins => 21 edges
binCenters = binEdges(1:end-1) + diff(binEdges)/2;

% Initialize a container for binned responses
binnedResponseMap = containers.Map('KeyType', 'double', 'ValueType', 'any');

% Group responses into bins for conditions 5 to 8
for cond = 5:8
    for neuron = 1:numNeurons
        neuron_times = squeeze(tim(cond, neuron, :));       % 1D time array
        neuron_responses = squeeze(PEON_responses_ordered(neuron, :, cond)); % 1D responses array
        
        for t = 1:numTimePoints
            current_time = neuron_times(t);
            current_response = neuron_responses(t);
            % Find the bin index for the current time
            binIdx = find(current_time >= binEdges(1:end-1) & current_time < binEdges(2:end), 1);
            if ~isempty(binIdx)
                binCenter = binCenters(binIdx);
                if isKey(binnedResponseMap, binCenter)
                    binnedResponseMap(binCenter) = [binnedResponseMap(binCenter), current_response];
                else
                    binnedResponseMap(binCenter) = current_response;
                end
            end
        end
    end
end

% Calculate mean and SEM for each bin
binnedMeanResponses = zeros(size(binCenters));
binnedSEMResponses = zeros(size(binCenters));
for b = 1:length(binCenters)
    if isKey(binnedResponseMap, binCenters(b))
        responses = binnedResponseMap(binCenters(b));
        binnedMeanResponses(b) = mean(responses);
        binnedSEMResponses(b) = std(responses) / sqrt(length(responses));
    else
        binnedMeanResponses(b) = NaN;
        binnedSEMResponses(b) = NaN;
    end
end

% Remove bins with NaN values
validIndices = ~isnan(binnedMeanResponses);
timeBins = binCenters(validIndices);
meansBins = binnedMeanResponses(validIndices);
stderrBins = binnedSEMResponses(validIndices);

% Define logistic function model and initial parameters
logistic_func1 = @(b, x) b(1) ./ (1 + exp(-b(2) * (x - b(3))));
initial_params1 = [mean(meansBins), 0.1, mean(timeBins)];

% Fit logistic model to binned data
[beta_hat1, ~, ~] = nlinfit(timeBins, meansBins, logistic_func1, initial_params1);

% Calculate R-squared value
yfit1 = logistic_func1(beta_hat1, timeBins);
SSresid1 = sum((meansBins - yfit1).^2);
SStotal1 = (length(meansBins) - 1) * var(meansBins);
rsq1 = 1 - SSresid1 / SStotal1;

% Generate high-resolution fitted curve
fitted_curve_time = linspace(min(timeBins), max(timeBins), 100);
fitted_curve1 = logistic_func1(beta_hat1, fitted_curve_time);

% Create first figure (polished style)
figure;
hold on;
% Use consistent marker and line styles
errorbar(timeBins, meansBins, stderrBins, 'ko', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
plot(fitted_curve_time, fitted_curve1, 'b-', 'LineWidth', 2);
title('Logistic Fit to Binned Data', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Number of Standards', 'FontSize', 12);
ylabel('Mean Response', 'FontSize', 12);
legend('Binned Data', 'Logistic Fit', 'Location', 'Best');
grid on;

% Annotate with logistic equation and R-squared
equation_text1 = sprintf('%.4f / (1 + exp(-%.4f * (t - %.4f)))\nR^2 = %.2f', ...
    beta_hat1(1), beta_hat1(2), beta_hat1(3), rsq1);
text(min(timeBins) + 0.1*range(timeBins), max(meansBins) - 0.1*range(meansBins), ...
    equation_text1, 'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k');
hold off;

% Set renderer and export as PDF
set(gcf, 'Renderer', 'painters');
exportgraphics(gcf, 'logistic_fit_binned_data.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%% --- PART 2: Aggregate Neuronal Response with Logistic Fit ---
% Compute aggregate data
dataAgg = mean(allAmat(PEONs_training,:,1), 3) + mean(allBmat(PEONs_training,:,7), 3);
timeAgg = linspace(1, 50, 50); 

% Compute mean and SEM at each timepoint
meansAgg = mean(dataAgg, 1);
stderrAgg = std(dataAgg, 0, 1) / sqrt(size(dataAgg, 1));

% Define logistic function and initial parameters
logistic_func2 = @(b, x) b(1) + (b(2) - b(1)) * exp(-b(3) * x);
initial_params2 = [1, 1, 0.1]; 

% Fit logistic model to aggregate data
[beta_hat2, ~, ~] = nlinfit(timeAgg, meansAgg, logistic_func2, initial_params2);

% Calculate R-squared value for the aggregate fit
yfit2 = logistic_func2(beta_hat2, timeAgg);
SSresid2 = sum((meansAgg - yfit2).^2);
SStotal2 = (length(meansAgg) - 1) * var(meansAgg);
rsq2 = 1 - SSresid2 / SStotal2;

% Generate high-resolution fitted curve for aggregate data
t_highres = linspace(min(timeAgg), max(timeAgg), 100);
fitted_curve2 = logistic_func2(beta_hat2, t_highres);

% Create second figure (polished style, matching first figure)
figure;
hold on;
errorbar(timeAgg, meansAgg, stderrAgg, 'ko', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
plot(t_highres, fitted_curve2, 'b-', 'LineWidth', 2);
title('Aggregate Neuronal Response with Logistic Fit', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time', 'FontSize', 12);
ylabel('Response', 'FontSize', 12);
legend('Mean Data', 'Logistic Fit', 'Location', 'Best');
grid on;

% Annotate with logistic equation and R-squared
A = beta_hat2(1);
B = beta_hat2(2);
k = beta_hat2(3);
equation_text2 = sprintf('y(t) = %.4f + %.4f * exp(-%.4f * t)\nR^2 = %.2f', A, B - A, k, rsq2);
text(min(timeAgg) + 0.05*range(timeAgg), max(meansAgg) - 0.1*range(meansAgg), ...
    equation_text2, 'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k');
hold off;

% Set renderer and export as PDF
set(gcf, 'Renderer', 'painters');
exportgraphics(gcf, 'aggregate_neuronal_response_logistic_fit.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
