clear all;
tic;

% Parameters of the generated image
pixelsize = 1;
Nx = 1000; % Number of pixels in x-direction
Ny = 1000; % Number of pixels in y-direction
meandgpx = 80.0; 
meandgpy = 80.0;
sddgpx = 5.0; 
sddgpy = 5.0; 
meanwgpx = 10.0; 
meanwgpy = 10.0;
sdwgpx = 0.00; 
sdwgpy = 0.00;
meanhgp = 10.0; % Keep meanhgp constant
qexp = 2;

% Define the values for sdhgp
sdhgp_values = [1, 2, 3, 4, 5]; % Different sdhgp values to test

% Number of iterations
numIterations = 10;

% Threshold values from 0.2 to 0.95
thresholds = 0.05:0.05:0.95;
numThresholds = length(thresholds);

% Initialize arrays to store results for all sdhgp values
num_sdhgp_values = length(sdhgp_values);
meanNNIs_all = zeros(num_sdhgp_values, numThresholds);
stdNNIs_all = zeros(num_sdhgp_values, numThresholds);
meanNumObjects_all = zeros(num_sdhgp_values, numThresholds);
stdNumObjects_all = zeros(num_sdhgp_values, numThresholds);

% Area of the region (for NNI calculation)
A = 800 * 800;

% Create a folder to save images if it doesn't already exist
if ~exist('iteration_images_periodic', 'dir')
   mkdir('iteration_images_periodic');
end

% Loop over different sdhgp values
for sdhgp_idx = 1:num_sdhgp_values
    sdhgp = sdhgp_values(sdhgp_idx);
    fprintf('sdhgp: %d\n', sdhgp);

    ngpx = floor(Nx / meandgpx) - 1; % Number of mounds in x-direction
    ngpy = floor(Ny / meandgpy) - 1; % Number of mounds in y-direction

    % To store Nearest Neighbor Index (NNI) and object counts for all iterations and thresholds
    NNIs = zeros(numIterations, numThresholds);
    numObjectsPerIteration = zeros(numIterations, numThresholds);

    % Run multiple iterations of point pattern generation and analysis
    for iter = 1:numIterations
        fprintf('Iteration: %d\n', iter);

        % Initializing the surface/image to be generated
        h = zeros(Ny, Nx); % Surface height map
        hgp = zeros(ngpy, ngpx); 
        wgpx = zeros(ngpy, ngpx); 
        wgpy = zeros(ngpy, ngpx);

        % Step 1: Generate periodic mound-like surface/image
        for igpx = 1:ngpx
            for igpy = 1:ngpy
                % Generate mound positions and heights
                xgp(igpy, igpx) = igpx * meandgpx + round(randn * sddgpx);  
                ygp(igpy, igpx) = igpy * meandgpy + round(randn * sddgpy); 
                hgp(igpy, igpx) = meanhgp + sdhgp * randn; % Use varying sdhgp
                wgpx(igpy, igpx) = round(meanwgpx + sdwgpx * randn);
                wgpy(igpy, igpx) = wgpx(igpy, igpx); % Isotropy assumption

                % Ensure width values are positive
                if (wgpx(igpy, igpx) <= 0) wgpx(igpy, igpx) = 1; end 
                if (wgpy(igpy, igpx) <= 0) wgpy(igpy, igpx) = 1; end

                % Boundaries for the mound region
                kgpminx = max(xgp(igpy, igpx) - 12 * wgpx(igpy, igpx), 1);
                kgpmaxx = min(xgp(igpy, igpx) + 12 * wgpx(igpy, igpx), Nx);
                kgpminy = max(ygp(igpy, igpx) - 12 * wgpy(igpy, igpx), 1);
                kgpmaxy = min(ygp(igpy, igpx) + 12 * wgpy(igpy, igpx), Ny);

                % Add each mound to the height map, keeping intensity of one object when merging
                for kgpx = kgpminx:kgpmaxx
                    for kgpy = kgpminy:kgpmaxy
                        newIntensity = hgp(igpy, igpx) * ...
                            exp(-(sqrt(((kgpx - xgp(igpy, igpx)) / (2 * wgpx(igpy, igpx)))^2 + ...
                                       ((kgpy - ygp(igpy, igpx)) / (2 * wgpy(igpy, igpx)))^2))^qexp);
                        
                        % Update pixel intensity only if the new intensity is greater
                        if h(kgpy, kgpx) < newIntensity
                            h(kgpy, kgpx) = newIntensity;
                        end
                    end
                end
            end
        end

        % Step 2: Normalize the height map to create an image with black peaks
        Is = mat2gray(h(101:900, 101:900));  % Crop the image to an 800x800 section
        Is = 1 - Is;  % Invert the grayscale image to make peaks black

        % Display the image
        imshow(Is);

        % Analyze image across all threshold levels
        for idx = 1:numThresholds
            % Step 3: Apply threshold to binarize the image
            threshold = thresholds(idx);
            binaryImage = Is < threshold;  % Threshold at each value (black peaks)
            binaryImage = imclearborder(binaryImage);  % Remove objects touching the border


            % Step 4: Calculate the centroids of the binary objects
            props = regionprops(binaryImage, 'Centroid');
            centroids = cat(1, props.Centroid);
            numObjects = size(centroids, 1);
            numObjectsPerIteration(iter, idx) = numObjects; % Save the number of objects

            if numObjects > 1
                % Step 5: Calculate observed nearest neighbor distances
                distances = pdist2(centroids, centroids);
                nearestDistances = zeros(numObjects, 1);

                for i = 1:numObjects
                    distances(i, i) = inf;  % Exclude self
                    nearestDistances(i) = min(distances(i, :));  % Find the nearest neighbor
                end

                % Mean observed nearest neighbor distance
                d_obs = mean(nearestDistances);

                % Step 6: Calculate expected nearest neighbor distance for a random pattern
                lambda = numObjects / A;  % Point density
                d_exp = 1 / (2 * sqrt(lambda));  % Expected nearest neighbor distance

                % Step 7: Calculate the Nearest Neighbor Index (NNI)
                NNIs(iter, idx) = d_obs / d_exp;
            else
                NNIs(iter, idx) = NaN;  % Set NNI to NaN if there are no or only one object
            end

            % Save the binary image for each iteration and threshold
            imageFileName = sprintf('iteration_images_periodic/Iteration_%d_Threshold_%.2f_sdhgp_%d.png', iter, threshold, sdhgp);
            imwrite(binaryImage, imageFileName);
        end
    end

    % Compute mean and standard deviation of NNIs across iterations
    meanNNIs = mean(NNIs, 'omitnan');
    stdNNIs = std(NNIs, 'omitnan');

    % Compute mean and standard deviation of object counts across iterations
    meanNumObjects = mean(numObjectsPerIteration, 'omitnan');
    stdNumObjects = std(numObjectsPerIteration, 'omitnan');

    % Store the results for the current sdhgp value
    meanNNIs_all(sdhgp_idx, :) = meanNNIs;
    stdNNIs_all(sdhgp_idx, :) = stdNNIs;
    meanNumObjects_all(sdhgp_idx, :) = meanNumObjects;
    stdNumObjects_all(sdhgp_idx, :) = stdNumObjects;
end


% Plot mean NNI vs Threshold for all sdhgp values with error bars for standard deviation
figure;
hold on;
colors = lines(num_sdhgp_values); % Get different colors for plotting
for sdhgp_idx = 1:num_sdhgp_values
    errorbar(thresholds, meanNNIs_all(sdhgp_idx, :), stdNNIs_all(sdhgp_idx, :), ...
        '-o', 'Color', colors(sdhgp_idx, :), 'LineWidth', 2);
end
xlabel('Threshold', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Mean Nearest Neighbor Index (NNI)', 'FontSize', 16, 'FontWeight', 'bold');
%title('Mean NNI vs Threshold for Different sdhgp Values', 'FontSize', 18, 'FontWeight', 'bold');
legendCell = cellstr(num2str(sdhgp_values', 'std depth = %.1f'));
legend(legendCell, 'Location', 'Best');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');  % Set axis font size and weight
grid on;
hold off;

% Save the figure
saveas(gcf, 'NNI_vs_Threshold_All_sdhgp_c10std10.png');

% Plot mean number of objects vs Threshold for all sdhgp values with error bars for standard deviation
figure;
hold on;
for sdhgp_idx = 1:num_sdhgp_values
    errorbar(thresholds, meanNumObjects_all(sdhgp_idx, :), stdNumObjects_all(sdhgp_idx, :), ...
        '-o', 'Color', colors(sdhgp_idx, :), 'LineWidth', 2);
end
xlabel('Threshold', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Mean Number of Objects', 'FontSize', 16, 'FontWeight', 'bold');
%title('Mean Number of Objects vs Threshold for Different sdhgp Values', 'FontSize', 18, 'FontWeight', 'bold');
legend(legendCell, 'Location', 'Best');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
grid on;
hold off;

% Save the figure
saveas(gcf, 'MeanNumObjects_vs_Threshold_All_sdhgp_c10std10.png');

% New: Plot Mean NNI vs Number of Objects for all sdhgp values
figure;
hold on;
for sdhgp_idx = 1:num_sdhgp_values
    % Extract mean number of objects and mean NNI for this sdhgp value
    meanNumObjects = meanNumObjects_all(sdhgp_idx, :);
    meanNNIs = meanNNIs_all(sdhgp_idx, :);

    % Plot NNI vs Number of Objects
    plot(meanNumObjects, meanNNIs, '-o', 'Color', colors(sdhgp_idx, :), 'LineWidth', 2, ...
        'DisplayName', sprintf('std depth= %.1f', sdhgp_values(sdhgp_idx)));
end
xlabel('Mean Number of Objects', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Mean Nearest Neighbor Index (NNI)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'Best');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
grid on;
hold off;

% Save the figure
saveas(gcf, 'NNI_vs_NumberOfObjects_All_sdhgpc10std10.png');

toc;


