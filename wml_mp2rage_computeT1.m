% This script reads in FA, MD, AD, and RD measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated and concatenates
% those measurements with behavioral measurements.

clear all; close all; clc
format shortG

% Set working directories and add needed paths.
rootDir = '/Volumes/Seagate/';
blprojectid = 'proj-5f4d3d2a42c1728f7d6440b0';
addpath(genpath(fullfile(rootDir, 'MP2RAGE-related-scripts')));
binsize = 10000;

% subjects who do not have mp2rage at either session
exclude = [001 002];

% subjects who have their inv2 and uni images in reverse order (in the series outputted from the scanner)
unusual_sub = [004 005 006 008 008 016 017];
unusual_ses = [2 1 1 1 2 2 1];

%% Construct B1 sensitivity plot for our specific acquisition parameters.

% Assign constants related to the mp2rage acquisition protocol.
MP2RAGE.B0=3;           % in Tesla
MP2RAGE.TR=5;           % MP2RAGE TR in seconds
MP2RAGE.TRFLASH=7.5e-3; % TR of the GRE readout, per Hu: (TI2-TI1)/(number of phase encoding steps) = (2500-700)/(256*93.8%) ms = 7.5 ms.
MP2RAGE.TIs=[700e-3 2500e-3];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
MP2RAGE.NZslices=[88 88];% Slices Per Slab * [PartialFourierInSlice-0.5  0.5] per Hu: PartialFourierInSlice = 1, so 176*[1-0.5 0.5] = 176*[0.5 0.5] = [88 88].
% NOTE: Franco thinks that NZslices should identify a subset of slices within the
% slab that does not contain empty space(i.e., non-sero), so maybe something like [70 90];
MP2RAGE.FlipDegrees=[4 5];% Flip angle of the two readouts in degrees

figure(1)
plotMP2RAGEproperties(MP2RAGE)
saveas(gcf, fullfile(rootDir, 'data_mp2rage', 'MP2RAGE_curve.png'));

%% Calculate T1map and R1map for each subject/session.

% Get contents of the directory where the images for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Initialize included subject counter.
count = 0;

% Read in denoised UNI images, i.e., mp2rage images, and brain mask images for each subject/session one at a time.
for i = 1:size(grp_contents, 1)
    
    % Grab subID.
    sub(i) = str2num(grp_contents(i).name(5:7));
    
    % Grab session.
    ses(i) = str2num(grp_contents(i).name(end));
    
    % Display current sub ID and session.
    disp(grp_contents(i).name)
    
    % only include the subjects that we care about.
    if sum(sub(i) == exclude) == 0
        
        % Update included subject counter.
        count = count + 1;
        
        %% Read in the denoised UNI images, i.e., the mp2rage image.
        
        % Get contents of the directory where the denoised mp2rage for this subject is stored.
        sub_mp2rage = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, '*mp2rage*/t1.nii.gz'));
        
        % Remove the '.' and '..' files.
        sub_mp2rage = sub_mp2rage(arrayfun(@(x) x.name(1), sub_mp2rage) ~= '.');
        
        % Assign file paths according to how the mp2rage code expects them to be assigned.
        MP2RAGE.filename = fullfile(sub_mp2rage(1).folder, sub_mp2rage(1).name); % standard MP2RAGE T1w image;
        MP2RAGE.filenameOUT = fullfile(rootDir, 'data_mp2rage', ['sub' num2str(sub(i), '%03.f') '-' num2str(ses(i)) '-MP2RAGE_curve.nii.gz']);% image without background noise;
        
        %% Read in the brain mask for the mp2rage image.
        
        % Get contents of the directory where the denoised mp2rage for this subject is stored.
        sub_mask = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, '*brain_extracted*/t1.nii.gz'));
        
        % Remove the '.' and '..' files.
        sub_mask = sub_mask(arrayfun(@(x) x.name(1), sub_mask) ~= '.');
        
        clear mask;
        mask = niftiRead(fullfile(sub_mask(1).folder, sub_mask(1).name));
        
        % Convert mask.data class from int32 to double for compatability with MP2RAGEimg.
        mask.data = double(mask.data);
        
        %% Calculate T1map and R1map.
        
        % load the MP2RAGE data - it can be either the SIEMENS one scaled from
        % 0 4095 or the standard -0.5 to 0.5
        clear MP2RAGEimg;
        MP2RAGEimg = load_untouch_nii(MP2RAGE.filename);
        
        % Caluclate the T1 and R1 maps for this subject.
        clear T1map R1map
        [T1map, R1map] = T1estimateMP2RAGE(MP2RAGEimg, MP2RAGE, 0.96);
        
        % Convert MP2RAGEimg.img class from int16 to double for compatability with mask.data.
        MP2RAGEimg.img = double(MP2RAGEimg.img);
        
        %% Plot histogram fit for this subject/session.
        
        % Open figure.
        if count == 1; figure(2); hold on; end
        
        % Set line properties for session 1 and session 2 differently.
        if ses(i) == 1; linecolor = 'r'; elseif ses(i) == 2; linecolor = 'b'; end
        
        % Plot histogram for this subject/session.
        subplot(3, 1, 1)
        toplot = MP2RAGEimg.img(:).*mask.data(:);
        toplot(find(toplot==0)) = NaN;
        temp = toplot(:);
        toplot2 = temp(~isnan(temp)); clear toplot;
%         histogram(toplot2, 'Normalization', 'pdf')
        [N, edges] = histcounts(toplot2, 'Normalization', 'probability');
        edges = edges(2:end) - (edges(2)-edges(1))/2;
        plot(edges, N, linecolor);
        pbaspect([1 1 1]);
        if count == 1; ylabel('Frequency'); xlabel('Signal Intensity'); hold on; end
        clear N edges;
%         if i == size(grp_contents, 1); pbaspect([5 3 1]); hold on; end
        
        subplot(3, 1, 2)
        toplot = T1map.img(:).*mask.data(:);
        toplot(find(toplot==0)) = NaN;
        temp = toplot(:);
        toplot2 = temp(~isnan(temp)); clear toplot;
%         histogram(toplot2, 'Normalization', 'pdf')
        [N, edges] = histcounts(toplot2, 'Normalization', 'probability');
        edges = edges(2:end) - (edges(2)-edges(1))/2;
        plot(edges, N, linecolor);
        axis([0 4 0 0.04]);
        pbaspect([1 1 1]);
        if count == 1; ylabel('Frequency'); xlabel('T1, seconds'); hold on; end
        clear N edges;
%         if i == size(grp_contents, 1); pbaspect([5 3 1]); hold on; end
        
        subplot(3, 1, 3)
        toplot = R1map.img(:).*mask.data(:);
        toplot(find(toplot==0)) = NaN;
        temp = toplot(:);
        toplot2 = temp(~isnan(temp)); clear toplot;
%         histogram(toplot2, 'Normalization', 'pdf')
        [N, edges] = histcounts(toplot2, 'Normalization', 'probability');
        edges = edges(2:end) - (edges(2)-edges(1))/2;
        plot(edges, N, linecolor);
        axis([0 2 0 0.03]);
        pbaspect([1 1 1]);
        if count == 1; ylabel('Frequency'); xlabel('R1, (seconds^-^1)'); hold on; end
        clear N edges;
%         if i == size(grp_contents, 1); pbaspect([5 3 1]); hold on; end

        %% Concatenate the T1, R1, and mp2rage intensity maps --for only brain tissue-- for this subject with previous subjects.
        T1(:, count) = T1map.img(:).*mask.data(:);
        R1(:, count) = R1map.img(:).*mask.data(:);
        mp2rage_intensity(:, count) = MP2RAGEimg.img(:).*mask.data(:);
        
    end %exclude
    
end % end i

hold off;
saveas(gcf, fullfile(rootDir, 'data_mp2rage', 'MP2RAGE_intersubjectreproducibility.png'));

% Convert zeros to NaN.
T1(find(T1==0)) = NaN;
R1(find(R1==0)) = NaN;
mp2rage_intensity(find(mp2rage_intensity==0)) = NaN;

% Remove NaN.
temp = T1(:);
T1_out = temp(~isnan(temp)); clear temp;
temp = R1(:);
R1_out = temp(~isnan(temp)); clear temp;
temp = mp2rage_intensity(:);
mp2rage_intensity_out = temp(~isnan(temp)); clear temp;

%% Plot summary figure, i.e., all subjects and all sessions together.
figure(3)

hold on;
facecolor = [0.4660 0.6740 0.1880]; %green
facealpha = 0.4;
edgecolor = 'none';

linecolor = [0.4940 0.1840 0.5560]; %purple
linewidth = 3;

subplot(3, 1, 1)
clear toplot;
toplot = double(mp2rage_intensity_out(:));
% histfit(toplot, 1000, 'kernel')
histogram(toplot, 'Normalization', 'probability', 'EdgeColor', edgecolor, 'FaceColor', facecolor, 'FaceAlpha', facealpha);
hold on;
[N, edges] = histcounts(toplot, 'Normalization', 'probability');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N, 'Color', linecolor, 'LineWidth', linewidth);
% axis([0 max(toplot) 0 13000])
pbaspect([1 1 1]);
ylabel('Frequency');
xlabel('Signal Intensity');
hold off;

subplot(3, 1, 2)
toplot = T1_out(:);
% histfit(toplot, 1000, 'kernel')
histogram(toplot, 'Normalization', 'probability', 'EdgeColor', edgecolor, 'FaceColor', facecolor, 'FaceAlpha', facealpha);
hold on;
[N, edges] = histcounts(toplot, 'Normalization', 'probability');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N, 'Color', linecolor, 'LineWidth', linewidth);
% axis([0 max(T1_out) 0 13000])
pbaspect([1 1 1]);
ylabel('Frequency');
xlabel('T1, seconds');
hold off;

subplot(3, 1, 3)
toplot = R1_out(:);
% % histfit(toplot, 1000, 'kernel')
histogram(toplot, 'Normalization', 'probability', 'EdgeColor', edgecolor, 'FaceColor', facecolor, 'FaceAlpha', facealpha);
hold on;
[N, edges] = histcounts(toplot, 'Normalization', 'probability');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N, 'Color', linecolor, 'LineWidth', linewidth);
% axis([0 max(R1_out) 0 13000])
pbaspect([1 1 1]);
ylabel('Frequency');
xlabel('R1, (seconds^-^1)');
hold off;

saveas(gcf, fullfile(rootDir, 'data_mp2rage', 'MP2RAGE_summary.png'));
