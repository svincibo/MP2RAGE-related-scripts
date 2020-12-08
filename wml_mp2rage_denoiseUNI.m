% This script reads in FA, MD, AD, and RD measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated and concatenates
% those measurements with behavioral measurements.

clear all; close all; clc
format shortG

prefix = ''; % 'streamlinecount25920-'

% Set working directories.
rootDir = '/Volumes/Seagate/';

% Selected regularization parameter specifically for this
% dataset using the interactive feature of the
% DemoRemoveBackgroundNoise.m code.
regularization = 96;
            
% subjects who do not have mp2rage at either session
exclude = [001 002];

% subjects who have their inv2 and uni images in reverse order (in the series outputted from the scanner)
unusual_sub = [004 005 006 008 008 016 017];
unusual_ses = [2 1 1 1 2 2 1];

addpath(genpath(fullfile(rootDir, 'MP2RAGE-related-scripts')));

count = 0;

%% TRACTOGRAPHY.

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, 'WML-NIFTIS'));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
% grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's' & arrayfun(@(x) x.name(end), grp_contents) == num2str(s));

% Load in each tract's tractography measures for this subject.
for i = 1:size(grp_contents, 1)
    
    % Grab subID.
    sub(i) = str2num(grp_contents(i).name(4:6));
    
    % Grab session.
    ses(i) = str2num(grp_contents(i).name(end));
    
    % Display current sub ID and session.
    disp(grp_contents(i).name)
    
    % only include the subjects that we care about.
    if sum(sub(i) == exclude) == 0
        
        % Get contents of the directory where the tract measures for this subject are stored.
        sub_contents = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, '*mp2rage*.nii.gz'));
        
        % Remove the '.' and '..' files.
        sub_contents = sub_contents(arrayfun(@(x) x.name(1), sub_contents) ~= '.');
        
        if ~isempty(sub_contents)
            
            count = count + 1;
            
            %% DENOISE.
            
            % Account for inconsistency in the output order of the three mp2rage images.
            if ((sum(sub(i) == unusual_sub) > 0) && (ses(i) == unusual_ses(min(find(sub(i)==unusual_sub))))) ...
                    || ((sum(sub(i) == unusual_sub) > 0) && (ses(i) == unusual_ses(max(find(sub(i)==unusual_sub)))))
                
                a = 1; b = 3; c = 2;
                
            else
                
                a = 1; b = 2; c = 3;
                
            end
            
            disp([a b c]);
            disp(size(sub_contents, 1));
            
            % Remove background noise.
            MP2RAGE.filenameUNI=fullfile(sub_contents(c).folder, sub_contents(c).name); % standard MP2RAGE T1w image;
            MP2RAGE.filenameINV1=fullfile(sub_contents(a).folder, sub_contents(a).name);% Inversion Time 1 MP2RAGE T1w image;
            MP2RAGE.filenameINV2=fullfile(sub_contents(b).folder, sub_contents(b).name);% Inversion Time 2 MP2RAGE T1w image;
            MP2RAGE.filenameOUT=fullfile(rootDir, ['data_mp2rage_reg' num2str(regularization)], ['sub' num2str(sub(i), '%03.f') '-' num2str(ses(i)) '-MP2RAGE_denoised.nii.gz']);% image without background noise;
            
            [MP2RAGEimgRobustPhaseSensitive]=RobustCombination(MP2RAGE, regularization);
            
            saveas(gcf, fullfile(rootDir, ['data_mp2rage_reg' num2str(regularization)], ['sub' num2str(sub(i), '%03.f') '-' num2str(ses(i)) '-MP2RAGE_denoised.png']));
            close all
            
        end % if exist
        
    end %exclude
    
end % end i

