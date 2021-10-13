function [Track_raw,Track_interp,varargout] = PALA_multiULM(IQ,listAlgo,ULM,PData,varargin)
%% [Track_raw,Track_interp,varargout] = PALA_multiULM(IQ_in,listAlgo,ULM,PData,varargin)
% performs ULM of 1 bloc of IQ for all ULM algorithms and store in cells.
%
% INPUTS:
%       - IQ: bloc of image where bubbles have to be detected
%       - listAlgo: list of algorithm s to be compared
%       - ULM: struct with various parameters
%       - PData: pixel informations of IQ images
%       - tracking (optional): 1/0
%       - savingFileName (optional): save data in 'savingFileName.mat'
% OUPUTS:
%       - Track_raw: list of non interpolated tracks. It returns the position of bubbles
%       and the frame index. It will be used later for pairing algorithm.
%       - Track_interp: list of interpolated tracks for density rendering.
%       - varargout: ProcessingTime
%
% Created by Arthur Chavignon 9/12/2019
%
% DATE 2020.07.22 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
% Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy, Nature Biomedical Engineering, 2021.
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018


if ~isfield(ULM,'parameters')% Create an empty structure for parameters hosting
    ULM.parameters = struct();
end
if ~isfield(ULM.parameters,'NLocalMax')
    ULM.parameters.NLocalMax = 3; % safeguard
end

tmp = strcmpi(varargin,'tracking');
if any(tmp),tracking= varargin{find(tmp)+1};else, tracking=1;end

tmp = strcmpi(varargin,'savingfilename');
if any(tmp),savingFileName = varargin{find(tmp)+1};SaveData=true;else, SaveData=false;end

% Convert the max_linking_distance into a appropriate scaling, the same as MatTracking
ULM.max_linking_distance = ULM.max_linking_distance*PData.PDelta(3);

%% Start localization for each compared algorithms.
IQ = abs(IQ);
Nalgo = numel(listAlgo);
ProcessingTime = zeros(Nalgo,1);

fprintf('Processing algo: ')
for ialgo = 1:Nalgo
    fprintf([num2str(ialgo) ' '])
    %% Detection and localization algorithm
    t0 = tic;

    switch lower(listAlgo{ialgo})
        case {'wa','weighted_average'}
            ULM.LocMethod = 'wa';

        case {'radial','radial_vivo','radial_silicio'}
            ULM.LocMethod = 'radial';

        case {'radial_sg'}
            ULM.LocMethod = 'radial_sg';

        case 'interp_cubic'
            ULM.LocMethod = 'interp';
            ULM.parameters.InterpMethod = 'cubic';

        case 'interp_lanczos'
            ULM.LocMethod = 'interp';
            ULM.parameters.InterpMethod = 'lanczos3';

        case {'interp_spline'}
            ULM.LocMethod = 'interp';
            ULM.parameters.InterpMethod = 'spline';

        case 'gaussian_fit'
            ULM.LocMethod = 'curvefitting';

        case {'interp_bilinear','no_localization','no_shift'}
            ULM.LocMethod = 'nolocalization';
        otherwise
            error('Wrong method selected')
    end
    [MatTracking] = ULM_localization2D(IQ,ULM);

    ProcessingTime(ialgo) = toc(t0);
    % MatTracking is the table that stores particles' values and positions
    % convert MatTracking from pixel to \lambda
    MatTracking(:,2:3) = (MatTracking(:,2:3) - [1 1]).*PData.PDelta([3 1]) + [PData.Origin(3) PData.Origin(1)]; % good origin

    %% Tracking algorithm
    if tracking
        [Track_raw{ialgo},Track_interp{ialgo}] = ULM_tracking2D(double(MatTracking),ULM,'pala');
    else
        Track_raw{ialgo} = single(MatTracking);
        Track_interp{ialgo} = [];
    end

    Track_interp{ialgo} = cellfun(@single,Track_interp{ialgo},'UniformOutput',false);
    Track_raw{ialgo} = cellfun(@single,Track_raw{ialgo},'UniformOutput',false);
end

if nargout ==3
    varargout{1} = ProcessingTime;
end

%% Save data in .mat file
if SaveData
    fprintf('saving... ')
    ProTime = ProcessingTime;
    save(savingFileName,'Track_raw','Track_interp','ProTime','ULM','PData','listAlgo','Nalgo','-v6')
end
fprintf('end.\n')
end
