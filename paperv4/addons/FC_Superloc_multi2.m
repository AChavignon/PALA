function [Track_tot,Track_tot_interp,varargout] = FC_Superloc_multi2(IQ,listAlgo,ULM,PData,FrameRate,varargin)
%% [Track_tot,Track_tot_interp,varargout] = FC_Superloc_multi2(IQ_in,listAlgo,ULM,PData,FrameRate,varargin)
% Perfom ULM of 1 bloc of IQ for all ULM algorithms and store in cells.
%
% INPUTS:
%       - IQ: bloc of image where bubbles have to be detected
%       - listAlgo: list of algorithm s to be compared
%       - ULM: struct with various parameters
%       - PData: pixel informations of IQ images
%       - FrameRate: framerate. 
% OUPUTS:
%       - Track_tot: list of tracks non interpolated. It returns the position of bubbles
%       and the frame index. It will be used later for pairing algorihm.
%       - Track_tot_interp: list of interpolated tracks for density rendering.
%       - varargout: ProcessingTime
%
% Created by Arthur Chavignon 9/12/2019
%
% DATE 2020.03.30 - VERSION 1.0.0
% AUTHROS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris  
% Code Available under Creative Commons Non-Commercial 4.0
% ACADEMIC REFERENCES TO BE CITED
% Details of the code published in 2020 article by Heiles, Chavignon, Hingot and Couture.
% Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018 

res = ULM.res;
ProcessingTime = [];

ULM.parameters.NLocalMax = 3; % safeguard

tmp = strcmpi(varargin,'tracking'); 
if any(tmp),tracking= varargin{find(tmp)+1};else, tracking=1;end

%% Start localization for each compared algorihms.
IQ = abs(IQ);

for ialgo = 1:numel(listAlgo)
    %% Detection and localization algorithm
    t0 = tic;
    
    switch lower(listAlgo{ialgo})
        case {'saf','saf_secure','saf_silicio','wa'}
            ULM.LocMethod = 'WA';
            [MatTracking] = ULM_Superloc2D(IQ,ULM);
            
        case {'radial','radial_vivo','radial_silicio'}
            ULM.LocMethod = 'Radial';
            [MatTracking] = ULM_Superloc2D(IQ,ULM);
            
        case {'interp_cubic','interp_cubic_silicio'}
            ULM.LocMethod = 'Interp';
            ULM.parameters.InterpMethod = 'cubic';
            [MatTracking] = ULM_Superloc2D(IQ,ULM);

        case {'interp_lanczos','interp_lanczos_silicio'}
            ULM.LocMethod = 'Interp';
            ULM.parameters.InterpMethod = 'lanczos3';
            [MatTracking] = ULM_Superloc2D(IQ,ULM);
            
         case {'interp_spline','interp_spline_silicio'}
            ULM.LocMethod = 'Interp';
            ULM.parameters.InterpMethod = 'spline';
            [MatTracking] = ULM_Superloc2D(IQ,ULM);  
            
        case 'gaussian_fit'
            ULM.LocMethod = 'CurveFitting';
            [MatTracking] = ULM_Superloc2D(IQ,ULM);

        case {'interp_bilinear','interp_bilinear_silicio','no_localization','no_shift'}
            ULM.LocMethod = 'NoLocalization';
            [MatTracking] = ULM_Superloc2D(IQ,ULM);    
        otherwise
            error('Wrong method selected')
    end
    ProcessingTime(ialgo) = toc(t0);
    % MatTracking is the table that stores the paricles values and position
%     MatTracking(:,2:3) = MatTracking(:,2:3)/ULM.res - [1 1]; % convert into \lambda
    % convert MatTracking from pixel to \lambda
    MatTracking(:,2:3) = (MatTracking(:,2:3) - [1 1]).*PData.PDelta([3 1]) + [PData.Origin(3) PData.Origin(1)]; % good origin
    
    %% Tracking algotihm
    if tracking
        max_linking_distance = ULM.max_linking_distance/res/PData.PDelta(3);
        [Track_tot{ialgo},Track_tot_interp{ialgo}] = ULM_tracking2D(double(MatTracking),ULM.size(3),FrameRate,max_linking_distance,ULM.max_gap_closing,ULM.min_length,'FightClub'); % FightClub noInterp
    else
        Track_tot{ialgo} = MatTracking;
        Track_tot_interp{ialgo} = [];
    end
        
end
    

if nargout ==3
    varargout{1} = ProcessingTime;
end

return 

% %% Debugg display position IQ
% % ifr = 15;
% ifr = ifr+1
% figure(14),
% hold off
% imagesc(abs(IQ_in(:,:,ifr)))
% ap  = MatTracking((1:ULM.numberOfParticles)+(ifr-1)*ULM.numberOfParticles,2:3)/ULM.res;
% hold on
% plot(ap(:,2),ap(:,1),'r.')

