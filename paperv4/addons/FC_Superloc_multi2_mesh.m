function MatTrackingTot = FC_Superloc_multi2_mesh(IQ_in,listAlgo,ULM,PData)
%% function MatTrackingTot = FC_Superloc_multi2_mesh(IQ_in,listAlgo,ULM,PData)
% Perfom ULM of 1 bloc of IQ for all ULM algorithms and store in cells.
%
% INPUTS:
%       - IQ: bloc of image where bubbles are located in the center
%       - listAlgo: list of algorithm s to be compared
%       - ULM: struct with various parameters
%       - PData: pixel informations of IQ images
% OUPUTS:
%       - MatTrackingTot: scatter position [z x] for each frame, each algo
%           Dim 1 : position [z x]
%           Dim 2 : algo [1:numel(listAlgo)]
%           Dim 3 : frame
%
% Created by Arthur Chavignon 9/12/2019

MatTrackingTot = {};

for ialgo = 1:numel(listAlgo)
    
    switch lower(listAlgo{ialgo})
        case {'saf','saf_secure','saf_silicio','wa'}
            ULM.LocMethod = 'WA';
            [MatTracking] = ULM_Superloc2D_mesh(IQ_in,ULM);    
            
        case {'radial','radial_vivo','radial_silicio'}
            ULM.LocMethod = 'Radial';
            [MatTracking] = ULM_Superloc2D_mesh(IQ_in,ULM);    
            
        case {'interp_cubic','interp_cubic_silicio'}
            ULM.LocMethod = 'Interp';
            ULM.parameters.InterpMethod = 'cubic';
            [MatTracking] = ULM_Superloc2D_mesh(IQ_in,ULM);    

        case {'interp_lanczos','interp_lanczos_silicio'}
            ULM.LocMethod = 'Interp';
            ULM.parameters.InterpMethod = 'lanczos3';
            [MatTracking] = ULM_Superloc2D_mesh(IQ_in,ULM);    
            
         case {'interp_spline','interp_spline_silicio'}
            ULM.LocMethod = 'Interp';
            ULM.parameters.InterpMethod = 'spline';
            [MatTracking] = ULM_Superloc2D_mesh(IQ_in,ULM);    
            
        case 'gaussian_fit'
            ULM.LocMethod = 'CurveFitting';
            [MatTracking] = ULM_Superloc2D_mesh(IQ_in,ULM);    

        case {'interp_bilinear','interp_bilinear_silicio','no_localization','no_shift'}
            ULM.LocMethod = 'NoLocalization';
            [MatTracking] = ULM_Superloc2D_mesh(IQ_in,ULM);    
        otherwise
            error('Wrong method selected')
    end
    
    % convert MatLTrakcing into \lambda
    MatTracking = (MatTracking-[1 1]).*PData.PDelta([3 1]) + [PData.Origin(3) PData.Origin(1)]; % good origin
    MatTrackingTot{ialgo} = MatTracking;  % [zpos xpos]
    clear MatTracking
end
disp('Localization perfomed') 

%% Cat and permute
% Dim 1 : position [z x]
% Dim 2 : algo [1:numel(listAlgo)]
% Dim 3 : frame
MatTrackingTot = cat(3,MatTrackingTot{:}); % dim : [frame, pos, ialgo]
MatTrackingTot = permute(MatTrackingTot,[2 3 1]); % [pos, ialgo, frame]

return 


