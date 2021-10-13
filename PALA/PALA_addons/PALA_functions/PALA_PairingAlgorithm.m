function [Stat_classification,ErrList,FinalPairs,MissingPoint,WrongLoc] = PALA_PairingAlgorithm(ListPos_ref,MatTrackedLoc,PData,Threshold_pairing,varargin)
%% function [Stat_classification,ErrList,FinalPairs,MissingPoint,WrongLoc] = PALA_PairingAlgorithm(ListPos_ref,MatTrackedLoc,PData,Threshold_pairing,varargin)
% This function paired the detected bubbles (MatTrackedLoc) with the reel simulated
% bubbles (ListPos_ref) by calculating the interdistance matrix for each frame. Bubbles
% are paired when the RMSE in less than Threshold_pairing.
% INPUTS:
%       - ListPos_ref: list of target simulated bubbles [x y z reflectivity]
%       - MatTrackedLoc: list of detected bubbles
%       - PData: parameters of IQ.
%       - Threshold_pairing: max distance between reel and detected point to consider a pair assignment. (\lambda/2)
%       - Threshold_TruePos: max distance between reel and detected point to consider a true positive localization. (varargin{1})
% OUTPUTS:
%       - Stat_classification: list of counting [Npos_in x Npos_loc x T_pos x F_neg x F_pos]
%           Npos_in: number of simulated bubbles, Npos_loc nb of localized bubbles
%           T_pos true position, F_neg false negative, F_pos false positive detection
%       - ErrList: list of localization [dz,dx,norm(dz,dx)]
%       - FinalPairs: list of good pairings, index ListPos_ref, index in MatTrackedLoc, error RMSE
%           [index ref, index loc, RMSE]
%       - MissingPoint: index of missed points (in ListPos_ref)
%       - WrongLoc: index of wrong detection (in MatTrackedLoc)
%
% Created by Arthur Chavignon 9/12/2019
%
% DATE 2020.08.10 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
% Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy, Nature Biomedical Engineering, 2021.
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018

if ~isempty(varargin)
    Threshold_TruePos = varargin{1};
    if Threshold_TruePos>Threshold_pairing
        error('Threshold_truepos must be smaller than Threshold_pairing')
    end
else
    Threshold_TruePos = Threshold_pairing;
end

ErrList = [];
for ifr = 1:size(ListPos_ref,3)
    %% for 1 Frame, Get positions of simulated points
    pos_ref = ListPos_ref(:,:,ifr);
    pos_ref = pos_ref(isfinite(pos_ref(:,1)),:);
    % remove points outside the image
    outLiers = true(size(pos_ref,1),1);
    outLiers(pos_ref(:,1)<PData.Origin(1))=0;
    outLiers(pos_ref(:,1)>PData.Origin(1)+(PData.Size(2)+1)*PData.PDelta(1))=0;
    outLiers(pos_ref(:,3)<PData.Origin(3))=0;
    outLiers(pos_ref(:,3)>PData.Origin(3)+(PData.Size(1)+1)*PData.PDelta(3))=0;
    pos_ref = pos_ref(outLiers,:);

    %% change the columns organization [z x y reflectivity]
    pos_ref = pos_ref(:,[3 1 2 4]);
    Points_ref{ifr} = pos_ref;
    Nb_ref = size(pos_ref,1); %number of target points

    %% get localized points after ULM
    idframe = find(MatTrackedLoc(:,3)==ifr);
    pos_loc = MatTrackedLoc(idframe,:); clear idframe
    Points_loc{ifr} = pos_loc;
    Nb_loc =  size(pos_loc,1); % number of localized points

    %% Compute inter-distances matrix
    InterDistMat = zeros(Nb_ref,Nb_loc);
    for i_0=1:Nb_ref
        InterDistMat(i_0,:) = vecnorm(pos_ref(i_0,[1 2])-pos_loc(:,[1 2]),2,2);
    end
    % Sorting and pairing
    InterDistMat_0 = InterDistMat;

    [minErr_loc,~] = min(InterDistMat_0,[],1);
    [~,indin_loc] = sort(minErr_loc(:));
    InterDistMat_0 = InterDistMat_0(:,indin_loc);

    [minErr_ref,~] = min(InterDistMat_0,[],2);
    [~,indin_ref] = sort(minErr_ref(:));
    InterDistMat_0 = InterDistMat_0(indin_ref,:);

    DD_mat_00 = InterDistMat_0;
    imin_last = 0;

    MatPairings = zeros(0,3); % pairing matrix [indice_ref,indice_loc,distance]
    for ii=1:size(InterDistMat_0,2)
        if isempty(DD_mat_00)
           break
        end
        [err,imin] = min(DD_mat_00(:,1));
        MatPairings = cat(1,MatPairings,[indin_ref(imin_last+imin),indin_loc(ii),err]);
        DD_mat_00 = DD_mat_00((imin+1):end,2:end);
        imin_last = imin_last+imin;
    end

    %% Store data in outputs
    GoodPairings        = MatPairings(MatPairings(:,3)<Threshold_pairing,:);
    FinalPairs{ifr}     = GoodPairings;     % [index ref, index loc, RMSE]
    MissingPoint{ifr}   = setdiff(1:Nb_ref,GoodPairings(:,1)); % [index ref]
    WrongLoc{ifr}       = setdiff(1:Nb_loc,GoodPairings(:,2)); % [index loc]
    ErrList             = cat(1,ErrList,[GoodPairings(:,3) pos_ref(GoodPairings(:,1),[1 2]) - pos_loc(GoodPairings(:,2),[1 2])]);

    TruePos_detections  = MatPairings(:,3)<Threshold_TruePos;
    T_pos(ifr)          = nnz(TruePos_detections);  % True positive points: good points
    F_neg(ifr)          = Nb_ref-T_pos(ifr);        % False negative: missing points
    F_pos(ifr)          = Nb_loc - T_pos(ifr);
%     F_pos(ifr)          = sum(~ismember(1:Nb_loc,GoodPairings(:,2))); % old)

    Npos_in(ifr)        = Nb_ref;
    Npos_loc(ifr)       = Nb_loc;

end
FinalPairs = FinalPairs';
MissingPoint = MissingPoint';
WrongLoc = WrongLoc';
Stat_classification = cat(2,Npos_in',Npos_loc',T_pos',F_neg',F_pos');

return
