function [Stat_classification,ErrList,FinalPairs,MissingPoint,WrongLoc] = PairingAlgorithm(ListPos_ref,MatTrackedLoc,PData,Threshold_pairing)
%% function [Stat_classification,ErrList,FinalPairs,MissingPoint,WrongLoc] = PairingAlgorithm(ListPos_ref,MatTrackedLoc,PData,Threshold_pairing)
% This function paired the detected bubbles (MatTrackedLoc) with the reel simulated
% bubbles (ListPos_ref) by calculating the interdistance matrix for each frame. Bubbles
% are paired when the RMSE in less than Threshold_pairing.
% INPUTS: 
%       - ListPos_ref: list of target simulated bubbles [x y z reflectivity]
%       - MatTrackedLoc: list of detected bubbles
%       - PData: parameters of IQ.
%       - Threshold_pairing: max distance between reel and detected point to consider a good localization. (\lambda/4)
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
% DATE 2020.03.30 - VERSION 1.0.0
% AUTHROS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris  
% Code Available under Creative Commons Non-Commercial 4.0
% ACADEMIC REFERENCES TO BE CITED
% Details of the code published in 2020 article by Heiles, Chavignon, Hingot and Couture.
% Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018 

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
    
    %% change the colums organisation [z x y reflectivity]
    pos_ref = pos_ref(:,[3 1 2 4]);  
    Points_ref{ifr} = pos_ref;
    N_ref = size(pos_ref,1); %number of target points

    %% get localized points after ULM
    idframe = find(MatTrackedLoc(:,3)==ifr);
    pos_loc = MatTrackedLoc(idframe,:); clear idframe

    Points_loc{ifr} = pos_loc;

    N_loc =  size(pos_loc,1); % number of localized points

    if 0
        % displaying
        lx = PData.Origin(1) + [0:PData.Size(2)-1]*PData.PDelta(1);
        lz = PData.Origin(3) + [0:PData.Size(1)-1]*PData.PDelta(3);
        figure(1);clf
        imagesc(lx,lz,20*log10(abs(temp.IQ(:,:,ifr))))
        hold on
        plot(pos_ref(:,2),pos_ref(:,1),'rx');
        plot(pos_loc(:,2),pos_loc(:,1),'wo');
    end

    %% Compute interdistance matrix
    DD_mat = zeros(N_ref,N_loc);
    for i_0=1:N_ref
        DD_mat(i_0,:) = vecnorm(pos_ref(i_0,[1 2])-pos_loc(:,[1 2]),2,2);  
    end
    % Sorting and pairing
    DD_mat_0 = DD_mat;

    [minErr_loc,~] = min(DD_mat_0,[],1);
    [~,indin_loc] = sort(minErr_loc(:));
    DD_mat_0 = DD_mat_0(:,indin_loc);

    [minErr_ref,~] = min(DD_mat_0,[],2);
    [~,indin_ref] = sort(minErr_ref(:));
    DD_mat_0 = DD_mat_0(indin_ref,:);

    DD_mat_00 = DD_mat_0;
    imin_last = 0;

    Pairings = zeros(0,3); % pairing matrix [indice_ref,indice_loc,distance]
    for ii=1:size(DD_mat_0,2)
        if isempty(DD_mat_00)
           break 
        end
        [err,imin] = min(DD_mat_00(:,1));
        Pairings = cat(1,Pairings,[indin_ref(imin_last+imin),indin_loc(ii),err]);
        DD_mat_00 = DD_mat_00((imin+1):end,2:end);
        imin_last = imin_last+imin;
    end

    %% Store data in outputs
    GoodPairings        = Pairings(Pairings(:,3)<Threshold_pairing,:);

    FinalPairs{ifr}     = GoodPairings;     % [index ref, index loc, RMSE]
    MissingPoint{ifr}   = setdiff(1:N_ref,GoodPairings(:,1)); % [index ref]
    WrongLoc{ifr}       = setdiff(1:N_loc,GoodPairings(:,2)); % [index loc]

    T_pos(ifr)          = size(GoodPairings,1);                      % true positive
    F_neg(ifr)          = N_ref-T_pos(ifr);
    F_pos(ifr)          = sum(~ismember(1:N_loc,GoodPairings(:,2))); % false positive (wrong loc)
    
    Npos_in(ifr)        = N_ref;
    Npos_loc(ifr)       = N_loc;
    
    ErrList             = cat(1,ErrList,[GoodPairings(:,3) pos_ref(GoodPairings(:,1),[1 2]) - pos_loc(GoodPairings(:,2),[1 2])]);
    
end
FinalPairs = FinalPairs';
MissingPoint = MissingPoint';
WrongLoc = WrongLoc';

Stat_classification = cat(2,Npos_in',Npos_loc',T_pos',F_neg',F_pos');

end






