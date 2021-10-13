function [ErrList_out,FinalPairs_out,MissingPoint_out,WrongLoc_out,Stat_class] = PALA_Stat_multi(myfilepath,Track_in,listAlgo,PData,Threshold_pairing,Threshold_TruePos)
%% function [ErrList_out,FinalPairs_out,MissingPoint_out,WrongLoc_out,Stat_class] = PALA_Stat_multi(myfilepath,Track_in,listAlgo,PData,Threshold_pairing,Threshold_TruePos)
% This function runs the PairingAlgorithm for all blocs of IQ. It requires all tracks and
% needs to load raw 'ListPos' from simulation.
% INPUTS:
%       - myfilepath: filepath of simulated IQ that contains ListPos
%       - Track_in: list of tracks for each algorithm and bloc of IQ. (non interpolated tracks !)
%       - listAlgo: list of used ULM algorithms
%       - PData: parameters of IQ.
%       - Threshold_pairing: max distance between reel and detected point to consider a pair assignee. (\lambda/2)
%       - Threshold_TruePos: max distance between reel and detected point to consider a true positive localization. (\lambda/4)
% OUTPUTS:
%       - ErrList_out: for each algorithm, list of localization [dz,dx,norm(dz,dx)]
%       - FinalPairs_out: for each algorithm, list of pairs (index in ref, index indetection, RMSE)
%       - MissingPoint_out: for each algorithm, index of missing points in ListPos
%       - WrongLoc_out: for each algorithm, index of wrong detections in MatTrackedLoc_algo
%       - Stat_class: [Npos_in x Npos_loc x T_pos x F_neg x F_pos x Jaccard]
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

ErrList = {}; FinalPairs = {}; MissingPoint = {};WrongLoc = {}; Dist2Mean = {};

parfor hhh = 1:size(Track_in,1) % can be replace by a for loop.
    %% For each bloc of data
    fprintf([num2str(hhh) ' '])
    temp = load([myfilepath '_IQ' num2str(hhh,'%03.0f') '.mat'],'ListPos');

    ErrList_i = {}; FinalPairs_i = {}; MissingPoint_i = {}; WrongLoc_i = {};

    for ialgo = 1:numel(listAlgo)
        MatTrackedLoc_algo = cell2mat(Track_in{hhh,ialgo}); % [z,x,idframe] in lambda
        [Stat_class_i,ErrList_i{ialgo},FinalPairs_i{ialgo},MissingPoint_i{ialgo},WrongLoc_i{ialgo}] = PALA_PairingAlgorithm(temp.ListPos,MatTrackedLoc_algo,PData,Threshold_pairing,Threshold_TruePos);
        % Stat_classification : [Npos_in x Npos_loc x T_pos x F_neg x F_pos]
        Stat_classification{hhh}(:,ialgo) = sum(Stat_class_i,1);

    end
    ErrList{hhh} = ErrList_i;
    FinalPairs{hhh} = FinalPairs_i;
    MissingPoint{hhh} = MissingPoint_i;
    WrongLoc{hhh} = WrongLoc_i;
end
fprintf('\n')

%% Compute Jaccard index
if 0 % pull over repetitions
    Stat_class = sum(cat(3,Stat_classification{:}),3);
    Ja = Stat_class(3,:)./(Stat_class(4,:)+Stat_class(5,:) + Stat_class(3,:));
    Stat_class(end+1,:)=Ja;
else
    Stat_class = cat(3,Stat_classification{:});
    Ja = Stat_class(3,:,:)./(Stat_class(4,:,:)+Stat_class(5,:,:) + Stat_class(3,:,:));
    Stat_class(end+1,:,:)=Ja;
end
%% Concatenate results
ErrList = cat(1,ErrList{:});
FinalPairs = cat(1,FinalPairs{:});
MissingPoint = cat(1,MissingPoint{:});
WrongLoc = cat(1,WrongLoc{:});

for ii=1:numel(listAlgo)
    ErrList_out{ii} = cat(1,ErrList{:,ii});
    FinalPairs_out{ii} = cat(1,FinalPairs{:,ii});
    MissingPoint_out{ii} = cat(1,MissingPoint{:,ii});
    WrongLoc_out{ii} = cat(1,WrongLoc{:,ii});
end

end

