function [ListAlgoName,ListColor,ListMarker,ListShortName] = PALA_GetFormat
%% function [ListAlgoName,ListColor,ListMarker,ListShortName] = PALA_GetFormat
%
% DATE 2020.07.22 - VERSION 1.1
% Created by Arthur Chavignon for ULM PALA 20/03/2020

ListAlgoName = {'No Shift','WA','Cub-Interp','Lz-Interp','Sp-Interp','Gauss-Fit','RS'};
ListShortName = {'NS','WA','CI','LI','SI','GF','RS'};

if 0
    ListColor = linspecer(8,'qualitative');
else
    ListColor =[0.9047    0.1918    0.1988;...
        0.2941    0.5447    0.7494;...
        0.3718    0.7176    0.3612;...
        1.0000    0.5482    0.1000;...
        0.8650    0.8110    0.4330;...
        0.6859    0.4035    0.2412;...
        0.9718    0.5553    0.7741;...
        0.6400    0.6400    0.6400];
end

ListColor(7,:)=[];
ListColor([1 7],:)=ListColor([7 1],:);
% ListColor = ListColor(1:8,:);
ListMarker = {'.','o','s','^','d','p','h'};

return

%% RGB color for Adobe Illustrator
% No shit :          163   163   163
% Weighed Average :   75   139   191
% Interp Cubic :      95   183    92
% Interp Lanczos :   255   140    25
% Interp Spline :    221   207   110
% Gaussian Fitting : 175   103    61
% Radial :           231    49    51

% print(gcf,[figname '.eps'],'-depsc2','-tiff','-r300','-painters')