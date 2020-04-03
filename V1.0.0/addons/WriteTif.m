function WriteTif(MatIn,OutPutColorMap,tif_filename,varargin)
%% function WriteTif(MatIn,OutPutColorMap,tif_filename)
% created by Arthur Chavignon 18/12/19
% example WriteTiff(MatOutConv(:,:,1:10:end),hot(256),'tiftest.tif')

if ~strcmp(tif_filename(end-3:end),'.tif')
    error('OutPutName should end by .tif ')
end

tmp = strcmpi(varargin,'overwrite'); 
if any(tmp),ForceOverWirte= varargin{find(tmp)+1};else, ForceOverWirte=0;end

if exist(tif_filename,'file')==2
    if ForceOverWirte==1
        choice = 'Overwrite';
    else
        % Ask the user if (s)he wants to overwrite the existing file:
        choice = questdlg(['The file  ',tif_filename,' already exists. Overwrite it?'], ...
            'The file already exists.','Overwrite','Cancel','Cancel');
    end
    % Overwriting basically means deleting and starting from scratch:
    
    if strcmp(choice,'Overwrite')
        delete(tif_filename)
    else
        clear tif_filename firstframe DelayTime DitherOption LoopCount frame
        error('The tiffing has been canceled.')
    end
end

c0 = [0 max(MatIn(:))];
tmp = strcmpi(varargin,'caxis'); 
if any(tmp)
    InputCaxis = varargin{find(tmp)+1}; 
    if numel(InputCaxis)==1
        c0 = c0.*InputCaxis; 
    elseif numel(InputCaxis)==2
        c0 = InputCaxis; 
    end
end

%% Convert
Ndigit = max(size(OutPutColorMap));

MatInSat = MatIn;
MatInSat(MatIn>c0(2)) = c0(2);
MatInSat(MatInSat<c0(1)) = c0(1);
MatInNorm = [MatInSat-c0(1)]/diff(c0);

MatInNorm = round(MatInNorm*Ndigit);
MatInNorm(MatInNorm>Ndigit)=Ndigit;

%% Write data

for ii=1:length(MatIn(1, 1, :))
    RGB = ind2rgb(MatInNorm(:,:,ii),OutPutColorMap); % convert to RGB
    imwrite(RGB, tif_filename, 'WriteMode', 'append');
end
% disp('Done')

return