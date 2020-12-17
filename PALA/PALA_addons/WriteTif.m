function WriteTif(MatIn,OutPutColorMap,tif_filename,varargin)
%% function WriteTif(MatIn,OutPutColorMap,tif_filename,varargin)
% created by Arthur Chavignon 18/12/19
% example WriteTif(MatOutConv(:,:,1:10:end),hot(256),'tiftest.tif')

if ~strcmp(tif_filename(end-3:end),'.tif')
    error('OutPutName should end by .tif ')
end

tmp = strcmpi(varargin,'overwrite');
if any(tmp),ForceOverWirte= varargin{find(tmp)+1};else, ForceOverWirte=0;end

if exist(tif_filename,'file')==2
    if ForceOverWirte==1
        choice = 'Overwrite';
    else
        % Ask the user whether to overwrite the existing file:
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

%% Convert into digits
if ~isempty(OutPutColorMap)
    Ndigit = max(size(OutPutColorMap));
    
    MatInSat = MatIn;
    MatInSat(MatIn>c0(2)) = c0(2);
    MatInSat(MatInSat<c0(1)) = c0(1);
    MatInNorm = [MatInSat-c0(1)]/diff(c0);
    
    MatInNorm = round(MatInNorm*Ndigit);
    MatInNorm(MatInNorm>Ndigit)=Ndigit;
else
    % if MatIn is already RGB
    MatInNorm = MatIn;
end
% Write RGB data

for ii=1:size(MatInNorm,3)
    if ~isempty(OutPutColorMap)
        RGB = ind2rgb(MatInNorm(:,:,ii),OutPutColorMap); % convert to RGB
    else
        RGB = squeeze(MatInNorm(:,:,ii,:));
    end
    imwrite(RGB, tif_filename,'WriteMode','append','Compression','lzw','Description','ArthurChavignon');
end

tmp = strcmpi(varargin,'VoxelSizeMm');
if any(tmp),VoxelSizeMm= varargin{find(tmp)+1};
    t_tif = Tiff(tif_filename,'r+');
    VoxSize_cm = VoxelSizeMm*1e-3;
    setTag(t_tif,'XResolution',1/VoxSize_cm(1))
    setTag(t_tif,'YResolution',1/VoxSize_cm(2))
    setTag(t_tif,'ResolutionUnit',Tiff.ResolutionUnit.Centimeter)
    close(t_tif);
end

return