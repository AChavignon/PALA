function IQf = SVDfilter(IQ,seuil)
%% function IQf = SVDfilter(IQ,seuil)
% SVD of a n-dim matrix
% the temporal dimension has to be the last dimension
% Arthur Chavignon 2019

initsize = size(IQ);
if seuil(end)>initsize(end)
    seuil = seuil(1):initsize(end);
end
if numel(seuil)==1
    seuil = seuil(1):initsize(end);
elseif numel(seuil)==2
    seuil = seuil(1):seuil(2);
end
    
if or(isequal(seuil,1:size(IQ,3)),seuil(1)<2)
    IQf = IQ;
    return
end

X = reshape(IQ,prod(initsize(1:end-1)),initsize(end)); % Reshape into Casorati matric
[U,~] = svd(X'*X);
% [U,~] = eig(X'*X);
V=X*U;
Reconst=V(:,seuil)*U(:,seuil)'; % Singular value decomposition

IQf = reshape(Reconst,initsize); % Reconstruction

end
