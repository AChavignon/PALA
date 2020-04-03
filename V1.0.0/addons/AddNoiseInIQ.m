function IQ_speckle = AddNoiseInIQ(IQ,NoiseParam)
%% function IQ_speckle = AddNoiseInIQ(IQ,NoiseParam)
% Takes raw IQ in input and returns a noised IQ simulating a clutter noise.
% requires Communications Toolbox.
%
% Created by Arthur Chavignon 25/02/2020

% Noise Param example
%       - NoiseParam.Power        = -2; %[dBW]
%       - NoiseParam.Impedance    = .2; %[ohms]
%       - NoiseParam.SigmaGauss   = 1.5; % gaussian filtering
%       - NoiseParam.clutterdB    = -20; % clutter level in dB
%       - NoiseParam.amplCullerdB = 10; % dB amplitude of clutter

IQ = abs(IQ);
IQ_speckle = IQ+imgaussfilt(max(IQ(:))*10^(NoiseParam.clutterdB/20)+reshape(wgn(numel(IQ),1,NoiseParam.Power,NoiseParam.Impedance),size(IQ,1),size(IQ,2),[])*max(IQ(:))*10^((NoiseParam.amplCullerdB+NoiseParam.clutterdB)/20),NoiseParam.SigmaGauss);

end