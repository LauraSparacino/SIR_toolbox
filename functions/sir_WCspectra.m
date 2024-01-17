%% Performs non-parametric estimation of the Power Spectral Density of the process
% with Weighted covariance estimator

%%% INPUT
% X: data matrix M*N (N data points, M processes)
% M: number of truncation points for cross-correlation analysis 
% typecorrest: 'unbiased' for BT estimator, 'biased' for JW estimator
% winname: 'rectwin' for rectangular, 'bartlett' for triangular, or window name ('blackman','hann','hamming',...)
% fs: sampling freq
% nfft: number of points of frequency axis (positive frequencies)

%%% OUTPUT
% S: power spectral density matrix 
% f: frequency axis

function [S,f] = sir_WCspectra(X,m,typecorrest,winname,fs,nfft)

narginchk(1,6); % from 1 to 6 input arguments
if nargin < 6, nfft=1000; end 
if nargin < 5, fs=1; end 
if nargin < 4, winname='parzenwin'; end 
if nargin < 3, typecorrest='biased'; end 
if nargin < 2, m=round((1.273*Fs)/(2.5)); end 

X=X-mean(X,2); % always work with zero-mean data
Q=size(X,1);

S=zeros(Q,Q,nfft);

for i=1:Q

    rx=xcorr(X(i,:),m,typecorrest); % correlation estimate, truncated at lag m
    eval(strcat('wi=',winname,'(',int2str(2*m+1),');'));
    rxw=rx.*wi'; % truncated and windowed correlation
    Px=fft(rxw,2*nfft); % fft estimate
    Px=(2/fs)*Px(1:nfft); 
    f=(0:fs/(2*(nfft-1)):fs/2)'; % frequency axis
    
    for j=1:Q
        if i~=j
            rxy=xcorr(X(i,:),X(j,:),m,typecorrest); % correlation estimate, truncated at lag m
            eval(strcat('wi=',winname,'(',int2str(2*m+1),');'));
            rxyw=rxy.*wi'; % truncated and windowed correlation
            Pxy=fft(rxyw,2*nfft); % fft estimate
            Pxy=(2/fs)*Pxy(1:nfft); 
            S(i,j,:) = Pxy;
        end
    end

    S(i,i,:) = Px;

end

end