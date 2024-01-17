%% Computation of spectral and time-domain entropy rate (ER) and mutual information rate (MIR)
% for two blocks of processes given the spectral matrix

%%% INPUT
% S: power spectral density matrix
% Mv: vector with the number of processes in each block
% i_1: first block index 
% i_2: second block index
% nrange: matrix with rows=number of ranges, columns=2 (limits of the bands in terms of frequency bins)
% --- example: if nfft=1000, nrange=[80 300] 
% --- nrange can be specified or not

%%% OUTPUT
% h12: spectral joint ER
% h1, h2: spectral ERs
% i12: spectral MIR
% H12: time domain joint ER
% H1, H2: time domain ERs
% I12: time domain MIR
% H12_band: integrated joint ER in the selected spectral bands
% H1_band, H2_band: integrated ERs in the selected spectral bands
% I_band: integrated MIR in the selected spectral bands

function out = sir_mir(S,Mv,i_1,i_2,nrange)

% check inputs
if nargin < 4
    error('Not enough input arguments')
end

[i1,i2] = sir_subindexes(Mv,i_1,i_2);

d1=length(i1); d2=length(i2);
nfft=size(S,3);

h12=nan*ones(nfft,1); h1=nan*ones(nfft,1); h2=nan*ones(nfft,1);
i12=nan*ones(nfft,1);
for n=1:nfft
    h12(n)=0.5*log((2*pi*exp(1))^(d1+d2)*abs(det(S([i1 i2],[i1 i2],n)))); % spectral ER (X1,X2)
    h1(n)=0.5*log((2*pi*exp(1))^d1*abs(det(S(i1,i1,n)))); % spectral ER X1
    h2(n)=0.5*log((2*pi*exp(1))^d2*abs(det(S(i2,i2,n)))); % spectral ER X2
    i12(n)=h1(n)+h2(n)-h12(n); % spectral MIR = ER (X1) + ER (X2) - joint ER (X1,X2)
end

% time domain
H12=sum(h12)/nfft; 
H1=sum(h1)/nfft; 
H2=sum(h2)/nfft; 
I12=sum(i12)/nfft; 

% put nan on the main diagonal
if i_1 == i_2
    i12(:)=nan; I12=nan;
    h12(:)=nan; H12=nan;
end

%%% integral of ER, MIR inside spectral bands
if nargin == 5
    nr = size(nrange,1);
    I_band = zeros(nr,1); 
    H12_band = I_band; H1_band = I_band; H2_band = I_band; 
    for ir=1:nr
        I_band(ir,1) = sum(i12(nrange(ir,1):nrange(ir,2)))/nfft;
        H12_band(ir,1) = sum(h12(nrange(ir,1):nrange(ir,2)))/nfft;
        H1_band(ir,1) = sum(h1(nrange(ir,1):nrange(ir,2)))/nfft;
        H2_band(ir,1) = sum(h2(nrange(ir,1):nrange(ir,2)))/nfft;
    end
    % integrated values
    out.I_band=I_band; 
    out.H12_band=H12_band; 
    out.H1_band=H1_band; 
    out.H2_band=H2_band; 
end

%% OUTPUT
out.h12=h12;
out.h1=h1;
out.h2=h2;
out.i12=i12;
out.H12=H12;
out.H1=H1;
out.H2=H2;
out.I12=I12;

end