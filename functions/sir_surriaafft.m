%% Iterative Amplitude Adjusted Fourier Tranform (iAAFT) surrogates 
% Schreiber & Schmitz algorithm - Physical Review Letters 1996

%%% INPUT
% y: input series
% nit: number of desired iterations (default 7)
% stop: 'spe' to save spectrum, 'dis' to save distribution

%%% OUTPUT
% ys: surrogate series

function ys=sir_surriaafft(y,nit,stop)

% check inputs
error(nargchk(1,3,nargin)); % min, max input arguments
if nargin < 3, stop='spe'; end 
if nargin < 2, nit=7; end % default 7 iterations

[ysorted,~]=sort(y); % from the lowest to the highest value
my=abs(fft(y));
ys=sir_surrshuf(y); % shuffling

for i=1:nit
    % step 1: set the spectrum
    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+1i*sin(faseys));
    ys=ifft(fys); ys=real(ys);
    ys=ys-mean(ys);

    % step 2: set the distribution
    [~,ysindice]=sort(ys);
    ypermuted=zeros(length(y),1);
    for ii=1:length(y)
        ypermuted(ysindice(ii))=ysorted(ii);
    end
    ys=ypermuted;

end

if stop == 'spe'
    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+1i*sin(faseys));
    ys=ifft(fys);ys=real(ys);
    ys=ys-mean(ys);
end




