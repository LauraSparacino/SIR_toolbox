%% Model Order Selection for identification of the strictly causal VAR model

%%% INPUT
% Y, M*N matrix of time series (each time series is in a row)
% pmax, maximum tested model order
% idMode, determines estimation algorithm (0: builtin least squares, else other methods [see mvar.m from biosig package])

%%% OUTPUT
% pottaic: model order optimized with multichannel Akaike Information Criterion (AIC)
% pottaic: model order optimized with  Minimum Description Length (MDL) criterion
% aic: values of AIC index as a function of the model order
% mdl: values of MDL index as a function of the model order

function [pottaic,pottmdl,aic,mdl] = sir_mos_idMVAR(Y,pmax,idMode)

if nargin < 3; idMode=0; end % default: builtin least squares

[M,N]=size(Y); 

% figures of merit
aic=NaN*ones(pmax,1); mdl=aic;

for p=1:pmax
    
    [~,S,~]=sir_idMVAR(Y,p,idMode); % S is the covariance matrix

    % multivariate AIC 
    aic(p)=N*log(det(S))+2*M*M*p; 
    
    % multivariate MDL
    mdl(p)=N*log(det(S))+log(N)*M*M*p;
    
        
end

pottaic=find(aic == min(aic));
pottmdl=find(mdl == min(mdl));
