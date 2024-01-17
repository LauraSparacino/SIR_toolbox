%% Computation of the gradient of the O-information rate when the block X_N is added to the group X^(N-1)

%%% INPUT
% S: power spectral density matrix
% Mv: vector with the number of processes in each block
% ij: complete vector of indexes
% j: index of the target to analyze within ij

%%% OUTPUT
% dO: delta OIR
% dOf: spectral delta OIR

function out = sir_deltaO(S,Mv,ij,j)

nfft=size(S,3);

%%% computation of dO(X-j;Xj)
assert(ismember(j,ij)); % verify target belongs to group
ii=setdiff(ij,j); % driver indexes
N=length(ij); % order of the OIR to compute

out1=sir_mir(S,Mv,ii,j); % MIR between Xj and X-j

i_cs=nchoosek(ii,N-2); % number of combinations to compute the sum of the deltaO 
dO=0; dOf=zeros(nfft,1); % init
for cnt=1:N-1
    outtmp=sir_mir(S,Mv,i_cs(cnt,:),j); % MIR between Xj and X-ij
    % time domain measure
    dO=dO+outtmp.I12;
    % spectral function
    dOf=dOf+outtmp.i12;
end

% time domain measures: last term 
dO=dO+(2-N)*out1.I12; % dO(X-j;Xj)

% spectral functions: last term
dOf=dOf+(2-N)*out1.i12; 

%% output
out.dO=dO;
out.dOf=dOf;

end

