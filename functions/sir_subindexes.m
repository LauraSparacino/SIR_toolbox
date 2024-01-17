%% Extraction of indexes of the two blocks to analyze

%%% INPUT
% Mv: vector with the number of processes in each block
% i_1, i_2: vector indexes of the variables to put in the two blocks to analyze

%%% OUTPUT
% i1, i2: vector indexes of the columns to select for the two blocks (before pruning)
% j1, j2: vector indexes of the columns of the two blocks after pruning

function [i1,i2,j1,j2]=sir_subindexes(Mv,i_1,i_2)

M=length(Mv);
Q=sum(Mv);
iQ=1:Q; im=cell(1,M);
for m=1:M
    istart=sum(Mv(1:m)); iend=istart-Mv(m)+1;
    im{m}=iQ(iend:istart);
end
i1=[]; i2=[];
for i=1:length(i_1), i1=[i1 im{i_1(i)}]; end % columns of 1st block
for i=1:length(i_2), i2=[i2 im{i_2(i)}]; end % columns of 2nd block

% indexes of the two blocks to analyze
Mr1=length(i1); Mr2=length(i2); Mr=Mr1+Mr2;
j1=1:Mr1;
j2=Mr1+1:Mr;

end