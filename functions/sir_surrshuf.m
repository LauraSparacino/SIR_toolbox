%% Generation of random shuffled (IID) surrogate data 

%%% INPUT
% x: original series

%%% OUTPUT
% xs: surrogate series

function xs = sir_surrshuf(x)

sx=size(x);
p=randperm(sx(1));
xs=zeros(sx(1),1);
for k = 1:sx(1)
	xs(k)=x(p(k));
end

end
