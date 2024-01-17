%% SIMULATION OF INTERACTING GAUSSIAN PROCESSES: STAR STRUCTURES
close all; clear; clc;
addpath([pwd '\..\functions\']);
warning('off');

% change these
spectrum='VAR'; % 'VAR': parametric spectral analysis,'noVAR': non parametric spectral analysis
sel_sim='out'; % 'in', 'out' star structure - switch simulation

%% theoretical simulation parameters
X=[1 2 3 4 5]; % processes
Q=length(X);
Mv=ones(1,Q); % blocks of processes
M=length(Mv); % Mv contains the dimension of each block

par.poles{X(1)}=([0.95 0.3]); % central hub 
for iout=2:Q
    par.poles{X(iout)}=([0.95 0.1]); % arms of the star
end

par.Su=ones(1,Q); % variance of innovation processes

% ingoing coefficients: c_in (from outside to inside)
% ongoing coefficients: c_out (from inside to outside)
c_in(2:Q)=0.5; c_out(2:Q)=0.5; 

switch sel_sim
    case 'in'
        par.coup=[X(2) X(1) 1 c_in(2)];
        for q=3:Q
            par.coup=[par.coup; ...
                X(q) X(1) 1 c_in(q)];
        end
    case 'out'
        par.coup=[X(1) X(2) 1 c_out(2)];
        for q=3:Q
            par.coup=[par.coup; ...
                X(1) X(q) 1 c_out(q)];
        end
end

[Am,Su]=theoreticalVAR(Q,par); %%% VAR parameters
Am=Am';

%% other parameters
Fs=1; % sampling frequency
nfft=1000; % number of points of freq axis
winname='parzenwin'; % 'rectwin' for rectangular, 'bartlett' for triangular, or window name ('blackman','hann','hamming',...)
typecorrest='biased'; % 'unbiased' for BT estimator, 'biased' for JW estimator
bw=20;
N=1000; % length of generated time series
pmax=8; % maximum scanned model order
alpha=0.01; % confidence interval
DimFont=14;

% frequency bands
rangeLF=[0.05 0.15]; rangeHF=[0.25 0.35];
nrangeLF=round((nfft*2/Fs)*rangeLF);
if rangeHF(2) < Fs/2
    nrangeHF=round((nfft*2/Fs)*rangeHF)+[1 0];
else
    nrangeHF(1)=round((nfft*2/Fs)*rangeHF(1))+1;
    nrangeHF(2)=round((nfft*2/Fs)*Fs/2);
end
nrange = [nrangeLF; nrangeHF];

%% from model parameters or data, compute the PSD matrix of the process
switch spectrum
    case 'VAR'
        %%%% PARAMETRIC APPROACH %%%%
        [S,~,f] = sir_VARspectra(Am,Su,nfft,Fs);
    case 'noVAR'
        %%%% NON-PARAMETRIC APPROACH %%%%
        % Generate time series
        U=mvnrnd(zeros(1,Q),Su,N);
        data=MVARfilter(Am,U');
        % compute WC spectrum
        [S,f] = sir_WCspectra(data,bw,typecorrest,winname,Fs,nfft);
end

%% compute information-theoretic measures: ER, MIR
% 1) for each process, compute the (spectral) entropy rate (ER)
% 2) for each pair of processes, compute the (spectral) mutual information rate (MIR)
i=nan(M,M,nfft); I=nan(M,M); I_band=nan(M,M,size(nrange,1));

% i_1, i_2: the pairs of blocks to analyze
% i_1: index of first (vector) process Z1
% i_2: index of second (vector) process Z2
for i_1=1:M % 1st block
    for i_2=1:M % 2nd block
        out = sir_mir(S,Mv,i_1,i_2,nrange);
        if i_1 ~= i_2 % MIR outside the main diagonal
            I(i_1,i_2)=out.I12; % time domain measures
            i(i_1,i_2,:)=out.i12;  % spectral profiles 
            I_band(i_1,i_2,:)=out.I_band; % integrated measures
        else % ER on the main diagonal
            I(i_1,i_2)=out.H1; % time domain measures
            i(i_1,i_2,:)=out.h1; % spectral profiles 
            I_band(i_1,i_2,:)=out.H1_band; % integrated measures
        end
    end
end

%% compute information-theoretic measures: deltaOIR, OIR
% 3) for each subset of order N of processes (N = 3,...,M), compute the (spectral) OIR

out = sir_oir(S,Mv,nrange);
dO=out.dO; dOf=out.dOf; % delta OIR
OIR=out.OIR;% time domain OIR
OIRf=out.OIRf; % spectral OIR
oir=out.OIR_band; % OIR integrated values in freq bands

%% plot ER, MIR 
nMIR=size(nchoosek(1:M,2),1);
nrow=(Q+nMIR)/Q;
figure('Color','w','WindowState','maximized')
for ii=1:M
    subplot(nrow+1,M,ii) % ERs
    plot(f,squeeze(i(ii,ii,:)),'r','LineWidth',1.8); hold on; 
    xlim([0 Fs/2]); xticks([]);
    ylim([0 1.2*max(squeeze(i(ii,ii,:)))])
    if ii~=1; yticks([]); end
    if mean(squeeze(i(ii,ii,:))) < 10^(-3)
        ylim([0 1]);
    end
    legend(horzcat('h_{X_{',num2str(ii),'}}'))
    legend box off
    set(gca,'FontSize',DimFont)
end
subpl=M+1;
for ii=1:M
    for jj=ii+1:M
        subplot(nrow+1,M,subpl) % MIRs
        plot(f,squeeze(i(ii,jj,:)),'k','LineWidth',1.8); hold on; 
        xlim([0 Fs/2]);
        if subpl <= (nMIR); xticks([]); else
            xticks(0:0.1:Fs/2);
            xlabel('f[Hz]');
        end
        ylim([0 1.2*max(squeeze(i(ii,jj,:)))])
        if mean(squeeze(i(ii,jj,:))) < 10^(-3)
            ylim([0 1]);
        end
        ylabel('[nats/Hz]')
        legend(horzcat('i_{X_{',num2str(ii),'};X_{',num2str(jj),'}}'))
        legend box off
        set(gca,'FontSize',DimFont)
        subpl=subpl+1;
    end
end

% plot bars for ER, MIR
n_ER=Q;
x_time=1:5:n_ER*4+n_ER; x{1}=2:5:n_ER*4+n_ER; x{2}=3:5:n_ER*4+n_ER;
subplot(nrow+1,M,[16 17]) % ERs
k=1;
for ii=1:M
    for jj=1:M
        if ii== jj
            bar(x_time(k),I(ii,jj),'r'); hold on;
            for irange=1:size(nrange,1)
                bar(x{irange}(k),I_band(ii,jj,irange),'r'); hold on;
            end
            k=k+1;
        end
    end
end
ylabel('[nats]')
xticks(x{1}); xticklabels({'H_{X_1}','H_{X_2}','H_{X_3}','H_{X_4}','H_{X_5}'});
set(gca,'FontSize',DimFont)

n_MIR=size(nchoosek(1:Q,2),1); 
x_time=1:5:n_MIR*4+n_MIR; x{1}=2:5:n_MIR*4+n_MIR; x{2}=3:5:n_MIR*4+n_MIR;
subplot(nrow+1,M,[19 20]) % MIRs
k=1;
for ii=1:M
    for jj=ii+1:M
        if ii ~= jj
            bar(x_time(k),I(ii,jj),'k'); hold on;
            for irange=1:size(nrange,1)
                bar(x{irange}(k),I_band(ii,jj,irange),'k'); hold on;
            end
            k=k+1;
        end
    end
end
ylabel('[nats]')
xticks(x{1}); xticklabels({'I_{X_1;X_2}','I_{X_1;X_3}','I_{X_1;X_4}','I_{X_1;X_5}',...
    'I_{X_2;X_3}','I_{X_2;X_4}','I_{X_2;X_5}','I_{X_3;X_4}','I_{X_3;X_5}','I_{X_4;X_5}'});
set(gca,'FontSize',DimFont)

%% plot OIR
% some parameters for barplots
ntot=0;
for in=3:M
    n(in)=size(nchoosek(1:Q,in),1); %#ok
    ntot=ntot+n(in);
end
x_time=1:5:ntot*4+ntot; x{1}=2:5:ntot*4+ntot; x{2}=3:5:ntot*4+ntot;

% colormap for multiplets (this is up to order 5)
load('colmap.mat');
colrange=VRVmap(1:256/ntot:256,:);
colors{3}=colrange(1:size(nchoosek(1:Q,3),1),:);
colors{4}=colrange(size(nchoosek(1:Q,3),1)+1:size(nchoosek(1:Q,3),1)+size(nchoosek(1:Q,4),1),:);
colors{5}=colrange(size(nchoosek(1:Q,3),1)+size(nchoosek(1:Q,4),1)+1:...
    size(nchoosek(1:Q,3),1)+size(nchoosek(1:Q,4),1)+size(nchoosek(1:Q,5),1),:);

%%%%%%%%%%%%%%%%%%%%%% spectral profiles
figure('Color','w','WindowState','maximized')
for m=3:M
    subplot(1,M-1,m-2)
    tmp=OIRf{m}; % OIR of order m
    comb=nchoosek(1:M,m);
    ncomb=size(comb,1);
    for icomb=1:ncomb
        plot(f,tmp{icomb},'Color',colors{m}(icomb,:), ...
            'LineWidth',1.8); hold on;
        legend(num2str(comb))
    end
    legend boxoff
    ylabel('[nats/Hz]')
    switch sel_sim
        case 'in'
            if m==M
                ylim([-6 0.05])
            else
                ylim([-1 0.05])
            end
        case 'out'
            if m==M
                ylim([0 5])
            else
                ylim([0 3])
            end
    end
    xlim([0 Fs/2]); xlabel('f[Hz]');
    title(horzcat('\nu_{X^',num2str(m),'}'))
    line([0 Fs/2],[0 0],'Color',[0.5 0.5 0.5],'LineStyle','--',...
        'HandleVisibility','off')
    set(gca,'FontSize',DimFont)
end

%%%%%%%%%%%%%%%%%%%%%% time domain barplot
subplot(1,M-1,4);
k=1;
for m=3:M
    for icomb=1:length(OIR{m})
        bar(x_time(k),OIR{m}(icomb),'EdgeColor',colors{m}(icomb,:),...
            'FaceColor',colors{m}(icomb,:)); hold on;
        for irange=1:size(nrange,1)
            bar(x{irange}(k),oir{m,irange}(icomb),'EdgeColor',colors{m}(icomb,:),...
            'FaceColor',colors{m}(icomb,:)); hold on;
        end
        k=k+1;
    end
end
xticks([]); yticks([]);
line([0 Fs/2],[0 0],'Color','k','LineStyle','-','HandleVisibility','off')
switch sel_sim
    case 'in'
        ylim([-0.8 0.05])
end
ylabel('[nats]')
xticks(x{1});
xticklabels({'\Omega_{1;2;3}','\Omega_{1;2;4}','\Omega_{1;2;5}','\Omega_{1;3;4}','\Omega_{1;3;5}',...
    '\Omega_{1;4;5}','\Omega_{2;3;4}','\Omega_{2;3;5}','\Omega_{2;4;5}','\Omega_{3;4;5}',...
    '\Omega_{1;2;3;4}','\Omega_{1;2;3;5}','\Omega_{1;2;4;5}','\Omega_{1;3;4;5}','\Omega_{2;3;4;5}',...
    '\Omega_{1;2;3;4;5}'});
line([0 Fs/2],[0 0],'Color','k','LineStyle','-','HandleVisibility','off')
title('\Omega')
set(gca,'FontSize',DimFont)
