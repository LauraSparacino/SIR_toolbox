%% theoretical simulation with 4 processes, showing coexistence of synergy and redundancy

clear; close all; clc
addpath([pwd '\..\functions\'])

%%% parameters
fs=100; % sampling frequency
nfft=1001; % number of points on frequency axis (total)
DimFont=16;
cpar=0.5; %% fixed coupling parameter to allow the presence of both syn and red

range_alpha=[8 12]; range_beta=[30 40];
nrange_alpha=round((nfft*2/fs)*range_alpha);
if range_beta(2) < fs/2
    nrange_beta=round((nfft*2/fs)*range_beta)+[1 0];
else
    nrange_beta(1)=round((nfft*2/fs)*range_beta(1))+1;
    nrange_beta(2)=round((nfft*2/fs)*fs/2);
end
nrange = [nrange_alpha; nrange_beta];

%% Simulation - Theoretical VAR process 
M=6; %%% number of processes
Mv=[2,2,2]; % structure of blocks
M1=length(Mv); % Note: Mv contains the dimension of each block
% iY=1; iX1=[2]; iX2=[3 4]; %% target and sources

%%% set poles and self-oscillations
par.poles{3}=([0.85 0.1]);
par.poles{4}=([0.95 0.35]);
par.poles{5}=([0.85 0.1]);
par.poles{6}=([0.95 0.2]);

%%% set the interactions as in Figure 1
c=cpar;
C_51=c;
C_31=1-c;
par.coup=[1 2 2 0.3;5 1 1 C_51;3 1 2 C_31;4 5 1 0.5;4 3 3 0.5;5 6 2 0.3];
par.Su=ones(1,M); %variance of innovation processes

[Am,Su]=theoreticalVAR(M,par); %% VAR parameters

% cross-spectral analysis: compute spectral matrices from VAR parameters
[S,~,f] = sir_VARspectra(Am',Su,nfft,fs);

%% compute information-theoretic measures: ER, MIR
% 1) for each process, compute the (spectral) entropy rate (ER)
% 2) for each pair of processes, compute the (spectral) mutual information rate (MIR)
i=nan(M1,M1,nfft); I=nan(M1,M1);
I_band=nan(M1,M1,size(nrange,1));

% i_1, i_2: the pairs of blocks to analyze
% i_1: index of first (vector) process Z1
% i_2: index of second (vector) process Z2
for i_1=1:M1 % 1st block
    for i_2=1:M1 % 2nd block
        out = sir_mir(S,Mv,i_1,i_2,nrange);
        if i_1 ~= i_2 % MIR outside the main diagonal
            i(i_1,i_2,:)=out.i12;  % freq measures
            I(i_1,i_2)=out.I12; % time domain measures
            I_band(i_1,i_2,:)=out.I_band; 
        else % ER on the main diagonal
            i(i_1,i_2,:)=out.h1; % freq measures
            I(i_1,i_2)=out.H1; % time domain measures
            I_band(i_1,i_2,:)=out.H1_band;
        end
    end
end

%% compute information-theoretic measures: deltaOIR, OIR
% 3) for each subset of order N of processes (N = 3,...,M), compute the (spectral) OIR
out = sir_oir(S,Mv,nrange);
dO=out.dO; dOf=out.dOf;
OIR=out.OIR; OIRf=out.OIRf;
oir=out.OIR_band;

%% plot ER, MIR 
nMIR=size(nchoosek(1:M1,2),1);
nrow=(M1+nMIR)/M1;
nrow=nrow+1;

figure('Color','w','WindowState','maximized')
for ii=1:M1
    subplot(nrow,M1,ii) % ERs
    plot(f,squeeze(i(ii,ii,:)),'r','LineWidth',1.8); hold on; 
    xlim([0 fs/2]); xticks([]);
    ylim([0 1.2*max(squeeze(i(ii,ii,:)))])
    if ii~=1
        yticks([])
    end
    if mean(squeeze(i(ii,ii,:))) < 10^(-3)
        ylim([0 1]);
    end
    legend(horzcat('h_{X_{',num2str(ii),'}}'))
    legend box off
    set(gca,'FontSize',DimFont)
end
subpl=M1+1;
for ii=1:M1
    for jj=ii+1:M1
        subplot(nrow,M1,subpl) % MIRs
        plot(f,squeeze(i(ii,jj,:)),'k','LineWidth',1.8); hold on; 
        xlim([0 fs/2]);
        if subpl <= (nMIR)
            xticks([]);
        else
            xticks(0:10:fs/2);
            xlabel('f[Hz]');
        end
        % if subpl~=M+1; yticks([]); elseif subpl~=nMIR+1; yticks([]); end
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

hold on

%% plot OIR
% some parameters for barplots
ntot=0;
for in=3:M1
    n(in)=size(nchoosek(1:3,in),1); %#ok
    ntot=ntot+n(in);
end
x_time=1:5:ntot*4+ntot; x{1}=2:5:ntot*4+ntot; x{2}=3:5:ntot*4+ntot;
load('colmap.mat');
colrange=VRVmap(1:256/ntot:256,:);
colors{3}=colrange(1:size(nchoosek(1:M1,3),1),:);
comb=nchoosek(1:M1,3);

subplot(3,3,7)
tmp=OIRf{3}; % OIR of order 3
icomb=1;
plot(f,tmp{icomb},'Color',colors{3}(1,:), ...
    'LineWidth',1.8); hold on;
xlim([0 fs/2]); xlabel('f[Hz]');
line([0 fs/2],[0 0],'Color',[0.5 0.5 0.5],'LineStyle','--',...
    'HandleVisibility','off')
title('\nu_{X^3}')
set(gca,'FontSize',DimFont)

m=3;
icomb=1;
k=1;
hold on
% integrated ER, MIR and OIR
ER_1=[I(1,1),I_band(1,1,1),I_band(1,1,2)];
ER_2=[I(2,2),I_band(2,2,1),I_band(2,2,2)];
ER_3=[I(3,3),I_band(3,3,1),I_band(3,3,2)];
ER=[ER_1;ER_2;ER_3];

MIR12_I=[I(1,2),I_band(1,2,1),I_band(1,2,2)];
MIR13_I=[I(1,3),I_band(1,3,1),I_band(1,3,2)];
MIR23_I=[I(2,3),I_band(2,3,1),I_band(2,3,2)];
MIR=[MIR12_I;MIR13_I;MIR23_I];

OIR_I=[OIR{3,1},oir{3,1},oir{3,2}];
TOT=[ER;MIR;OIR_I];

subplot(3,3,[8 9])
bar(TOT,'k');
xticks(1:7);
xticklabels({'H_{X_1}','H_{X_2}','H_{X_3}',...
    'I_{X_1;X_2}','I_{X_1;X_3}','I_{X_2;X_3}','\Omega_{X_1;X_2;X_3}'})
ylim([-0.2 4])
set(gca,'FontSize',DimFont)

