%% APPLICATION TO PHYSIOLOGICAL NETWORK OF CARDIOVASCULAR, RESPIRATORY & CEREBROVASCULAR TIME SERIES
close all; clear; clc;
addpath([pwd '\..\functions\']);
addpath([pwd '\data\']);

do_surr=1; nsurr=100; 

%% load data
load('data_physio.mat'); % HP, SAP, MAP, CBFV, RESP
Yo=Yo'; 
[N,Q]=size(Yo);
Mv=ones(1,Q); % blocks of data
M=length(Mv);

Yo(:,1)=Yo(:,1)*1000; % sec to msec

names={'T', 'S', 'M', 'F', 'R'}; % HP, SAP, MAP, CBFV, RESP

%% parameters
pmax=12; pfilter=0.94;
nfft=1000;
chunklength=round(N/5);
bw=20; % number of truncation point for cross-correlation analysis
winname='parzenwin'; % 'rectwin' for rectangular, 'bartlett' for triangular, or window name ('blackman','hann','hamming',...)
typecorrest='biased'; % 'unbiased' for BT estimator, 'biased' for JW estimator
rangeLF=[0.04 0.15]; rangeHF=[0.15 0.4];
alpha=0.05;
DimFont=14;

%% analysis
fs=1000/mean(Yo(:,1)); % sampling frequency

% compute bands in frequency bins
nrangeLF=round((nfft*2/fs)*rangeLF);
if rangeHF(2) < fs/2
    nrangeHF=round((nfft*2/fs)*rangeHF)+[1 0];
else
    nrangeHF(1)=round((nfft*2/fs)*rangeHF(1))+1;
    nrangeHF(2)=round((nfft*2/fs)*fs/2);
end
nrange=[nrangeLF; nrangeHF];

% pre-processing: AR filtering and removal of mean value
Yf=nan*ones(N,Q); Y=nan*ones(N,Q);
for m=1:Q
    Yf(:,m)=AR_filter(Yo,m,pfilter); % AR highpass filtered series
    Y(:,m)=Yo(:,m)-mean(Yo(:,m)); % zero-mean series
    % Y(:,m)=zscore(Yo(:,m));
end

%%%% NON-PARAMETRIC APPROACH %%%%
[S,f] = sir_WCspectra(Y',bw,typecorrest,winname,fs,nfft);

%% MIR
I=zeros(M,M);
I_band=zeros(M,M,size(nrange,1)); 
i=zeros(M,M,nfft);

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

%% OIR
out = sir_oir(S,Mv,nrange);
dO=out.dO; dOf=out.dOf; % delta OIR
OIR=out.OIR; OIRf=out.OIRf; % OIR (time domain, spectral profiles)
oir=out.OIR_band; % OIR integrated in freq band

%% surrogates

if do_surr == 1

    % init
    I_surr=nan(nsurr,M,M); % time domain
    i_surr=nan(nsurr,M,M,size(nrange,1)); % frequency domain
    i_surr_f=nan(nsurr,M,M,nfft); % spectral profiles
    I_perc_f=cell(M,M);
    OIR_sign=OIR; % save new variables
    oir_sign=oir;

    %% ER, MIR statistical significance
    for ns=1:nsurr
    
        %%%%%%%%%%%%%%%%%%%%% statistical validation of ER 
        for q=1:Q
            data_surr_ER(q,:) = sir_surrshuf(Y(:,q)); %#ok
        end
        % non parametric spectral analysis
        S_est = sir_WCspectra(data_surr_ER,bw,typecorrest,winname,fs,nfft);
        % compute ER
        for i_1=1:M % 1st block
            for i_2=1:M % 2nd block
                out = sir_mir(S_est,Mv,i_1,i_2,nrange);
                if i_1 == i_2 % ER on the main diagonal
                    i_surr(ns,i_1,i_2,:)=out.H1_band; % freq measures
                    I_surr(ns,i_1,i_2)=out.H1; % time domain measures
                    i_surr_f(ns,i_1,i_2,:)=out.h1; % spectral profiles
                end
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%% statistical validation of MIR 
        for q=1:Q
            data_surr_MIR(q,:) = sir_surriaafft(Y(:,q)); %#ok
        end
        % non parametric spectral analysis
        S_est = sir_WCspectra(data_surr_MIR,bw,typecorrest,winname,fs,nfft);
        % compute MIR 
        for i_1=1:M % 1st block
            for i_2=1:M % 2nd block
                out = sir_mir(S_est,Mv,i_1,i_2,nrange);
                if i_1 ~= i_2 % MIR outside the main diagonal
                    i_surr(ns,i_1,i_2,:)=out.I_band;  % freq measures
                    I_surr(ns,i_1,i_2)=out.I12; % time domain measures
                    i_surr_f(ns,i_1,i_2,:)=out.i12; % spectral profiles
                end
            end
        end
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% test significance ER, MIR     
    for i_1=1:M % 1st block
        for i_2=1:M % 2nd block
    
            %%%%%%%%%%%%%%%% spectral profiles
            if i_1 ~= i_2 % MIR
                I_perc_f{i_1,i_2}=prctile(squeeze(i_surr_f(:,i_1,i_2,:)),...
                    [50 100-alpha*100],1); % < 95° perc
            else % ER
                I_perc_f{i_1,i_2}=prctile(squeeze(i_surr_f(:,i_1,i_2,:)),...
                    [alpha/2*100 50 100-alpha/2*100],1); % > 5° perc
            end
    
            %%%%%%%%%%%%%%%% time domain
            if i_1 ~= i_2 % MIR
                I_perc=prctile(I_surr(:,i_1,i_2),100-alpha*100);
                if I(i_1,i_2) < I_perc % non significant (< 95° perc)
                    I(i_1,i_2) = nan;
                end
            else % ER
                I_perc=prctile(I_surr(:,i_1,i_2),alpha*100);
                if I(i_1,i_2) > I_perc % non significant (> 5° perc)
                    I(i_1,i_2) = nan;
                end
            end
    
            %%%%%%%%%%%%%%%% frequency domain
            for irange=1:size(nrange,1)
                if i_1 ~= i_2 % MIR
                    i_perc=prctile(i_surr(:,i_1,i_2,irange),100-alpha*100);
                    if I_band(i_1,i_2,irange) < i_perc % non significant (< 95° perc)
                        I_band(i_1,i_2,irange) = nan;
                    end
                else % ER
                    i_perc=prctile(i_surr(:,i_1,i_2,irange),[alpha/2*100 100-alpha/2*100]);
                    if I_band(i_1,i_2,irange) < i_perc(1) % non significant (< 2.5°)
                        I_band(i_1,i_2,irange) = nan;
                    end
                end
            end
    
        end
    end

    %% OIR statistical significance
    % build bootstrap time series
    data_boot = sir_block_bootstrap(Y,chunklength,nsurr);

    for ns=1:nsurr

        Yb=data_boot{ns}';

        % spectral analysis
        [S,f] = sir_WCspectra(Yb,bw,typecorrest,winname,fs,nfft);
    
        % compute information-theoretic measures: deltaOIR, OIR
        out = sir_oir(S,Mv,nrange);
        for m=3:M
            tmp = out.OIRf(m); % spectral OIR
            for icomb=1:size(nchoosek(1:M,m),1)
                tmp_OIR_b=cell2mat(out.OIR(m));
                OIR_b{m}{icomb}(ns)=tmp_OIR_b(icomb); %#ok
                for irange=1:size(nrange,1)
                    tmp_oir_b=cell2mat(out.OIR_band(m,irange));
                    oir_b{irange}{m,icomb}(ns)=tmp_oir_b(icomb); %#ok
                end
                OIRf_b{m,icomb}(ns,:) = tmp{1}{icomb}; %#ok
            end
        end
       
    end
    
    % compute percentiles and check significance
    % OIR is not significant if the bootstrap distribution comprises the
    % zero value
    for m=3:M
        for icomb=1:size(nchoosek(1:M,m),1)
            OIR_perc_f{m,icomb}=prctile(OIRf_b{m,icomb},[alpha/2*100 50 100-alpha/2*100],1); %#ok
            % time domain
            OIR_perc=prctile(OIR_b{m}{icomb},[alpha/2*100 50 100-alpha/2*100]); 
            if OIR_perc(3) > 0 && OIR_perc(1) < 0
                OIR_sign{m}(icomb)=nan;
            end
            % integrated measures
            for irange=1:size(nrange,1)
                oir_perc=prctile(oir_b{irange}{m,icomb},[alpha/2*100 50 100-alpha/2*100]); 
                if oir_perc(3) > 0 && oir_perc(1) < 0
                    oir_sign{m,irange}(icomb)=nan;
                end
            end
        end
    end

end

%% plot representative ER, MIR
nMIR=size(nchoosek(1:M,2),1);
nrow=(Q+nMIR)/Q;
figure('Color','w','WindowState','maximized')
for ii=1:M
    subplot(nrow+1,M,ii) % ERs
    if do_surr==1
        area(f,I_perc_f{ii,ii}(3,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8],...
            'HandleVisibility','off'); hold on;
        area(f,I_perc_f{ii,ii}(1,:),'FaceColor',[1 1 1],'EdgeColor',[1 1 1],...
            'HandleVisibility','off');
        plot(f,I_perc_f{ii,ii}(1,:),'Color',[0.2 0.2 0.2],'LineWidth',1.2,...
            'HandleVisibility','off'); hold on; % 2.5° prctile
        plot(f,I_perc_f{ii,ii}(3,:),'Color',[0.2 0.2 0.2],'LineWidth',1.2,...
            'HandleVisibility','off'); hold on; % 97.5° prctile
        plot(f,I_perc_f{ii,ii}(2,:),'Color',[0.5 0.5 0.5],'LineWidth',1.2,...
            'HandleVisibility','off'); hold on; % 50° prctile
    end
    plot(f,squeeze(i(ii,ii,:)),'r','LineWidth',2); hold on; % original
    xlim([0 fs/2]); xticks([]);
    if squeeze(i(ii,ii,:)) < 0; ylim([1.1*min(squeeze(i(ii,ii,:))) 0]);
    else; ylim([0 1.1*max(squeeze(i(ii,ii,:)))]); end
    legend(horzcat('h_{',names{ii},'}'))
    legend box off
    set(gca,'FontSize',DimFont)
end
subpl=M+1;
for ii=1:M
    for jj=ii+1:M
        subplot(nrow+1,M,subpl) % MIRs
        ylabel('[nats]')
        if do_surr==1
            area(f,I_perc_f{ii,jj}(2,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8],...
                'HandleVisibility','off'); hold on;
            plot(f,I_perc_f{ii,jj}(1,:),'Color',[0.5 0.5 0.5],'LineWidth',1.2,...
                'HandleVisibility','off'); hold on; % 50° prctile
            plot(f,I_perc_f{ii,jj}(2,:),'Color',[0.2 0.2 0.2],'LineWidth',1.2,...
                'HandleVisibility','off'); hold on; % 95° prctile
        end
        plot(f,squeeze(i(ii,jj,:)),'k','LineWidth',2); hold on; % original
        xlim([0 fs/2]);
        if subpl <= (nMIR); xticks([]);
        else
            xticks(0:0.1:fs/2);
            xlabel('f[Hz]');
        end
        if subpl~=M+1 && subpl~=2*M+1
            yticks([])
        end
        ylabel('[nats/Hz]')
        ylim([0 1.5])
        legend(horzcat('i_{',names{ii},';',names{jj},'}'))
        legend box off
        set(gca,'FontSize',DimFont)
        subpl=subpl+1;
    end
end

% plot bars for ER, MIR
n_ER=Q;
x_time=1:5:n_ER*4+n_ER; x{1}=2:5:n_ER*4+n_ER; x{2}=3:5:n_ER*4+n_ER;
subplot(nrow+1,M,[16 17])
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
xticks(x{1}); xticklabels({'H_{T}','H_{S}','H_{M}','H_{F}','H_{R}'});
set(gca,'FontSize',DimFont)
line([0 fs/2],[0 0],'Color','k','LineStyle','-','HandleVisibility','off')

n_MIR=size(nchoosek(1:Q,2),1); 
x_time=1:5:n_MIR*4+n_MIR; x{1}=2:5:n_MIR*4+n_MIR; x{2}=3:5:n_MIR*4+n_MIR;
subplot(nrow+1,M,[19 20])
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
xticks(x{1}); xticklabels({'I_{T;S}','I_{T;M}','I_{T;F}','I_{T;R}',...
    'I_{S;M}','I_{S;F}','I_{S;R}','I_{M;F}','I_{M;R}','I_{F;R}'});
set(gca,'FontSize',DimFont)
line([0 fs/2],[0 0],'Color','k','LineStyle','-','HandleVisibility','off')

%% plot representative OIR
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
kk=1;
for m=3:M
    tmp=OIRf{m}; % OIR of order m
    comb=nchoosek(1:M,m);
    ncomb=size(comb,1);
    for icomb=1:ncomb
        subplot(3,6,kk)
        if do_surr==1
            area(f,OIR_perc_f{m,icomb}(3,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8],...
                'HandleVisibility','off'); hold on;
            if OIR_perc_f{m,icomb}(1,:) < 0
                area(f,OIR_perc_f{m,icomb}(1,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8],...
                    'HandleVisibility','off'); hold on;
            else
                for n=1:nfft
                    if OIR_perc_f{m,icomb}(1,n) > 0
                        line([f(n) f(n)],[0 OIR_perc_f{m,icomb}(1,n)],...
                            'Color',[1 1 1],'HandleVisibility','off'); hold on;
                    else
                        line([f(n) f(n)],[OIR_perc_f{m,icomb}(1,n) 0],...
                            'Color',[0.8 0.8 0.8],'HandleVisibility','off'); hold on;
                    end
                end
            end
            plot(f,OIR_perc_f{m,icomb}(1,:),'Color',[0.2 0.2 0.2],'LineWidth',1.2,...
                'HandleVisibility','off'); hold on; % 2.5° prctile
            plot(f,OIR_perc_f{m,icomb}(2,:),'Color',[0.5 0.5 0.5],'LineWidth',1.2,...
                'HandleVisibility','off'); hold on; % 50° prctile
            plot(f,OIR_perc_f{m,icomb}(3,:),'Color',[0.2 0.2 0.2],'LineWidth',1.2,...
                'HandleVisibility','off'); hold on; % 97.5° prctile
        end
        plot(f,tmp{icomb},'Color',colors{m}(icomb,:), ...
                'LineWidth',2); hold on; % original
        xlim([0 fs/2]); xlabel('f[Hz]');
        legend(horzcat('\nu^{',names{comb(icomb,:)},'}'))
        legend box off
        ylabel('[nats/Hz]')
        line([0 fs/2],[0 0],'Color',[0.5 0.5 0.5],'LineStyle','--',...
            'HandleVisibility','off')
        set(gca,'FontSize',DimFont)

        kk=kk+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%% time domain barplot
subplot(3,6,[17 18])
k=1;
for m=3:M
    comb=nchoosek(1:M,m);
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
title('\Omega')
ylabel('[nats]')
xticks(x{1});
xticklabels({'\Omega_{T;S;M}','\Omega_{T;S;F}','\Omega_{T;S;R}','\Omega_{T;M;F}','\Omega_{T;M;R}',...
    '\Omega_{T;F;R}','\Omega_{S;M;F}','\Omega_{S;M;R}','\Omega_{S;F;R}','\Omega_{M;F;R}',...
    '\Omega_{T;S;M;F}','\Omega_{T;S;M;R}','\Omega_{T;S;F;R}','\Omega_{T;M;F;R}','\Omega_{S;M;F;R}',...
    '\Omega_{T;S;M;F;R}'});
line([0 fs/2],[0 0],'Color','k','LineStyle','-','HandleVisibility','off')
set(gca,'FontSize',DimFont)

