%% Application to EEG data (motor execution task)
close all; clear; clc
addpath([pwd '\data\']);
addpath([pwd '\..\functions\']);

load('EEG_data.mat')
data=EEG.samp;

%%% parameters
nfft=1024;
Fs=EEG.fsamp;
IP=10; % fixed model order

Mv=[2 2 2];
M=length(Mv);
Q=6; % number of processes
nsurr=100;
alpha=0.05;

range_alpha=[7 15]; range_beta=[18 26];
nrange_alpha=round((nfft*2/Fs)*range_alpha);
if range_beta(2) < Fs/2
    nrange_beta=round((nfft*2/Fs)*range_beta)+[1 0];
else
    nrangeHF(1)=round((nfft*2/Fs)*rangeHF(1))+1;
    nrangeHF(2)=round((nfft*2/Fs)*Fs/2);
end
nrange = [nrange_alpha; nrange_beta];

% init
i=zeros(M,M,nfft);
I=zeros(M,M);
I_band=zeros(M,M,size(nrange,1));
i_surr_f=zeros(nsurr,M,M,nfft);
I_surr=zeros(nsurr,M,M);
i_surr=zeros(nsurr,M,M,size(nrange,1));
I_perc_f=cell(M,M);

for uu=1:size(data,3) % trial

    tic
    disp(['Trial: ' num2str(uu)])

    % Model Identification 
    [Am,Su]=sir_idMVAR(data(:,:,uu),IP);

    % Estimation of the spectral matrix
    [S,~,f] = sir_VARspectra(Am,Su,nfft,Fs);
    
    % compute ER and MIR
    for i_1=1:M % 1st block
        for i_2=1:M % 2nd block
            out= sir_mir(S,Mv,i_1,i_2,nrange);
            if i_1 ~= i_2 % MIR outside the main diagonal
                i(i_1,i_2,:)=out.i12;  % freq measures
                I(i_1,i_2)=out.I12; % time domain measures
                I_band(i_1,i_2,:)=out.I_band; % integrated measures
            else % ER on the main diagonal
                i(i_1,i_2,:)=out.h1; % freq measures
                I(i_1,i_2)=out.H1; % time domain measures
                I_band(i_1,i_2,:)=out.H1_band; % integrated measures
            end
        end
    end

    % compute delta OIR and OIR
    clear out
    out= sir_oir(S,Mv,nrange);
    dO=out.dO; dOf=out.dOf;
    OIR_band=out.OIR_band;
    OIR=out.OIR; OIRf=out.OIRf;
    
    % save data in a struct
    RES.i{uu}=i;
    RES.I{uu}=I;
    RES.I_band{uu}=I_band;
    RES.OIR_band{uu}=OIR_band;
    RES.dO{uu}=dO;
    RES.dOf{uu}=dOf;
    RES.OIR{uu}=OIR;
    RES.OIRf{uu}=OIRf;
   
    % surrogate data analysis
    for ns=1:nsurr

        %%%%%%%%%%%%%%%%%%%%% statistical validation of ER
        for q=1:Q
            data_surr_ER(q,:) = sir_surrshuf(data(q,:,uu)'); %#ok
        end
        % parametric spectral analysis
        [Am_est,Su_est]=sir_idMVAR(data_surr_ER,IP);
        % Estimation of the spectral matrix
        S_est = sir_VARspectra(Am_est,Su_est,nfft,Fs);
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
            data_surr_MIR(q,:) = sir_surriaafft(data(q,:,uu)'); %#ok
        end
        % parametric spectral analysis
        [Am_est,Su_est]=sir_idMVAR(data_surr_MIR,IP);
        % Estimation of the spectral matrix
        [S_est,~,f] = sir_VARspectra(Am_est,Su_est,nfft,Fs);
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

    RES.i_surr{uu}=i_surr;
    RES.I_surr{uu}=I_surr;
    RES.i_surr_f{uu}=i_surr_f;
        
    %%%%%%%%%%%%%%%%%%%%%%% test significance ER, MIR     
    for i_1=1:M % 1st block
        for i_2=1:M % 2nd block

            %%%%%%%%%%%%%%%% spectral profiles
            if i_1 ~= i_2 % MIR
                I_perc_f{i_1,i_2}=prctile(squeeze(i_surr_f(:,i_1,i_2,:)),...
                    [50 100-alpha*100],1); 
            else % ER
                I_perc_f{i_1,i_2}=prctile(squeeze(i_surr_f(:,i_1,i_2,:)),...
                    [alpha/2*100 50 100-alpha/2*100],1); 
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

    RES.i_perc{uu}=i_perc;
    RES.I_perc{uu}=I_perc;
    RES.I_perc_f{uu}=I_perc_f;
    
    %% OIR statistical significance
    % build bootstrap time series
    chunklength=round(size(data(:,:,uu)',1)/5);
    data_boot = sir_block_bootstrap(data(:,:,uu)',chunklength,nsurr);
    
    for ns=1:nsurr

        Yb=data_boot{ns}';
        % parametric spectral analysis
        [Am_est,Su_est]=sir_idMVAR(Yb,IP);
        % Estimation of the spectral matrix
        [S_est,~,f] = sir_VARspectra(Am_est,Su_est,nfft,Fs);
        % compute information-theoretic measures: deltaOIR, OIR
        out = sir_oir(S_est,Mv,nrange);
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
    RES.OIRf_b{uu}=OIRf_b;

    % compute percentiles and check significance
    % OIR is not significant if the bootstrap distribution comprises the zero value
    for m=3:M
        for icomb=1:size(nchoosek(1:M,m),1)
            OIR_perc_f{m,icomb}=prctile(OIRf_b{m,icomb},[alpha/2*100 50 100-alpha/2*100],1); %#ok
            % time domain
            OIR_perc=prctile(OIR_b{m}{icomb},[alpha/2*100 50 100-alpha/2*100]); 
            if OIR_perc(3) > 0 && OIR_perc(1) < 0
                OIR_sign{m}(icomb)=nan; %#ok
            end
            % integrated measures
            for irange=1:size(nrange,1)
                oir_perc=prctile(oir_b{irange}{m,icomb},[alpha/2*100 50 100-alpha/2*100]); 
                if oir_perc(3) > 0 && oir_perc(1) < 0
                    oir_sign{m,irange}(icomb)=nan; %#ok
                end
            end
        end
    end
   
    RES.OIR_perc_f{uu}=OIR_perc_f;
    RES.OIR_perc{uu}=OIR_perc;
    RES.oir_perc{uu}=oir_perc;

    toc

end

clearvars -except RES f data

%% averaging over trials
for t=1:20
    ER_f(:,:,:,t)=RES.i{t}; %#ok
    ER_whole(:,:,t)=RES.I{t}; %#ok
    ER_band(:,:,:,t)=RES.I_band{t}; %#ok
    % OIR
    OIR_f(:,t)=RES.OIRf{t}{3,1}{1,1}; %#ok
    OIR_whole(t)=RES.OIR{t}{3,1}; %#ok
    OIR_band(1,1,t)=RES.OIR_band{t}{3,1};
    OIR_band(1,2,t)=RES.OIR_band{t}{3,2};
    % surrogates
    % f domain
    ER11_perc_f(:,:,t)=RES.I_perc_f{t}{1,1}; %#ok
    ER22_perc_f(:,:,t)=RES.I_perc_f{t}{2,2}; %#ok
    ER33_perc_f(:,:,t)=RES.I_perc_f{t}{3,3}; %#ok
    MIR12_perc_f(:,:,t)=RES.I_perc_f{t}{1,2}; %#ok
    MIR13_perc_f(:,:,t)=RES.I_perc_f{t}{1,3}; %#ok
    MIR23_perc_f(:,:,t)=RES.I_perc_f{t}{2,3}; %#ok
    % OIR surr
    OIR_f_perc(:,:,t)=RES.OIR_perc_f{t}{3,1}; %#ok
end

% averaging over trials
ER_f=squeeze(mean(ER_f,4));
ER_whole=squeeze(mean(ER_whole,3));
ER_band=squeeze(mean(ER_band,4));
OIR_f=squeeze(mean(OIR_f,2));
OIR_whole=squeeze(mean(OIR_whole,2));
OIR_alpha=squeeze(mean(OIR_band(1,1,:,:),3));
OIR_beta=squeeze(mean(OIR_band(1,2,:,:),3));
% surrogates 
ER11_perc=squeeze(mean(ER11_perc_f,3));
ER22_perc=squeeze(mean(ER22_perc_f,3));
ER33_perc=squeeze(mean(ER33_perc_f,3));
MIR12_perc=squeeze(mean(MIR12_perc_f,3));
MIR13_perc=squeeze(mean(MIR13_perc_f,3));
MIR23_perc=squeeze(mean(MIR23_perc_f,3));
% surrogate OIR_f
OIR_perc=squeeze(mean(OIR_f_perc,3));

I_perc{1,1}{1,1}=ER11_perc;
I_perc{1,1}{2,2}=ER22_perc;
I_perc{1,1}{3,3}=ER33_perc;
I_perc{1,1}{1,2}=MIR12_perc;
I_perc{1,1}{1,3}=MIR13_perc;
I_perc{1,1}{2,3}=MIR23_perc;
OIR_perc_T{1,1}=OIR_perc;

%% %% plot representative ER, MIR
DimFont=14;
do_surr=1;
M=3; Q=3;
nMIR=size(nchoosek(1:M,2),1);
nrow=(Q+nMIR)/Q;

figure('Color','w','WindowState','maximized')
for ii=1:M
    subplot(nrow+1,M,ii) % ERs
    if do_surr==1
        area(f,I_perc{1}{ii,ii}(3,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]); hold on;
        area(f,I_perc{1}{ii,ii}(1,:),'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
        plot(f,I_perc{1}{ii,ii}(1,:),'Color',[0.2 0.2 0.2],'LineWidth',1.2); hold on; % 2.5° prctile
        plot(f,I_perc{1}{ii,ii}(3,:),'Color',[0.2 0.2 0.2],'LineWidth',1.2); hold on; % 97.5° prctile
        plot(f,I_perc{1}{ii,ii}(2,:),'Color',[0.5 0.5 0.5],'LineWidth',1.2); hold on; % 50° prctile
    end
    plot(f,squeeze(ER_f(ii,ii,:)),'r','LineWidth',2); hold on; % original
    xlim([0 30]); xticks([]);
    if squeeze(ER_f(ii,ii,:)) < 0; ylim([1.1*min(squeeze(ER_f(ii,ii,:))) 0]);
    else; ylim([0 1.1*max(squeeze(ER_f(ii,ii,:)))]); end
    set(gca,'FontSize',DimFont)
end
subpl=M+1;
for ii=1:M
    for jj=ii+1:M
        subplot(nrow+1,M,subpl) % MIRs
        if do_surr==1
            area(f,I_perc{1}{ii,jj}(2,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]); hold on;
            plot(f,I_perc{1,1}{ii,jj}(1,:),'Color',[0.5 0.5 0.5],'LineWidth',1.2); hold on; % 50° prctile
            plot(f,I_perc{1,1}{ii,jj}(2,:),'Color',[0.2 0.2 0.2],'LineWidth',1.2); hold on; % 95° prctile
        end
        plot(f,squeeze(ER_f(ii,jj,:)),'k','LineWidth',2); hold on; % original
        xlim([0 30]);
        if subpl <= (nMIR); xticks([]);
        else
            xticks(0:5:30);
            xlabel('f[Hz]');
        end
        if subpl~=M+1 && subpl~=2*M+1
            yticks([])
        end
        ylim([0 2])
        set(gca,'FontSize',DimFont)
        subpl=subpl+1;
    end
end

% plot OIR with surrogates distribution
hold on
subplot(3,3,7)
area(f,OIR_perc_T{1,1}(3,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]); hold on;
area(f,OIR_perc_T{1,1}(1,:),'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
plot(f,OIR_perc_T{1,1}(1,:),'Color',[0.2 0.2 0.2],'LineWidth',1.2); hold on; % 2.5° prctile
plot(f,OIR_perc_T{1,1}(3,:),'Color',[0.2 0.2 0.2],'LineWidth',1.2); hold on; % 97.5° prctile
plot(f,OIR_perc_T{1,1}(2,:),'Color',[0.5 0.5 0.5],'LineWidth',1.2); hold on; % 50° prctile
plot(f,OIR_f,'r','LineWidth',2); hold on; % original
xlim([0 30]); xticks([]);
xticks(0:5:30);xlabel('f[Hz]');xlabel('f[Hz]');
set(gca,'FontSize',DimFont)
hold on
% barplot ER-MIR-OIR
% ER
ER_1=[ER_whole(1,1),ER_band(1,1,1),ER_band(1,1,2)];
ER_2=[ER_whole(2,2),ER_band(2,2,1),ER_band(2,2,2)];
ER_3=[ER_whole(3,3),ER_band(3,3,1),ER_band(3,3,2)];
ER=[ER_1;ER_2;ER_3];

MIR12_I=[ER_whole(1,2),ER_band(1,2,1),ER_band(1,2,2)];
MIR13_I=[ER_whole(1,3),ER_band(1,3,1),ER_band(1,3,2)];
MIR23_I=[ER_whole(2,3),ER_band(2,3,1),ER_band(2,3,2)];
MIR=[MIR12_I;MIR13_I;MIR23_I];
OIR_I=[OIR_whole,OIR_alpha,OIR_beta];

TOT=[ER;MIR;OIR_I];
subplot(3,3,[8 9])
AX=bar([TOT;NaN,NaN,NaN],'FaceColor','flat');
AX(1).CData(1:3,:)=repmat([1 0 0],3,1);
AX(2).CData(1:3,:)=repmat([1 0 0],3,1);
AX(3).CData(1:3,:)=repmat([160/255 160/255 160/255],3,1);
AX(1).CData(4:6,:)=repmat([0 0 0],3,1);
AX(2).CData(4:6,:)=repmat([0 0 0],3,1);
AX(3).CData(4:6,:)=repmat([0 0 0],3,1);
AX(1).CData(7,:)=[0 255/255 128/255];
AX(2).CData(7,:)=[0 255/255 128/255];
AX(3).CData(7,:)=[0 255/255 128/255];
set(gca,'XTickLabel',{'H_{X_1}','H_{X_2}','H_{X_3}','I_{X_1;X_2}','I_{X_1;X_3}','I_{X_2;X_3}','\Omega_{X^3}',''})
ylim([-1 1])
