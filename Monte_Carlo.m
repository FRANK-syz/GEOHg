clear;
clc;
%% 常量
tic;
Times=10000;

[Area,R]=readgeoraster("Monte Carlo/0.25x0.25base.tif"); 
Lc=readgeoraster("Monte Carlo/2015Lc025major.tif");

load("Monte Carlo/GFED5_month.mat");

AGBwoody=readgeoraster("Monte Carlo/live woody carbon 025.tif")/0.49;
AGBwoody=AGBwoody(:,:,11:20);
AGBwoody(isnan(AGBwoody))=0;
AGBwoody_uncertainvalue=readgeoraster("Monte Carlo/AGBwoody_uncertainty.tif");
AGBnonwoody=readgeoraster("Monte Carlo/2010AGB025.tif");
AGBnonwoody(AGBnonwoody<0)=0;
AGBnonwoody_uncertaintrate=readgeoraster("Monte Carlo/2010AGBuncertainty025.tif");
AGBnonwoody_uncertaintrate(AGBnonwoody_uncertaintrate<0)=0;
SLA=readgeoraster("Monte Carlo/Global_SLA.tif");
load("Monte Carlo/LAI.mat");
load("Monte Carlo/NDVI01-22.mat");
load("Monte Carlo/Tc.mat")

Foliar_Mc=ncread("Monte Carlo/GlobalFoliarHg_MachineLearning.nc","FoliarHg_concentrations");
Foliar_Mc(isnan(Foliar_Mc))=0;
Foliar_P5=ncread("Monte Carlo/GlobalFoliarHg_MachineLearning.nc",'FoliarHg_concentrations_P5');
Foliar_P5(isnan(Foliar_P5))=0;
Foliar_P95=ncread("Monte Carlo/GlobalFoliarHg_MachineLearning.nc",'FoliarHg_concentrations_P95');
Foliar_P95(isnan(Foliar_P95))=0;
% 计算对数正态分布参数mu和sigma
Foliar_sigma = (log(Foliar_P95) - log(Foliar_P5)) / (norminv(0.95) - norminv(0.05));
Foliar_sigma(isnan(Foliar_sigma))=0;
Foliar_sigma(Foliar_sigma==Inf)=0;
Foliar_mu = log(Foliar_P5) - norminv(0.05) * Foliar_sigma;
Foliar_mu(Foliar_mu<0)=0;

clear Foliar_P5 Foliar_P95

litterHgflux=readgeoraster("Monte Carlo/LitterHgflux.tif");
littersigma=readgeoraster("Monte Carlo/Littersigma.tif");

PeatCover=readgeoraster("Monte Carlo/Peatland025.tif");
% PeatHg=readgeoraster("Monte Carlo/PeatHg.tif");
load("Monte Carlo/PeatHglogP50.mat")
load("Monte Carlo/PeatHglogsigma.mat")
load("Monte Carlo/SW.mat")

t1=toc;
disp(['Initial Time:',num2str(t1,'%.1f'),'s'])
%% 循环
diary('Monte_log.txt'); % 开始记录
tic;

EmissionP50=zeros(720,1440,10);
EmissionP5=zeros(720,1440,10);
EmissionP95=zeros(720,1440,10);
EmissionMean=zeros(720,1440,10);
EmissionSigma=zeros(720,1440,10);

for lat=1:720
    for lon=1:1440
% for lat=431
%     t2=toc;
%     for lon=1431
        
        %判断是否有火
        pixelBA=squeeze(GFED5_month(lat,lon,:,:));
        if sum(pixelBA(:))==0
            continue
        else
            pixelemission_Times=zeros(Times,10);
            pixelLc=Lc(lat,lon);
            pixelAGBwoody=squeeze(AGBwoody(lat,lon,:));
            pixelAGBwoody_uncertainty=AGBwoody_uncertainvalue(lat,lon);
            pixelAGBnonwoody=AGBnonwoody(lat,lon);
            pixelAGBnonwoody_uncertainrate=AGBnonwoody_uncertaintrate(lat,lon);
            pixelLAI=squeeze(LAI(lat,lon,:));
            pixelSLA=squeeze(SLA(lat,lon,:));
            alpha=3.83;

            pixelNDVI01_22=squeeze(NDVI(lat,lon,:,:));
            pixelTc=squeeze(Tc(lat,lon,:,:));

            pixelFoliar_mu=Foliar_mu(lat,lon);
            pixelFoliar_sigma=Foliar_sigma(lat,lon);
            
            pixelLitterflux=litterHgflux(lat,lon);
            pixelLittersigma=littersigma(lat,lon);
            
            pixelPeatCover=PeatCover(lat,lon);
            % pixelPeatHg=PeatHg(lat,lon);
            pixelPeatHglogP50=logP50(lat,lon);
            pixelPeatHglogsigma=logsigma(lat,lon);
            pixelSW=squeeze(SW(lat,lon,:,:));

            % if pixelLc~=2
            %     pixelLc=0;
            % end

            switch pixelLc
                case 1
                    for times=1:Times
                        %BA(normal distribution 2sigma=50%)
                        BA_uncertainty=normrnd(pixelBA,pixelBA*0.25);

                        %AGB
                        NDVI01_22_uncertainty=unifrnd(0.9*pixelNDVI01_22,1.1*pixelNDVI01_22);
                        AGBwoody_uncertainty=normrnd(pixelAGBwoody,pixelAGBwoody_uncertainty);
    
                        Leaf_AGB_uncertainty=unifrnd(0.9*pixelLAI,1.1*pixelLAI)./unifrnd(0.9*pixelSLA,1.1*pixelSLA)*10*alpha;
                        Leaf_AGB_uncertainty(isnan(Leaf_AGB_uncertainty) | Leaf_AGB_uncertainty==Inf)=0;
                        Stem_AGB_uncertainty=AGBwoody_uncertainty-Leaf_AGB_uncertainty;
                        Branch_AGB_uncertainty=Stem_AGB_uncertainty*0.13;
                        Bolewood_AGB_uncertainty=(Stem_AGB_uncertainty-Branch_AGB_uncertainty).*0.85;
                        Bark_AGB_uncertainty=(Stem_AGB_uncertainty-Branch_AGB_uncertainty).*0.15;
                        Canopy_AGB_uncertainty=Leaf_AGB_uncertainty+Branch_AGB_uncertainty;

                        Leaf_AGB_uncertainty=transpose(Leaf_AGB_uncertainty);
                        Branch_AGB_uncertainty=transpose(Branch_AGB_uncertainty);
                        Bolewood_AGB_uncertainty=transpose(Bolewood_AGB_uncertainty);
                        Bark_AGB_uncertainty=transpose(Bark_AGB_uncertainty);
                        Canopy_AGB_uncertainty=transpose(Canopy_AGB_uncertainty);

                        %CF
                        Tc_uncertainty=unifrnd(0.9*pixelTc,1.1*pixelTc);
                        CF_uncertainty=calculate_CF(NDVI01_22_uncertainty,Tc_uncertainty,10,19);
                        CFleaf_uncertainty=Canopy_AGB_uncertainty.*CF_uncertainty./(Leaf_AGB_uncertainty+3/8*Branch_AGB_uncertainty);
                        CFleaf_uncertainty(CFleaf_uncertainty>1)=1;CFleaf_uncertainty(isnan(CFleaf_uncertainty))=0;
                        CFbranch_uncertainty=(Canopy_AGB_uncertainty.*CF_uncertainty-Leaf_AGB_uncertainty.*CFleaf_uncertainty)./Branch_AGB_uncertainty;
                        CFbranch_uncertainty(isnan(CFbranch_uncertainty))=0;CFbranch_uncertainty(CFbranch_uncertainty<0)=0;CFbranch_uncertainty(CFbranch_uncertainty>1)=0;
                        CFbolewood_uncertainty=CFbranch_uncertainty;CFbark_uncertainty=CFbranch_uncertainty;

                        %Mc
                        foliar_Mc_uncertainty=lognrnd(pixelFoliar_mu,pixelFoliar_sigma);
                        Leaf_Mc_uncertainty=lognrnd(logn_mu(27,14),logn_sigma(27,14));
                        Leaf_Mc_uncertainty(foliar_Mc_uncertainty~=0)=foliar_Mc_uncertainty(foliar_Mc_uncertainty~=0);
                        Branch_Mc_uncertainty=lognrnd(logn_mu(16,7),logn_sigma(16,7));
                        Bolewood_Mc_uncertainty=lognrnd(logn_mu(2,1),logn_sigma(2,1));
                        Bark_Mc_uncertainty=lognrnd(logn_mu(17,7),logn_sigma(17,7));
                        pixelLeaf_Ash=9.1;pixelBranch_Ash=4.6;pixelBolewood_Ash=1.8;pixelBark_Ash=4.1;%需要不确定性吗？？？？？？？？？

                        pixelLeafemission=BA_uncertainty.*Leaf_AGB_uncertainty.*CFleaf_uncertainty.*Leaf_Mc_uncertainty.*(1-pixelLeaf_Ash/100)./10^4;%kg
                        pixelBranchemission=BA_uncertainty.*Branch_AGB_uncertainty.*CFbranch_uncertainty.*Branch_Mc_uncertainty.*(1-pixelBranch_Ash/100)./10^4;%kg
                        pixelBolewoodemission=BA_uncertainty.*Bolewood_AGB_uncertainty.*CFbolewood_uncertainty.*Bolewood_Mc_uncertainty.*(1-pixelBolewood_Ash/100)./10^4;%kg
                        pixelBarkemission=BA_uncertainty.*Bark_AGB_uncertainty.*CFbark_uncertainty.*Bark_Mc_uncertainty.*(1-pixelBark_Ash/100)./10^4;%kg
                        
                        [pixelLitteremission,pixelPeatemission]=calculate_litter_peat(lat,BA_uncertainty,pixelLitterflux,pixelLittersigma,pixelPeatCover,pixelSW,pixelPeatHglogP50,pixelPeatHglogsigma);
                        
                        pixelTotalemission=pixelLeafemission+pixelBranchemission+pixelBolewoodemission+pixelBarkemission+pixelLitteremission+pixelPeatemission;
                        pixelannualemission=squeeze(sum(pixelTotalemission,1));
                        pixelemission_Times(times,:)=pixelannualemission;
                    end
                    %prctile for 10000 times
                    EmissionP95(lat,lon,:)=prctile(pixelemission_Times,95,1);
                    EmissionP5(lat,lon,:)=prctile(pixelemission_Times,5,1);
                    EmissionP50(lat,lon,:)=prctile(pixelemission_Times,50,1);
                    % EmissionMean(lat,lon,:)=mean(pixelemission_Times,1);
                    % EmissionSigma(lat,lon,:)=std(pixelemission_Times,0,1);

                case 2
                    for times=1:Times
                        %BA(normal distribution 2sigma=50%)
                        BA_uncertainty=normrnd(pixelBA,pixelBA*0.25);

                        %AGB
                        NDVI01_22_uncertainty=unifrnd(0.9*pixelNDVI01_22,1.1*pixelNDVI01_22);
                        AGBwoody_uncertainty=normrnd(pixelAGBwoody,pixelAGBwoody_uncertainty);
    
                        Leaf_AGB_uncertainty=unifrnd(0.9*pixelLAI,1.1*pixelLAI)./unifrnd(0.9*pixelSLA,1.1*pixelSLA)*10*alpha;
                        Leaf_AGB_uncertainty(isnan(Leaf_AGB_uncertainty) | Leaf_AGB_uncertainty==Inf)=0;
                        Stem_AGB_uncertainty=AGBwoody_uncertainty-Leaf_AGB_uncertainty;
                        Branch_AGB_uncertainty=Stem_AGB_uncertainty*0.28;
                        Bolewood_AGB_uncertainty=(Stem_AGB_uncertainty-Branch_AGB_uncertainty).*0.85;
                        Bark_AGB_uncertainty=(Stem_AGB_uncertainty-Branch_AGB_uncertainty).*0.15;
                        Canopy_AGB_uncertainty=Leaf_AGB_uncertainty+Branch_AGB_uncertainty;

                        Leaf_AGB_uncertainty=transpose(Leaf_AGB_uncertainty);
                        Branch_AGB_uncertainty=transpose(Branch_AGB_uncertainty);
                        Bolewood_AGB_uncertainty=transpose(Bolewood_AGB_uncertainty);
                        Bark_AGB_uncertainty=transpose(Bark_AGB_uncertainty);
                        Canopy_AGB_uncertainty=transpose(Canopy_AGB_uncertainty);

                        %CF
                        Tc_uncertainty=unifrnd(0.9*pixelTc,1.1*pixelTc);
                        CF_uncertainty=calculate_CF(NDVI01_22_uncertainty,Tc_uncertainty,10,19);
                        CFleaf_uncertainty=Canopy_AGB_uncertainty.*CF_uncertainty./(Leaf_AGB_uncertainty+3/9*Branch_AGB_uncertainty);
                        CFleaf_uncertainty(CFleaf_uncertainty>1)=1;CFleaf_uncertainty(isnan(CFleaf_uncertainty))=0;
                        CFbranch_uncertainty=(Canopy_AGB_uncertainty.*CF_uncertainty-Leaf_AGB_uncertainty.*CFleaf_uncertainty)./Branch_AGB_uncertainty;
                        CFbranch_uncertainty(isnan(CFbranch_uncertainty))=0;CFbranch_uncertainty(CFbranch_uncertainty<0)=0;CFbranch_uncertainty(CFbranch_uncertainty>1)=0;
                        CFbolewood_uncertainty=CFbranch_uncertainty;CFbark_uncertainty=CFbranch_uncertainty;

                        %Mc
                        foliar_Mc_uncertainty=lognrnd(pixelFoliar_mu,pixelFoliar_sigma);
                        Leaf_Mc_uncertainty=lognrnd(logn_mu(53,24),logn_sigma(53,24));
                        Leaf_Mc_uncertainty(foliar_Mc_uncertainty~=0)=foliar_Mc_uncertainty(foliar_Mc_uncertainty~=0);
                        Branch_Mc_uncertainty=lognrnd(logn_mu(12,11),logn_sigma(12,11));
                        Bolewood_Mc_uncertainty=lognrnd(logn_mu(2,1),logn_sigma(2,1));
                        Bark_Mc_uncertainty=lognrnd(logn_mu(4,3),logn_sigma(4,3));
                        pixelLeaf_Ash=9.1;pixelBranch_Ash=4.6;pixelBolewood_Ash=1.8;pixelBark_Ash=4.1;%需要不确定性吗？？？？？？？？？

                        pixelLeafemission=BA_uncertainty.*Leaf_AGB_uncertainty.*CFleaf_uncertainty.*Leaf_Mc_uncertainty.*(1-pixelLeaf_Ash/100)./10^4;%kg
                        pixelBranchemission=BA_uncertainty.*Branch_AGB_uncertainty.*CFbranch_uncertainty.*Branch_Mc_uncertainty.*(1-pixelBranch_Ash/100)./10^4;%kg
                        pixelBolewoodemission=BA_uncertainty.*Bolewood_AGB_uncertainty.*CFbolewood_uncertainty.*Bolewood_Mc_uncertainty.*(1-pixelBolewood_Ash/100)./10^4;%kg
                        pixelBarkemission=BA_uncertainty.*Bark_AGB_uncertainty.*CFbark_uncertainty.*Bark_Mc_uncertainty.*(1-pixelBark_Ash/100)./10^4;%kg

                        [pixelLitteremission,pixelPeatemission]=calculate_litter_peat(lat,BA_uncertainty,pixelLitterflux,pixelLittersigma,pixelPeatCover,pixelSW,pixelPeatHglogP50,pixelPeatHglogsigma);
                        
                        pixelTotalemission=pixelLeafemission+pixelBranchemission+pixelBolewoodemission+pixelBarkemission+pixelLitteremission+pixelPeatemission;
                        pixelannualemission=squeeze(sum(pixelTotalemission,1));
                        pixelemission_Times(times,:)=pixelannualemission;
                    end
                    %prctile for 10000 times
                    EmissionP95(lat,lon,:)=prctile(pixelemission_Times,95,1);
                    EmissionP5(lat,lon,:)=prctile(pixelemission_Times,5,1);
                    EmissionP50(lat,lon,:)=prctile(pixelemission_Times,50,1);
                    % EmissionMean(lat,lon,:)=mean(pixelemission_Times,1);
                    % EmissionSigma(lat,lon,:)=std(pixelemission_Times,0,1);

                case 3
                    for times=1:Times
                        %BA(normal distribution 2sigma=50%)
                        BA_uncertainty=normrnd(pixelBA,pixelBA*0.25);

                        %AGB
                        NDVI01_22_uncertainty=unifrnd(0.9*pixelNDVI01_22,1.1*pixelNDVI01_22);
                        AGBwoody_uncertainty=normrnd(pixelAGBwoody,pixelAGBwoody_uncertainty);
    
                        Leaf_AGB_uncertainty=unifrnd(0.9*pixelLAI,1.1*pixelLAI)./unifrnd(0.9*pixelSLA,1.1*pixelSLA)*10*alpha;
                        Leaf_AGB_uncertainty(isnan(Leaf_AGB_uncertainty) | Leaf_AGB_uncertainty==Inf)=0;
                        Stem_AGB_uncertainty=AGBwoody_uncertainty-Leaf_AGB_uncertainty;
                        Branch_AGB_uncertainty=Stem_AGB_uncertainty*0.08;
                        Bolewood_AGB_uncertainty=(Stem_AGB_uncertainty-Branch_AGB_uncertainty).*0.85;
                        Bark_AGB_uncertainty=(Stem_AGB_uncertainty-Branch_AGB_uncertainty).*0.15;
                        Canopy_AGB_uncertainty=Leaf_AGB_uncertainty+Branch_AGB_uncertainty;
                        
                        Leaf_AGB_uncertainty=transpose(Leaf_AGB_uncertainty);
                        Branch_AGB_uncertainty=transpose(Branch_AGB_uncertainty);
                        Bolewood_AGB_uncertainty=transpose(Bolewood_AGB_uncertainty);
                        Bark_AGB_uncertainty=transpose(Bark_AGB_uncertainty);
                        Canopy_AGB_uncertainty=transpose(Canopy_AGB_uncertainty);

                        %CF
                        Tc_uncertainty=unifrnd(0.9*pixelTc,1.1*pixelTc);
                        CF_uncertainty=calculate_CF(NDVI01_22_uncertainty,Tc_uncertainty,10,19);
                        CFleaf_uncertainty=Canopy_AGB_uncertainty.*CF_uncertainty./(Leaf_AGB_uncertainty+3/8*Branch_AGB_uncertainty);
                        CFleaf_uncertainty(CFleaf_uncertainty>1)=1;CFleaf_uncertainty(isnan(CFleaf_uncertainty))=0;
                        CFbranch_uncertainty=(Canopy_AGB_uncertainty.*CF_uncertainty-Leaf_AGB_uncertainty.*CFleaf_uncertainty)./Branch_AGB_uncertainty;
                        CFbranch_uncertainty(isnan(CFbranch_uncertainty))=0;CFbranch_uncertainty(CFbranch_uncertainty<0)=0;CFbranch_uncertainty(CFbranch_uncertainty>1)=0;
                        CFbolewood_uncertainty=CFbranch_uncertainty;CFbark_uncertainty=CFbranch_uncertainty;

                        %Mc
                        foliar_Mc_uncertainty=lognrnd(pixelFoliar_mu,pixelFoliar_sigma);
                        Leaf_Mc_uncertainty=lognrnd(logn_mu(49,16),logn_sigma(49,16));
                        Leaf_Mc_uncertainty(foliar_Mc_uncertainty~=0)=foliar_Mc_uncertainty(foliar_Mc_uncertainty~=0);
                        Branch_Mc_uncertainty=lognrnd(logn_mu(19,2),logn_sigma(19,2));
                        Bolewood_Mc_uncertainty=lognrnd(logn_mu(3,0),logn_sigma(3,0));
                        Bark_Mc_uncertainty=lognrnd(logn_mu(19,2),logn_sigma(19,2));
                        pixelLeaf_Ash=9.1;pixelBranch_Ash=4.6;pixelBolewood_Ash=1.8;pixelBark_Ash=4.1;%需要不确定性吗？？？？？？？？？

                        pixelLeafemission=BA_uncertainty.*Leaf_AGB_uncertainty.*CFleaf_uncertainty.*Leaf_Mc_uncertainty.*(1-pixelLeaf_Ash/100)./10^4;%kg
                        pixelBranchemission=BA_uncertainty.*Branch_AGB_uncertainty.*CFbranch_uncertainty.*Branch_Mc_uncertainty.*(1-pixelBranch_Ash/100)./10^4;%kg
                        pixelBolewoodemission=BA_uncertainty.*Bolewood_AGB_uncertainty.*CFbolewood_uncertainty.*Bolewood_Mc_uncertainty.*(1-pixelBolewood_Ash/100)./10^4;%kg
                        pixelBarkemission=BA_uncertainty.*Bark_AGB_uncertainty.*CFbark_uncertainty.*Bark_Mc_uncertainty.*(1-pixelBark_Ash/100)./10^4;%kg

                        [pixelLitteremission,pixelPeatemission]=calculate_litter_peat(lat,BA_uncertainty,pixelLitterflux,pixelLittersigma,pixelPeatCover,pixelSW,pixelPeatHglogP50,pixelPeatHglogsigma);
                        
                        pixelTotalemission=pixelLeafemission+pixelBranchemission+pixelBolewoodemission+pixelBarkemission+pixelLitteremission+pixelPeatemission;
                        pixelannualemission=squeeze(sum(pixelTotalemission,1));
                        pixelemission_Times(times,:)=pixelannualemission;
                    end
                    %prctile for 10000 times
                    EmissionP95(lat,lon,:)=prctile(pixelemission_Times,95,1);
                    EmissionP5(lat,lon,:)=prctile(pixelemission_Times,5,1);
                    EmissionP50(lat,lon,:)=prctile(pixelemission_Times,50,1);
                    % EmissionMean(lat,lon,:)=mean(pixelemission_Times,1);
                    % EmissionSigma(lat,lon,:)=std(pixelemission_Times,0,1);

                case 4
                    for times=1:Times
                        %BA(normal distribution 2sigma=50%)
                        BA_uncertainty=normrnd(pixelBA,pixelBA*0.25);

                        %AGB
                        NDVI01_22_uncertainty=unifrnd(0.9*pixelNDVI01_22,1.1*pixelNDVI01_22);
                        AGBwoody_uncertainty=normrnd(pixelAGBwoody,pixelAGBwoody_uncertainty);
    
                        Leaf_AGB_uncertainty=unifrnd(0.9*pixelLAI,1.1*pixelLAI)./unifrnd(0.9*pixelSLA,1.1*pixelSLA)*10*alpha;
                        Leaf_AGB_uncertainty(isnan(Leaf_AGB_uncertainty) | Leaf_AGB_uncertainty==Inf)=0;
                        Stem_AGB_uncertainty=AGBwoody_uncertainty-Leaf_AGB_uncertainty;
                        Branch_AGB_uncertainty=Stem_AGB_uncertainty*0.23;
                        Bolewood_AGB_uncertainty=(Stem_AGB_uncertainty-Branch_AGB_uncertainty).*0.85;
                        Bark_AGB_uncertainty=(Stem_AGB_uncertainty-Branch_AGB_uncertainty).*0.15;
                        Canopy_AGB_uncertainty=Leaf_AGB_uncertainty+Branch_AGB_uncertainty;

                        Leaf_AGB_uncertainty=transpose(Leaf_AGB_uncertainty);
                        Branch_AGB_uncertainty=transpose(Branch_AGB_uncertainty);
                        Bolewood_AGB_uncertainty=transpose(Bolewood_AGB_uncertainty);
                        Bark_AGB_uncertainty=transpose(Bark_AGB_uncertainty);
                        Canopy_AGB_uncertainty=transpose(Canopy_AGB_uncertainty);
                        
                        %CF
                        Tc_uncertainty=unifrnd(0.9*pixelTc,1.1*pixelTc);
                        CF_uncertainty=calculate_CF(NDVI01_22_uncertainty,Tc_uncertainty,10,19);
                        CFleaf_uncertainty=Canopy_AGB_uncertainty.*CF_uncertainty./(Leaf_AGB_uncertainty+5/9*Branch_AGB_uncertainty);
                        CFleaf_uncertainty(CFleaf_uncertainty>1)=1;CFleaf_uncertainty(isnan(CFleaf_uncertainty))=0;
                        CFbranch_uncertainty=(Canopy_AGB_uncertainty.*CF_uncertainty-Leaf_AGB_uncertainty.*CFleaf_uncertainty)./Branch_AGB_uncertainty;
                        CFbranch_uncertainty(isnan(CFbranch_uncertainty))=0;CFbranch_uncertainty(CFbranch_uncertainty<0)=0;CFbranch_uncertainty(CFbranch_uncertainty>1)=0;
                        CFbolewood_uncertainty=CFbranch_uncertainty;CFbark_uncertainty=CFbranch_uncertainty;
                        
                        %Mc
                        foliar_Mc_uncertainty=lognrnd(pixelFoliar_mu,pixelFoliar_sigma);
                        Leaf_Mc_uncertainty=lognrnd(logn_mu(41,15),logn_sigma(41,15));
                        Leaf_Mc_uncertainty(foliar_Mc_uncertainty~=0)=foliar_Mc_uncertainty(foliar_Mc_uncertainty~=0);
                        Branch_Mc_uncertainty=lognrnd(logn_mu(12,5),logn_sigma(12,5));
                        Bolewood_Mc_uncertainty=lognrnd(logn_mu(2,1),logn_sigma(2,1));
                        Bark_Mc_uncertainty=lognrnd(logn_mu(9,6),logn_sigma(9,6));
                        pixelLeaf_Ash=9.1;pixelBranch_Ash=4.6;pixelBolewood_Ash=1.8;pixelBark_Ash=4.1;%需要不确定性吗？？？？？？？？？

                        pixelLeafemission=BA_uncertainty.*Leaf_AGB_uncertainty.*CFleaf_uncertainty.*Leaf_Mc_uncertainty.*(1-pixelLeaf_Ash/100)./10^4;%kg
                        pixelBranchemission=BA_uncertainty.*Branch_AGB_uncertainty.*CFbranch_uncertainty.*Branch_Mc_uncertainty.*(1-pixelBranch_Ash/100)./10^4;%kg
                        pixelBolewoodemission=BA_uncertainty.*Bolewood_AGB_uncertainty.*CFbolewood_uncertainty.*Bolewood_Mc_uncertainty.*(1-pixelBolewood_Ash/100)./10^4;%kg
                        pixelBarkemission=BA_uncertainty.*Bark_AGB_uncertainty.*CFbark_uncertainty.*Bark_Mc_uncertainty.*(1-pixelBark_Ash/100)./10^4;%kg

                        [pixelLitteremission,pixelPeatemission]=calculate_litter_peat(lat,BA_uncertainty,pixelLitterflux,pixelLittersigma,pixelPeatCover,pixelSW,pixelPeatHglogP50,pixelPeatHglogsigma);
                        
                        pixelTotalemission=pixelLeafemission+pixelBranchemission+pixelBolewoodemission+pixelBarkemission+pixelLitteremission+pixelPeatemission;
                        pixelannualemission=squeeze(sum(pixelTotalemission,1));
                        pixelemission_Times(times,:)=pixelannualemission;
                    end
                    %prctile for 10000 times
                    EmissionP95(lat,lon,:)=prctile(pixelemission_Times,95,1);
                    EmissionP5(lat,lon,:)=prctile(pixelemission_Times,5,1);
                    EmissionP50(lat,lon,:)=prctile(pixelemission_Times,50,1);
                    % EmissionMean(lat,lon,:)=mean(pixelemission_Times,1);
                    % EmissionSigma(lat,lon,:)=std(pixelemission_Times,0,1);

                case 5
                    for times=1:Times
                        %BA(normal distribution 2sigma=50%)
                        BA_uncertainty=normrnd(pixelBA,pixelBA*0.25);

                        %AGB
                        NDVI01_22_uncertainty=unifrnd(0.9*pixelNDVI01_22,1.1*pixelNDVI01_22);
                        AGBwoody_uncertainty=normrnd(pixelAGBwoody,pixelAGBwoody_uncertainty);
    
                        Leaf_AGB_uncertainty=unifrnd(0.9*pixelLAI,1.1*pixelLAI)./unifrnd(0.9*pixelSLA,1.1*pixelSLA)*10*alpha;
                        Leaf_AGB_uncertainty(isnan(Leaf_AGB_uncertainty) | Leaf_AGB_uncertainty==Inf)=0;
                        Stem_AGB_uncertainty=AGBwoody_uncertainty-Leaf_AGB_uncertainty;
                        Branch_AGB_uncertainty=Stem_AGB_uncertainty*0.18;
                        Bolewood_AGB_uncertainty=(Stem_AGB_uncertainty-Branch_AGB_uncertainty).*0.85;
                        Bark_AGB_uncertainty=(Stem_AGB_uncertainty-Branch_AGB_uncertainty).*0.15;
                        Canopy_AGB_uncertainty=Leaf_AGB_uncertainty+Branch_AGB_uncertainty;

                        Leaf_AGB_uncertainty=transpose(Leaf_AGB_uncertainty);
                        Branch_AGB_uncertainty=transpose(Branch_AGB_uncertainty);
                        Bolewood_AGB_uncertainty=transpose(Bolewood_AGB_uncertainty);
                        Bark_AGB_uncertainty=transpose(Bark_AGB_uncertainty);
                        Canopy_AGB_uncertainty=transpose(Canopy_AGB_uncertainty);
                        
                        %CF
                        Tc_uncertainty=unifrnd(0.9*pixelTc,1.1*pixelTc);
                        CF_uncertainty=calculate_CF(NDVI01_22_uncertainty,Tc_uncertainty,10,19);
                        CFleaf_uncertainty=Canopy_AGB_uncertainty.*CF_uncertainty./(Leaf_AGB_uncertainty+5/8*Branch_AGB_uncertainty);
                        CFleaf_uncertainty(CFleaf_uncertainty>1)=1;CFleaf_uncertainty(isnan(CFleaf_uncertainty))=0;
                        CFbranch_uncertainty=(Canopy_AGB_uncertainty.*CF_uncertainty-Leaf_AGB_uncertainty.*CFleaf_uncertainty)./Branch_AGB_uncertainty;
                        CFbranch_uncertainty(isnan(CFbranch_uncertainty))=0;CFbranch_uncertainty(CFbranch_uncertainty<0)=0;CFbranch_uncertainty(CFbranch_uncertainty>1)=0;
                        CFbolewood_uncertainty=CFbranch_uncertainty;CFbark_uncertainty=CFbranch_uncertainty;

                        %Mc
                        foliar_Mc_uncertainty=lognrnd(pixelFoliar_mu,pixelFoliar_sigma);
                        Leaf_Mc_uncertainty=lognrnd(logn_mu(25,15),logn_sigma(25,15));
                        Leaf_Mc_uncertainty(foliar_Mc_uncertainty~=0)=foliar_Mc_uncertainty(foliar_Mc_uncertainty~=0);
                        Branch_Mc_uncertainty=lognrnd(logn_mu(5,2),logn_sigma(5,2));
                        Bolewood_Mc_uncertainty=lognrnd(logn_mu(3,1),logn_sigma(3,1));
                        Bark_Mc_uncertainty=lognrnd(logn_mu(13,5),logn_sigma(13,5));
                        pixelLeaf_Ash=9.1;pixelBranch_Ash=4.6;pixelBolewood_Ash=1.8;pixelBark_Ash=4.1;%需要不确定性吗？？？？？？？？？

                        pixelLeafemission=BA_uncertainty.*Leaf_AGB_uncertainty.*CFleaf_uncertainty.*Leaf_Mc_uncertainty.*(1-pixelLeaf_Ash/100)./10^4;%kg
                        pixelBranchemission=BA_uncertainty.*Branch_AGB_uncertainty.*CFbranch_uncertainty.*Branch_Mc_uncertainty.*(1-pixelBranch_Ash/100)./10^4;%kg
                        pixelBolewoodemission=BA_uncertainty.*Bolewood_AGB_uncertainty.*CFbolewood_uncertainty.*Bolewood_Mc_uncertainty.*(1-pixelBolewood_Ash/100)./10^4;%kg
                        pixelBarkemission=BA_uncertainty.*Bark_AGB_uncertainty.*CFbark_uncertainty.*Bark_Mc_uncertainty.*(1-pixelBark_Ash/100)./10^4;%kg

                        [pixelLitteremission,pixelPeatemission]=calculate_litter_peat(lat,BA_uncertainty,pixelLitterflux,pixelLittersigma,pixelPeatCover,pixelSW,pixelPeatHglogP50,pixelPeatHglogsigma);
                        
                        pixelTotalemission=pixelLeafemission+pixelBranchemission+pixelBolewoodemission+pixelBarkemission+pixelLitteremission+pixelPeatemission;
                        pixelannualemission=squeeze(sum(pixelTotalemission,1));
                        pixelemission_Times(times,:)=pixelannualemission;
                    end
                    %prctile for 10000 times
                    EmissionP95(lat,lon,:)=prctile(pixelemission_Times,95,1);
                    EmissionP5(lat,lon,:)=prctile(pixelemission_Times,5,1);
                    EmissionP50(lat,lon,:)=prctile(pixelemission_Times,50,1);
                    % EmissionMean(lat,lon,:)=mean(pixelemission_Times,1);
                    % EmissionSigma(lat,lon,:)=std(pixelemission_Times,0,1);

                case {6,7,8}
                    for times=1:Times
                        %BA(normal distribution 2sigma=50%)
                        BA_uncertainty=normrnd(pixelBA,pixelBA*0.25);

                        %AGB
                        NDVI01_22_uncertainty=unifrnd(0.9*pixelNDVI01_22,1.1*pixelNDVI01_22);
                        AGBwoody_uncertainty=normrnd(pixelAGBwoody,pixelAGBwoody_uncertainty);
    
                        Leaf_AGB_uncertainty=unifrnd(0.9*pixelLAI,1.1*pixelLAI)./unifrnd(0.9*pixelSLA,1.1*pixelSLA)*10*alpha;
                        Leaf_AGB_uncertainty(isnan(Leaf_AGB_uncertainty) | Leaf_AGB_uncertainty==Inf)=0;
                        Branch_AGB_uncertainty=AGBwoody_uncertainty-Leaf_AGB_uncertainty;

                        AGBwoody_uncertainty=transpose(AGBwoody_uncertainty);
                        Leaf_AGB_uncertainty=transpose(Leaf_AGB_uncertainty);
                        Branch_AGB_uncertainty=transpose(Branch_AGB_uncertainty);

                        %CF
                        Tc_uncertainty=unifrnd(0.9*pixelTc,1.1*pixelTc);
                        CF_uncertainty=calculate_CF(NDVI01_22_uncertainty,Tc_uncertainty,10,19);
                        CFleaf_uncertainty=AGBwoody_uncertainty.*CF_uncertainty./(Leaf_AGB_uncertainty+6/9*Branch_AGB_uncertainty);
                        CFleaf_uncertainty(CFleaf_uncertainty>1)=1;CFleaf_uncertainty(isnan(CFleaf_uncertainty))=0;
                        CFbranch_uncertainty=(AGBwoody_uncertainty.*CF_uncertainty-Leaf_AGB_uncertainty.*CFleaf_uncertainty)./Branch_AGB_uncertainty;
                        CFbranch_uncertainty(isnan(CFbranch_uncertainty))=0;CFbranch_uncertainty(CFbranch_uncertainty<0)=0;CFbranch_uncertainty(CFbranch_uncertainty>1)=0;
                        
                        %Mc
                        foliar_Mc_uncertainty=lognrnd(pixelFoliar_mu,pixelFoliar_sigma);
                        Leaf_Mc_uncertainty=lognrnd(logn_mu(19,17),logn_sigma(19,17));
                        Leaf_Mc_uncertainty(foliar_Mc_uncertainty~=0)=foliar_Mc_uncertainty(foliar_Mc_uncertainty~=0);
                        Branch_Mc_uncertainty=lognrnd(logn_mu(6,7),logn_sigma(6,7));
                        pixelLeaf_Ash=8.75;pixelBranch_Ash=4.6;%需要不确定性吗？？？？？？？？？

                        pixelLeafemission=BA_uncertainty.*Leaf_AGB_uncertainty.*CFleaf_uncertainty.*Leaf_Mc_uncertainty.*(1-pixelLeaf_Ash/100)./10^4;%kg
                        pixelBranchemission=BA_uncertainty.*Branch_AGB_uncertainty.*CFbranch_uncertainty.*Branch_Mc_uncertainty.*(1-pixelBranch_Ash/100)./10^4;%kg
                        
                        [pixelLitteremission,pixelPeatemission]=calculate_litter_peat(lat,BA_uncertainty,pixelLitterflux,pixelLittersigma,pixelPeatCover,pixelSW,pixelPeatHglogP50,pixelPeatHglogsigma);
                        
                        pixelTotalemission=pixelLeafemission+pixelBranchemission+pixelLitteremission+pixelPeatemission;
                        pixelannualemission=squeeze(sum(pixelTotalemission,1));
                        pixelemission_Times(times,:)=pixelannualemission;
                    end
                    %prctile for 10000 times
                    EmissionP95(lat,lon,:)=prctile(pixelemission_Times,95,1);
                    EmissionP5(lat,lon,:)=prctile(pixelemission_Times,5,1);
                    EmissionP50(lat,lon,:)=prctile(pixelemission_Times,50,1);
                    % EmissionMean(lat,lon,:)=mean(pixelemission_Times,1);
                    % EmissionSigma(lat,lon,:)=std(pixelemission_Times,0,1);

                case {9,10,11,14,16}
                    for times=1:Times
                        %BA(normal distribution 2sigma=50%)
                        BA_uncertainty=normrnd(pixelBA,pixelBA*0.25);

                        %AGB
                        NDVI01_22_uncertainty=unifrnd(0.9*pixelNDVI01_22,1.1*pixelNDVI01_22);
                        AnnualNDVI_uncertainty=transpose(squeeze(max(NDVI01_22_uncertainty(:,10:19),[],1)));
                        AGBnonwoody_uncertainty=normrnd(pixelAGBnonwoody,pixelAGBnonwoody.*pixelAGBnonwoody_uncertainrate);
                        AGBnonwoody_uncertainty=AnnualNDVI_uncertainty./AnnualNDVI_uncertainty(1).*AGBnonwoody_uncertainty;
                        AGBnonwoody_uncertainty(isnan(AGBnonwoody_uncertainty))=0;

                        Leaf_AGB_uncertainty=unifrnd(0.9*pixelLAI,1.1*pixelLAI)./unifrnd(0.9*pixelSLA,1.1*pixelSLA)*10*alpha;
                        Leaf_AGB_uncertainty(isnan(Leaf_AGB_uncertainty) | Leaf_AGB_uncertainty==Inf)=0;
                        Leaf_AGB_uncertainty=transpose(Leaf_AGB_uncertainty);

                        %CF
                        Tc_uncertainty=unifrnd(0.9*pixelTc,1.1*pixelTc);
                        CFleaf_uncertainty=calculate_CF(NDVI01_22_uncertainty,Tc_uncertainty,10,19);

                        %Mc
                        foliar_Mc_uncertainty=lognrnd(pixelFoliar_mu,pixelFoliar_sigma);
                        Leaf_Mc_uncertainty=lognrnd(logn_mu(20,10),logn_sigma(20,10));
                        Leaf_Mc_uncertainty(foliar_Mc_uncertainty~=0)=foliar_Mc_uncertainty(foliar_Mc_uncertainty~=0);
                        pixelLeaf_Ash=8.4;%需要不确定性吗？？？？？？？？？

                        pixelLeafemission=BA_uncertainty.*Leaf_AGB_uncertainty.*CFleaf_uncertainty.*Leaf_Mc_uncertainty.*(1-pixelLeaf_Ash/100)./10^4;%kg
                        
                        [pixelLitteremission,pixelPeatemission]=calculate_litter_peat(lat,BA_uncertainty,pixelLitterflux,pixelLittersigma,pixelPeatCover,pixelSW,pixelPeatHglogP50,pixelPeatHglogsigma);
                        
                        pixelTotalemission=pixelLeafemission+pixelLitteremission+pixelPeatemission;
                        pixelannualemission=squeeze(sum(pixelTotalemission,1));
                        pixelemission_Times(times,:)=pixelannualemission;
                    end
                    %prctile for 10000 times
                    EmissionP95(lat,lon,:)=prctile(pixelemission_Times,95,1);
                    EmissionP5(lat,lon,:)=prctile(pixelemission_Times,5,1);
                    EmissionP50(lat,lon,:)=prctile(pixelemission_Times,50,1);
                    % EmissionMean(lat,lon,:)=mean(pixelemission_Times,1);
                    % EmissionSigma(lat,lon,:)=std(pixelemission_Times,0,1);

                case 12
                    for times=1:Times
                        %BA(normal distribution 2sigma=50%)
                        BA_uncertainty=normrnd(pixelBA,pixelBA*0.25);
                        
                        %AGB/CF
                        NDVI01_22_uncertainty=unifrnd(0.9*pixelNDVI01_22,1.1*pixelNDVI01_22);
                        AnnualNDVI_uncertainty=squeeze(max(NDVI01_22_uncertainty(:,10:19),[],1));
                        AGBnonwoody_uncertainty=normrnd(pixelAGBnonwoody,pixelAGBnonwoody.*pixelAGBnonwoody_uncertainrate);
                        AGBnonwoody_uncertainty=(AnnualNDVI_uncertainty)/AnnualNDVI_uncertainty(1)*AGBnonwoody_uncertainty;
                        AGBnonwoody_uncertainty(isnan(AGBnonwoody_uncertainty))=0;
                        
                        Leaf_AGB_uncertainty=0.45*AGBnonwoody_uncertainty;Branch_AGB_uncertainty=0.55*AGBnonwoody_uncertainty;
                        CFleaf_uncertainty=0.9;CFbranch_uncertainty=0.9;
                        
                        %Mc
                        foliar_Mc_uncertainty=lognrnd(pixelFoliar_mu,pixelFoliar_sigma);
                        Leaf_Mc_uncertainty=lognrnd(logn_mu(28,16),logn_sigma(28,16));
                        Branch_Mc_uncertainty=lognrnd(logn_mu(21,13),logn_sigma(21,13));
                        Leaf_Mc_uncertainty(foliar_Mc_uncertainty~=0)=foliar_Mc_uncertainty(foliar_Mc_uncertainty~=0);
                        pixelLeaf_Ash=5.9;pixelBranch_Ash=5.9;%需要不确定性吗？？？？？？？？？

                        pixelLeafemission=BA_uncertainty.*Leaf_AGB_uncertainty.*CFleaf_uncertainty.*Leaf_Mc_uncertainty.*(1-pixelLeaf_Ash/100)./10^4;%kg
                        pixelBranchemission=BA_uncertainty.*Branch_AGB_uncertainty.*CFbranch_uncertainty.*Branch_Mc_uncertainty.*(1-pixelBranch_Ash/100)./10^4;%kg
                        
                        [pixelLitteremission,pixelPeatemission]=calculate_litter_peat(lat,BA_uncertainty,pixelLitterflux,pixelLittersigma,pixelPeatCover,pixelSW,pixelPeatHglogP50,pixelPeatHglogsigma);
                        
                        pixelTotalemission=pixelLeafemission+pixelBranchemission+pixelLitteremission+pixelPeatemission;
                        pixelannualemission=squeeze(sum(pixelTotalemission,1));
                        pixelemission_Times(times,:)=pixelannualemission;
                    end
                    %prctile for 10000 times
                    EmissionP95(lat,lon,:)=prctile(pixelemission_Times,95,1);
                    EmissionP5(lat,lon,:)=prctile(pixelemission_Times,5,1);
                    EmissionP50(lat,lon,:)=prctile(pixelemission_Times,50,1);
                    % EmissionMean(lat,lon,:)=mean(pixelemission_Times,1);
                    EmissionSigma(lat,lon,:)=std(pixelemission_Times,0,1);
            end

        end
    end

    % 进度显示
    interval=20;
    current_iter=lat;
    max_iter=720;
    tcurrent=round(toc);
    if mod(current_iter, round(max_iter/interval)) == 0 || current_iter == max_iter
        progress = current_iter/max_iter * 100;
        fprintf('progress %.1f%% (%d/%d)\n', progress, current_iter, max_iter);
        disp(['Spend Time:',num2str(tcurrent),'s'])
        
        % [usr, sys] = memory;
        % disp(['当前 MATLAB 占用内存: ', num2str(usr.MemUsedMATLAB/1e6), ' MB']);
        % disp(['系统可用内存: ', num2str(sys.PhysicalMemory.Available/1e6), ' MB']);
    end
end

tEnd=round(toc);
disp(['TotalTime:',num2str(tEnd),'s'])

%% 结算
Multiyear_P50=mean(EmissionP50,3);
Multiyear_P5=mean(EmissionP5,3);
Multiyear_P95=mean(EmissionP95,3);
% Multiyear_Mean=mean(EmissionMean,3);
Multiyear_Sigma=mean(EmissionSigma,3);
% Multiyear_CV=Multiyear_Sigma./Multiyear_Mean;
Multiyear_uncertainty=(Multiyear_P95-Multiyear_P5)./Multiyear_P50;

AnnualP50=squeeze(sum(EmissionP50,[1 2]))./10^3;%t
AnnualP5=squeeze(sum(EmissionP5,[1 2]))./10^3;%t
AnnualP95=squeeze(sum(EmissionP95,[1 2]))./10^3;%t
AnnualMean=squeeze(sum(EmissionMean,[1 2]))./10^3;%t
AnnualSigma=squeeze(sum(EmissionSigma,[1 2]))./10^3;%t
% figure(1)
% img1=Multiyear_P50;
% imagesc(img1);colorbar;
% clim([0 20]);
figure(2)
bar(1:10, AnnualP50, 'FaceColor', [0.7 0.7 0.9]);
hold on;
errorbar(1:10,AnnualP50,AnnualP50-AnnualSigma,AnnualP50+AnnualSigma,'LineStyle', 'none')
% figure(3)
% img3=Multiyear_uncertainty;
% imagesc(img3);colorbar
% clim([0 2]);

save("Monte Carlo/Multiyear_P5.mat",'Multiyear_P5');
save("Monte Carlo/Multiyear_P50.mat",'Multiyear_P50');
save("Monte Carlo/Multiyear_P95.mat",'Multiyear_P95');
save("Monte Carlo/Multiyear_uncertainty.mat",'Multiyear_uncertainty');
save("Monte Carlo/Multiyear_Sigma.mat",'Multiyear_Sigma');

disp('P5:')
disp(transpose(AnnualP5))
disp('P50:')
disp(transpose(AnnualP50))
disp('P95:')
disp(transpose(AnnualP95))
disp(mean(AnnualP5))
disp(mean(AnnualP50))
disp(mean(AnnualP95))
diary off
%% 

function mu=logn_mu(m,v)
    mu = log((m^2)/sqrt(v+m^2));
end

function sigma=logn_sigma(m,v)
    sigma = sqrt(log(v/(m^2)+1));
end

function CF=calculate_CF(NDVI,Tc,Start,End)
    CF=zeros(size(Tc));
    maxNDVI=max(NDVI(:));minNDVI=min(NDVI(:));
    if maxNDVI-minNDVI~=0
        VCI=(NDVI(:,Start:End)-minNDVI)/(maxNDVI-minNDVI);
        mcf=VCI;
        mcf(mcf>0 & mcf<=(100/6))=0.33;
        mcf(mcf>(100/6) & mcf<=(200/6))=0.5;
        mcf(mcf>(200/6) & mcf<=(300/6))=1;
        mcf(mcf>(300/6) & mcf<=(400/6))=2;
        mcf(mcf>(400/6) & mcf<=(500/6))=4;
        mcf(mcf>(500/6) & mcf<=(100))=5;
    
        CF(Tc>0 & Tc<=40)=-2.13.*VCI(Tc>0 & Tc<=40)./100+1.38;
        CF(Tc>0 & Tc<=40 & CF<0.44)=0.44;
        CF(Tc>0 & Tc<=40 & CF>0.98)=0.98;
        CF(Tc>40 & Tc<=60)=exp(-0.013.*Tc(Tc>40 & Tc<=60));
        CF(Tc>60 & Tc<100)=(1-1/exp(1)).^mcf(Tc>60 & Tc<100);  
    end
end

function [pixelLitteremission,pixelPeatemission]=calculate_litter_peat(lat,BA_uncertainty,pixelLitterflux,pixelLittersigma,pixelPeatCover,pixelSW,pixelPeatHglogP50,pixelPeatHglogsigma)
    %Litterfall
    LitterCF=0.8;
    Litter_Ash=9.1;
    litterflux_uncertainty=normrnd(pixelLitterflux,pixelLittersigma);
    pixelLitteremission=BA_uncertainty.*LitterCF.*litterflux_uncertainty*(1-Litter_Ash/100)./10^3;%kg;
    
    %Peat
    pixelPeatemission=0;
    PeatCF=0.8;
    Peat_Ash=4.9;
    BD=0.1;%g/cm3
    if pixelPeatCover>5
        PeatCover_uncertainty=unifrnd(0.9*pixelPeatCover,1.1*pixelPeatCover);
        SW_uncertainty=unifrnd(0.9*pixelSW,1.1*pixelSW);
        DOB_uncertainty=zeros(size(pixelSW));
        if lat>200 && lat<=540
            DOB_uncertainty=-51.*SW_uncertainty./0.45+57;
            DOB_uncertainty(DOB_uncertainty<0)=0;
        elseif lat<=200 || lat>540
            BorealSW_uncertainty=SW_uncertainty;
            BorealSW_uncertainty(BorealSW_uncertainty<=0.2)=NaN;
            DOB_uncertainty=236.*(42*10*BorealSW_uncertainty./(10*BorealSW_uncertainty+1)-28).^(-0.38)-100;
            DOB_uncertainty(DOB_uncertainty<0)=0;DOB_uncertainty(DOB_uncertainty>40)=40;
        end
        DOB_uncertainty(isnan(DOB_uncertainty))=0;
        logPeatHg_uncertainty=normrnd(pixelPeatHglogP50,pixelPeatHglogsigma);
        PeatHg_uncertainty=exp(1).^logPeatHg_uncertainty;
        PeatHg_uncertainty(logPeatHg_uncertainty==0)=0;
        % PeatHg_uncertainty(PeatHg_uncertainty>800)=0;
        pixelPeatemission=BA_uncertainty*PeatCover_uncertainty/100*BD.*DOB_uncertainty*PeatCF*PeatHg_uncertainty*(1-Peat_Ash./100)./10^2;%km2/cm2*ng=10*12;
    end
end