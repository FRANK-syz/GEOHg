clear;
clc;
tic;
[Area,R]=readgeoraster("Data\WORLD\0.25x0.25base.tif"); 
continent=readgeoraster("Data\WORLD\continent025.tif");
Lc=readgeoraster("Data\MODIS_LC\2015\2015Lc025major.tif");
Lc3d=repmat(Lc,[1 1 10]); 
Lc4d=repmat(Lc,[1 1 12 10]);

[Forest,Shrub,Grass]=deal(Lc);
Forest(Lc>5)=0;Forests=repmat(Forest,[1 1 12 10]);
Shrub(Lc<6 | Lc>8)=0;Shrubs=repmat(Shrub,[1 1 12 10]);
Grass(Lc<9 | Lc==12 | Lc==13 | Lc==15 | Lc==17)=0;Grasses=repmat(Grass,[1 1 12 10]);

%% BA
% GFED5_month=read_GFED5(2010,2019,"Total");

load("Data\GFED5.1\GFED5month2002-2024.mat");
GFED5_month=single(BAmonth02_24(:,:,:,2010-2001:2019-2001));

% imagesc(BAmonth(:,:,7,14))
% GFED4s_2015=zeros(720,1440,12);
% for month=1:12
%     GFED4s_month=transpose(h5read("Data/GFED4s/GFED4.1s_2015.hdf5",strcat('/burned_area/',num2str(month,'%02d'),'/burned_fraction'))).*Area/10^6;%km2
%     GFED4s_2015(:,:,month)=GFED4s_month;
% end
% MCD64A1_2015=readgeoraster("Data\MCD64A1\BA2015.tif");
% GFED5_month(:,:,:,6)=MCD64A1_2015;
disp("BA done")

%% F
Tc=readTc("Data\MODIS_VCF");
Tc=Tc(:,:,1:10);
meanTc=mean(Tc,3);
AGB10_19=readgeoraster("Data\AGB\AGB10-19.tif");
LAI=zeros(720,1440,10);
LAIdir=dir("Data\LAI\LAIsets*");
for year=1:10
    currentLAIs=ncread(strcat(LAIdir(year).folder,"\",LAIdir(year).name),'lai');
    currentLAI=max(currentLAIs,[],3);
    LAI(:,:,year)=rot90(currentLAI);
end
SLA=readgeoraster("Data\SLA\Global_SLA.tif");
% testSLA=ncread("Data\SLA\Global_Maps_SLA.nc","variable");
% SLA=imresize(transpose(testSLA(:,:,4)),[720 1440],'bilinear');

Leaf_AGB=LAI./SLA*10;%kg/m2->t/ha
Leaf_AGB(isnan(Leaf_AGB) | Leaf_AGB==Inf)=0;

meanLeaf_AGB=mean(Leaf_AGB,3,"omitmissing");meanAGB=mean(AGB10_19,3,"omitmissing");
meanLeaf_AGB(continent==4)=0;meanAGB(continent==4)=0;
meanLeaf_AGB(meanTc>=10)=0;meanAGB(meanTc>=10)=0;
Leaf_AGB_series=meanLeaf_AGB(meanLeaf_AGB>0 & meanAGB>0 & meanAGB<10);
AGB_series=meanAGB(meanLeaf_AGB>0 & meanAGB>0 & meanAGB<10);
ratios=AGB_series./Leaf_AGB_series;
alpha=round(median(ratios(ratios>0)),2);
Leaf_AGB=Leaf_AGB.*alpha;

% Leaf_AGB(Leaf_AGB>AGB10_19)=AGB10_19(Leaf_AGB>AGB10_19);
% Leaf_AGB(repmat(Grass,[1 1 10])~=0)=AGB10_19(repmat(Grass,[1 1 10])~=0);
% figure(11)
% histogram(ratios,0:1:20)
% figure(12)
% scatter(Leaf_AGB_series*alpha,AGB_series,1);xlim([0 10]);ylim([0 10])

% ratio_Leaf_AGB=Leaf_AGB./AGB10_19.*100;         
% ratio_Leaf_AGB(ratio_Leaf_AGB>100)=0;
% imagesc(ratio_Leaf_AGB(:,:,1));colorbar

clear meanTc LAIdir currentLAI currentLAIs meanLeaf_AGB meanAGB AGB_series Leaf_AGB_series ratios ratio_Leaf_AGB year SLA LAI
%%
Stem_AGB=AGB10_19-Leaf_AGB;
Branch_AGB=Stem_AGB;%总木质13-23%
Branch_AGB(Lc3d==1)=Stem_AGB(Lc3d==1).*0.13;
Branch_AGB(Lc3d==2)=Stem_AGB(Lc3d==2).*0.28;
Branch_AGB(Lc3d==3)=Stem_AGB(Lc3d==3).*0.08;
Branch_AGB(Lc3d==4)=Stem_AGB(Lc3d==4).*0.23;
Branch_AGB(Lc3d==5)=Stem_AGB(Lc3d==5).*0.18;

%15树干=树皮
Bolewood_AGB=(Stem_AGB-Branch_AGB).*0.85;
Bark_AGB=(Stem_AGB-Branch_AGB).*0.15;

%农田
Leaf_AGB(Lc3d==12)=0.45.*AGB10_19(Lc3d==12);
Branch_AGB(Lc3d==12)=0.55.*AGB10_19(Lc3d==12);

clear Stem_AGB
disp("AGB done")

%% CF
load("Data\PKU_NDVI\NDVI01-22.mat");
maxNDVI=max(NDVI,[],[3 4]);maxNDVIs=repmat(maxNDVI,[1 1 12 10]);
minNDVI=min(NDVI,[],[3 4]);minNDVIs=repmat(minNDVI,[1 1 12 10]);
NDVI=NDVI(:,:,:,10:19);
VCI=(NDVI-minNDVIs)./(maxNDVIs-minNDVIs).*100;
VCI(isnan(VCI))=0;

mcf=VCI;
mcf(mcf>0 & mcf<=(100/6))=0.33;
mcf(mcf>(100/6) & mcf<=(200/6))=0.5;
mcf(mcf>(200/6) & mcf<=(300/6))=1;
mcf(mcf>(300/6) & mcf<=(400/6))=2;
mcf(mcf>(400/6) & mcf<=(500/6))=4;
mcf(mcf>(500/6) & mcf<=(100))=5;

% Tc=readTc("Data\MODIS_VCF");
Tc=To4D(Tc);

CFforest=(1-1/exp(1)).^mcf;
CFforest(Tc<=60)=0;

CFwood=exp(-0.013.*Tc);
CFwood(Tc<=40 | Tc>60)=0;

CFgrass=-2.13.*VCI./100+1.38;
CFgrass(VCI==0)=0;
CFgrass(CFgrass<0.44 & CFgrass~=0)=0.44;
CFgrass(CFgrass>0.98)=0.98;
CFgrass(Tc>40)=0;

CF=CFforest+CFwood+CFgrass;
CF(Lc4d==12)=0.9;%crop

clear NDVI VCI mcf maxNDVI maxNDVIs minNDVI minNDVIs Tc
%%
ratio=[3/8;3/9;3/8;5/9;5/8;6/9];

Canopy_AGB=Leaf_AGB+Branch_AGB;
CFleaf1=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(1)));
CFleaf1(Lc4d~=1)=0;
CFleaf2=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(2)));
CFleaf2(Lc4d~=2)=0;
CFleaf3=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(3)));
CFleaf3(Lc4d~=3)=0;
CFleaf4=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(4)));
CFleaf4(Lc4d~=4)=0;
CFleaf5=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(5)));
CFleaf5(Lc4d~=5)=0;
CFleaf6=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(6)));%shrub
CFleaf6(Shrubs==0)=0;

CFleaf7=CFgrass;
CFleaf7(Grasses==0)=0;%grass

CFleaf=CFleaf1+CFleaf2+CFleaf3+CFleaf4+CFleaf5+CFleaf6+CFleaf7;
CFleaf(CFleaf>1)=1;CFleaf(To4D(Leaf_AGB)==0)=0;

CFbranch=(To4D(Canopy_AGB).*CF-To4D(Leaf_AGB).*CFleaf)./To4D(Branch_AGB);
CFbranch(Forests==0 & Shrubs==0)=0;
CFbranch(To4D(Branch_AGB)==0)=0;

CFleaf(Lc4d==12)=0.9;
CFbranch(Lc4d==12)=0.9;

CFbolewood=CFbranch;CFbark=CFbranch;
CFbolewood(Forests==0)=0;CFbark(Forests==0)=0;

clear CFleaf1 CFleaf2 CFleaf3 CFleaf4 CFleaf5 CFleaf6 CFleaf7
disp("CF done")

%% MC
Foliar_Mc=ncread("Data\FoliarHg\GlobalFoliarHg_MachineLearning.nc","FoliarHg_concentrations");
[Leaf_Mc,Branch_Mc,Bolewood_Mc,Bark_Mc]=deal(double(zeros(720,1440)));

Leaf_Mc(Lc==1)=27;
Leaf_Mc(Lc==2)=53;
Leaf_Mc(Lc==3)=49;
Leaf_Mc(Lc==4)=41;
Leaf_Mc(Lc==5)=25;
Leaf_Mc(Lc==6)=19;
Leaf_Mc(Lc==7)=19;
Leaf_Mc(Lc==8)=19;
Leaf_Mc(Lc==9)=20;
Leaf_Mc(Lc==10)=20;
Leaf_Mc(Lc==11)=20;
Leaf_Mc(Lc==12)=28;
Leaf_Mc(Lc==14)=20;
Leaf_Mc(Lc==16)=20;

Branch_Mc(Lc==1)=16;
Branch_Mc(Lc==2)=12;
Branch_Mc(Lc==3)=19;
Branch_Mc(Lc==4)=12;
Branch_Mc(Lc==5)=5;
Branch_Mc(Lc==6)=6;
Branch_Mc(Lc==7)=6;
Branch_Mc(Lc==8)=6;
Branch_Mc(Lc==12)=21;

Bolewood_Mc(Lc==1)=2;
Bolewood_Mc(Lc==2)=2;
Bolewood_Mc(Lc==3)=3;
Bolewood_Mc(Lc==4)=2;
Bolewood_Mc(Lc==5)=3;

Bark_Mc(Lc==1)=17;
Bark_Mc(Lc==2)=4;
Bark_Mc(Lc==3)=19;
Bark_Mc(Lc==4)=9;
Bark_Mc(Lc==5)=13;

Leaf_Mc(Foliar_Mc~=0)=Foliar_Mc(Foliar_Mc~=0);%填补空值
%% Ash
[Leaf_Ash,Branch_Ash,Bark_Ash,Bolewood_Ash]=deal(double(zeros(720,1440)));

Leaf_Ash(Forest~=0)=9.1;%树叶
Leaf_Ash(Shrub~=0)=8.75;%灌木叶
Leaf_Ash(Grass~=0)=8.4;%草叶
Leaf_Ash(Lc==12)=5.9;

Branch_Ash(Forest | Shrub~=0)=4.6;%树枝
Branch_Ash(Lc==12)=5.9;%农

Bolewood_Ash(Forest~=0)=1.8;    
Bark_Ash(Forest~=0)=4.1;

%% Litter
litterHgflux=readgeoraster("Data\Litterfall\LitterHgflux.tif");
litterfallCF=0.8;
litter_Ash=9.1;
%% Peat
PeatCover=readgeoraster("Data\Peatland\Peatland025.tif");
PeatCover(PeatCover<5)=0;
BD=0.1;%g/cm3
load("Data\Soil\DOB.mat");
PeatCF=0.8;
Peat_Ash=4.9;
PeatHg=readgeoraster("Data\Peatland\PeatHg.tif");

%Peat burned ratio
PeatBA=repmat(PeatCover,[1 1 12 10]).*GFED5_month/100;
sumPeatBA=mean(squeeze(sum(PeatBA,[1 2 3])))
sumGlobalBA=mean(squeeze(sum(GFED5_month,[1 2 3])))
sumPeatBA/sumGlobalBA*100
%% Calc

Quick_emission=GFED5_month.*To4D(AGB10_19).*CF.*20./10^4;
Quick_series=squeeze(sum(Quick_emission,[1 2 3]))./10^3;%t

Leaf_emission=GFED5_month.*To4D(Leaf_AGB).*CFleaf.*repmat(Leaf_Mc,[1 1 12 10]).*(1-repmat(Leaf_Ash,[1 1 12 10])./100)./10^4;%km2*t/ha*ng/g=10^6m2*10^6g/10^4m2*10^-9=10^-1g
Branch_emission=GFED5_month.*To4D(Branch_AGB).*CFbranch.*repmat(Branch_Mc,[1 1 12 10]).*(1-repmat(Branch_Ash,[1 1 12 10])./100)./10^4;%kg
Bolewood_emission=GFED5_month.*To4D(Bolewood_AGB).*CFbolewood.*repmat(Bolewood_Mc,[1 1 12 10]).*(1-repmat(Bolewood_Ash,[1 1 12 10])./100)./10^4;%kg
Bark_emission=GFED5_month.*To4D(Bark_AGB).*CFbark.*repmat(Bark_Mc,[1 1 12 10]).*(1-repmat(Bark_Ash,[1 1 12 10])./100)./10^4;%kg
Litterfall_emission=GFED5_month.*litterfallCF.*repmat(litterHgflux,[1 1 12 10])*(1-litter_Ash/100)./10^3;%kg
Peat_emission=GFED5_month.*repmat(PeatCover,[1 1 12 10])./100.*BD.*DOB.*PeatCF.*repmat(PeatHg,[1 1 12 10]).*(1-Peat_Ash./100)./10^2;%km2/cm2*ng=10*12
Norm_emission=Leaf_emission+Branch_emission+Bolewood_emission+Bark_emission;
Total_emission=Norm_emission+Litterfall_emission+Peat_emission;

%Interannual
Leaf_emission_series=squeeze(sum(Leaf_emission,[1 2 3]))./10^3;%t
Branch_emission_series=squeeze(sum(Branch_emission,[1 2 3]))./10^3;%t
Bolewood_emission_series=squeeze(sum(Bolewood_emission,[1 2 3]))./10^3;%t
Bark_emission_series=squeeze(sum(Bark_emission,[1 2 3]))./10^3;%t
Litterfall_emission_series=squeeze(sum(Litterfall_emission,[1 2 3]))./10^3;%t
Peat_emission_series=squeeze(sum(Peat_emission,[1 2 3]))./10^3;%t
Norm_emission_series=squeeze(sum(Norm_emission,[1 2 3]))./10^3;%t
Total_emission_series=squeeze(sum(Total_emission,[1 2 3]))./10^3;%t

%Trend analysis
[Z_score, p_val, ~] = MK_Trend_Test(Total_emission_series, 0.05);
% fprintf('Z统计量: %.4f\nP值: %.4e\n趋势: %s\n', Z_score, p_val, trend_direction);

%MultiYear Average
Annual_Total_emission=squeeze(sum(Total_emission,3));
MultiYear_Norm_emission=mean(squeeze(sum(Norm_emission,3)),3);
MultiYear_Total_emission=mean(Annual_Total_emission,3);

% MultiYear_Leaf_flux=mean(squeeze(sum(Leaf_emission,3)),3)./Area*10^9;
% MultiYear_Branch_flux=mean(squeeze(sum(Branch_emission,3)),3)./Area*10^9;
% MultiYear_Bolewood_flux=mean(squeeze(sum(Bolewood_emission,3)),3)./Area*10^9;
% MultiYear_Bark_flux=mean(squeeze(sum(Bark_emission,3)),3)./Area*10^9;
% MultiYear_Litterfall_flux=mean(squeeze(sum(Litterfall_emission,3)),3)./Area*10^9;
% MultiYear_Peat_flux=mean(squeeze(sum(Peat_emission,3)),3)./Area*10^9;
MultiYear_Total_flux=MultiYear_Total_emission./Area*10^9; %ug/m2

% figure(1)
% imagesc(MultiYear_Total_flux)
%% Dry Matter burned
% Leaf_DM=GFED5_month.*To4D(Leaf_AGB).*CFleaf*10^2;%km2*t/ha=10^6m2*10^3kg/10^4m2=10^2t
% Branch_DM=GFED5_month.*To4D(Branch_AGB).*CFbranch.*10^2;%t
% Bolewood_DM=GFED5_month.*To4D(Bolewood_AGB).*CFbolewood.*10^2;%t
% Bark_DM=GFED5_month.*To4D(Bark_AGB).*CFbark.*10^2;%t
% DMB=Leaf_DM+Branch_DM+Bolewood_DM+Bark_DM;
% Total_DMB=squeeze(sum(DMB,[1 2 3]))/10^9*0.49;%PgC
% MeanSOC=readgeoraster("Data\Soil\SOC\MeanSOC.tif");
% PeatC=GFED5_month.*repmat(PeatCover,[1 1 12 10])./100.*DOB.*PeatCF.*MeanSOC;%km2*cm*t/ha=10^6*10^-2*/10^4=1t
% PeatC=squeeze(sum(PeatC,[1 2 3]))/10^6;%TgC

disp("calculation done")
%%
%Horizontal,Zonal
lon=-179.875:0.25:179.875;
lat=89.875:-0.25:-89.875;
Horizontal=[transpose(lon),transpose(squeeze(sum(MultiYear_Total_emission,1)))];
Zonal=[transpose(lat),sum(MultiYear_Total_emission,2)];%kg'

%Tropic,Temperate,Boreal
[Tropical_emission,Boreal_emission]=deal(squeeze(sum(Total_emission,[2 3]))/10^3);
Tropical_emission=sum(Tropical_emission(241:480,:),1);
Boreal_emission=sum(Boreal_emission(1:160,:),1);
Temperate_emission=transpose(Total_emission_series)-Tropical_emission-Boreal_emission;

[Tropical_Peat,Boreal_Peat]=deal(squeeze(sum(Peat_emission,[2 3]))/10^3);
Tropical_Peat=sum(Tropical_Peat(241:480,:),1);
Boreal_Peat=sum(Boreal_Peat(1:160,:),1);
Temperate_Peat=transpose(Peat_emission_series)-Tropical_Peat-Boreal_Peat;

%% 10a trend
% Trendamount=zeros(720,1440);
% Trendrate=zeros(720,1440);
% for i=1:720
%     for j=1:1440
%         Y=squeeze(Annual_Total_emission(i,j,:));
%         X=transpose(1:1:10);
%         b=regress(Y,[ones(size(X)),X]);
%         annual_change=b(2);
%         annual_rate=b(2)/(b(1)+b(2))*100;
% 
%         if b(1)+b(2)==0 || b(2)==0
%             Trendamount(i,j)=0;
%             Trendrate(i,j)=0;
%         else
%            Trendamount(i,j)=annual_change; 
%            Trendrate(i,j)=annual_rate; 
%         end
%     end
% end
% Trendrate(Trendrate>20)=20;
% Trendrate(Trendrate<-20)=-20;
% h = fspecial('average', [3 3]);  % 3×3均值滤波器
% Trendamount=imfilter(Trendamount,h);
% Trendrate=imfilter(Trendrate,h);
% 
% geotiffwrite("OBBHg results\Trendamount.tif",Trendamount,R)
% geotiffwrite("OBBHg results\Trendrate.tif",Trendrate,R)

%%
%plant tissue
Tissue_list=[mean(Leaf_emission_series),mean(Branch_emission_series),mean(Bolewood_emission_series),mean(Bark_emission_series),mean(Litterfall_emission_series),mean(Peat_emission_series)];

%land type
[ENF,EBF,DNF,DBF,MF,CS,OS,WS,SV,GS,PW,Crop,Mosaic,Barren]=deal(MultiYear_Total_emission);
ENF(Lc~=1)=0;
EBF(Lc~=2)=0;
DNF(Lc~=3)=0;
DBF(Lc~=4)=0;
MF(Lc~=5)=0;
CS(Lc~=6)=0;
OS(Lc~=7)=0;
WS(Lc~=8)=0;
SV(Lc~=9)=0;
GS(Lc~=10)=0;
PW(Lc~=11)=0;
Crop(Lc~=12)=0;
Mosaic(Lc~=14)=0;
Barren(Lc~=16)=0;
Landtype_list=transpose([sum(ENF(:)),sum(EBF(:)),sum(DNF(:)),sum(DBF(:)),sum(MF(:)),sum(CS(:)),sum(OS(:)),sum(WS(:)), ...
    sum(SV(:)),sum(GS(:)),sum(PW(:)),sum(Crop(:)),sum(Mosaic(:)),sum(Barren(:))])/10^3;

%land area
[Forest_emission,Shrub_emission,Grass_emission,Agri_emission]=deal(Norm_emission);
Forest_emission(Forests==0)=0;Forest_emission_series=squeeze(sum(Forest_emission,[1 2 3]))/10^3;
Shrub_emission(Shrubs==0)=0;Shrub_emission_series=squeeze(sum(Shrub_emission,[1 2 3]))/10^3;
Grass_emission(Grasses==0)=0;Grass_emission_series=squeeze(sum(Grass_emission,[1 2 3]))/10^3;
Agri_emission(Lc4d~=12)=0;Agri_emission_series=squeeze(sum(Agri_emission,[1 2 3]))/10^3;

[ForestLitter_emission,ShrubLitter_emission,GrassLitter_emission]=deal(Litterfall_emission);
ForestLitter_emission(Forests==0)=0;Forestlitter_series=squeeze(sum(ForestLitter_emission,[1 2 3]))/10^3;
ShrubLitter_emission(Shrubs==0)=0;Shrublitter_series=squeeze(sum(ShrubLitter_emission,[1 2 3]))/10^3;
GrassLitter_emission(Grasses==0 & Lc4d~=12)=0;Grasslitter_series=squeeze(sum(GrassLitter_emission,[1 2 3]))/10^3;

[ForestPeat_emission,ShrubPeat_emission,GrassPeat_emission]=deal(Peat_emission);
ForestPeat_emission(Forests==0)=0;ForestPeat_series=squeeze(sum(ForestPeat_emission,[1 2 3]))/10^3;
ShrubPeat_emission(Shrubs==0)=0;ShrubPeat_series=squeeze(sum(ShrubPeat_emission,[1 2 3]))/10^3;
GrassPeat_emission(Grasses==0 & Lc4d~=12)=0;GrassPeat_series=squeeze(sum(GrassPeat_emission,[1 2 3]))/10^3;
Veg_list=[Forest_emission_series,Forestlitter_series,ForestPeat_series,Shrub_emission_series,Shrublitter_series,ShrubPeat_series,Grass_emission_series,Grasslitter_series,GrassPeat_series,Agri_emission_series];

%continent
[AS_emission,NA_emission,EU_emission,AF_emission,SA_emission,OA_emission]=deal(Annual_Total_emission);
continent3d=repmat(continent,[1 1 10]);
AS_emission(continent3d~=1)=0;
NA_emission(continent3d~=2)=0;
EU_emission(continent3d~=3)=0;
AF_emission(continent3d~=4)=0;
SA_emission(continent3d~=5)=0;
OA_emission(continent3d~=6)=0;
Continent_list=transpose([squeeze(sum(AS_emission,[1 2])),squeeze(sum(NA_emission,[1 2])),squeeze(sum(EU_emission,[1 2])), ...
    squeeze(sum(AF_emission,[1 2])),squeeze(sum(SA_emission,[1 2])),squeeze(sum(OA_emission,[1 2]))])./10^3;
Continent_list=mean(Continent_list,2);

%month variation
Monthly_emission=mean(Total_emission,4);
Global_month=squeeze(sum(Monthly_emission,[1 2])./10^3);
[AS_month,NA_month,AF_month,SA_month]=deal(Monthly_emission);
continent3d=repmat(continent,[1 1 12]);
AS_month(continent3d~=1)=0;AS_month=squeeze(sum(AS_month,[1 2])./10^3);
NA_month(continent3d~=2)=0;NA_month=squeeze(sum(NA_month,[1 2])./10^3);
AF_month(continent3d~=4)=0;AF_month=squeeze(sum(AF_month,[1 2])./10^3);
SA_month(continent3d~=5)=0;SA_month=squeeze(sum(SA_month,[1 2])./10^3);
Month_list=[Global_month,AF_month,AS_month,NA_month,SA_month];
clear AF_month AS_month NA_month SA_month

%GFEDregions
GFEDregion=readgeoraster("Data\WORLD\GFEDregions.tif");
GFEDregion3d=repmat(GFEDregion,[1 1 10]);
[BONA,TENA,CEAM,NHSA,SHSA,EURO,MIDE,NHAF,SHAF,BOAS,CEAS,SEAS,EQAS,AUST]=deal(Annual_Total_emission);
BONA(GFEDregion3d~=1)=0;
TENA(GFEDregion3d~=2)=0;
CEAM(GFEDregion3d~=3)=0;
NHSA(GFEDregion3d~=4)=0;
SHSA(GFEDregion3d~=5)=0;
EURO(GFEDregion3d~=6)=0;
MIDE(GFEDregion3d~=7)=0;
NHAF(GFEDregion3d~=8)=0;
SHAF(GFEDregion3d~=9)=0;
BOAS(GFEDregion3d~=10)=0;
CEAS(GFEDregion3d~=11)=0;
SEAS(GFEDregion3d~=12)=0;
EQAS(GFEDregion3d~=13)=0;
AUST(GFEDregion3d~=14)=0;
GFEDregion_list=[squeeze(sum(BONA,[1 2])),squeeze(sum(TENA,[1 2])),squeeze(sum(CEAM,[1 2])),squeeze(sum(NHSA,[1 2])), ...
    squeeze(sum(SHSA,[1 2])),squeeze(sum(EURO,[1 2])),squeeze(sum(MIDE,[1 2])),squeeze(sum(NHAF,[1 2])),squeeze(sum(SHAF,[1 2])), ...
    squeeze(sum(BOAS,[1 2])),squeeze(sum(CEAS,[1 2])),squeeze(sum(SEAS,[1 2])),squeeze(sum(EQAS,[1 2])),squeeze(sum(AUST,[1 2]))]./10^3;
meanGFEDregion=mean(GFEDregion_list,1);

%%
% SESA season
Global_month=Monthly_emission;Global_month=squeeze(sum(Global_month,[1 2])/10^3);
Global_season=[sum(Global_month(3:5)),sum(Global_month(6:8)),sum(Global_month(9:11)),sum(Global_month(1:2))+Global_month(12)];
Forest_Global_month=Monthly_emission;Forest_Global_month(repmat(Forest,[1 1 12])==0)=0;Forest_Global_month=squeeze(sum(Forest_Global_month,[1 2])/10^3);
Forest_Global_season=[sum(Forest_Global_month(3:5)),sum(Forest_Global_month(6:8)),sum(Forest_Global_month(9:11)),sum(Forest_Global_month(1:2))+Forest_Global_month(12)];
Crop_Global_month=Monthly_emission;Crop_Global_month(repmat(Lc,[1 1 12])~=12)=0;Crop_Global_month=squeeze(sum(Crop_Global_month,[1 2])/10^3);
Crop_Global_season=[sum(Crop_Global_month(3:5)),sum(Crop_Global_month(6:8)),sum(Crop_Global_month(9:11)),sum(Crop_Global_month(1:2))+Crop_Global_month(12)];

SESA_month=Monthly_emission;SESA_month(repmat(GFEDregion,[1 1 12])~=12)=0;SESA_month=squeeze(sum(SESA_month,[1 2])/10^3);
SESA_season=[sum(SESA_month(3:5)),sum(SESA_month(6:8)),sum(SESA_month(9:11)),sum(SESA_month(1:2))+SESA_month(12)];
Forest_SESA_month=Monthly_emission;Forest_SESA_month(repmat(GFEDregion,[1 1 12])~=12 | repmat(Forest,[1 1 12])==0)=0;Forest_SESA_month=squeeze(sum(Forest_SESA_month,[1 2])/10^3);
Forest_SESA_season=[sum(Forest_SESA_month(3:5)),sum(Forest_SESA_month(6:8)),sum(Forest_SESA_month(9:11)),sum(Forest_SESA_month(1:2))+Forest_SESA_month(12)];
Crop_SESA_month=Monthly_emission;Crop_SESA_month(repmat(GFEDregion,[1 1 12])~=12 | repmat(Lc,[1 1 12])~=12)=0;Crop_SESA_month=squeeze(sum(Crop_SESA_month,[1 2])/10^3);
Crop_SESA_season=[sum(Crop_SESA_month(3:5)),sum(Crop_SESA_month(6:8)),sum(Crop_SESA_month(9:11)),sum(Crop_SESA_month(1:2))+Crop_SESA_month(12)];

Season_list=transpose([Global_season;Forest_Global_season;Crop_Global_season;SESA_season;Forest_SESA_season;Crop_SESA_season]);
%% Events
Brazil025=readgeoraster("Data\WORLD\Brazil.tif");
Indonesia025=readgeoraster("Data\WORLD\Indonesia.tif");
Australia025=readgeoraster("Data\WORLD\Australia.tif");
China025=readgeoraster("Data\WORLD\China.tif");
India025=readgeoraster("Data\WORLD\India.tif");
Myanmar025=readgeoraster("Data\WORLD\Myanmar.tif");
Russia025=readgeoraster("Data\WORLD\Russia.tif");
Canada025=readgeoraster("Data\WORLD\Canada.tif");
America025=readgeoraster("Data\WORLD\USA.tif");
Mexico025=readgeoraster("Data\WORLD\Mexico.tif");

Brazil_2010=Annual_Total_emission(:,:,1);Brazil_2010(Brazil025~=0)=0;Brazil_2010=sum(Brazil_2010,[1 2])/10^3;%2010
Indonesia_2015=Annual_Total_emission(:,:,6);Indonesia_2015(Indonesia025~=0)=0;Indonesia_2015=sum(Indonesia_2015,[1 2])/10^3;%2015
Australia_2012=Annual_Total_emission(:,:,3);Australia_2012(Australia025~=0)=0;Australia_2012=sum(Australia_2012,[1 2])/10^3;%2012
Australia_2019=Annual_Total_emission(:,:,10);Australia_2019(Australia025~=0)=0;Australia_2019=sum(Australia_2019,[1 2])/10^3;%2019
China_2014=Annual_Total_emission(:,:,5);China_2014(China025~=0)=0;China_2014=sum(China_2014,[1 2])/10^3;
Myanmar_2010=Annual_Total_emission(:,:,1);Myanmar_2010(Myanmar025~=0)=0;Myanmar_2010=sum(Myanmar_2010,[1 2])/10^3;
Russia_2012=Annual_Total_emission(:,:,3);Russia_2012(Russia025~=0)=0;Russia_2012=sum(Russia_2012,[1 2])/10^3;
Canada_2013=Annual_Total_emission(:,:,4);Canada_2013(Canada025~=0)=0;Canada_2013=sum(Canada_2013,[1 2])/10^3;
America_2015=Annual_Total_emission(:,:,6);America_2015(America025~=0)=0;America_2015=sum(America_2015,[1 2])/10^3;
Mexico_2011=Annual_Total_emission(:,:,2);Mexico_2011(Mexico025~=0)=0;Mexico_2011=sum(Mexico_2011,[1 2])/10^3;

Brazil_average=MultiYear_Total_emission;Brazil_average(Brazil025~=0)=0;Brazil_average=sum(Brazil_average,[1 2])/10^3;
Indonesia_average=MultiYear_Total_emission;Indonesia_average(Indonesia025~=0)=0;Indonesia_average=sum(Indonesia_average,[1 2])/10^3;
Australia_average=MultiYear_Total_emission;Australia_average(Australia025~=0)=0;Australia_average=sum(Australia_average,[1 2])/10^3;
China_average=MultiYear_Total_emission;China_average(China025~=0)=0;China_average=sum(China_average,[1 2])/10^3;
India_average=MultiYear_Total_emission;India_average(India025~=0)=0;India_average=sum(India_average,[1 2])/10^3;
Myanmar_average=MultiYear_Total_emission;Myanmar_average(Myanmar025~=0)=0;Myanmar_average=sum(Myanmar_average,[1 2])/10^3;
Russia_average=MultiYear_Total_emission;Russia_average(Russia025~=0)=0;Russia_average=sum(Russia_average,[1 2])/10^3;
Canada_average=MultiYear_Total_emission;Canada_average(Canada025~=0)=0;Canada_average=sum(Canada_average,[1 2])/10^3;
America_average=MultiYear_Total_emission;America_average(America025~=0)=0;America_average=sum(America_average,[1 2])/10^3;
Mexico_average=MultiYear_Total_emission;Mexico_average(Mexico025~=0)=0;Mexico_average=sum(Mexico_average,[1 2])/10^3;

Eventlist=[Brazil_2010,Brazil_average;
    Indonesia_2015,Indonesia_average;
    Australia_2012,Australia_average;
    Australia_2019,Australia_average;
    China_2014,China_average;
    Myanmar_2010,Myanmar_average;
    Russia_2012,Russia_average;
    Canada_2013,Canada_average;
    America_2015,America_average;
    Mexico_2011,Mexico_average];
toc;
%% Output
% imagesc(Annual_Total_emission(:,:,1)-MultiYear_Total_emission);
% clim([-10,10])

% GFED4s_2015=Annual_Total_emission(:,:,6);
% save("OBBHg results\GFED4Semission2015.mat","GFED4s_2015");
% MCD64_2015=Annual_Total_emission(:,:,6);
% save("OBBHg results\MCD64emission2015.mat","MCD64_2015");

% geotiffwrite("OBBHg results\MultiYear_Total_emission.tif",MultiYear_Total_emission,R);
% geotiffwrite("OBBHg results\MultiYear_Total_flux.tif",MultiYear_Total_flux,R);
% geotiffwrite("OBBHg results\MultiYear_Leaf_flux.tif",MultiYear_Leaf_flux,R);
% geotiffwrite("OBBHg results\MultiYear_Branch_flux.tif",MultiYear_Branch_flux,R);
% geotiffwrite("OBBHg results\MultiYear_Bolewood_flux.tif",MultiYear_Bolewood_flux,R);
% geotiffwrite("OBBHg results\MultiYear_Bark_flux.tif",MultiYear_Bark_flux,R);
% geotiffwrite("OBBHg results\MultiYear_Litterfall_flux.tif",MultiYear_Litterfall_flux,R);
% geotiffwrite("OBBHg results\MultiYear_Peat_flux.tif",MultiYear_Peat_flux,R);
% geotiffwrite("OBBHg results\Annual_Total_emission.tif",Annual_Total_emission,R);

% meanLeaf_AGB=mean(Leaf_AGB,3);
% meanBranch_AGB=mean(Branch_AGB,3);meanBranch_AGB(Forest==0 & Shrub==0 & Crop==0)=0;
% meanBolewood_AGB=mean(Bolewood_AGB,3);
% meanBark_AGB=mean(Bark_AGB,3);
% geotiffwrite("OBBHg results\Leaf_AGB.tif",meanLeaf_AGB,R);
% geotiffwrite("OBBHg results\Branch_AGB.tif",meanBranch_AGB,R);
% geotiffwrite("OBBHg results\Bolewood_AGB.tif",meanBolewood_AGB,R);
% geotiffwrite("OBBHg results\Bark_AGB.tif",meanBark_AGB,R);

% CF(CF==0)=NaN;meanCF=mean(CF,[3,4],'omitmissing');meanCF(Lc==12)=0;
% CFleaf(CFleaf==0)=NaN;meanCFleaf=mean(CFleaf,[3,4],'omitmissing');meanCFleaf(Lc==12)=0;
% CFbranch(CFbranch==0)=NaN;meanCFbranch=mean(CFbranch,[3,4],'omitmissing');meanCFbranch(Lc==12)=0;
% CFbolewood(CFbolewood==0)=NaN;meanCFbolewood=mean(CFbolewood,[3,4],'omitmissing');meanCFbolewood(Lc==12)=0;
% CFbark(CFbark==0)=NaN;meanCFbark=mean(CFbark,[3,4],'omitmissing');meanCFbark(Lc==12)=0;
% geotiffwrite("OBBHg results\CF.tif",meanCF,R);
% geotiffwrite("OBBHg results\CFleaf.tif",meanCFleaf,R);
% geotiffwrite("OBBHg results\CFbranch.tif",meanCFbranch,R);
% geotiffwrite("OBBHg results\CFbolewood.tif",meanCFbolewood,R);
% geotiffwrite("OBBHg results\CFbark.tif",meanCFbark,R);
%%
function BAmonth=read_GFED5(start_year,end_year,fire_type)
    duration=end_year-start_year+1;
    BAmonth=zeros(720,1440,12,duration);
    for currentYear=start_year:end_year
        year=currentYear-start_year+1;
        Annual=zeros(720,1440);
        dir_GFED5=dir("Data\GFEDv5\BA\BA"+(currentYear)+"*.nc"); 
        for month=1:12
            currentMonth=rot90(ncread("Data\GFEDv5\BA\"+dir_GFED5(month).name,fire_type));
            Annual=Annual+currentMonth;
            BAmonth(:,:,month,year)=currentMonth;
        end
    end
end

function Tc=readTc(dirname)
    Tcdir=dir(fullfile(dirname,"*Tc025*"));
    Tc=zeros(720,1440,length(Tcdir));
    for i=1:10
        currentTc=readgeoraster(fullfile(Tcdir(i).folder,Tcdir(i).name));
        Tc(:,:,i)=currentTc;
    end
end

function Mat4D=To4D(Mat3D)
    Mat4D=permute(repmat(Mat3D,[1 1 1 12]),[1 2 4 3]);
end

function [Z, p_value, trend] = MK_Trend_Test(data, alpha)
% 输入:
%   data: 待检验的一维时间序列数据向量
%   alpha: 显著性水平 (默认0.05)
% 输出:
%   Z: 标准化统计量Z值
%   p_value: 检验的P值
%   trend: 趋势判断 ('increasing', 'decreasing', 'no trend')

    n = length(data);
    S = 0; % 初始化统计量S

    % 计算符号函数总和 (核心的双重循环)
    for i = 1:n-1
        for j = i+1:n
            S = S + sign(data(j) - data(i));
        end
    end

    % 计算方差 (考虑可能存在的结值，此处为简化版)
    VarS = (n*(n-1)*(2*n+5)) / 18;

    % 计算标准化统计量Z
    if S > 0
        Z = (S - 1) / sqrt(VarS);
    elseif S < 0
        Z = (S + 1) / sqrt(VarS);
    else
        Z = 0;
    end

    % 计算P值 (双尾检验)
    p_value = 2 * (1 - normcdf(abs(Z))); % normcdf为标准正态分布累积函数

    % 判断趋势
    if abs(Z) >= norminv(1 - alpha/2) % 与显著性水平对应的临界值比较
        if Z > 0
            trend = 'increasing';
        else
            trend = 'decreasing';
        end
    else
        trend = 'no trend';
    end
end