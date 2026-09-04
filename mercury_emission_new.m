clear;
clc;
tic;

startyear=2010;
endyear=2024;
duration=endyear-startyear+1;

[Area,R]=readgeoraster("Data\WORLD\0.25x0.25base.tif"); 
continent=readgeoraster("Data\WORLD\continent025.tif");
Lc=readgeoraster("Data\MODIS_LC\2015\2015Lc025major.tif");
Lc3d=repmat(Lc,[1 1 duration]); 
Lc4d=repmat(Lc,[1 1 12 duration]);

[Forest,Shrub,Grass]=deal(Lc);
Forest(Lc>5)=0;Forests=repmat(Forest,[1 1 12 duration]);
Shrub(Lc<6 | Lc>8)=0;Shrubs=repmat(Shrub,[1 1 12 duration]);
Grass(Lc<9 | Lc==12 | Lc==13 | Lc==15 | Lc==17)=0;Grasses=repmat(Grass,[1 1 12 duration]);

toc;
disp("初始化完成")
%% BA
% BAinfo=ncinfo("Data\GFED5.1\GFED5.1_ecosystem_2002.nc");
% files=dir("Data\GFED5.1\GFED5.1*");
% BAmonth02_24=zeros(720,1440,12,size(files,1));
% for i=1:size(files,1)
%     GFED5month=ncread(files(i).folder+"\"+files(i).name,"burned_area"); 
%     for currentmonth=1:12
%         currentmonthBA=GFED5month(:,:,currentmonth);
%         currentmonthBA=rot90(currentmonthBA)/10^6; %km2
%         BAmonth02_24(:,:,currentmonth,i)=currentmonthBA;
%     end
% end

load("Data\GFED5.1\GFED5month2002-2024.mat");
BAmonth=single(BAmonth02_24(:,:,:,startyear-2001:endyear-2001));
% imagesc(BAmonth(:,:,7,14))

toc;
disp("BA完成")
clear BAmonth02_24

%% 
% BorealBA=sum(BAmonth(1:140,:,:,:),3)/10^5;%mha
% BorealBA_series=squeeze(sum(BorealBA,[1 2]));
% plot(startyear:1:endyear,BorealBA_series)
% xticks(startyear:1:endyear)
% clear BAmonth02_24 BorealBA

%% AGB
load("Data\GLASS\GLASSmonthNDVI.mat")
load("Data\GLASS\GLASSmonthLAI.mat")
NDVI=GLASSNDVI02_24(:,:,:,startyear-2001:endyear-2001);
AMNDVI=squeeze(max(NDVI,[],3));

Tcdir=dir(fullfile("Data\MODIS_VCF","*Tc025.tif"));
Tc10_20=single(zeros(720,1440,length(Tcdir)));
for i=1:length(Tcdir)
    currentTc=readgeoraster(fullfile(Tcdir(i).folder,Tcdir(i).name));
    Tc10_20(:,:,i)=currentTc;
end
Tc10_24=single(zeros(720, 1440, length(Tcdir)+4));
Tc10_24(:, :, 1:length(Tcdir)) = Tc10_20;
Tc10_24(:, :, length(Tcdir)+1:length(Tcdir)+4) = repmat(Tc10_20(:, :, length(Tcdir)), [1, 1, 4]);
Tc=Tc10_24(:,:,startyear-2009:endyear-2009);
Tc4d=To4D(Tc);
meanTc=mean(Tc,3);

AGB2010=readgeoraster("Data\AGB\AGBC\2010AGB025.tif");
AGB2010(AMNDVI(:,:,9)==0)=0;
AGB=AGB2010.*(AMNDVI./AMNDVI(:,:,9));
AGB(AGB>500)=0;AGB(isnan(AGB))=0;

LAI=GLASSmonthLAI(:,:,:,startyear-2001:endyear-2001);
AMLAI=squeeze(max(LAI,[],3,'omitmissing'));
SLA=single(readgeoraster("Data\SLA\Global_SLA.tif"));

Leaf_AGB=AMLAI./SLA*10;%t/ha
Leaf_AGB(isnan(Leaf_AGB) | Leaf_AGB==Inf)=0;

meanLeaf_AGB=mean(Leaf_AGB,3,"omitmissing");meanAGB=mean(AGB,3,"omitmissing");
meanLeaf_AGB(continent==4)=0;meanAGB(continent==4)=0;
meanLeaf_AGB(meanTc>=10)=0;meanAGB(meanTc>=10)=0;
Leaf_AGB_series=meanLeaf_AGB(meanLeaf_AGB>0 & meanAGB>0 & meanAGB<10);
AGB_series=meanAGB(meanLeaf_AGB>0 & meanAGB>0 & meanAGB<10);
ratios=AGB_series./Leaf_AGB_series;
alpha=round(median(ratios(ratios>0)),2);
Leaf_AGB=Leaf_AGB.*alpha;

Stem_AGB=AGB-Leaf_AGB;
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
Leaf_AGB(Lc3d==12)=0.45.*AGB(Lc3d==12);
Branch_AGB(Lc3d==12)=0.55.*AGB(Lc3d==12);

Leaf_AGB_4d=To4D(Leaf_AGB);
Branch_AGB_4d=To4D(Branch_AGB);
Bolewood_AGB_4d=To4D(Bolewood_AGB);
Bark_AGB_4d=To4D(Bark_AGB);
Canopy_AGB_4d=Leaf_AGB_4d+Branch_AGB_4d;

clear SLA LAI AMLAI AGB2010 Stem_AGB GLASSmonthLAI GLASSNDVI02_24 Tc10_24 Tc10_20 AGB_series Leaf_AGB_series currentTc Tcdir ratios meanAGB meanLeaf_AGB meanTc
toc;
disp("AGB done")

%% CF
maxNDVI=max(NDVI,[],[3 4]);
minNDVI=min(NDVI,[],[3 4]);
VCI=(NDVI-minNDVI)./(maxNDVI-minNDVI)*100;
VCI(isnan(VCI))=0;

thresholds=[0,100/6,200/6,300/6,400/6,500/6,100];
values=[0,0.33,0.5,1,2,4,5];
mcf=zeros(size(VCI));
for i = 1:length(thresholds)-1
    mask=VCI>thresholds(i)&VCI<=thresholds(i+1);
    mcf(mask)=values(i+1);
end

CFforest=(1-1/exp(1)).^mcf;
CFforest(Tc4d<=60)=0;

CFwood=exp(-0.013.*Tc4d);
CFwood(Tc4d<=40 | Tc4d>60)=0;

CFgrass = max(0.44,min(0.98,-2.13.*VCI./100+1.38));
CFgrass(VCI==0|Tc4d>40)=0;

CF=CFforest+CFwood+CFgrass;
CF(Lc4d==12)=0.9;%crop

clear maxNDVI minNDVI NDVI VCI mcf Tc Tc4d mask thresholds values
%%
ratio=[3/8;3/9;3/8;5/9;5/8;6/9];

CFleaf_Forest=single(zeros(720,1440,12,duration));
for i=1:5
    CFleafi=Canopy_AGB_4d.*CF./(Leaf_AGB_4d+Branch_AGB_4d*ratio(i));
    CFleafi(Lc4d~=i)=0;
    CFleaf_Forest=CFleaf_Forest+CFleafi;
end
CFleaf_Shrub=Canopy_AGB_4d.*CF./(Leaf_AGB_4d+Branch_AGB_4d*ratio(6));
CFleaf_Shrub(Shrubs==0)=0;
CFleaf_Grass=CFgrass;
CFleaf_Grass(Grasses==0)=0;%grass

CFleaf=CFleaf_Forest+CFleaf_Shrub+CFleaf_Grass;
CFleaf(CFleaf>1)=1;CFleaf(Leaf_AGB_4d==0)=0;

% CFleaf1=Canopy_AGB_4d.*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(1)));
% CFleaf1(Lc4d~=1)=0;
% CFleaf2=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(2)));
% CFleaf2(Lc4d~=2)=0;
% CFleaf3=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(3)));
% CFleaf3(Lc4d~=3)=0;
% CFleaf4=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(4)));
% CFleaf4(Lc4d~=4)=0;
% CFleaf5=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(5)));
% CFleaf5(Lc4d~=5)=0;
% CFleaf6=To4D(Canopy_AGB).*CF./(To4D(Leaf_AGB)+To4D(Branch_AGB*ratio(6)));%shrub
% CFleaf6(Shrubs==0)=0;

CFbranch=(Canopy_AGB_4d.*CF-Leaf_AGB_4d.*CFleaf)./Branch_AGB_4d;
CFbranch(Forests==0 & Shrubs==0)=0;
CFbranch(Branch_AGB_4d==0)=0;

CFleaf(Lc4d==12)=0.9;
CFbranch(Lc4d==12)=0.9;

CFbolewood=CFbranch;CFbark=CFbranch;
CFbolewood(Forests==0)=0;CFbark(Forests==0)=0;

clear CFleaf_Forest CFleaf_Shrub CFleaf_Grass CFleafi
disp("CF done")

%% MC
Foliar_Mc=single(ncread("Data\FoliarHg\GlobalFoliarHg_MachineLearning.nc","FoliarHg_concentrations"));
[Leaf_Mc,Branch_Mc,Bolewood_Mc,Bark_Mc]=deal(single(zeros(720,1440)));

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
[Leaf_Ash,Branch_Ash,Bark_Ash,Bolewood_Ash]=deal(single(zeros(720,1440)));

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

% DOB=readDOB("Data\Soil\ERA5_SW_2010-2024.nc",startyear,endyear);
load("Data\Soil\DOB_2010-2024.mat");
DOB=single(DOB(:,:,:,startyear-2009:endyear-2009));
PeatCF=0.8;
Peat_Ash=4.9;
PeatHg=readgeoraster("Data\Peatland\PeatHg.tif");

%% Calc
Leaf_emission=BAmonth.*Leaf_AGB_4d.*CFleaf.*repmat(Leaf_Mc,[1 1 12 duration]).*(1-repmat(Leaf_Ash,[1 1 12 duration])./100)./10^4;%km2*t/ha*ng/g=10^6m2*10^6g/10^4m2*10^-9=10^-1g
Branch_emission=BAmonth.*Branch_AGB_4d.*CFbranch.*repmat(Branch_Mc,[1 1 12 duration]).*(1-repmat(Branch_Ash,[1 1 12 duration])./100)./10^4;%kg
Bolewood_emission=BAmonth.*Bolewood_AGB_4d.*CFbolewood.*repmat(Bolewood_Mc,[1 1 12 duration]).*(1-repmat(Bolewood_Ash,[1 1 12 duration])./100)./10^4;%kg
Bark_emission=BAmonth.*Bark_AGB_4d.*CFbark.*repmat(Bark_Mc,[1 1 12 duration]).*(1-repmat(Bark_Ash,[1 1 12 duration])./100)./10^4;%kg
Litterfall_emission=BAmonth.*litterfallCF.*repmat(litterHgflux,[1 1 12 duration])*(1-litter_Ash/100)./10^3;%kg
Peat_emission=BAmonth.*repmat(PeatCover,[1 1 12 duration])./100.*BD.*DOB.*PeatCF.*repmat(PeatHg,[1 1 12 duration]).*(1-Peat_Ash./100)./10^2;%km2/cm2*ng=10*12
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

%MultiYear Average
Annual_Total_emission=squeeze(sum(Total_emission,3));
MultiYear_Norm_emission=mean(squeeze(sum(Norm_emission,3)),3);
MultiYear_Total_emission=mean(Annual_Total_emission,3);
%%
%Soil emission
Soil_year_emission=single(zeros(720,1440,duration));
Legacy_dir=dir("Data\Topsoil_Hgemissions\*Mg*");
for i=1:size(Legacy_dir,1)
    currentSoil=single(readgeoraster(Legacy_dir(i).folder+"\"+Legacy_dir(i).name));
    Soil_year_emission(:,:,i)=currentSoil*1000;
end
Soil_year_emission=Soil_year_emission(:,:,startyear-2009:endyear-2009);
Soil_year_emission(PeatCover>=5)=0;%去泥炭地
% Soil_year_emission=Soil_year_emission*0;disp("无土壤排放！！！")
BAyear=squeeze(sum(BAmonth,3));
Soil_emission_series=squeeze(sum(Soil_year_emission,[1 2]));
Soil_emission=To4D(Soil_year_emission).*BAmonth./To4D(BAyear);%排放比例系数
Soil_emission(To4D(BAyear)==0)=0;%去掉不重合部分

All_emission=Total_emission+Soil_emission;%%没有月
Annual_All_emission=Annual_Total_emission+Soil_year_emission;
MultiYear_All_emission=mean(Annual_All_emission,3);
All_emission_series=squeeze(sum(All_emission,[1 2 3]))./10^3;%t

clear currentSoil Soil_dir

%% Plot
%Boreal emission
Boreal_BA=BAmonth(1:160,:,:,:);%55° 1:140 50° 1:160
Boreal_BA_series=squeeze(sum(Boreal_BA,[1 2 3]))./10^3;
Boreal_emission=All_emission(1:160,:,:,:);
Boreal_emission_series=squeeze(sum(Boreal_emission,[1 2 3]))./10^3;
Boreal_soil_emission=Soil_emission(1:160,:,:,:);
Boreal_soil_emission_series=squeeze(sum(Boreal_soil_emission,[1 2 3]))./10^3;%t
BONA_emission=Boreal_emission(:,1:600,:,:);%-30% 1:600
BONA_emission_series=squeeze(sum(BONA_emission,[1 2 3]))./10^3;
BOAS_emission=Boreal_emission(:,601:1440,:,:);%-30% 600:1440
BOAS_emission_series=squeeze(sum(BOAS_emission,[1 2 3]))./10^3;

%Canada emission
Canada025=readgeoraster("Data\WORLD\Canada.tif");
Canada025_4d=repmat(Canada025,[1 1 12 duration]);
Canada_BA=BAmonth;
Canada_BA(Canada025_4d~=0)=0;
Canada_BA_series=squeeze(sum(Canada_BA,[1 2 3]))./10^3;
Canada_emission=All_emission;
Canada_emission(Canada025_4d~=0)=0;
Canada_emission_series=squeeze(sum(Canada_emission,[1 2 3]))./10^3;

%Russia emission
Russia025=readgeoraster("Data\WORLD\Russia.tif");
Russia025_4d=repmat(Russia025,[1 1 12 duration]);
Russia_BA=BAmonth;
Russia_BA(Russia025_4d~=0)=0;
Russia_BA_series=squeeze(sum(Russia_BA,[1 2 3]))./10^3;
Russia_emission=All_emission;
Russia_emission(Russia025_4d~=0)=0;
Russia_emission_series=squeeze(sum(Russia_emission,[1 2 3]))./10^3;

%% Analyze
%Horizontal,Zonal
lon=-179.875:0.25:179.875;
lat=89.875:-0.25:-89.875;
% Horizontal=[transpose(lon),transpose(squeeze(sum(MultiYear_All_emission,1)))];
Zonal=[transpose(lat),sum(MultiYear_All_emission,2)];%kg
% Zonal225=sum(MultiYear_All_emission,2);
% Zonal225=transpose(mean(reshape(Zonal225,[8 90]),1)*8);
Zonal_Boreal_front=mean(sum(Annual_All_emission(1:160,:,6:10),2),3);
Zonal_Boreal_front=transpose(sum(reshape(Zonal_Boreal_front,[4 40]),1));
Zonal_Boreal_behind=mean(sum(Annual_All_emission(1:160,:,11:15),2),3);
Zonal_Boreal_behind=transpose(sum(reshape(Zonal_Boreal_behind,[4 40]),1));
Zonal_Boreal2021=sum(Annual_All_emission(1:160,:,12),2);
Zonal_Boreal2021=transpose(sum(reshape(Zonal_Boreal2021,[4 40]),1));
Zonal_Boreal2023=sum(Annual_All_emission(1:160,:,14),2);
Zonal_Boreal2023=transpose(sum(reshape(Zonal_Boreal2023,[4 40]),1));
%%
% figure(1)
% imagesc(Annual_All_emission(1:160,:,12))
% clim([0 10])

% figure(2)
% Anomaly_2021=Annual_All_emission(:,:,12)-mean(Annual_All_emission(:,:,1:10),3);
% imagesc(Anomaly_2021)
% clim([-5 5])

figure(3)
hold on
% plot(startyear:1:endyear,Total_emission_series)
plot(startyear:1:endyear,Boreal_emission_series,"r",startyear:1:endyear,BONA_emission_series,"g",startyear:1:endyear,BOAS_emission_series,"b",startyear:1:endyear,Boreal_soil_emission_series,"--",startyear:1:endyear,Boreal_emission_series-Boreal_soil_emission_series,"--")  
xticks(startyear:1:endyear)
hold off

figure(4)
hold on
plot(startyear:1:endyear,Canada_emission_series,"g")
plot(startyear:1:endyear,Russia_emission_series,"b")
xticks(startyear:1:endyear)
hold off

figure(5)
hold on
plot(startyear:1:endyear,Boreal_BA_series,"r")
plot(startyear:1:endyear,Canada_BA_series,"g")
plot(startyear:1:endyear,Russia_BA_series,"b")
xticks(startyear:1:endyear)
hold off
toc;

clear Canada025 Canada025_4d Russia025 Russia025_4d
%%
function Mat4D=To4D(Mat3D)
    Mat4D=permute(repmat(Mat3D,[1 1 1 12]),[1 2 4 3]);
end

function DOB=readDOB(file,startyear,endyear)
    duration=endyear-startyear+1;
    PeatCover=readgeoraster("Data\Peatland\Peatland025.tif");
    PeatCover=repmat(PeatCover,[1,1,12*duration]);

    SW0_7=ncread(file,"swvl1");
    SW7_28=ncread(file,"swvl2");
    SW0_7=SW0_7(:,:,((startyear-2010)*12+1):((endyear-2010)*12+12));
    SW7_28=SW7_28(:,:,((startyear-2010)*12+1):((endyear-2010)*12+12));

    SW1=zeros(720,1440,12*duration);
    SW2=zeros(720,1440,12*duration);
    for i=1:12*duration
        SW1(:,:,i)=transpose(imresize(SW0_7(:,:,i),[1440 720],"nearest"));
        SW2(:,:,i)=transpose(imresize(SW7_28(:,:,i),[1440 720],"nearest"));
    end
    SW=0.2.*SW1+0.8.*SW2;
    SW=[SW(:,721:1440,:),SW(:,1:720,:)];

    TropicalDOB=-51.*SW./0.45+57;% the volumetric soil water content = SM/0.45 m3 m−3       0-57cm
    TropicalDOB(TropicalDOB<0)=0;
    TropicalDOB(PeatCover<5)=0;TropicalDOB(1:200,:,:)=0;TropicalDOB(541:720,:,:)=0;% 30N/S

    BorealSW=SW;BorealSW(SW<=0.2)=NaN;
    BorealDOB=236.*(42*10.*BorealSW./(10.*BorealSW+1)-28).^(-0.38)-100; %soil gravimetric moisture content >67% <90%  10-100cm
    BorealDOB(PeatCover<5)=0;
    BorealDOB(isnan(BorealDOB))=0;BorealDOB(BorealDOB>40)=40;
    BorealDOB(201:540,:,:)=0;
    DOB=TropicalDOB+BorealDOB;
    DOB=reshape(DOB,[720 1440 12 duration]);
end