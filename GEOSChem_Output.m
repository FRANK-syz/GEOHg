clear;
clc;
tic;

%%
DryDep_info=ncinfo("Data\GEOSChem\BaseGFED\GEOSChem.DryDep.20170101_0000z.nc4");
HgEmis_info=ncinfo("Data\GEOSChem\BaseGFED\GEOSChem.MercuryEmis.20170101_0000z.nc4");
HgConc_info=ncinfo("Data\GEOSChem\BaseGFED\GEOSChem.SpeciesConc.20170101_0000z.nc4");
WetLossConv_info=ncinfo("Data\GEOSChem\BaseGFED\GEOSChem.WetLossConv.20170101_0000z.nc4");
WetLossLS_info=ncinfo("Data\GEOSChem\BaseGFED\GEOSChem.WetLossLS.20170101_0000z.nc4");
Area225=single(transpose(ncread("Data\GEOSChem\BaseGFED\GEOSChem.WetLossLS.20170101_0000z.nc4","AREA")));
[~,R]=readgeoraster("Data\WORLD\2x2.5zero.tif");

%% 读取各方案[DD:ug/m2,HE:kg,CC:ng/m3,WD:ug/m2]
Basedir="Data\GEOSChem\BaseGFED";
[BaseHg0dep,BaseHE,BaseCC,BaseHg2dep]=getGCdir(Basedir,2015,2023);
GFEDoffdir="Data\GEOSChem\NoEF";
[GFEDoffHg0dep,GFEDoffHE,GFEDoffCC,GFEDoffHg2dep]=getGCdir(GFEDoffdir,2015,2024);
MBOBBdir="Data\GEOSChem\MBOBBnew";
[MBOBBHg0dep,MBOBBHE,MBOBBCC,MBOBBHg2dep]=getGCdir(MBOBBdir,2015,2024);

% Borealdir="Data\GEOSChem\BorealOBB";
% [BorealDD,BorealHE,BorealCC,BorealWD]=getGCdir(Borealdir,2015,2024);
toc;

%% 异常分析
Hg0depfromGFED =BaseHg0dep -GFEDoffHg0dep(:,:,1:9);
Hg0depfromMBOBB=MBOBBHg0dep-GFEDoffHg0dep;
AnnualHg0depfromGFED =squeeze(sum(Hg0depfromGFED,3));
AnnualHg0depfromMBOBB=squeeze(sum(Hg0depfromMBOBB,3));

CCfromGFED=BaseCC-GFEDoffCC(:,:,1:9);
CCfromMBOBB=MBOBBCC-GFEDoffCC;
MeanCCfromGFED=squeeze(mean(CCfromGFED,3));
MeanCCfromMBOBB=squeeze(mean(CCfromMBOBB,3));

Hg2depfromGFED=BaseHg2dep-GFEDoffHg2dep(:,:,1:9);
Hg2depfromMBOBB=MBOBBHg2dep-GFEDoffHg2dep;
AnnualHg2depfromGFED=squeeze(sum(Hg2depfromGFED,3));
AnnualHg2depfromMBOBB=squeeze(sum(Hg2depfromMBOBB,3));

Anomaly2021Hg0dep=AnnualHg0depfromMBOBB(:,:,7)-mean(AnnualHg0depfromMBOBB(:,:,1:5),3);
Anomaly2023Hg0dep=AnnualHg0depfromMBOBB(:,:,9)-mean(AnnualHg0depfromMBOBB(:,:,1:5),3);
Anomaly2020sHg0dep=mean(AnnualHg0depfromMBOBB(:,:,6:9),3)-mean(AnnualHg0depfromMBOBB(:,:,1:5),3);
Anomaly2021CC=MeanCCfromMBOBB(:,:,7)-mean(MeanCCfromMBOBB(:,:,1:5),3);
Anomaly2023CC=MeanCCfromMBOBB(:,:,9)-mean(MeanCCfromMBOBB(:,:,1:5),3);
Anomaly2020sCC=mean(MeanCCfromMBOBB(:,:,6:9),3)-mean(MeanCCfromMBOBB(:,:,1:5),3);
Anomaly2021Hg2dep=AnnualHg2depfromMBOBB(:,:,7)-mean(AnnualHg2depfromMBOBB(:,:,1:5),3);
Anomaly2023Hg2dep=AnnualHg2depfromMBOBB(:,:,9)-mean(AnnualHg2depfromMBOBB(:,:,1:5),3);
Anomaly2020sHg2dep=mean(AnnualHg2depfromMBOBB(:,:,6:9),3)-mean(AnnualHg2depfromMBOBB(:,:,1:5),3);

%多年平均
Multiyear_Hg0depGFED=squeeze(mean(AnnualHg0depfromGFED,3));
Multiyear_GFEDCC=squeeze(mean(MeanCCfromGFED,3));
Multiyear_Hg2depGFED=squeeze(mean(AnnualHg2depfromGFED,3));
Multiyear_Hg0depMBOBB=squeeze(mean(AnnualHg0depfromMBOBB,3));
Multiyear_MBOBBCC=squeeze(mean(MeanCCfromMBOBB,3));
Multiyear_Hg2depMBOBB=squeeze(mean(AnnualHg2depfromMBOBB,3));
Multiyear_TDGFED =Multiyear_Hg0depGFED+Multiyear_Hg2depGFED;
Multiyear_TDMBOBB=Multiyear_Hg0depMBOBB+Multiyear_Hg2depMBOBB;

%% 趋势（注意geos-chem输出纬度不规则，需用Area225做面积加权）
GFEDHgTD_series=squeeze(sum((Hg0depfromGFED+Hg2depfromGFED).*Area225,[1 2 3]))/10^12;
MBOBBHgTD_series=squeeze(sum((Hg0depfromMBOBB+Hg2depfromMBOBB).*Area225,[1 2 3]))/10^12;
MBOBBHg0dep_series=squeeze(sum((Hg0depfromMBOBB).*Area225,[1 2 3]))/10^12;
MBOBBHg2dep_series=squeeze(sum((Hg2depfromMBOBB).*Area225,[1 2 3]))/10^12;

GFED4_Emis_series=squeeze(sum(BaseHE,[1 2 3])/10^3);%kg->t
MBOBB_Emis_series=squeeze(sum(MBOBBHE,[1 2 3])/10^3);%kg->t
Hg0depGFED_series=squeeze(sum(AnnualHg0depfromGFED.*Area225,[1 2]))/10^12;
Hg2depGFED_series=squeeze(sum(AnnualHg2depfromGFED.*Area225,[1 2]))/10^12;
Hg0depMBOBB_series=squeeze(sum(AnnualHg0depfromMBOBB.*Area225,[1 2]))/10^12;
Hg2depMBOBB_series=squeeze(sum(AnnualHg2depfromMBOBB.*Area225,[1 2]))/10^12;

%Horizontal,Zonal
lon=-178.75:2.5:178.75;
lat=89:-2:-89;

ZonalHg0dep=sum(Multiyear_Hg0depMBOBB.*Area225,2)/10^9;%kg
ZonalHg0dep=Resample91to90(ZonalHg0dep);
ZonalHg2dep=sum(Multiyear_Hg2depMBOBB.*Area225,2)/10^9;%kg
ZonalHg2dep=Resample91to90(ZonalHg2dep);
ZonalDep=[transpose(lat),ZonalHg0dep,ZonalHg2dep,ZonalHg0dep+ZonalHg2dep];%kg

%% 贡献
DDratioMBOBB=mean(AnnualDDfromMBOBB(:,:,1:5),3)./mean(AnnualMBOBBDD(:,:,1:5),3)*100;
CCratioMBOBB=mean(MeanCCfromMBOBB(:,:,1:5),3)./mean(MeanMBOBBCC(:,:,1:5),3)*100;
WDratioMBOBB=mean(AnnualWDfromMBOBB(:,:,1:5),3)./mean(AnnualMBOBBWD(:,:,1:5),3)*100;

DDratioMBOBB2023=AnnualDDfromMBOBB(:,:,9)./AnnualMBOBBDD(:,:,9)*100;
CCratioMBOBB2023=MeanCCfromMBOBB(:,:,9)./MeanMBOBBCC(:,:,9)*100;
WDratioMBOBB2023=AnnualWDfromMBOBB(:,:,9)./AnnualMBOBBWD(:,:,9)*100;
%% 北方
NorthEmis_series=squeeze(sum(MBOBBHE(1:20,:,:,:),[1 2 3])/10^3);

NorthZonalDep=ZonalDep(1:20,:);
NorthHg0dep=AnnualHg0depfromMBOBB(1:20,:,:);
NorthHg0dep_series=squeeze(sum(NorthHg0dep.*Area225(1:20,:,:),[1 2]))/10^12;
NorthHg2dep=AnnualHg2depfromMBOBB(1:20,:,:);
NorthHg2dep_series=squeeze(sum(NorthHg2dep.*Area225(1:20,:,:),[1 2]))/10^12;
NorthHgTD_series=NorthHg0dep_series+NorthHg2dep_series;
%% 画图
Draw1=0;
Draw2=0;
Draw3=1;
if Draw1==1
    year=9;
    figure(1);imagesc(AnnualDDfromGFED(:,:,year));colorbar
    clim([0,10])
    figure(2);imagesc(AnnualDDfromMBOBB(:,:,year));colorbar
    clim([0,10])
    figure(3);imagesc(MeanCCfromGFED(:,:,year));colorbar
    clim([0,0.2])
    figure(4);imagesc(MeanCCfromMBOBB(:,:,year));colorbar
    clim([0,0.2])
    figure(5);imagesc(AnnualWDfromGFED(:,:,year));colorbar
    clim([0,0.4])
    figure(6);imagesc(AnnualWDfromMBOBB(:,:,year));colorbar
    clim([0,0.4])
end

if Draw2==1
    figure(7);imagesc(Anomaly2021DD);colorbar
    figure(8);imagesc(Anomaly2021CC);colorbar
    figure(9);imagesc(Anomaly2021WD);colorbar
    clim([-0.1 0.1])
    figure(10);imagesc(Anomaly2023DD);colorbar
    figure(11);imagesc(Anomaly2023CC);colorbar
    figure(12);imagesc(Anomaly2023WD);colorbar
    clim([-0.1 0.1])
end

if Draw3==1
    figure(13);imagesc(DDratioMBOBB);colorbar
    clim([0 25])
    figure(14);imagesc(CCratioMBOBB);colorbar
    clim([0 25])
    figure(15);imagesc(WDratioMBOBB);colorbar
    clim([0 10])
    figure(16);imagesc(DDratioMBOBB2023);colorbar
    clim([0 25])
    figure(17);imagesc(CCratioMBOBB2023);colorbar
    clim([0 25])
    figure(18);imagesc(WDratioMBOBB2023);colorbar
    clim([0 10])
end

%% 输出
% Hg0dep=Resample91to90(Multiyear_Hg0depGFED);
% CC=Resample91to90(Multiyear_GFEDCC);
% Hg2dep=Resample91to90(Multiyear_Hg2depGFED);
% geotiffwrite("BorealHg results\GFED_Hg0dep.tif",Hg0dep,R)
% geotiffwrite("BorealHg results\GFED_CC.tif",CC,R)
% geotiffwrite("BorealHg results\GFED_Hg2dep.tif",Hg2dep,R)
% 
% Hg0dep=Resample91to90(Multiyear_Hg0depMBOBB);
% CC=Resample91to90(Multiyear_MBOBBCC);
% Hg2dep=Resample91to90(Multiyear_Hg2depMBOBB);
% geotiffwrite("BorealHg results\MBOBB_Hg0dep.tif",Hg0dep,R)
% geotiffwrite("BorealHg results\MBOBB_CC.tif",CC,R)
% geotiffwrite("BorealHg results\MBOBB_Hg2dep.tif",Hg2dep,R)
% 
% Anomaly2021Hg0dep=Resample91to90(Anomaly2021Hg0dep);
% Anomaly2021CC=Resample91to90(Anomaly2021CC);
% Anomaly2021Hg2dep=Resample91to90(Anomaly2021Hg2dep);
% Anomaly2023Hg0dep=Resample91to90(Anomaly2023Hg0dep);
% Anomaly2023CC=Resample91to90(Anomaly2023CC);
% Anomaly2023Hg2dep=Resample91to90(Anomaly2023Hg2dep);
% Anomaly2020sHg0dep=Resample91to90(Anomaly2020sHg0dep);
% Anomaly2020sCC=Resample91to90(Anomaly2020sCC);
% Anomaly2020sHg2dep=Resample91to90(Anomaly2020sHg2dep);
% geotiffwrite("BorealHg results\Anomaly2021Hg0dep.tif",Anomaly2021Hg0dep,R)
% geotiffwrite("BorealHg results\Anomaly2021CC.tif",Anomaly2021CC,R)
% geotiffwrite("BorealHg results\Anomaly2021Hg2dep.tif",Anomaly2021Hg2dep,R)
% geotiffwrite("BorealHg results\Anomaly2023Hg0dep.tif",Anomaly2023Hg0dep,R)
% geotiffwrite("BorealHg results\Anomaly2023CC.tif",Anomaly2023CC,R)
% geotiffwrite("BorealHg results\Anomaly2023Hg2dep.tif",Anomaly2023Hg2dep,R)
% geotiffwrite("BorealHg results\Anomaly2020sHg0dep.tif",Anomaly2020sHg0dep,R)
% geotiffwrite("BorealHg results\Anomaly2020sCC.tif",Anomaly2020sCC,R)
% geotiffwrite("BorealHg results\Anomaly2020sHg2dep.tif",Anomaly2020sHg2dep,R)
% 
% Hg0dep2021=Resample91to90(AnnualHg0depfromMBOBB(:,:,7));
% CC2021=Resample91to90(MeanCCfromMBOBB(:,:,7));
% Hg2dep2021=Resample91to90(AnnualHg2depfromMBOBB(:,:,7));
% Hg0dep2023=Resample91to90(AnnualHg0depfromMBOBB(:,:,9));
% CC2023=Resample91to90(MeanCCfromMBOBB(:,:,9));
% Hg2dep2023=Resample91to90(AnnualHg2depfromMBOBB(:,:,9));
% Hg0dep15_19=Resample91to90(mean(AnnualHg0depfromMBOBB(:,:,1:5),3));
% CC15_19=Resample91to90(mean(MeanCCfromMBOBB(:,:,1:5),3));
% Hg2dep15_19=Resample91to90(mean(AnnualHg2depfromMBOBB(:,:,1:5),3));
% geotiffwrite("BorealHg results\2021Hg0dep.tif",Hg0dep2021,R)
% geotiffwrite("BorealHg results\2021CC.tif",CC2021,R)
% geotiffwrite("BorealHg results\2021Hg2dep.tif",Hg2dep2021,R)
% geotiffwrite("BorealHg results\2023Hg0dep.tif",Hg0dep2023,R)
% geotiffwrite("BorealHg results\2023CC.tif",CC2023,R)
% geotiffwrite("BorealHg results\2023Hg2dep.tif",Hg2dep2023,R)
% geotiffwrite("BorealHg results\15-19Hg0dep.tif",Hg0dep15_19,R)
% geotiffwrite("BorealHg results\15-19CC.tif",CC15_19,R)
% geotiffwrite("BorealHg results\15-19Hg2dep.tif",Hg2dep15_19,R)

%%
function [DryDep,DryDepHg0,HgEmis,HgConc,WetDep]=getGCdir(targetdir,startyear,endyear)%收集各要素文件夹
    DryDepdir=dir(strcat(targetdir+"\*DryDep*"));% molec cm s
    HgEmisdir=dir(strcat(targetdir+"\*MercuryEmis*"));%kg s
    HgConcdir=dir(strcat(targetdir+"\*SpeciesCon*"));%mol mol dry
    WetConvdir=dir(strcat(targetdir+"\*WetLossConv*"));
    WetLossLSdir=dir(strcat(targetdir+"\*WetLossLS*"));%kg s

    [DryDep, DryDepHg0]=getDryDep(DryDepdir,startyear,endyear);
    HgEmis=getHgEmis(HgEmisdir,startyear,endyear);
    HgConc=getConc(HgConcdir,startyear,endyear);
    WetDep=getWetDep(WetConvdir,WetLossLSdir,startyear,endyear);
end

function [DryDep, DryDep_Hg0]=getDryDep(dir,startyear,endyear)
    num=length(dir);
    nyear=endyear-startyear+1;
    DryDep_Hg0=single(zeros(91,144,12,nyear));
    DryDep    =single(zeros(91,144,12,nyear));%总干沉降（Hg0+氧化汞）
    [years,months]=parseYM(dir);

    %氧化汞物种列表（与DryDep文件中的变量名一一对应）
    HgSp={'DryDep_Hg2ORGP','DryDep_Hg2ClP','DryDep_HgCl2','DryDep_HgOHOH',...
          'DryDep_HgOHBrO','DryDep_HgOHClO','DryDep_HgOHHO2','DryDep_HgOHNO2',...
          'DryDep_HgClOH','DryDep_HgClBr','DryDep_HgClBrO','DryDep_HgClClO',...
          'DryDep_HgClHO2','DryDep_HgClNO2','DryDep_HgBr2','DryDep_HgBrOH',...
          'DryDep_HgBrClO','DryDep_HgBrBrO','DryDep_HgBrHO2','DryDep_HgBrNO2'};
    nSp=numel(HgSp);

    for i=1:num
        if years(i)>=startyear && years(i)<=endyear
            im=months(i);
            iy=years(i)-startyear+1;
            month_days=eomday(years(i),im);%当月天数
            conv_factor=(6.02*10^23)/200.59/10^6/10^4/(month_days*24*3600);%molec cm-2 s-1->ug/m2/month
            fullpath=fullfile(dir(i).folder,dir(i).name);

            %Hg0
            DD0=rot90(single(ncread(fullpath,'DryDep_Hg0')))/conv_factor;
            DryDep_Hg0(:,:,im,iy)=DD0;
            %氧化汞求和（共用同一conv_factor，单位与Hg0一致）
            tmp=single(zeros(size(DD0)));
            for s=1:nSp
                tmp=tmp+rot90(single(ncread(fullpath,HgSp{s})));
            end
            DryDep(:,:,im,iy)=DD0+tmp/conv_factor;
            if im==12
                disp(['已处理',num2str(years(i)),'年干沉降']);
            end
        end
    end
    toc;
end

function HgEmis=getHgEmis(dir,startyear,endyear)
    num=length(dir);
    nyear=endyear-startyear+1;
    HgEmis=single(zeros(91,144,12,nyear));
    [years,months]=parseYM(dir);
    for i=1:num
        if years(i)<startyear || years(i)>endyear
            continue;
        end
        im=months(i);
        iy=years(i)-startyear+1;
        fullpath=fullfile(dir(i).folder,dir(i).name);
        %检查文件是否含EmisHg0biomass变量，无则该月留0
        info=ncinfo(fullpath);
        varNames=arrayfun(@(v) v.Name,info.Variables,'UniformOutput',false);
        if ~any(strcmp(varNames,'EmisHg0biomass'))
            continue;
        end
        month_days=eomday(years(i),im);
        HgEmis(:,:,im,iy)=rot90(single(ncread(fullpath,'EmisHg0biomass')))*month_days*24*3600;
        if im==12
            disp(['已处理',num2str(years(i)),'年汞排放'])
        end
    end
    toc;
end

function HgConc=getConc(dir,startyear,endyear)
    num=length(dir);
    nyear=endyear-startyear+1;
    HgConc=single(zeros(91,144,12,nyear));
    % mol/mol dry -> ng/m3
    CONV_FACTOR=200.59/22.4*10^12;
    [years,months]=parseYM(dir);
    for i=1:num
        if years(i)<startyear || years(i)>endyear
            continue;
        end
        im=months(i);
        iy=years(i)-startyear+1;
        ncfile=rot90(single(ncread(fullfile(dir(i).folder,dir(i).name),'SpeciesConcVV_Hg0')))*CONV_FACTOR;
        HgConc(:,:,im,iy)=ncfile(:,:,1);
        if im==12
            disp(['已处理',num2str(years(i)),'年汞浓度'])
        end
    end
    toc;
end

function WetDep=getWetDep(WCdir,WLdir,startyear,endyear)
    num=length(WCdir);
    nyear=endyear-startyear+1;
    WetDep=single(zeros(91,144,12,nyear));

    %变量名对所有文件一致，仅读取一次
    info_WC=ncinfo(fullfile(WCdir(1).folder,WCdir(1).name));
    info_WL=ncinfo(fullfile(WLdir(1).folder,WLdir(1).name));
    varNames_WC=arrayfun(@(j) info_WC.Variables(j).Name, 34:53, 'UniformOutput', false);
    varNames_WL=arrayfun(@(k) info_WL.Variables(k).Name, 14:33, 'UniformOutput', false);
    nWC=numel(varNames_WC);
    nWL=numel(varNames_WL);

    Area=rot90(ncread(fullfile(WCdir(1).folder,WCdir(1).name),'AREA'));
    CONST_BASE=3600*24*10^9./Area;%kg/s->ug/m2的常数部分

    [years,months]=parseYM(WCdir);
    for i=1:num
        if years(i)<startyear || years(i)>endyear
            continue;
        end
        month_days=eomday(years(i),months(i));
        im=months(i);
        iy=years(i)-startyear+1;
        conv_factor=CONST_BASE*month_days;%kg/s->ug/m2/month

        %== WetLossConv: 打开文件一次，读取所有变量，再关闭 ==
        ncid=netcdf.open(fullfile(WCdir(i).folder,WCdir(i).name),'NC_NOWRITE');
        acc_conv=single(zeros(144,91));
        for j=1:nWC
            varid=netcdf.inqVarID(ncid,varNames_WC{j});
            data=netcdf.getVar(ncid,varid,'single');
            acc_conv=acc_conv+squeeze(sum(data,3));
        end
        netcdf.close(ncid);
        conv_sum=rot90(acc_conv);

        %== WetLossLS: 打开文件一次，读取所有变量，再关闭 ==
        ncid=netcdf.open(fullfile(WLdir(i).folder,WLdir(i).name),'NC_NOWRITE');
        acc_ls=single(zeros(144,91));
        for k=1:nWL
            varid=netcdf.inqVarID(ncid,varNames_WL{k});
            data=netcdf.getVar(ncid,varid,'single');
            acc_ls=acc_ls+squeeze(sum(data,3));
        end
        netcdf.close(ncid);
        ls_sum=rot90(acc_ls);
        WetDep(:,:,im,iy)=(conv_sum+ls_sum).*conv_factor;
        if im==12
            disp(['已处理',num2str(years(i)),'年湿沉降'])
        end
    end
    toc;
end

function [years,months]=parseYM(filenames)
%从文件名解析年月，避免在各读取函数中重复正则解析
    n=numel(filenames);
    years=zeros(1,n);
    months=zeros(1,n);
    for i=1:n
        token=regexp(filenames(i).name,'\.(\d{4})(\d{2})\d{2}_','match');
        years(i)=str2double(token{1}(2:5));
        months(i)=str2double(token{1}(6:7));
    end
end

function Outputfile=Resample91to90(GCfile)
    if size(GCfile,1)==91
        Outputfile=(GCfile(1:end-1,:) + GCfile(2:end,:))/2;
    else
        Outputfile=GCfile;
    end
end