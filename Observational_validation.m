clear;
clc;
tic;

%%
[~,R]=readgeoraster("Data\WORLD\2x2.5zero.tif");
% DryDep_info=ncinfo("Data\GEOSChem\BaseGFED\GEOSChem.DryDep.20170101_0000z.nc4");

%% 读取各方案[DD:ug/m2,HE:kg,CC:ng/m3,WD:ug/m2]
Basedir="Data\GEOSChem\BaseGFED";
[BaseHg0dep,~,BaseCC]=getGCdir(Basedir,2019,2023);
MBOBBdir="Data\GEOSChem\MBOBB";
[MBOBBHg0dep,~,MBOBBCC]=getGCdir(MBOBBdir,2019,2024);%浓度垂直层！

%% 生长季提取
% GrowSeasonHg0dep_GFED=BaseHg0dep(:,:,4:9,:);%生长季4-9月，2019-2023年
GrowSeasonHg0dep_MBOBB=MBOBBHg0dep(:,:,4:9,:);
% GrowSeasonCC_GFED=BaseCC(:,:,4:9,:);
GrowSeasonCC_MBOBB=MBOBBCC(:,:,4:9,:);

%% 站点序列提取
Station=readtable("Data\Hg观测\MLN\MLN_Mercury_Litterfall.xls");
StationLat=table2array(Station(:,2));
StationLon=table2array(Station(:,3));
StationElv=table2array(Station(:,4));
StationHgC=table2array(Station(:,5:end));

%模拟点位提取
SimuHg0dep=zeros(size(StationHgC));
SimuCC    =zeros(size(StationHgC));
for num=1:length(StationElv)
    [X,Y]=latlon2grid(StationLat(num),StationLon(num));%改
    SimuHg0dep(num,:)=transpose(squeeze(sum(GrowSeasonHg0dep_MBOBB(X,Y,:,:),3)));
    SimuCC(num,:)=transpose(squeeze(mean(GrowSeasonCC_MBOBB(X,Y,:,:),3)));
end
% imagesc(squeeze(sum(GrowSeasonHg0dep_MBOBB(:,:,:,5),3)))
%%
figure(1)
hold on 
for i=1:length(StationLon)
    plot(SimuHg0dep(i,:))
end
hold off

figure(2)
hold on 
for i=1:length(StationLon)
    plot(SimuCC(i,:))
end
hold off

%% 2023年干沉降异常分布
Hg0dep2023=squeeze(sum(GrowSeasonHg0dep_MBOBB(:,:,:,5),3));
Hg0depNorm=(squeeze(sum(GrowSeasonHg0dep_MBOBB(:,:,:,:),[3,4]))-Hg0dep2023)/(size(GrowSeasonHg0dep_MBOBB,4)-1);

figure(3)
imagesc(Hg0dep2023-Hg0depNorm);
Diff2023Hg0dep=Resample91to90(Hg0dep2023-Hg0depNorm);
% geotiffwrite("BorealHg results\Diff2023Hg0dep.tif",Diff2023Hg0dep,R)

%%
function [Hg0dep,HgEmis,HgConc]=getGCdir(targetdir,startyear,endyear)%收集各要素文件夹
    DryDepdir=dir(strcat(targetdir+"\*DryDep*"));% molec cm s
    HgEmisdir=dir(strcat(targetdir+"\*MercuryEmis*"));%kg s
    HgConcdir=dir(strcat(targetdir+"\*SpeciesCon*"));%mol mol dry

    Hg0dep=getDryDep(DryDepdir,startyear,endyear);
    HgEmis=getHgEmis(HgEmisdir,startyear,endyear);
    HgConc=getConc(HgConcdir,startyear,endyear);
end

function [DryDep_Hg0]=getDryDep(dir,startyear,endyear)
    num=length(dir);
    nyear=endyear-startyear+1;
    DryDep_Hg0 =single(zeros(91,144,12,nyear));
    [years,months]=parseYM(dir);

    for i=1:num
        if years(i)>=startyear && years(i)<=endyear
            im=months(i);
            iy=years(i)-startyear+1;
            month_days=eomday(years(i),im);%当月天数
            conv_factor=(6.02*10^23)/200.59/10^6/10^4/(month_days*24*3600);%molec cm-2 s-1->ug/m2/month
            fullpath=fullfile(dir(i).folder,dir(i).name);

            %Hg0
            Hg0dep0=rot90(single(ncread(fullpath,'DryDep_Hg0')))/conv_factor;
            DryDep_Hg0(:,:,im,iy)=Hg0dep0;
            
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

function [lat_idx, lon_idx] = latlon2grid(lat,lon)
    lat_idx=round((92-lat)/2); %91-89°对应1，89-87对应2，0°对应46,，-89--91°对应91
    lon_idx=round((lon+181.25)/2.5);%180-178.75对应1
end
