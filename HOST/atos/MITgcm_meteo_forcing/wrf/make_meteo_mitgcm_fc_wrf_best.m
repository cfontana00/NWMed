function make_meteo_mitgcm_fc_wrf(date)

% Check OS
%
if isunix
    slash = '/';
else
    slash = '\';                                                                                               
end
%glob_path=strcat(slash,'home',slash,'innocenti',slash,'MITgcm_BFM',slash,...
%          'NWMed',slash,'HOST',slash,'hpc-fe1',slash,'MITgcm_meteo_forcing');
glob_path=strcat(slash,'home',slash,'itai',slash,'MITgcm_BFM',slash,...
          'NWMed',slash,'HOST',slash,'atos',slash,'MITgcm_meteo_forcing');
glob_path2=strcat(slash,'ec',slash,'res4',slash,'hpcperm',slash,...
          'itai');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bathymetry and Land sea Mask

fileinN=strcat(glob_path2,slash,'bathymetry',slash,'E5_2018.dtm');
fidN=netcdf.open(fileinN);
Ni=ncinfo(fileinN);
ND={Ni.Dimensions.Name};
NV={Ni.Variables.Name};
batN=ncread(fileinN,'DEPTH');
latN=ncread(fileinN,'LINES');
lonN=ncread(fileinN,'COLUMNS');
% figure;pcolor(batW');shading flat;colorbar
netcdf.close(fidN);
[LAN,LON]=meshgrid(latN,lonN);

% North side
fileinS=strcat(glob_path2,slash,'bathymetry',slash,'F5_2018.dtm');
fidS=netcdf.open(fileinS);
Si=ncinfo(fileinS);
SD={Si.Dimensions.Name};
SV={Si.Variables.Name};
batS=ncread(fileinS,'DEPTH');
latS=ncread(fileinS,'LINES');
lonS=ncread(fileinS,'COLUMNS');
% figure;pcolor(batE');shading flat;colorbar
netcdf.close(fidS);
[LAS,LOS]=meshgrid(latS,lonS);

% 1/128 degrees grid
nx=784; ny=336; % (98x8) x (48x7)  or (49x16) x (48x7)  or  (56x14) x (28x12)
res=1/128;      %      2 nodi              4 nodi                6 nodi

yc0=41.88125+(1/256); 
yce=yc0+(ny-1)*res;  % 42.1 - (336-308)/128;
xc0=6.0625+(1/256); 
xce=xc0+(nx-1)*res;     % 6 + 8/128; (altrimenti il bordo sud est viene tagliato prima della costa)
[Y,X]=meshgrid(yc0:res:yce,xc0:res:xce);

% zoom on the Italian coast of the EMODnet dataset
lonm=5; lonM=13; latm=41.5; latM=45;
lonzN=lonN(lonN>lonm & lonN<lonM);
latzN=latN(latN>latm & latN<latM);
lonzS=lonS(lonS>lonm & lonS<lonM);
latzS=latS(latS>latm & latS<latM);
% WARNING!!! 24 points overlap in the x direction -> start from 25
lonz=lonzN; %[lonzN; lonzS(25:end)];
latz=[latzS; latzN(25:end)]; %latzN; % =latzE;
[LAZ,LOZ]=meshgrid(latz,lonz);

nxN=size(lonzN,1);
nyN=size(latzN,1);
nxS=size(lonzS,1);
nyS=size(latzS,1);

batN01=batN(LON>lonm & LON<lonM & LAN>latm & LAN<latM);
batS01=batS(LOS>lonm & LOS<lonM & LAS>latm & LAS<latM);
batN02=reshape(batN01,nxN,nyN);
batS02=reshape(batS01,nxS,nyS);

bat=[batS02 batN02(:,25:end)];
% land points
bat(bat>=0.0)=0.0;
batNaN=bat;
batNaN(batNaN==0.0)=NaN;
% figure%('position',[10 10 1000 1000])
% pcolor(LOZ,LAZ,batNaN);colorbar;shading flat;colormap(jet);
% colorbar; colormap(jet); caxis([-3000 0])
% XLIM=get(gca,'Xlim'); YLIM=get(gca,'Ylim');
% maxdepth=-min(min(bat));
% text(LOZ(40,1),LAZ(1,70),['max depth= ',num2str(maxdepth,'%15.3f'), ' m'],...
%     'FontSize',30,'color','w');
% set(gcf,'renderer','zbuffer','PaperPositionMode','auto')
%print -dpng -r300 test_EMODnet_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2 - zoom on the Italian coast (at 1/128 deg.) of the EMODnet dataset
batz=interp2(LAZ,LOZ,bat,Y,X);
% land points
batz(batz>=0.0)=0.0;
% "flatten" the shallowest areas to the depth of the first level (1.5 m)
batz(batz>=-3.001 & batz<0.0)=-3.001;

% % manual adjustments
for i=2:length(batz(:,1))-1
    for j=2:length(batz(1,:))-1
        if ((batz(i,j) < 0.0) && (batz(i+1,j) == 0.0) && (batz(i-1,j) == 0.0) ...
               && (batz(i,j-1) == 0.0) && (batz(i,j+1) == 0.0))
           batz(i,j) = 0.0;
        end
    end
end

i=-9; j=27;

% Correzioni per griglia con vertice + 1/256

batz(7,148) = 0.0;
batz(9,148) = 0.0;
batz(11,146) = -3.001;
batz(22,145) = 0.0;
batz(56,149) = 0.0;
batz(57,150) = -3.001;
batz(103,199) = 0.0;
batz(136,215) = 0.0;
batz(137,217) = 0.0;
batz(162,233) = 0.0;
batz(163,233) = 0.0;
batz(212,246) = 0.0;
%batz(277,274) = 0.0;
batz(427,306) = 0.0;
batz(444,298) = 0.0;
batz(454,293) = 0.0;
batz(456,290) = 0.0;
batz(466,287) = 0.0;
batz(484,278) = 0.0;
batz(485,278) = 0.0;
batz(486,277) = -3.001;
batz(507,277) = 0.0;
batz(543,217) = 0.0;
batz(547,203) = 0.0;
batz(571,176) = 0.0;
batz(566,142) = -3.001;
batz(600,134) = 0.0;
batz(624,111) = 0.0;
batz(696,62) = 0.0;
batz(700,60) = 0.0;
batz(703,59) = 0.0;
batz(727,37) = 0.0;

batz(560,120) = 0.0;
batz(543,112) = 0.0;
batz(541,121) = 0.0;
batz(546,119) = 0.0;
batz(446,58) = 0.0;
batz(393,109) = -3.001;
batz(384,103) = 0.0;
batz(319,62) = 0.0;
batz(328,48) = -3.001;
batz(321,45) = 0.0;
batz(321,35) = -3.001;
batz(339,28) = 0.0;
batz(338,19) = 0.0;
batz(337,4) = 0.0;
batz(351,2) = 0.0;

batz(39,145) = -3.001;
batz(75,165) = -3.001;
batz(78,185) = 0.0;
batz(103,198) = 0.0;
batz(113,209) = 0.0;
batz(161,232) = 0.0;
batz(161,233) = 0.0;
batz(170,236) = 0.0;
batz(404,310) = -60.;
batz(426,307) = 0.0;
batz(438,302) = 0.0;
batz(454,292) = 0.0;
batz(458,289) = -3.001;
batz(461,289) = -3.001;
batz(499,276) = -3.001;
batz(547,203) = -3.001;
batz(551,202) = 0.0;
batz(558,195) = 0.0;
batz(568,136) = 0.0;
batz(782,1) = 0.0;
batz(422,134) = 0.0;
batz(414,102) = 0.0;
batz(384,101) = 0.0;
batz(355,93) = 0.0;
batz(342,88) = 0.0;
batz(335,76) = 0.0;
batz(334,73) = 0.0;
batz(320,46) = 0.0;
batz(323,39) = 0.0;
batz(324,37) = 0.0;
batz(324,34) = 0.0;
batz(344,22) = 0.0;
batz(334,14) = 0.0;
batz(328,7) = 0.0;
batz(543,111) = 0.0;
batz(547,119) = 0.0;

for l=2:length(batz(:,1))-1
    for m=2:length(batz(1,:))-1
        if ((batz(l,m) < 0.0) && (batz(l,m) >= -3.001) && ((batz(l+1,m) == 0.0) || (batz(l-1,m) == 0.0) ...
               || (batz(l,m-1) == 0.0) || (batz(l,m+1) == 0.0)))
           batz(l,m) = -3.001;
        end
        if ((batz(l,m) < -3.001) && (batz(l,m) >= -6.235) && ((batz(l+1,m) == 0.0) || (batz(l-1,m) == 0.0) ...
               || (batz(l,m-1) == 0.0) || (batz(l,m+1) == 0.0)))
           batz(l,m) = -6.235;
        end
    end
end


% Arno:
batz(549+i,203+j)=-3.001;
batz(549+i:553+i,203+j)=-3.001;
batz(553+i,204+j)=-3.001;
% Fiora:
batz(714+i,30+j)=-3.001;
batz(714+i,30+j:34+j)=-3.001;
% Ombrone:
batz(643+i,73+j)=-3.001;
batz(643+i:647+i,73+j)=-3.001;
batz(647+i,74+j)=-3.001;
% Serchio:
batz(547+i,217+j)=-3.001;
batz(547+i:551+i,217+j)=-3.001;
batz(551+i,218+j)=-3.001;
% Magra:
batz(512+i,250+j)=-3.001;
batz(512+i,250+j:254+j)=-3.001;
% Roia:
batz(206+i,217+j)=-3.001;
batz(206+i,217+j:221+j)=-3.001;
% Var:
batz(154+i,200+j)=-3.001;
batz(154+i,200+j:204+j)=-3.001;
% Argens:
batz(95+i,168+j)=-3.001;
batz(91+i:95+i,168+j)=-3.001;
batz(91+i,169+j)=-3.001;

batzNaN=batz;
batzNaN(batzNaN==0.0)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
date_year=datestr(datenum(date,'yyyymmdd'),'yyyy');
pathinC=strcat(glob_path2,slash,'meteo_nc',slash,date_year,slash); % path wrf NetCDF data
pathinM=strcat(glob_path2,slash,'meteo_txt',slash,date,slash); % path output

% mask wrf
%lsm=ncread([pathinM, 'const.nc'],'LAND_GDS10_SFC'); % land-sea mask
%crot=ncread([pathinM, 'const.nc'],'g10_rot_2'); crot=double(crot); % rotation (due to the different grid projection)
% per forecast delle 00
%clat=ncread([pathinC,date, '00_arw_ecm_3km.nc'],'latitude'); 
%clon=ncread([pathinC,date, '00_arw_ecm_3km.nc'],'longitude'); 
% per forecast delle 12
daypath=datestr(datenum(date,'yyyymmdd')-1,'yyyymmdd');
if exist(strcat(pathinC,daypath,'12_arw_ecm_3km.nc'), 'file') == 2
    fileexist = strcat(pathinC,daypath,'12_arw_ecm_3km.nc');
elseif exist(strcat(pathinC,daypath,'00_arw_ecm_3km.nc'), 'file') == 2
    fileexist = strcat(pathinC,daypath,'00_arw_ecm_3km.nc');
else
    daypath = datestr(datenum(date,'yyyymmdd')-2,'yyyymmdd');
    if exist(strcat(pathinC,daypath,'12_arw_ecm_3km.nc'), 'file') == 2
        fileexist = strcat(pathinC,daypath,'12_arw_ecm_3km.nc');
    else
        disp('file ', fileexist, ' not found');
        return;
    end
end
        
%clat=ncread([pathinC,daypath, '12_arw_ecm_3km.nc'],'latitude'); 
%clon=ncread([pathinC,daypath, '12_arw_ecm_3km.nc'],'longitude');
clat=ncread(fileexist,'latitude'); 
clon=ncread(fileexist,'longitude');
clat=double(clat);
clon=double(clon);
[CLAT,CLON]=meshgrid(clat,clon);%double(clat);
Cmask = batzNaN;
Cmask(batzNaN<0) = 1.0;

% variables and output files
aqh=1;lwdown=2;swdown=3;precip=4;apress=5;...
    uwind=6;vwind=7;atemp=8;lwflux=9;swflux=10;
NAMEVAROUT={'aqh','lwdown','swdown','precip','apress',...
    'uwind','vwind','atemp'};%,'lwflux','swflux'}; % runoff a parte
VAROUT=[aqh,lwdown,swdown,precip,apress,...
    uwind,vwind,atemp];%,lwflux,swflux];
% size output
% Adriatic
% mlat=39.0:0.02:46.0;
% mlon=12.0:0.02:21.0;
% NADRI
% mlat=43.0:0.02:46.0;
% mlon=12.0:0.02:16.5;
% AZA Lazio
%mlat=40.5:0.02:42.5;
%mlon=11.0:0.02:14.5;
%mlat=41.88125:1/128:44.4984;
mlat=linspace(41.88125+1/256,44.4984+1/256,336);
mlon=6.0625+1/256:1/128:12.1797+1/256;
[MLAT,MLON]=meshgrid(mlat,mlon);
%[MLON,MLAT]=meshgrid(mlon,mlat);

% time parameters (2018)
%dstart=datenum('20220106','yyyymmdd');
dend=datenum(date,'yyyymmdd')-1;
dstart = dend-7;  %%%% Cambiare se si cambia il numero di giorni di analisi

varm=zeros([size(CLON),24]); % mean variables (24 records/hours)
varc=zeros([size(CLON),25]); % cumulated variables (25 records/hours)

% cycle on variables
for v=1:size(VAROUT,2)
    VAR=VAROUT(v);
    % open output file
    eval(['fidout' int2str(v) '=fopen(''' pathinM 'BC_' NAMEVAROUT{VAR} ''',''w'',''l'');']);
    disp(NAMEVAROUT{VAR})
    
    % cycle on days
    for day=dstart:dend % dstart=20180101; dend=20181231
        disp(day)
        fileexist = strcat(pathinC,datestr(day,'yyyymmdd'),'12_arw_ecm_3km.nc');
        day_i = 0;
        if exist(fileexist, 'file') ~= 2
            if exist(strcat(pathinC,datestr(day+1,'yyyymmdd'),'00_arw_ecm_3km.nc'),'file') == 2
                dayex=day+1;
                run_start='00';
                day_i = 1;
                disp('file non presente, uso giorno +1 alle 00')
            elseif exist(strcat(pathinC,datestr(day,'yyyymmdd'),'00_arw_ecm_3km.nc'),'file') == 2
                dayex=day;
                run_start='00';
                day_i = 0;
                disp('file non presente, uso giorno +0 alle 00')
            elseif exist(strcat(pathinC,datestr(day-1,'yyyymmdd'),'12_arw_ecm_3km.nc'),'file') == 2
                dayex=day-1;
                run_start='12';
                day_i = -1;
                disp('file non presente, uso giorno -1 alle 12')
            elseif exist(strcat(pathinC,datestr(day-1,'yyyymmdd'),'00_arw_ecm_3km.nc'),'file') == 2
                dayex=day-1;
                run_start='00';
                day_i = -1;
                disp('file non presente, uso giorno -1 alle 00')
            end
        else
            dayex=day;
            run_start='12';
        end
        if strcmp(NAMEVAROUT{VAR},'aqh')
            clear var
            eval(['d2m=ncread(''' pathinC ...
                datestr(dayex,'yyyymmdd') run_start '_arw_ecm_3km.nc'',''t2m'');'])
            eval(['msl=ncread(''' pathinC ...
                datestr(dayex,'yyyymmdd') run_start '_arw_ecm_3km.nc'',''msl'');'])
            Td=d2m-273.15; % K to deg C
            p=msl./100; % Pa to mb
            
            e = 6.112.*exp((17.67.*Td)./(Td + 243.5));
            q = (0.622 .* e)./(p - (0.378 .* e));
            
            % https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
            % (Bolton 1980)
            % where:
            %    e = vapor pressure in mb;
            %    Td = dew point in deg C;
            %    p = surface pressure in mb;
            %    q = specific humidity in kg/kg.
            %
            % (Note the final specific humidity units are in g/kg = (kg/kg)*1000.0)
            
            if (day ~= dend) && (day_i == 0) && strcmp(run_start,'12')
                var=q(:,:,1+12:24+12);
            elseif (day ~= dend) && (day_i == 1) && strcmp(run_start,'00')
                var=q(:,:,1:24);
            elseif (day ~= dend) && (day_i == 0) && strcmp(run_start,'00')
                var=q(:,:,1+24:24+24);
            elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'12')
                var=q(:,:,1+36:24+36);
            elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'00')
                var=q(:,:,1+48:24+48);
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
                var=q(:,:,1+12:72+12);
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'00')
                var=q(:,:,1+24:72);
            elseif (day == dend) && (day_i == -1) && strcmp(run_start,'12')
                var=q(:,:,1+36:48+36);
            end
            
        elseif strcmp(NAMEVAROUT{VAR},'lwdown') % average value (since the beginning of the forecast)
            %  cycle on 24 records
            clear var
            %for md=1:24
                eval(['varh=ncread(''' pathinC ...
                    datestr(dayex,'yyyymmdd') run_start '_arw_ecm_3km.nc'',''p3115'');'])
                %varm(:,:,md)=varh;
                
                if (day ~= dend) && (day_i == 0) && strcmp(run_start,'12')
                    var=varh(:,:,1+12:24+12);
                elseif (day ~= dend) && (day_i == 1) && strcmp(run_start,'00')
                    var=varh(:,:,1:24);
                elseif (day ~= dend) && (day_i == 0) && strcmp(run_start,'00')
                    var=varh(:,:,1+24:24+24);
                elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'12')
                    var=varh(:,:,1+36:24+36);
                elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'00')
                    var=varh(:,:,1+48:24+48);
                elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
                    var=varh(:,:,1+12:72+12);
                elseif (day == dend) && (day_i == 0) && strcmp(run_start,'00')
                    var=varh(:,:,1+24:72);
                elseif (day == dend) && (day_i == -1) && strcmp(run_start,'12')
                    var=varh(:,:,1+36:48+36);
                end
                % de-averaging
%                 if md==1
%                     var(:,:,md)=varm(:,:,md); %#ok<SAGROW>
%                 else
%                     var(:,:,md)=varm(:,:,md)*md-varm(:,:,md-1)*(md-1); %#ok<SAGROW>
%                 end
            %end
            
        elseif strcmp(NAMEVAROUT{VAR},'swdown') % average value (since the beginning of the forecast)
            %  cycle on 24 records
            clear var
%            for md=1:24
                % sum of the direct and diffused component of sw radiation
                eval(['swdr=ncread(''' pathinC ...
                    datestr(dayex,'yyyymmdd') run_start '_arw_ecm_3km.nc'',''p3116'');'])
%                 eval(['swdf=ncread(''' pathinC 'ci2_ogs_'...
%                     datestr(day,'yyyymmdd') '.nc'',''VAR_23_GDS10_SFC_ave' int2str(md) 'h'');'])
%                varm(:,:,md)=swdr+swdf;
                
                if (day ~= dend) && (day_i == 0) && strcmp(run_start,'12')
                    var=swdr(:,:,1+12:24+12);
                elseif (day ~= dend) && (day_i == 1) && strcmp(run_start,'00')
                    var=swdr(:,:,1:24);
                elseif (day ~= dend) && (day_i == 0) && strcmp(run_start,'00')
                    var=swdr(:,:,1+24:24+24);
                elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'12')
                    var=swdr(:,:,1+36:24+36);
                elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'00')
                    var=swdr(:,:,1+48:24+48);
                elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
                    var=swdr(:,:,1+12:72+12);
                elseif (day == dend) && (day_i == 0) && strcmp(run_start,'00')
                    var=swdr(:,:,1+24:72);
                elseif (day == dend) && (day_i == -1) && strcmp(run_start,'12')
                    var=swdr(:,:,1+36:48+36);
                end
                % de-averaging
%                 if md==1
%                     var(:,:,md)=varm(:,:,md); %#ok<SAGROW>
%                 else
%                     var(:,:,md)=varm(:,:,md)*md-varm(:,:,md-1)*(md-1); %#ok<SAGROW>
%                 end
%             end
            
        elseif strcmp(NAMEVAROUT{VAR},'precip') % cumulated value (since the beginning of the forecast)
            %  cycle on 25 records
            clear var
 %           for md=1:25
 %               if md==1
                    eval(['varh=ncread(''' pathinC ...
                        datestr(dayex,'yyyymmdd') run_start '_arw_ecm_3km.nc'',''p3059'');'])
%                     varc(:,:,md)=varh;
                    
                    if (day ~= dend) && (day_i == 0) && strcmp(run_start,'12')
                        var=varh(:,:,1+12:24+12);
                    elseif (day ~= dend) && (day_i == 1) && strcmp(run_start,'00')
                        var=varh(:,:,1:24);
                    elseif (day ~= dend) && (day_i == 0) && strcmp(run_start,'00')
                        var=varh(:,:,1+24:24+24);
                    elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'12')
                        var=varh(:,:,1+36:24+36);
                    elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'00')
                        var=varh(:,:,1+48:24+48);
                    elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
                        var=varh(:,:,1+12:72+12);
                    elseif (day == dend) && (day_i == 0) && strcmp(run_start,'00')
                        var=varh(:,:,1+24:72);
                    elseif (day == dend) && (day_i == -1) && strcmp(run_start,'12')
                        var=varh(:,:,1+36:48+36);
                    end
%                 else
%                     eval(['varh=ncread(''' pathinC 'arw_ecm_3km_'...
%                         datestr(day,'yyyymmdd') '00.nc'',''tp' int2str(md-1) 'h'');'])
%                     varc(:,:,md)=varh;
%                 end
 %           end
            % de-cumulation
%             for md=1:24
%                 var(:,:,md)=varc(:,:,md+1)-varc(:,:,md);
%             end
            var=var/1000/3600; % from kg/m^2 (i.e., mm/hr) to m/s
            
        elseif strcmp(NAMEVAROUT{VAR},'apress')
            clear var
            eval(['var=ncread(''' pathinC datestr(dayex,'yyyymmdd') run_start '_arw_ecm_3km.nc'',''msl'');'])

            if (day ~= dend) && (day_i == 0) && strcmp(run_start,'12')
                var=var(:,:,1+12:24+12);
            elseif (day ~= dend) && (day_i == 1) && strcmp(run_start,'00')
                var=var(:,:,1:24);
            elseif (day ~= dend) && (day_i == 0) && strcmp(run_start,'00')
                var=var(:,:,1+24:24+24);
            elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'12')
                var=var(:,:,1+36:24+36);
            elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'00')
                var=var(:,:,1+48:24+48);
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
                var=var(:,:,1+12:72+12);
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'00')
                var=var(:,:,1+24:72);
            elseif (day == dend) && (day_i == -1) && strcmp(run_start,'12')
                var=var(:,:,1+36:48+36);
            end
            
        elseif strcmp(NAMEVAROUT{VAR},'uwind')
            clear var
            eval(['uw=ncread(''' pathinC datestr(dayex,'yyyymmdd') run_start '_arw_ecm_3km.nc'',''u10'');'])
%            uw=uw(:,:,1:24);
            
            if (day ~= dend) && (day_i == 0) && strcmp(run_start,'12')
                uw=uw(:,:,1+12:24+12);
            elseif (day ~= dend) && (day_i == 1) && strcmp(run_start,'00')
                uw=uw(:,:,1:24);
            elseif (day ~= dend) && (day_i == 0) && strcmp(run_start,'00')
                uw=uw(:,:,1+24:24+24);
            elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'12')
                uw=uw(:,:,1+36:24+36);
            elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'00')
                uw=uw(:,:,1+48:24+48);
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
                uw=uw(:,:,1+12:72+12);
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'00')
                uw=uw(:,:,1+24:72);
            elseif (day == dend) && (day_i == -1) && strcmp(run_start,'12')
                uw=uw(:,:,1+36:48+36);
            end
%             eval(['vw=ncread(''' pathinC 'arw_ecm_3km_' datestr(day,'yyyymmdd') '.nc'',''V_GDS10_HTGL'');'])
%             vw=vw(:,:,1:24);
            if (day ~= dend)
                for hr=1:24
                    var(:,:,hr)=uw(:,:,hr);%sin(crot).*vw(:,:,hr)+cos(crot).*uw(:,:,hr);
                end
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
                for hr=1:72
                    var(:,:,hr)=uw(:,:,hr);%sin(crot).*vw(:,:,hr)+cos(crot).*uw(:,:,hr);
                end
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'00')
                for hr=1:48
                    var(:,:,hr)=uw(:,:,hr);%sin(crot).*vw(:,:,hr)+cos(crot).*uw(:,:,hr);
                end
            elseif (day == dend) && (day_i == -1) && strcmp(run_start,'12')
                for hr=1:48
                    var(:,:,hr)=uw(:,:,hr);%sin(crot).*vw(:,:,hr)+cos(crot).*uw(:,:,hr);
                end
            end
%             for hr=1:24
%                 var(:,:,hr)=uw(:,:,hr);%sin(crot).*vw(:,:,hr)+cos(crot).*uw(:,:,hr);
%             end
        elseif strcmp(NAMEVAROUT{VAR},'vwind')
            clear var
%             eval(['uw=ncread(''' pathinC 'ci2_ogs_' datestr(day,'yyyymmdd') '.nc'',''U_GDS10_HTGL'');'])
%             uw=uw(:,:,1:24);
            eval(['vw=ncread(''' pathinC datestr(dayex,'yyyymmdd') run_start '_arw_ecm_3km.nc'',''v10'');'])
%            vw=vw(:,:,1:24);
            
            if (day ~= dend) && (day_i == 0) && strcmp(run_start,'12')
                vw=vw(:,:,1+12:24+12);
            elseif (day ~= dend) && (day_i == 1) && strcmp(run_start,'00')
                vw=vw(:,:,1:24);
            elseif (day ~= dend) && (day_i == 0) && strcmp(run_start,'00')
                vw=vw(:,:,1+24:24+24);
            elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'12')
                vw=vw(:,:,1+36:24+36);
            elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'00')
                vw=vw(:,:,1+48:24+48);
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
                vw=vw(:,:,1+12:72+12);
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'00')
                vw=vw(:,:,1+24:72);
            elseif (day == dend) && (day_i == -1) && strcmp(run_start,'12')
                vw=vw(:,:,1+36:48+36);
            end
            
            if (day ~= dend)
                for hr=1:24
                    var(:,:,hr)=vw(:,:,hr);%sin(crot).*vw(:,:,hr)+cos(crot).*uw(:,:,hr);
                end
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
                for hr=1:72
                    var(:,:,hr)=vw(:,:,hr);%sin(crot).*vw(:,:,hr)+cos(crot).*uw(:,:,hr);
                end
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'00')
                for hr=1:48
                    var(:,:,hr)=vw(:,:,hr);%sin(crot).*vw(:,:,hr)+cos(crot).*uw(:,:,hr);
                end
            elseif (day == dend) && (day_i == -1) && strcmp(run_start,'12')
                for hr=1:48
                    var(:,:,hr)=vw(:,:,hr);%sin(crot).*vw(:,:,hr)+cos(crot).*uw(:,:,hr);
                end
            end
%             for hr=1:24
%                 var(:,:,hr)=vw(:,:,hr);%cos(crot).*vw(:,:,hr)-sin(crot).*uw(:,:,hr);
%             end
        elseif strcmp(NAMEVAROUT{VAR},'atemp')
            clear var
            eval(['var=ncread(''' pathinC datestr(dayex,'yyyymmdd') run_start '_arw_ecm_3km.nc'',''t2m'');'])
%            var=var(:,:,1:24);
            
            if (day ~= dend) && (day_i == 0) && strcmp(run_start,'12')
                var=var(:,:,1+12:24+12);
            elseif (day ~= dend) && (day_i == 1) && strcmp(run_start,'00')
                var=var(:,:,1:24);
            elseif (day ~= dend) && (day_i == 0) && strcmp(run_start,'00')
                var=var(:,:,1+24:24+24);
            elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'12')
                var=var(:,:,1+36:24+36);
            elseif (day ~= dend) && (day_i == -1) && strcmp(run_start,'00')
                var=var(:,:,1+48:24+48);
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
                var=var(:,:,1+12:72+12);
            elseif (day == dend) && (day_i == 0) && strcmp(run_start,'00')
                var=var(:,:,1+24:72);
            elseif (day == dend) && (day_i == -1) && strcmp(run_start,'12')
                var=var(:,:,1+36:48+36);
            end
            
%         elseif strcmp(NAMEVAROUT{VAR},'lwflux') % average value (since the beginning of the forecast)
%             %  cycle on 24 records
%             clear var
%             for md=1:24
%                 eval(['varh=ncread(''' pathinC 'ci2_ogs_'...
%                     datestr(day,'yyyymmdd') '.nc'',''ATHB_S_GDS10_SFC_ave' int2str(md) 'h'');'])
%                 varm(:,:,md)=varh;
%                 % de-averaging
%                 if md==1
%                     var(:,:,md)=varm(:,:,md); %#ok<SAGROW>
%                 else
%                     var(:,:,md)=varm(:,:,md)*md-varm(:,:,md-1)*(md-1); %#ok<SAGROW>
%                 end
%             end
%             % different sign convention
%             var=-var;
            
%         elseif strcmp(NAMEVAROUT{VAR},'swflux') % average value (since the beginning of the forecast)
%             %  cycle on 24 records
%             clear var
%             for md=1:24
%                 eval(['varh=ncread(''' pathinC 'ci2_ogs_'...
%                     datestr(day,'yyyymmdd') '.nc'',''ASOB_S_GDS10_SFC_ave' int2str(md) 'h'');'])
%                 varm(:,:,md)=varh;
%                 % de-averaging
%                 if md==1
%                     var(:,:,md)=varm(:,:,md); %#ok<SAGROW>
%                 else
%                     var(:,:,md)=varm(:,:,md)*md-varm(:,:,md-1)*(md-1); %#ok<SAGROW>
%                 end
%             end
%             % different sign convention
%             var=-var;
         end
        % %{
        if (day ~= dend)
        for md=1:24
            disp(md)
            vart=squeeze(var(:,:,md));
            % spatial interpolation
            % var2Dint=griddata(CLON(~isnan(Cmask)),CLAT(~isnan(Cmask)),...
            %     vart(~isnan(Cmask)),MLON,MLAT,'nearest'); %#ok<GRIDD>
            % reduced grid to reduce computational costs
            %mx=160;Mx=290;my=190;My=290;
            %C2LON=CLON(mx:Mx,my:My); C2LAT=CLAT(mx:Mx,my:My);
            %C2mask=Cmask(mx:Mx,my:My); 
            %vart2=vart(mx:Mx,my:My);
            %F=scatteredInterpolant(CLON(~isnan(Cmask)),CLAT(~isnan(Cmask)),...
            %    vart(~isnan(Cmask)),'nearest');
            %F=scatteredInterpolant(CLON,CLAT,vart,'nearest');
            %var2Dint=F(MLON,MLAT);
            var2Dint=interp2(CLAT,CLON,vart,MLAT,MLON);
            var2DintM = var2Dint.*Cmask;
            var2DintM(isnan(var2DintM)) = 0.0;
            %var2Dint=vart2;
            % further "filterings"
%             if strcmp(NAMEVAROUT{VAR},'precip') || strcmp(NAMEVAROUT{VAR},'swdown')
%                 var2Dint(var2Dint<0.0)=0.0;
%             elseif strcmp(NAMEVAROUT{VAR},'swflux')
%                 var2Dint(var2Dint>0.0)=0.0;
%             end
            
            eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
            % "extra" records
            if day==dend && md==24
                eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
                eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
                eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
            end
        end
        
        elseif (day == dend) && (day_i == 0) && strcmp(run_start,'12')
        
        for md=1:72
            disp(md)
            vart=squeeze(var(:,:,md));
            % spatial interpolation
            % var2Dint=griddata(CLON(~isnan(Cmask)),CLAT(~isnan(Cmask)),...
            %     vart(~isnan(Cmask)),MLON,MLAT,'nearest'); %#ok<GRIDD>
            % reduced grid to reduce computational costs
            %mx=160;Mx=290;my=190;My=290;
            %C2LON=CLON(mx:Mx,my:My); C2LAT=CLAT(mx:Mx,my:My);
            %C2mask=Cmask(mx:Mx,my:My); 
            %vart2=vart(mx:Mx,my:My);
            %F=scatteredInterpolant(CLON(~isnan(Cmask)),CLAT(~isnan(Cmask)),...
            %    vart(~isnan(Cmask)),'nearest');
            %F=scatteredInterpolant(CLON,CLAT,vart,'nearest');
            %var2Dint=F(MLON,MLAT);
            var2Dint=interp2(CLAT,CLON,vart,MLAT,MLON);
            var2DintM = var2Dint.*Cmask;
            var2DintM(isnan(var2DintM)) = 0.0;
            %var2Dint=vart2;
            % further "filterings"
%             if strcmp(NAMEVAROUT{VAR},'precip') || strcmp(NAMEVAROUT{VAR},'swdown')
%                 var2Dint(var2Dint<0.0)=0.0;
%             elseif strcmp(NAMEVAROUT{VAR},'swflux')
%                 var2Dint(var2Dint>0.0)=0.0;
%             end
            
            eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
            % "extra" records
            if day==dend && md==72
                eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
                eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
                eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
            end
        end
        
        elseif (day == dend)
        for md=1:48
            disp(md)
            vart=squeeze(var(:,:,md));
            % spatial interpolation
            % var2Dint=griddata(CLON(~isnan(Cmask)),CLAT(~isnan(Cmask)),...
            %     vart(~isnan(Cmask)),MLON,MLAT,'nearest'); %#ok<GRIDD>
            % reduced grid to reduce computational costs
            %mx=160;Mx=290;my=190;My=290;
            %C2LON=CLON(mx:Mx,my:My); C2LAT=CLAT(mx:Mx,my:My);
            %C2mask=Cmask(mx:Mx,my:My); 
            %vart2=vart(mx:Mx,my:My);
            %F=scatteredInterpolant(CLON(~isnan(Cmask)),CLAT(~isnan(Cmask)),...
            %    vart(~isnan(Cmask)),'nearest');
            %F=scatteredInterpolant(CLON,CLAT,vart,'nearest');
            %var2Dint=F(MLON,MLAT);
            var2Dint=interp2(CLAT,CLON,vart,MLAT,MLON);
            var2DintM = var2Dint.*Cmask;
            var2DintM(isnan(var2DintM)) = 0.0;
            %var2Dint=vart2;
            % further "filterings"
%             if strcmp(NAMEVAROUT{VAR},'precip') || strcmp(NAMEVAROUT{VAR},'swdown')
%                 var2Dint(var2Dint<0.0)=0.0;
%             elseif strcmp(NAMEVAROUT{VAR},'swflux')
%                 var2Dint(var2Dint>0.0)=0.0;
%             end
            
            eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
            % "extra" records
            if day==dend && md==48
                eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
                eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
                eval(['fwrite(fidout' int2str(v) ',var2DintM,''float32'');']);
            end
        end
        
        end
        % %}
    end
%    figure;pcolor(var2DintM');shading flat;colorbar
    clear var vart
end

