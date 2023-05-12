%{
cross-correlation 
Timing
For Cluser 5VPS
By M. M. Wang
%}
clear all;
clc;
PathIn = 'D:\ClusterData\';
pathOut = 'D:\figure\';
path_cis='D:\CDF_data\';
%
begin_time='2005-03-26T23:43:00';
end_time='2005-03-26T23:45:00';
%%%time window for calculating the time delay
begin_time00='2005-03-26T23:43:24';%%%the begin time00 should be [timeAhead] s  after the begin time
end_time00='2005-03-26T23:43:32';
%%%time window for calculating the time delay
%%%%%%%
timeAhead = 4;%%%time search start before the start of C1, in second
timeLagM = 8;%%%max timeLag in second
%%%%%%
beginYear=str2num(begin_time(1:4));
beginMonth=str2num(begin_time(6:7));
beginDate=str2num(begin_time(9:10));
beginHour=str2num(begin_time(12:13));
beginMinute=str2num(begin_time(15:16));
beginSecond=str2num(begin_time(18:19));
endYear=str2num(end_time(1:4));
endMonth=str2num(end_time(6:7));
endDate=str2num(end_time(9:10));
endHour=str2num(end_time(12:13));
endMinute=str2num(end_time(15:16));
endSecond=str2num(end_time(18:19));
StimeInDay=datenum([beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond]); %%%compared
EtimeInDay=datenum([endYear,endMonth,endDate,endHour,endMinute,endSecond]);
clNum = 1;
[clusterVtimeC1,posVariable_out_C1,fileVersion] = Fun_getCAA_FGM_5VPS_dataV04pathset(PathIn,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond,['sc_pos_xyz_gse__C',num2str(clNum,'%1.1d'),'_CP_FGM_5VPS'],clNum);
[clusterVtimeC1,bMagVariable_out_C1,fileVersion] = Fun_getCAA_FGM_5VPS_dataV04pathset(PathIn,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond,['B_mag__C',num2str(clNum,'%1.1d'),'_CP_FGM_5VPS'],clNum);
clNum = 2;
[clusterVtimeC2,posVariable_out_C2,fileVersion] = Fun_getCAA_FGM_5VPS_dataV04pathset(PathIn,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond,['sc_pos_xyz_gse__C',num2str(clNum,'%1.1d'),'_CP_FGM_5VPS'],clNum);
[clusterVtimeC2,bMagVariable_out_C2,fileVersion] = Fun_getCAA_FGM_5VPS_dataV04pathset(PathIn,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond,['B_mag__C',num2str(clNum,'%1.1d'),'_CP_FGM_5VPS'],clNum);
clNum = 3;
[clusterVtimeC3,posVariable_out_C3,fileVersion] = Fun_getCAA_FGM_5VPS_dataV04pathset(PathIn,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond,['sc_pos_xyz_gse__C',num2str(clNum,'%1.1d'),'_CP_FGM_5VPS'],clNum);
[clusterVtimeC3,bMagVariable_out_C3,fileVersion] = Fun_getCAA_FGM_5VPS_dataV04pathset(PathIn,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond,['B_mag__C',num2str(clNum,'%1.1d'),'_CP_FGM_5VPS'],clNum);
clNum = 4;
[clusterVtimeC4,posVariable_out_C4,fileVersion] = Fun_getCAA_FGM_5VPS_dataV04pathset(PathIn,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond,['sc_pos_xyz_gse__C',num2str(clNum,'%1.1d'),'_CP_FGM_5VPS'],clNum);
[clusterVtimeC4,bMagVariable_out_C4,fileVersion] = Fun_getCAA_FGM_5VPS_dataV04pathset(PathIn,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond,['B_mag__C',num2str(clNum,'%1.1d'),'_CP_FGM_5VPS'],clNum);
%
clNum = 1;
beginYear00=str2num(begin_time00(1:4));
beginMonth00=str2num(begin_time00(6:7));
beginDate00=str2num(begin_time00(9:10));
beginHour00=str2num(begin_time00(12:13));
beginMinute00=str2num(begin_time00(15:16));
beginSecond00=str2num(begin_time00(18:19));
endYear00=str2num(end_time00(1:4));
endMonth00=str2num(end_time00(6:7));
endDate00=str2num(end_time00(9:10));
endHour00=str2num(end_time00(12:13));
endMinute00=str2num(end_time00(15:16));
endSecond00=str2num(end_time00(18:19));
StimeInDay00=datenum([beginYear00,beginMonth00,beginDate00,beginHour00,beginMinute00,beginSecond00-3]); %%%compared
EtimeInDay00=datenum([endYear00,endMonth00,endDate00,endHour00,endMinute00,endSecond00+3]);
[clusterVtimeC100,bMagVariable_out_C100,fileVersion] = Fun_getCAA_FGM_5VPS_dataV04pathset(PathIn,beginYear00,beginMonth00,beginDate00,beginHour00,beginMinute00,beginSecond00,...
    endYear00,endMonth00,endDate00,endHour00,endMinute00,endSecond00,['B_mag__C',num2str(clNum,'%1.1d'),'_CP_FGM_5VPS'],clNum);
[clusterVtimeC100,posVariable_out_C100,fileVersion] = Fun_getCAA_FGM_5VPS_dataV04pathset(PathIn,beginYear00,beginMonth00,beginDate00,beginHour00,beginMinute00,beginSecond00,...
    endYear00,endMonth00,endDate00,endHour00,endMinute00,endSecond00,['sc_pos_xyz_gse__C',num2str(clNum,'%1.1d'),'_CP_FGM_5VPS'],clNum);
%
bStartC1 = bMagVariable_out_C100(1);
tStartC1 = clusterVtimeC100(1);
inBStartC1 = find(clusterVtimeC1 == tStartC1);
%%%calculate the time delay betweent C1 and C2
fprintf('calculate the time delay betweent C1 and C2...');
for tt = 1 : 5*timeLagM
    fprintf('tt %d\n',tt);
    inBStartC2 = (inBStartC1 - timeAhead*5)  + tt;
    inBEndC2 = inBStartC2 + length(bMagVariable_out_C100) - 1;
    bMagVariable_out_C200 = bMagVariable_out_C2(inBStartC2: inBEndC2);
    if length(bMagVariable_out_C100) == length(bMagVariable_out_C200)
    A = corrcoef(bMagVariable_out_C100,bMagVariable_out_C200);%%%correlation coefficient
    B = cov(bMagVariable_out_C100,bMagVariable_out_C200);%%%covariance
    corr21(tt) = A(1,2);
    cov21(tt) = B(1,2);
    end
end%%tt = 1 : 5*timeLagM
corr21Max = max(corr21);
inCorr21Max = find(corr21 == corr21Max);
timeLag21 = (inCorr21Max - timeAhead*5)/5;
s2 = inBStartC1 - timeAhead*5 + inCorr21Max;
e2 = s2 + length(bMagVariable_out_C100) - 1;
posDiff21 = posVariable_out_C2(s2,:) - posVariable_out_C1(inBStartC1,:);
%%%calculate the time delay betweent C1 and C2
%
%%%calculate the time delay betweent C1 and C3
fprintf('calculate the time delay betweent C1 and C3...');
for tt = 1 : 5*timeLagM
    fprintf('tt %d\n',tt);
    inBStartC3 = (inBStartC1 - timeAhead*5)  + tt;
    inBEndC3 = inBStartC3 + length(bMagVariable_out_C100) - 1;
    bMagVariable_out_C300 = bMagVariable_out_C3(inBStartC3: inBEndC3);
    if length(bMagVariable_out_C100) == length(bMagVariable_out_C300)
    A = corrcoef(bMagVariable_out_C100,bMagVariable_out_C300);
    corr31(tt) = A(1,2);
    end
end%%tt = 1 : 5*timeLagM
corr31Max = max(corr31);
inCorr31Max = find(corr31 == corr31Max);
timeLag31 = (inCorr31Max - timeAhead*5)/5;
s3 = inBStartC1 - timeAhead*5 + inCorr31Max;
e3 = s3 + length(bMagVariable_out_C100) - 1;
posDiff31 = posVariable_out_C3(s3,:) - posVariable_out_C1(inBStartC1,:);
%%%calculate the time delay betweent C1 and C3
%
%%%calculate the time delay betweent C1 and C4
fprintf('calculate the time delay betweent C1 and C4...');
for tt = 1 : 5*timeLagM
    fprintf('tt %d\n',tt);
    inBStartC4 = (inBStartC1 - timeAhead*5)  + tt;
    inBEndC4 = inBStartC4 + length(bMagVariable_out_C100) - 1;
    bMagVariable_out_C400 = bMagVariable_out_C4(inBStartC4: inBEndC4);
    if length(bMagVariable_out_C100) == length(bMagVariable_out_C400)
    A = corrcoef(bMagVariable_out_C100,bMagVariable_out_C400);
    corr41(tt) = A(1,2);
    end
end%%tt = 1 : 5*timeLagM
corr41Max = max(corr41);
inCorr41Max = find(corr41 == corr41Max);
timeLag41 = (inCorr41Max - timeAhead*5)/5;
s4 = inBStartC1 - timeAhead*5 + inCorr41Max;
e4 = s4 + length(bMagVariable_out_C100) - 1;
posDiff41 = posVariable_out_C4(s4,:) - posVariable_out_C1(inBStartC1,:);
%%%calculate the time delay betweent C1 and C4
%
[vTiming,nor]=TIMING2(timeLag21,timeLag31,timeLag41,posDiff21,posDiff31,posDiff41)

clNum=1;
[Vtime_Vx1,Variable_out_Vx1,fileVersion]=Fun_getCluster_CISdataV07(path_cis,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond-300,...
    beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,['V_HIA_x_gse__C',num2str(clNum,'%1.1d'),'_PP_CIS'],clNum);
[Vtime_Vy1,Variable_out_Vy1,fileVersion]=Fun_getCluster_CISdataV07(path_cis,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond-300,...
    beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,['V_HIA_y_gse__C',num2str(clNum,'%1.1d'),'_PP_CIS'],clNum);
[Vtime_Vz1,Variable_out_Vz1,fileVersion]=Fun_getCluster_CISdataV07(path_cis,beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond-300,...
    beginYear,beginMonth,beginDate,beginHour,beginMinute,beginSecond,['V_HIA_z_gse__C',num2str(clNum,'%1.1d'),'_PP_CIS'],clNum);

[Vtime_Vx2,Variable_out_Vx2,fileVersion]=Fun_getCluster_CISdataV07(path_cis,endYear,endMonth,endDate,endHour,endMinute,endSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond + 300,['V_HIA_x_gse__C',num2str(clNum,'%1.1d'),'_PP_CIS'],clNum);
[Vtime_Vy2,Variable_out_Vy2,fileVersion]=Fun_getCluster_CISdataV07(path_cis,endYear,endMonth,endDate,endHour,endMinute,endSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond + 300,['V_HIA_y_gse__C',num2str(clNum,'%1.1d'),'_PP_CIS'],clNum);
[Vtime_Vz2,Variable_out_Vz2,fileVersion]=Fun_getCluster_CISdataV07(path_cis,endYear,endMonth,endDate,endHour,endMinute,endSecond,...
    endYear,endMonth,endDate,endHour,endMinute,endSecond + 300,['V_HIA_z_gse__C',num2str(clNum,'%1.1d'),'_PP_CIS'],clNum);

Variable_out_Vx = [Variable_out_Vx1, Variable_out_Vx2];
Variable_out_Vy = [Variable_out_Vy1, Variable_out_Vy2];
Variable_out_Vz = [Variable_out_Vz1, Variable_out_Vz2];

Vx=nanmean(Variable_out_Vx);
Vy=nanmean(Variable_out_Vy);
Vz=nanmean(Variable_out_Vz);
V_cis=[Vx; Vy; Vz];
%
vSWFrame1 = vTiming - dot(nor,V_cis);
%
U = [];
U(1) = vTiming * nor(1,:);
U(2) = vTiming * nor(2,:);
U(3) = vTiming * nor(3,:);
mVSW = sqrt(Vx^2 + Vy^2 + Vz^2);
Vcis = [Vx Vy Vz];
vSW = [];
vSW(1) = Vx/mVSW;
vSW(2) = Vy/mVSW;
vSW(3) = Vz/mVSW;
vSWFrame2 = dot(U,vSW) .* vSW - Vcis;
mVswFrame2 = sqrt(vSWFrame2(1)^2 + vSWFrame2(2)^2 + vSWFrame2(3)^2);
%%%figure
clColor = {'k','r','g','b'};
numSub = 3;
figure('Position', [1 1 1200 4000],'Units','inches');
set(gcf, 'Units', 'pixels');
%
iSub = 1;
h(iSub) = subplot(numSub,1,iSub);
tickduration=[1/60,1/12,1/6,0.25,0.5,1,5,10,15,30,60,120,180,240,300,360,720]/60/24; %minute/60/24
CompIntl=10^10;
for L=1:length(tickduration)
    if abs((EtimeInDay-StimeInDay)/tickduration(L)-8)<=CompIntl
        CompIntl=abs((EtimeInDay-StimeInDay)/tickduration(L)-8);
        itickD=L;
    end
end
plot(clusterVtimeC1, bMagVariable_out_C1,'k', clusterVtimeC2, bMagVariable_out_C2,'r',...
    clusterVtimeC3, bMagVariable_out_C3, 'g', clusterVtimeC4, bMagVariable_out_C4, 'b');
ylabel({'Bt', '(nT)'}, 'FontSize',10);
tis = sprintf('Cluster FGM 5VPS: %4.4d-%2.2d-%2.2d by CC', beginYear, beginMonth, beginDate);
tit1 = title(tis);
set(tit1,'fontsize',14);
axis tight;
ymax_min=get(h(iSub),'YLim');
if ymax_min(1)<=0&ymax_min(2)>=0
    plot([StimeInDay,EtimeInDay],[0,0],':m');
end
xlim([StimeInDay,EtimeInDay]);
x_tick=StimeInDay:tickduration(itickD):EtimeInDay;
set(gca,'xtick',x_tick,'fontsize',7);
datetick('x',13,'keeplimits','keepticks');
%
%
iSub = 2;
h(iSub) = subplot(numSub,1,iSub);
tickduration=[1/60,1/12,1/6,0.25,0.5,1,5,10,15,30,60,120,180,240,300,360,720]/60/24; %minute/60/24
CompIntl=10^10;
for L=1:length(tickduration)
    if abs((EtimeInDay00-StimeInDay00)/tickduration(L)-8)<=CompIntl
        CompIntl=abs((EtimeInDay00-StimeInDay00)/tickduration(L)-8);
        itickD=L;
    end
end
%
plot(clusterVtimeC100,bMagVariable_out_C100,'k',clusterVtimeC2(s2:e2),bMagVariable_out_C2(s2:e2),'r');
hold on;
plot(clusterVtimeC3(s3:e3),bMagVariable_out_C3(s3:e3),'g');
plot(clusterVtimeC4(s4:e4),bMagVariable_out_C4(s4:e4),'b');
ylabel({'Bt', '(nT)'}, 'FontSize',10);
%
axis tight;
ymax_min=get(h(iSub),'YLim');
if ymax_min(1)<=0&ymax_min(2)>=0
    plot([StimeInDay00,EtimeInDay00],[0,0],':m');
end
xlim([StimeInDay00,EtimeInDay00]);
x_tick=StimeInDay00:tickduration(itickD):EtimeInDay00;
set(gca,'xtick',x_tick,'fontsize',7);
datetick('x',13,'keeplimits','keepticks');
%
%
y2=get(gca,'ylim');
x2=get(gca,'xlim');
ymin2=min(min(y2));
ymax2=max(max(y2));
xmin2=min(min(x2));
xmax2=max(max(x2));
text(EtimeInDay00,ymin2+0.8*(ymax2-ymin2),'  C1','Color','k','fontsize',14);
text(EtimeInDay00,ymin2+0.6*(ymax2-ymin2),'  C2','Color','r','fontsize',14);
text(EtimeInDay00,ymin2+0.4*(ymax2-ymin2),'  C3','Color','g','fontsize',14);
text(EtimeInDay00,ymin2+0.2*(ymax2-ymin2),'  C4','Color','b','fontsize',14);
%
%
iSub = 3;
h(iSub) = subplot(numSub,1,iSub);
box on;
set(gca,'xtick',[],'fontsize',7);
set(gca,'ytick',[],'fontsize',7);
y3=get(gca,'ylim');
x3=get(gca,'xlim');
ymin3=min(min(y3));
ymax3=max(max(y3));
xmin3=min(min(x3));
xmax3=max(max(x3));
text(xmin3+0.015*(xmax3-xmin3),ymin3+0.8*(ymax3-ymin3),['Time delay (C2 - C1): ', num2str(timeLag21, '%9.1f'),'s; correlation coefficient: ', num2str(corr21Max,'%9.4f')],'fontsize',12);
text(xmin3+0.015*(xmax3-xmin3),ymin3+0.5*(ymax3-ymin3),['Time delay (C3 - C1): ', num2str(timeLag31, '%9.1f'),'s; correlation coefficient: ', num2str(corr31Max,'%9.4f')],'fontsize',12);
text(xmin3+0.015*(xmax3-xmin3),ymin3+0.2*(ymax3-ymin3),['Time delay (C4 - C1): ', num2str(timeLag41, '%9.1f'),'s; correlation coefficient: ', num2str(corr41Max,'%9.4f')],'fontsize',12);
text(xmin3+0.55*(xmax3-xmin3),ymin3+0.85*(ymax3-ymin3),['Velocity in the satellite frame: ', num2str(vTiming, '%9.2f'),' km/s'],'fontsize',12);
text(xmin3+0.55*(xmax3-xmin3),ymin3+0.65*(ymax3-ymin3),['Normal: (', num2str(nor(1,:), '%9.2f'),',',num2str(nor(2,:), '%9.2f'),',',num2str(nor(3,:), '%9.2f'),')'],'fontsize',12);
text(xmin3+0.55*(xmax3-xmin3),ymin3+0.45*(ymax3-ymin3),['Velocity in the solar wind frame: ', num2str(vSWFrame1, '%9.2f'),' km/s'],'fontsize',12);
text(xmin3+0.55*(xmax3-xmin3),ymin3+0.25*(ymax3-ymin3),['V_{SW}: (', num2str(V_cis(1,:), '%9.2f'),',',num2str(V_cis(2,:), '%9.2f'),',',num2str(V_cis(3,:), '%9.2f'),'); '...
    ,'n_{SW}: (', num2str(vSW(1), '%9.2f'), ',', num2str(vSW(2), '%9.2f'), ',', num2str(vSW(3), '%9.2f'), ')'], 'fontsize',12);
text(xmin3+0.55*(xmax3-xmin3),ymin3+0.05*(ymax3-ymin3),['Velocity in the solar wind frame: ', num2str(mVswFrame2, '%9.2f'),' km/s'],'fontsize',12);

%
dir_name=sprintf('TimeDelayTiming_C1234CC');
ifexist=dir([pathOut,dir_name]);
if length(ifexist)==0
    mkdir(pathOut,dir_name);
end
picname=sprintf('%s%s\\CC5VPSDataCavitonVelocityInSolarWindFrame_%s%s%s_GSE_%s%s%s-%s%s%s',...
    pathOut,dir_name,begin_time(1:4),begin_time(6:7),begin_time(9:10),begin_time(12:13),...
    begin_time(15:16),begin_time(18:19),begin_time(12:13),begin_time(15:16),begin_time(18:19));%可以修改文件名

saveas(gcf, picname, 'bmp');
% print('-djpeg', '-r300', picname);
% print('-dbmp', '-r500', picname);
