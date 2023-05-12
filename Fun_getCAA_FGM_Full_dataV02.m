%{
By W. J. Sun
For reading Cluster Full Resolution FGM and Pos
%}
function [Vtime,Variable_out,fileVersion]=Fun_getCAA_FGM_Full_dataV02(path,Syear,Smonth,Sdate,Shour,Sminute,Ssecond,Eyear,Emonth,Edate,Ehour,Eminute,Esecond,NameVariable_out,num_cl)%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% path='D:\MyWork\CAA_data\FGM_FULL\';
%
file_name=Fun_FindFilename_V01(path,Syear,Smonth,Sdate,Shour,Sminute,Ssecond,Eyear,Emonth,Edate,Ehour,Eminute,Esecond,NameVariable_out);
%
fileVersion=[];
%
if isempty(file_name)
    fprintf('The data file %s does not exist in the given directory. Please download the file!!!\n',file_name);
    Vtime(1)=NaN;
    Variable_out(1,1:3)=NaN;
else
    file_read=[path,'C',num2str(num_cl,'%1.1d'),'\C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL','\',file_name];
    %
    fprintf('    Reading Variable <%s> ................\n',NameVariable_out);
    Vtime(1)=NaN;
    Variable_out(1,1:3)=NaN;
    num_out=1;
    fileVersion=[];
    Stime=datenum([Syear,Smonth,Sdate,Shour,Sminute,Ssecond]);
    Etime=datenum([Eyear,Emonth,Edate,Ehour,Eminute,Esecond]);
    Iyear=Syear;
    Imonth=Smonth;
    for Idate=Sdate:1:Edate
        Epoch_cell=cdfread(file_read,'variables',['time_tags__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL'],'CombineRecords', true, 'ConvertEpochtoDatenum', true); %%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!time_tags__C1_CP_FGM_FULL
        %                 if flag_Vname~=0
        flag_exist=0;
        Len_Epoch=length(Epoch_cell);
        if Len_Epoch>1
            num_out=1;
            for k=1:Len_Epoch
                timek=Epoch_cell(k);
                if timek>=Stime && timek<=Etime
                    Vtime(num_out)=timek;
                    record_k(num_out)=k;
                    num_out=num_out+1;
                    flag_exist=1; %data existing in this duration
                end
            end
            k_start=record_k(1);
            k_end=record_k(end);
            clear  record_k;
            %         end %if flag_Vname~=0
            if flag_exist==1;
                fprintf('Successfully reading Epoch from file <%s>!\n',file_read);
                %                 flag_Vname=1;
                if isequal(NameVariable_out,['B_vec_xyz_gse__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL'])%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!B_vec_xyz_gse__C1_CP_FGM_FULL
                    Variable_cell=cdfread(file_read,'variables',['B_vec_xyz_gse__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL']);%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    num_out=1;
                    for m=k_start:k_end
                        Variable_out(num_out,:)=double(Variable_cell{m})';
                        num_out=num_out+1;
                    end
                elseif isequal(NameVariable_out,['B_mag__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL'])%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!B_mag__C1_CP_FGM_FULL
                    Variable_cell=cdfread(file_read,'variables',['B_mag__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL']);%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    num_out=1;
                    for m=k_start:k_end
                        Variable_out(num_out)=double(Variable_cell{m});
                        num_out=num_out+1;
                    end
                    %                     num_dim=1;%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                elseif isequal(NameVariable_out,['B_vec_xyz_gsm__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL'])
                    %
                    Variable_cell=cdfread(file_read,'variables',['B_vec_xyz_gse__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL']);%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    num_out=1;
                    for m=k_start:k_end
                        %
                        timek=Epoch_cell(m);
                        time_1=datevec(timek);
                        A_GSE2GSM=GSE2GSM(time_1(1),time_1(2),time_1(3),time_1(4),time_1(5),time_1(6));
                        Variable_double(m,:)=double(Variable_cell{m})';
                        B_GSM=A_GSE2GSM*[Variable_double(m,1);Variable_double(m,2);Variable_double(m,3)];
                        Variable_out(num_out,1)=B_GSM(1);
                        Variable_out(num_out,2)=B_GSM(2);
                        Variable_out(num_out,3)=B_GSM(3);
                        num_out=num_out+1;
                    end
                elseif isequal(NameVariable_out,['sc_pos_xyz_gse__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL'])%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!B_vec_xyz_gse__C1_CP_FGM_FULL
                    Variable_cell=cdfread(file_read,'variables',['sc_pos_xyz_gse__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL']);%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    num_out=1;
                    for m=k_start:k_end
                        Variable_out(num_out,:)=double(Variable_cell{m})';
                        num_out=num_out+1;
                    end
                elseif isequal(NameVariable_out,['sc_pos_xyz_gsm__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL'])
                    %
                    Variable_cell=cdfread(file_read,'variables',['sc_pos_xyz_gse__C',num2str(num_cl,'%1.1d'),'_CP_FGM_FULL']);%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    num_out=1;
                    for m=k_start:k_end
                        %
                        timek=Epoch_cell(m);
                        time_1=datevec(timek);
                        A_GSE2GSM=GSE2GSM(time_1(1),time_1(2),time_1(3),time_1(4),time_1(5),time_1(6));
                        Variable_double(m,:)=double(Variable_cell{m})';
                        Pos_GSM=A_GSE2GSM*[Variable_double(m,1);Variable_double(m,2);Variable_double(m,3)];
                        Variable_out(num_out,1)=Pos_GSM(1);
                        Variable_out(num_out,2)=Pos_GSM(2);
                        Variable_out(num_out,3)=Pos_GSM(3);
                        num_out=num_out+1;
                    end
                else
                    fprintf('!!The Variable name <%s> is not supported by this (get) function for the given file <%s> !!\n',NameVariable_out,file_read);
                    %                     flag_Vname=0;
                    break;% can not be delete if there are two days in the total duration
                end
                clear Variable_cell;
                ttt=find(Variable_out>=1E20);
                Variable_out(ttt)=NaN;
            else
                fprintf('!!No data in the required duration from file <%s>!\n',file_read);
            end
        end %        if Len_Epoch>1
        %
        
    end %if length(ifexist)==0
    % end %for Idate=Sdate:1:Edate
    % fprintf('\n');
end


%*********************************************************************************************
%******  The calculation of the transformation matrix from GSE to GSM coordinates   **********
%******  by Q. Q. SHI    2003-09-20                                                 **********
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put this file into the same dictionary in which you use some program to call it.  **********
% We should call this program like this,for exmample:   **************************************
%exmample line 1:    BGSE=[40; -10; 20];
%exmample line 2:    A=GSE2GSM(2001,10,1,5,40,30);
%exmample line 3:    BGSM=A*BGSE;
% Then we get the vector in GSM coordinate system!      **************************************
%*********************************************************************************************
function A_GSE2GSM=GSE2GSM(year,month,date,hour,minute,second)
     %Calculating the day of year & the sencond of day
      month_days=[31,28,31,30,31,30,31,31,30,31,30,31];
      if mod(year,4)==0
         month_days(2)=29;
      end
      mdsum(1)=0;
      for k=2:12
          mdsum(k)=mdsum(k-1)+month_days(k-1);
      end
      nday=date+mdsum(month);
      nsecond=hour*3600+minute*60+second;
     %Calculating the Position of the Sun in GEI Coordinates
      [S(1),S(2),S(3),GST]=SUN_DIR(year,nday,nsecond);
     %Calculating the transformation matrix from GSE to GEI =======================================
      y=cross([0,-0.3978,0.9175],S);
      A_GEI2GSE=[S(1),S(2),S(3);   y(1),y(2),y(3);   0,-0.3978,0.9175];
      A_GSE2GEI=A_GEI2GSE';
     %Calculating the transformation matrix from GEI to GSM =======================================
      A_GEO2GEI=[cos(GST),-sin(GST),0;   sin(GST),cos(GST),0;   0,0,1];
      pp1=3.1415926/180.; 
      %D=[0.06252; -0.18060; 0.98157];  %%%%%%%%%%%%% IGRF 1985.0 
      %phi=79.30*pp1;  lambda=288.59*pp1;  %%%%%%% IGRF 1995 
      phi=79.54*pp1;  lambda=288.43*pp1;   %%%%%%% IGRF 2000
      DD1=cos(phi)*cos(lambda);  DD2=cos(phi)*sin(lambda);DD3=sin(phi);  %%%%%%% IGRF 
      D=[DD1; DD2; DD3];  %%%%%%% IGRF 
      D1=A_GEO2GEI*D;
      Y0=cross(D1,S);
      aY0=sqrt(Y0(1)^2+Y0(2)^2+Y0(3)^2);
      Y=Y0./aY0;
      Z=cross(S,Y);
      A_GEI2GSM=[S(1),S(2),S(3);   Y(1),Y(2),Y(3);   Z(1),Z(2),Z(3);];
      %Calculating the transformation matrix from GSE to GSM ====OK!!
      A_GSE2GSM=A_GEI2GSM*A_GSE2GEI;
return

%***********************************************************************************
%******  The Calculation of the Position of the Sun in GEI Coordinates   ***********
%******  By G.D. Mead  (modified by sqq in 2003)                         ***********
%***********************************************************************************
function [S1,S2,S3,GST]=SUN_DIR(IYR, IDAY, SECS)
% PROGRAM TO CALCULATE SIDEREAL TIME AND POSITION OF THE SUN. 
% GOOD FOR YEARS 1901 THROUGH 2099. ACCURACY 0.006 DEGREE.
% INPUT IS IYR, IDAY (INTEGERS), AND SECS, DEFINING UN. TIME. 
% OUTPUT IS GREENWICH MEAN SIDEREAL TIME (GST) IN DEGREES,
% LONGITUDE ALONG ECLIPTIC (SLONG), AND APPARENT RIGHT ASCENSION
% AND DECLINATION (SRASN, SDEC) OF THE SUN, ALL IN DEGREES 
format long;
  RAD=57.29578;  %DATA RAD /57.29578/ 
  %DOUBLE PRECISION DJ, FDAY 
  %IF(IYR. LT. 1901. OR. IYR. GT. 2099) RETURN
  FDAY = SECS/86400;
  DJ = 365* (IYR-1900) + (IYR-1901)/4 + IDAY + FDAY -0.5; 
  T = DJ / 36525 ;
  VL = mod(279.696678 + 0.9856473354*DJ, 360.) ;
  GST = mod(279.690983 + 0.9856473354*DJ + 360.*FDAY + 180., 360.);
  G = mod(358.475845 + 0.985600267*DJ, 360.) / RAD ;
  SLONG = VL + (1.91946 -0.004789*T)*sin(G) + 0.020094*sin (2.*G); 
  OBLIQ = (23.45229 -0.0130125*T) / RAD ;
  SLP = (SLONG -0.005686) / RAD ;
  SIND = sin (OBLIQ)*sin(SLP) ;
  COSD = sqrt(1.-SIND^2);
  SDEC = RAD * atan(SIND/COSD); 
  SRASN = 180. -RAD*atan2(cot(OBLIQ)*SIND/COSD, -cos(SLP)/COSD); 
  %
  pp=3.1415926/180.;
  S1 = cos(SRASN*pp)*cos(SDEC*pp);
  S2 = sin(SRASN*pp)*cos(SDEC*pp);
  S3 = sin(SDEC*pp);
  GST = GST*pp;
return