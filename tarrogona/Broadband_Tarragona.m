function [Raw_data,Results, filename]=Broadband_Tarragona(Input_ref,p1,p2,Input_insult,p3,p4, lambda, bAutoPath, res, bSave, bVis)
%
% This programme reads the tarragona data file and output pathlength840 and
% pathlength740 together with concentrations using user defined wavelength
% range (resolution of 1nm).
%
% Created by Tingting Zhu 18/08/2011
%
% Inputs:
% Input_ref is the Tarragona's reference intensity file including its path
% eg: 'C:\....\Input_ref.tdf'
% Input_insult is the Tarragona's insult intensity file including its path
% eg: 'C:\....\Input_insult.tdf'
% p1 is the start time stamp point of the reference intensity file
% (default: p1=1)
% p2 is the end time stemp point of the reference intensity file
% p3 is the start time stamp point of the insult intensity file
% (default: p3=1)
% p4 is the end time stamp point of the insult intensity file
%
% ----------------------------------------------------------------------------------
% Please note that p2 and p4 can be obtained through running the function
% [result,stindex,endindex]=tarragona_intensityread(filename,average,stindex,endindex)
% where p2 can be obtained as
%     [result]=tarragona_intensityread(Input_ref,0,1,2);
%          p2 = result.num_spectra ;
% and p4 can be obtained as
%     [result]=tarragona_intensityread(Input_insult,0,1,2);
%          p4 = result.num_spectra ;
% ----------------------------------------------------------------------------------
%
% Outputs:
% Raw_data contains 3 structures, each structure has fields:
% - Intensity contains 2 fields:
%   -  mean_Ref indicates the mean refence intensity
%   -  Insult indicates the insult intensity
% - Abs_Attenuation indicates the absolute attenuation
% - Change_Attenuation indicates the change attenuation
% Results contains 9 structures, each structure has fields:
% - Change_attenuation indicates the fitted  and smoothed change attenuation
%   according to the user defined wavelength range
% - Abs_attenuation indicates the fitted absolute attenuation according to
%   the user defined wavelength range
% - pathlength840 indicates the pathelength estimated at 840nm by fitting
%   the 2nd differential of water (Matcher et al. 1994)
% - pathlength740 indicates the pathelength estimated at 740nm by fitting
%   the 2nd differential of water (Matcher et al. 1994)
% - Abs_Hb760 indicates the absolute concentration of Hb estimated
%   using pathlength obtained at 760nm (Matcher et al. 1994)
% - Time indicates the time stamp for the insult (p3:1:p4)
% - Concentration indicates the concentration changes, each line corresponds
%   to Hb,HbO2,CtOx,HbT,HbDiff respectively
% - Residuals contains  3 fields:
%   -  ThreeFit indicates the residuals between estimated and measured
%      change attenuation during a 3-component fit (Matcher et al. 1995)
%   -  TwoFit indicates the residuals between estimated and measure change
%      attenuation during a 2-component fit (Matcher et al. 1995)
%   -  Conc2 indicates the concentration changes of Hb, HbO2,HbT
% - Selected_Residual contains 3 fiels:
%   -  ThreeFit indicates the residuals of change attenuation in selected
%      time stamps during a 3-component fit (Matcher et al. 1995)
%   -  TwoFit indicates the residuals of change attenuation in selected
%      time stamps during a 2-component fit (Matcher et al. 1995)



set(0, 'DefaultTextInterpreter', 'tex');
global pth;
global bSave;
global visString;
[pth, name, ext] = fileparts(Input_insult);
filename = [pth, filesep, name];

if bVis
    visString = 'on';
else
    visString = 'off';
end
%--------------------------------------------------------------------------
% Read reference data
% Input_ref= input('\nEnter the reference file name:','s');
[result_ref]=tarragona_intensityread(Input_ref,0,p1,p2);
% convert from log(intensity) to intensity
result_ref.Block.intensity =10.^(result_ref.Block.intensity);

% figure;plot(result_ref.wavelength,result_ref.Block.intensity);grid;
% xlabel('Wavelength(nm)');ylabel('Intensity(counts)');
% title(['Reference Intensity for ',int2str(result_ref.num_spectra),' samples'],'fontsize',12,'fontweight','b');

% %Average all the samples of the reference data
I0=mean(result_ref.Block.intensity,ndims(result_ref.Block.intensity));
% figure;plot(result_ref.wavelength,I0);grid;
% xlabel('Wavelength(nm)');ylabel('Intensity(counts)');
% title('Averaged Reference Intensity','fontsize',12,'fontweight','b');

%--------------------------------------------------------------------------
% Read the insult data
% Input_insult= input('\nEnter the insult file name:','s');

[result_insult,stindex,endindex]=tarragona_intensityread(Input_insult,0,p3,p4);
result_insult.Block.intensity1= 10.^(result_insult.Block.intensity); % convert from log(intensity) to intensity
% figure;plot(result_insult.wavelength,result_insult.Block.intensity);grid;
% xlabel('Wavelength(nm)');ylabel('Intensity(counts)');
% title(['Measured Intensity for ',int2str(result_insult.num_spectra),' samples'],'fontsize',12,'fontweight','b');

%--------------------------------------------------------------------------

i=1:60;j=1;

if endindex >=5000
    %    I0 = I0(1:60:end);
    result_insult.Block.intensity=zeros(1022,fix(size(result_insult.Block.intensity1,2))./60);
    while max(i)<=size(result_insult.Block.intensity1,2)
        result_insult.Block.intensity(:,j) = mean(result_insult.Block.intensity1(:,i),2);
        i=i+60;
        j=j+1;
    end
    endindex=size(result_insult.Block.intensity,2);
else result_insult.Block.intensity = result_insult.Block.intensity1;
end

Raw_data.Intensity.mean_Ref=I0;
Raw_data.Intensity.Insult=result_insult.Block.intensity;
%--------------------------------------------------------------------------
% Calculate Absolute Attenuation
Abs_A=-log10(result_insult.Block.intensity)+log10(repmat(I0,1,size(result_insult.Block.intensity,2)));
% figure;plot(result_insult.wavelength,Abs_A);grid;
% xlabel('Wavelength(nm)');ylabel('Absolute Attenuation');
% title(['Absolute Attenuation for ',int2str(result_insult.num_spectra),' samples'],'fontsize',12,'fontweight','b');
Raw_data.Abs_Attenuation= Abs_A;
%--------------------------------------------------------------------------
% Calculate Change attenuation
A = result_insult.Block.intensity;
OD = log10(repmat(A(:,1),1,size(A,2))./A);
% figure;
% subplot(2,1,1);plot(result_insult.wavelength,OD);grid;
% xlabel('Wavelength(nm)');ylabel('Attenuation(OD)');
% title(['Change Attenuation for ',int2str(result_insult.num_spectra),' samples'],'fontsize',12,'fontweight','b');
Raw_data.Change_Attenuation= OD;

%--------------------------------------------------------------------------
% Fit Tarragona's wavelengths
Wavelengths= result_insult.wavelength;
%--------------------------------------------------------------------
% Fit changes in attenuation with respect to the selected wavelengths
if exist('lambda', 'var')
    lambda1 = lambda(1);
    lambda2 = lambda(2);
else
    
    I_c = menu('Please select the wavelength range from 700 to 900nm:',...
        '780 to 900nm  (default)','your own wavelength range');
    switch I_c
        case I_c==1
            lambda1= 780;
            lambda2= 900;
            flag=0;
        otherwise
            lambda1= input('Enter the minimum lambda: ');
            lambda2= input('Enter the maximum lambda: ');
            if lambda1 <740
                warning('The specific extinction coefficient will not be corrected by wavelength dependent factor!'); %#ok<WNTAG>
                flag = 1;
            else
                flag=0;
            end
    end
end

if lambda1 < 740
    flag = 1;
else
    flag = 0;
end
lambda_w = lambda1:1:lambda2;
% I_w='n';
% % I_w=input('Would you like to input the specific wavelength combinations?(y/n)','s');
% if I_w=='n'
% %       lambda1= input('Enter the minimum lambda: ');
% %       lambda2= input('Enter the maximum lambda: ');
% %       res = input('Enter resolution (nm):(default 1nm) ');
%       lambda_w = lambda1:res:lambda2;
% elseif I_w=='y'
%       lambda_w = input('Enter all the wavelengths: ');
%       lambda1= min(lambda_w);
%       lambda2= max(lambda_w);
% else
%       error('Answer eith y or n only');
% end
n2=size(lambda_w,2);

%--------------------------------------------------------------------------
% Interpolate Change Attenuation to correspond the correct wavelengths
q2=1; OD_Intp = zeros(length(lambda_w),length(stindex:endindex));
for j=stindex:endindex;
    OD_Intp(:,q2)=interp1(Wavelengths,OD(:,j),lambda_w(1:n2),'spline');
    q2=q2+1;
end
OD_fit = squeeze(OD_Intp); % OD for corresponded wavelengths
% Smooth Change Attenuation
for i =stindex:1:endindex
    OD_fit(:,i) = smooth(lambda_w,OD_fit(:,i),0.2,'rloess');
end

Results.Change_attenuation = OD_fit;
%--------------------------------------------------------------------------
% Interpolate Absolute attenuation with wavelegnths range from 700 to 900 nm
q2_A=1;Abs_A_Intp = zeros(201,length(stindex:endindex));
for j=stindex:endindex;
    Abs_A_Intp(:,q2_A)=interp1(Wavelengths,Abs_A(:,j),700:900,'spline');
    q2_A=q2_A+1;
end
% Absolute attenuation for corresponded wavelengths
Abs_A_fit = squeeze(Abs_A_Intp);
Results.Abs_attenuation = Abs_A_fit;
%------------------------------------------------------------------------------
% Calculation of pathlength & cocnentration of hhb
if exist('bAutoPath', 'var')
    choice = 1;
else
    choice = menu('Please select the following options for pathlength estimation:',...
        '740nm and 840nm (recommended)','Manually input pathlength yourself');
end
if choice==2
    pl_m= input('Enter the pathlength(cm): ');
    pl=repmat(pl_m,1,p4-p3+1);
    [~,Abs_Hb760]=pathlength(Abs_A_fit,840);close(gcf);
    Results.pathlength =pl;
elseif choice==1
    figure;[pl_840,Abs_Hb760]=pathlength(Abs_A_fit,840);
    figure;[pl_740,~]=pathlength(Abs_A_fit,740);
    Results.pathlength840=pl_840;
    Results.pathlength740=pl_740;
    pl = pl_840;
else error('You need the pathlength to estimate concentration changes!');
end
Results.Abs_Hb760 =Abs_Hb760;
Results.Time = p3:p4;
close(gcf);
%-------------------------------------------------------------------------
WavelengthsFit = lambda_w';
% subplot(2,1,2);plot(WavelengthsFit,OD_fit);grid;
% xlabel('wavelengths(nm)');ylabel('Change Attenuation')
% figure;plot(WavelengthsFit,Abs_A_fit);grid; xlabel('wavelengths(nm)');ylabel('Absolute Attenuation')

Nwl=length(WavelengthsFit);
OD_fit_results = OD_fit';
% lambda = WavelengthsFit';
K=Nwl; % select all wavelengths
[Results1]=Conc_n_residual(stindex,endindex,OD_fit_results,lambda_w,Abs_A_fit,K,lambda1,lambda2,pl,flag);
Results.Concentrations = cat(1,Results1.Conc3,Results1.HBD);
Results.Residuals.ThreeFit = Results1.err3;
Results.Residuals.TwoFit = Results1.err2;
Results.Residuals.Conc2 = Results1.Conc2;
Results.Selected_Residual = Results1.Selected_Residual;

% Save to text file
m1 = strfind(Input_insult, '\');
m2 = strfind(Input_insult, '.tdf');


output2txt(Results,filename)
end


function [PL,Abs_Hb]=pathlength(Abs_A_fit,w_selected,width,order,width_s)
% Pathlength  and absolute concentration of deoxygenated haemoglobin are
% obtained at 840nm/740nm by fitting the second differential of the absolute
% attenuation spectra to the second differential of the water and Hb
% absorption spectra between 700 to 900 nm
% (assuming an average cerebral water content of 85%)
% For default: [PL]=pathlength(Abs_A_fit,w_selected)
% Inputs:
% Abs_A_fit is the absolute attenuation that is fitted for 700 to 900nm with resolution of 1nm
% w_selected is the wavelength selected to estimate the pathlength (either
% 740nm or 840nm)
% width is the window length for the Savitzky-Golay filter
% order is the polynomial order for the Savitzky-Golay filter
% width_s is the window size for a moving average filter for the purpose of smoothing the signal
% Outputs:
% PL is the estimated pathlength using 840nm second differential of water absorption.
% Abs_Hb is the absolute concentration of Hb estimated using pathlength obtained at 760nm/840nm

% Note: If the 740nm is chosen for estimating the pathlength, the width_s
% should be as small as possible (eg: width_s = 5) since there are
% different chromophores (eg: water,Hb,HbO)that contribute to the
% second differential of the absolute attenuation
%
% Created by Tingting Zhu 27/06/11


if nargin==2
    width=11; % Matcher et al. 1994
    order=3;  % Matcher et al. 1994
    width_s=21;
end;
global pth;
global bSave;
global visString;

% Specific	extinction	coeffs	measured by	Mark Cope
c1=abs_Hb; % Measured in OD/mm/uM
indx11= find(c1(:,1)==700);
indx21= find(c1(:,1)==900);
c_hb = c1(indx11:indx21,2)*10000; % Measured in OD/cm/mMolar
c_hb =smooth(c_hb,width_s);  % Smoothing

% Specific	extinction	coeffs of water	measured by	Matcher
c=abs_water;
indx1= find(c(:,1)==700);
indx2= find(c(:,1)==900);
c_water = c(indx1:indx2,2);% measured in OD/cm


Att = Abs_A_fit';
od_range= Att(:,1:end);
[i,j]=size(od_range);
N=order; % the polynomial order
F=width;  % window length
[~,g]=sgolay(N,F);% Calculate S-G coefficients
% dx=1;
% figure;
c_water =smooth(c_water,width_s);  % Smoothing
HalfWin  = ((F+1)/2) -1;
w_range=701:895;
k=1;
while k <=i
    
    od = od_range(k,:);
    od_pl=smooth(od,width_s); % Smoothing for 840nm 2nd differential of water
    od_Hb=smooth(od,21); % Smoothing for 760nm 2nd differential of Hb
    
    for n = (F+1)/2:j-(F+1)/2
        % Second differential
        zw2(n)=2*dot(g(:,3)',c_water(n - HalfWin: n + HalfWin)); %#ok<*AGROW> %water
        zh2(n)=2*dot(g(:,3)',c_hb(n - HalfWin: n + HalfWin)); %Hb
        d2_pl(n)=2*dot(g(:,3)',od_pl(n - HalfWin: n + HalfWin)');% Abs. OD
        d2_Hb(n)=2*dot(g(:,3)',od_Hb(n - HalfWin: n + HalfWin)');% Abs. OD
    end
    %    % Turn differential into second derivative
    %    zw2= zw2/(dx*dx);
    %    zh2= zh2/(dx*dx);
    %    d2= d2/(dx*dx);
   
    figure('Visible', visString);
    plot(w_range,zw2,'r',w_range,d2_pl,'k',w_range,zh2,'b');grid;
    legend('water(ODcm^-^1nm^-^2)','Absolute Attenuation(ODnm^-^2)','HHb(ODcm^-^1mM^-^1nm^-^2)');
    xlabel('Wavelength(nm)');
    title('Second differential');
    
    % Interpolates the 2nd differential for selected wavelength
    %    w_selected=840;
    zw2_8(:,k)=interp1(w_range,zw2,w_selected);
    d2_8(:,k)=interp1(w_range,d2_pl,w_selected);
    % Estimate the pathlength @840nm/740nm
    PL(:,k)=d2_8(:,k)*pinv((zw2_8(:,k))')'/0.85;
    
    w_selected2=760;
    zw2_Hb(:,k)=interp1(w_range,zw2,w_selected2);
    d2_Hb_760(:,k)=interp1(w_range,d2_Hb,w_selected2);
    % Estimate the pathlength @760nm
    PL_Hb(:,k)=d2_Hb_760(:,k)*pinv((zw2_Hb(:,k))')'/0.85;
    
    %estimate the required amount of pure 2nd differential chromophore
    abs2= [zw2;zh2];
    conc1(:,k)=(d2_Hb*pinv(abs2)*PL_Hb(:,k))';% measured in mMolar %760nm PL
    conc(:,k)=(d2_Hb*pinv(abs2)*PL(:,k))';%#ok<NASGU> % measured in mMolar   %840nm PL
    k=k+1;
end
if bSave
    figSaveEPS(gcf,[pth, filesep, 'abs_atten_', int2str(w_selected),'.eps'],'fixeps', true);
end
close(gcf);

% Estimate the absolute Hb
H2O = conc1(1,:);
Hb = conc1(2,:);
Abs_Hb=Hb*0.85./H2O*1000;  % measured in uMolar
% mean_Hb= mean(Abs_Hb);
% sd_Hb=std(Abs_Hb);
% figure;plot(Abs_Hb);xlabel('Time(min)');ylabel('Absolute Hb(uM)');grid;



%Absolute HG Figure
figure('Visible', visString);
plot(0:i-1, Abs_Hb);xlabel('Time(mins)');
ylabel('Absolute Hb (uMolar)');grid;
title(['Absolute Hb estimated using ',int2str(w_selected2),'nm']);
if bSave
    figSaveEPS(gcf,[pth, filesep, 'absolute_Hb_', int2str(w_selected),'.eps'],'fixeps', true);
end
close(gcf);

%Pathlength estimate
figure('Visible', visString);
plot(0:i-1,PL);xlabel('Time(mins)');
ylabel('Pathlength(cm)');
grid;
title(['Pathlength estimated using ',int2str(w_selected),'nm']);
%----------------------------------------------------------------------------------
%smooth pathlength
PL=smooth(1:i,PL,0.05,'rloess')';
hold on;
plot(0:i-1,PL,'r','LineWidth',2);legend('Before smoothing','After smoothing');
if bSave
    figSaveEPS(gcf,[pth, filesep, 'absolute_pathlength_', int2str(w_selected),'.eps'],'fixeps', true);
end
close(gcf);
end


function [Results]=Conc_n_residual(stindex,endindex,OD_fit_results,lambda,Abs_A_fit,K,lambda1,lambda2,pl,flag)
global pth;
global bSave;
global visString;
% calculation of Concentration & Residuals
t=stindex:endindex;
R_fit2=zeros(1,size(OD_fit_results,3));
R_fit3=R_fit2; %#ok<NASGU>

[Conc_2,err_OD_2,Hb_T_2,~,R_fit2]=OD2Conc(K,t,2,lambda,OD_fit_results',lambda1,lambda2,stindex,endindex,Abs_A_fit,pl,flag);
[Conc_3,err_OD_3,Hb_T_3,~,R_fit3]=OD2Conc(K,t,3,lambda,OD_fit_results',lambda1,lambda2,stindex,endindex,Abs_A_fit,pl,flag);
Conc_2=cat(1,Conc_2,Hb_T_2); % hhb;hbo;ctOx;hbT

Results.Conc2=Conc_2;
% Results.Conc3=Conc_3;
Results.err2=err_OD_2;
Results.err3=err_OD_3;
Results.R_fit2=R_fit2;
Results.R_fit3=R_fit3;
% Results.OD=OD_fit_results;
% Results.pl = pl;
%---------------------------------------------------------------------------------------------------
%Estimate HBO-HHB
if size(OD_fit_results,3)==1
    HBD=Conc_3(2,:)-Conc_3(1,:);
    figure('Visible', visString);
    plot(t,HBD);
    grid;%title(nae);
    Results.HBD= HBD;
    xlabel('Time(min)');
    ylabel('Concentration(uMolar)');
    legend('HbDiff');
    if bSave
        figSaveEPS(gcf,[pth, filesep, 'HbDiff.eps'],'fixeps', true);
    end
    close(gcf);
    
end

Results.Conc3 = cat(1,Conc_3,Hb_T_3);

%---------------------------------------------------------------------------------------------------
fig2=figure('Visible', visString);
Conc_3_r=Results.Conc3(:,t);
plot(t,Conc_3_r(1,:),'b',t,Conc_3_r(2,:),'r',t,Conc_3_r(3,:),'g',t,Conc_3_r(4,:),'k');
grid;
legend('Hb','HbO_2','CtOx','HbT');
xlabel('Time(min)');
ylabel('Concentration(uMolar)');%title(Input_insult)
if bSave
    figSaveEPS(gcf,[pth, filesep, 'concentrations.eps'],'fixeps', true);
end
close(gcf);

[error_sel_3fit,S]= sel_error(err_OD_3(:,t),fig2); % errors of 1st combination of wavelengths
if error_sel_3fit==0
    Selected_Residual.ThreeFit.ErrorLambda= error_sel_3fit;
    disp('No points are selected for calculating residual!');
else % No wavelength combination is selected or all wavelengths are selected
    error_lambda3 = cat(2,lambda',error_sel_3fit(2:end,:));
    Selected_Residual.ThreeFit.ErrorLambda= error_lambda3;
end
Selected_Residual.ThreeFit.Time= error_sel_3fit(1,:);

%--------------------------------------------------------------------------------------------------
%Selection for residual calculation

Conc_2_r=Conc_2(:,t);
fig1=figure('Visible', visString);%title(['Change Concentration for wavelength ',int2str(lambda1),' to ',int2str(lambda2), ' nm'],'fontsize',12,'fontweight','b');
plot(t,Conc_2_r(1,:),'b',t,Conc_2_r(2,:),'r',t,Conc_2_r(3,:),'k');
grid;
legend('Hb','HbO_2','HbT');
xlabel('Time(min)');
ylabel('Concentration(uMolar)');%title(Input_insult)
if bSave
    figSaveEPS(gcf,[pth, filesep, 'conc_wavelength.eps'],'fixeps', true);
end
close(gcf);
[error_sel_2fit]= sel_error(err_OD_2(:,t),fig1,error_sel_3fit(1,:),S);% errors of 1st combination of wavelengths

if error_sel_2fit==0
    Selected_Residual.TwoFit.ErrorLambda= error_sel_2fit;
    disp('No points are selected for calculating residual!');
else % No wavelength combination is selected or all wavelengths are selected
    error_lambda2 = cat(2,lambda',error_sel_2fit(2:end,:));
    Selected_Residual.TwoFit.ErrorLambda= error_lambda2;
end
Selected_Residual.TwoFit.Time= error_sel_2fit(1,:);

Results.Selected_Residual=  Selected_Residual;
toc

end


function [error_sel,S]= sel_error(err,fig,x,d,res)
if nargin <3 %#ok<ALIGN>
    
    %S=input('Please enter the no. of datapoint you would like to select for calculating residual:(enter 0 to skip)  ');
    S=0; %for automation
    error_sel= zeros(size(err,1)+1,S);
    if S>=1 %#ok<ALIGN>
        % Save residuals for the user's selected datapoints
        for x= 1:1:S
            
            [sel_x,sel_y]=User_GUI(fig);%#ok<NASGU> % read data cursor
            error_sel(1,x)= round(sel_x);% Time stamp
            error_sel(2:end,x)= err(:,round(sel_x));
        end
    else error_sel=0;end
elseif d>=1
    error_sel= zeros(size(err,1)+1,d);
    error_sel(1,:)= x;% Time stamp
    error_sel(2:end,:)= err(:,x);
    S = d;
else  error_sel=0;end
end

function [x,y]=User_GUI(fig)
refreshdata
dcm_obj = datacursormode(fig);
set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')
disp('Please select the point(s), then press SPACE to continue.')
pause                            % Wait while the user does this.
c_info = getCursorInfo(dcm_obj);
set(c_info.Target,'LineWidth',2)  % Make selected line wider
x=c_info.Position(1);
y=c_info.Position(2);

end

function [conc_fit,err_OD,Hb_T,OD_est,R_fit]=OD2Conc(K,~,component,WavelengthsFit,OD_fit,~,~,stindex,endindex,~,pl,flag)

if flag==0
    c=abscoeff;     % extract absorption coefficient from function named 'abscoeff'
    c_lambda = c(:,1);
    kDPF=c(:,8);
    ABSMAT1_Intp= zeros(K,3);
    for i = 1:K
        ABSMAT1_Intp(i,1)=interp1(c_lambda,c(:,2),WavelengthsFit(i))*interp1(c_lambda,kDPF,WavelengthsFit(i));
        ABSMAT1_Intp(i,2)=interp1(c_lambda,c(:,3),WavelengthsFit(i))*interp1(c_lambda,kDPF,WavelengthsFit(i));
        ABSMAT1_Intp(i,3)=interp1(c_lambda,c(:,4),WavelengthsFit(i))*interp1(c_lambda,kDPF,WavelengthsFit(i));
    end
elseif flag ==1 % no wavelength dependenced factor included
    c=abscoeff1;     % extract absorption coefficient from function named 'abscoeff1'
    c_lambda = c(:,1);
    ABSMAT1_Intp= zeros(K,3);
    for i = 1:K
        ABSMAT1_Intp(i,1)=interp1(c_lambda,c(:,2),WavelengthsFit(i));
        ABSMAT1_Intp(i,2)=interp1(c_lambda,c(:,3),WavelengthsFit(i));
        ABSMAT1_Intp(i,3)=interp1(c_lambda,c(:,4),WavelengthsFit(i));
    end
end
%  figure;plot(WavelengthsFit,ABSMAT1,'xr',WavelengthsFit,ABSMAT1_Intp,'xb');
ABSMAT1=ABSMAT1_Intp;

if component == 2
    ABSMAT = ABSMAT1(:,1:2);
elseif  component == 3
    ABSMAT = ABSMAT1;
end
%PIMAT	=inv(ABSMAT'*ABSMAT)*ABSMAT';	    % pseudo inverse matrix     (mM . cm / OD)
PIMAT	=(ABSMAT'*ABSMAT)\ (ABSMAT');

conc_fit=PIMAT*OD_fit;% row1:Hb row2:HbO2 rwo3:Cyt  (mMpathlength)

if size(pl,2)==2
    conc_fit= conc_fit./repmat(pl(stindex:endindex,2)',size(ABSMAT,2),1);   % convert mMpathlength tp mM
elseif size(pl,2)==3 && component == 2
    conc_fit= conc_fit./pl(:,1:2)';
else %conc_fit= conc_fit./pl';
    conc_fit= conc_fit./repmat(pl,component,1);
end
conc_fit = conc_fit*1000;  % convert mM tp uM
Hb_T = conc_fit(1,:)+conc_fit(2,:);
%   if component == 2
%     figure;plot(t,conc_fit(1,:),'b',t,conc_fit(2,:),'r',t,Hb_T,'k');grid;
%     legend('HHB','HBO_2','HBT');
%   elseif component == 3
%     figure;plot(t,conc_fit(1,:),'b',t,conc_fit(2,:),'r',t,conc_fit(3,:),'g',t,Hb_T,'k');grid;
%     legend('HHB','HBO_2','CtOx','HBT');
%   end
%   xlabel('Time(min)');ylabel('Concentration(uMolar)');
%   title(['Change Concentration for wavelength ',int2str(lambda1),' to ',int2str(lambda2), ' nm'],'fontsize',12,'fontweight','b');

%--------------------------------------------------------------------------
% Estimation of attenuation
%   OD_est = ABSMAT*conc_fit.*repmat(pl(stindex:endindex,2)',K,1)/1000;
OD_est = ABSMAT*conc_fit.*repmat(pl(1,stindex:endindex),K,1)/1000;
%--------------------------------------------------------------------------
%   figure;subplot(2,1,1);plot(WavelengthsFit,OD_est);grid;xlabel('Wavelength(nm)');ylabel('Attenuation(OD)');
%   title(['Estimated Change Attenuation of wavelength ',int2str(lambda1),' to ',int2str(lambda2), ' nm for all samples'],'fontsize',12,'fontweight','b');
err_OD= OD_est - OD_fit;
R_fit = 1-sum(sum(err_OD.^2))./sum(sum((OD_fit-repmat(mean(OD_fit),size(OD_fit,1),1)).^2));
%   subplot(2,1,2); plot(WavelengthsFit,err_OD);grid;
%   title(['Residual of Change ODs considering ',int2str(component),' chromophores for all samples'],'fontsize',12,'fontweight','b');
%   xlabel('Wavelength(nm)');ylabel('Absolute Residual');
% err_OD_mean = mean(err_OD,1);err_OD_std = std(err_OD,0,1);
% figure;errorbar(t,err_OD_mean,err_OD_std,'-*r'); grid; xlabel('Time(min)');ylabel('error')
end

function isFieldResult = myIsField (inStruct, fieldName)
% inStruct is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
isFieldResult = 0;
f = fieldnames(inStruct(1));
for i=1:length(f)
    if(strcmp(f{i},strtrim(fieldName)))
        isFieldResult = 1;
        return;
    elseif isstruct(inStruct(1).(f{i}))
        isFieldResult = myIsField(inStruct(1).(f{i}), fieldName);
        if isFieldResult
            return;
        end
    end
end
end

function output2txt(Results,filename)
% This program generates text file.

% open the file with write permission
k = filename;
fid = fopen(k, 'w');
if myIsField(Results, 'pathlength840')==0
    Results.pathlength840 = zeros(1,size(Results.Time,2));
    Results.pathlength740 = Results.pathlength840;
    str4 = ['Pathlength =  ',int2str(Results.pathlength(1)),'cm'];
else
    str4 = [];
end

str2 = ['Time ','    Hb    ','   HbO2    ','    CtOx      ',...
    '   HbT      ',' HbDiff    ','Pathlength840   ',...
    'Pathlength740     ','Absolute Hb    '];

str3 = ['(min) ','(uMolar) ',' (uMolar) ','   (uMolar) ',...
    '   (uMolar) ','    (uMolar) ','       (cm)    ',...
    '      (cm)     ','      (uMolar)     '];

fprintf(fid,'\r\n');
fprintf(fid, '%s\r\n', str2);
fprintf(fid, '%s\r\n', str3);
fprintf(fid, '%s\r\n', str4);
fprintf(fid,'\r\n');
OUTPUT= cat(2,Results.Time',Results.Concentrations',...
    Results.pathlength840',Results.pathlength740',...
    Results.Abs_Hb760');
fprintf(fid, '%2i  %9.4f %10.4f%11.4f%12.4f%13.4f%14.4f%15.4f%16.4f\r\n', OUTPUT');
fprintf(fid,'\r\n');
fclose(fid);
end






