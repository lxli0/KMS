% -------------------  Calculation of KMS/b ratio -----------------------
% -----------------------------------------------------------------------

% This code calculates the KMS/b ratio under a given catalogsize.
% When the sample size is large, it can be time-consuming. 
% We are building a parallel version.


% INPUT
% catalog_size: the earthquake sizes of interest
% catime: the number of random catalogs generated for a given catalog size

%
% Reference: Linxuan Li, Gang Luo, and Mian Liu (in revision). The K–M
% slope: a potential supplement for b-value.
% CopyRight: Linxuan Li, Wuhan University
% Tested under Matlab_R2020b 
%

clc,clear
catalog_size=[3860 2750];
catime=[20 20];
for xtime=1:length(catalog_size)
    L=catalog_size(xtime);
    time_L=catime(xtime);
    b=1;
    Result_Poisson=[];
    for time=1:time_L
        Mmin=0;
        meandt=1;
        td=zeros(1,L);
        mM=10; %magnitude upper limit
        M=-1/b.*log10(1-rand(1,L)*(1-10^(-b*mM))); %random magnitude
        dt_Poisson=exprnd(meandt,1,L-1); %random time interval
        t_Poisson=cumsum(dt_Poisson);
        t_Poisson=[0,t_Poisson];
        Result_Poisson(time)=VGA(t_Poisson,M,L);
    end
    mout(xtime)=mean(Result_Poisson); %mean value
    sout(xtime)=std(Result_Poisson); %standard deviation
end


%% VGA
function [Slope] =VGA(Dt,Dm,L)
    VG=zeros(1,L);
    for i=1:L
        for j=1:L
            if i==j
                VG(i)=VG(i);
            elseif abs(i-j)==1
                VG(i)=VG(i)+1;
            elseif i-j>=2
                t=0;
                for k=j+1:i-1
                    if Dm(k)<Dm(j)+(Dm(i)-Dm(j))*(Dt(k)-Dt(j))/(Dt(i)-Dt(j))
                        t=t+1;
                    end
                end
                if t==(i-j-1)
                VG(i)=VG(i)+1;
                end     
            elseif j-i>=2
                t=0;
                for k=i+1:j-1
                    if Dm(k)<Dm(i)+(Dm(j)-Dm(i))*(Dt(k)-Dt(i))/(Dt(j)-Dt(i))
                        t=t+1;
                    end
                end
                if t==(j-1-i)
                VG(i)=VG(i)+1;
                end
            end
        end
    end
    p=polyfit(Dm,VG,1);
    Slope=p(1);
end
