% -----------  B-value estimation by the improved K-M slope (KMS) method  ------------
% ------------------------------------------------------------------------------------

% This script is the function to calculate the b-value by the improved K-M slope (KMS) method
% The calculation is based on the traditional definition of a visibility graph. 
% The input is a series of Magnitude filted by Mc
% Important note: Magnitude should be a row vector

% Reference: 
% 1. Linxuan Li, Gang Luo, and Mian Liu (2023). The K–M slope: a potential
% supplement for b-value. Seismological Research Letters.
% 2. Linxuan Li and Gang Luo (submitted). Can we obtain reliable seismic b-values for real-time catalogs?

% CopyRight: Linxuan Li, Wuhan University (lucas.linxuan.li@gmail.com)
% Tested under Matlab_R2020b

%%
clc,clear
% replace the part by your data
    L=500;
    b=1;
    Pcm=rand(1,L);
    Magnitude=-1/b.*log10(Pcm);
    b_KMS=KMS(Magnitude);
%
function [b_kms] =KMS(Magn)
    size=length(Magn);
    kms0=-15.15*(log10(size))^(-2.14)+11.85;
    time_random=10;% the number of synthetic catalogs for the improved KMS method
    kms=zeros(1,time_random);
    for l=1:time_random
        dt_Poisson=exprnd(1,1,size-1);
        t_Poisson=cumsum(dt_Poisson);
        t_Poisson=[0,t_Poisson];
        randIndex_A = randperm(size);
        Magn_cal= Magn(randIndex_A);
        kms(l)=VGA(t_Poisson,Magn_cal);
    end
    im_kms=mean(kms);
    b_kms=im_kms/kms0;
end

function [Slope] =VGA(Dt,Dm)
    L=length(Dm);
    VG=2*ones(1,L);
    VG(1)=1;
    VG(end)=1;
    for i=1:L-2
        for j=i+2:L
            tem_Dm=Dm(i+1:j-1);
            tem_Dt=Dt(i+1:j-1);
            cri_Dm=Dm(i)+(Dm(j)-Dm(i))*(tem_Dt-Dt(i))/(Dt(j)-Dt(i));
            cri=tem_Dm-cri_Dm;
            if max(cri)<0
                VG(i)=VG(i)+1;
                VG(j)=VG(j)+1;
            end
        end
    end
    p=polyfit(Dm,VG,1);
    Slope=p(1);
end
