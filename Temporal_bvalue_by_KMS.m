% -------  Image the Temporal Evolution of b-value by K-M slope (KMS)  -----------
% --------------------------------------------------------------------------------

% This code calculates the temporal evolution of b-value by least-squares
% regression (LSR), maximum likelihood estimation (MLE), K-M slope (KMS),
% and the improved KMS. Note: the magnitude of completeness (Mc) is considered
% constant for the entire catalog. Once consider the temporal evolution of
% Mc, parallel computing is required to obtain the time-varying KMS0.

% INPUT
% Catalog: the earthquake catalog used for analysis, the first column
% represents the magnitude, and the second column represents the time.
% bin: the events in a window used for calculation
% step: the number of events moved forward for each calculation

%
% Reference: Linxuan Li, Gang Luo, and Mian Liu (under review). The K–M
% slope: a potential supplement for b-value.
% CopyRight: Linxuan Li, Wuhan University
% Tested under Matlab_R2020b 
%

clc,clear
Catalog=load('Catalog_out.txt');
M=Catalog(:,4);
t=Catalog(:,1);
M=0.001*ones(length(M),1)+M; %To deal with the accuracy of data, use as appropriate
Mrange=min(M):0.1:max(M);
%% Mc
cntnumb=hist(M,Mrange);
[~,ind]=max(cntnumb);
Mmin=Mrange(ind)+0.2;            
jkf=find(M>=Mmin);
MM=M(jkf);
tt=t(jkf);
L0=length(MM);
tm=[tt,MM];
tm=sortrows(tm,1);

%% Divide time window
bin=200;
step=10;
num_bin=floor((L0-bin)/step)+1;
for i=1:100
    Pcm=[];
    t0=[];
    Pcm=rand(1,bin);
    M0=-log10(Pcm);
    dt0=exprnd(1,1,bin-1);
    t0=cumsum(dt0);
    t0=[0,t0];
    KMS0_0(i)=VGA(t0, M0,0);
end
KMS0=mean(KMS0_0);

%% Calculation
bleast=[];
bmax=[];
KMS_origin=[];
KMS_Poisson=[];
t_cal=[];

for j=1:num_bin
    event_begin=(j-1)*step+1;
    event_end=(j-1)*step+bin;
    t_origin_bin=tm(event_begin:event_end,1);
    m_bin=tm(event_begin:event_end,end);

    [bleast(j),bmax(j)]=GR(m_bin,Mmin);
    KMS_origin(j)=VGA(t_origin_bin,m_bin,Mmin);  

    for k=1:10
        mean_dt=mean(diff(t_origin_bin));
        dt_Poisson=exprnd(mean_dt,1,bin-1);
        t_Poisson0=cumsum(dt_Poisson); 
        t_Poisson_bin=[t_origin_bin(1),t_origin_bin(1)+t_Poisson0];
        randIndex_A = randperm(bin);
        m_Poisson_bin = m_bin(randIndex_A);
        KMS_Poisson0(k)=VGA(t_Poisson_bin, m_Poisson_bin,Mmin);
    end

    KMS_Poisson(j)=mean(KMS_Poisson0);
    t_cal(j)=t_origin_bin(end);     
end

bleast_normal=bleast;
bmax_normal=bmax;
KMS_origin_normal=KMS_origin/KMS0;
KMS_Poisson_normal=KMS_Poisson/KMS0;


%% Plot
color=1/255*[117 114 181;91 183 205; 197 86 89;203 180 123];
figure('units','normalized','position',[0.1,0.1,0.7,0.9])
subplot(5,1,1)
plot(t_cal,bleast_normal,'color',color(1,:),'linewidth',1.5);hold on;
plot(t_cal,bmax_normal,'color',color(2,:),'linewidth',1.5);hold on;
plot(t_cal,KMS_origin_normal,'color',color(3,:),'linewidth',1.5);hold on;
plot(t_cal,KMS_Poisson_normal,'color',color(4,:),'linewidth',1.5);hold on;
legend('LSR','MLE','KMS','KMS (Improved)','NumColumns',4,'Location','Northwest')
set(gca,'position',[0.15 0.82 0.75 0.16])
set(gca,'xlim',[1980 2021]);
set(gca,'xticklabel',[]);
ylabel('b value')
ylim([0.2,1.8])
box on;
grid on;
set(gca,'fontsize',16)

subplot(5,1,2)
plot(t_cal,bleast_normal,'color',color(1,:),'linewidth',1.5);
set(gca,'xticklabel',[]);
ylabel('b value')
ylim([0.2,1.8])
box on;
grid on;
set(gca,'fontsize',16)
set(gca,'position',[0.15 0.64 0.75 0.16])
set(gca,'xlim',[1980 2021]);
legend('LSR','Location','Northwest')

subplot(5,1,3)
plot(t_cal,bmax_normal,'color',color(2,:),'linewidth',1.5);
set(gca,'xticklabel',[]);
ylim([0.2,1.8])
ylabel('b value')
box on;
grid on;
set(gca,'fontsize',16)
set(gca,'position',[0.15 0.46 0.75 0.16])
set(gca,'xlim',[1980 2021]);
legend('MLE','Location','Northwest')

subplot(5,1,4)
plot(t_cal,KMS_origin_normal,'color',color(3,:),'linewidth',1.5);
ylabel('b value')
legend('KMS','Location','Northwest')
ylim([0.2,1.8])
box on;
grid on;
set(gca,'fontsize',16)
set(gca,'position',[0.15 0.28 0.75 0.16])
set(gca,'xlim',[1980 2021]);
set(gca,'xticklabel',[]);

subplot(5,1,5)
plot(t_cal,KMS_Poisson_normal,'color',color(4,:),'linewidth',1.5);
ylim([0.2,1.8])
box on;
grid on;
set(gca,'fontsize',16)
set(gca,'position',[0.15 0.1 0.75 0.16])
set(gca,'xlim',[1980 2021]);
xlabel('Year');
ylabel('b value')
legend('KMS (Improved)','Location','Northwest')
    
    
%% traditional b-value estimation
function [b_least,b_max]=GR(Dm,Mmin)
    L=length(Dm);
    dm=0.1;
    m=0:dm:floor(max(Dm)/dm)*dm;
    for i=1:floor(max(Dm)/dm)+1
        n0(i)=length(find(Dm>=dm*(i-1)&Dm<dm*i)); 
    end
    
    cn0(1)=L;
    for i=2:floor(max(Dm)/dm)+1
        cn0(i)=L-sum(n0(1:i-1));
    end
    n=log10(n0);
    cn=log10(cn0);
 
    jkf1=find(m>=Mmin);
    Mgemin=m(jkf1);
    LogNgemin=cn(jkf1);
    p=polyfit(Mgemin,LogNgemin,1);
    b_least=-p(1);
    
    jkf2=find(Dm>=Mmin);
    Mgemin2=Dm(jkf2);
    b_max=1/(log(10)*(mean(Mgemin2)-Mmin+0.05));
end

%% KMS estimation
function [Slope] =VGA(Dt0,Dm0,Mc)
    ij=find(Dm0>=Mc);
    Dt=Dt0(ij);
    Dm=Dm0(ij);
    L=length(Dm);
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
