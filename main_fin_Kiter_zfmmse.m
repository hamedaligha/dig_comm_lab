clear all
% close all
clc
SNR_list=0:5:20;
phi1_2=pi/3;  h=2.5;  Pt=1e5;
% U_locs=[-1.5 1.5 0;-1.5 1.4 0;1.5 2.2 0;0.5 2 0];
% U_locs=[0.5 1.5 0;1.5 0.5 0;2 1 0];
K_all = [2,4,6];
meanh=0;
iter_num_F=20; iter_num=20;
T=10;
user_x=-1.9:.1:1.9;
user_y=user_x;
[ux0,uy0] = meshgrid(user_x);
user_loc_all=[ux0(:), uy0(:)];
user_loc_all(:,3)=0;

dx=0.05;  drx=0.05;
Nt=36; 
Nrx=10;  Nrz=6;
N=Nrz*Nrx;    % K=size(U_locs,1);  
Bsx=-((sqrt(Nt)-1)/2*dx):dx:((sqrt(Nt)-1)/2*dx);
Bsy=Bsx;
[sx0,sy0] = meshgrid(Bsx);
BS_loc=[sx0(:), sy0(:) h*ones(Nt,1)];
RISx=-(((Nrx)-1)/2*drx):drx:(((Nrx)-1)/2*drx);
RISz=-(((Nrz)-1)/2*drx):drx:(((Nrz)-1)/2*drx);

RISy=2;
[sx1,sz1] = meshgrid(RISx,RISz);
RIS_loc=[sx1(:), RISy*ones(N,1), sz1(:)+h/2];
for k_iter = 1:length(K_all)
    K = K_all(k_iter)
    for iter_user=1:1000
%         tic
        user_index = randperm(size(user_loc_all,1),K);        
        U_locs=user_loc_all(user_index,:);        
        for k=1:K
            U_loc=U_locs(k,:);
            H_los(:,k)=Pt*h_fun_los(BS_loc,U_loc,phi1_2);
            H_nlos(:,:,k)=Pt*h_fun_nlos(BS_loc,RIS_loc,U_loc,phi1_2);
        end
        for iter_snr=1:length(SNR_list)
            
            SNRdb=SNR_list(iter_snr);
%             [iter_user SNRdb k_iter]
            SNR=10^(SNRdb/10);
            sigm2_n=1/SNR/Nt;
            for iter=1:iter_num_F
                ranp=reshape(randperm(N),K,[]);
                F=zeros(N,K);
                for k=1:K
                    F(ranp(k,:),k)=1;
                    hf2(:,k)=H_nlos(:,:,k)'*F(:,k);
                end
                Htot=H_los+hf2;
                pd=sigm2_n*K;
                Wzf=Htot*inv(Htot'*Htot);
                Wzf=manifold_projection(Wzf);
                SINRzf_all(iter)=sinr_fun_allk2(Wzf,Htot,sigm2_n);
                WzfMMSE=Htot*inv(Htot'*Htot+pd*eye(K));
                WzfMMSE=manifold_projection(WzfMMSE);
                SINRMMSE_all(iter)=sinr_fun_allk2(WzfMMSE,Htot,sigm2_n);
%                 FF=F;
%                 mu0=0.0005;
%                 W0=2*(rand(Nt,K)-0.5);
%                 Wn=manifold_projection(W0);
%                 for iterwf=1:10
%                     W0=2*(rand(Nt,K)-0.5);
%                     Wn=manifold_projection(W0);
%                     Htot=Htot_fun(H_los,H_nlos,FF);
%                     [Wn SINR_al objf]=W_finder4(Wn,Htot,sigm2_n,mu0);
%                     SINR_all(iterwf)=SINR_al(end);
%                     FF=F_finder(Wn,H_los,H_nlos,sigm2_n);
%                 end
%                 SINRnew_all(iter)=max(SINR_all);
%                 [Wn0 SINR_al0 objf0]=W_finder4(W0,H_los,sigm2_n,mu0);
%                 SINRnew0_all(iter)=SINR_al0(end);
            end
            SINRzf_all_fin(iter_snr)=mean(SINRzf_all);
            SINRMMSE_all_fin(iter_snr)=mean(SINRMMSE_all);
%             SINRnew_all_fin(iter_snr)=mean(SINRnew_all);
%             SINRnew0_all_fin(iter_snr)=mean(SINRnew0_all);
        end
        SINR_all_fin(:,:,iter_user)=[SINRzf_all_fin;SINRMMSE_all_fin];
%         if rem(iter_user,10)==0
%             disp(log10(mean(SINR_all_fin,3)))
%         end
%         toc
    end
    SINR_all_fin_mean=mean(SINR_all_fin,3);
    SINR_all_fin_mean_kall(:,:,k_iter) = SINR_all_fin_mean
    save('SINR_all_fin_mean_kall_zfmmse','SINR_all_fin_mean_kall');
    figure
    semilogy(SNR_list,SINR_all_fin_mean(1,:),'-*r',SNR_list,SINR_all_fin_mean(2,:),'-+b')
    grid on
    legend('ZF','MMSE')
end

