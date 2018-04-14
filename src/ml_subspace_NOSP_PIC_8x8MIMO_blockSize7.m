% Test the subspace detection performance
addpath('/user/ufmc/xwang/MatlabWorks/MIMO-SD/')
clear all

Nt=8; % Tx Antenna
Nr=8; % Rx Antenna
N_ch = 1e4;
nb =2;
GrayLabeling=1;
mapPamQam =1; % QAM
[symlabel,QAM_Symbols] = get_mapping_and_labeling(nb,GrayLabeling,mapPamQam);
rx_active=[0 0 0 1 0]; % [ZF MMSE fullML, blockML fullAPP] 1 yes, 0 no

SNRs=-10:2.5:30;
Niter=10;
withPIC = 1;
useDistMetric = 1;
MLorAPP = 0;
maxLog=1;
La=0;
BlockSize=7;

ber_full_ml = zeros(length(SNRs),1);
ber_block_ml = zeros(length(SNRs),Niter+1);
ber_zf = zeros(length(SNRs),1);
ber_mmse = zeros(length(SNRs),1);
ber_full_App = zeros(length(SNRs),1);
MI_full_ml = zeros(length(SNRs),1);
MI_block_ml = zeros(length(SNRs),Niter+1);
MI_zf_hard = zeros(length(SNRs),1);
MI_mmse_hard = zeros(length(SNRs),1);
MI_full_App = zeros(length(SNRs),1);
for iSNR=1:length(SNRs)
    noise_var = 10^(-SNRs(iSNR)/10);
    brrNr_fullMl =0;
    brrNr_blockMl =zeros(Niter+1,1);
    brrNr_ZF =0;
    brrNr_MMSE =0;
    brrNr_fullAPP =0;
    tx_bits_all=-ones(nb*N_ch,Nt);
    rx_bits_zf_hard_all =-ones(nb*N_ch,Nt);
    rx_bits_mmse_hard_all = -ones(nb*N_ch,Nt);
    rx_bits_full_ML_all =-ones(nb*N_ch,Nt);
    rx_bits_block_ML_all =-ones(nb*N_ch,Nt,Niter+1);
    rx_bits_full_APP_all =-ones(nb*N_ch,Nt);
    LLR_full_APP_all =-ones(nb*N_ch,Nt);
    for iCh=1:N_ch
        %symbol generation
        sym_idx = randi(4,Nt,1) -1; % generate QPSK symbols
        tx_bits = de2bi(sym_idx,nb);
        tx_bits_all((1+(iCh-1)*nb):(iCh*nb),:) =tx_bits';
        s = mapping(symlabel,QAM_Symbols,sym_idx);
        
        % channel matrix
        H=1/sqrt(2)/sqrt(Nt) *(randn(Nr,Nt)+1i*randn(Nr,Nt));
        n=1/sqrt(2) *sqrt(noise_var)*(randn(Nr,1)+1i*randn(Nr,1));
        
        r = H*s+n; % received signal vector
        
        if(rx_active(1))
            %ZF rx
            s_hat = H\r;
            rx_bits_zf_hard  =hardDemapping(symlabel,QAM_Symbols,s_hat);
            rx_bits_zf_hard_all((1+(iCh-1)*nb):(iCh*nb),:) =rx_bits_zf_hard';
            brrNr_ZF = brrNr_ZF+sum(sum(double((tx_bits~=rx_bits_zf_hard))));
        end
        
        if(rx_active(2))
            %MMSE rx
            s_hat_mmse = (((H')*H+noise_var*eye(Nt))\(H'))*r;
            rx_bits_mmse_hard  =hardDemapping(symlabel,QAM_Symbols,s_hat_mmse);
            rx_bits_mmse_hard_all((1+(iCh-1)*nb):(iCh*nb),:) =rx_bits_mmse_hard';
            brrNr_MMSE = brrNr_MMSE+sum(sum(double((tx_bits~=rx_bits_mmse_hard))));
        end
        
        if(rx_active(3))
            %full ML
            s_hat_ml = ML_rx(r,H,QAM_Symbols);
            rx_bits_ml_hard  =hardDemapping(symlabel,QAM_Symbols,s_hat_ml);
            rx_bits_full_ML_all((1+(iCh-1)*nb):(iCh*nb),:) =rx_bits_ml_hard';
            brrNr_fullMl = brrNr_fullMl + sum(sum(double((tx_bits~=rx_bits_ml_hard))));
        end
        
        %block ML search
        if(rx_active(4))
            s_block_all_iterations = SubspaceMIMORx( r,H,BlockSize,QAM_Symbols,withPIC,Niter,useDistMetric,MLorAPP,symlabel,noise_var,maxLog,La,0,0);
            for iter = 1:(Niter+1)
               rx_bits_BlockML_hard  =hardDemapping(symlabel,QAM_Symbols,s_block_all_iterations(:,iter));
               rx_bits_block_ML_all((1+(iCh-1)*nb):(iCh*nb),:,iter) =rx_bits_BlockML_hard';
               brrNr_blockMl(iter,1) = brrNr_blockMl(iter,1)+sum(sum(double((tx_bits~=rx_bits_BlockML_hard))));
            end
        end
        % full APP
        if(rx_active(5))
            LLR_full_APP_vec = softDemapping(symlabel,QAM_Symbols,r,noise_var,2,H,0);
            LLR_full_APP = reshape(LLR_full_APP_vec,[nb Nt]);
            LLR_full_APP_all((1+(iCh-1)*nb):(iCh*nb),:) =LLR_full_APP;
            rx_bits_app_hard  =((LLR_full_APP>0)+0)';
            rx_bits_full_APP_all((1+(iCh-1)*nb):(iCh*nb),:) =rx_bits_app_hard';
            brrNr_fullAPP = brrNr_fullAPP + sum(sum(double((tx_bits~=rx_bits_app_hard))));
        end
    end
    % add up mutual info for all streams
    for iStream=1:Nt
        if(rx_active(1))
            MI_zf_hard(iSNR) = MI_zf_hard(iSNR) + getMutualInformationHardInput(tx_bits_all(:,iStream),rx_bits_zf_hard_all(:,iStream),nb);
        end
        if(rx_active(2))
            MI_mmse_hard(iSNR) = MI_mmse_hard(iSNR) + getMutualInformationHardInput(tx_bits_all(:,iStream),rx_bits_mmse_hard_all(:,iStream),nb);
        end
        if(rx_active(3))
            MI_full_ml(iSNR) = MI_full_ml(iSNR) + getMutualInformationHardInput(tx_bits_all(:,iStream),rx_bits_full_ML_all(:,iStream),nb);
        end
        if(rx_active(4))
            for iter=1:(Niter+1)
                MI_block_ml(iSNR,iter) = MI_block_ml(iSNR,iter) + getMutualInformationHardInput(tx_bits_all(:,iStream),rx_bits_block_ML_all(:,iStream,iter),nb);
            end
        end
        if(rx_active(5))
            MI_full_App(iSNR) = MI_full_App(iSNR) + getMutualInformationSoftInput(tx_bits_all(:,iStream),LLR_full_APP_all(:,iStream),nb);
        end
    end
    iSNR
    if(rx_active(1))
        ber_zf(iSNR) = brrNr_ZF/(nb*Nt*N_ch);
    end
    if(rx_active(2))
        ber_mmse(iSNR) = brrNr_MMSE/(nb*Nt*N_ch);
    end
    if(rx_active(3))
        ber_full_ml(iSNR) = brrNr_fullMl/(nb*Nt*N_ch);
    end
    if(rx_active(4))
        ber_block_ml(iSNR,:) = brrNr_blockMl/(nb*Nt*N_ch);
    end
    if(rx_active(5))
        ber_full_App(iSNR) = brrNr_fullAPP/(nb*Nt*N_ch);
    end
end

save(['RES_BER_Block_ML_NOSP_iterative_PIC_',num2str(Nr),'x',num2str(Nt),'MIMO','_blockSize',num2str(BlockSize),'.mat'],'ber_full_ml','ber_block_ml','ber_zf','ber_full_App','MI_zf_hard','MI_full_ml','MI_block_ml','MI_full_App')