% Test the subspace detection performance
% addpath('C:\Users\vivianma\Documents\MATLAB\Forschungsarbeit')
% 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 19 20
clear all
% ------------------------------------------------------------------------------
% some configuration for debuging
ON = true;
OFF = false;
DEBUG = ON;
OUTPUT = OFF;
TEST = ON;

% ------------------------------------------------------------------------------
% global variables
Nt = 2; % Tx Antenna
Nr = 2; % Rx Antenna
N_ch = 10;
nb = 1;
GrayLabeling=1; % PSK
% mapPamQam =1; % QAM
mapPamQam = 0;
[symlabel,QAM_Symbols] = get_mapping_and_labeling(nb,GrayLabeling,mapPamQam);
% rx_active=[ 0 1 0]; % [ZF MMSE ML] 1 yes, 0 no
%% Element definition

SNRs = -10:2.5:14; % use array not for loop
Niter = 5;

ber_zf = zeros(length(SNRs),1);
ser_zf = zeros(length(SNRs),1);
ber_zf_soft = zeros(length(SNRs),1);

ber_mmse = zeros(length(SNRs),1);
ber_mmse_soft = zeros(length(SNRs),1);
ser_mmse= zeros(length(SNRs),1);

ber_ml= zeros(length(SNRs),1);
ser_ml= zeros(length(SNRs),1);

ber_sic= zeros(length(SNRs),1);
ser_sic= zeros(length(SNRs),1);

ber_sphere = zeros(length(SNRs),1);
ber_full_App = zeros(length(SNRs),1);

ber_BF_MRF = zeros(length(SNRs),1);
ber_BF_GAMP = zeros(length(SNRs),1);

MI_zf_bit_hard = zeros(length(SNRs),1);% Mutual Information Matrix
MI_zf_symbol_hard = zeros(length(SNRs),1);
MI_zf_soft_all = zeros(length(SNRs),1);
MI_zf_soft_noscaling = zeros(length(SNRs),N_ch);
MI_zf_soft_noscaling_all = zeros(length(SNRs),1);

MI_mmse_hard = zeros(length(SNRs),1);
MI_mmse_soft_all = zeros(length(SNRs),1);
MI_mmse_soft_scaled = zeros(length(SNRs),N_ch);
MI_mmse_soft_scaled_all = zeros(length(SNRs),1);
MI_mmse_soft_unbiased = zeros(length(SNRs),N_ch);
MI_mmse_soft_unbiased_all = zeros(length(SNRs),1);
MI_mmse_soft_noscaling = zeros(length(SNRs),N_ch);
MI_mmse_soft_noscaling_all = zeros(length(SNRs),1);

MI_ml = zeros(length(SNRs),1);
MI_sic_hard = zeros(length(SNRs),1);
MI_sphere_hard = zeros(length(SNRs),1);

MI_full_App = zeros(length(SNRs),1);
MI_bf_bit_MRF = zeros(length(SNRs),1);
MI_bf_bit_GAMP = zeros(length(SNRs),1);
% MI_CM_ZF_all = zeros(length(SNRs),1);
% MI_CM_ZF = zeros(length(SNRs),1);
% MI_CM_MMSE_all = zeros(length(SNRs),1);

% dmin_BPSK = zeros(length(SNRs),1);
ber_siso_all = zeros(length(SNRs),1);

% ------------------------------------------------------------------------------
%% Main part
for iSNR=1:length(SNRs)
	% some local variables
	noise_var = 10^(-SNRs(iSNR)/10);
	brrNr_siso = 0;
	symerrNr_siso = 0;
	dmin_BPSK = 0;

	brrNr_ZF =0;% number of bit error
	symerrNr_ZF=0;% number of symbol error
	brrNr_zf_soft = 0;
	brrNr_ZF_Sphere = 0;

	brrNr_MMSE =0;
	brrNr_MMSE_soft = 0;
	brrNr_mmse_soft = 0;
	symerrNr_MMSE=0;

	brrNr_ML=0;
	brrNr_ml_soft = 0;
	symerrNr_ML=0;

	brrNr_SIC=0;
	symerrNr_SIC=0;

	brrNr_fullAPP = 0;

	brrNr_BF_MRF = 0;
	brrNr_BF_GAMP = 0;
	N_sym_per_Ch = 1;
	channel_weight = 1;

	ZForMMSE=1;
	ordering=1;
	algorithm=1;
	La=0;

	tx_bits_SISO_all = -ones(N_ch,nb);
	tx_bits_all_soft = -ones(0,1);

	tx_symbols_all = -ones(N_ch,Nt);
	rx_symbol_zf_hard_all = -ones(N_ch,Nt);
	rx_symbol_mmse_hard_all = -ones(N_ch,Nt);
	rx_symbol_ml_all = -ones(N_ch,Nt);
	rx_symbol_sic_hard_all = -ones(N_ch,Nt);
	rx_sphere_hard_all = -ones(N_ch,Nt);

	tx_bits_all=-ones(nb*N_ch,Nt);
	rx_bits_zf_hard_all =-ones(nb*N_ch,Nt);
	rx_bits_mmse_hard_all = -ones(nb*N_ch,Nt);
	rx_bits_ml_all =-ones(nb*N_ch,Nt);
	rx_bits_sic_hard_all = -ones(nb*N_ch,Nt);
	rx_bits_bf_hard_all = -ones(nb*N_ch,Nt);
	rx_bits_bf_GAMP_all = -ones(nb*N_ch,Nt);

	rx_bits_zf_soft_all_channel = -ones(0,1);
	rx_bits_zf_soft_all_noscaling_channel = -ones(0,1);
	rx_bits_mmse_soft_unbiased_all_channel = -ones(0,1);
	rx_bits_mmse_soft_scaled_all_channel = -ones(0,1);
	rx_bits_mmse_soft_all_noscaling_channel = -ones(0,1);

	dmin_BPSK= 2/noise_var;
	ber_siso = 1/2*erfc(dmin_BPSK);
	ber_siso_all(iSNR,1) = ber_siso;

	% --------------------------------------------------------------------------
	% one SNR value corresponding to 1e4 times sending
	for iCh = 1:N_ch
		% symbol generation MIMO
		%       sym_idx = randi(4,Nt,1) -1; % generate QPSK symbols 4 symbol from 0 to 3
		sym_idx = randi(2,Nt,1) -1;
		tx_bits = de2bi(sym_idx,nb); % converts a nonnegative decimal integer to a binary row vector
		tx_bits_all((1+(iCh-1)*nb):(iCh*nb),:) = tx_bits'; %'
		% put tx_bits in tx_bits_all
		s = mapping(symlabel,QAM_Symbols,sym_idx);% original symbol in constellation map
		tx_symbols_all(iCh,:) =s.'; %'

		% ----------------------------------------------------------------------
		%% symbol generation SISO
		sym_SISO_idx = randi(2,1,1) -1;
		tx_SISO_bits = de2bi(sym_SISO_idx,nb);
		tx_bits_SISO_all(iCh,:) =  tx_SISO_bits;
		s_SISO = mapping(symlabel,QAM_Symbols, sym_SISO_idx);
		tx_symbols_SISO_all(iCh,:) = s_SISO;

		% ----------------------------------------------------------------------
		% MIMO channel matrix
		% generate Nt*Nr H randn create mean=0,delta=1 matrix
		H = 1/sqrt(2)/sqrt(Nt) *(randn(Nr,Nt)+1i*randn(Nr,Nt));
		n = 1/sqrt(2) *sqrt(noise_var)*(randn(Nr,1)+1i*randn(Nr,1));
		r = H*s+n; % received signal vector

		% ZF rx
		% ----------------------------------------------------------------------
		% SISO channel matrix
		% n_SISO = 1/sqrt(2)*sqrt(noise_var)*(randn(1)+1i*randn(1));
		% r_SISO = s_SISO + n_SISO;
		%
		% s_est_SISO = r_SISO;
		% [rx_bits_siso_hard,rx_symbol_siso_hard] = hardDemapping(symlabel,QAM_Symbols,s_est_SISO);
		% rx_bits_siso_hard_all(iCh,:) = rx_bits_siso_hard;
		% rx_symbol_siso_hard_all(iCh,:)= rx_symbol_siso_hard;
		% brrNr_siso = brrNr_siso + sum(sum(double(tx_SISO_bits~= rx_bits_siso_hard)));% sum of bit error
		% symerrNr_siso = symerrNr_siso + sum(sum(double(s_SISO~= rx_symbol_siso_hard)));
		% ber_siso(iSNR) = brrNr_siso/(nb*N_ch);
		% ser_siso(iSNR) = symerrNr_siso/(N_ch);

		% ZF RECEIVER
		% ----------------------------------------------------------------------
		% harddemapping
		s_est_zf = H\r; % (inverser H)*r
		[rx_bits_zf_hard,rx_symbol_zf_hard] = hardDemapping(symlabel,QAM_Symbols,s_est_zf);
		rx_bits_zf_hard_all((1+(iCh-1)*nb):(iCh*nb),:) = rx_bits_zf_hard'; %'
		rx_symbol_zf_hard_all(iCh,:) = rx_symbol_zf_hard.'; %'
		% sum of bit error
		brrNr_ZF = brrNr_ZF+sum(sum(double(tx_bits~= rx_bits_zf_hard)));
		symerrNr_ZF = symerrNr_ZF+sum(sum(double(s ~= rx_symbol_zf_hard.'))); %'

		% ----------------------------------------------------------------------
		% softdemapping
		s_est_zf = reshape(s_est_zf,Nt*N_sym_per_Ch,1);
		rx_bits_zf_soft  = softDemappingSISO(symlabel,QAM_Symbols,s_est_zf,1,2,ones(Nt*N_sym_per_Ch,1),0);
		rx_bits_zf_soft_all_noscaling = rx_bits_zf_soft;
		noise_scaling_zf = abs(diag((H'*H)\eye(Nt)))/Nt; %'
		% noise_scaling_factor_all(:,iCh) = noise_scaling;
		rx_bits_zf_soft = rx_bits_zf_soft./(repmat(noise_scaling_zf,[N_sym_per_Ch nb]));
		rx_bits_zf_soft_all = rx_bits_zf_soft;
		% rx_bits_zf_soft_all_channel = [rx_bits_zf_soft_all_channel;rx_bits_zf_soft_all(:)];
		% rx_bits_zf_soft_all_noscaling_channel = [rx_bits_zf_soft_all_noscaling_channel;rx_bits_zf_soft_all_noscaling(:)];
		rx_bits_zf_hard = logical((sign(rx_bits_zf_soft)+1)/2);
		brrNr_zf_soft = brrNr_zf_soft+sum(sum(double((tx_bits~=rx_bits_zf_hard))));
		%  noise_scaling_zf_all((1+(iCh-1)*nb):(iCh*nb),:)  = noise_scaling_zf*noise_var;

		% ----------------------------------------------------------------------
		% ZF Sphere Detection
		r_sphere = sphdec(H,r,QAM_Symbols);
		rx_sphere_hard_all(iCh,:) = r_sphere.'; %'
		[rx_bits_sphere_hard,rx_symbol_sphere_hard] = hardDemapping(symlabel,QAM_Symbols,r_sphere);
		rx_bits_sphere_hard_all((1+(iCh-1)*nb):(iCh*nb),:) = rx_bits_sphere_hard'; %'
		brrNr_ZF_Sphere = brrNr_ZF_Sphere + sum(sum(double(tx_bits~= rx_bits_sphere_hard)));

		% if (rx_active(2))
			% ------------------------------------------------------------------
			%  MMSE
			W_mmse = (((H')*H+noise_var*eye(Nt))\(H'));
			s_est_mmse = W_mmse*r;
			[rx_bits_mmse_hard,rx_symbol_mmse_hard] = hardDemapping(symlabel,QAM_Symbols,s_est_mmse);
			rx_bits_mmse_hard_all((1+(iCh-1)*nb):(iCh*nb),:) = rx_bits_mmse_hard'; %'
			rx_symbol_mmse_hard_all(iCh,:)= rx_symbol_mmse_hard;
			brrNr_MMSE = brrNr_MMSE+sum(sum(double((tx_bits~= rx_bits_mmse_hard))));% sum of bit error
			% symerrNr_MMSE = symerrNr_MMSE+sum(sum(double((s.'~= rx_symbol_mmse_hard)))); %'

			W_mmse = (((H')*H+noise_var*eye(Nt))\(H'));
			s_hat_mmse = W_mmse*r;
			rx_post_mat = W_mmse*H;
			unbiased_scaling = diag(rx_post_mat);
			s_est_mmse = reshape(s_est_mmse,Nt*N_sym_per_Ch,1);
			rx_bits_mmse_soft = softDemappingSISO(symlabel,QAM_Symbols,s_est_mmse,1,2,ones(Nt*N_sym_per_Ch,1),0);
			rx_bits_mmse_soft_all_noscaling = rx_bits_mmse_soft;
			s_est_mmse = s_est_mmse./repmat(unbiased_scaling,[N_sym_per_Ch 1]);
			rx_bits_mmse_soft = softDemappingSISO(symlabel,QAM_Symbols,s_est_mmse,1,2,ones(Nt*N_sym_per_Ch,1),0);
			rx_bits_mmse_soft_unbiased = rx_bits_mmse_soft;
			Interference_power = sum((abs(rx_post_mat - diag(unbiased_scaling))).^2,2);
			noise_power = abs(diag(W_mmse*W_mmse'))*noise_var; %'
			noise_scaling_mmse = (Interference_power+noise_power)./(abs(unbiased_scaling).^2);
			rx_bits_mmse_soft_scaled = rx_bits_mmse_soft./(repmat(noise_scaling_mmse,[N_sym_per_Ch nb]));
			rx_bits_mmse_hard = logical((sign(rx_bits_mmse_soft)+1)/2);
			brrNr_mmse_soft = brrNr_mmse_soft+sum(sum(double((tx_bits~=rx_bits_mmse_hard))));

		% end % if (rx_active(2))

		% if(rx_active(3))
			% full ML
			s_est_ml = ML_rx(r,H,QAM_Symbols);
			[rx_bits_ml,rx_symbol_ml] = hardDemapping(symlabel,QAM_Symbols,s_est_ml);
			rx_bits_ml_all((1+(iCh-1)*nb):(iCh*nb),:) =rx_bits_ml'; %'
			rx_symbol_ml_all(iCh,:) = rx_symbol_ml;
			brrNr_ML = brrNr_ML + sum(sum(double((tx_bits~=rx_bits_ml))));
			symerrNr_ML = symerrNr_ML+sum(sum(double((s.'~=rx_symbol_ml)))); %'
		% end % if(rx_active(3))
		% {
		%    if (rx_active(2))
		%        ZForMMSE=1;
		%    else
		%        ZForMMSE=0;
		%    end % if (rx_active(2))
		% }

		% ----------------------------------------------------------------------
		% APP
		LLR_full_APP_vec = softDemapping(symlabel,QAM_Symbols,r,noise_var,2,H,0);
		LLR_full_APP = reshape(LLR_full_APP_vec,[nb Nt]);
		LLR_full_APP_all((1+(iCh-1)*nb):(iCh*nb),:) =LLR_full_APP;
		rx_bits_app_hard = ((LLR_full_APP>0)+0)'; %'
		rx_bits_full_APP_all((1+(iCh-1)*nb):(iCh*nb),:) =rx_bits_app_hard'; %'
		brrNr_fullAPP = brrNr_fullAPP + sum(sum(double((tx_bits~=rx_bits_app_hard))));

		% ----------------------------------------------------------------------
		% ML softdemapping
		LLR_ml_bits = softDemapping(symlabel,QAM_Symbols,s_est_zf,noise_var,algorithm,channel_weight,La);
		LLR_ml_bits_all((1+(iCh-1)*nb):(iCh*nb),:) = LLR_ml_bits;
		rx_bits_ml_llr_hard =((LLR_ml_bits>0)+0)'; %'
		rx_bits_ml_llr_hard_all((1+(iCh-1)*nb):(iCh*nb),:) = rx_bits_ml_llr_hard'; %'
		brrNr_ml_soft = brrNr_ml_soft + sum(sum(double((tx_bits~=rx_bits_ml_llr_hard))));

		% ----------------------------------------------------------------------
		% SIC
		[rx_symbol_sic_hard,~,rx_bits_sic_hard] = mimoSICrx(H,r,noise_var,ZForMMSE,ordering,QAM_Symbols,symlabel);
		rx_symbol_sic_hard_all(iCh,:)= rx_symbol_sic_hard.'; %'
		rx_bits_sic_hard_all((1+(iCh-1)*nb):(iCh*nb),:) = rx_bits_sic_hard'; %'
		symerrNr_SIC=symerrNr_SIC+sum(sum(double((s~=rx_symbol_sic_hard))));
		brrNr_SIC = brrNr_SIC + sum(sum(double((tx_bits~=rx_bits_sic_hard))));

		tx_bits_all_soft = [tx_bits_all_soft;tx_bits(:)];
		rx_bits_zf_soft_all_channel = [rx_bits_zf_soft_all_channel;rx_bits_zf_soft_all(:)];
		rx_bits_zf_soft_all_noscaling_channel = [rx_bits_zf_soft_all_noscaling_channel;rx_bits_zf_soft_all_noscaling(:)];
		rx_bits_mmse_soft_unbiased_all_channel = [rx_bits_mmse_soft_unbiased_all_channel;rx_bits_mmse_soft_unbiased(:)];
		rx_bits_mmse_soft_scaled_all_channel = [rx_bits_mmse_soft_scaled_all_channel;rx_bits_mmse_soft_scaled(:)];
		rx_bits_mmse_soft_all_noscaling_channel = [rx_bits_mmse_soft_all_noscaling_channel;rx_bits_mmse_soft_all_noscaling(:)];

		% ----------------------------------------------------------------------
		% belief propagation
		% Markov random field
		LLR_BF_MRF = mimo_BP_MRF_detector_v02(r,H,noise_var,QAM_Symbols,symlabel,Niter,0.2);
		if TEST
			% r,H,noise_var,QAM_Symbols,symlabel,Niter,0.2
			LLR_BF_MRF
			TEST = OFF;
		end
		rx_bits_bf_hard = logical((1+sign(LLR_BF_MRF(:,:,Niter+1)))/2);
		rx_bits_bf_hard_all((1+(iCh-1)*nb):(iCh*nb),:) = rx_bits_bf_hard'; %'
		brrNr_BF_MRF = brrNr_BF_MRF + sum(sum(double((tx_bits~= rx_bits_bf_hard))));
		% brrNr_BF_MRF  = brrNr_BF_MRF (:,:); %??ά
		% semilogy(1:Niter+1,ber_BF_MRF (1,:)),

		% ----------------------------------------------------------------------
		% Gaussian Approximation Message Passing
		% LLR_BF_GAMP = mimo_BP_GAMP_detector(r,H,noise_var,QAM_Symbols,symlabel,Niter);
		% rx_bits_bf_GAMP = logical((1+sign(LLR_BF_GAMP(:,:,Niter+1)))/2);
		% rx_bits_bf_GAMP_all((1+(iCh-1)*nb):(iCh*nb),:) = rx_bits_bf_Gaus'; %'
		% brrNr_BF_GAMP = brrNr_BF_GAMP + sum(sum(double((tx_bits~= rx_bits_bf_GAMP))));

	end

	% --------------------------------------------------------------------------
	% Soft Mutual information
	MI_zf_soft_all(iSNR,1) = MI_zf_soft_all(iSNR,1) + getMutualInformationSoftInput(tx_bits_all_soft(:),rx_bits_zf_soft_all_channel(:),nb*Nt);
	MI_zf_soft_noscaling_all(iSNR,1) = MI_zf_soft_noscaling_all(iSNR,1) + getMutualInformationSoftInput(tx_bits_all_soft(:),rx_bits_zf_soft_all_noscaling_channel(:),nb*Nt);
	MI_mmse_soft_unbiased_all(iSNR,1) = MI_mmse_soft_unbiased_all(iSNR,1) + getMutualInformationSoftInput(tx_bits_all_soft(:),rx_bits_mmse_soft_unbiased_all_channel(:),nb*Nt);
	MI_mmse_soft_scaled_all(iSNR,1) = MI_mmse_soft_scaled_all(iSNR,1) + getMutualInformationSoftInput(tx_bits_all_soft(:),rx_bits_mmse_soft_scaled_all_channel(:),nb*Nt);
	MI_mmse_soft_noscaling_all(iSNR,1) = MI_mmse_soft_noscaling_all(iSNR,1) + getMutualInformationSoftInput(tx_bits_all_soft(:),rx_bits_mmse_soft_all_noscaling_channel(:),nb*Nt);

	% --------------------------------------------------------------------------
	% add up mutual info for all streams
	for iStream=1:Nt

		MI_zf_bit_hard(iSNR) = MI_zf_bit_hard(iSNR) + getMutualInformationHardInput(tx_bits_all(:,iStream),rx_bits_zf_hard_all(:,iStream),nb);

		MI_mmse_hard(iSNR) = MI_mmse_hard(iSNR) + getMutualInformationHardInput(tx_bits_all(:,iStream),rx_bits_mmse_hard_all(:,iStream),nb);
		%
		MI_ml(iSNR) = MI_ml(iSNR) + getMutualInformationHardInput(tx_bits_all(:,iStream),rx_bits_ml_all(:,iStream),nb);

		MI_sic_hard(iSNR) = MI_sic_hard(iSNR) + getMutualInformationHardInput(tx_bits_all(:,iStream),rx_bits_sic_hard_all(:,iStream),nb);

		MI_sphere_hard(iSNR) =  MI_sphere_hard(iSNR)+ getMutualInformationHardInput(tx_bits_all(:,iStream),rx_bits_sphere_hard_all(:,iStream),nb);

		MI_full_App(iSNR) = MI_full_App(iSNR) + getMutualInformationSoftInput(tx_bits_all(:,iStream),LLR_full_APP_all(:,iStream),nb);

		MI_bf_bit_MRF(iSNR) = MI_bf_bit_MRF(iSNR) + getMutualInformationHardInput(tx_bits_all(:,iStream), rx_bits_bf_hard_all(:,iStream),nb);

		% MI_bf_bit_GAMP(iSNR) = MI_bf_bit_MRF(iSNR) + getMutualInformationHardInput(tx_bits_all(:,iStream), rx_bits_bf_GAMP_all(:,iStream),nb);
	end

	ber_zf(iSNR) = brrNr_ZF/(nb*Nt*N_ch);
	ber_zf_soft(iSNR) = brrNr_zf_soft/(nb*Nt*N_ch);
	ser_zf(iSNR) = symerrNr_ZF/(Nt*N_ch);

	ber_sphere(iSNR) = brrNr_ZF_Sphere/(nb*Nt*N_ch);

	ber_mmse(iSNR) = brrNr_MMSE/(nb*Nt*N_ch);

	ber_mmse_soft(iSNR) = brrNr_mmse_soft/(nb*Nt*N_ch);
	ser_mmse(iSNR) = symerrNr_MMSE/(Nt*N_ch);

	ber_ml(iSNR) = brrNr_ML/(nb*Nt*N_ch);
	ser_ml(iSNR) = symerrNr_ML/(Nt*N_ch);

	ber_full_App(iSNR) = brrNr_fullAPP/(nb*Nt*N_ch);

	ber_sic(iSNR) = brrNr_SIC/(nb*Nt*N_ch);
	ser_sic(iSNR) = symerrNr_SIC/(Nt*N_ch);

	ber_BF_MRF(iSNR) = brrNr_BF_MRF/(nb*Nt*N_ch);

	if DEBUG
		iSNR
	end

end % for iSNR=1:length(SNRs)

% ------------------------------------------------------------------------------
% debug and output area
if OUTPUT
	% save output to files
	save('data.mat','ber_siso','ber_zf','ber_zf_soft','ber_sphere','ber_mmse','ber_mmse_soft','ber_ml','ber_full_App','ber_sic','ber_BF_MRF', 'MI_zf_bit_hard','MI_zf_soft_all','MI_ml','MI_mmse_hard','MI_mmse_soft_scaled_all','MI_full_App','MI_sic_hard','MI_sphere_hard','MI_bf_bit_MRF');

	grid on;
	title('bit error rate');
	figure(1);
	hold on;
	loglog(SNRs,ber_siso_all,'*-','LineWidth',1.5),
	loglog(SNRs,ber_zf),
	loglog(SNRs,ber_zf_soft,'--'),
	loglog(SNRs,ber_mmse),
	loglog(SNRs,ber_mmse_soft,'--'),
	loglog(SNRs,ber_ml,'.-'),
	loglog(SNRs,ber_sic),
	loglog(SNRs,ber_sphere),
	loglog(SNRs,ber_full_App),
	loglog(SNRs,ber_BF_MRF),
	xlabel('SNR/dB'),ylabel('bit error rate'),
	legend('SISO','ZF','ZF Soft','MMSE hard','MMSE soft','ML','SIC','Sphere detection','fullApp','BF MRF');
	% legend('SISO');

	figure(2);
	grid on;
	hold on;
	semilogx(SNRs, MI_zf_bit_hard),
	semilogx(SNRs, MI_zf_soft_all,'--'),
	semilogx(SNRs, MI_mmse_hard),
	semilogx(SNRs, MI_mmse_soft_scaled_all,'--'),
	semilogx(SNRs, MI_ml,'.-'),
	semilogy(SNRs, MI_sic_hard),
	semilogx(SNRs, MI_full_App),
	semilogy(SNRs, MI_sphere_hard),
	semilogx(SNRs, MI_bf_bit_MRF),
	xlabel('SNR/dB'),ylabel('Mutual Information'),
	legend('ZF hard','ZF soft','MMSE hard','MMSE_soft','ML ','SIC','Full App','Sphere hard','BF MRF');
end
