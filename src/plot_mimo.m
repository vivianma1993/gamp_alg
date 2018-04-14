clear all;
SNRs = -10:2.5:10;
grid on;
title('symbol error rate');
figure(1);
hold on;
% loglog(SNRs,ber_siso),
% loglog(SNRs,ber_zf,'b'),
% loglog(SNRs,ber_mmse,'r'),
loglog(SNRs,ber_ml,'k'),
% loglog(SNRs,ber_sic,'g'),
loglog(SNRs,ber_mmse_soft,'--'),
% loglog(SNRs,ber_zf_soft,'--'),
% loglog(SNRs,ber_sphere),
% loglog(SNRs,ber_full_App),
loglog(SNRs,ber_BF_MRF),
xlabel('SNR/dB'),ylabel('bit error rate'),
legend('MMSE','ZF','ML','SIC','SISO','MMSE Soft','ZF Soft','Sphere detection','fullApp','BF MRF');

figure(2);
hold on;
% semilogx(SNRs, MI_mmse_hard),
% semilogx(SNRs, MI_zf_bit_hard),
semilogx(SNRs, MI_ml),
% semilogy(SNRs, MI_sic_hard),
% semilogx(SNRs, MI_full_App),
% semilogy(SNRs, MI_sphere_hard),
semilogx(SNRs, MI_bf_bit_hard),
% semilogx(SNRs, MI_zf_soft_all,'--'),
%semilogx(SNRs, MI_mmse_soft_scaled_all,'--'),
xlabel('SNR/dB'),ylabel('Mutual Information'),
legend('MMSE hard','ZF hard','ML hard','SIC hard','Sphere hard','ZF soft','MMSE soft','fullApp','MI bf bit hard');
