function [s_Hat, s_Dist,s_DistVec] = ML_rx(y,H,QAM_Symbols)

% input: y received symbols in Nr * Nsym form 
%        H the channel matrix  Nr * Nt*Nsym form 
%        QAM_Symbols possible QAM symbols
% output: estimated symbol vectors 

[Nr, Nt, Nsym] = size(H); % number of tx antennas
M_lattice = length(QAM_Symbols);



M_b_qam = ceil(log2(M_lattice));
Nr_Hypothesis = (M_lattice)^Nt;

% construct the set of all Hypotheses 
bi_index = de2bi(0:(Nr_Hypothesis-1));
qam_index = zeros(Nr_Hypothesis,Nt);
for iNt  =1:Nt
qam_index(:,iNt) = 1+bi2de(bi_index(:,(1+(iNt-1)*M_b_qam):(iNt*M_b_qam)));
end
if(Nt>1)
S_hyp_set = QAM_Symbols(qam_index).';
else
    S_hyp_set = QAM_Symbols(qam_index);
end

s_Hat = zeros(Nt,Nsym);
s_Dist = zeros(1,Nsym);
s_DistVec = zeros(Nr,Nsym);
for idxSym = 1:Nsym
    %ML detector 
    r_mat =repmat(y(:,idxSym),1,Nr_Hypothesis);
    err_vec = (r_mat-H(:,:,idxSym)*S_hyp_set);
    
    if(Nr>1)
      euc_dist = sum(abs(err_vec).^2);
    else
      euc_dist =(abs(err_vec).^2);
    end
    
    [min_err,min_idx]=min(euc_dist);
    s_Hat(:,idxSym)= S_hyp_set(:,min_idx);
    s_Dist(:,idxSym) = min_err;
    s_DistVec(:,idxSym) = abs(err_vec(:,min_idx)).^2;
end



