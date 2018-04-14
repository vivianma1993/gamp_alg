function[LLRs] = mimo_BP_MRF_detector(y,H_channel,noise_var,const_pointsSymbol,const_pointsLabel,Niter)
% BP detector based Markov Random Field
[Nr,Nt] = size(H_channel);
n_VN = Nt;
n_FN = n_VN*(n_VN-1)/2;
FN_edges = nchoosek(1:n_VN,2);
nb = log2(numel(const_pointsSymbol));
LLRs = zeros(Nt,nb,Niter+1);

% create variable nodes
for iVN = 1:n_VN
    VN{iVN,:} =  variable_node_GAMP(iVN);
end

% create Function nodes
for iFN = 1:n_FN
    nFNIdx = n_VN+iFN;
    FN{iFN,:} =  function_node_GAMP(nFNIdx);
end

% establish connections between VN and FN
for iVN = 1:n_VN
    idxes = (FN_edges(:,1) == iVN) | (FN_edges(:,2) == iVN);
    VN{iVN}.setup_link(FN(idxes,:)');
end

% Link connection Function node 建立FN与VN之间的对应关系
for iFN = 1:n_FN
    FN{iFN}.setup_link(VN(FN_edges(iFN,:))');
end

%setup intial message =0.5
Pa = 0.5;
%La = log(Pa/(1-Pa));
y_MF = H_channel'*y/noise_var; %等效y
R_cov = H_channel'*H_channel/noise_var;
% channel_weight = diag(R_cov);
% algorithm = 2;
% initMsg = La + softDemappingSISO(const_pointsLabel,const_pointsSymbol,y_MF,noise_var,algorithm,channel_weight,La);
nb = log2(numel(const_pointsSymbol));
ConstLabels_bits = de2bi(const_pointsLabel,nb);
initMsg = zeros(n_VN,nb);

% compute the marginal prob. using joint pdf
for inb= 1:nb
    l_o = const_pointsSymbol(ConstLabels_bits(:,inb)==1);
    l_u = const_pointsSymbol(ConstLabels_bits(:,inb)==0);
    prob_o= 0;
    prob_u =0;
    for iTest = 1:2^(nb-1)
        prob_o = prob_o + exp(real(conj(y_MF)*l_o(iTest)))*Pa;
        prob_u = prob_u + exp(real(conj(y_MF)*l_u(iTest)))*(1-Pa);
    end
    initMsg(:,inb) = log(prob_o./prob_u);
end

for  iVN = 1:n_VN
    VN{iVN}.setup_init_msg(initMsg(iVN,:));
end


% %setup first block of receive signal
% for  iFN = 1:n_FN
%     FN{iFN}.ReceivedSignal = {};
% end


LLRs(:,:,1) = initMsg;
for iter = 1:Niter
    %Update Function node
    for Ln = 1:n_FN
        FN{Ln}.update_node({R_cov,noise_var,const_pointsSymbol,const_pointsLabel});
    end
    
    %update variable node
    for Ln = 1:n_VN
        VN{Ln}.update_node({y_MF,noise_var,const_pointsSymbol,const_pointsLabel});
    end
    
    %SaveLLR for each iteration
    
    for nVN = 1:n_VN
        LLRs(nVN,:,iter+1) = initMsg(nVN,:)+ sum(cell2mat(VN{nVN}.inbound_msg));
    end
end

end