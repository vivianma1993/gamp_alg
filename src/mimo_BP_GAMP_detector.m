function[LLRs] = mimo_BP_GAMP_detector(y,H_channel,noise_var,const_pointsSymbol,const_pointsLabel,Niter)
% BP detector based Markov Random Field
[Nr,Nt] = size(H_channel);
n_VN = Nt;
n_FN = Nr;
FN_edges = nchoosek(1:n_VN,n_VN);
nb = log2(numel(const_pointsSymbol));
% LLRs = zeros(Nt,nb,Niter+1);% Number of iteration
 LLRs = zeros(Nt,Niter+1);

% create variable nodes
for iVN = 1:n_VN
    VN{iVN,:} =  variable_node_GAMP(iVN);
end

% create Function nodes
for iFN = 1:n_FN
    nFNIdx = n_VN+iFN; % 给与指定VN相连的FN编号
    FN{iFN,:} =  function_node_GAMP(nFNIdx);
end

% establish connections between VN and FN
for iVN = 1:n_VN
    VN{iVN}.setup_link(FN(iVN,:)');
end

% Link connection Function node
for iFN = 1:n_FN
    FN{iFN}.setup_link(VN(FN_edges(1,:))');
end

nb = log2(numel(const_pointsSymbol));
ConstLabels_bits = de2bi(const_pointsLabel,nb);
initMsg = zeros(n_VN,1);
mean_0 = zeros(n_VN,1);
var_GAMP_0 = ones(n_VN,1); %?
s_init = 0;
% pro_x = ones(n_VN,1);

% compute the marginal prob. using joint pdf
for i = 1 : n_FN
    for j = 1: n_VN
        if (i == j)
            continue ;
        end
        mean_0(i,:) = mean_0(i,:) + H_channel(i,j) * mean_0(i,:);
        var_GAMP_0(i,:) = abs(H_channel(i,j)'*H_channel(i,j)) * var_GAMP_0(i,:);
    end
    var_GAMP = var_GAMP_0 + noise_var ;
    initMsg(i,:) = 4 / var_GAMP(i,:) * (real(conj( H_channel(i,i))*(y(i)-mean_0(i,:))));
end
% for i = 1 : n_FN
%     for j = 1: n_VN
%         if (i == j)
%             continue ;           
%         end
%         s_init = s_init + initMsg(j,:);
%     end

for  iVN = 1:n_VN
    VN{iVN}.setup_init_msg(initMsg(iVN,:));
end


% %setup first block of receive signal
% for  iFN = 1:n_FN
%     FN{iFN}.ReceivedSignal = {};
% end


LLRs(:,1) = initMsg;
for iter = 1:Niter
    %Update Function node
    for Ln = 1:n_FN
        FN{Ln}.update_node({y,H_channel,mean_0,var_GAMP});
    end
   
    %update variable node
    for Ln = 1:n_VN
        VN{Ln}.update_node({y,H_channel,mean_0,var_GAMP});
    end
    
    %SaveLLR for each iteration
    
    for nVN = 1:n_VN
        LLRs(nVN,:,iter+1) = initMsg(nVN,:)+ sum(cell2mat(VN{nVN}.inbound_msg));
    end
end

end