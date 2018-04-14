function[DecodedSignal, FGExportData] = FactorGraphDecode_2(x,h,varience,const_pointsSymbol,const_pointsLabel,data_user)
data_user = data_user';
[n_FN,n_VN] = size(h);
% 
for nVN = 1:n_VN
    VN{nVN,:} =  variable_node(nVN);
end

%Function nodes
for nFN = 1:n_FN
    nFNIdx = n_VN+nFN;
    FN{nFN,:} =  function_node(nFNIdx);
end

%Link connection Variable Node
for Ln = 1:size(h,2)
    pos_ones = find(h(:,Ln));
    VN{Ln}.setup_link(FN(pos_ones)');
end

%Link connection Function node
for Ln = 1:size(h,1)
    pos_ones = find(h(Ln,:));
    FN{Ln}.setup_link(VN(pos_ones)');
end

%setup intial message =0.5
Pa = 0.5;
La = log(Pa/(1-Pa));
for  Ln = 1:size(h,2)
    VN{Ln}.setup_init_msg(La);
    
end


%setup first block of receive signal
for  Ln = 1:size(h,1)
    FN{Ln}.ReceivedSignal = x(Ln,:);
end



for i = 1:10
    LLRValueInEachIteration = zeros(size(data_user));
    %Update Function node
    for Ln = 1:n_FN
        FN{Ln}.update_node({h,varience,const_pointsSymbol,const_pointsLabel});
    end
    %update variable node
    for Ln = 1:n_VN
        VN{Ln}.update_node({h,varience,const_pointsSymbol,const_pointsLabel});
    end
    
    %SaveLLR for each iteration
    
    for nVN = 1:n_VN
        LLRValueInEachIteration(nVN,:) = sum(cell2mat(VN{nVN}.inbound_msg));
        LLRValues = LLRValueInEachIteration(nVN,:);
        UserTxData = data_user(nVN,:);
        FGExportData.Iteration(i).User(nVN).MutualInfo = measure_mutual_information_histogram(-LLRValues(:),UserTxData(:));
        
        LLRValues(LLRValues >=0) = 0;
        LLRValues(LLRValues <0) = 1;        
        BER = UserTxData-LLRValues;
        BER(abs(BER)>10^-4) = 1;        
        FGExportData.Iteration(i).User(nVN).BER= nnz(BER)/numel(UserTxData);
    end 
    FGExportData.Iteration(i).count = i;
    Iteration = i
end

for nVN = 1:n_VN
    LValue(nVN,:) =  sum(cell2mat(VN{nVN}.inbound_msg));
end

LValue(LValue>=0) = 0;
LValue(LValue<0) = 1;
LValue = LValue+1;
for UserNo = 1:size(LValue,1)
    UserConstSymb = (const_pointsSymbol(:,UserNo));
    %HardDemapping
    DecodedSignal(UserNo,:) = UserConstSymb(LValue(UserNo,:));
end



% DecodedSignal = LValue;
%__________________________________________________________________________
    
end