classdef function_node_GAMP < MotherNode_GAMP
    %FUNCTION_NODE_GAMP 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        init_msg
        ReceivedSignal
    end
    
    methods
        function s = function_node_GAMP(id)
            s = s@MotherNode_GAMP(id);
            s.ReceivedSignal = {};
        end
        function LLRValue = factor_fun(s, in_msg, from_id, to_id,default_msg)
            y = default_msg{1};
            H_channel = default_msg{2};
            mean_0 = default_msg{3};
            var_GAMP = default_msg{4};
            ConstSymbols = default_msg{5};
            ConstLabels = default_msg{6};
            nb = log2(numel(ConstSymbols));
            ConstLabels_bits = de2bi(ConstLabels,nb);
%           prob_syms = zeros(2^nb,1);
            LLRValue = zeros(1,nb);
            
            
            for i = 1 : n_FN
                for j = 1: n_VN
                    mean_0(i,:) = mean_0(i,:) + H_channel(i,j) * mean_0(i,:);
                    var_GAMP_0(i,:) = abs(H_channel(i,j)'*H_channel(i,j)) * var_GAMP_0(i,:);
                end
                var_GAMP = var_GAMP_0 + noise_var ;
                LLRValue(i,:) = 4 / var_GAMP(i,:) * (real(conj( H_channel(i,j))*(y(i)-mean_0(i,:))));
                
            end
        end
    end
end


