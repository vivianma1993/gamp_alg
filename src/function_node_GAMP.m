classdef function_node_GAMP < MotherNode_GAMP
    % FUNCTION_NODE_GAMP 此处显示有关此类的摘要
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
        function msg_FN = factor_fun(s, in_msg, from_id, to_id,default_msg)
            y = default_msg{1};
            H_channel = default_msg{2};
            mean = default_msg{3};
            noise_var = default_msg{4};
            var_GAMP = default_msg{5};
            est_x = default_msg{6};
            LLRValue = zeros(1,1);
            [Nr,Nt] = size(H_channel);
            n_VN = Nt;
            n_FN = Nr;
            var_x = ones(n_VN,1);
            est_x = zeros (n_VN,1);
            msg_FN = zeros(n_FN,n_VN);
            mean = zeros(n_FN,n_VN);
            var_x = ones (n_VN,1);
            var_GAMP = ones (n_FN,n_VN);

            for i = 1 : n_FN
                for j = 1: n_VN
                    for k = 1:n_FN
                        if (k == i)
                            continue;
                        end
                        mean (i,j) = mean (i,j) + H_channel(k,j) * est_x(k);
                        var_x(i) = 1 - est_x(k).^2;
                        var_GAMP(i,j) = abs(H_channel(k,j)'*H_channel(k,j)) * var_x(k);
                    end
                end
                var_GAMP (i,j)= var_GAMP(i,j) + noise_var ;
                for s = 1:n_VN
                    msg_FN(i,s) = 4 ./ var_GAMP(i,s) * (real(conj( H_channel(i,s))*(y(i)-mean(i,s))));
                end
            end
        end

    end
end
