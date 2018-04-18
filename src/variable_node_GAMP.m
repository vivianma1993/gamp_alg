classdef variable_node_GAMP < MotherNode_GAMP
    % VARIABLE_NODE_GAMP 此处显示有关此类的摘要
    % 此处显示详细说明

    properties
        msg_VN
        ReceivedSignal
        n_VN
        n_FN
        DEBUG
        TAG = 'variable_node_GAMP'
    end

    methods
        function s = variable_node_GAMP(id, vn, fn)
            s = s@MotherNode_GAMP(id);
            s.msg_VN = {};
            s.ReceivedSignal = {};
            s.DEBUG = true;
            s.n_VN = vn;
            s.n_FN = fn;
        end

        function setup_init_msg(s, msg_VN)
            if size(msg_VN,1) ~= 1
                sprintf('Node %0.0f ERROR from setup_init_msg: message must be row vector !!!\n', s.id)

            else
                s.init_msg{1} = msg_VN;
                construct message
                snd_msg = s.factor_fun(s.init_msg,{},{},{});
                snd_msg = msg_VN;

                % send message
                for linkNo = 1:numel(s.linklist)
                    s.linklist{linkNo}.rx_msg(s, snd_msg);
                    % update outbound record
                    s.to_node{linkNo} = s.linklist{linkNo};
                    s.outbound_msg{linkNo} = snd_msg;
                    s.to_id(linkNo) = [s.linklist{linkNo}.id];
                end
            end
        end

        function reset(s)
            reset@MotherNode_GAMP(s);
            s.init_msg = {};
        end

        % 这边如何传两个数过去
        function [msg_VN, est_x] = factor_fun(s, in_msg, from_id, to_id,default_msg)
            in_msg
            msg = []; % 这个从哪里获得
            msg_FN = []; % 这个从哪里得到
            for j = 1:s.n_VN % n_VN is not defined
                for i = 1: s.n_FN % n_FN is not defined
                    for k = 1:s.n_FN % n_FN is not defined
                        if (j == k)
                            continue;
                        end
                        % 这里不能这么写，没看明白这里要做什么？
                        msg(j) = msg(j) + msg_FN(i,j);
                    end
                end

                prob_o(j) = exp(msg(j))./(1+exp(msg(j)));
                prob_u(j)= 1- prob_o(j);
                msg_VN(j) = log(prob_o(j)/prob_u(j));
                est_x(j) = tanh(msg_VN(j)./2);
            end
        end

    end

    methods (Access = protected)
        function rx_msg(s, from_node, msg_FN)
            from_nodeID = from_node.id;
            from_nodeIndx = find(s.link_id == from_nodeID);
            fprintf('%s: from_nodeID: %d, from_nodeIndx: %d \n', s.TAG, from_nodeID, from_nodeIndx);
            if 0 == length(from_nodeIndx)
                disp([s.TAG ': length of from_nodeIndx is equal to 0']);
            else
                disp([s.TAG ': add to s.from_node'])
                s.from_node{from_nodeIndx}      = from_node;
                s.inbound_msg{from_nodeIndx}    = msg_FN;
                s.from_id(from_nodeIndx)        = [from_node.id];
            end
        end

    end
end


