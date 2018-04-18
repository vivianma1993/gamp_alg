classdef variable_node_GAMP < MotherNode_GAMP
    % VARIABLE_NODE_GAMP 姝ゅ鏄剧ず鏈夊叧姝ょ被鐨勬憳瑕?
    % 姝ゅ鏄剧ず璇︾粏璇存槑
    
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
%             disp("variable_node_GAMP ---")
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
            s.msg_VN = {};
        end
        
        % 杩欒竟濡備綍浼犱袱涓暟杩囧幓
        function [msg_VN, est_x] = factor_fun(s, in_msg, from_id, to_id,default_msg)
            %             in_msg
            %             msg = []; % 杩欎釜浠庡摢閲岃幏寰?
            %             msg_FN = []; % 杩欎釜浠庡摢閲屽緱鍒?
            %             for j = 1:s.n_VN % n_VN is not defined
            %                 for i = 1: s.n_FN % n_FN is not defined
            %                     for k = 1:s.n_FN % n_FN is not defined
            %                         if (j == k)
            %                             continue;
            %                         end
            %                         % 杩欓噷涓嶈兘杩欎箞鍐欙紝娌＄湅鏄庣櫧杩欓噷瑕佸仛浠?箞锛?
            %                         msg(j) = msg(j) + msg_FN(i,j);
            %                     end
            %                 end
            %
            %                 prob_o(j) = exp(msg(j))./(1+exp(msg(j)));
            %                 prob_u(j)= 1- prob_o(j);
            %                 msg_VN(j) = log(prob_o(j)/prob_u(j));
            %                 est_x(j) = tanh(msg_VN(j)./2);
            %             end
            for i = 1 : s.n_FN
%                 for i = 1:length(s.linklist)
%                     if (s.linklist{i}.id == s.link_id)
%                         FN = s.linklist{i};
%                         break;
%                     end
%                 end
%                 in_msg(i) = FN.outboundmsg;
                msg_VN(i) = sum(sum(in_msg));
                prob_o(i) = exp(msg_VN(i))./(1+exp(msg_VN(i)));
                prob_u(i)= 1- prob_o(i);
                msg_VN(i) = log(prob_o(i)/prob_u(i));
                est_x(i) = tanh(msg_VN(i)./2);
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


