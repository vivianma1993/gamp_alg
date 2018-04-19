classdef variable_node_GAMP< MotherNode_GAMP
    %VARIABLE_NODE_GAMP
    
    properties
        msg_VN
        ReceivedSignal
        n_VN
        n_FN
        DEBUG
        TAG = 'variable_node_GAMP'
    end
    
    methods
        %         function s = variable_node_GAMP(id, vn, fn)
        % %             disp("variable_node_GAMP ---")
        %             s = s@MotherNode_GAMP(id);
        %             s.msg_VN = {};
        %             s.ReceivedSignal = {};
        %             s.DEBUG = true;
        %             s.n_VN = vn;
        %             s.n_FN = fn;
        %         end
        
        function s = variable_node_GAMP(id)
            s = s@MotherNode_GAMP(id);
            s.msg_VN = {};
            s.ReceivedSignal = {};
            s.DEBUG = true;
        end
        
        
        function reset(s)
            reset@MotherNode_GAMP(s);
            s.msg_VN = {};
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
        
        
        function msg_VN = factor_fun(s, in_msg, from_id, to_id,default_msg) %% 这边如何传两个数过去
            %    VN{Ln}.update_node({y,H_channel,mean,var_GAMP,msg_FN});
            y = default_msg{1};
            H_channel = default_msg{2};
            mean = default_msg{3};
            var_GAMP = default_msg{4};
            msg_FN = default_msg{5};
            [Nr,Nt] = size(H_channel);
            n_VN = Nt;
            n_FN = Nr;
            msg_VN = zeros (n_VN,1);

            for i = 1 : n_VN
                %msg_FN(i,:) = zeros(1,n_VN);%% 什么好方法
                tmp = msg_FN(i,:);
                msg_FN(i,:) = zeros(1, length(msg_FN(i,:)));
                msg_VN(i) = sum(sum(msg_FN(:,:)));
                prob_o(i) = exp(msg_VN(i))./(1+exp(msg_VN(i)));
                prob_u(i)= 1- prob_o(i);
                msg_VN(i) = log(prob_o(i)/prob_u(i));
                est_x(i) = tanh(msg_VN(i)./2);
                msg_FN(i,:) = tmp;
            end
            msg_VN = est_x;
                          
%                              for i = 1 : n_VN
%                                 msg = sum(msg_FN.');
%                                 msg(i) = [];
%                                 msg_VN = sum(msg(:,:));
%                                 prob_o(i) = exp(msg_VN(i))./(1+exp(msg_VN(i)));
%                                 prob_u(i)= 1- prob_o(i);
%                                 msg_VN(i) = log(prob_o(i)/prob_u(i));
%                                 est_x(i) = tanh(msg_VN(i)./2);
%                           end
                              
                              



        end
    end
 
        
        
        methods (Access = protected)
            function rx_msg(s, from_node, msg_FN)
                from_nodeID = from_node.id;
                from_nodeIndx = find(s.link_id == from_nodeID);
%                 fprintf('%s: from_nodeID: %d, from_nodeIndx: %d \n', s.TAG, from_nodeID, from_nodeIndx);
                if 0 ~= length(from_nodeIndx)
%                     disp([s.TAG ': length of from_nodeIndx is equal to 0']);
%                 else
%                     disp([s.TAG ': add to s.from_node'])
                    s.from_node{from_nodeIndx}      = from_node;
                    s.inbound_msg{from_nodeIndx}    = msg_FN;
                    s.from_id(from_nodeIndx)        = [from_node.id];
                end
            end
        end
    end


