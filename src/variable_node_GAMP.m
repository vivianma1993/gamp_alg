classdef variable_node_GAMP< MotherNode_GAMP
    %VARIABLE_NODE_GAMP 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
         init_msg
        ReceivedSignal
    end
    
    methods
        function s = variable_node_GAMP(id)
            s = s@MotherNode_GAMP(id);
            s.init_msg = {};
            s.ReceivedSignal = {};   
        end
        
      function setup_init_msg(s, msg)
            if size(msg,1) ~= 1
                sprintf('Node %0.0f ERROR from setup_init_msg: message must be row vector !!!\n',s.id)
            else
                s.init_msg{1} = msg;
                % construct message
                %snd_msg = s.factor_fun(s.init_msg,{},{},{});
                snd_msg = msg;
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
        
        %        function setup_link(s,linklist)
        %            if size(linklist,1)*size(linklist,2) ~= 1
        %                sprintf('Node %0.0f ERROR from setup_link: number of links of an evident node must be 1 !!!\n',s.id)
        %            else
        %                setup_link@MotherNode(s,linklist);
        %            end
        %        end
        
        function msg = factor_fun(s, in_msg, from_id, to_id,default_msg)
            if(~isempty(in_msg))
                msg = s.init_msg{1};
                for priorN = 1:size(in_msg,1)
                   % msg = msg+in_msg{priorN};
                    msg = msg +in_msg{priorN};
                end
                                
                for i = nVN
                    prob_o = 0;
                    prob_u = 0;
                   for j= nFN
                        prob_o = prob_o + exp(initMsg(j,:) );
                        prob_u = 1-pro_o;
            
                   end 
                end
                LLRValue = log(prob_o/prob_u);
                est_x = tanh(2./LLRValue);
            else
                %nb = log2(numel(default_msg{3}));
                %msg =zeros(1,nb);
                msg=s.init_msg{1};
            end
            
        end
    end
    
    methods (Access = protected)
        function rx_msg(s, from_node, msg)
            from_nodeID = from_node.id;
            from_nodeIndx = find(s.link_id == from_nodeID);
            s.from_node{from_nodeIndx} = from_node;
            s.inbound_msg{from_nodeIndx} = msg;
            s.from_id(from_nodeIndx) = [from_node.id];
        end
    end
end


