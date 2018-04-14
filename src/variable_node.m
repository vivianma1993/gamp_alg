classdef variable_node < MotherNode
    
    properties % class当中的公共属性
        init_msg
        ReceivedSignal
    end
    
    methods
        function s = variable_node(id)
            s = s@MotherNode(id);
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
            reset@MotherNode(s);
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
%                 if (~isempty(default_msg))
%                     y_MF = default_msg{1};
%                     ConstSymbols = default_msg{3};
%                     ConstLabels = default_msg{4};
%                     prob_in_msg = 1- 1./(1+exp(msg)); % prob of being 1 in the corresponding position
%                     nb = log2(numel(ConstSymbols));
%                     ConstLabels_bits = de2bi(ConstLabels,nb);
%                     prob_syms = ones(2^nb,1);
%                     
%                     for inb = 1:nb
%                         prob_syms(ConstLabels_bits(:,inb)==1)=prob_syms(ConstLabels_bits(:,inb)==1)*prob_in_msg(inb);
%                         prob_syms(ConstLabels_bits(:,inb)==0)=prob_syms(ConstLabels_bits(:,inb)==0)*(1-prob_in_msg(inb));
%                     end
%                     
%                     % compute the marginal prob. using joint pdf
%                     for inb= 1:nb
%                         l_o = ConstSymbols(ConstLabels_bits(:,inb)==1);
%                         l_u = ConstSymbols(ConstLabels_bits(:,inb)==0);
%                         prob_o= 0;
%                         prob_u =0;
%                         for iTest = 1:2^(nb-1)
%                             prob_o = prob_o + exp(real(conj(y_MF(s.id))*l_o(iTest)))*prob_syms(iTest);
%                             prob_u = prob_u + exp(real(conj(y_MF(s.id))*l_u(iTest)))*prob_syms(iTest);
%                         end
%                         msg(1,inb) = log(prob_o/prob_u);
%                     end
%                 end
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
