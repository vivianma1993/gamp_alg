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

        function setup_init_msg(s, msg) % 把VN初始化的数值传进来 此处S为VN
            if size(msg,1) ~= 1 % 保证初始化数据是个行向量
                sprintf('Node %0.0f ERROR from setup_init_msg: message must be row vector !!!\n',s.id)
            else
                s.init_msg{1} = msg;
                snd_msg = msg; % 初始化的数据作为要发送的数据
                % send message
                for linkNo = 1:numel(s.linklist) % s.linklist 是与同一个VN连接的FN
                    s.linklist{linkNo}.rx_msg(s, snd_msg); %把VN的初始值传给每个和他连接的FN
                    % update outbound record
                    s.to_node{linkNo} = s.linklist{linkNo};
                    s.outbound_msg{linkNo} = snd_msg; % VN outboundmsg
                    s.to_id(linkNo) = [s.linklist{linkNo}.id];
                end
            end
        end

        function reset(s)
            reset@MotherNode(s);
            s.init_msg = {};
        end

        function msg = factor_fun(s, in_msg, from_id, to_id,default_msg)
            if(~isempty(in_msg))
                msg = s.init_msg{1};
                for priorN = 1:size(in_msg,1)
                   % msg = msg+in_msg{priorN};
                    msg = msg +in_msg{priorN}; %最终的 msg=msg+initial_msg
                end

            else
                msg=s.init_msg{1};
            end

        end
    end

    methods (Access = protected)
        % 确认哪些代码可以调用此方法 protected - 从类或子类的方法进行访问
        function rx_msg(s, from_node, msg) % S是VN
            from_nodeID = from_node.id;
            from_nodeIndx = find(s.link_id == from_nodeID);
            % s.link_id from_node.id 是与VN 相连的FN 的id
            % 返回一个包含数组 X 中每个非零元素的线性索引的矢量
            s.from_node{from_nodeIndx} = from_node; %from_node 是一个function node
            s.inbound_msg{from_nodeIndx} = msg;
            s.from_id(from_nodeIndx) = [from_node.id];
        end
    end
end
