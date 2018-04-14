classdef function_node < MotherNode

	properties
		init_msg
		ReceivedSignal
	end

	methods
		function s = function_node(id)
			s = s@MotherNode(id);
			s.ReceivedSignal = {};
		end

		function LLRValue = factor_fun(s, in_msg, from_id, to_id,default_msg)
			R_cov = default_msg{1};
			NoiseVar = default_msg{2};
			ConstSymbols = default_msg{3};
			ConstLabels = default_msg{4};
			%in_msg{1,1} = zeros(1,2);
			%prob_in_msg = 1- 1./(1+exp(in_msg{1,1})); % prob of being 1 in the corresponding position
			nb = log2(numel(ConstSymbols));
			ConstLabels_bits = de2bi(ConstLabels,nb);
			prob_syms = zeros(2^nb,1);
			LLRValue = zeros(1,nb);

			for inb = 1:nb
				prob_syms(ConstLabels_bits(:,inb)==1)=prob_syms(ConstLabels_bits(:,inb)==1) + in_msg{1,1}(inb);
				prob_syms(ConstLabels_bits(:,inb)==0)=prob_syms(ConstLabels_bits(:,inb)==0) + 0;
			end

			% compute the marginal prob. using joint pdf
			for inb= 1:nb
				l_o = ConstSymbols(ConstLabels_bits(:,inb)==1);
				l_u = ConstSymbols(ConstLabels_bits(:,inb)==0);
				X_o = real(-R_cov(to_id,from_id)*(l_o')*ConstSymbols)+ repmat(prob_syms',2^(nb-1),1);
				X_u = real(-R_cov(to_id,from_id)*(l_u')*ConstSymbols)+ repmat(prob_syms',2^(nb-1),1);
				max_o = max(max(X_o));
				max_u = max(max(X_u));
				prob_o= 0;
				prob_u =0;
				for iTest = 1:2^(nb-1)
					for iTest2 = 1:2^(nb)
						prob_o = prob_o + exp(X_o(iTest,iTest2)-max_o);
						prob_u = prob_u + exp(X_u(iTest,iTest2)-max_u);
					end
				end
				LLRValue(1,inb) = max_o - max_u +  log(prob_o/prob_u);
			end
		end

	end % methods

end
