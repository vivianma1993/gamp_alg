function[LLRs] = mimo_BP_MRF_detector_v02(y,H_channel,noise_var,const_pointsSymbol,const_pointsLabel,Niter,DampingFactor)
% name:
%		mimo_BP_MRF_detector_v02
%
% describe:
%		Markov random field
%
% Inputs:
%		- y: r,
%		- H_channel: H,
%		- noise_var: noise_var,
%		- const_pointsSymbol: QAM_Symbols,
%		- const_pointsLabel: symlabel,
%		- Niter: Niter,
%		- DampingFactor: 0.2
%
% Outputs:
%		- LLRs:

	y,H_channel,noise_var,const_pointsSymbol,const_pointsLabel,Niter,DampingFactor
	% BP detector based Markov Random Field
	[Nr,Nt] = size(H_channel);
	n_VN = Nt;
	n_FN = n_VN*(n_VN-1)/2;
	FN_edges = nchoosek(1:n_VN,2);
	nb = log2(numel(const_pointsSymbol));
	LLRs = zeros(Nt,nb,Niter+1);

	% --------------------------------------------------------------------------
	% create variable nodes
	for iVN = 1:n_VN
		VN{iVN,:} =  variable_node(iVN);
	end

	% --------------------------------------------------------------------------
	% create Function nodes
	for iFN = 1:n_FN
		nFNIdx = n_VN+iFN;
		FN{iFN,:} =  function_node(nFNIdx);
	end

	% --------------------------------------------------------------------------
	% establish connections between VN and FN
	for iVN = 1:n_VN
		idxes = (FN_edges(:,1) == iVN) | (FN_edges(:,2) == iVN);
		VN{iVN}.setup_link(FN(idxes,:)'); %'
	end

	% --------------------------------------------------------------------------
	% Link connection Function node
	for iFN = 1:n_FN
		FN{iFN}.setup_link(VN(FN_edges(iFN,:))'); %'
	end

	% --------------------------------------------------------------------------
	% setup intial message = 0.5
	Pa = 0.5;
	La = log(Pa/(1-Pa));
	y_MF = H_channel' * y / noise_var; %' z
	R_cov = H_channel' * H_channel / noise_var; %'
	% channel_weight = diag(R_cov);
	% algorithm = 2;
	% initMsg = La + softDemappingSISO(const_pointsLabel,const_pointsSymbol,y_MF,noise_var,algorithm,channel_weight,La);
	nb = log2(numel(const_pointsSymbol));
	ConstLabels_bits = de2bi(const_pointsLabel,nb);
	initMsg = zeros(n_VN,nb);

	% --------------------------------------------------------------------------
	% compute the marginal prob. using joint pdf
	for inb = 1:nb
		l_o = const_pointsSymbol(ConstLabels_bits(:,inb)==1);
		l_u = const_pointsSymbol(ConstLabels_bits(:,inb)==0);
		prob_o = 0;
		prob_u = 0;
		X_o = real(conj(y_MF)*l_o);
		X_u = real(conj(y_MF)*l_u);
		max_o = max(X_o,[],2);
		max_u = max(X_u,[],2);
		for iTest = 1:2^(nb-1)
			prob_o = prob_o + exp(X_o(:,iTest)-max_o);
			prob_u = prob_u + exp(X_u(:,iTest)-max_u);
		end
		initMsg(:,inb) = La + max_o - max_u + log(prob_o./prob_u);
	end

	for  iVN = 1:n_VN
		VN{iVN}.setup_init_msg(initMsg(iVN,:));
	end

	% --------------------------------------------------------------------------
	% setup first block of receive signal
	% for  iFN = 1:n_FN
	%     FN{iFN}.ReceivedSignal = {};
	% end

	LLRs(:,:,1) = initMsg;
	message_old = zeros(n_VN-1,nb,n_VN);
	for iter = 1:Niter
		% ----------------------------------------------------------------------
		% Update Function node
		for Ln = 1:n_FN
			FN{Ln}.update_node({R_cov,noise_var,const_pointsSymbol,const_pointsLabel});
		end

		% ----------------------------------------------------------------------
		% message damping
		if(iter~=1)
			for nVN = 1:n_VN
				for idx =1:length(VN{nVN}.inbound_msg)
					msg_old = message_old(idx,:,nVN);
					msg_new = VN{nVN}.inbound_msg{idx};
					prob_in_msg_old = 1./(1+exp(-msg_old));
					prob_in_msg_new = 1./(1+exp(-msg_new));
					for inb = 1:nb
						if(~(prob_in_msg_old(inb)==1 && prob_in_msg_new(inb)==1))
							 VN{nVN}.inbound_msg{idx}(inb) = log((DampingFactor*prob_in_msg_old(inb) + (1-DampingFactor)*prob_in_msg_new(inb)))-log(DampingFactor*(1-prob_in_msg_old(inb)) + (1-DampingFactor)*(1-prob_in_msg_new(inb)));
						end
					end
					%max_o = max([msg_old + log(DampingFactor); msg_new + log(1-DampingFactor); msg_old + msg_new]);
					%max_u = max([msg_old + log(DampingFactor); msg_new + log(1-DampingFactor); zeros(1,nb)]);
					%VN{nVN}.inbound_msg{idx} = max_o - max_u + log(exp(msg_old + log(DampingFactor)-max_o) + exp(msg_new + log(1-DampingFactor)-max_o) + exp(msg_old + msg_new -max_o)) - log(exp(msg_old + log(DampingFactor)-max_u) + exp(msg_new + log(1-DampingFactor)-max_u) + exp(zeros(1,nb) -max_u));
				end
			end

		end

		% ----------------------------------------------------------------------
		% update variable node
		for Ln = 1:n_VN
			VN{Ln}.update_node({y_MF,noise_var,const_pointsSymbol,const_pointsLabel});
		end

		% ----------------------------------------------------------------------
		% SaveLLR for each iteration
		for nVN = 1:n_VN
			message_old(:,:,nVN) = cell2mat(VN{nVN}.inbound_msg);
			LLRs(nVN,:,iter+1) = initMsg(nVN,:)+ sum(cell2mat(VN{nVN}.inbound_msg));
		end
	end

end
