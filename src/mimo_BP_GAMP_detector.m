function[LLRs] = mimo_BP_GAMP_detector(y,H_channel,noise_var,const_pointsSymbol,const_pointsLabel,Niter)
% BP detector based Markov Random Field
%   1. miu_z_i_k <- sum(h_i_j * E_x_i)
%   2. o_z_i_k <- sum(|h_i_j|^2 Var(x_j)) + o^2
%   3. v_i_k <- 4/o_z_i_k^2 * R(h_i_k* (y_i - miu_z_i_k))
%   4. p_i_k_plus <- exp( sum(l!=i v_l_k) ) / ( 1+ exp( sum(l!=i v_l_k) ))
%
% Inputs:
%   x: s
%   y: r
%   H_channel: H_channel
%   E: 0
%   ni: n
%   o: noise_var
% -----------------------------
%   y,
%   H_channel,
%   noise_var,
%   const_pointsSymbol,
%   const_pointsLabel,
%   Niter
%
% Outputs:
%   xxxxx
%

    [Nr,Nt] = size(H_channel);
    n_VN = Nt;
    n_FN = Nr;
    FN_edges = nchoosek(1:n_VN,n_VN);
    nb = log2(numel(const_pointsSymbol));
    % LLRs = zeros(Nt,nb,Niter+1);% Number of iteration
    LLRs = zeros(Nt,Nr,Niter+1);
    est_x = zeros (n_VN,1);
    msg = zeros(n_VN,1);
    msg_VN = zeros(n_VN,1);
    prob_o = zeros(n_VN,1);
    prob_u = zeros(n_VN,1);

    % create variable nodes
    for iVN = 1:n_VN
        VN{iVN,:} =  variable_node_GAMP(iVN, n_VN, n_FN);
    end

    % create Function nodes
    for iFN = 1:n_FN
        nFNIdx = n_VN+iFN; % 给与指定VN相连的FN编号
        FN{iFN,:} =  function_node_GAMP(nFNIdx);
    end

    % establish connections between VN and FN
    for iVN = 1:n_VN
        VN{iVN}.setup_link(FN(iVN,:)');
    end

    % Link connection Function node
    for iFN = 1:n_FN
        FN{iFN}.setup_link(VN(FN_edges(1,:))');
    end

    nb = log2(numel(const_pointsSymbol));
    ConstLabels_bits = de2bi(const_pointsLabel,nb);
    msg_FN = zeros(n_FN,n_VN);
    mean = zeros(n_FN,n_VN);
    var_x = ones (n_VN,1);
    var_GAMP = ones (n_FN,n_VN);
    est_x = zeros (n_VN,1);

    % compute the marginal prob. using joint pdf
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
            var_GAMP (i,j)= var_GAMP(i,j) + noise_var ;
            for iVn = 1:n_VN
                msg_FN(i,iVN) = 4 ./ var_GAMP(i,iVN) * (real(conj( H_channel(i,iVN))*(y(i)-mean(i,iVN))));
            end
        end
    end

    for j = 1:n_VN
        for i = 1: n_FN
            for k = 1:n_FN
                if (j == k)
                    continue;
                end
            msg(j) = msg(j) + msg_FN(i,j);
            end
        end
            prob_o(j) = exp(msg(j))./(1+exp(msg(j)));
            prob_u(j)= 1- prob_o(j);
            msg_VN(j) = log(prob_o(j)/prob_u(j));
            est_x(j) = tanh(msg_VN(j)./2);
    end


    % for  iVN = 1:n_VN
    %     VN{iVN}.setup_init_msg(initMsg(iVN,:));
    % end

    LLRs(:,:,1) = msg_FN;
    for iter = 1:Niter
        %Update Function node
        for Ln = 1:n_FN
            FN{Ln}.update_node({y,H_channel,mean,noise_var,var_GAMP,est_x}); % 输入多少个值就有多少default msg
        end

        %update variable node
        for Ln = 1:n_VN
            VN{Ln}.update_node({y,H_channel,mean,var_GAMP,msg_FN});
        end

        %SaveLLR for each iteration

        for nVN = 1:n_VN
            LLRs(nVN,:,iter+1) = initMsg(nVN,:)+ sum(cell2mat(VN{nVN}.inbound_msg));
        end
    end

end
