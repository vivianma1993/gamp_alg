function [ s_est_hard_SIC,s_est_soft_SIC,b_est_hard_SIC] = mimoSICrx(channelH,rxSig,noise_var,ZForMMSE,ordering,QAM_Symbols,symLabeling)

    [Nr, Nt] = size(channelH); % get the dimension of MIMO systems

    if(ZForMMSE) % 0 ZF, 1 MMSE
        mmse_enable = 1;
    else
        mmse_enable = 0;
    end

    w_rx =  cell(Nt,1);
    w_rx{1} = (channelH'*channelH+noise_var*eye(Nt)*mmse_enable)\channelH';
    %sort the channel according to their power level
    if(ordering)
       streamIdx = get_first_SIC_stream_todetect(w_rx{1},channelH,noise_var);
    else
       streamIdx = 1;
    end

    w1_1 = w_rx{1}(streamIdx,:);
    y_1 = w1_1*rxSig;
    s_est_soft = y_1;

    idxes_tobeDetected = 1:Nt;
    [b_est,s_est] = hardDemapping(symLabeling,QAM_Symbols,y_1);
    detected_stream_index = streamIdx;
    idxes_tobeDetected(streamIdx)=[];

    for iStream =2:Nt
        H_prev = channelH(:,detected_stream_index);
        H_current = channelH(:,idxes_tobeDetected);
        w_rx{iStream} = (H_current'*H_current+noise_var*eye(Nt-iStream+1)*mmse_enable)\H_current';
        if(ordering)
           streamIdx = get_first_SIC_stream_todetect(w_rx{iStream},H_current,noise_var);
        else
            streamIdx = 1;
        end

        %w1_1 = w_rx{1}(streamIdx,:);
        %y_1 = w1_1*rxSig;
        w_rx_sic = w_rx{iStream}(streamIdx,:);
        r_sic = w_rx_sic*(rxSig - H_prev*s_est);
        [b_est_current, s_est_current] = hardDemapping(symLabeling,QAM_Symbols,r_sic);
        s_est = [s_est;s_est_current];

        b_est=[b_est;b_est_current];
        s_est_soft = [s_est_soft;r_sic];

        detected_stream_index = [detected_stream_index;idxes_tobeDetected(streamIdx)];
        idxes_tobeDetected(streamIdx) = [];
    end

    if(ordering)
        s_est_hard_SIC(detected_stream_index,1)=s_est;
        b_est_hard_SIC(detected_stream_index,:)=b_est;
        s_est_soft_SIC(detected_stream_index,1)=s_est_soft;
    else
        s_est_hard_SIC = s_est;
        b_est_hard_SIC = b_est;
        s_est_soft_SIC = s_est_soft;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunktion of SIC algorithm
function streamIdx = get_first_SIC_stream_todetect(w_rx,channelH,noise_var)
    G = w_rx*channelH;
    P_Sig = abs(diag(G)).^2;
    P_ICI = (abs(diag((G-diag(diag(G)))*((G-diag(diag(G)))')))).^2;
    P_noise = noise_var*diag((w_rx*(w_rx')));
    SINR = P_Sig./(P_ICI+P_noise);
    [~,streamIdxSorted] = sort(SINR,'descend');
    streamIdx = streamIdxSorted(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
