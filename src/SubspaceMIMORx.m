function [ s_block_all_iterations ,s_block_all_iterations_LLR ] = SubspaceMIMORx( rxSig,channelH,BlockSize,QAM_Symbols,withPIC,NiterMax,useDistMetric,MLorAPP,sym_labeling,NoiseVar,maxLog,La,givenInitialEstimates,iniEstimates)
%SUBSPACEMIMORX Summary of this function goes here
%   Detailed explanation goes here
% input received signal, channel matrix, size of block, QAM_Symbols,apply
% PIC?, max Iteration of PIC, using DistMetric for updating ?, MLorAPP (0 ML/1 APP)

[Nr,Nt]=size(channelH);

Nr_Blocks = ceil(Nt/BlockSize);
lastBlockSize = Nt-(Nr_Blocks-1)*BlockSize;

% construct block channel matrices for subspace detection
H_block=zeros(Nr,Nt,Nr_Blocks);
Q_block=zeros(Nr,Nr,Nr_Blocks);
R_block_real=zeros(Nr,Nt,Nr_Blocks);
R_block=zeros(Nt,Nt,Nr_Blocks);
y_block=zeros(Nr,Nr_Blocks);
s_block_all = zeros(Nt,1);
s_block_all_iterations = zeros(Nt,NiterMax+1);
nb= log2(length(QAM_Symbols));
s_block_all_iterations_LLR = zeros(Nt,nb,NiterMax+1);
s_block_all_LLR = zeros(Nt,nb);
for iBlock= 1:(Nr_Blocks-1)
    H_block(:,:,(Nr_Blocks-iBlock+1))=circshift(channelH,[0,(iBlock-1)*BlockSize]);
    [Q_block(:,:,(Nr_Blocks-iBlock+1)), R_block_real(:,:,(Nr_Blocks-iBlock+1))] = qr(H_block(:,:,(Nr_Blocks-iBlock+1)));
    R_block(:,:,(Nr_Blocks-iBlock+1))=R_block_real(1:Nt,:,(Nr_Blocks-iBlock+1));
    y_block(:,(Nr_Blocks-iBlock+1)) = Q_block(:,:,(Nr_Blocks-iBlock+1))'*rxSig;
    if(givenInitialEstimates)
        s_block_all = iniEstimates;
    else
        if(MLorAPP)
            LLR_bits_vec = softDemapping(sym_labeling,QAM_Symbols,y_block((end-Nr+Nt-BlockSize+1):(end-Nr+Nt),(Nr_Blocks-iBlock+1)),NoiseVar,maxLog,R_block((end-BlockSize+1):end,(end-BlockSize+1):end,(Nr_Blocks-iBlock+1)),La);
            LLR_sym_vec = reshape(LLR_bits_vec,[nb BlockSize]);
            s_block_all_LLR((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),:) = LLR_sym_vec.';
            s_block_bits = logical((sign(LLR_sym_vec.')+1)/2);
            sym_idxes = bi2de(s_block_bits); 
            s_block_all((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),1) = mapping(sym_labeling,QAM_Symbols,sym_idxes);
        else
            [s_block_all((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),1), ~] = ML_rx(y_block((end-Nr+Nt-BlockSize+1):(end-Nr+Nt),(Nr_Blocks-iBlock+1)),R_block((end-BlockSize+1):end,(end-BlockSize+1):end,(Nr_Blocks-iBlock+1)),QAM_Symbols);
        end
    end
end

% last block, first few streams
iBlock=Nr_Blocks;
H_block(:,:,(Nr_Blocks-iBlock+1))=circshift(channelH,[0,(iBlock-1)*BlockSize]);
[Q_block(:,:,(Nr_Blocks-iBlock+1)), R_block_real(:,:,(Nr_Blocks-iBlock+1))] = qr(H_block(:,:,(Nr_Blocks-iBlock+1)));
R_block(:,:,(Nr_Blocks-iBlock+1))=R_block_real(1:Nt,:,(Nr_Blocks-iBlock+1));
y_block(:,(Nr_Blocks-iBlock+1)) = Q_block(:,:,(Nr_Blocks-iBlock+1))'*rxSig;
if(~givenInitialEstimates)
    if(MLorAPP)
        LLR_bits_vec = softDemapping(sym_labeling,QAM_Symbols,y_block((end-Nr+Nt-lastBlockSize+1):(end-Nr+Nt),1),NoiseVar,maxLog,R_block((end-lastBlockSize+1):end,(end-lastBlockSize+1):end,(Nr_Blocks-iBlock+1)),La);
            LLR_sym_vec = reshape(LLR_bits_vec,[nb lastBlockSize]);
            s_block_all_LLR(1:lastBlockSize,:) = LLR_sym_vec.';
            s_block_bits = logical((sign(LLR_sym_vec.')+1)/2);
            sym_idxes = bi2de(s_block_bits); 
            s_block_all(1:lastBlockSize,1) = mapping(sym_labeling,QAM_Symbols,sym_idxes);
    else
        [s_block_all(1:lastBlockSize,1), ~] = ML_rx(y_block((end-Nr+Nt-lastBlockSize+1):(end-Nr+Nt),1),R_block((end-lastBlockSize+1):end,(end-lastBlockSize+1):end,1),QAM_Symbols);
    end
end

s_block_all_iterations(:,1) = s_block_all;
s_block_all_iterations_LLR(:,:,1) = s_block_all_LLR;
% if subspace PIC used
if(withPIC)
    % iterative subspace with PIC
    y_block_new = zeros(Nr,Nr_Blocks);
    s_block_all_new = zeros(Nt,1);
    dist = zeros(Nr_Blocks,1);
    dist_vec = zeros(Nr_Blocks-1,Nr);
    dist_vec_last = zeros(1,Nr);
    dist_new = zeros(Nr_Blocks,1);
    dist_vec_new = zeros(Nr_Blocks-1,Nr);
    dist_vec_last_new = zeros(1,Nr);
    exitIter = 1;
    s_block_all_LLR_new = zeros(size(s_block_all_LLR));
    for iter =1:NiterMax
        % perform interference cancellation
        for iBlock= 1:(Nr_Blocks-1)
            s_block_interference = circshift(s_block_all,[(iBlock-1)*BlockSize, 0]);
            s_block_interference((end-BlockSize+1):end)=[];
            y_block_new(:,(Nr_Blocks-iBlock+1))=y_block(:,(Nr_Blocks-iBlock+1)) - R_block_real(:,1:((end-BlockSize)),(Nr_Blocks-iBlock+1))*s_block_interference;
            if(MLorAPP)
                 LLR_bits_vec = softDemapping(sym_labeling,QAM_Symbols,y_block_new(:,(Nr_Blocks-iBlock+1)),NoiseVar,maxLog,R_block_real(:,(end-BlockSize+1):end,(Nr_Blocks-iBlock+1)),La);
                 LLR_sym_vec = reshape(LLR_bits_vec,[nb BlockSize]);
                 s_block_all_LLR_new((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),:) = LLR_sym_vec.';
                 s_block_bits = logical((sign(LLR_sym_vec.')+1)/2);
                 sym_idxes = bi2de(s_block_bits); 
                 s_block_all_new((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),1) = mapping(sym_labeling,QAM_Symbols,sym_idxes);
                 dist_new(Nr_Blocks-iBlock+1) = sum(abs(y_block_new(:,(Nr_Blocks-iBlock+1))-R_block_real(:,(end-BlockSize+1):end,Nr_Blocks-iBlock+1)*s_block_all_new((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize))).^2);
            else
                [s_block_all_new((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),1), dist_new(Nr_Blocks-iBlock+1),dist_vec_new((Nr_Blocks-iBlock+1),:)] = ML_rx(y_block_new(:,(Nr_Blocks-iBlock+1)),R_block_real(:,(end-BlockSize+1):end,(Nr_Blocks-iBlock+1)),QAM_Symbols);
            end
            if(iter==1)
                dist(Nr_Blocks-iBlock+1) = sum(abs(y_block_new(:,(Nr_Blocks-iBlock+1))-R_block_real(:,(end-BlockSize+1):end,Nr_Blocks-iBlock+1)*s_block_all((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize))).^2);
                dist_vec(Nr_Blocks-iBlock+1,:) = (abs(y_block_new(:,(Nr_Blocks-iBlock+1))-R_block_real(:,(end-BlockSize+1):end,Nr_Blocks-iBlock+1)*s_block_all((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize))).^2);
            end
        end
        % last block, first few streams
        s_block_interference = s_block_all;
        s_block_interference(1:lastBlockSize)=[];
        y_block_new(:,1)=y_block(:,1) - R_block_real(:,1:((end-lastBlockSize)),1)*s_block_interference;
        if(MLorAPP)
            LLR_bits_vec = softDemapping(sym_labeling,QAM_Symbols,y_block_new(:,1),NoiseVar,maxLog,R_block_real(:,(end-lastBlockSize+1):end,1),La);
             LLR_sym_vec = reshape(LLR_bits_vec,[nb lastBlockSize]);
                 s_block_all_LLR_new(1:lastBlockSize,:) = LLR_sym_vec.';
                 s_block_bits = logical((sign(LLR_sym_vec.')+1)/2);
                 sym_idxes = bi2de(s_block_bits); 
                 s_block_all_new(1:lastBlockSize,1) = mapping(sym_labeling,QAM_Symbols,sym_idxes);
                 dist_new(1) = sum(abs(y_block_new(:,1)-R_block_real(:,(end-lastBlockSize+1):end,1)*s_block_all_new(1:lastBlockSize)).^2);
        else
            [s_block_all_new(1:lastBlockSize,1), dist_new(1),dist_vec_last_new(1,:)] = ML_rx(y_block_new(:,1),R_block_real(:,(end-lastBlockSize+1):end,1),QAM_Symbols);
        end
        if(iter==1)
            dist(1) = sum(abs(y_block_new(:,1)-R_block_real(:,(end-lastBlockSize+1):end,1)*s_block_all(1:lastBlockSize)).^2);
            dist_vec_last(1,:) = (abs(y_block_new(:,1)-R_block_real(:,(end-lastBlockSize+1):end,1)*s_block_all(1:lastBlockSize)).^2);
        end
        
        %update estimates
        if(useDistMetric==1) % use the distance metric for each subspace
            for iBlock= 1:(Nr_Blocks-1)
                if(dist_new(Nr_Blocks-iBlock+1)<dist(Nr_Blocks-iBlock+1))
                    s_block_all((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),1) = s_block_all_new((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),1);
                    s_block_all_LLR((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),:) = s_block_all_LLR_new((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),:);
                    dist(Nr_Blocks-iBlock+1) = dist_new(Nr_Blocks-iBlock+1);
                    exitIter = 0;
                end
            end
            if(dist_new(1)<dist(1))
                s_block_all(1:lastBlockSize,1) = s_block_all_new(1:lastBlockSize,1);
                s_block_all_LLR(1:lastBlockSize,:) = s_block_all_LLR_new(1:lastBlockSize,:);
                dist(1) = dist_new(1);
                exitIter=0;
            end
            s_block_all_iterations(:,iter+1) = s_block_all;
            s_block_all_iterations_LLR(:,:,iter+1) = s_block_all_LLR;
            if(exitIter)
                s_block_all_iterations(:,(iter+2):end) = repmat(s_block_all,[1 NiterMax-iter]);
                s_block_all_iterations_LLR(:,:,(iter+2):end) = repmat(s_block_all_LLR,[1 1 NiterMax-iter]);
                break;
            end
            exitIter =1;
        elseif(useDistMetric==2) % use APP for each bits
            for iBlock= 1:(Nr_Blocks-1)
                if(sum(sum(abs(s_block_all_LLR_new((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),:))))>sum(sum(abs(s_block_all_LLR((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),:)))))
                    s_block_all_LLR((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),:) = s_block_all_LLR_new((Nt-iBlock*BlockSize+1):(Nt-(iBlock-1)*BlockSize),:);
                    exitIter = 0;
                end
            end
            if(sum(sum(abs(s_block_all_LLR_new(1:lastBlockSize,:))))>sum(sum(abs(s_block_all_LLR(1:lastBlockSize,:)))))
                    s_block_all_LLR(1:lastBlockSize,:) = s_block_all_LLR_new(1:lastBlockSize,:);
                    exitIter = 0;
            end
              s_block_bits = logical((sign(s_block_all_LLR)+1)/2);
                 sym_idxes = bi2de(s_block_bits); 
                 s_block_all = mapping(sym_labeling,QAM_Symbols,sym_idxes);
              s_block_all_iterations(:,iter+1) = s_block_all;
              s_block_all_iterations_LLR(:,:,iter+1) = s_block_all_LLR;

            if(exitIter)
                s_block_all_iterations(:,(iter+2):end) = repmat(s_block_all,[1 NiterMax-iter]);
                s_block_all_iterations_LLR(:,:,(iter+2):end) = repmat(s_block_all_LLR,[1 1 NiterMax-iter]);
                break;
            end
            exitIter =1;
        else
            s_block_all = s_block_all_new;
            s_block_all_LLR = s_block_all_LLR_new;
            s_block_all_iterations(:,iter+1) = s_block_all;
            s_block_all_iterations_LLR(:,:,iter+1) = s_block_all_LLR;
        end 
    end
end
end

