function [ I_txrx ] = getMutualInformationHardInput( tx_bits , rx_bits ,n_b)
%   Calculation of the mutual information of bitvectors tx_bit and rx_bit
%   under the assumption of beeing n_b multiplexed bitvectors
%
%   @return I_txrx scalar real double in range 0 to n_b
%
%   tx_bits = [ bitvector1(1),
%               bitvector2(1),
%               ...
%               bitvectornb(1),
%               bitvector1(2),
%               bitvector2(2),
%               ...
%               bitvectornb(2),
%               bitvector1(1),
%               ...
%               ]
%
%   same for rx_bits

    % ensure column vectors
    tx_bit_vec = tx_bits(:);
    rx_bit_vec = rx_bits(:);

    % trivial case
    if numel(tx_bit_vec) < 2*n_b
        I_txrx = 0;
        return;
    end

    % demultiplex vectors
    if n_b > 1
        tx_bit_mat = reshape(tx_bit_vec,n_b,[]).';
        rx_bit_mat = reshape(rx_bit_vec,n_b,[]).';
    else
        tx_bit_mat = tx_bit_vec;
        rx_bit_mat = rx_bit_vec;
    end

    % calculate mutual information for individual bitsteams
    I = -ones(n_b,1);
    for idx = 1:n_b
        I(idx) = getMutualInformationHardBitInput(tx_bit_mat(:,idx),...
            rx_bit_mat(:,idx));
    end

    % OPTIONAL: return mutual information of individual bitstreams
    %
    % final calculation
    I_txrx = sum(I(:));

end

function [ I_txrx ] = getMutualInformationHardBitInput( tx_bits ,rx_bits)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%   Calculation of the mutual information of bitvectors tx_bit and rx_bit
%
%   @return I_txrx scalar real double in range 0 to 1
%

    %esnsure column vectors
    tx_bit_vec = tx_bits(:);
    rx_bit_vec = rx_bits(:);

    % trivial cases
    if numel(tx_bit_vec) < 2
        I_txrx = 0;
        return;
    end
    if isequal(tx_bit_vec,rx_bit_vec) || isequal(tx_bit_vec,1 - rx_bit_vec)
        I_txrx = 1;
        return;
    end

    % avoid divisions by zero for short bit sequences by adding some
    % (indepentent) bits (decreases mutual information)
    % to decrease inpairments: duplication of input bits to sz bits
    tare = 1;
    sz = 1e5;
    if numel(tx_bit_vec) < sz
        dup_factor = floor(sz / numel(tx_bit_vec));
        tx_bit_vec = repmat(tx_bit_vec,[dup_factor,1]);
        rx_bit_vec = repmat(rx_bit_vec,[dup_factor,1]);
    end

    counter_total = numel(tx_bit_vec) + 4*tare;

    % count number of occourences
    counter_tx_0 = nnz(tx_bit_vec == 0) + 2*tare;
    counter_tx_1 = nnz(tx_bit_vec == 1) + 2*tare;
    counter_rx_0 = nnz(rx_bit_vec == 0) + 2*tare;
    counter_rx_1 = nnz(rx_bit_vec == 1) + 2*tare;

    counter_0_0 = nnz(tx_bit_vec == 0 & rx_bit_vec == 0) + tare;
    counter_1_0 = nnz(tx_bit_vec == 1 & rx_bit_vec == 0) + tare;
    counter_0_1 = nnz(tx_bit_vec == 0 & rx_bit_vec == 1) + tare;
    counter_1_1 = nnz(tx_bit_vec == 1 & rx_bit_vec == 1) + tare;

    % calculate probabilities of occurence
    p_tx_0 = counter_tx_0 / counter_total;
    p_tx_1 = counter_tx_1 / counter_total;
    p_rx_0 = counter_rx_0 / counter_total;
    p_rx_1 = counter_rx_1 / counter_total;

    p_0_0 = counter_0_0 / counter_total;
    p_1_0 = counter_1_0 / counter_total;
    p_0_1 = counter_0_1 / counter_total;
    p_1_1 = counter_1_1 / counter_total;

    I_0_0 = p_0_0 * log2(p_0_0 / p_tx_0 / p_rx_0);
    I_1_0 = p_1_0 * log2(p_1_0 / p_tx_1 / p_rx_0);
    I_0_1 = p_0_1 * log2(p_0_1 / p_tx_0 / p_rx_1);
    I_1_1 = p_1_1 * log2(p_1_1 / p_tx_1 / p_rx_1);

    % final calculation
    I_txrx = I_0_0 + I_1_0 + I_0_1 + I_1_1;

    if isempty(I_txrx)
        disp('Changing I from [] to 0');
        I_txrx = 0;
    end
    if isnan(I_txrx)
        disp('I_txrx = NaN');
    end

end % getMutualInformationHardBitInput
