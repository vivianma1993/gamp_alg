function [ I_txrx ] = getMutualInformationSoftInput( tx_bits , LLR_bits ,n_b)
%   Calculation of the mutual information of bitvectors tx_bit and rx_bit
%   under the assumption of beeing n_b multiplexed bitvectors
%
%   @return I_txrx scalar real double in range 0 to n_b
%
%
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

%%

% ensure column vectors
tx_bit_vec = tx_bits(:);
LLR_bits_vec = LLR_bits(:);


% trivial case
if numel(tx_bit_vec) < 2*n_b
    I_txrx = 0;
    return;
end


% demultiplex vectors
if n_b > 1
    tx_bit_mat = reshape(tx_bit_vec,n_b,[]).';
    LLR_bits_mat = reshape(LLR_bits_vec,n_b,[]).';
else
    tx_bit_mat = tx_bit_vec;
    LLR_bits_mat = LLR_bits_vec;
end

% calculate mutual information for individual bitsteams
I = -ones(n_b,1);
for idx = 1:n_b
    I(idx) = measure_mutual_information_histogram(...
        LLR_bits_mat(:,idx),tx_bit_mat(:,idx));
end

% OPTIONAL: return mutual information of individual bitstreams
%
% final calculation
I_txrx = sum(I(:));


end
