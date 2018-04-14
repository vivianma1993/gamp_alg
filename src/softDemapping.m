function LLR_bits = softDemapping(sym_labeling,symbol_set,rxSig,NoiseVar,algorithm,channel_weight,La)
% input: the symbols and its labeling and La (app of bits)
%        received signal, noise vairance, and algorithm (1 exact, 2 max log, 3 max-log jacob)
%        the weighting factor of channel (flat rayleight fading coeff.)

	%[m, l] = get_mapping_and_labeling(nb,labeling,mapPamQam);
	%Nstream = size(channel_weight,2);

	[Nr,Nt]=size(channel_weight);

	ns=length(symbol_set); % number of symbols
	nb=log2(ns);
	n_bits_per_vec = nb*Nt;
	nr_llr_hyp_per_bit = (1/2 *ns^Nt);
	% preallocation of sets
	sset_0 = zeros(Nt,nr_llr_hyp_per_bit,n_bits_per_vec);
	sset_1 = zeros(Nt,nr_llr_hyp_per_bit,n_bits_per_vec);

	bit_label=zeros(ns,nb);
	for k=1:ns
		bit_label(k,:)=bitget(sym_labeling(k),1:1:nb);
	end

	% iteration over all bit indices
	for m = 1:n_bits_per_vec
		% calculation of sets of vector symbols with x_m = {0,1}
		[sset_0(:,:,m),sset_1(:,:,m)] = get_vec_symbols_by_bit_index(m,Nt,bit_label,symbol_set);
	end
	% preallocation of bit vectors
	LLR_bits = -ones(n_bits_per_vec,1);

	% iteration over all bits in this tranmitted vector symbol
	for m = 1:n_bits_per_vec

		% get set of vector symbols with x_m = {0,1}
		s_0 = sset_0(:,:,m);
		s_1 = sset_1(:,:,m);


		% matrix containing the y-H*s for all s in s_0
		sum_temp_0 = reshape(sum((repmat(channel_weight,[1,1,nr_llr_hyp_per_bit]) .* ...
			 permute(repmat(s_0,[1,1,Nr]),[3,1,2])),2),[Nr,nr_llr_hyp_per_bit]);
		diff_0_vec = repmat(rxSig,[1,nr_llr_hyp_per_bit]) - sum_temp_0;

%         diff_0_vec = repmat(rxSig,[1,nr_llr_hyp_per_bit]) - ...
%             reshape(sum((repmat(channel_weight,[1,1,nr_llr_hyp_per_bit]) .* ...
%             permute(repmat(s_0,[1,1,Nr]),[3,1,2])),2),[Nr,nr_llr_hyp_per_bit]);
		% matrix containing the y-H*s for all s in s_1
		sum_temp_1 = reshape(sum((repmat(channel_weight,[1,1,nr_llr_hyp_per_bit]) .* ...
			 permute(repmat(s_1,[1,1,Nr]),[3,1,2])),2),[Nr,nr_llr_hyp_per_bit]);
		diff_1_vec = repmat(rxSig,[1,nr_llr_hyp_per_bit]) - sum_temp_1;

%         diff_1_vec = repmat(rxSig,[1,nr_llr_hyp_per_bit]) - ...
%             reshape(sum((repmat(channel_weight,[1,1,nr_llr_hyp_per_bit]) .* ...
%             permute(repmat(s_1,[1,1,Nr]),[3,1,2])),2),[Nr,nr_llr_hyp_per_bit]);


		% vector containing all |y-H*s|^2 for all s in s_0, s_1
		diff_0 = sum(abs(diff_0_vec).^2,1);
		diff_1 = sum(abs(diff_1_vec).^2,1);

		if(algorithm==1)
			LLR_bits(m) = log(sum(exp(- diff_1 / NoiseVar)) / sum(exp(- diff_0 / NoiseVar)));
		else
			LLR_bits(m) = (min(diff_0) - min(diff_1)) / NoiseVar;
		end
	end
	LLR_bits = LLR_bits +La;
	% end bitwise soft demapping -------------------------------------
end

function [ vec_symbol_mat_0 , vec_symbol_mat_1 ] = get_vec_symbols_by_bit_index( bit_no , tx_total ,bit_label, qam_symbol_vec)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%
% returns two matrices whose columns are vector symbols that could result
% from a mapping of bit vectors that have the value 0 or 1 at the
% position bit_no (using MSBF) using tx_total transmit antennas
%
% bit_no range is from 1 to (log2(numel(qam_symbol_vec))*tx_total)
%
% tx_no range is from 1 to tx_total
%
% size of vec_symbol_matrix_0 and vec_symbol_matrix_1
%                     is [tx_total, (1/2 * numel(qam_symbol_vec)^tx_total)]
%
%   JPF 2015-02-19

	% number of constelation symbols
	n_symbols = numel(qam_symbol_vec);
	% bit per symbol
	M = log2(n_symbols);

	% get index number of corresponding transmit antenna
	tx_no    = floor((bit_no-1)/M)+1;
	% get index of current bit bit_no at its antenna
	bit_no_a = mod(bit_no-1,M)+1;

	% mapping indices --> symbols
	symbols_a_0 = qam_symbol_vec(bit_label(:,bit_no_a)==0);
	symbols_a_1 = qam_symbol_vec(bit_label(:,bit_no_a)==1);


	% preallocation of output matrix
	vec_symbol_mat_0 = zeros(tx_total,1/2 * n_symbols^tx_total);
	vec_symbol_mat_1 = zeros(tx_total,1/2 * n_symbols^tx_total);

	% add symbols of the antenna of interrest to matrix multiple times
	vec_symbol_mat_0(tx_no,:) = repmat(symbols_a_0,[1,n_symbols^(tx_total-1)]);
	vec_symbol_mat_1(tx_no,:) = repmat(symbols_a_1,[1,n_symbols^(tx_total-1)]);

	% fill matrix with all possible permutations of symbols -------------------

	% all possible symbol indices
	all_ind = 1:n_symbols;

	% preallocation of padding matrix
	fill_matrix = zeros(tx_total-1,1/2 * n_symbols^tx_total);

	% iteration over all other transmit antennas
	for t = 1:(tx_total-1)
		tmp      = repmat(all_ind,1/2*n_symbols^t,n_symbols^(tx_total-1-t));
		curr_row = (tmp(:)).'; %'
		% mapping indices --> symbols
		fill_matrix(t,:) = qam_symbol_vec(curr_row);
	end

	vec_symbol_mat_0([1:(tx_no-1),(tx_no+1):end],:) = fill_matrix;
	vec_symbol_mat_1([1:(tx_no-1),(tx_no+1):end],:) = fill_matrix;

end
