function tx_sym =  mapping(sym_labeling,symbol_set,tx_sym_idx)
% QAM/PAM/PSK Mapper
% Inputs:  - nb: number of bits per symbol, max supported to 8
%          - labeling: Gray, anti-Gray, natural, random
%          - mapPamQam: PAM/QAM/PSK
%          - labArray: data index
% Outputs: the transmitting symbols 

l_max=length(symbol_set); % max index of symbol
%[m, l] = get_mapping_and_labeling(nb, labeling, mapPamQam);

tx_sym = tx_sym_idx;

for n = 1:l_max
    tx_sym(tx_sym==sym_labeling(n)) = symbol_set(n);
end

end
