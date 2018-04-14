function [ soft_syms,second_moment] = convert_LLR_to_soft_symbols(LLR,symLabeling, symbols)
%CONVERT_LLR_TO_SOFT_SYMBOLS

% this function converts the LLR to soft symbols based on used
% constellation and bit labeling
%   INPUT: LLR is matrix
%          [ sym1_bit1 sym1_bit2 ... ; MSB right most
%            sym2_bit1 sym2_bit2 ... ;
%            ...                        ]
%           LLR_bit = ln(P_0/P_1)

N_bit_symbol = log2(length(symbols));
Nsym=size(LLR,1);
soft_syms = zeros(Nsym,1);
second_moment = zeros(Nsym,1);
% convert LLR for bits to prob.
P_0 = 1./(1+exp(LLR));
P_1 = 1 - P_0;

for iSymConst = 0:(2^N_bit_symbol-1)
    label = de2bi(iSymConst,N_bit_symbol);
    if(label(1)==0)
        prob = P_0(:,1);
    else
        prob = P_1(:,1);
    end
    for iBit = 2:N_bit_symbol
        if(label(iBit)==0)
            prob =prob.*P_0(:,iBit);
        else
            prob =prob.*P_1(:,iBit);
        end
    end
    soft_syms = soft_syms + prob*mapping(symLabeling, symbols, iSymConst);
    second_moment = second_moment + prob*(abs(mapping(symLabeling, symbols, iSymConst)))^2;
end

end

