function LLR_bits = softDemappingSISO( sym_labeling,symbol_set,rxSig,NoiseVar,algorithm,channel_weight,La)
% input: the symbols and its labeling and La (app of bits)
%        received signal, noise vairance, and algorithm (1 exact, 2 max log, 3 max-log jacob)
%        the weighting factor of channel (flat rayleight fading coeff.)

%[m, l] = get_mapping_and_labeling(nb,labeling,mapPamQam);
ns=length(symbol_set); % number of symbols
nb=log2(ns);

bit=zeros(ns,nb);
for k=1:ns
    bit(k,:)=bitget(sym_labeling(k),1:1:nb);
end

% compute L-values
LLR_bits = zeros(length(rxSig),nb);

% kL: Nummer von L-Wert bei QPSK z.B. 2 L-Werte zu berechnen.
for kL=1:nb
    
    l_o = symbol_set(bit(:,kL)==1);
    l_u = symbol_set(bit(:,kL)==0);
    
    % Lch
    variance = NoiseVar/2;
    c=-1./(2*variance.*abs(1./(channel_weight.^2)));
    
    if(algorithm==1)
        exp_o = 0;
        exp_u = 0;
        for n=1:(ns/2)
            exp_o = exp_o + exp(c.*(abs(rxSig-l_o(n))).^2);
            exp_u = exp_u + exp(c.*(abs(rxSig-l_u(n))).^2);
        end
        y_Lch = log( exp_o ./ exp_u );
        
    else
        
        max_o = -1000;
        max_u = -1000;
        for n=1:(ns/2)
            if(algorithm==2)
                r_o = log(1 + exp(-abs( max_o - c.*(rxSig-l_o(n)).^2 )));
                r_u = log(1 + exp(-abs( max_u - c.*(rxSig-l_u(n)).^2 )));
            end
            max_o = max ( max_o, c.*abs(rxSig-l_o(n)).^2 );
            max_u = max ( max_u, c.*abs(rxSig-l_u(n)).^2 );
            if(algorithm==2)
                max_o = max_o + r_o;
                max_u = max_u + r_u;
            end
        end
        y_Lch = max_o - max_u;
    end
    
    y_L = y_Lch + La;
    
    LLR_bits(1:end,kL) = y_L;
    
end

end