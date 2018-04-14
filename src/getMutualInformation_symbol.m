
function [ MI_symbol_txrx ] = getMutualInformation_symbol( tx_symbols,rx_symbols )
 tx_symbol_vec = tx_symbols(:);
 rx_symbol_vec = rx_symbols(:);
  if isequal(tx_symbol_vec,rx_symbol_vec)
    MI_symbol_txrx = 1;
    return;
  end
 tare = 1;
 sz = 1e5;
 

 if numel(tx_symbol_vec) < sz
    dup_factor = floor(sz / numel(tx_symbol_vec));
    tx_symbol_vec = repmat(tx_symbol_vec,[dup_factor,1]);
    rx_symbol_vec = repmat(rx_symbol_vec,[dup_factor,1]);
end


counter_total = numel(tx_symbol_vec) + 16*tare;
% count number of occurance
counter_tx_0 = nnz(tx_symbol_vec == -0.7071-0.7071i) + 4*tare;
counter_tx_1 = nnz(tx_symbol_vec == 0.7071-0.7071i) + 4*tare;
counter_tx_2 = nnz(tx_symbol_vec == -0.7071+0.7071i) + 4*tare;
counter_tx_3 = nnz(tx_symbol_vec == 0.7071+0.7071i) + 4*tare;

counter_rx_0 = nnz(rx_symbol_vec == -0.7071-0.7071i) + 4*tare;
counter_rx_1 = nnz(rx_symbol_vec== 0.7071-0.7071i) + 4*tare;
counter_rx_2 = nnz(rx_symbol_vec == -0.7071+0.7071i) + 4*tare;
counter_rx_3 = nnz(rx_symbol_vec == 0.7071+0.7071i) + 4*tare;

counter_0_0 = nnz(tx_symbol_vec==-0.7071-0.7071i & rx_symbol_vec==-0.7071-0.7071i)+tare;
counter_0_1 = nnz(tx_symbol_vec==-0.7071-0.7071i & rx_symbol_vec== 0.7071-0.7071i)+tare;
counter_0_2 = nnz(tx_symbol_vec==-0.7071-0.7071i & rx_symbol_vec==-0.7071+0.7071i)+tare;
counter_0_3 = nnz(tx_symbol_vec==-0.7071-0.7071i & rx_symbol_vec== 0.7071+0.7071i)+tare;

counter_1_0 = nnz(tx_symbol_vec== 0.7071-0.7071i & rx_symbol_vec==-0.7071-0.7071i)+tare;
counter_1_1 = nnz(tx_symbol_vec== 0.7071-0.7071i & rx_symbol_vec== 0.7071-0.7071i)+tare;
counter_1_2 = nnz(tx_symbol_vec== 0.7071-0.7071i & rx_symbol_vec==-0.7071+0.7071i)+tare;
counter_1_3 = nnz(tx_symbol_vec== 0.7071-0.7071i & rx_symbol_vec== 0.7071+0.7071i)+tare;

counter_2_0 = nnz(tx_symbol_vec==-0.7071+0.7071i & rx_symbol_vec==-0.7071-0.7071i)+tare;
counter_2_1 = nnz(tx_symbol_vec==-0.7071+0.7071i & rx_symbol_vec== 0.7071-0.7071i)+tare;
counter_2_2 = nnz(tx_symbol_vec==-0.7071+0.7071i & rx_symbol_vec==-0.7071+0.7071i)+tare;
counter_2_3 = nnz(tx_symbol_vec==-0.7071+0.7071i & rx_symbol_vec==0.7071+0.7071i)+tare;

counter_3_0 = nnz(tx_symbol_vec== 0.7071+0.7071i & rx_symbol_vec==-0.7071-0.7071i)+tare;
counter_3_1 = nnz(tx_symbol_vec== 0.7071+0.7071i & rx_symbol_vec== 0.7071-0.7071i)+tare;
counter_3_2 = nnz(tx_symbol_vec== 0.7071+0.7071i & rx_symbol_vec==-0.7071+0.7071i)+tare;
counter_3_3 = nnz(tx_symbol_vec== 0.7071+0.7071i & rx_symbol_vec== 0.7071+0.7071i)+tare;

%calculate probabilities of accurance
p_tx_0 = counter_tx_0 / counter_total;
p_tx_1 = counter_tx_1 / counter_total;
p_tx_2 = counter_tx_2 / counter_total;
p_tx_3 = counter_tx_3 / counter_total;

p_rx_0 = counter_rx_0 / counter_total;
p_rx_1 = counter_rx_0 / counter_total;
p_rx_2 = counter_rx_0 / counter_total;
p_rx_3 = counter_rx_0 / counter_total;

p_0_0 = counter_0_0 / counter_total;
p_0_1 = counter_0_1 / counter_total;
p_0_2 = counter_0_2 / counter_total;
p_0_3 = counter_0_3 / counter_total;

p_1_0 = counter_1_0 / counter_total;
p_1_1 = counter_1_1 / counter_total;
p_1_2 = counter_1_2 / counter_total;
p_1_3 = counter_1_3 / counter_total;

p_2_0 = counter_2_0 / counter_total;
p_2_1 = counter_2_1 / counter_total;
p_2_2 = counter_2_2 / counter_total;
p_2_3 = counter_2_3 / counter_total;

p_3_0 = counter_3_0 / counter_total;
p_3_1 = counter_3_1 / counter_total;
p_3_2 = counter_3_2 / counter_total;
p_3_3 = counter_3_3 / counter_total;

I_0_0 = p_0_0 * log2(p_0_0/p_tx_0/p_rx_0);
I_0_1 = p_0_1 * log2(p_0_1/p_tx_0/p_rx_1);
I_0_2 = p_0_2 * log2(p_0_2/p_tx_0/p_rx_2);
I_0_3 = p_0_3 * log2(p_0_3/p_tx_0/p_rx_3);

I_1_0 = p_0_0 * log2(p_1_0/p_tx_1/p_rx_0);
I_1_1 = p_0_1 * log2(p_1_1/p_tx_1/p_rx_1);
I_1_2 = p_0_2 * log2(p_1_2/p_tx_1/p_rx_2);
I_1_3 = p_0_3 * log2(p_1_3/p_tx_1/p_rx_3);

I_2_0 = p_0_0 * log2(p_2_0/p_tx_2/p_rx_0);
I_2_1 = p_0_1 * log2(p_2_1/p_tx_2/p_rx_1);
I_2_2 = p_0_2 * log2(p_2_2/p_tx_2/p_rx_2);
I_2_3 = p_0_3 * log2(p_2_3/p_tx_2/p_rx_3);

I_3_0 = p_3_0 * log2(p_3_0/p_tx_3/p_rx_0);
I_3_1 = p_3_1 * log2(p_3_1/p_tx_3/p_rx_1);
I_3_2 = p_3_2 * log2(p_3_2/p_tx_3/p_rx_2);
I_3_3 = p_3_3 * log2(p_3_3/p_tx_3/p_rx_3);

%final calculation
MI_txrx= I_0_0 + I_0_1 + I_0_2 + I_0_3 + I_1_0 + I_1_1 + I_1_2 + I_1_3 + ...
         I_2_0 + I_2_1 + I_2_2 + I_2_3 + I_3_0 + I_3_1 + I_3_2 + I_3_3;  

 if isempty(MI_txrx)
     disp('changing I from [] to 0');
     MI_txrx = 0;
 end
 
 if isnan(MI_txrx)
     disp('MI_txrx=NaN');
 end 
  MI_symbol_txrx = sum(MI_txrx(:));

end
