function [ muEfficiency,SINR,SNR ] = caculateMultiUserEfficiency(signal,mui,noise_term)
% Caculates Multiuser effciency = SINR/SNR
% input: signal  
%        mui multiuser interference 
%        noise tern
% output: multiuser efficiency

noise_power = sum(abs(noise_term).^2);
interference_power = sum(abs(mui).^2);
signal_power = sum(abs(signal).^2);
SNR=signal_power/noise_power;
SINR= signal_power/(noise_power+interference_power);
muEfficiency =SINR/SNR;
end

