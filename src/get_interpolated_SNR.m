function [ requiredSNR ] = get_interpolated_SNR(information_rate,SNRs,maxRate,codeRate )
% this function caculates the require SNR for achieving target information
% rate given as coderate, linear interpolation is used 

targetRate = maxRate*codeRate;

if(targetRate>max(information_rate))
    requiredSNR=1000; % unachievable 
else
    idx = find(information_rate>targetRate,1);
    if(idx==1)
       requiredSNR=-1000; % not enough data  
    else
       slope_line = (information_rate(idx)-information_rate(idx-1))/(SNRs(idx)-SNRs(idx-1));
       offset_c = information_rate(idx) - slope_line*SNRs(idx);
       requiredSNR = (targetRate - offset_c)/slope_line;
    end
end
end

