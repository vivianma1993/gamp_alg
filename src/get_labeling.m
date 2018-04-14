function labelArray = get_labeling(nb,grey,qam)

nb_codes = 2^nb;
Gray_tab = zeros(nb_codes,nb);
Gray_tab(1,:)=-1;

if (grey == 1)
    for i = 0:nb-1
        for j = 0:(2^i)-1
            for s = 0:nb-1
                Gray_tab(2^(i+1)-1-j+1,s+1) = Gray_tab(j+1,s+1);
            end
            Gray_tab(2^(i+1)-1-j+1,i+1) = 1;
        end
    end
    
else
    % anti-grey
    for i = 0:nb-2
        for j = 0:(2^i)-1
            for s = 0:nb-2
                Gray_tab(2^(i+1)-1-j+1,s+1) = Gray_tab(j+1,s+1);
            end
            Gray_tab(2^(i+1)-1-j+1,nb) = -1;
            Gray_tab(2^(i+1)-1-j+1,i+1) = 1;
        end
    end
    
    % second: rearrange and insert inverted versions (see Hamming-book)
    for i = nb_codes/2-1:-1:0
        for s = 0:nb-1
            Gray_tab(2*i+2,s+1) = -Gray_tab(i+1,s+1);
            Gray_tab(2*i+1,s+1) = Gray_tab(i+1,s+1);
        end
    end
    
end

Gray_tab(Gray_tab==-1)=0;

Gray_dec = zeros(nb_codes,nb);
for i=1:nb
    Gray_dec(:,i)=2^(i-1);
end
Gray_dec = cumsum(Gray_dec.*Gray_tab,2);
labelArray = Gray_dec(:,nb)';

if(qam==1&&nb>2)
    labelArray = reshape(labelArray,sqrt(nb_codes),[])';
    labelArray(2:2:end) = fliplr(labelArray(2:2:end,:));
    labelArray = labelArray(:);
    labelArray =  labelArray';
end

end

