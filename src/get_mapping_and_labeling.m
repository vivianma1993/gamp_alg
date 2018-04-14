function [m, l] = get_mapping_and_labeling(nb,labeling,mapPamQam)
% QAM/PAM/PSK Mapper 
% Inputs:  - nb: number of bits per symbol 
%          - labeling: Gray, anti-Gray, natural, random
%          - mapPamQam: PAM/QAM/PSK
%          - labArray: data index 
% Outputs: m : the bits 
%          l : the symbols 

if( nb == 1 )
    m = [0, 1]; % Gray, anti-Gray, natural, random 
    l = [-1, 1];
elseif( nb == 2 ) 

   mArray= [ get_labeling(nb,1,1); % Gray, QAM
             get_labeling(nb,0,1); % anti-Gray
                     0:3;  % natural
                     randperm(2^nb)-ones(1,2^nb); % random
                     get_labeling(nb,1,0); % Gray, PAM
                     get_labeling(nb,0,0); % anti-Gray
                     0:3;  % natural
                     randperm(2^nb)-ones(1,2^nb); % random
                     get_labeling(nb,1,0); % Gray, PSK
                     get_labeling(nb,0,0); % anti-Gray
                     0:3;  % natural
                     randperm(2^nb)-ones(1,2^nb)]; % random
                     
    lArray =  [ (1/sqrt(2))*[-1-1i, 1-1i, 1+1i, -1+1i ]; ... % QAM
                    (1/sqrt(5))*[-3+0i, -1+0i, 1+0i, 3+0]; ...  % PAM
                    (1/sqrt(2))*[-1-1i, 1-1i, 1+1i, -1+1i]]; % PSK
                    
    m = mArray(labeling+(mapPamQam-1)*4,1:end);
    l = lArray(mapPamQam,1:end);
    
elseif( nb == 3 )

    mArray = [ 0, 2, 6, 4, 1, 3, 7, 5; ...
                      0, 6, 4, 2, 1, 7, 5, 3; ...
                      0, 2, 4, 6, 1, 3, 5, 7;
                      randperm(2^nb)-ones(1,2^nb);
                      0, 1, 3, 2, 6, 7, 5, 4;
                      0, 7, 1, 6, 3, 4, 2, 5;
                      0, 1, 2, 3, 4, 5, 6, 7;
                      randperm(2^nb)-ones(1,2^nb);
                      0, 1, 3, 2, 6, 7, 5, 4;
                      0, 7, 1, 6, 3, 4, 2, 5;
                      0, 1, 2, 3, 4, 5, 6, 7;
                      randperm(2^nb)-ones(1,2^nb)];
                      
       lArray = [ (0.4597)*[-1-1i,1-1i,1+1i,-1+1i, -(1+sqrt(3))*1i,1+sqrt(3),(1+sqrt(3))*1i,-1-sqrt(3)];
                      (1/sqrt(21))*[-7, -5, -3, -1, 1, 3, 5, 7];
                      (1/sqrt(2))*(-1-1i), -1i, (1/sqrt(2))*(1-1i), 1, (1/sqrt(2))*(1+1i), 1i,(1/sqrt(2))*( -1+1i), -1 ];
                      
    m = mArray(labeling+(mapPamQam-1)*4,1:end);
    l = lArray(mapPamQam,1:end);
    
elseif( nb == 4 )
        mArray = [ get_labeling(nb,1,1); 
                          get_labeling(nb,0,1); 
                          0:15;
                          randperm(2^nb)-ones(1,2^nb);
                          get_labeling(nb,1,0); %0, 1, 3, 2, 6, 7, 5, 4, 12, 13, 15, 14, 10, 11, 9, 8; 
                          get_labeling(nb,0,0); 
                          0:15;
                          randperm(2^nb)-ones(1,2^nb);
                          get_labeling(nb,1,0); 
                          get_labeling(nb,0,0); 
                          0:15;
                          randperm(2^nb)-ones(1,2^nb)];
                          
       lArray =  [ (1/sqrt(10))*[-3-3i,-1-3i,1-3i,3-3i, -3-1i,-1-1i,1-1i,3-1i, -3+1i,-1+1i,1+1i,3+1i, -3+3i,-1+3i,1+3i,3+3i]; 
                     (1/sqrt(85))*[-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15];
                     exp(10i*pi/8), exp(11i*pi/8), -1i, exp(13i*pi/8), exp(14i*pi/8), exp(15i*pi/8), 1, exp(1i*pi/8), exp(2i*pi/8), exp(3i*pi/8), 1i, exp(5i*pi/8), exp(6i*pi/8), exp(7i*pi/8), -1, exp(9i*pi/8)];
                     
    m = mArray(labeling+(mapPamQam-1)*4,1:end);
    l = lArray(mapPamQam,1:end);
    
elseif( nb == 6 )

        mArray = [ get_labeling(nb,1,1);
                   get_labeling(nb,0,1); 
                   0:63;
                   randperm(2^nb)-ones(1,2^nb);
		               get_labeling(nb,1,0);
                   get_labeling(nb,0,0);
                   0:63;
                   randperm(2^nb)-ones(1,2^nb); 
		               get_labeling(nb,1,0);
                   get_labeling(nb,0,0);
                   0:63;
                   randperm(2^nb)-ones(1,2^nb)];
                          
       lArray =  [ (1/sqrt(42))*[ reshape(ones(8)*(-9-9i)+cumsum(ones(8)*2i,2)+cumsum(ones(8)*2),1,[]) ]; 
                   (1/sqrt(1365))*[-63:2:63];
                   cumprod(ones(1,64)*exp(i*pi/32)) ];
               
    m = mArray(labeling+(mapPamQam-1)*4,1:end);
    l = lArray(mapPamQam,1:end); 
 
 else % nb == 8
 
        mArray = [ get_labeling(nb,1,1);
                   get_labeling(nb,0,1);
                   0:255;
                   randperm(2^nb)-ones(1,2^nb);
                   get_labeling(nb,1,0);
                   get_labeling(nb,0,0);
                   0:255;
                   randperm(2^nb)-ones(1,2^nb);
                   get_labeling(nb,1,0);
                   get_labeling(nb,0,0);
                   0:255;
                   randperm(2^nb)-ones(1,2^nb)];
                          
       lArray =  [ (1/sqrt(170))*[ reshape(ones(16)*(-17-17i)+cumsum(ones(16)*2i,2)+cumsum(ones(16)*2),1,[]) ]; 
                   (1/sqrt(21845))*[-255:2:255];
                   cumprod(ones(1,256)*exp(i*pi/128)) ];
                          
    m = mArray(labeling+(mapPamQam-1)*4,1:end);
    l = lArray(mapPamQam,1:end);
end

end


