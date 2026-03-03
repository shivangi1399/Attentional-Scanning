function [met] = subfunc_pairwisecosine(inp)

% See Landau ... Fries 2015 Current Biology and Ni ... Fries 2017 Neuron

v1 = inp.v1;
v2 = inp.v2;

clear inp

met = 0;

ntrl1 = size(v1,1);
ntrl2 = size(v2,1);

nsiz = numel(size(v1));

if nsiz == 2
    
    for k = 1:ntrl1
        
        for l = 1:ntrl2
            
            met = met+cos(angle(v1(k,:).*conj(v2(l,:))));
            
        end
        
    end
    
elseif nsiz == 3
    
    for k = 1:ntrl1
        
        for l = 1:ntrl2
            
            met = met+cos(angle(squeeze(v1(k,:,:)).*conj(squeeze(v2(l,:,:)))));
            
        end
        
    end
    
end

met = met./(ntrl1*ntrl2);