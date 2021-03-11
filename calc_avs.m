function [avs, phis, moment_arms]= calc_avs(W, D, L, f_contour, f_rho)
    phis = 1:180;
    moment_arms = zeros(size(phis));
    
    for phi = 1 : length(phis)
        moment_arms(phi) = moment_arm_func(phi, W, D, L, f_contour, f_rho);
    end
    
    moment_arms_cropped = moment_arms(10:170);
    
    % Get the AVS
    % The AVS satisfies two conditions: 
    % 1. The sign changes. 
    % 2. The derivative is negative
    avs = phis(find(diff(sign(moment_arms_cropped)) ...
                & diff(moment_arms_cropped) < 0)); 
    
    if isempty(avs)
        avs = 0;
        return
    end
    
    % Sometimes the first element of phis is mistakenly included
    if avs(1) == phis(1)
        avs = avs(2:end);
    end
    
    if length(avs) > 1
        error("multiple AVS's found with params [%f, %f]", W, D)
    end
    
end