function [avs, phis, moment_arms]= calc_avs(W, D, L, f_contour, f_rho)
    phis = 1:180;
    moment_arms = zeros(size(phis));
    
    for phi = 1 : length(phis)
        moment_arms(phi) = moment_arm_func(phi, W, D, L, f_contour, f_rho);
    end
    
    crop_min = 10;
    crop_max = 170;
    moment_arms_cropped = moment_arms(10:170);
    
    % Get the AVS
    % The AVS satisfies two conditions: 
    % 1. The sign changes. 
    % 2. The derivative is negative
    found_indices = crop_min - 1 + ...
            find(diff(sign(moment_arms_cropped)) &...
                 diff(moment_arms_cropped) < 0);

    v_avs = phis(found_indices); 
    
    if isempty(v_avs)
        avs = 0;
        return
    end
    
    % Sometimes the first element of phis is mistakenly included
    if v_avs(1) == phis(1)
        v_avs = v_avs(2:end);
    end
    
    % Return the first (minimal) avs
    avs = v_avs(1);
    
end