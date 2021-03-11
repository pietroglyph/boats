function [center_of_mass, center_of_buoyancy, mass_boat, mass_water] = boat_sim(heel_angle, W, D, L, contour_func, rho_func, should_draw)
    samples_per = 100;
    rho_water = 1000; % kg/m^3
    
    % Check if optional argument exists
    if ~exist("should_draw", "var")
        should_draw = true;
    end
    
    %% Define boat
    % Base meshgrid
    [x, y] = meshgrid(linspace(-W, W, samples_per), linspace(0, D, samples_per));
    p = [x(:), y(:)];
    
    x_vals = unique(p(:, 1));
    y_vals = unique(p(:, 2));
    dx = abs(x_vals(2) - x_vals(1));
    dy = abs(y_vals(2) - y_vals(1));
    dA = dx * dy;

    % Roate boat
    p_rotated = p * [cosd(heel_angle) -sind(heel_angle); sind(heel_angle) cosd(heel_angle)]';
    
    % Check if in boat
    inside_boat = contour_func(p(:, 1), p(:, 2));
    
    %% Integrate conditional masses
    function [mass_difference, masses_boat, masses_water, inside_boat_under_water] = find_masses(d)
        inside_boat_under_water = inside_boat & p_rotated(:, 2) <= min(p_rotated(inside_boat, 2)) + d;
        
        masses_boat = inside_boat .* rho_func(p(:, 2)) * dA * L;
        masses_water = inside_boat_under_water * rho_water * dA * L;
        
        mass_difference = masses_boat - masses_water;
    end

    %% Solve for water level
    % Define possible range of depths for fzero
    y_min = min(p(:, 2));
    y_max = max(p(:, 2));
    
    % Solve for buoyant water level
    try
        water_level = fzero(@(d)( sum(find_masses(d)) ), [y_min, y_max]);
    catch exception
        if strcmp(exception.identifier, 'MATLAB:fzero:ValuesAtEndPtsSameSign')
            % TODO: Check this BEFORE the fzero try/catch to save time
            if sum(find_masses(y_max)) > 0
                %fprintf("Boatsim saturating waterlevel heavy phi=%f, [W, D] = [%f, %f]\n", heel_angle, W, D);
                water_level = y_max;
            end
            
        else
            % Otherwise throw error
            fprintf("Boatsim error at phi=%f, [W, D] = [%f, %f]\n", heel_angle, W, D);
            throw(exception)
        end  
    end
    [~, masses_boat, masses_water, inside_boat_under_water] = find_masses(water_level);
    
    
    %% Draw output
    if should_draw
      plot(p_rotated(inside_boat, 1), p_rotated(inside_boat, 2), "ro", p_rotated(inside_boat_under_water, 1), p_rotated(inside_boat_under_water, 2), "bo");
      axis equal;
    end
    
    %% Set return outputs
    mass_boat = sum(masses_boat);
    mass_water = sum(masses_water);
    center_of_mass = (masses_boat' * p_rotated) ./ mass_boat;
    center_of_buoyancy = (masses_water' * p_rotated) ./ mass_water;
end