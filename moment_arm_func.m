function moment_arm = moment_arm_func(phi, W, D, L, f_contour, f_rho)
    [com, cob, ~, ~] = boat_sim(phi, W, D, L, f_contour, f_rho, false);
    moment_arm = com(1) - cob(1);
end