function [y_meas, x_next] = LPV_step_fast(x_curr, u_curr, Wf_eq, models)
    % Interpolates the LPV plant model to the current Wf_eq scheduling variable,
    % evaluates the output, and updates the state.
    % We do this natively here instead of calling Simulink for a massive speedup!

    % Find the closest operating point
    wf_values = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8];
    [~, k_plant] = min(abs(wf_values - Wf_eq));
    
    A_p = models.A_dt(:,:,k_plant); 
    B_p = models.B_dt(:,:,k_plant);
    C_p = models.C_dt(:,:,k_plant); 
    D_p = models.D_dt(:,:,k_plant);
    u0_p = models.offset.u(:,k_plant); 
    z0_p = models.offset.z(:,k_plant);
    
    idx_track = find(contains(lower(string(models.sys_eng.OutputName)), "n_h"), 1);
    if isempty(idx_track); idx_track = 2; end
    
    y_meas = z0_p(idx_track) + C_p(idx_track,:)*x_curr + D_p(idx_track,:)*(u_curr - u0_p);
    x_next = A_p*x_curr + B_p*(u_curr - u0_p);
end
