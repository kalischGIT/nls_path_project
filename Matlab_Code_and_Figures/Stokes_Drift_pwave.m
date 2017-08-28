function ud = Stokes_Drift_pwave(eta1,s,Om,om,k0)

    q1 = -s*Om*eta1;
          
    qt = 1i*Om*eta1;
    
    % Compute phix
    R1 = q1;
    
    % Compute phiz    
    R1t = qt;
    ud = -2*k0/Om*(abs(R1)^2 + abs(R1t)^2 - om/Om*imag(R1.*conj(R1t)));
    
end

