function rhs = nonlinearity(lhs, a3)

% lhs is in spectral/frequency space.  Return to physical space 
    
    phys_lhs = ifft(lhs); 
    
% compute nonlinearity and then return to frequency space
    
    rhs = fft( a3.*phys_lhs.*real(phys_lhs.*conj(phys_lhs)) );
    