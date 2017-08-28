function rhs = phi_eval_pwave_ho_1d_theta(t,lhs,s,A,anl,Om,om,n0,n2,k0,ep)
    
    Omt = Om + anl*A^2*ep^2;
    
    rhs = Omt + ep^2*k0*(om*n0-2*Om*k0)*A^2 ...
          + ep*k0*(om - s*Om)*(2*A*cos(lhs)) ...
          + ep^2*k0*((om-2*s*Om)*n2- abs(k0)*om)*(2*A^2*cos(2*lhs));
    
end

