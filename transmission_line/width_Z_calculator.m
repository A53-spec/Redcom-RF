function W = width_Z_calculator (Z0, e_r, d)
  %% CALCULATIONS
  A = Z0/60 *sqrt((e_r+1)/2) + (e_r-1)/(e_r+1)*(0.23 + 0.11/e_r);
  B = 377*pi/(2*Z0*sqrt(e_r));
  W_over_d = 8*exp(A)/(exp(2*A)-2);

  if W_over_d < 2
    W = W_over_d * d;
  elseif W_over_d > 2
    W =  d * (2/pi * (B - 1 -log(2*B-1) + (e_r-1)/(2*e_r)*(log(B-1) + 0.39 - 0.61/e_r)));
  end
endfunction
