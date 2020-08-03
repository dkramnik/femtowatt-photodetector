function [ z_p ] = prl( z_1, z_2 )
% Parallel impedance of two transfer functions.
    
    z_p = z_1 * z_2 / (z_1 + z_2);

end

