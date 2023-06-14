
% [aos, aoa] = sphere_fibonacci_grid_points(1000);

function [aos, aoa] = sphere_fibonacci_grid_points ( ng )

%*****************************************************************************80
%
%% sphere_fibonacci_grid_points(): Fibonacci spiral gridpoints on a sphere.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    27 April 2015
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Richard Swinbank, James Purser,
%    Fibonacci grids: A novel approach to global modelling,
%    Quarterly Journal of the Royal Meteorological Society,
%    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.
%
%  Input:
%
%    integer NG, the number of points.
%
%  Output:
%
%    real XG(N,3), the grid points.
%
  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;

  i = ( - ( ng - 1 ) : 2 : ( ng - 1 ) )';
  theta = 2 * pi * i / phi;
  sphi = i / ng;
  cphi = sqrt ( ( ng + i ) .* ( ng - i ) ) / ng;

  xg = zeros ( ng, 3 );

  xg(1:ng,1) = cphi .* sin ( theta );
  xg(1:ng,2) = cphi .* cos ( theta );
  xg(1:ng,3) = sphi;
    
  [aos, aoa] = aos_aoa_converter(xg);
   
end


function [aos, aoa] = aos_aoa_converter(xg)
    % two circles
    % s_vector = [-1, 0, 0]  % xyz
    aos = zeros(1, length(xg));
    aoa = zeros(1, length(xg));

    for i = 1:length(xg)
        [azimuth, elevation, ~] = cart2sph(xg(i, 1), xg(i, 2), xg(i, 3));
        aos(i) = azimuth;
        aoa(i) = elevation;
    end
    
end
