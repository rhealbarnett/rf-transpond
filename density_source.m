%-------------------------------------------------------------------------%
% Calculate density source term for RHS of continuity equation            %
%
% density_source(rate_coeff, fact, nxax, vxax, npts, neut_max, vx_new,
%                   n_new)
%
% INPUTS
% rate_coeff: ionsation rate coefficient 
%
% fact: percentage of the domain where the neutral density profile is
% non-zero for one boundary ONLY e.g. if you enter 0.05, this will be 5%
% of the total domain on either side (for 10% total across the whole
% domain). 
%
% nxax, vxax, npts: density and velocity grids, number of points.
%
% neut_max: maximum value for the neutral density profile.
%
% vx_new, n_new: velocity and density profiles used to calculate the
% neutral source term, fluxes at boundaries. 
%-------------------------------------------------------------------------%

function [n_source] = density_source(rate_coeff,fact,nxax,vxax,npts,neut_max,vx_new,n_new)

    %--
    % Get maximum grid value from the density grid. 
    xmax = max(vxax);
    xmin = min(vxax);
    dx = (xmax - xmin)/(npts-1);

    %--
    % Approximate size of non-zero portion of neutral profile.
    decay_loc = xmax - fact*xmax;
    
    %--
    % Determine the index for this density location. 
    a = find(nxax >= decay_loc);
    
    %--
    % Determine the number of points included in the non-zero region of the
    % neutral density. 
    decay_index = npts - a(1);
    
    %--
    % Initialise the enutral density array. 
    n_neut = zeros(1,npts);

    %-- 
    % Calculate shape of neutral profile. Approximated as half a cosine
    % wave.
    cosax = linspace(pi,2*pi,decay_index);
    n_neut(end-decay_index+1:end) = neut_max*((cos(cosax) + 1)/2);

    %--
    % Fill the remaining half of the array with the inner most calculated 
    % value for the neutral density.
    n_neut(end/2 + 1:end-decay_index+1) = neut_max*((cos(pi) + 1)/2);
    
    %--
    % For two metallic boundaries, mirror around the centre. 
    n_neut(1,1:end/2) = fliplr(n_neut(end/2 + 1:end));
    
    %--
    % Interpolate onto the variable density grid.
    n_neut = interp1(linspace(xmin-0.5*dx,xmax+0.5*dx,npts),n_neut,nxax,'linear');

    %--
    % Calculate the density source.
    n_source = n_new.*(n_neut)*(rate_coeff);

    %--
    % Interpolate source onto velocity grid to calculate the flux at the
    % boundary.
    source_avg = interp1(nxax,n_source,vxax,'linear');
    
    %--
    % Interpolate the densities onto the velocity grid. 
    n_avg = interp1(nxax,n_new,vxax);
    
    %--
    % Calculate the integral of the density source over the velocity grid.
    source_int = trapz(vxax,source_avg);

    %--
    % Calculate the fluxes at both boundaries. 
    rflux = n_avg(end)*vx_new(end);
    lflux = n_avg(1)*vx_new(1);
    
    %--
    % Calculate the multiplier to match density out = density in. 
    ns_mult = (rflux-lflux)/source_int;
    
    %--
    % Multiply the density source term by the constant ns_mult. 
    n_source = (n_source*ns_mult)*1.0e-2;
    
    %--
    % Calculate the new 'normalised' source term on the velocity grid. 
    nv_source = source_avg*ns_mult*1.0e-2;
    
    %--
    % Set the ghost point values for the density source to zero, as they
    % are not used in the calculation. 
    n_source(1,1) = 0.0; n_source(1,end) = 0.0;

end