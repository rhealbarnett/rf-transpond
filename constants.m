%-------------------------------------------------------%
% physical constants                                    %
%-------------------------------------------------------%

% e = 1.6022e-19;
% mp = 1.67e-27;
% me = 9.11e-31;
% kb = 1.38e-23;
% mu0 = 4.0*pi*1.0e-7;
% eps0 = 8.85e-12;
% c0 = 1.0/sqrt(eps0*mu0);
function const = constants()
    const.e = 1.6022e-19;
    const.mp = 1.67e-27;
    const.me = 9.11e-31;
    const.kb = 1.38e-23;
    const.mu0 = 4.0*pi*1.0e-7;
    const.eps0 = 8.85e-12;
    const.c0 = 1.0/sqrt(const.eps0*const.mu0);
end