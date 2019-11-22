function [ans] = pond_source(me,mp,omega,q,Efield,dx)
    pond_const = (1.0/4.0)*((q^2)/(me*omega^2));
    for ii = 2:length(dx)
        Ediff(1,ii-1) = -(dx(1,ii)/(dx(1,ii-1)*(dx(1,ii) + dx(1,ii-1))))*Efield(1,ii-1) +...
            ((dx(1,ii) - dx(1,ii-1))/(dx(1,ii)*dx(1,ii-1)))*Efield(1,ii) +...
            (dx(1,ii-1)/(dx(1,ii)*(dx(1,ii) + dx(1,ii-1))))*Efield(1,ii+1);
    end
    ans = (1.0/mp)*pond_const*Ediff;
end