
function [phase_uw] = unwrapper(phase,npts)

    phase_uw = phase;
    for ii=2:npts
        diff = phase(ii) - phase(ii-1);
        if diff > (pi - 0.02*pi)
            phase_uw(ii:end) = phase_uw(ii:end) - pi;
        elseif diff < (-pi + 0.02*pi)
            phase_uw(ii:end) = phase_uw(ii:end) + pi;
        end
    end
end