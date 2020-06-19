

function ans = coupled_rf_transp(hh,ww)

%     kx = (10 + (hh-1)*5)*1i;
    kx = hh*1i;
    Nmax = ww;
   
    if hh==10
        source_mult = 1.0e5;
    elseif hh==15
        source_mult = 1.3e5;
    elseif hh==20
        source_mult = 1.6e5;
    elseif hh==25
        source_mult = 2.0e5;
    end

    transport_1d

    [density_pert,mean_pert] = den_pert(n_init, n_new, nxax, 0);

    transport.density_pert = density_pert;
    transport.mean_pert = mean_pert;

    filename = strcat('/home/c3149416/coupled_results/coupled_',num2str(ii),'_',num2str(source_mult),'_',num2str(kx(1)),...
        '_',num2str(Nmax),'.mat');
    save(filename,'-struct','transport','-append');

end

%'coupled_transport_',num2str(ii),'.mat'
