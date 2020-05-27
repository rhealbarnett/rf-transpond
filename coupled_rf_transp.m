

function ans = coupled_rf_transp(hh)

%     kx = (10 + (hh-1)*5)*1i;
    kx = hh*1i;

    transport_1d

    [density_pert,mean_pert] = den_pert(n_init, n_new, nxax, 0);

    transport.density_pert = density_pert;
    transport.mean_pert = mean_pert;

    filename = strcat('/Volumes/DATA/LAPD/matlab/results_test/coupled_',num2str(source_mult),'_',num2str(kx),'_',...
        num2str(ky),'_',num2str(Nmax),'.mat');
    save(filename,'-struct','transport');

end

%'coupled_transport_',num2str(ii),'.mat'