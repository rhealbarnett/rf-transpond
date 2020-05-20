%--------------------------------%
% FDFD wave solver pml function
% rlbarnett c3149416
% 20/05/07
%--------------------------------%

function cpml = wave_pml(a_max,kap_max,np_pml,dz,eps0,eta0,cpdt,omega,m,ma,npts)

    sig_pml_s = zeros(1,np_pml);
    sig_pml_p = zeros(1,np_pml);
    a_pml = zeros(1,np_pml);
    kap_pml = zeros(1,np_pml);
    sz_pml_s = zeros(1,np_pml);
    sz_pml_p = zeros(1,np_pml);
    
    epss = (cpdt(1,1,1));
    epsp = (cpdt(3,3,1));
    
    sig_max_s = 25*(0.8*(m+1)/(eta0*dz*sqrt(epss)));
    sig_max_p = 25*(0.8*(m+1)/(eta0*dz*sqrt(epsp)));
    
%     Sr_d = 1.0;
%     Sr_dd = -1.0;
    
%     test_pml = zeros(1,np_pml);
    
    for kk=1:np_pml
    
        sig_pml_s(1,kk) = sig_max_s*((np_pml - kk)/(np_pml - 1.0))^m;
        sig_pml_p(1,kk) = sig_max_p*((np_pml - kk)/(np_pml - 1.0))^m;
        a_pml(1,kk) = a_max*((kk)/(np_pml - 1.0))^ma;
        kap_pml(1,kk) = 1.0 + (kap_max - 1.0)*((np_pml - kk)/(np_pml - 1.0))^m;
        sz_pml_s(1,kk) = kap_pml(1,kk) + sig_pml_s(1,kk)/(a_pml(1,kk) - 1i*omega*eps0);
        sz_pml_p(1,kk) = kap_pml(1,kk) + sig_pml_p(1,kk)/(a_pml(1,kk) - 1i*omega*eps0);
        
%         test_pml(1,kk) = 1.0 + (Sr_d + 1i*Sr_dd)*((np_pml - kk)/(np_pml - 1.0))^2;
        
    end
    
    cpml = ones(3,3,npts);
    
    for kk=1:npts
        
        cpml(1,1,1:np_pml) = sz_pml_s;
        cpml(1,1,end-np_pml+1:end) = fliplr(sz_pml_s);
        cpml(2,2,:) = cpml(1,1,:);
        cpml(3,3,1:np_pml) = sz_pml_p;
        cpml(3,3,end-np_pml+1:end) = fliplr(sz_pml_p);

%         cpml(1,1,1:np_pml) = test_pml;
%         cpml(1,1,end-np_pml+1:end) = fliplr(test_pml);
%         cpml(2,2,:) = cpml(1,1,:);
%         cpml(3,3,:) = cpml(1,1,:);
        
    end

end