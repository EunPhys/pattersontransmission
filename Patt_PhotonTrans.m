 %%%%%%%%%%%%%%%%%%%% Patterson transmission %%%%%%%%%%%%%%%%%%%%%
function Phi = diff(t)

    c = (3*10^8)/1.3460;
    % PER MM-1
   
mu_s = 1.25;
    mu_a = 0.0950;
    g = 0.9;
    d = 0.831;




    
    D = 1./(3*(mu_s+mu_a));
    z0 = 1./((1-g)*mu_s); 

            T1 = ((4*pi*D*c).^(-1/2)) .* t.^(-3/2) .* exp(-mu_a*c .*t);

            E1 = (((d-z0).^2) ./ (4*D*c.*t)); 
            T2 = (d-z0) * exp(-E1);

            E2 = (((d+z0).^2) ./ (4*D*c.*t)); 
            T3 = (d+z0) * exp(-E2);

            E3 = (((3*d-z0).^2) ./ (4*D*c.*t)); 
            T4 = (3*d-z0) * exp(-E3);

            E4 = (((3*d+z0).^2) ./ (4*D*c.*t)); 
            T5 = (3*d+z0) * exp(-E4);

    TRANS = T1.* ((T2 - T3 + T4 - T5));
    Phi = TRANS;
    
end