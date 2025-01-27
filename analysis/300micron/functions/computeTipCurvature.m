function rk = computeTipCurvature(p,rr1,opt_deg);
    px = p(opt_deg); pxx = 0;
    for m=1:opt_deg-1
        px = px+(opt_deg-m+1)*p(m)*rr1.^(opt_deg-m);
        pxx = pxx+(opt_deg-m+1)*(opt_deg-m)*p(m)*rr1.^(opt_deg-m-1);
    end
    % curvature
    k0 = max(abs(pxx)./(1+px.^2).^(3/2)); 
    % radius of curvature
    rk = 1/k0; 
end