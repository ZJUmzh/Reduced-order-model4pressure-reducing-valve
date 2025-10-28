function An_new = poly_fit(ph_new,ph,An,nmodes)

for i=1:nmodes
    a_fit = fit(ph,An(:,i),'poly33','Normalize','on','Robust','Bisquare');
    An_new(:,i) = a_fit(ph_new);
end

end