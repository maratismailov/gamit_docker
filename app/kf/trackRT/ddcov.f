% Now do analytic
covd = -ones(np);
for ks = 1:ns                       ! Number of satellites
    for kl = 1:nl                   ! Number of stations
        k = (kl-1)*ns + ks;         ! Parameter number
        r1 = (kl-1)*ns + nr(kl);    ! New ref parameter number
        for ms = 1: ns                     ! Number of satellites
            for ml = 1:nl                  ! Number of stations
                m = (ml-1)*ns + ms;        ! Parameter number
                r2 = (ml-1)*ns + nr(ml);   ! New ref parameter number at this site 
                if kl == ml                ! Same site block 
                    if k == r1 && m == r1 
                        covd(k,m) = covo(r1,r1);
                    elseif k == r1
                        covd(k,m) = (covo(r1,r1)-covo(m,r1));
                    elseif m == r1
                        covd(k,m) = (covo(r1,r1)-covo(r1,k));                        
                    else
                        covd(k,m) = (covo(m,k)-covo(r1,m))-(covo(k,r1)-covo(r1,r1));
                    end
                else                       ! Upper diagonal site block (lower generated at same time).
                    if m > k 
                        if k == r1 && m == r2
                            covd(k,m) = covo(k,m);
                            covd(m,k) = covo(k,m);
                        elseif k == r1 && m ~= r2
                            covd(k,m) = covo(r1,r2)-covo(k,m);
                            covd(m,k) = covd(k,m);
                        elseif m == r2
                            covd(k,m) = covo(r1,r2)-covo(k,m);
                            covd(m,k) = covd(k,m);
                        else
                            covd(k,m) = (covo(k,m)-covo(r1,m))-(covo(r2,k)-covo(r1,r2));
                            covd(m,k) = covd(k,m);
                        end
                    end

                end

            end
        end
    end
end
