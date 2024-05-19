%
%  This Matlab function produces the inverse Kostka matrix,
%  which is the matrix that maps monomial symmetric polynomials
%  to a linear combination of Schur polynomials.  The algorithm
%  is based on the recursion given by Egecioglu and Remmel (1990).
%  This version uses double precision.
%
function C = kostkad_inv(k)
%  When k<=28, it is faster to compute
%  K and take its inverse.
   if k<=28
      [~,C] = kostkad(k);
      return
   end
   F = ones(k);
   for j=2:k
       F(:,j) = F(:,j-1);
       F(j,j) = F(j,j)+1;
       for i=j+1:k
           F(i,j) = F(i,j)+F(i-j,j);
       end
   end
   C = cell(k,1);
   C{1} = 1;
   C{2} = [1 -1; 0 1];
   f = gammaln(2:k+1);
   for kk=3:k
       s = ip_desc(kk);
       s1 = cumsum(s,2);
       m = size(s,1);
       l = sum(s>0,2);   % l is a vector of number of parts for integer partition of kk
       C{kk} = eye(m);
%  Egecioglu and Remmel (1990), Corollary 1 (iii)
       for j=2:m
           if s(j,2)==1
              C{kk}(1,j) = (-1)^(l(j)-1);
           end
       end
       for i=2:m-1
           lam = s(i,1:l(i));
%  Egecioglu and Remmel (1990), Corollary 3
           if lam(1)==2
              l0 = sum(lam==2);
              k0 = sum(lam==1);
              for j=i+1:m
                  mu = s(j,1:l(j));
                  p0 = l0-sum(mu==2);
                  C{kk}(i,j) = (-1)^p0*nchoosek(k0+p0,p0);
              end
%  Egecioglu and Remmel (1990), Corollary 2, but with corrections
           elseif lam(2)==1   
              s0 = lam(1);
              C{kk}(i,m) = (-1)^(s0-1)*(kk-s0+1);
              for j=i+1:m-1
                  mu = s(j,1:l(j));
                  r0 = sum(mu(2:end)==1);
                  k0 = sum(mu(2:end)==2);
                  if mu(2)<=2&&kk-s0-r0<=k0&&k0<=kk-s0
                     C{kk}(i,j) = (-1)^(s0-mu(1)); 
                  end
              end
%  Egecioglu and Remmel (1990), Corollary 4 (i)     
           elseif l(i)==2
              if lam(1)==lam(2)
                 for j=i+1:m
                     mu = s(j,1:l(j));
                     if mu(2)==1
                        l0 = lam(1)-mu(1);
                        C{kk}(i,j) = (-1)^(kk-mu(1)-1);
                     else
                        r0 = sum(mu(3:end)==1);
                        k0 = sum(mu(3:end)==2);
                        if l(j)>2&&r0==mu(1)-mu(2)&&k0==lam(1)-mu(1)
                           C{kk}(i,j) = (-1)^(mu(1)-mu(2));
                        end
                     end
                 end
              else
%  Egecioglu and Remmel (1990), Corollary 4 (ii)     
                 for j=i+1:m
                     mu = s(j,1:l(j));
                     if mu(2)==1
                        if mu(1)<=lam(2)
                           C{kk}(i,j) = 2*(-1)^(kk-mu(1)-1);
                        else
                           C{kk}(i,j) = (-1)^(kk-mu(1)-1);
                        end
                     else
                        r0 = sum(mu(3:end)==1);
                        k0 = sum(mu(3:end)==2);
                        if mu(1)==lam(1)-k0&&mu(2)==lam(2)-r0-k0
                           C{kk}(i,j) = (-1)^r0;
                        elseif mu(1)==lam(1)-k0-r0-1&&mu(2)==lam(2)-k0+1
                           C{kk}(i,j) = (-1)^(r0+1);
                        elseif mu(1)==lam(2)-k0&&mu(2)==lam(1)-r0-k0
                           C{kk}(i,j) = (-1)^r0;
                        end
                     end
                end
              end
           else             
              cc0 = f(l(i)-1);
              for jj=1:lam(1)
                  v = sum(lam==jj);
                  if v>1
                     cc0 = cc0-f(v);
                  end
              end
              cc0 = exp(cc0);
%  Egecioglu and Remmel (1990), Corollary 1 (iv)
              C{kk}(i,m) = (-1)^(kk-l(i))*round(cc0*l(i));
              for j=i+1:m
                  if l(j)>=l(i)
                     mu = s(j,1:l(j));
%  Duan (2003), Lemma 6
                     if mu(2)==1
                        cc1 = (-1)^(l(j)-l(i))*sum(lam>=mu(1));
                        C{kk}(i,j) = round(cc1*cc0);
%  Duan (2003), Corollary 2
                     elseif lam(1)==mu(1)
                        kk1 = kk-lam(1);
                        ii = 2;
                        while lam(ii)==mu(ii)
                           kk1 = kk1-lam(ii);
                           ii = ii+1;
                        end
                        C{kk}(i,j) = C{kk1}(part_ind(lam(ii:end)),part_ind(mu(ii:end)));
%  A new reduction formula
                     elseif l(i)==l(j)&&lam(end)==mu(end)
                        kk1 = kk-lam(end);
                        ii = l(j)-1;
                        while lam(ii)==mu(ii)
                           kk1 = kk1-lam(ii);
                           ii = ii-1;
                        end
                        C{kk}(i,j) = C{kk1}(part_ind(lam(1:ii)),part_ind(mu(1:ii)));
%  Duan (2003), Theorem 3, but only for the case that \lambda_1=\mu_1+1
                     elseif lam(1)==mu(1)+1  
                        kk1 = kk-lam(1);
                        ind1 = part_ind(lam(2:l(i)));
                        mu1 = mu(2:l(j));
                        for ii=find(diff([mu1 0])<0)   
                            omega = mu1;
                            omega(ii) = omega(ii)-1;
                            C{kk}(i,j) = C{kk}(i,j)-C{kk1}(ind1,part_ind(omega));
                        end     
%  Find the second largest element of lam and see if it is equal to mu_1
                        for ii=2:l(i)
                            if lam(ii)<lam(1)
                               if lam(ii)==mu(1)
                                  kk1 = kk-lam(ii);
                                  ind1 = part_ind(lam([1:ii-1 ii+1:l(i)]));
                                  C{kk}(i,j) = C{kk}(i,j)+C{kk1}(ind1,part_ind(mu1));
                               end
                               break
                            end
                        end       
                     elseif all(s1(i,1:l(i))>=s1(j,1:l(i)))   % check if lam dominates mu                
%  Recursion based on Egecioglu and Remmel (1990)
                        for ii=1:l(j)
                            mui = mu(ii)+l(j)-ii;
                            jj = find(lam==mui,1);
                            if ~isempty(jj)
                               kk1 = kk-mui;
                               lam1 = lam;
                               lam1(jj) = [];
                               mu1 = mu;
                               mu1(ii) = [];
                               mu1(ii:end) = mu1(ii:end)-1;
                               C{kk}(i,j) = C{kk}(i,j)+(-1)^(l(j)-ii)*C{kk1}(part_ind(lam1),part_ind(mu1));
                            end
                        end
                     end
                  end
              end
           end
       end
   end
   C = C{k};
function y = part_ind(nu)
%
%  This function takes nu, an integer partition of k, and returns the location of 
%  nu within the list of integer partitions with |kappa|=k.
%
   n1 = kk1;
   y = F(n1,n1);
   for i2=1:length(nu)
       k1 = nu(i2);
       if k1>1
          y = y-F(n1,k1-1);
          n1 = n1-k1;
       else
          break
       end
   end
end
end