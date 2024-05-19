%
%  This Matlab function produces the mapping of 
%  integer partitions of k to its conjugate partition.
%
function ind = conjind(k)
   F = ones(k);
   for j=2:k
       F(:,j) = F(:,j-1);
       F(j,j) = F(j,j)+1;
       for i=j+1:k
           F(i,j) = F(i,j)+F(i-j,j);
       end
   end
   if k<=5
      ind = [F(k,k):-1:1]';
      return
   end
   s = ip_desc(k);
   m = F(k,k);
   l = sum(s>0,2);   % l is a vector of number of parts for each integer partition of k
   ind = ones(m,1)*m;
   ind(1:3) = [m m-1 m-2]';
   ind(m-2:m) = [3 2 1]';
   for i=4:m-3
       m1 = l(i);
       kappa = s(i,1:m1);
       lam = zeros(1,kappa(1));
       lam(1) = m1;
       n1 = k-m1;
       ind(i) = ind(i)-F(k,m1-1);
       for j=2:kappa(1)
           lam(j) = sum(kappa>=j);
           if lam(j)>1
              k2 = lam(j);
              ind(i) = ind(i)-F(n1,k2-1);
              n1 = n1-k2;
           else
              break
           end
       end
   end