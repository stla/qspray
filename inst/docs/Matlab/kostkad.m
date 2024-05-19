%
%  This Matlab function produces the Kostka matrix,
%  which is the matrix that maps Schur polynomials
%  to a linear combination of monomial symmetric
%  polynomials.  Optionally, it also returns the
%  inverse of the Kostka matrix. This version uses 
%  double precision.  Caveat: For C, the matrix 
%  is not fully accurate when k>=31.  
%
function [C,Ci] = kostkad(k)
   s = ip_desc(k);
   s1 = cumsum(s,2);
   m = size(s,1);
   l = sum(s>0,2);   % l is a vector of number of parts for each integer partition of k
   rho = zeros(m,1);
   rho(1) = k*(k-1);
   for i=2:m
       rho(i) = sum(s(i,1:l(i)).*(s(i,1:l(i))-1-2*(0:l(i)-1)));
   end
   d = ones(m,1);
   v = zeros(k,1);
   if m<3
      m1 = fix(m*(m-1)/4);
   else
      m1 = round(m^(7/5)+6);  % An approximate length of vector a
   end
%  a is a vector which is used to store a_{\lambda}(\mu) for (k)>=\mu>\lambda
%  alist is used to store the index for \mu
%  count is used to store the number of admissible \mu's for a given \lambda
   a = zeros(m1,1);
   alist = zeros(m1,1);
   count = zeros(m,1);
   for jj=2:m
       count(jj) = count(jj-1);
       m1 = l(jj);
       kappa = s(jj,1:m1);
%  Convert kappa=(kappa_1,kappa_2,...) to [1^{v_1}2^{v_2}...]
       for i=1:kappa(1)
           v(i) = sum(kappa==i);
       end
       for kk=1:jj-1
           if check2(kappa,s(kk,1:l(kk)))
              for q2=m1:-1:2
                  if kappa(q2)~=s(kk,q2)
                     break
                  end
              end
              t = s(kk,q2);
              while kappa(q2)==s(kk,q2-1)
                 q2 = q2-1;
              end
              t = kappa(q2)-t;
              for q1=q2-1:-1:1
                  if s(kk,q1)==kappa(q1)+t||(kappa(q1)~=s(kk,q1)&&kappa(q1-1)==s(kk,q1))
                     break
                  end
              end
              if kappa(q1)==kappa(q2)
                 w = v(kappa(q1));
                 w = w*(w-1)/2;
              else
                 w = v(kappa(q1))*v(kappa(q2));
              end
              count(jj) = count(jj)+1;
              a(count(jj)) = 2*(kappa(q1)-kappa(q2)+2*t)*w;
              alist(count(jj)) = kk;
           end
       end
   end
%  C is the transistion matrix from Schur to monomials
   C = eye(m);
   C(1,:)= 1;
   for ii=2:m-1
       for jj=ii+1:m
           if all(s1(ii,:)>=s1(jj,:))
              alist1 = alist(count(jj-1)+1:count(jj));
              ind = find(alist1>=ii);
              C(ii,jj) = round(C(ii,alist1(ind))*a(count(jj-1)+ind)/(rho(ii)-rho(jj)));
           end
       end
   end
   if nargout==2
%  When k>=29, taking the inverse of C is not fully accurate
      if k>=29
         Ci = kostkad_inv(k);
      else
         Ci = C;
         for j=2:m
             for i=1:j-1
                 Ci(i,j) = -Ci(i,i:j-1)*Ci(i:j-1,j);
             end
         end
      end
   end
   