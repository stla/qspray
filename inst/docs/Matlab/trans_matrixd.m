%
%   This Matlab program produces the transition matrix
%   between six different bases of symmetric polynomials.
%   The six bases are
%   e: elementary symmetric function
%   h: complete homogeneous symmetric function
%   m: monomial symmetric function
%   f: forgotten symmetric function
%   s: Schur function
%   p: power-sum symmetric function
%   Input: 
%   k: order of symmetric polynomials
%   b1: base 1 
%   b2: base 2
%   bases 1 and 2 are either 'e', 'h', 'm', 'f', 's' or 'p'
%   Output
%   M: transition matrix that expresses symmetric 
%      polynomials in base 1 in terms of linear 
%      combination of symmetric polynomials in base 2.
%   M1: When the transition matrix consists of rational
%       numbers (i.e., base 2 is p), we can specifiy two 
%       output arguments.  In that case, both M and M1 
%       are int64 matrices, and the transition matrix can 
%       be computed by using sym(M)./sym(M1)
%   Note 1: This version uses double precision.
%   Note 2: All the transition matrices have integer
%           elements, with the exception that when base 2 is p.
%   Reference: Symmetric Functions and Hall Polynomials,
%              Macdonald, 2nd edition, 1995, p.101
%
function [M,M1] = trans_matrixd(k,b1,b2)
   F = ones(k);
   for j=2:k
       F(:,j) = F(:,j-1);
       F(j,j) = F(j,j)+1;
       for i=j+1:k
           F(i,j) = F(i,j)+F(i-j,j);
       end
   end
   m = F(k,k);
   if b1==b2
      M = eye(m);
      return
   end
   if b1=='p'||b2=='p'
      f = cumprod(1:k);
      s = ip_desc(k);
      l = sum(s>0,2);   % l is a vector of number of parts for each integer partition of k
      zi = ones(m,1);
      zi1 = ones(m,1);
      for jj=1:m
          kappa = s(jj,1:l(jj));
%  Convert kappa=(kappa_1,kappa_2,...) to [1^{v_1}2^{v_2}...]
          for i=1:kappa(1)
              v = sum(kappa==i);
              if v>0
                 zi(jj) = zi(jj)*f(v)*i^v;
                 zi1(jj) = zi1(jj)*f(v);
              end
          end
      end
      T1 = mtop2(k);
   end   
   if b1=='e'
      if b2=='h'
         ind = conjind(k);
         [K,Ki] = kostkad(k);         
         M = triu(K'*Ki(:,ind)');
      elseif b2=='m'
         ind = conjind(k);
         K = kostkad(k);
         M = K'*K(ind,:);      
      elseif b2=='f'
         K = kostkad(k);
         M = K'*K;
      elseif b2=='s'
         ind = conjind(k);
         K = kostkad(k);
         M = K(ind,:)';
      elseif b2=='p'
         ei = (-1).^(k-l);
         if nargout<2
            ei = (-1).^(k-l);
            M = tril_inv(T1)'.*(zi1./(zi.*ei)');
         else
            M = int64(tril_inv(T1)'.*zi1);
            M1 = int64(ones(m,1)*(zi.*ei)');
            for i=1:m
                for j=i:m
                    if M(i,j)==0
                       M1(i,j) = 1;
                    else
                       g = gcd(M(i,j),M1(i,j))*sign(M1(i,j));
                       M(i,j) = M(i,j)/g;
                       M1(i,j) = M1(i,j)/g;
                    end
                end
            end
         end
      end
   elseif b1=='h'
      if b2=='e'
         ind = conjind(k);
         [K,Ki] = kostkad(k);
         M = triu(K'*Ki(:,ind)');
      elseif b2=='m'
         K = kostkad(k);
         M = K'*K;      
      elseif b2=='f'
         ind = conjind(k);
         K = kostkad(k);
         M = K'*K(ind,:);
      elseif b2=='s'
         M = kostkad(k)';
      elseif b2=='p'
         if nargout<2
            M = tril_inv(T1)'.*(zi1./zi');
         else
            M = int64(tril_inv(T1)'.*zi1);
            M1 = int64(ones(m,1)*zi');
            for i=1:m
                for j=i:m
                    if M(i,j)==0
                       M1(i,j) = 1;
                    else
                       g = gcd(M(i,j),M1(i,j));
                       M(i,j) = M(i,j)/g;
                       M1(i,j) = M1(i,j)/g;
                    end
                end
            end
         end
      end
   elseif b1=='m'
      if b2=='e'
         ind = conjind(k);
         Ki = kostkad_inv(k);
         M = triu(Ki(:,ind)*Ki(ind,:)');
         M = M(:,ind);
      elseif b2=='h'
         Ki = kostkad_inv(k);
         M = Ki*Ki';      
      elseif b2=='f'
         ind = conjind(k);
         [K,Ki] = kostkad(k);
         M = Ki*K(ind,:);
      elseif b2=='s'
         M = kostkad_inv(k);
      elseif b2=='p'
         if nargout<2
            M1 = T1./zi1;
         else
            M = int64(T1);
            M1 = int64(zi1*ones(1,m));
            for i=1:m
                for j=1:i
                    if M(i,j)==0
                       M1(i,j) = 1;
                    else
                       g = gcd(M(i,j),M1(i,j));
                       M(i,j) = M(i,j)/g;
                       M1(i,j) = M1(i,j)/g;
                    end
                end
            end
         end
      end
   elseif b1=='f'
      if b2=='e'
         Ki = kostkad_inv(k);
         M = Ki*Ki';
      elseif b2=='h'
         ind = conjind(k);
         Ki = kostkad_inv(k);
         M = triu(Ki(:,ind)*Ki(ind,:)');
         M = M(:,ind);
      elseif b2=='m'
         ind = conjind(k);
         [K,Ki] = kostkad(k);
         M = Ki*K(ind,:);
      elseif b2=='s'
         ind = conjind(k);
         Ki = kostkad_inv(k);
         M = Ki(:,ind);
      elseif b2=='p'
         ei = (-1).^(k-l);
         if nargout<2
            M = T1.*(1./zi1*ei');
         else
            M = int64(T1);
            M1 = int64(zi1*ei');
            for i=1:m
                for j=1:i
                    if M(i,j)==0
                       M1(i,j) = 1;
                    else
                       g = gcd(M(i,j),M1(i,j))*sign(M1(i,j));                       
                       M(i,j) = M(i,j)/g;
                       M1(i,j) = M1(i,j)/g;
                    end
                end
            end
         end
      end
   elseif b1=='s'
      if b2=='e'
         ind = conjind(k);
         Ki = kostkad_inv(k);
         M = Ki(:,ind)';
      elseif b2=='h'
         M = kostkad_inv(k)';
      elseif b2=='m'
         M = kostkad(k);
      elseif b2=='f'
         ind = conjind(k);
         K = kostkad(k);
         M = K(ind,:);
      elseif b2=='p'
         K = kostkad(k);
         if nargout<2
            M = K*diag(1./zi1)*T1;
         else
            zi2 = zi1(end);
            zi1 = zi2./zi1;
            zi2 = int64(zi2);
            M = int64(K*diag(zi1)*T1);
            M1 = zeros(m,'int64');
            for i=1:m
                for j=1:m
                    if M(i,j)==0
                       M1(i,j) = 1;
                    else
                       g = gcd(M(i,j),zi2);
                       M(i,j) = M(i,j)/g;
                       M1(i,j) = zi2/g;
                    end
                end
            end
         end
      end
   elseif b1=='p'
      if b2=='e'
         ei = (-1).^(k-l);
         M = round(T1'.*((ei.*zi)./zi1'));
      elseif b2=='h'
         M = round(T1'.*(zi./zi1'));
      elseif b2=='m'
         M = tril_inv(T1).*zi1';
      elseif b2=='f'
         ei = (-1).^(k-l);
         M = tril_inv(T1).*(ei*zi1');
      elseif b2=='s'
         Ki = kostkad_inv(k);
         M = tril_inv(T1)*diag(zi1)*Ki;
      end
   end