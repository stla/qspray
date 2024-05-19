function S = ip_desc(n)
%
% Adapted from "Algorithm ZS1" in
% ANTOINE ZOGHBI and IVAN STOJMENOVIC (1998), "FAST ALGORITHMS FOR 
% GENERATING INTEGER PARTITIONS", International Journal of Computer 
% Mathematics, Volume 70, Issue 2, pages 319-332
% DOI: 10.1080/00207169808804755
%
   x = ones(1,n);
   x(1) = n;
   m = 1;
   h = 1;
   S = [n zeros(1,n-m)];
   while x(1)~=1
      if x(h)==2 
         m = m+1;
         x(h) = 1;
         h = h-1;
      else
         r = x(h)-1;
         t = m-h+1;
         x(h) = r;
         while t>=r
            h = h+1;
            x(h) = r;
            t = t-r;
         end
         if t==0
            m = h;
         else
            m = h+1;
            if t>1
               h = h+1;
               x(h) = t;
            end
         end
      end
      S = cat(1,S,[x(1:m) zeros(1,n-m)]);
   end
