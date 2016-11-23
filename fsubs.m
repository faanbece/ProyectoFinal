function p = fsubs(x,c)
m = length(c);
n = length(x);
p = zeros(1,m-1);

for i = 1:m
   
   ci = c(i);
   
   if (i<m)
      cf = c(i+1)-1;
   else
       cf = n;
   end
   
   p(i) = prod(x(ci:cf));
end
end