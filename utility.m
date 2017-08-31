function U=utility(c)

global sigma;

if sigma==1;
    U=log(c);
else
   U=(c.^(1-sigma)-1)/(1-sigma);
end

