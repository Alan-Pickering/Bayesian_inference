function gout = gfxn(gin, ulimits, choice)

s=size(ulimits,1);
if s>0
   Mx=ulimits(1);
end
if s==2
   Mn=ulimits(2);
end
if choice==-1
    gout = (Mx/(Mx+1)).*(Mx + 1 - gin);
elseif choice==-3
    gout = Mx + 1 - gin;
elseif choice==1
    gout = gin; %a simple linear function
elseif choice==2
    gout = gin.^2; %non-linear function used in Bogacz example
elseif choice==3
    gout=Mx./gin;
else
    %@@@ add other function choices
end

%preserving the min value of u
if s==2
  gout=max(gout,Mn);
  %testgo=gout<Mn;
  %disp(sum(testgo));
end

endfunction
