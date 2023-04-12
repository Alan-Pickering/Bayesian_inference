function gout = gfxn(gin, M, choice)

if choice==-1
    %meanu_fromv= (MAXV/(MAXV+1)).*(MAXV + 1 - gfxn(vrange,fxntype));
    gout = (M/(M+1)).*(M + 1 - gin);
elseif choice==1
    gout = gin; %a simple linear function
elseif choice==2
    gout = gin.^2; %non-linear function used in Bogacz example
elseif choice==3
    gout=M./gin;
else
    %@@@ add other function choices
end

endfunction
