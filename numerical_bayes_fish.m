%numerical_bayes_fish.m
%v1; 28.03.23
%
%Bogacz (2017) exercise 1 extended to fishing
%
%This does a numerical integration to implement an "exact" Bayesian
%estimation.
%v is the value being inferred/estimated (the number of fish that can be caught in a specific
%length of time at a place on a specific lake or river).
%From past fishing experince you have a prior distribution of probabilities for v = p(v)
%But you also judge how long in minutes will be the wait until your first bite. 
%This is parameter u (with a certain amount of noise in your estimate of u).
%You can infer a conditional posterior distribution for v given u, p(v|u) using Bayes'
%theorem.
%
%p(v|u) = [p(v)*p(u|v)]/p(u)
%
%where p(u|v) is the likelihood of u given possible values of v
%i.e the likelihood that the fishing conditions are estimated to be of
%observed value u for each potential number of fish available.
%This assumes a specific function relating v to u, and error variance in the estimate u.
%And p(u) is the evidence (or marginal likelihood) and acts as a normalizsation
%term to ensure p(u) is a proper density (values sum to 1). 
%p(u) is computed by integrating p(v)*p(u|v) over all values of v.
%Here we do the integration numerically over a range of values of v in small
%steps, DV.
%

clc;
clear variables;

%values for integration
%for Bogacz model 1 set bog_or_fish = 1; =2 for fish model
bog_or_fish=2;
if bog_or_fish==1
    %food size example from Bogacz 2017
    MINV = 0.01; % minimum value of v for which posterior computed
    DV = 0.01; % interval between values of v for which posterior found
    MAXV = 6; %can vary for plotting; maximum value of v for which posterior computed 
elseif bog_or_fish==2
    %fishing example
    MAXV=45;
    MINV=0;
    DV=0.01;
end

vrange = MINV:DV:MAXV; %range of possible v values
if bog_or_fish==1
    u = 2;  %observed light intensity
elseif bog_or_fish==2
    u = 1.5; %observed wait to first bite in minutes, try diff values in range 1-45
end

fxntype=-1; %this chooses the fxn g(v) that drives the value of u
if bog_or_fish==1
    fxntype=2; %this is the function that Bogacz used
end
if fxntype==-1 
    %linearly decreasing function g(v)=v
    %needs a bit of rescaling
    meanu_fromv= (MAXV/(MAXV+1)).*(MAXV + 1 - gfxn(vrange,fxntype));
elseif fxntype==1
    %linearly increasing function g(v)=v
    meanu_fromv= gfxn(vrange,fxntype);
elseif fxntype==2
    %non-linearly increasing function, used by Bogacz
    %g(v)=v*v;
    meanu_fromv=gfxn(vrange,fxntype);
elseif fxntype==3
    %a non-linearly decreasing function based around u=1/v; 
    vrange=MINV:DV:MAXV;
    meanu_fromv=MAXV.*gfxn(vrange+1,fxntype); %+1, computed based on 1-46 to avoid infinite values
end

fxnplot=1; %1=plot the function first, 0 do not
if fxnplot==1
    figure;
    plot(vrange,meanu_fromv,'-b','LineWidth', 4);
    ax = gca; %get number of current axis
    ax.FontSize = 12; %set axis font size
    if bog_or_fish==1
        xlabel('v=Radius of spherical food pellet','FontSize',16); %make label font bigger
        ylabel('u=Estimate of light intensity (arbitray units)','FontSize',16);
    elseif bog_or_fish==2
        xlabel('v=Number of fish available to be caught','FontSize',16); %make label font bigger
        ylabel('u=Estimate of time to first bite (mins)','FontSize',16);
    end
end

pdftype=1;
if pdftype==1
    %normally distribution likelihood for u
    if bog_or_fish==1
        sigma_u = 1; % standard deviation of estimate of light intensity
    elseif bog_or_fish==2
        sigma_u = 1; % standard deviation of noise affecting time to first bite
    end
    likelihd=normpdf(u,meanu_fromv, sigma_u);
    %normal PDF for v and 
    if bog_or_fish==1
        v_p = 3; % mean of prior distribution of size (radius) of food pellet
        sigma_p = 1; % standard deviation of prior
    elseif bog_or_fish==2
        v_p = 10; % mean of prior distribution of roach available to be caught in 45 minutes
        sigma_p = 3; % standard deviation of prior
    end
    vpriorpdf=normpdf(vrange,v_p,sigma_p);
elseif pdftype==2
    %normally distributed likelihood for u
    sigma_u = 1; % standard deviation of noise time to first bite
    likelihd=normpdf(u,meanu_fromv, sigma_u); 
    %and compute this funky prior distribution for v
    v_p = 0.5*MAXV;
    vpriorpdf=sqrt(vrange.*(MAXV-vrange));
elseif pdftype==3
    %add others perhaps
end

%normalise to turn the pdfs into proper densities
priornorm=sum(vpriorpdf);
vpriorpdf=vpriorpdf./priornorm;
%@@we need to work out sd of prior here when non-standard pdf

%now we do the integration
%computing the marginal likelihoods
if pdftype==1 || pdftype==2
    numerator =  vpriorpdf.* likelihd;
else
    %not sure there is an else here?!
end
normalization = sum(numerator); %(numerator.*DV); %@@@changed from Bogacz to be more intutitive IMO
prior_p= vpriorpdf; %prior probability p(v)
post_p = numerator / normalization ; %this is the posterior probability p(v|u) by Bayes 
figure;
plot (vrange , prior_p, 'b','LineWidth', 4);
hold on;
plot (vrange , post_p, 'r','LineWidth', 4);
ax = gca; %get number of current axis
ax.FontSize = 12; %set axis font size
if bog_or_fish==1
    xlabel ('Size of food item (v)','FontSize',16);
elseif bog_or_fish==2
    xlabel ('No. of fish available to be caught (v)','FontSize',16);
end
ylabel ('Probability density');
legend('Prior p(v)','Post. p(v|u)','Location','best')
%work out how to convert density back into actual summary numbers of fish
%the posterior probability is a PDF as the values sum to 1
%because the values of post_p sum to 1 then the mean is just the sum of the values
mean_v=sum(post_p.*vrange); 
[max_post, max_idx]=max(post_p); %max_idx is the index of the array for max posterior
mode_v=vrange(max_idx); %work out the v value for the max posterior index
%what is the weighted sd of the distribution?
sumsqdev=sum(post_p.*(vrange - mean_v).*(vrange - mean_v)); %sum of the squared deviations
M=sum(post_p~=0); %number of non-zero weights
sumwts=sum(post_p);
sd_v=sqrt((M*sumsqdev)/((M-1)*sumwts));
%@@@maybe add details of prior mean and sd here and point estimate based on the
%function using the observed value
%pointval_v = (MAXV-1).*gfxn(u,fxntype);
%disp([v_p,pointval_v]);
if bog_or_fish==1
    disp(['The mean of posterior distribution of size of food pellet = ' num2str(mean_v)]);
    disp(['The mode of posterior distribution of size of food pellet = ' num2str(mode_v)]);
    disp(['The s.d. of posterior distribution of size of food pellet = ' num2str(sd_v)]);
elseif bog_or_fish==2
    disp(['The mean of posterior distribution of # of fish available = ' num2str(mean_v)]);
    disp(['The mode of posterior distribution of # of fish available = ' num2str(mode_v)]);
    disp(['The s.d. of posterior distribution of # of fish available = ' num2str(sd_v)]);
end

%@@@compute the surprisal of the distribution cf the prior

function gout = gfxn(gin, choice)

if choice==-1 || choice==1
    gout = gin; %a simple linear function
elseif choice==2
    gout = gin.^2; %non-linear function used in Bogacz example
elseif choice==3
    gout=1./gin;
else
    %@@@ add other function choices
end

end
