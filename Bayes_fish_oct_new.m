%Bayes_fish_oct_new.m
%v1; 12.04.23
%version for Octave or Matlab
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
%Important information
%*********************
%i)   This is the Octave compatible version, which runs under Octave or Matlab
%     Under Octave needs statistics package
%     Search Matlab in code for small edits to optimise for Matlab
%ii)  This requires the bespoke function gfxn.m
%iii) search @@@ for features which could be added
%
%Usage Notes
%***********
%you can select values of variables below to run variations:-
%bog_or_fish chooses the Bogacz original or fishing
%fxntype controls the nature of the function relating u to v; i.e. u=g(v)
%fxnplot controls whether the function is plotted
%pdftype controls the likelihood used for u|v and pdf used for v
%
%and you can vary the values of the model paremeters etc
%MINV and MAXV control the range of v values over which integration is done
%u is the estimate of light intensity (Bogacz) or time to first bite (fishing)
%sigma_u is the s.d. of u
%Maxu controls the maximum value of u
%v is the thing you are trying to estimate (food size, number of fish)
%v_p is the mean of the prior for v
%sigma_p is the s.d. of the prior for v

clc;
clear variables;

%values for integration
%for Bogacz model 1 set bog_or_fish = 1; =2 for fish model
bog_or_fish=1; %2
if bog_or_fish==1
    %food size example from Bogacz 2017
    MINV = 0.01; % minimum value of v for which posterior computed
    DV = 0.01; % interval between values of v for which posterior found
    %MAXV=maximum value of v for which posterior computed
    %can vary for plotting;
    MAXV = 6;
    %Maxu is not relevant for the fxn used by Bogacz
    Maxu=MAXV.^2;
    Minu=MINV.^2;
elseif bog_or_fish==2
    %fishing example
    Minu=1;
    Maxu=45;
end
urange=[Maxu; Minu];

if bog_or_fish==1
    u = 2;  %observed light intensity
elseif bog_or_fish==2
    u = 1.5; %observed wait to first bite in minutes, try diff values in range 1-45
end

%for the fishing example use fxntype -3 or 3
fxntype=-3; %3; %this chooses the fxn g(v) that drives the value of u
if bog_or_fish==1
    fxntype=2; %this is the function that Bogacz used
end

if fxntype==3
    %these values work well for the non-linearly decreasing function
    MAXV=20;
    MINV=0;
    DV=0.01;
elseif fxntype==-3
    %these values work well for the linearly decreasing function
    MAXV=50;
    MINV=30;
    DV=0.01;
end;

vrange = MINV:DV:MAXV; %range of possible v values

%could easily remove all the if elseif below
%and just write
%meanu_fromv= gfxn(vrange, urange, fxntype);
%for all cases, although urange not relevant/undefined for some
if fxntype==-1
    %linearly decreasing function g(v) = (M/(M+1))*(M + 1 -v)
    meanu_fromv= gfxn(vrange, urange, fxntype);
elseif fxntype==-3
    %linearly decreasing function g(v) = M + 1 -v
    meanu_fromv= gfxn(vrange, urange, fxntype);
elseif fxntype==1
    %linearly increasing function g(v)=v
    meanu_fromv= gfxn(vrange, [ ], fxntype);
elseif fxntype==2
    %non-linearly increasing function, used by Bogacz
    %g(v)=v*v;
    %could pass urange to the function
    meanu_fromv=gfxn(vrange, [ ], fxntype);
elseif fxntype==3
    %a non-linearly decreasing function u=M/(v+1);
    %+1 used to avoid infinite values when v=0
    meanu_fromv=gfxn(vrange, urange, fxntype);
end

fxnplot=1; %1=plot the function first, 0 do not
if fxnplot==1
    figure;
    plot(vrange,meanu_fromv,'-b','LineWidth', 4);
    ax = gca; %get number of current axis
    %the following line works in Matlab but not in Octave
    %ax.FontSize = 12; %set axis font size
    if bog_or_fish==1
        xlabel('v=Radius of spherical food pellet','FontSize',16); %make label font bigger
        ylabel('u=Estimate of light intensity (arbitrary units)','FontSize',16);
    elseif bog_or_fish==2
        xlabel('v=Number of fish available (arbitrary units)','FontSize',16); %make label font bigger
        ylabel('u=Estimate of time to first bite (mins)','FontSize',16);
    end
end

pdftype=1; %^stick to using 1 at the moment
if pdftype==1
    %normally distribution likelihood for u
    if bog_or_fish==1
        sigma_u = 1; % standard deviation of estimate of light intensity
    elseif bog_or_fish==2
        sigma_u = 1; % standard deviation of noise affecting time to first bite
    end
    likelihd=normpdf(u,meanu_fromv, sigma_u);
    if bog_or_fish==1
        v_p = 3; % mean of prior distribution of size (radius) of food pellet
        sigma_p = 1; % standard deviation of prior
    elseif bog_or_fish==2
        if fxntype==-3
            %use these values for the linear decreasing model
            v_p = 41.5; %mean of prior distribution of roach available to be caught in 45 minutes
            sigma_p = 2; %standard deviation of prior
        elseif fxntype==3
            %use these values for the non-linear decreaing model
            v_p = 9;
            sigma_p = 3;
        else
            %these are the values used by Bogacz
            v_p = 3;
            sigma_p = 1;
        end
    end
    %normal PDF for v
    vpriorpdf=normpdf(vrange,v_p,sigma_p);
else
    %add others perhaps
end

%normalise to turn the pdfs into proper densities
priornorm=sum(vpriorpdf);
vpriorpdf=vpriorpdf./priornorm;
%@@@we need to work out sd of prior here when non-standard pdf

%what are the u values expected based on the prior dist from v
%95% CLs and mean
pre_v_vals=[v_p-1.96*sigma_p; v_p; v_p+1.96*sigma_p];
pre_u_vals=gfxn(pre_v_vals,urange,fxntype);
disp('Based on the prior distribution of v, we predict that ')
disp(['the expected value of u = ' num2str(pre_u_vals(2))]);
disp(['and it will lie between 95% CLs of ' num2str(pre_u_vals(1)) ' and ' num2str(pre_u_vals(3))]);
disp(' ');
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
%the following line works in Matlab but not in Octave
%ax.FontSize = 12; %set axis font size
if bog_or_fish==1
    xlabel ('Size of food item (v; arbitrary units)','FontSize',16);
elseif bog_or_fish==2
    xlabel ('No. of fish available (v; arbitrary units)','FontSize',16);
end
ylabel ('Probability density','FontSize',16);
legend('Prior p(v)','Post. p(v|u)','Location','northeast')
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
post_v_vals=[mean_v-1.96*sd_v; mean_v; mean_v+1.96*sd_v];
post_u_vals=gfxn(post_v_vals,urange,fxntype);
if bog_or_fish==1
    disp(['The mean of posterior distribution of size of food pellet = ' num2str(mean_v)]);
    disp(['The mode of posterior distribution of size of food pellet = ' num2str(mode_v)]);
    disp(['The s.d. of posterior distribution of size of food pellet = ' num2str(sd_v)]);
elseif bog_or_fish==2
    disp(['The mean of posterior distribution of # of fish available = ' num2str(mean_v)]);
    disp(['The mode of posterior distribution of # of fish available = ' num2str(mode_v)]);
    disp(['The s.d. of posterior distribution of # of fish available = ' num2str(sd_v)]);
end
disp(' ');
disp('Based on the posterior distribution of v, we predict that ')
disp(['the expected value of u = ' num2str(post_u_vals(2))]);
disp(['and it will lie between 95% CLs of ' num2str(post_u_vals(1)) ' and ' num2str(post_u_vals(3))]);


