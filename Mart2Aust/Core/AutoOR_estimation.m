function [final_OR,final_halfwidth,OR_metadata] = AutoOR_estimation(ebsd,options)
% Given a Martensite EBSD scan, estimate the Orientation Relationship

% Pre-Flight stuff to skip the OR determination if unnecessary
assert(ismatrix(options.OR_ksi) && length (options.OR_ksi) == 3,...
    'options.OR_ksi must be a length 3 matrix of positive angles')
assert(isfloat(options.OR_noise) && length (options.OR_noise) == 1,...
    'options.OR_noise must be a length 1 float')
options.OR_noise = abs(options.OR_noise);
if options.OR_ksi(1) <=0 || options.OR_noise <=0
    disp("Calculating the Orientation Relationship ...")
     [final_OR,final_halfwidth,OR_metadata] = Calculate_OR(ebsd,options);
else
    disp("skipping OR calculation, using saved values for OR and Martensite background noise")
    final_OR = options.OR_ksi;
    final_halfwidth = options.OR_noise;
    OR_metadata = [];
end
end


    function [final_OR,final_halfwidth,OR_metadata] = Calculate_OR(ebsd,options)
% Given a Martensite EBSD scan, estimate the Orientation Relationship
% This is done in steps:
%   1) use MTEX + graph cut to find a possible Aus Grain (Aus_Grn_guess)
%   2) attempt to find a reasonable Aus OR within allowed tolerance
%   3) if fails, repeat using scan that wasn't in previous assumed grain

tic;
%% Useful variables for calculations
CS_HT = ebsd.CSList{1}; % equivilant to CS_R
CS_LT = ebsd.CSList{2}; % equivilant to CS_T
SS = ebsd.opt.SS;

%% search the low temp phase (Phase 2) for acceptably close PAG
% for first iteration, make the entire scan searchable
searchable_ebsd = ebsd;
for iteration = 1:10
    % NOTE TO FUTURE USERS: The first step in this code is to use the MTEX
    % calcgrain tool to find likely martensite grains, then merge grains
    % that share boundaries that are reasonably close to those expected
    % between martensite varients. It then takes the largest region from
    % this and discards the rest. The result, ideally, should be a region
    % containing at least one (but likely a few) Prior Austenite grain(s).
    % the next step is to convert_to_fundamental_zone, find the strongest
    % ODF concentration (which should be equal to the orientation of the
    % largest PAG), then do a graph cut on the map, which ideally will pull
    % out one and ONLY one PAG. This is all done inside the function
    % "HT_grain_guess.m"
    % Now, IN THE ORIGINAL VERSION, the function gave back two ebsds: the
    % first was the Austenite grain guess, the second was the LEFTOVER from
    % the MTEX grain, NOT the leftover of the entire ebsd. If the MCMC
    % Auto_OR portion failed, the next attempt would only search the
    % leftovers, not the whole map. This seems like a mistake:
    %   1). since it is possible to actually just get 1 grain from MTEX,
    %   this makes it very easy to end up in situations where you run out
    %   of iterations before finding a good grain that passes the following
    %   MCMC OR guess in this function.
    %   2) If you have a small PAG size, you can get forced quickly into a
    %   situation where you are making these estimations on just a few
    %   dozen pixels, which seems bad.
    %   3) why limit your data pool?
    %ALL THAT SAID THOUGH, there might be something happening here I don't
    %fully understand (after all, old version had a 90% success rate,) so I
    %am switching it over, BUT also leaving a coded out version to return
    %to the leftovers-only search, JIC. Feel free to try either way.
    %% Method 1: Search all unused areas
    % find indices of a possible grain
    [Aus_Grain_Ids,~,AusOr] = HT_grain_guess_Alex(searchable_ebsd);
    % extract the expected grain and leave the remainder for future searches
    grain_ebsd = searchable_ebsd(Aus_Grain_Ids);
    searchable_ebsd(Aus_Grain_Ids) = [];
    %% Method 2: Search only MTEX leftovers on subsequent searches
    %    % find indices of a possible grain
    %    [Aus_Grain_Ids,grain_remainder_Ids,AusOr] = HT_grain_guess(searchable_ebsd);
    %    % extract the expected grain and leave the remainder for future searches
    %    grain_ebsd = searchable_ebsd(Aus_Grain_Ids);
    %    searchable_ebsd = searchable_ebsd(grain_remainder_Ids)
    %%
    % remove everything from this grain that isnt the LT phase (phase 2)
    %    grain_ebsd = grain_ebsd(grain_ebsd.phase == 2);
    % We don't need the whole grain, just a representative sample of the
    % Martensite from a single Austenite grain. Plus, a 2x larger sample
    % causes an 8x slowdown in the graph cut, so time to downsample
    martensite=grain_ebsd.orientations;
    martensite.CS = CS_LT;
    keep=randperm(length(martensite));
    if length(martensite) > options.OR_sampling_size
        martensite=martensite(keep(1:options.OR_sampling_size));
    else
    end
    
    % plot the martensite pole figure and grain if desired
    if options.OR_plot_PAG_Mart_Guess == 1
        figure();
        plotPDF(martensite,Miller({0,0,1},martensite.CS),'antipodal','points','all','marker','.');
        figure();
        plot(grain_ebsd("Martensite"),grain_ebsd("Martensite").orientations)
    end
    if options.OR_plot_PAG_Mart_Guess == 1
    end
    
    %% Generate initial guess for parameters to be estimated
    % initial guess for ksi and halfwidth.
    ksi_initial=[5.26 10.3 10.53];      % KS
    halfwidth_in=2.5*degree;
    % make sure suggested ksi values make sense
    [~,flag]=calc_T2R(ksi_initial,CS_HT,CS_LT);
    if flag
        warning('non-physical initial ksi angles!');
    end
    % Other options for guesses (these should be moved to a help document
    % ksi_initial=[5.26 10.3 10.53];      % KS
    % ksi_initial=[0.00,9.74,9.74];       % NW
    % ksi_initial=[3.3,8.5,8.9];        % KS-Like
    % ksi_initial=[4.37,7.79,9.24];
    % ksi_initial=[2.3114    8.9269    9.1332];  % 1st Samp Exp Obs
    % ksi_initial = [2.986,8.228,8.584];      % Crop3 Ksi average for good fits
    
    % clear unneeded variables before next section
    clear Aus_Grain_Ids keep flag
    %% set parameters for inference
    
    % Set mu and sigma for prior probability on ksi. Values based on rough
    % Estimates from Yardley Payton 2014 conference paper
    % ksi_prior_mu=[3,8,9];
    ksi_prior_mu = [5,9,10];
    % ksi_prior_sigma=[1.2,1.2,1.2];
    ksi_prior_sigma = [2,2,2];
    
    % Noise is modeled as a unimodal odf centered on the cube orientation.
    % Parameter to be estimated by Bayesian inference is the halfwidth of the
    % odf kernel. Halfwidth distribution is assumed folded Gaussian. This
    % approximates uniform from 0->~1 degree then decaying at larger noise values
    halfwidth_prior_mu=1;
    halfwidth_prior_sigma=2;
    
    %austenite - twp prior considered for austenite based on cases. If ksi
    %angles can be carefully measured then austenite prior will be a unimodal
    %odf about the global pole figure modal ODF. If initial ksi values are
    %approximate then austenite prior is uniform over the fundamental zone
    
    %leftover can ignore
    austenite_prior_odf=uniformODF(CS_HT,SS);
    
    prior_pars=struct;
    prior_pars.ksi_prior_mu=ksi_prior_mu;
    prior_pars.ksi_prior_sigma=ksi_prior_sigma;
    prior_pars.halfwidth_prior_mu=halfwidth_prior_mu;
    prior_pars.halfwidth_prior_sigma=halfwidth_prior_sigma;
    prior_pars.austenite_prior_odf=austenite_prior_odf;
    prior_pars.CS_A=CS_HT;
    prior_pars.CS_M=CS_LT;
    prior_pars.SS=SS;
    
    %% MAP estimate of parameters by optimization
    
    MAP_options=optimset('fminsearch');
    MAP_options=optimset(MAP_options,'display','iter');
    MAP_options.TolFun = options.OR_fit_TolFun;
    MAP_options.TolX = options.OR_fit_TolX;
    optimfunc= @(samples) -posterior_pdf_fminunc(samples,prior_pars,martensite);
    x0=[ksi_initial,halfwidth_in/degree];
    
    % Ensure the constraint that ksi_1 < ksi_2 and ksi_3. Since we
    % don't know for certain if ksi_2 SHOULD be < ksi_3, force this to
    % occur if the first constraint is met.
    init_guess=0;
    count = 0;
    while init_guess == 0
        count = count+1;
        % Optimization function outside of ML add-on
        [MAPpars,~,~,~]=fminsearch(optimfunc,x0,MAP_options);
        keyboard
        if (MAPpars(1) < MAPpars(2) && MAPpars(3)) || count > 2
            init_guess = 1;
            % Enforce constraint if need be
            if(MAPpars(2) > MAPpars(3))
                MAPpars(2) = (MAPpars(3)-1e-3);
            end
        end
    end
    ksi_1st_guess=MAPpars(1:3);
    halfwidth_in=MAPpars(4)*degree;
    
    %% Set MCMC sampler parameters
    
    %burnin is number of samples to disregard at the beginning maybe try 100
    %for paper
    burnin=0;
    num_samples=1;
    
    %ksi values update parameters
    scale=0.15;
    ksi1_width=.2*scale;
    ksi2_width=.08*scale;
    ksi3_width=.08*scale;
    
    %leftover needs to be cleaned up
    austenite_current=AusOr;
    
    %% componentwise MCMC sampler
    
    ksi1_current=ksi_1st_guess(1);
    ksi2_current=ksi_1st_guess(2);
    ksi3_current=ksi_1st_guess(3);
    HW=halfwidth_in;
    
    %austenite_posterior=austenite_current;
    ksi1_posterior=ksi1_current;
    ksi2_posterior=ksi2_current;
    ksi3_posterior=ksi3_current;
    
    
    M=martensite;
    num_mart_Ors=100;%;num_martensite_orientations;
    for ii=1:num_samples+burnin
        
        mm=randperm(length(M));
        martensite=M(mm<=num_mart_Ors);
        
        p_current=martensite_posterior_log_likelihood(martensite,[ksi1_current,ksi2_current,ksi3_current],HW,austenite_current,prior_pars);
        
        flag=1;
        flag_counter=0;
        skip=0;
        while flag
            flag_counter=flag_counter+1;
            ksi3_proposal=abs(randn*ksi3_width+ksi3_current);
            [~, flag]=YardleyVariants([ksi1_current,ksi2_current,ksi3_proposal]);
            if flag_counter>=10
                flag=0;
                ksi3_proposal=ksi3_current;
                skip=1;
            end
        end
        
        [T2R,~]=calc_T2R([ksi1_current,ksi2_current,ksi3_proposal],CS_HT,CS_LT);
        
        austenite=symmetrise(martensite)*T2R;
        
        [austenite_proposal] = global_pole_figure_estimation(austenite,CS_HT,SS,1);
        
        for kk=1:length(austenite_proposal)
            temp(kk)=martensite_posterior_log_likelihood(martensite,[ksi1_current,ksi2_current,ksi3_proposal],HW,austenite_proposal(kk),prior_pars);
        end
        id=find(temp==max(temp),1);
        p_proposal=temp(id);
        clear temp
        
        %stricly the log probability
        p_accept=p_proposal-p_current;
        
        if skip
            accept = false;
        else
            accept=log(rand)<p_accept;
        end
        
        if accept
            austenite_current=austenite_proposal(id);
            ksi3_current=ksi3_proposal;
        end
        
        ksi3_posterior=[ksi3_posterior;ksi3_current];
        
        % ksi 2 sample in order ksi3, ksi2, ksi1, halfwidth
        
        p_current=martensite_posterior_log_likelihood(martensite,[ksi1_current,ksi2_current,ksi3_current],HW,austenite_current,prior_pars);
        
        flag=1;
        flag_counter=0;
        skip=0;
        while flag
            flag_counter=flag_counter+1;
            ksi2_proposal=abs(randn*ksi2_width+ksi2_current);
            [~, flag]=YardleyVariants([ksi1_current,ksi2_proposal,ksi3_current]);
            if flag_counter>=10
                flag=0;
                %ksi2_proposal=ksi2_currrent;
                skip=1;
            end
        end
        
        [T2R,~]=calc_T2R([ksi1_current,ksi2_proposal,ksi3_current],CS_HT,CS_LT);
        
        austenite=symmetrise(martensite)*T2R;
        
        [austenite_proposal] = global_pole_figure_estimation(austenite,CS_HT,SS,1);
        
        for kk=1:length(austenite_proposal)
            temp(kk)=martensite_posterior_log_likelihood(martensite,[ksi1_current,ksi2_proposal,ksi3_current],HW,austenite_proposal(kk),prior_pars);
        end
        id=find(temp==max(temp),1);
        p_proposal=temp(id);
        clear temp
        
        p_accept=p_proposal-p_current;
        
        if skip
            accept = false;
        else
            accept=log(rand)<p_accept;
        end
        
        if accept
            austenite_current=austenite_proposal(id);
            ksi2_current=ksi2_proposal;
        end
        
        ksi2_posterior=[ksi2_posterior;ksi2_current];
        
        % ksi 1 sample in order ksi3, ksi2, ksi1, halfwidth
        
        p_current=martensite_posterior_log_likelihood(martensite,[ksi1_current,ksi2_current,ksi3_current],HW,austenite_current,prior_pars);
        
        flag=1;
        flag_counter=0;
        skip=0;
        while flag
            flag_counter=flag_counter+1;
            ksi1_proposal=abs(randn*ksi1_width+ksi1_current);
            [~, flag]=YardleyVariants([ksi1_proposal,ksi2_current,ksi3_current]);
            if flag_counter>=10
                flag=0;
                ksi1_proposal=ksi1_current;
                skip=1;
            end
        end
        
        [T2R,~]=calc_T2R([ksi1_proposal,ksi2_current,ksi3_current],CS_HT,CS_LT);
        
        austenite=symmetrise(martensite)*T2R;
        
        [austenite_proposal] = global_pole_figure_estimation(austenite,CS_HT,SS,1);
        
        for kk=1:length(austenite_proposal)
            temp(kk)=martensite_posterior_log_likelihood(martensite,[ksi1_proposal,ksi2_current,ksi3_current],HW,austenite_proposal(kk),prior_pars);
        end
        id=find(temp==max(temp),1);
        p_proposal=temp(id);
        clear temp
        
        p_accept=p_proposal-p_current;
        
        if skip
            accept = false;
        else
            accept=log(rand)<p_accept;
        end
        
        if accept
            austenite_current=austenite_proposal(id);
            ksi1_current=ksi1_proposal;
        end
        ksi1_posterior=[ksi1_posterior;ksi1_current];
    end
    
    % If the likelihood value is large enough, conclude with this OR.
    % If not, and we've gone through 10 iterations with no successful
    % OR, return error and make user input manually.
    if p_proposal > 4e2
        break
    elseif iteration == 10
        error('Automatic Determination of OR Failed. Try Manual Segmentation of Potential PAG!')
    end
end
%% Consolidate ksi data to one term
ksi_out=[ksi1_posterior(burnin+1:end),ksi2_posterior(burnin+1:end),ksi3_posterior(burnin+1:end)];
ksi_curr = mean(ksi_out);
% clear out variables no longer needed
clear searchable_ebsd

%%    Plot if desired
% Plot the Discrete ODFs
% [R2T,~]=calc_R2T(ksi_curr',CS_HT,CS_LT);
% curr_mart = symmetrise(austenite_current)*R2T;
% curr_mart_odf = calcODF(curr_mart,CS_LT,SS,CS_HT);
% figure; plotPDF(curr_mart_odf,Miller({0,0,1},{1,1,0},{1,1,1},CS_M),'antipodal','points','all','marker','.');
if options.OR_plot_ODF_of_PAG_OR_guess == 1
    martensite=generate_simulated_data(austenite_proposal,[ksi_curr(1),ksi_curr(2),ksi_curr(3)],HW,2000,CS_LT);
    figure; plotPDF(martensite,Miller({0,0,1},martensite.CS),'antipodal','points','all','marker','.');
end
if options.OR_plot_ksi_spread == 1
    for i = 1:3
        figure; histogram(ksi_out(:,i),20)
        if i == 1
            title('\xi_1')
        elseif i == 2
            title('\xi_2')
        else
            title('\xi_3')
        end
        xlim([min(ksi_out(:,i))-1,max(ksi_out(:,i)+1)])
    end
end
%    figure()
%    aaa = searchable_ebsd(searchable_ebsd.phaseId == 2);
%    plot(aaa("Martensite"),aaa("Martensite").orientations)
%    figure()
%    aaa = grain_ebsd(grain_ebsd.phaseId == 2);
%    plot(aaa("Martensite"),aaa("Martensite").orientations)

%% Output result
% actual ksi values and kernel halfwidth
final_OR = ksi_curr;
final_halfwidth = HW;
OR_metadata.ORLikelihood = p_proposal;
OR_metadata.ksi.ksi1 = ksi_out(:,1);
OR_metadata.ksi.ksi2 = ksi_out(:,2);
OR_metadata.ksi.ksi3 = ksi_out(:,3);
OR_metadata.iterations = iteration;
OR_metadata.solve_time = toc;
end
%%

% at this point, we have identical scans with the following phase IDs:
%-1 : Unindexed
% 1 : Untransformed Parent (High temperature, or HT)
% 2 : Transformed Child (Low Temperature, or LR)
% 3 : Reconstructed Parent (starts empty) (Reconstructed, or R)