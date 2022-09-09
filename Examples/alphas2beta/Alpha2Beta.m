




%%
% ======================================================================= %
%       ab_ti.m
%       Proof of concept for Example for reconstructing multiple EBSD 
%       scans where the orientation relationship is not known a priori
% ======================================================================= %

% Set some options
options = load_options("default");
m = options.RGC_in_plane_m;
b = options.RGC_in_plane_b;

tic;
% load in some files
%ebsd = EBSD.load('small_ti.ctf');
mtexdata alphaBetaTitanium
% manual phase map fix

% high temp phase
beta = ebsd(ebsd.phase==1);
csbeta = beta.CS;
% low temp phase
alpha = ebsd(ebsd.phase==2);
csalpha = alpha.CS;

% Orientation Relationships
beta2alpha = orientation('map',Miller(1,1,0,beta.CS),Miller(0,0,0,1,alpha.CS),...
        Miller(-1,1,-1,beta.CS),Miller(2,-1,-1,0,alpha.CS));
alpha2beta=inv(beta2alpha);
tocA = toc;

tic;
% building the MDF
oriParent=orientation.id(csbeta);
oriVariants=unique(oriParent.symmetrise * inv(beta2alpha));
high_low_misos = inv(orientation.id(csalpha)).*(oriVariants);
low_low_misos = inv(oriVariants).*(oriVariants);
pc_kernel = deLaValleePoussinKernel('halfwidth',2*degree,'bandwidth',96);
cc_kernel = deLaValleePoussinKernel('halfwidth',2*degree,'bandwidth',96);
parent_child_MDF = calcDensity(high_low_misos,'kernel',pc_kernel);
variants_MDF = calcDensity(high_low_misos,'kernel',cc_kernel);


% precalculate the  neighbor list
fprintf("finding friends...\n")
neigh_list = find_neigh_ti(alpha,2);
tocB = toc;

tic;
% precalculate In-plane weights
fprintf("ranking those friends...\n")
IP_wts = get_In_plane_weights(neigh_list,alpha,variants_MDF,m,b);

% Now take all the precalculated data required for doing the reconstruction
% and attach it to the EBSD structure under "ebsd.opt"
ebsd.opt.parent_child_MDF = parent_child_MDF;
ebsd.opt.variants_MDF = variants_MDF;
ebsd.opt.pc_kernel = pc_kernel;
ebsd.opt.LT_MDF = variants_MDF;
ebsd.opt.alpha2beta = alpha2beta;
ebsd.opt.beta2alpha = beta2alpha;
% NOTE: this is the point where users should change inputs if they deire,
% such as kernel sizes, graph weights, etc. Everything up to here is the
% setup for the actual reconstruction
ebsd.opt.max_recon_attempts = 500;
ebsd.opt.min_cut_size = 5;

% Do the reconstruction
[recon_ebsd, likelihoods] = Reconstruct_HT_grains_ti(ebsd, neigh_list, IP_wts);


%% ============================= %%

% load options, cause Austin hard coded things in like a noob
options = load_options("default");
% precalculate the  neighbor list
fprintf("finding friends...\n")
neigh_list = find_neigh_ti(alpha);
tocA = toc;
% find the in-plane weights, which never change
fprintf("calculating In-plane weights...\n")
m = options.RGC_in_plane_m;
b = options.RGC_in_plane_b;
IP_wts = get_In_plane_weights(neigh_list,alpha,variant_MDF,m,b);
tocB = toc;



ebsd.opt.LT_MDF = variant_MDF;
ebsd.opt.alpha2beta = alpha2beta;
ebsd.opt.beta2alpha = beta2alpha;
ebsd.opt.psi = kernel;
[recon_ebsd, LT_likelihoods] = Reconstruct_HT_grains_ti(ebsd,...
    neigh_list,...
    IP_wts,...
    options);

%% ============================= %%
% % Now copy the original ebsd and replace the old values for the LT phases
% % with the new ones.
% temp_ebsd = recon_ebsd;
% temp_ebsd.phase = 2;
% HT_ebsd = orig_ebsd;
% HT_ebsd(orig_ebsd.phaseId == 2).orientations = temp_ebsd.orientations;
% HT_ebsd(orig_ebsd.phaseId == 2).phaseId = recon_ebsd.phaseId;
% likelihoods(orig_ebsd.phaseId == 2) = LT_likelihoods;
% likelihoods = transpose(likelihoods);
% end






    % Segment the variants using the low temp and high temp maps
    Tasks(i).variant_int_map = variants(Tasks(i).ebsd, ...
        Tasks(i).Recon_ebsd,Tasks(i).options);
% ================================================================ %
%                              Plotting                            %
% ================================================================ %
% These scans have a lot going on, so for ease we have provided
% some examples of how to plot this data.
O = Tasks(9).ebsd;
R = Tasks(9).Recon_ebsd;
options = Tasks(9).options;

MO = O(O.phaseId == 2);
[gmap, MO.grainId] = calcGrains(R(R.phaseId == 3));

v_map = variants(O, R, options);
key = ipfHSVKey(MO);
mp = key.orientation2color(MO.orientations);
steel_cmap = get_subgrain_colormaps('Steel');

figure()
for i = 1:24
    mask = rem(v_map,24)==(rem(i,24));
%    mask = logical(mask.*transpose(MO.grainId ==mode(MO.grainId)));
    plot(MO(mask),mp(mask),'FaceColor',steel_cmap.Variant(i,:))
    hold on
end










%% Steve example below here

clear all
close all

%load OR_stuff.mat
csAlpha = crystalSymmetry('622',[3 3 4.7],'mineral','Ti (alpha)');
csBeta = crystalSymmetry('432',[3.3 3.3 3.3],'mineral','Ti (beta)');


beta2alpha = orientation('map',Miller(1,1,0,csBeta),Miller(0,0,0,1,csAlpha),...
        Miller(-1,1,-1,csBeta),Miller(2,-1,-1,0,csAlpha));

alpha2beta=inv(beta2alpha);

%%
%oriParent=orientation.id(csBeta);
oriParent=orientation.rand(csBeta);
oriChild=oriParent * inv(beta2alpha);

oriVariants=unique(oriParent.symmetrise*inv(beta2alpha));



%%
%stupid test # 1

possible_parents=(oriVariants.symmetrise*inv(alpha2beta));

psi=calcKernel(possible_parents);

global_odf=calcODF(possible_parents,'kernel',psi);

predicted_beta=calcModes(global_odf);

error_test_1=angle(predicted_beta,oriParent)
%%
%id_alpha=orientation('euler',1e-6,1e-6,0,csAlpha);
%end_alpha=orientation('euler',1e-6,1e-6,60*degree,csAlpha);

%f=fibre(id_alpha,end_alpha);

%fodf = fibreODF(f,'halfwidth',10*degree);

%k=calcKernel(oriVariants);
%vodf=calcODF(oriVariants,'kernel',k);

%% calculate alpha variant fibers and variant ODF



for ii=1:12
    phi1=oriVariants.phi1(ii);
    Phi=oriVariants.Phi(ii);
    
    start_orientation=orientation('euler',phi1,Phi,0,csAlpha);
    end_orientation=orientation('euler',phi1,Phi,60*degree,csAlpha);
    f=fibre(start_orientation,end_orientation);
    
    temp_odf=fibreODF(f,'halfwidth',0.5*degree);
    
    if ii==1
        variant_odf=temp_odf;
    else
        variant_odf=variant_odf+temp_odf;
    end
    
end

variant_orientations=calcOrientations(variant_odf,10000);
%% Project variant ODF into beta space


possible_betas=(variant_orientations.symmetrise*inv(alpha2beta));

%k=calcKernel(possible_betas);

possible_beta_odf=calcODF(possible_betas);

%% see if max is at true beta orientation

predicted_beta=calcModes(possible_beta_odf);

error_test_2=angle(predicted_beta,oriParent)