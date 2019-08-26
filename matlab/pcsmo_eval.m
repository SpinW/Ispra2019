% SpinW / Horace integration example
%
% This tutorial illustrates how to evaluate a (fairly complex) SpinW
% model on Horace ( http://horace.isis.rl.ac.uk ) INS data.
%
% The dataset are sets of measurements of Pr(Ca0.9Sr0.1)2Mn2O7 using
% the MAPS spectrometer at ISIS at different neutron energies.

% Loop through the different energy data and make 2D cuts along (h00) vs E
ei = [25 35 50 70 140];
qstp = [0.02 0.02 0.025 0.025 0.05];
estp = ei ./ 100;
lint = [0.5 0.2 0.2 0.2 0.2];
proj = projaxes([1 0 0], [0 1 0]);
for ii = 1:numel(ei);
    sqw_file = sprintf('C:/HoraceWorkshop/data/pcsmo_ei%d_base.sqw', ei(ii));
    ws_cut(ii) = cut_sqw(sqw_file, proj, [-5,qstp(ii),5], [-1,1]*lint(ii), [-10,10], [estp(ii)]);
    % Symmetrise about h=0
    ws_cut(ii) = symmetrise_sqw(ws_cut(ii), [0 0 1], [0 1 0], [0 0 0]);
    plot(ws_cut(ii))
    lz 0 10
    keep_figure;
end

%%
% Make cuts at high Q to use as a phonon background for subtraction
for ii = 1:numel(ei)
    idx = find(sum(ws_cut(ii).data.s,2)>0);
    qmax = ws_cut(ii).data.p{1}(idx(end)) * 0.5;
    ws_bkg(ii) = cut_sqw(ws_cut(ii), [qmax-5*qstp(ii),qmax], []);
    % Uncomment the code below to generate plots of the q range
    % of the background cuts and what they look like
    %plot(ws_cut(ii)); hold all; plot([qmax-5*qstp(ii),qmax],[0 0],'-r'); keep_figure;
    %plot(ws_bkg(ii)); keep_figure;
    ws_sub(ii) = ws_cut(ii) - replicate(ws_bkg(ii), ws_cut(ii));
    plot(ws_sub(ii))
    lz 0 10
    lx -2 2
    keep_figure;
end

%%
% Try to reproduce figure 2 of Johnstone et al., PRL 109 237202 (2012)

% First make the required cuts
ws25 = cut_sqw(ws_sub(1), [0.2,qstp(1),1], [0,estp(1),15], '-nopix');
ws35 = cut_sqw(ws_sub(2), [0.2,qstp(2),1], [15,estp(2),30], '-nopix');
ws50 = cut_sqw(ws_sub(3), [0,qstp(3),0.2], [0,estp(3),30], '-nopix');
ws70 = cut_sqw(ws_sub(4), [0,qstp(4),2], [0,estp(4),45], '-nopix');
ws140 = cut_sqw(ws_sub(5), [0,qstp(5),2], [45,estp(5),100], '-nopix');
% Plot them
plot(ws25); lz 0 10; keep_figure;
plot(ws35); lz 0 10; keep_figure;
plot(ws50); lz 0 10; keep_figure;
plot(ws70); lz 0 10; keep_figure;
plot(ws140); lz 0 10; keep_figure;

% Now parse through them and extract the data matrix
wsarr = [ws25 ws35 ws50 ws70 ws140];
for ii = 1:numel(wsarr);
    % Horace stores bin boundaries rather bin centers in coordinate axes
    dat(ii).x = (wsarr(ii).p{1}(1:(end-1))+wsarr(ii).p{1}(2:end))/2;
    dat(ii).y = (wsarr(ii).p{2}(1:(end-1))+wsarr(ii).p{2}(2:end))/2;
    dat(ii).z = wsarr(ii).s;
end

% Now use normal Matlab plotting commands to make the composite plots
figure
subplot(121); hold all;
    pcolor(dat(1).x, dat(1).y, dat(1).z');
    pcolor(dat(2).x, dat(2).y, dat(2).z');
    pcolor(dat(3).x, dat(3).y, dat(3).z');
    shading flat;
    caxis([0 10])
    plot([1 1]*0.2, [0 30], '-r');
    plot([0.2 1], [1 1]*15, '-r');
    ylim([0 30]);
    xlabel('(h,0) (rlu)');
    ylabel('Energy (meV)');
subplot(122); hold all;
    pcolor(dat(4).x, dat(4).y, dat(4).z');
    pcolor(dat(5).x, dat(5).y, dat(5).z');
    shading flat;
    caxis([0 10]);
    plot([0 2], [45 45], '-r');
    xlabel('(h,0) (rlu)');
    ylabel('Energy (meV)');

% Why don't we see the same as in Johnstone's Fig 2?
    
%%
% Spin Hamiltonian parameters from Johnstone et al.
JF1 = -11.39;
JA = 1.5;
JF2 = -1.35;
JF3 = 1.5;
Jperp = 0.88;
D = 0.074;

% Define the SpinW model - same as the file prcasrmn2o7.m
% From the "real world" systems tutorial.
lat = [5.408 5.4599 19.266];
alf = [90 90 90];
SM4 = 7/4;   % Spin length for Mn4+
SM3 = 7/4;   % Spin length for Mn3+
pcsmo = spinw;
pcsmo.genlattice('lat_const', lat.*[2 2 1], 'angled', alf, 'spgr', 'x,y+1/2,-z'); 
[~,ffn3] = sw_mff('MMn3');
[~,ffn4] = sw_mff('MMn4');
myaddatom3 = @(x,y,z) pcsmo.addatom('label', x, 'r', y, 'S', SM3, 'color', z, ...
    'formfactn', ffn3, 'formfactx', 'MMn3', 'Z', 25, 'b', sw_nb('MMn3'));
myaddatom4 = @(x,y,z) pcsmo.addatom('label', x, 'r', y, 'S', SM4, 'color', z, ...
    'formfactn', ffn4, 'formfactx', 'MMn4', 'Z', 25, 'b', sw_nb('MMn4'));
myaddatom4('Mn4-up', [0 0 0.1], 'gold');
myaddatom4('Mn4-up', [0.5 0.5 0.1], 'gold');
myaddatom4('Mn4-dn', [0 0.5 0.1], 'gold');
myaddatom4('Mn4-dn', [0.5 0 0.1], 'gold');
myaddatom3('Mn3-up', [0.25 0.75 0.1], 'black');
myaddatom3('Mn3-up', [0.75 0.75 0.1], 'black');
myaddatom3('Mn3-dn', [0.25 0.25 0.1], 'black');
myaddatom3('Mn3-dn', [0.75 0.25 0.1], 'black');
% Generate the CE magnetic structure
S0 = [0; 1; 0];
spin_up = find(~cellfun(@isempty, strfind(pcsmo.table('matom').matom, 'up')));
spin_dn = find(~cellfun(@isempty, strfind(pcsmo.table('matom').matom, 'dn')));
SS = zeros(3, 16);
SS(:, spin_up) = repmat(S0, 1, numel(spin_up));
SS(:, spin_dn) = repmat(-S0, 1, numel(spin_dn));
pcsmo.genmagstr('mode', 'direct', 'S', SS)
% Generate the exchange interactions
pcsmo.gencoupling('forceNoSym', true)
pcsmo.addmatrix('label', 'JF1', 'value', JF1, 'color', 'green');
pcsmo.addmatrix('label', 'JA', 'value', JA, 'color', 'yellow');
pcsmo.addmatrix('label', 'JF2', 'value', JF2, 'color', 'white');
pcsmo.addmatrix('label', 'JF3', 'value', JF3, 'color', 'red');
pcsmo.addmatrix('label', 'Jperp', 'value', Jperp, 'color', 'blue');
pcsmo.addmatrix('label', 'D', 'value', diag([0 0 D]), 'color', 'white');
% The zig-zag chains couple Mn3-Mn4 with same spin.
pcsmo.addcoupling('mat', 'JF1', 'bond', 1, 'atom', {'Mn3-up', 'Mn4-up'})
pcsmo.addcoupling('mat', 'JF1', 'bond', 1, 'atom', {'Mn3-dn', 'Mn4-dn'})
% And vice-versa for the inter-chain interaction
pcsmo.addcoupling('mat', 'JA', 'bond', 1, 'atom', {'Mn3-up', 'Mn4-dn'})
pcsmo.addcoupling('mat', 'JA', 'bond', 1, 'atom', {'Mn3-dn', 'Mn4-up'})
pcsmo.addcoupling('mat', 'Jperp', 'bond', 2)
% JF3 couples Mn3 within the same zig-zag (same spin)
pcsmo.addcoupling('mat', 'JF3', 'bond', 3, 'atom', 'Mn3-up')
pcsmo.addcoupling('mat', 'JF3', 'bond', 3, 'atom', 'Mn3-dn')
% Find indexes of the Mn4+ atoms which have a=0.5:
idmid = find((~cellfun(@isempty, strfind(pcsmo.table('matom').matom, 'Mn4'))) ...
    .* (pcsmo.table('matom').pos(:,1)==0.5));
bond8 = pcsmo.table('bond', 8);
% Finds the bonds which start on one of these atoms and goes along +b
idstart = find(ismember(bond8.idx1, idmid) .* (bond8.dr(:,2)>0));
% Finds the bonds which ends on one of these atoms and goes along -b
idend = find(ismember(bond8.idx2, idmid) .* (bond8.dr(:,2)<0));
pcsmo.addcoupling('mat', 'JF2', 'bond', 8, 'subIdx', [idstart; idend]')
pcsmo.addaniso('D')
% Optimise structure and plot a picture of the model
res = pcsmo.optmagsteep()
plot(pcsmo, 'range', [0 1; 0 1; 0 0.2])

%%
% Now try to evaluate the SpinW on the 70meV and 140meV cuts
% Using a Gaussian resolution function
cpars = {'mat', {'JF1', 'JA', 'JF2', 'JF3', 'Jperp', 'D(3,3)'}, ...
    'hermit', false, 'optmem', 0, 'useFast', true, 'formfact', true, ...
    'resfun', 'gauss', 'coordtrans', diag([2 2 1 1])};
% We have to scale qh and qk by 2 because the unit cell we have to use
% in SpinW is 2x2 times the high temperature unit cell used in the data.
swpref.setpref('usemex',true);

idx = 5;  % EIs: [25 35 50 70 140], index 5 == 140meV
fwhm = ei(idx)/30;  % Typical resolution ~ 3% of Ei

% Do evaluation on DnD object - much fewer pixels, so faster computation
% - The sqw objects have ~100x more q points - so will take >100x longer!
% - with the danger that you will run out of memory too!
kk = multifit_sqw(d2d(ws_sub(idx)));
kk = kk.set_fun (@pcsmo.horace_sqw, {[JF1 JA JF2 JF3 Jperp D fwhm] cpars{:}});
kk = kk.set_free ([1, 1, 1, 1, 1, 1, 1]);
kk = kk.set_options ('list',2);

% Time a single iteration
tic
wsim = kk.simulate;
t_spinw_single = toc;

plot(ws_sub(idx));
lz 0 10
keep_figure;
plot(wsim);
keep_figure;

% Put the two figures on one plot
wss = symmetrise_sqw(ws_sub(idx),[0 1 0],[0 0 1],[0 0 0]);
spaghetti_plot([compact(d2d(wss)) compact(wsim)])
lz 0 2