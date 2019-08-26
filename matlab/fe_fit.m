% Fitting data to a SpinW model using Horace example 
%
% This tutorial will demonstrate how to fit a single crystal inelastic
% neutron scattering dataset to a SpinW using the Horace 
% ( http://horace.isis.rl.ac.uk ) package for data handling.
%
% The dataset is a measurement of the spinwaves of bcc-Fe using the MAPS
% spectrometer at the ISIS facility.
%
% The dataset, and other materials from the Horace course run periodically
% by ISIS may be found here:  
% ftp://ftp.nd.rl.ac.uk/scratch/Ewings/HoraceWorkshop2017.zip
%
% Edit this line with the location of the sqw file
sqw_file = 'C:/HoraceWorkshop/data/fe.sqw';

% If the data file doesn't exist, create a fake version
if ~exist(sqw_file, 'file')
    % Use the Horace "demo" generating scripts 
    % - please be patient: this will take about 10 min.
    cur_dir = pwd;
    demo_dir = [fileparts(which('fake_data')) '/demo'];
    cd(demo_dir);
    setup_demo_data;
    spes = arrayfun(@(x)sprintf('HoraceDemoDataFile%d.nxspe',x), 1:23, 'UniformOutput', 0);
    gen_sqw(spes, '', sqw_file, 787, 1, [1 1 1]*2.87, [1 1 1]*90, ...
        [1 0 0], [0 1 0], 0:4:90, 0, 0, 0, 0);
    delete('*.nxspe');
    cd(cur_dir);
end

% Make a series of 1D cuts of the data
proj.u  = [1,1,0]; 
proj.v  = [-1,1,0]; 
proj.type  = 'rrr';
energy_range = [80:20:160];
for i = 1:numel(energy_range)
    my_cuts(i) = cut_sqw(sqw_file, proj, [-3,0.05,3], [-1.05,-0.95], [-0.05,0.05], ...
        [-10 10]+energy_range(i));
end

%%
% Run the fitting with an analytical expression for the cross-section
%  - using an exact expression for the dispersion of a body-centered cubic 
%    FM combined with a damped harmonic oscillator for the intensity.
% The parameters of the function is: [JS Delta Gamma Temperature Amplitude]
p0 = [35 0 30 10 300];
% Expression for the dispersion
om = @(h,k,l,e,js,d) d + (8*js)*(1-cos(pi*h).*cos(pi*k).*cos(pi*l));
% The magnetic form factor of Fe2+
A=0.0706; a=35.008;  B=0.3589; b=15.358;  C=0.5819; c=5.561;  D=-0.0114;
q2 = @(h,k,l) (1/(2*2.87)).^2 .* (h(:).^2 + k(:).^2 + l(:).^2);
ff = @(h,k,l) sum(bsxfun(@times, [A B C], exp(bsxfun(@times, -[a b c], q2(h,k,l)))), 2) + D;
% Put it altogether
fe_sqw = @(h,k,l,e,p) ff(h,k,l).^2 .* (p(5)/pi) .* (e./(1-exp(-11.602.*e./p(4)))) .* ...
    4.*p(3).*om(h,k,l,e,p(1),p(2)) ./ ((e.^2-om(h,k,l,e,p(1),p(2)).^2).^2 + 4.*(p(3).*e).^2);

% Starting parameters for fit
J = 35;
D = 0;
gam = 30;
temp = 10;
amp = 300;

kk = multifit_sqw (my_cuts);
kk = kk.set_fun (fe_sqw, [J D gam temp amp]);
kk = kk.set_free ([1, 0, 1, 0, 1]);
kk = kk.set_bfun (@linear_bg, [0.1,0]);
kk = kk.set_bfree ([1,0]);
kk = kk.set_options ('list',2);

% Time it to see how long it takes to do the fit
tic
[wfit, fitdata] = kk.fit('comp');
t_ana = toc;

for i=1:numel(my_cuts)
    acolor black
    plot(my_cuts(i));
    acolor red
    pl(wfit.sum(i));
    pl(wfit.back(i));
    keep_figure;
end

%%
% Now setup the SpinW model and try to run it
a = 2.87;

fe = spinw;
fe.genlattice('lat_const', [a a a], 'angled', [90 90 90], 'spgr', 'I m -3 m')  % bcc Fe
fe.addatom('label', 'MFe3', 'r', [0 0 0], 'S', 5/2, 'color', 'gold')
fe.gencoupling()
fe.addmatrix('label', 'J1', 'value', 1, 'color', 'gray')
fe.addmatrix('label', 'D', 'value', diag([0 0 -1]), 'color', 'green')
fe.addcoupling('mat', 'J1', 'bond', 1)
fe.addaniso('D')
fe.genmagstr('mode', 'direct', 'S', [0 0 1; 0 0 1]');  % Ferromagnetic

plot(fe, 'range', [2 2 2])

% Constant parameters for SpinW model
% Note that we use the damped harmonic oscillator resolution model ('sho')
cpars = {'mat', {'J1', 'D(3,3)'}, 'hermit', false, 'optmem', 1, ...
    'useFast', true, 'resfun', 'sho', 'formfact', true};
swpref.setpref('usemex',true);

% Initial parameters:
J = -16;     % Exchange constant in meV - Note previous value was J*S (S=2.5)
D = -0.1;    % SIA constant in meV
gam = 66;
temp = 10;
amp = 131;

kk = multifit_sqw (my_cuts);
kk = kk.set_fun (@fe.horace_sqw, {[J D gam temp amp] cpars{:}});
kk = kk.set_free ([1, 0, 1, 0, 1]);
kk = kk.set_bfun (@linear_bg, [0.1,0]);
kk = kk.set_bfree ([1,0]);
kk = kk.set_options ('list',2);

% Time a single iteration
tic
wsim = kk.simulate('comp');
t_spinw_single = toc;

% Time the fit
tic
[wfit, fitdata] = kk.fit('comp');
t_spinw = toc;

for i=1:numel(my_cuts)
    acolor black
    plot(my_cuts(i));
    acolor red
    pl(wfit.sum(i));
    pl(wfit.back(i));
    keep_figure;
end

% Display how long it takes to fit
if exist('t_ana', 'var');
    fprintf('Time for analytical fitting = %f s\n', t_ana);
end
fprintf('Time for SpinW single iteration = %f s\n', t_spinw_single);
fprintf('Time for SpinW fit = %f s\n', t_spinw);