%% Install SpinW
% Download the latest version of SpinW from https://github.com/SpinW/spinw/releases
% Or use the Add-ons option on the MATLAB toolbar.
% Unzip the file and run install_spinw from the main directory
%
% Verify that spinw is installed:


%% Help on SpinW
% Use the help command to get help on any SpinW function. For any function
% that starts with 'sw_' use help swfiles, for spinw class methods use help
% spinw.function_name. For help on plotting commands, use help swplot.

help swfiles

% To find the location of the spinw library use
sw_rootdir

% To open any of the functions in the Matlab editor, use
edit spinw.plot

% To look at any of the spinw object properties, double click on the
% "Workspace" view on the name of the object. Also the data files
% (symmetry.dat, atom.dat, color.dat, magion.dat) can be easily edited, for
% example

edit symmetry.dat

%% Lattice
% we create a triangular lattice

tri = spinw;

% give lattice parameters and angles in degree
tri.genlattice('lat_const',[3 3 4],'angled',[90 90 120])

% plot the lattice
plot(tri)

% In the plot window, you can zoom with the mouse wheel, pan by pressing
% the Ctrl button while dragging. Change the plot range and view direction
% by pressing the corresponding button on the top.


%% Atoms
% we add 1 atom at the origin of the unit cell, spin-3/2 Cr3+ ion

tri.addatom('r',[0 0 0],'S',3/2,'label','MCr3')

swpref.setpref('fontsize',12)
plot(tri)


%% Spin-Hamiltonian
% we create an antiferromagnetic first neighbor Hamiltonian
% plus easy plane single ion anisotropy
% red ellipsoids represent the single ion anistropy on the plot
% (equienergetic surface)
% examine the plot and test different values of A0 with different signs

A0 = -0.1;

tri.addmatrix('label','J1','value',1)
tri.addmatrix('label','A','value',[0 0 0;0 0 0;0 0 A0])


tri.gencoupling

tri.addcoupling('mat','J1','bond',1)
tri.addaniso('A')

% make good quality plot
swpref.setpref('nmesh',0,'npatch',10)
plot(tri,'range',[3 3 1/2],'cellMode','inside')


%% Magnetic structure
% the ground state magnetic structure of the above Hamltonian is a spiral,
% with propagation vector of (1/3,1/3,0). We define the plane of the
% spiral as the ab plane
% Careful: the given spin vector is column vector!
% What are the angles between neares neighbor moments?

tri.genmagstr('mode','helical','S',[1;0;0],'k',[1/3 1/3 0],'n',[0 0 1],'nExt',[1 1 1])

plot(tri,'range',[3 3 1/2],'cellMode','inside','magColor','red')


%% Spin wave dispersion
% we calculate the spin wave dispersion along the (H,H,0) high symmetry
% direction
% How many number of modes are there and why?
% What does the red line mean?
%
% Did you got any warning?
%

spec = tri.spinwave({[0 0 0] [1 1 0] 500},'hermit',false);

figure
sw_plotspec(spec,'mode','disp','imag',true,'colormap',[0 0 0],'colorbar',false)
axis([0 1 0 5])

%% Spin-Spin correlation functions
% the spin-spin correlations is already calculated, however it contains 9
% numbers per Q-point per mode. It is not possible to show this on a single
% plot.
%
% BUT!
%
% 1. we  can calculate the neutron scattering cross section
% 2. we can select one of the components S_alpha_beta(Q,w)
% 3. we can sum up the diagonal S_alpha_alpha(Q,w)
%
% What do you see? Where is the largest intensity? How is it related to the
% magnetic propagation vector?
%
% Why ar esome modes gapped? Which correlations are gapped?
%
% Why do we have Szz?
%
% Should not the spins precess around M, thus shouldn't each spin wave
% mode corresponds to a precession in a specific plane not a motion along a
% specific axis?
%

%spec = sw_egrid(spec,'component','Sperp','Evect',0:0.01:5.5);
spec = sw_egrid(spec,'component',{'Sxx+Syy' 'Szz'},'Evect',0:0.01:5);
%spec = sw_egrid(spec,'component','Syz','Evect',0:0.01:5);

figure
sw_plotspec(spec,'mode','color','dE',0.2,'imag',false)
%hold on
%sw_plotspec(spec,'mode','disp','color','k','linestyle','-')
axis([0 1 0 5.5])
caxis([0 3])


%% k=0 magnetic structure
% We can generate the same magnetic structure with zero propagation vector
% on a magnetic supercell.
%
% Why would we do that? Not for a triangular lattice, but in more complex
% cases, we will have to!
%
% You can keep the previous structure plot by pressing the red circle on
% the figure to compare the new magnetic structure. Is there any
% difference?
%
% What are we doing here? Can you tell from the script?
%
% check out the tri.magstr command? What data is stored there? What are the
% dimensions of the different matrices?
%
% You can also compare the energy per spin of the old magnetic structure
% and the new magnetic structure using the spinw.energy() function. Is
% there any difference?
%
% You can also store the new magnetic structure in a separate SpinW object
% by first duplicating the original object using the copy() command!
%
% What happens when we use tri2 = triNew, without the copy command?

triNew = copy(tri);
tri2   = triNew;

triNew.genmagstr('mode','random','next',[3 3 1])
triNew.optmagsteep('nRun',1e4)
% Converged?

triNew.genmagstr('mode','rotate','n',[0 0 1])

phi1 = atan2(triNew.magstr.S(2,1),triNew.magstr.S(1,1));

triNew.genmagstr('mode','rotate','n',[0 0 1],'phi',-phi1)

plot(triNew,'range',[3 3 1])

% How does the magnetic structures compare? Are they the same? Why?
%

tri.magstr
tri2.magstr

%% Spin wave dispersion on the k=0 magnetic structure

% we calculate the spin wave dispersion along the (H,H,0) high symmetry
% direction
% How many modes do we have? Is there more than before?
% Which is the right one then?
%
% Why are there vertical lines in the dispersion? Is it a bug?
%

spec = tri.spinwave({[0 0 0] [1 1 0] 500},'hermit',false);

figure
subplot(2,1,1)
sw_plotspec(spec,'mode','disp','imag',true,'colormap',[0 0 0])
colorbar off
axis([0 1 0 5])

% spin-spin correlation functions
spec = sw_egrid(spec,'component','Sperp','Evect',0:0.01:5.5);

subplot(2,1,2)
sw_plotspec(spec,'mode','color','dE',0.2,'imag',false)
axis([0 1 0 5.5])
caxis([0 3])
colormap jet


