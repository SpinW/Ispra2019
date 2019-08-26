% SpinW Magnetic Structure tutorial 1, Horace/SpinW Workshop 2019
%
% The aim of this tutorial is to calculate the propagation vector
% of the cycloid of the multiferroic compound BiFeO3.
%
% An experimental description of the Hamiltonian may be found in
% the paper:
% Magnetic Dispersion and Anisotropy in Multiferroic BiFeO3,
%   M. Matsuda et al., Phys. Rev. Lett. 109 067205 (2012)
%   https://web.ornl.gov/~okapon/manybodytheory/publication/pdf/PhysRevLett_109_067205.pdf
%
% BiFeO3 is a famous room-temperature multiferroic. It is a slightly
% distorted perovskite structure. Its magnetic structure is a very
% long period cycloid which is stabilised by a Dzyaloshinskii-Moriya
% interaction. In the absence of the DM interaction, it would have
% a simple G-type AFM structure.
%
% We will apply the DM interaction and use SpinW to calculate the
% propagation vector k of the cycloid, using optmagk.

% First we define the model. The crystal structure may be obtained
% from ICSD, e.g. http://icsd.kisti.re.kr/icsd/icsd_view1.jsp?num=4067

bfo = spinw;
bfo.genlattice('lat_const', [5.58 5.58 13.86], 'angled', [90 90 120], 'spgr', 'R 3 c');
bfo.addatom('r', [0 0 0.2212], 'S', 2.5, 'label', 'MFe3', 'color', 'gold');

% The PRL paper (and most spin wave models of BiFeO3) considers two
% exchange interactions, a DM interaction and an axial single-ion
% anisotropy. We will use the parameters in the PRL paper.
%
% Note that the spin length used by SpinW is g*sqrt(S*(S+1)) where S
% is defined in the addatom function. We have used the free-ion
% ("theoretical") value of S=5/2 for Fe3+ here, whereas Matsuda
% et al. have used the measured moment at 200K which is 4.1uB
% Compared to the theoretical spin length this is smaller by 1.44
% as the mention in page 4 of the paper. Thus the exchange parameters
% should be scaled down by this amount

J1 = 4.5;  % == 6.48 / 1.44
J2 = 0.2;  % == 0.29 / 1.44
D = 0.1623;
K = -0.0068; % negative to indicate axial SIA rather than planar

bfo.gencoupling('maxDistance', 20)
bfo.addmatrix('label', 'J1', 'value', J1, 'color', 'green')
bfo.addmatrix('label', 'J2', 'value', J2, 'color', 'white')
bfo.addcoupling('mat', 'J1', 'bond', 1)
bfo.addcoupling('mat', 'J2', 'bond', 2)

% Plots the unit cell - check it resembles fig 1 of paper
plot(bfo)

% Now we add the DM term. SpinW has a short cut where if you give it
% a 3-vector, it would interpret this as a DM vector and construct the
% correct 3x3 exchange tensor.
% In the case of BiFeO3, the DM vector is along [1 -1 0] and connects
% atoms along the spiral direction [1 1 0] (e.g. 2nd neighbours).
%
% Note that this is superficially different to what the paper has
% (it says the spiral is along [1 -1 0]) because the paper uses the
% pseudo-cubic (rhombohedral) unit cell, whereas in SpinW, 'R 3 c'
% defaults to a hexagonal cell. The principal axes in the two cells
% are: [1 1 0]hex = [1 -1 0]pc, and [0 0 1]hex = [1 1 1]pc.

bfo.addmatrix('label', 'D', 'value', [1 -1 0]*D, 'color', [255 225 175])
bfo.addcoupling('mat', 'D', 'bond', 2)

% Print the exchange tensor for the DM interaction
id = find(~cellfun(@isempty, strfind(bfo.matrix.label, 'D')));
bfo.matrix.mat(:,:,id)

% Finally add the single-ion anisotropy
bfo.addmatrix('label', 'K', 'value', diag([0 0 K]), 'color', [185 135 0])
bfo.addaniso('K')

% Now optimise the propagation vector.
% We know that it should be qm = [delta delta 0] so we give optmagk
% a helping hand by setting this as a basis vector
res = bfo.optmagk('kbase', [1; 1; 0])
res.k

% What happens if we don't set the basis vector as [1; 1; 0]?
% E.g. if we used:
%res = bfo.optmagk()   % without giving 'kbase' option

% What about the energy of this ground state compared to the one with
% kbase=[1;1;0]? (Note, use:
%format long g
% To print more that 8 significant figures).

% What happens if we increase the value of D? (Or decrease J?)

%%
% You should find k = [0.0036 0.0036 0] or [0.9964 0.9964 0] because
% there is no energetic difference between k and 1-k.
% However, delta=0.0036 is quite different from delta=0.0045 that we
% were expecting.
%
% This is because in the harmonic approximation of treating the
% spin cycloid structure as S_i = S0.exp(ik.r_i), the propagation
% vector is independent of the single-ion anisotropy K
% (see "Origin of the long period magnetic ordering in BiFeO3",
% I. Sosnowska and AK. Zvezdin, J. Magn. Magn. Mat. 140-144
% (1995) 167-168, which also shows that for J2=0, the cycloid
% wavelength = 4J/D, or k = (D/J) / (2a) where a=5.58A
% is the lattice parameter, since cycloid wavelength = (a/2)/k.)
%
% The effect of the single-ion anisotropy is to make the cycloid
% anharmonic, so that S_i = S0.exp(i phi_i) where
% phi_i = k.r_i + psi_i, and psi_i is proportional to sin(2k.r_i).
%
% This anharmonicity results in higher order satelite peaks in
% the neutron diffraction, and more importantly makes the magnon
% modes dipole active - that is, makes them electromagnons.
% However, including these effects in SpinW is not currently
% possible, except through the use of a large supercell
% (as in the Matsuda paper). The computational cost, however,
% for BiFeO3 is prohibitive due to the long wavelength of the
% cycloid (the calculations in the Matsuda paper used an optimised
% code on a cluster).
%
% Calculations with SpinW are thus restricted to the K=0 limit
% in figure 5, e.g.

J1 = 4.5;  % == 6.48 / 1.44
J2 = 0.2;  % == 0.29 / 1.44
D = 0.185;
K = 0;

% Substitute in new parameter values
bfo.matparser('mat', {'J1', 'J2', 'K(3,3)'}, 'param', [J1 J2 K]);
% For the DM vector we have to put in the matrix explicitly
dmmat = [0 0 1; 0 0 1; -1 -1 0];
bfo.matparser('mat', {'D'}, 'selector', dmmat, 'param', [D]);

% Optimise the propagation vector again
res = bfo.optmagk('kbase', [1; 1; 0])
res.k

plot(bfo)

% Calculate the spin wave dispersion around the propagation vector
% with the ('hermit', false) option and see where the imaginary
% eigenvalues are (you will have to zoom in to below 2meV).

% Why are there still imaginary energies?

% Print the magnetic moment directions after the optmagk run:
%bfo.magstr.S

% Now run bfo.optmagsteep print the moments again.
% What has changed? Why?

% Now plot the spin wave dispersion again. Verify the imaginary
% energies are gone.

%%
% You can also try to verify the theoretical dependence of k on
% J1 and D

J1 = 0.1:0.1:10;
J2 = 0;
D = 0.2;
K = 0;
km = J1*0;
for ii = 1:numel(J1)
    bfo.matparser('mat', {'J1', 'J2', 'K(3,3)'}, 'param', [J1(ii) J2 K]);
    bfo.matparser('mat', {'D'}, 'selector', dmmat, 'param', [D]);
    res = bfo.optmagk('kbase', [1; 1; 0]);
    km(ii) = min([res.k(1) 1-res.k(1)]);
end
figure; hold all;
plot(J1, km, '.');
plot(J1, (D./J1)./(2*5.58), '-');
xlabel('J_1 (meV)'); ylabel('k (r.l.u.)');

J1 = 4;
J2 = 0;
D = 0.05:0.05:5;
K = 0;
km = D*0;
for ii = 1:numel(D)
    bfo.matparser('mat', {'J1', 'J2', 'K(3,3)'}, 'param', [J1 J2 K]);
    bfo.matparser('mat', {'D'}, 'selector', dmmat, 'param', [D(ii)]);
    res = bfo.optmagk('kbase', [1; 1; 0]);
    km(ii) = min([res.k(1) 1-res.k(1)]);
end
figure; hold all;
plot(D, km, '.');
plot(D, (D./J1)./(2*5.58), '-');
xlabel('D (meV)'); ylabel('k (r.l.u.)');

% Is the theoretical expression matched by the calculation?

% If not, why not?