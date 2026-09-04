# Symmetry-adapted shared displacement-response ensemble

Development derivation. Numerical validation status is tracked separately in
IMPLEMENTATION_PLAN.md. This construction combines established response
definitions and finite-displacement symmetry methods; it is not a claim of a
new physical definition of Born charges or force constants.

## 1. Common derivatives

Let kappa label atoms, alpha and beta Cartesian components, u atomic
displacements, and E the macroscopic electric field. In an insulating system
with a fixed cell and specified electrical boundary conditions, expand the
electric enthalpy near the reference configuration:

\[
\mathcal F(\mathbf u,\mathbf E)=\mathcal F_0
-\mathbf F_0^T\mathbf u+\tfrac12\mathbf u^T\Phi\mathbf u
-\boldsymbol\mu_0^T\mathbf E
-e\mathbf E^T\mathcal Z\mathbf u
-\tfrac12\mathbf E^T\alpha^{\rm el}\mathbf E+\cdots .
\]

Here \(\mathcal Z=(Z_1^*,\ldots,Z_N^*)\), and

\[
Z^*_{\kappa,\alpha\beta}
=\frac1e\frac{\partial\mu_\alpha}{\partial u_{\kappa\beta}}
=\frac1e\frac{\partial F_{\kappa\beta}}{\partial E_\alpha},
\qquad
\Phi_{\lambda\alpha,\kappa\beta}
=-\frac{\partial F_{\lambda\alpha}}{\partial u_{\kappa\beta}}.
\]

For a periodic bulk, replace dipole derivatives by
\(\partial\boldsymbol\mu=\Omega\partial\mathbf P\), with a continuous
Berry-polarization branch. The reference need not have exactly zero residual
forces for the derivative to exist: explicitly subtract \(\mathbf F_0\).
However, a phonon stability or equilibrium dielectric interpretation requires
a sufficiently relaxed reference, and a converged SCF is not proof of ionic
convergence.

One displaced SCF therefore supplies two compatible observations:

\[
\Delta\boldsymbol\mu/e=Z_\kappa^*\mathbf d+O(|\mathbf d|^2),
\qquad
-\Delta\mathbf F_\lambda=\Phi_{\lambda\kappa}\mathbf d
+O(|\mathbf d|^2).
\]

Force evaluation must be enabled in that SCF. Polarization still needs its
Berry or real-space postprocessing. Enabling force output does not itself
compute polarization, the electronic dielectric tensor, or Raman tensors.

## 2. Symmetry covariance

Let g be a symmetry of the reference Hamiltonian and boundary conditions.
Its orthogonal Cartesian matrix R_g and atomic permutation p_g satisfy

\[
Z^*_{p_g(\kappa)}=R_g Z^*_\kappa R_g^T,
\qquad
\Phi_{p_g(\lambda),p_g(\kappa)}
=R_g\Phi_{\lambda\kappa}R_g^T.
\]

For a fractional-coordinate lattice operation W_g and row-vector lattice L,
\(R_g=L^T W_g L^{-T}\). In a nonorthogonal cell W_g is not a Cartesian
rotation. Every atom map is checked for species consistency, uniqueness,
and a stated positional tolerance. Approximate symmetry must not be accepted
with a substantially nonorthogonal R_g.

For a site stabilizer, p_g(kappa)=kappa. A measured displacement/response pair
then supplies the transformed pair

\[
(\mathbf d,\Delta\boldsymbol\mu,\Delta\mathbf F_\lambda)
\mapsto(R_g\mathbf d,R_g\Delta\boldsymbol\mu,
R_g\Delta\mathbf F_\lambda),
\]

with force atom lambda relabeled p_g(lambda). The permutation is essential;
rotating force components without moving the atom labels is incorrect.
Atomic BEC tensors are not generally symmetric in their Cartesian indices.
Do not impose \(Z^*=Z^{*T}\).

## 3. Identifiability and general symmetry

For each inequivalent atom, collect seed displacements and all their site
stabilizer images as columns of X_kappa. Collect the corresponding dipole and
negative-force changes as columns of Y_kappa. Solve

\[
Y_\kappa=B_\kappa X_\kappa,\qquad
B_\kappa=Y_\kappa X_\kappa^+,
\qquad \operatorname{rank}X_\kappa=3.
\]

The top three rows of B give Z_kappa; successive three-row force blocks give
Phi_lambda,kappa. The same geometric rank criterion therefore supplies all
three displacement columns for both responses. Equivalent atoms are obtained
by covariance and averaged over equivalent mappings. Singular values,
condition numbers, and fitting residuals are recorded.

This is not a universal claim that two SCFs suffice. In a site with trivial
symmetry, three independent directions are necessary; a central scheme uses
both signs. In higher symmetry, one mixed direction can generate a spanning
orbit. Distinct inequivalent atoms need distinct seed sets. Phonopy supplies
these geometry-dependent seeds, and ZStar verifies their actual ranks rather
than assuming a fixed reduction factor. Additional physical constraints can
sometimes reduce the required observations further; the present scheme does
not use them to rescue a rank-deficient data set.

For a finite crystallographic symmetry group this follows directly from
linear algebra, independent of the space-group label. Validation over 32
crystallographic point groups is an analytic test of the reconstruction,
not numerical DFT validation of every possible material, magnetic group,
molecular continuous group, or electrical boundary condition.

## 4. Step length and truncation error

The ABACUS Phonopy default step of 0.02 bohr is
0.01058354421806 Angstrom. It can be described as approximately 0.01 Angstrom
in prose, but must not be replaced by 0.01 in a derivative denominator.
The implementation reads both serialized STRU files, computes their actual
minimum-image Cartesian difference, verifies that exactly one atom moved,
and uses the entire displacement vector in X. No individual Cartesian
component is treated as the total displacement length.

With an explicit opposite displacement, an odd difference cancels the
quadratic response term and yields an O(delta^2) derivative error. If a site
symmetry maps a seed onto its negative, the transformed response supplies the
opposite observation under the same Hamiltonian, so a second physical SCF is
unnecessary. Phonopy's automatic +/- selection retains explicit negatives
when required. Explicit forward-only sampling need not be second order.
Electronic convergence, force noise, real-space integration, Berry branch
selection, and printed precision remain independent sources of error.

The spatial group must also be respected to numerical accuracy by the response
evaluation. A small spurious transverse polarization induced by a mostly
normal displacement can contaminate the inferred in-plane derivative when
the in-plane seed component is small. Compare mixed and Cartesian controls,
converge the Berry integration mesh, and inspect symmetry-forbidden Cartesian
responses rather than assuming that full printed precision ensures physical
accuracy. A smaller displacement alone need not remove this sampling error.

## 5. Gamma scope and dimensions

A primitive-cell displacement repeats in every lattice cell. Its force
derivative is the lattice-summed Gamma matrix

\[
\Phi^\Gamma_{\lambda\alpha,\kappa\beta}
=\sum_{\mathbf R}\Phi_{0\lambda\alpha,\mathbf R\kappa\beta}.
\]

An identity-supercell shared ensemble yields this Gamma response, not the
individual real-space couplings needed for a general phonon dispersion.
A conventional input cell instead gives the corresponding folded Gamma
matrix; it is not silently interpreted as a primitive cell. Supercell and
finite-q reuse require separate mappings and convergence studies.

For dim=3, use Berry polarization in all directions. For dim=2, use in-plane
Berry changes and the localized, charge-neutral slab dipole normal to the
layer. For dim=1, use the periodic-direction Berry change and two transverse
localized dipoles. For dim=0, use a molecular dipole derivative/APT or a
validated large-cell Berry equivalent. Do not reduce symmetries that mix
periodic and open directions or are broken by fields, magnetic order, or
the electrostatic treatment. Current shared code restricts these unsupported
cases explicitly and uses a valid periodic-supercell subgroup for molecules,
which need not achieve the maximal molecular reduction.

Dipole correction and vacuum separation are part of the physical setup, not
optional changes made between the BEC and force-constant calculations.
Localized neutral real-space dipoles are not arbitrarily wrapped by a
bulk polarization quantum. Charged, delocalized, or boundary-crossing systems
need additional validation of the real-space integral.

## 6. Constraints and dielectric closure

Preserve raw tensors and report charge neutrality, force translational sums,
Hessian reciprocity, residual reference forces, and fit residuals before
applying any projection. For a neutral system, the explicit projected outputs
obey \(\sum_\kappa Z_\kappa^*=0\), Hessian reciprocity, and translational
invariance. Charge corrections are weighted by the actual full atom count,
not the number of inequivalent representatives. Molecular rotational
constraints are separate from this translational projection.

With stable optical modes, equilibrium displacement minimizes the quadratic
enthalpy and gives the phonon response

\[
\epsilon^{\rm ph}(0)
=\frac{e^2}{\epsilon_0\Omega}\mathcal Z
(\Phi^\Gamma)^+\mathcal Z^T,
\qquad
\epsilon(0)=\epsilon^\infty+\epsilon^{\rm ph}(0).
\]

The inverse excludes rigid translations and any other genuine zero modes.
Negative optical eigenvalues do not define a stable static susceptibility;
do not hide them by dropping modes or taking absolute frequencies. In SI,
convert the Hessian from eV/Angstrom^2 before evaluating this expression.
For low-dimensional systems, report the appropriately normalized molecular,
line, or sheet polarizability instead of interpreting a vacuum-dependent
supercell dielectric constant as an intrinsic bulk material quantity.

Benchmark the entire chain: raw/constrained Z, raw/constrained Phi, matched
mode frequencies (degenerate subspaces rather than arbitrary eigenvectors),
mode-effective-charge strengths, and static response where stable. Quote
actual allocated CPU core-hours and the SCF and PYATB components separately.
Displacement count reduction alone is not a measured timing speedup.

## References

- X. Gonze and C. Lee, Phys. Rev. B 55, 10355 (1997),
  https://doi.org/10.1103/PhysRevB.55.10355.
- Phonopy finite-displacement method and interfaces:
  https://phonopy.github.io/phonopy/formulation.html and
  https://phonopy.github.io/phonopy/abacus.html.
- W. Ding et al., Nature Communications 8, 14956 (2017),
  https://doi.org/10.1038/ncomms14956. New In2Se3 calculations match the stated
  PBE/k-point/vacuum/dipole/force settings, not the PAW/LCAO numerical basis.
