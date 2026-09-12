# Independent source audit: T-REVERSED

This audit reconstructs the OLD sector from the original eight-page handwritten PDF and the supplied OLD transcription. Previously computed sectors are not sources of the physics. The PDF was reviewed visually, including the original lattice on page 1, the Majorana construction on pages 2-3, the real-space Hamiltonian on pages 3-4, the arrow diagram on page 4, the gauge and Fourier convention on page 5, the expanded Fourier Hamiltonian on pages 5-7, and the Bloch coefficients, transformation and energy on pages 7-8. Source notation is consistently rendered as `u` below.

Primary snapshot: `sources/original_notes.pdf`. Secondary unmodified transcription: `sources/OLD_reference_unchanged.md`. The transcription explicitly preserves the handwritten coefficients; it is not an independently corrected derivation.

The snapshots have the same SHA-256 as their supplied originals: PDF `e8d763f2b40c0df4cbd764b1d467d278e348076f589922752da54a5df5761e1e`; OLD transcription `3c8dec124881e58a66e5f86cb95b07124bd82130a6aea79159bddae8d38b6156`.

## 1. Geometry reconstructed before momentum coefficients

[DERIVADO DO PDF ORIGINAL] Pages 1, 3 and 4 identify a four-sublattice unit cell with translation `2 a1` and internal positions

\[
r_1=0,\quad r_2=a_1,\quad r_3=a_1+a_2,\quad r_4=a_2.
\]

Let `theta_m(n)` mean the Majorana on sublattice `m` at `2 n a1 + r_m`. For visual orientation one may draw `a1=(1,0)` and `a2=(1/2,1)`; these Cartesian drawing coordinates do not change the lattice vectors in the Fourier derivation. The original triangular bonds form alternating x and w blocks.

The following list is reconstructed from the real-space Hamiltonian, Eq. (11), on pages 3-4. `cell_shift` is the endpoint cell minus the starting cell in units of `2 a1`; `(da1,da2)` is the complete displacement of the endpoint relative to the start, INCLUDING internal sublattice positions. The two quantities must not be confused.

| Stored variable | Start | End | OLD u | cell_shift | da1 | da2 | Full displacement |
|---|---:|---:|---:|---:|---:|---:|---|
| u12x | 1 | 2 | +1 | 0 | +1 | 0 | a1 |
| u34x | 3 | 4 | +1 | 0 | -1 | 0 | -a1 |
| u23y | 2 | 3 | +1 | 0 | 0 | +1 | a2 |
| u41y | 4 | 1 | -1 | 0 | 0 | -1 | -a2 |
| u42z | 4 | 2 | +1 | 0 | +1 | -1 | a1-a2 |
| u31z | 3 | 1 | +1 | +1 | +1 | -1 | a1-a2 |
| u43w | 4 | 3 | +1 | -1 | -1 | 0 | -a1 |
| u12w | 1 | 2 | -1 | -1 | -1 | 0 | -a1 |
| u13t | 1 | 3 | -1 | 0 | +1 | +1 | a1+a2 |
| u24t | 2 | 4 | -1 | +1 | +1 | +1 | a1+a2 |

[TESTE INDEPENDENTE] Each row must satisfy

\[
2\,\mathrm{cell\_shift}\,a_1+r_{\rm end}-r_{\rm start}
=da1\,a_1+da2\,a_2.
\]

The potentially misleading real-space t term in Eq. (11) is `theta_2(R-a1) theta_4(R+a2)`. Its start-to-end displacement is `a1+a2`; shifting the summation origin gives the canonical-cell form `theta_2(n) theta_4(n+1)`. The z term `theta_3(R+a2-a1) theta_1(R)` likewise becomes `theta_3(n) theta_1(n+1)`.

## 2. Gauge checked against arrows and Eq. (15)

[DERIVADO DO PDF ORIGINAL] Page 5 states that the arrow `i -> j` means `u_ij,gamma=+1`, and reversing it gives `u_ji,gamma=-u_ij,gamma`. The ten signs printed in Eq. (15) agree with the ten OLD entries above and with the page-4 arrows: x arrows point `1 -> 2` and `3 -> 4`; y arrows point `2 -> 3` and `1 -> 4`; z arrows point `4 -> 2` and `3 -> 1`; w arrows point `4 -> 3` and `2 -> 1`; t arrows point `3 -> 1` and `4 -> 2`, with the cell translations in the table.

[DEDUÇÃO MATEMÁTICA] T-REVERSED changes only `u13t` and `u24t` from -1 to +1. Thus its t arrows point `1_n -> 3_n` and `2_n -> 4_(n+1)`. All other coordinates, cell translations and arrows remain unchanged.

## 3. Complete anticlockwise plaquette audit

[DERIVADO DO PDF ORIGINAL] Eq. (13), page 4, defines `W_p` as the product of directed u variables around a triangular plaquette anticlockwise. Eq. (14) defines `W_p=exp(i Phi_p)`. These are the source's real Z2 flux conventions. No extra phase factor is inserted here.

[DEDUÇÃO MATEMÁTICA] The x block consists of `(1_n,2_n,3_n,4_n)`. The adjacent w block consists, in anticlockwise quadrilateral order, of `(2_n,1_(n+1),4_(n+1),3_n)`. This cell information is essential: the labels alone do not identify the w triangles.

The following directed factors follow from traversing the geometric triangles, not from assigning target fluxes. A minus sign preceding a stored variable is precisely the application of antisymmetry. Every stated cycle has positive signed Cartesian area for `a1=(1,0), a2=(1/2,1)`.

| Plaquette | Exact anticlockwise path | First factor | Second factor | Third factor | OLD product | T product |
|---|---|---|---|---|---|---|
| W124 | 1_n -> 2_n -> 4_n -> 1_n | u12x, same | -u42z, reverse | u41y, same | (+1)(-1)(-1)=+1 | (+1)(-1)(-1)=+1 |
| W234 | 2_n -> 3_n -> 4_n -> 2_n | u23y, same | u34x, same | u42z, same | (+1)(+1)(+1)=+1 | (+1)(+1)(+1)=+1 |
| W123 | 1_n -> 2_n -> 3_n -> 1_n | u12x, same | u23y, same | -u13t, reverse | (+1)(+1)(+1)=+1 | (+1)(+1)(-1)=-1 |
| W134 | 1_n -> 3_n -> 4_n -> 1_n | u13t, same | u34x, same | u41y, same | (-1)(+1)(-1)=+1 | (+1)(+1)(-1)=-1 |
| W213 | 2_n -> 1_(n+1) -> 3_n -> 2_n | -u12w, reverse | -u31z, reverse | -u23y, reverse | (+1)(-1)(-1)=+1 | (+1)(-1)(-1)=+1 |
| W143 | 1_(n+1) -> 4_(n+1) -> 3_n -> 1_(n+1) | -u41y, reverse | u43w, same | u31z, same | (+1)(+1)(+1)=+1 | (+1)(+1)(+1)=+1 |
| W214 | 2_n -> 1_(n+1) -> 4_(n+1) -> 2_n | -u12w, reverse | -u41y, reverse | -u24t, reverse | (+1)(+1)(+1)=+1 | (+1)(+1)(-1)=-1 |
| W243 | 2_n -> 4_(n+1) -> 3_n -> 2_n | u24t, same | u43w, same | -u23y, reverse | (-1)(+1)(-1)=+1 | (+1)(+1)(-1)=-1 |

The products in this note are the mathematical consequences of the reconstructed paths; the production pipeline must independently evaluate them from its graph. T reverses four of the eight listed triangular fluxes: W123, W134, W214 and W243. OLD is zero flux under the source convention.

[DEDUÇÃO MATEMÁTICA] For an arbitrary site gauge `u_ij -> s_i u_ij s_j`, each site factor appears twice around a closed path, so every plaquette product is invariant. The changed plaquettes prove inequivalence even allowing cell-dependent site gauges. Searching the 16 cell-periodic choices is a second, narrower test. A simple contradiction is already present: unchanged x, y bonds imply `s1=s2=s3=s4`, which cannot reverse a t bond.

## 4. Real-space Hamiltonian and direct difference

[DERIVADO DO PDF ORIGINAL] With all bonds counted once, Eq. (11) and Eq. (15) yield

\[
\begin{aligned}
H_{\rm OLD}=i\sum_n\{&K_x[\theta_{1n}\theta_{2n}+\theta_{3n}\theta_{4n}]
+K_y[\theta_{2n}\theta_{3n}-\theta_{4n}\theta_{1n}]\\
&+K_z[\theta_{4n}\theta_{2n}+\theta_{3n}\theta_{1,n+1}]
+K_w[\theta_{4n}\theta_{3,n-1}-\theta_{1n}\theta_{2,n-1}]\\
&-K_t[\theta_{1n}\theta_{3n}+\theta_{2n}\theta_{4,n+1}]\}.
\end{aligned}
\]

The superscript y on every theta is suppressed only for legibility. T has precisely the same first two lines and a plus sign before the final `K_t` bracket. Direct subtraction in real space, before any coefficient formula, gives

\[
\Delta H_T=2iK_t\sum_n[\theta_{1n}\theta_{3n}+\theta_{2n}\theta_{4,n+1}].
\]

## 5. Fourier RAW and Eq. (27)

[DERIVADO DO PDF ORIGINAL] Page 5 uses `sqrt(4/N)` and the half-BZ Majorana Fourier sum with a positive phase multiplying the annihilation operator and its conjugate phase multiplying the creation operator. The subsequent expanded real-to-Fourier expression includes its factor 2 over the half BZ before being rewritten using a full BZ. Writing `q1=k.a1`, `q2=k.a2`, the single-bond Hermitian matrix contribution is

\[
h_{ij}\mathrel{+}=iK_\gamma u_{ij,\gamma}e^{i(q1\,da1+q2\,da2)},\qquad
h_{ji}\mathrel{+}=h_{ij}^{*}.
\]

[DEDUÇÃO MATEMÁTICA] The OLD upper entries reconstructed from the displacement table are

\[
\begin{aligned}
f_{12}^{\rm raw}&=i(K_xe^{iq1}-K_we^{-iq1}),\\
f_{13}^{\rm raw}=f_{24}^{\rm raw}&=-i e^{iq2}(K_te^{iq1}+K_ze^{-iq1}),\\
f_{14}^{\rm raw}=f_{23}^{\rm raw}&=iK_y e^{iq2},\\
f_{34}^{\rm raw}&=i(K_xe^{-iq1}-K_we^{iq1}).
\end{aligned}
\]

These signs also agree with the explicitly expanded Eq. (17), pages 5-7. In particular, the z phase contributing to f13 is `exp[-i(q1-q2)]`, not `exp[+i(q1-q2)]`. A high-resolution inspection of the relevant page-6 line confirms the minus sign before `i(q1-q2)`.

For T the only changed RAW entries are

\[
f_{13,T}^{\rm raw}=f_{24,T}^{\rm raw}
=i e^{iq2}(K_te^{iq1}-K_ze^{-iq1}).
\]

Every cross-chain UPPER entry has a common `exp(+i q2)`. This is fixed by endpoint displacements, independently of the later basis transformation.

[DERIVADO DO PDF ORIGINAL] Eq. (27), page 8, displays the diagonal phase matrix `U=diag(1,1,exp(-i q2),exp(-i q2))`. Its arrow notation does not independently name the two spinors as RAW and Bloch.

[DEDUÇÃO MATEMÁTICA] The project interpretation `Theta_raw=U Theta_Bloch` is the one consistent with Eq. (17); substitution gives `H_Bloch=U† H_raw U`. Since each cross-chain upper entry is multiplied on the right by `exp(-i q2)`, q2 cancels. U is unitary, and the lower entries cancel by conjugation. Applying `U H_raw U†` instead leaves a phase `exp(+2i q2)`; reversing both the RAW q2 sign and the transformation would accidentally cancel two errors. An independent RAW test is therefore necessary.

After the transformation, set `q=q1`; f12, f14, f23 and f34 remain the OLD expressions. Only f13 and f24 have their Kt term reversed.

## 6. Explicit source inconsistency: Eq. (26)

[DERIVADO DO PDF ORIGINAL] The literal page-7 Eq. (26), preserved in the OLD transcription, is

\[
f_{34}^{(26),\,\mathrm{literal}}=i(K_xe^{iq}-K_we^{-iq}).
\]

[DEDUÇÃO MATEMÁTICA] The real-space bonds `3 -> 4` and `4 -> 3` have displacements `-a1` in their respective stored orientations. Both the geometric reconstruction and Eq. (17) instead give

\[
f_{34}^{\mathrm{geometry}}=i(K_xe^{-iq}-K_we^{iq})=-f_{12}^{*}.
\]

Eq. (27) cannot change f34 because sublattices 3 and 4 carry the same diagonal phase. Thus this is an inconsistency between the literal Eq. (26) and the preceding geometric derivation, not a change introduced by T. The pipeline should implement the geometric expression and state the discrepancy; it must not call that expression the literal Eq. (26), overwrite the PDF, or claim the OLD transcription supplied the correction.

## 7. Majorana normalization and double coverage

[DEDUÇÃO MATEMÁTICA] Fix the convention `{theta_i,theta_j}=2 delta_ij`. For the once-counted real-space bond sum, define a real antisymmetric matrix A by `A_ij=2 K_gamma u_ij` for each stored oriented bond, accumulating repeated endpoints if necessary. Then

\[
H={i\over4}\theta^T A\theta,\qquad h={iA\over2}.
\]

The matrix h has real eigenvalues in opposite pairs globally. If a positive eigenvalue is `e`, the corresponding complex-fermion term is `e(2 f†f-1)`, so its vacuum contribution is `-e` and its occupation excitation costs `2e`. Consequently

\[
E_0=\sum_{e_j(h)<0}e_j(h)
\]

has NO additional one-half factor for this definition of h. A two-Majorana dimer `iK theta1 theta2` has exact ground energy `-|K|`, while h has eigenvalues `±|K|`; it is an independent calibration of the factor. The matrix `iA`, sometimes used elsewhere, would have eigenvalues doubled and would require the corresponding half factor. These matrix conventions must not be mixed.

For M physical four-site cells there are `N=4M` matter Majoranas and M independent four-by-four Fourier blocks over any width-pi interval in `q=k.a1` because the cell translation is `2a1`. The momenta satisfy `exp(i 2 M q)=1`. With

\[
F(q)=\sum_{n=1}^{4}E_n(q)\,\Theta[-E_n(q)],
\]

the finite unprojected matter vacuum energy per site is `sum_{M blocks} F(q)/(4M)`. The four-component representation makes the spectrum pi periodic, with `G_pi=diag(1,-1,-1,1)` implementing the pi shift; F has the same period. Consequently

\[
\epsilon={1\over4\pi}\int_{q_0}^{q_0+\pi}F(q)dq
={1\over8\pi}\int_{-\pi}^{\pi}F(q)dq.
\]

On a doubled grid with 2M samples over `[-pi,pi)`, each sample represents half of a physical block. Thus `sum_{2M} F/(8M)` equals `sum_M F/(4M)`. The plotting interval has two equivalent spectral copies, not twice as many physical states. The factor 1/4 divides by sites per cell; the half weight removes duplicated blocks and must not be mistaken for an additional Majorana vacuum factor.

No same-q opposite-energy pairing assumption is needed: global particle-hole symmetry relates q and -q, and summing all negative eigenvalues handles all cases. The finite-chain comparison must assemble its matrix directly from translated real bonds, not inverse transform a closed Bloch formula. Comparing those independently assembled matrices and an exact dimer/many-body calibration checks both displacement and absolute normalization.

## 8. Physical limits that must accompany the conclusion

[HIPÓTESE / NÃO PROVADO] The source derives a Majorana representation and fixes link eigenvalues, but these eight pages do not specify or carry out the local projection from the enlarged Majorana Hilbert space, the allowed matter fermion parity in each finite gauge sector, or global Wilson-loop/boundary sectors. A periodic-chain vacuum calculation therefore tests the unprojected quadratic matter problem. In a finite physical spin system an allowed-parity restriction can change the lowest state in a fixed link configuration; its subextensive correction must be distinguished from the thermodynamic energy density. Do not present an unprojected finite-chain equality as an exact proof of projected finite-spin energies.

The requested `zero_energy_band_gap=2 min_abs_E` is a precisely defined diagnostic. The matrix eigenvalues E use the h convention above; occupation excitation energies are 2|E|. If a physical constraint only permits multi-particle excitations, the physical excitation gap requires a separate calculation. Gap size alone does not determine energetic preference.

Changed plaquette fluxes prove OLD and T are gauge inequivalent. Equal spectra or equal integrated energies, if later found, would not undo this proof: isospectrality is weaker than gauge equivalence. Likewise a gapped spectrum, a flux label, and a comparison of several translation-preserving gauges do not prove a topological invariant or global stability against all gauges, enlarged cells, boundary sectors, interactions or perturbations. These remain open until independently computed.

The audit approves only the source reconstruction and the stated mathematical consequences. Numerical convergence, random validation, actual script hashes, plot inspection and the final physical preference must be established from fresh T computations.

## 9. Independent explicit many-body normalization calibration

[TESTE INDEPENDENTE] Executed with `OPENBLAS_NUM_THREADS=1 python3 -W error` on 2026-09-11. Twelve matter Majoranas for three cells were represented by six Pauli factors on a 64-dimensional complex-fermion Hilbert space:

\[
\theta_{2j}=Z_0\cdots Z_{j-1}X_j,\qquad
\theta_{2j+1}=Z_0\cdots Z_{j-1}Y_j.
\]

The maximum matrix error in `{theta_i,theta_j}=2 delta_ij` was exactly zero. The 64-by-64 operator H was assembled by multiplying the Majorana matrices in the ten real bonds, using the literal source table of starts, ends, signs and cell shifts, independently of `core.GEOMETRY` and all momentum formulas. Direct `eigvalsh(H)[0]/12` was compared with the production finite-chain `-sum(abs(eigvalsh(iA)))/(4*12)`. This tests the absolute factor through the explicit many-body operator, not only agreement between two quadratic one-particle matrices.

| Couplings (Kx,Ky,Kz,Kt,Kw) | Sector | Explicit many-body E0/N | Absolute difference from finite iA formula |
|---|---|---:|---:|
| (1,1,1,1,1) | OLD | -1.062967019772246 | 1.33e-15 |
| (1,1,1,1,1) | T-REVERSED | -1.048583770354864 | 6.66e-16 |
| (1.13,0.79,1.31,0.92,0.67) | OLD | -1.060108523583230 | 4.44e-16 |
| (1.13,0.79,1.31,0.92,0.67) | T-REVERSED | -1.024274743110539 | 2.22e-16 |
| (-1.21,0.71,1.41,-0.93,0.53) | OLD | -1.070036983742864 | 1.55e-15 |
| (-1.21,0.71,1.41,-0.93,0.53) | T-REVERSED | -0.902506507679779 | 2.22e-15 |

All six cases passed, including positive and mixed-sign couplings. This is the unprojected matter Hilbert space; it does not implement the separate physical spin projection. Even if the isotropic OLD and T thermodynamic integrals agree, their periodic finite-size energies need not be identical at the same number of cells, as the three-cell values explicitly demonstrate.
