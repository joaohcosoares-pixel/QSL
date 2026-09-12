# Quadrupolar Spin Liquid in 1D — OLD reference configuration

> **Purpose.** This file establishes the **OLD/original gauge sector** directly from the handwritten PDF.  
> No algebraic correction has been applied to the equations below. Signs, phases, matrix entries and equation structure are kept as in the source.  
> The bond variable is written as \(u_{ij,\gamma}\) to match the notation used with the advisor; this is only a typographical relabeling of the same bond variable.

---

## 1. Original gauge convention

The arrows define the bond variable by

$$
u_{ij,\gamma}=+1
$$

when the arrow points from \(i\) to \(j\), with

$$
u_{ij,\gamma}=-u_{ji,\gamma}.
$$

The original four-sublattice gauge configuration is

$$
\boxed{
\begin{aligned}
u_{12,x}&=+1, &
u_{23,y}&=+1,\\
u_{34,x}&=+1, &
u_{41,y}&=-1,\\
u_{42,z}&=+1, &
u_{43,w}&=+1,\\
u_{31,z}&=+1, &
u_{12,w}&=-1,\\
u_{13,t}&=-1, &
u_{24,t}&=-1.
\end{aligned}}
$$

This is the configuration referred to throughout this file as

$$
\boxed{\text{OLD}}
$$

or the original gauge configuration.

---

## 2. Gauge flux definition

For a triangular plaquette \(p\),

$$
W_p:=\prod_{\langle ij\rangle_\gamma\in p}u_{ij,\gamma},
$$

where the product is taken by circulating the plaquette in the anticlockwise direction.

The gauge flux is defined by

$$
W_p=e^{i\Phi_p}.
$$

The original configuration shown in the PDF is introduced as the free-gauge-flux sector,

$$
\Phi_p=0.
$$

---

## 3. Majorana Hamiltonian for a fixed gauge

For a fixed configuration of the bond variables,

$$
\boxed{
\hat H
=
i
\sum_{\gamma\in\{x,y,z,t,w\}}
\sum_{\langle ij\rangle_\gamma}
K_q^\gamma
u_{\langle ij\rangle_\gamma}
\hat\theta_i^y
\hat\theta_j^y
}
$$

with \(u_{\langle ij\rangle_\gamma}=\pm1\).

The four sublattices are placed as

$$
1:\mathbf R,
\qquad
2:\mathbf R+\mathbf a_1,
\qquad
4:\mathbf R+\mathbf a_2,
\qquad
3:\mathbf R+\mathbf a_1+\mathbf a_2.
$$

---

## 4. Fourier transformation

The Fourier convention used in the original notes is

$$
\boxed{
\hat\theta_m^y(\mathbf R)
=
\sqrt{\frac{4}{N}}
\sum_{\mathbf k\in \mathrm{BZ}/2}
\left[
e^{i\mathbf k\cdot\mathbf R}
\hat\theta_m^y(\mathbf k)
+
e^{-i\mathbf k\cdot\mathbf R}
\hat\theta_m^{y\dagger}(\mathbf k)
\right].
}
$$

The notes also impose

$$
\hat{\mathbf k}=k\hat{\mathbf x}.
$$

---

## 5. Bloch spinor

The spinor is

$$
\boxed{
\hat\Theta(\mathbf k)
=
\left(
\hat\theta_1^y(\mathbf k),
\hat\theta_2^y(\mathbf k),
\hat\theta_3^y(\mathbf k),
\hat\theta_4^y(\mathbf k)
\right)^T.
}
$$

The Hamiltonian is written as

$$
\boxed{
\hat H_I
=
\sum_{\mathbf k\in\mathrm{BZ}}
\hat\Theta^\dagger(\mathbf k)
\hat{\mathcal H}_I(\mathbf k)
\hat\Theta(\mathbf k).
}
$$

with

$$
\boxed{
\hat{\mathcal H}_I(\mathbf k)
=
\begin{pmatrix}
0 & f_{12}(\mathbf k) & f_{13}(\mathbf k) & f_{14}(\mathbf k)\\
f_{12}^*(\mathbf k) & 0 & f_{23}(\mathbf k) & f_{24}(\mathbf k)\\
f_{13}^*(\mathbf k) & f_{23}^*(\mathbf k) & 0 & f_{34}(\mathbf k)\\
f_{14}^*(\mathbf k) & f_{24}^*(\mathbf k) & f_{34}^*(\mathbf k) & 0
\end{pmatrix}.
}
$$

---

## 6. OLD Bloch coefficients — transcribed without algebraic correction

The coefficients written in Eqs. (21)–(26) of the original PDF are

$$
\boxed{
f_{12}(\mathbf k)
=
i
\left(
K_q^x e^{i\mathbf k\cdot\mathbf a_1}
-
K_q^w e^{-i\mathbf k\cdot\mathbf a_1}
\right)
}
\tag{21}
$$

$$
\boxed{
f_{13}(\mathbf k)
=
-i
\left(
K_q^t e^{i\mathbf k\cdot\mathbf a_1}
+
K_q^z e^{-i\mathbf k\cdot\mathbf a_1}
\right)
}
\tag{22}
$$

$$
\boxed{
f_{14}(\mathbf k)=iK_q^y
}
\tag{23}
$$

$$
\boxed{
f_{23}(\mathbf k)=iK_q^y
}
\tag{24}
$$

$$
\boxed{
f_{24}(\mathbf k)
=
-i
\left(
K_q^t e^{i\mathbf k\cdot\mathbf a_1}
+
K_q^z e^{-i\mathbf k\cdot\mathbf a_1}
\right)
}
\tag{25}
$$

$$
\boxed{
f_{34}(\mathbf k)
=
i
\left(
K_q^x e^{i\mathbf k\cdot\mathbf a_1}
-
K_q^w e^{-i\mathbf k\cdot\mathbf a_1}
\right)
}
\tag{26}
$$

> **Important:** Eq. (26) above is intentionally reproduced as written in the PDF.  
> This OLD reference file does **not** apply the later audit/correction of that equation.

---

## 7. Gauge/Bloch transformation of Eq. (27)

The transformation displayed in the original notes is

$$
\boxed{
\begin{pmatrix}
\hat\theta_1^y(\mathbf k)\\
\hat\theta_2^y(\mathbf k)\\
\hat\theta_3^y(\mathbf k)\\
\hat\theta_4^y(\mathbf k)
\end{pmatrix}
\longmapsto
\begin{pmatrix}
1&0&0&0\\
0&1&0&0\\
0&0&e^{-i\mathbf k\cdot\mathbf a_2}&0\\
0&0&0&e^{-i\mathbf k\cdot\mathbf a_2}
\end{pmatrix}
\begin{pmatrix}
\hat\theta_1^y(\mathbf k)\\
\hat\theta_2^y(\mathbf k)\\
\hat\theta_3^y(\mathbf k)\\
\hat\theta_4^y(\mathbf k)
\end{pmatrix}.
}
\tag{27}
$$

No additional reinterpretation of Eq. (27) is introduced in this reference file.

---

## 8. Energy dispersions and ground-state energy

Diagonalizing \(\hat{\mathcal H}_I(\mathbf k)\) gives

$$
E_{I,n}(\mathbf k),
\qquad
n=1,2,3,4.
$$

The ground-state-energy expression written in the notes is

$$
\boxed{
E_I^{(0)}
=
\frac{1}{N}
\sum_{n=1}^{4}
\sum_{\mathbf k\in \mathrm{BZ}}
E_{I,n}(\mathbf k)\,
\Theta\!\left[-E_{I,n}(\mathbf k)\right].
}
\tag{28}
$$

where \(\Theta(x)\) is the Heaviside step function.

This quantity is used to determine which gauge-field configuration has the lowest ground-state energy.

---

# OLD reference summary

The OLD sector used as the reference configuration is therefore

$$
\boxed{
\begin{aligned}
u_{12,x}&=+1, &
u_{34,x}&=+1,\\
u_{23,y}&=+1, &
u_{41,y}&=-1,\\
u_{42,z}&=+1, &
u_{31,z}&=+1,\\
u_{43,w}&=+1, &
u_{12,w}&=-1,\\
u_{13,t}&=-1, &
u_{24,t}&=-1.
\end{aligned}}
$$

and the equations in this document are preserved as they appear in the supplied handwritten source.

**This file is a source-faithful OLD reference, not an audited/corrected derivation.**
