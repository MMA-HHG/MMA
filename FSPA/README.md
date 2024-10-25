# Full-Saddle-Point-Approximation dipoles

## Theoretical background

This routine provides the dipoles,
$$ d_{\text{S-P}}(\omega) \approx - \frac{i}{\sqrt{2\pi}} \left( \frac{-2\pi i}{\sqrt{\det \Phi''_\omega}} \right)^{\frac{5}{2}} \left[ \mathcal{E}(t_i) \cdot \mathbf{d}\left( \mathbf{k}_{tr,t_i}^{(0)} + \mathbf{A}(t_i) \right) e^{i \Phi_\omega (t_r,t_i, \mathbf{k}_{tr,t_i}^{(0)})} \right] \mathbf{d}^*\left( \mathbf{k}_{tr,t_i}^{(0)} + \mathbf{A}(t_r) \right) + c.c. \,,
 $$
 obtained by computing the saddle points in all the integration variables in the [SFA integral approach](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.49.2117) for a driving field of the form $A(t)=A_0\cos(\omega_0t)$, with the electric field $\mathcal{E}(t) = -\partial_t A(t) = A_0 \omega_0 \sin(\omega_0 t)$. The complex saddle points times of ionisation $t_i$ and recombination $t_r$ are obtained by the Newton-Raphson method applied to the set of equations
 $$\frac{[\mathbf{k}_s + \mathbf{A}(t_i)]^2}{2} + I_p = 0, \quad
\int_{t_i}^{t_r} [\mathbf{k}_s + \mathbf{A}(t')] dt' = 0, \quad
\frac{[\mathbf{k}_s + \mathbf{A}(t_r)]^2}{2} + I_p = \omega \,. $$
([See Eqs. (1.29-1.31) in this reference for more details](https://arxiv.org/pdf/1304.2413).) Note that saddles in the momentum space $\mathbf{k}_s$ are obtained analytically $\mathbf{k}_{t, t'}^{(\text{sp})} = - \int_{t'}^{t} \mathbf{A}(t'') \, dt''/(t - t')$ for the given driver's field $A(t)$. Effectivelly, only the equations for $t_i$ and $t_r$ remain. In order to compute the response for a given harmonic $H$, we set $\omega = H\omega_0$. The model of the transition-dipole-elemnts $\mathbf{d}$ follows the model of truncated-harmonic potential given by [Eq. (21) in Lewenstein *et al.* (1994)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.49.2117) (note that $\mathbf{d}$ does not affect the saddle points).

The routine provides the solution for short and long of trajectories determined by the initial conditions $(t_i,t_r)_{\text{short}} = ((0.3+2\mathrm{i})T_{0}, (0.44 + 2\mathrm{i}) T_0)$, $(t_i,t_r)_{\text{long}} = ((0.15+2\mathrm{i})T_{0}, (0.75 + 2\mathrm{i}) T_0)$, where $T_0 = 2\pi/\omega_0$. The solution is problematic in the cut-off, because the [uniform approximation](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.66.043413) is not applied. This usually leads to an unphysical branch of the solution for the short trajectories. A crude estimation is to interpolate this region to link the plateau short branch with the above-cut-off long branch.

## Code and compilation
The code consists from two simple files `main.f90` and `tools.f90` it does not depend on any external libraries and is directly compiled by the script `compile.sh`.

The inputs are in the in the form of a text file `param.inp` consting of (all in atomic units):
* *the ionization potential of the target* $I_p$,
* *the fundamental frequency of the driving field* $\omega_0$,
* *the harmonic order of the interest* $H$,
* *the minimal value of the vector potential $A_{0,\mathrm{min}}$*,
* *the discretisation in the vector potential $\Delta A$*,
* *the number of steps in the vector potential*.

The outputs are `times_long.dat`, `times_short.dat`, `phase_long.dat`, `phase_short.dat`.

The file `times_xxx.dat` contains:\
$I[\mathrm{a.u.}]$ | $\Re(t_i)$ | $\Im(t_i)$  | $\Re(t_r)$ | $\Im(t_r)$  | $\Re(k_s)$ | $\Im(k_s)$  \
where $I[\mathrm{a.u.}] = (\omega_0 A_0)^2$

The file `phase_xxx.dat` contains:\
$I[\mathrm{a.u.}]$ | $\Re(\Phi_{\omega})$ | $\Im(\Phi_{\omega})$ | $\Re((d_{\text{S-P}}))$ | $\Im((d_{\text{S-P}}))$ | $\mathrm{e}^{-\Im(\mathrm{Arg}(d_{\text{S-P}}))}$ | $|d_{\text{S-P}}|$ | $\mathrm{Arg}(d_{\text{S-P}})$ |