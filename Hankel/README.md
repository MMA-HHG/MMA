# Hankel transform module
The modul `Hankel_transform` provides the toolbox for computing Hankel transform in both short and long media. The main two computational routines, which may stand independently, are `HankelTransform` and `Hankel_long`, repsectively.  Addressing the long media, the module sources material constants for rare gases. It means that they are ready for a direct use.

This documentation provides a transition between the mathematical formulation in the main manuscript and the Pythonic implementation. The manuscript provides details about theoretical background, while technical details about the implementation are directly in the documentation within the module.

## Main ideas

The integral to compute is
$$
    \hat{\mathcal{E}}(z,\rho,\omega)
=
    \int_0^L
    \int_0^{\tilde{\rho}_{\text{max}}}

         \underbrace{
            \varrho(\tilde{z},\tilde{\rho})
            \mathrm{e}^{\mathbf{i} \Phi(\tilde{z},\tilde{\rho},\omega)}
        }_{\text{pre-factor}}

        \underbrace{
            \mathrm{e}^{-\frac{\mathbf{i}\omega \tilde{\rho}^2}{2c(z-\tilde{z})}}
        }_{\substack{\text{near-field} \\ \text{factor}}}

        \underbrace{
            [\widehat{\partial_t j}(\tilde{z},\tilde{\rho},\omega)]_{F_{v}}
        }_{\substack{\text{dipole acceleration}\\
                     v\text{-reference frame}}}

        J_0 \left( \frac{\omega \rho \tilde{\rho}}{c(z-\tilde{z})} \right)

        \underbrace{
            \tilde{\rho} \, \mathrm{d}\tilde{\rho}\mathrm{d}\tilde{z}
        }_{\substack{\text{volumetric} \\
                     \text{integration}}} \,,
$$
the phase, $\Phi$, in the pre-factor includes local density $\varrho$, linear complex dispersion (which includes absorption), and the adjustement to the proper reference frame.
$$
        \Phi(\tilde{z},\tilde{\rho},\omega)
    =
        \int_0^{\tilde{z}}
            \frac{\omega}{c}
            (n(\tilde{\tilde{z}},\tilde{\rho},\omega)
            -
            n_{\mathrm{eff}})
        \,\mathrm{d}\tilde{\tilde{z}}
        +
         \mathbf{i} \Phi_{\text{renorm}}
    \,.
$$

There are details about various parts of the integral and the phase:
* **pre-factor**: It handles 1) the longitudinal aspects of the model that are the phase-matching and absorption, 2) the macroscopic density profile as $\widehat{\partial_t j}(\tilde{z},\tilde{\rho},\omega)$ is the single-atom response ignoring the local macroscopic density.
* **near-field factor**: This is an optional pre-factor switching between the Fresnel (factor included) and Fraunhoffer (without the factor) diffraction. Fresnel diffraction is more accurate and allows the first-order description of converging beams. The default set-up includes this pre-factor.
* **the dipole acceleration $[\widehat{\partial_t j}(\tilde{z},\tilde{\rho},\omega)]_{F_{v}}$**: This is the microscopic response from the microscopic model. $F_v$ denotes the co-moving reference frame, with the velocity $v$, of $\widehat{\partial_t j}$. See $n_{\mathrm{eff}}$ for the synchronisation of the frames. (The reason for this is that raw data from CTDSE can be fed directly. And it can be adapted easily for user's convenience of arbitrary generated dipole responses.)
* **The refractive index $n$**: The complex refractive index as a function of $\omega$. It, therefore, includes the absorption. By default, it sources the tables either from `NIST` or `Henke`. So, it is not accessed directly except the table specification.
* **The effective refractive index $n_{\mathrm{eff}}$**: It adjusts the reference frame of $[\widehat{\partial_t j}(\tilde{z},\tilde{\rho},\omega)]_{F_{v}}$ by setting $n_{\mathrm{eff}}=c/v$.
* **The renormalisation factor $\mathbf{i} \Phi_{\text{renorm}}$**: This factor makes the absorption to act as attenuation. The last plane in the medium is then weighted by one. It then makes the generation comparable for varius medium lengths. Finally, this factor is applied continuously during the propagation for a commulative integral, in the case it is required on the output.  
There is the explanation of this factor. The absorption is modelled by the imaginary part of the refractive index $n=n_r-\mathbf{i}n_i$. Concenrning on the absorption, the signal is computed as $S_{\text{ampl.}}=\int_0^L f(z)\mathrm{e}^{F(z)} \,\mathrm{d}z$, where $F(z)=\int_0^z n_i(\tilde{z})\, \mathrm{d}\tilde{z}$. In fact, the model of imaginary refractive index thus treats the absorption as *amplification* along $z$. The model of attenuation is, instead, $S_{\text{att.}}=\int_0^L f(z)\mathrm{e}^{-\int_z^{L}n_i(\tilde{z})\,\mathrm{d}\tilde{z}} \,\mathrm{d}z$. By simple manipulation with the integral: $\int_z^{L}n_i(\tilde{z})\,\mathrm{d}\tilde{z}=\int_0^{L}n_i(\tilde{z})\,\mathrm{d}\tilde{z}-\int_0^{z}n_i(\tilde{z})\,\mathrm{d}\tilde{z} = F(L)-F(z)$. The *attenuated* signal is thus $S_{\text{att.}}=\int_0^L f(z)\mathrm{e}^{F(z)-F(L)} \,\mathrm{d}z = \mathrm{e}^{-F(L)}S_{\text{ampl.}}$ Finally, $\mathrm{e}^{-F(L)}$ is then exactly the "renormalisation" in the model.
* **The phase factor $\Phi$**: This quantity then describes the dephasing of the elementary emitters $\widehat{\partial_t j}$ along $z$. It is important to note that $\widehat{\partial_t j}$ *carries* its own phase (in our case the phase imprinted by the numerical field from CUPRAD + the linear phase given by the $v$-reference frame). So the total phase is the sum of the phase of $\widehat{\partial_t j}$ and $\Phi$.  
Let us clarify the phase issues, which are complex due to the two origins of the phase and the reference frame. Assume there is a homogeneous medium and only linear dispersion given by $n_{\mathrm{IR}}=1+\chi_{\mathrm{IR}}$ (for the driver) and $n_{\mathrm{XUV}}=1+\chi_{\mathrm{XUV}}$ (for the harmonics). Next, the phase given by the reference frame $[\cdot]_{F{v}}$ is $\omega z /v$. In total $$ \Phi_{\text{total}} = \Phi -\frac{\omega n_{\mathrm{IR}}}{c}z + \frac{\omega z}{v} = \int_0^z \left( \frac{\omega n_{\mathrm{XUV}}}{c} -\frac{\omega}{v} \right) \, \mathrm{d}z -\frac{\omega n_{\mathrm{IR}}}{c}z + \frac{\omega z}{v} = \\ = \frac{\omega}{c}(n_{\mathrm{XUV}}-n_{\mathrm{IR}})z = \frac{\omega}{c} (\chi_{\mathrm{XUV}} - \chi_{\mathrm{IR}}) \,. $$ This is the expected outcome where the signal is perfectly phase-matched (all the emitters are in phase and sum constructively) for $\chi_{\mathrm{XUV}} = \chi_{\mathrm{IR}}$.





# Development
The main API is the `HFn2.py`-module. There shall be some good examples of calling it. THe cannonical operation is pipelined after 1D-TDSE, it processes `*.hdf5` dipoles and compute their diffraction integral. Beside this, there is a few other uses of this code (*Maker fringes*, $\omega + 2\omega$ HHG, ...) The preparation shall consider following ideas:

## `Hfn2.py`
* Actual inplementation assumes all the dipoles are in the input array. This might be too memory-demanding. Maybe we can create a *class* providing the dipoles. This class could internally treat storing them or providing them on-the-fly. It would need some thinking about the efficiency.

## Calling the `Hfn2.py`
* The actual example for large-scale applications uses `multiprocessing`, i.e. parallelisation limited to multithreading. THe performace is thus limited by the # of threads per core... We can think about using `mpi4py`, it would require to develop a data treatment as they cannot be easily managed in MPI.

## General remarks
* The output is now in *arbitrary units*. Principally, all physical constants are included and we shall be able to retrieve true XUV intensity.[^1]
* Clear junk files. (This directory contains many junk files for testing etc.)


[^1]: This would require some testing and verifying. So far, we've been using *arbitrary units* for all comparisons with experiments etc. 