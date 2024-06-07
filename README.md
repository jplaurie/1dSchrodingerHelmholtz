# 1dOpticalWaveTurbulence
Numerical code for types of 1D optical wave turbulence equation for the 1D complex wave function $`\psi = \psi(x,t)\in \mathbb{C}`$.

## Non-local OWT Equation
```math
i\frac{\partial \psi}{\partial t}=-\frac{1}{2}\frac{\partial^2 \psi}{\partial x^2} -\frac{1}{2}\psi\left( 1- \frac{1}{g}\frac{\partial^2}{\partial x ^2}\right)^{-1}|\psi|^2 -i\alpha \left(-\frac{\partial^2}{\partial x^2}\right)^{\alpha_{power}}\psi -i\nu \left(-\frac{\partial^2}{\partial x^2}\right)^{\nu_{power}}\psi + iF
```

## Long-Wave Limit OWT Equation
```math
i\frac{\partial \psi}{\partial t}=-\frac{1}{2}\frac{\partial^2 \psi}{\partial x^2}  -\frac{1}{2}\psi|\psi|^2 - \frac{1}{2g}\psi\frac{\partial^2|\psi|^2}{\partial x ^2}  -i\alpha \left(-\frac{\partial^2}{\partial x^2}\right)^{\alpha_{power}}\psi -i\nu \left(-\frac{\partial^2}{\partial x^2}\right)^{\nu_{power}}\psi + iF
```

## Parameters

The main dynamical parameters are the chemical potential $`\mu`$ and level of nonlinearity $`g`$. Parameters $`\alpha`$ and $`\nu`$ correspond to the hypo- an hyper-viscosity coefficients with the order of the two dissipations given by $`\alpha_{power}`$ and $`\nu_{power}`$.

The structure of the additive forcing is such that
```math
F = \sum_{k} f_k \eta_k(t),
```
where $`\eta_k(t)`$ is a complex Gaussian white noise at mode $`k`$ and $`f_k\in \mathbb{R}`$ is the distribution of the forcing amplitudes in Fourier space.

## Boundary Conditions
The equation is solved with periodic boundary conditions with a domain size $`L`$. By default we set $`L=2\pi`$ so that the wavenumebers take integer values.

## Spatial Discretization
Pseudo-spectral method with full dealiasing with the 3/2-rule. We use the FFTW library that stores the wavenumbers in the for $`\mathbf{k} = 0,1,2,\dots, N/2, -N/2 +1, \dots, -2,-1`$, where $`N`$ is the spatial resolution.

## Temporal Discretization
4th-order exponential time differencing Runge-Kutta method (see Kassam 2005, and Cox 2002 for details.).




# How to Run the Code

Please see the `instructions.pdf` file in the folder `/docs/`

