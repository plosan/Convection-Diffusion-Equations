# Convection--Diffusion Equations

The convection--diffusion equations are the following set of conservation laws
<!-- $$
    \frac{\partial \rho}{\partial t} + \nabla \vdot (\rho \mathbf{v}) = 0
$$ --> 

<div align="center"><img style="background: white;" src="svg\IuGE4WMyWw.svg"></div>
<!-- $$
    \rho \frac{\partial \phi}{\partial t} + \rho \mathbf{v} \vdot \nabla \phi = \nabla \vdot \left( \Gamma_\phi \nabla \phi \right) + \dot{s}_\phi
$$ --> 

<div align="center"><img style="background: white;" src="svg\ZMvskRUSsx.svg"></div>
where <!-- $\rho$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\H1XCgAMfat.svg"> is the density of the fluid, $\mathbf{v}$ is the velocity field, $\phi$ is a scalar magnitude of the fluid (such as the temperature $T$ or the concentration of a pollutant $Y$), $\Gamma_\phi$ is the diffusion coefficient and $\dot{s}_\phi$ is the source term. When $\phi$ is a vector magnitude (for instance the velocity field $\mathbf{v}$), the equation is
$$
    \frac{\partial(\rho \phi)}{\partial t} + \nabla \vdot (\rho \mathbf{v} \otimes \mathbf{\phi}) = \nabla \vdot \left( \Gamma_\phi \nabla \phi \right) + \dot{s}_\phi
$$
where $\otimes$ is the exterior product of two vectors.








