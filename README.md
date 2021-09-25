# Convection-Diffusion Equations

The convection--diffusion equations is a set of conservation laws, constituted by the continuity equation <!-- $\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0$ --> 

<img style="transform: translateY(0.1em); background: white;" src="svg\EjZrnvyEsQ.svg"> 

and the general convection-diffusion equation <!-- $\rho \frac{\partial \phi}{\partial t} + \rho \mathbf{v} \cdot \nabla \phi = \nabla \cdot \left( \Gamma_\phi \nabla \phi \right) + \dot{s}_\phi$ --> 

<img style="transform: translateY(0.1em); background: white;" src="svg\96T6pB215v.svg">

Here $\rho$ is the density of the fluid, $\mathbf{v}$ is the velocity field, $\phi$ is a scalar magnitude of the fluid (such as the temperature or the concentration of a pollutant ), $\Gamma_\phi$ is the diffusion coefficient and $\dot{s}_\phi$ is the source term. When $\phi$ is a vector magnitude (for instance the velocity field $\mathbf{v}$), the equation is $$\frac{\partial(\rho \phi)}{\partial t} + \nabla \cdot (\rho \mathbf{v} \otimes \mathbf{\phi}) = \nabla \cdot \left( \Gamma_\phi \nabla \phi \right) + \dot{s}_\phi$$ where $\otimes$ denotes the exterior product of two vectors.








