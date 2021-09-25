# Convection-Diffusion Equations

The convection--diffusion equations is a set of conservation laws, constituted by the continuity equation 

<!-- $$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Cpartial%20%5Crho%7D%7B%5Cpartial%20t%7D%20%2B%20%5Cnabla%20%5Ccdot%20(%5Crho%20%5Cmathbf%7Bv%7D)%20%3D%200"></div>

and the general convection-diffusion equation 

<!-- $$
\rho \frac{\partial \phi}{\partial t} + \rho \mathbf{v} \cdot \nabla \phi = \nabla \cdot \left( \Gamma_\phi \nabla \phi \right) + \dot{s}_\phi
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Crho%20%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20t%7D%20%2B%20%5Crho%20%5Cmathbf%7Bv%7D%20%5Ccdot%20%5Cnabla%20%5Cphi%20%3D%20%5Cnabla%20%5Ccdot%20%5Cleft(%20%5CGamma_%5Cphi%20%5Cnabla%20%5Cphi%20%5Cright)%20%2B%20%5Cdot%7Bs%7D_%5Cphi"></div>

Here, <!-- $\rho$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Crho"> is the density of the fluid, <!-- $\mathbf{v}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bv%7D"> is the velocity field, <!-- $\phi$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cphi"> is a scalar magnitude of the fluid (such as the temperature or the concentration of a pollutant), $\Gamma_\phi$ is the diffusion coefficient and <!-- $\dot{s}_\phi$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cdot%7Bs%7D_%5Cphi"> is the source term. When <!-- $\phi$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cphi"> is a vector magnitude (for instance the velocity field <!-- $\mathbf{v}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bv%7D">), the equation is 

<!-- $$
\frac{\partial(\rho \phi)}{\partial t} + \nabla \cdot (\rho \mathbf{v} \otimes \mathbf{\phi}) = \nabla \cdot \left( \Gamma_\phi \nabla \phi \right) + \dot{s} _\phi
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Cpartial(%5Crho%20%5Cphi)%7D%7B%5Cpartial%20t%7D%20%2B%20%5Cnabla%20%5Ccdot%20(%5Crho%20%5Cmathbf%7Bv%7D%20%5Cotimes%20%5Cmathbf%7B%5Cphi%7D)%20%3D%20%5Cnabla%20%5Ccdot%20%5Cleft(%20%5CGamma_%5Cphi%20%5Cnabla%20%5Cphi%20%5Cright)%20%2B%20%5Cdot%7Bs%7D%20_%5Cphi"></div> 

where <!-- $\otimes$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cotimes"> denotes the exterior product of two vectors.








