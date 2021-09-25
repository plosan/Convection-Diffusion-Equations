# Convection-Diffusion Equations (CDEs)

## The equations

The convection--diffusion equations (CDEs) is a set of conservation laws, constituted by the continuity equation 

<!-- $$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Cpartial%20%5Crho%7D%7B%5Cpartial%20t%7D%20%2B%20%5Cnabla%20%5Ccdot%20(%5Crho%20%5Cmathbf%7Bv%7D)%20%3D%200"></div>

and the general convection-diffusion equation (GCDE)

<!-- $$
\rho \frac{\partial \phi}{\partial t} + \rho \mathbf{v} \cdot \nabla \phi = \nabla \cdot \left( \Gamma_\phi \nabla \phi \right) + \dot{s}_\phi
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Crho%20%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20t%7D%20%2B%20%5Crho%20%5Cmathbf%7Bv%7D%20%5Ccdot%20%5Cnabla%20%5Cphi%20%3D%20%5Cnabla%20%5Ccdot%20%5Cleft(%20%5CGamma_%5Cphi%20%5Cnabla%20%5Cphi%20%5Cright)%20%2B%20%5Cdot%7Bs%7D_%5Cphi"></div>

Here, <!-- $\rho$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Crho"> is the density of the fluid, <!-- $\mathbf{v}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bv%7D"> is the velocity field, <!-- $\phi$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cphi"> is a scalar magnitude of the fluid (such as the temperature or the concentration of a pollutant), <!-- $\Gamma_\phi$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5CGamma_%5Cphi"> is the diffusion coefficient and <!-- $\dot{s}_\phi$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cdot%7Bs%7D_%5Cphi"> is the source term. When <!-- $\phi$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cphi"> is a vector magnitude (for instance the velocity field <!-- $\mathbf{v}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bv%7D">), the equation is 

<!-- $$
\frac{\partial(\rho \phi)}{\partial t} + \nabla \cdot (\rho \mathbf{v} \otimes \mathbf{\phi}) = \nabla \cdot \left( \Gamma_\phi \nabla \phi \right) + \dot{s} _\phi
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Cpartial(%5Crho%20%5Cphi)%7D%7B%5Cpartial%20t%7D%20%2B%20%5Cnabla%20%5Ccdot%20(%5Crho%20%5Cmathbf%7Bv%7D%20%5Cotimes%20%5Cmathbf%7B%5Cphi%7D)%20%3D%20%5Cnabla%20%5Ccdot%20%5Cleft(%20%5CGamma_%5Cphi%20%5Cnabla%20%5Cphi%20%5Cright)%20%2B%20%5Cdot%7Bs%7D%20_%5Cphi"></div> 

where <!-- $\otimes$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cotimes"> denotes the exterior product of two vectors.

Notice that the GCDE for a scalar magnitude <!-- $\phi$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cphi"> is ''combination'' of the linear transport equation 

<!-- $$
\frac{\partial \phi}{\partial t} + \mathbf{v} \cdot \nabla \phi = \dot{s}_\phi
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20t%7D%20%2B%20%5Cmathbf%7Bv%7D%20%5Ccdot%20%5Cnabla%20%5Cphi%20%3D%20%5Cdot%7Bs%7D_%5Cphi"></div>

and the diffusion/heat equation

<!-- $$
\frac{\partial \phi}{\partial t} = \Delta \phi + \dot{s}_\phi
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20t%7D%20%3D%20%5CDelta%20%5Cphi%20%2B%20%5Cdot%7Bs%7D_%5Cphi"></div>

Therefore it is expected that the solution to a boundary-value problem involving the GCDE, when it exists, shares properties with the solution to a linear transport problem and a diffusion problem.

## Project outline

The project focuses on solving the CDEs numerically following a finite volume approach. The problems considered always take place in rectangular domains, hence a cartesian mesh is suitable in order to solve them. The report is structured in the following way:

1. Introduction: brief summary of the project.
2. Convection-diffusion equations: rigorous derivation of the CDEs.
3. Numerical study of the convection-diffusion equations: discretization of the CDEs in a rectangular domain discretized by means of a cartesian mesh. The algorithm to solve a general transient problem is also given.
4. Diagonal flow case: numerical solution to a CDEs steady state problem, taking place in a square domain with a flow in the diagonal diretion.
5. Smith-Hutton case: numerical solution to a CDEs steady state problem, taking place in a rectangular domain with a ''circular'' flow.
6. Appendices: quick reference for some facts in Measure Theory, Ordinary Differential Equations and Numerical resolution of linear systems.

A C++ code was developed in order to solve these problems numerically. 

Hereinafter, <!-- $\rho$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Crho">, <!-- $\Gamma$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5CGamma"> are known constants, <!-- $\dot{s}_\phi = 0$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cdot%7Bs%7D_%5Cphi%20%3D%200"> and <!-- $\mathbf{v}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bv%7D"> does not depend upon time. Under these hypothesis, the CDEs are

<!-- $$
\nabla \cdot \mathbf{v} = 0
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cnabla%20%5Ccdot%20%5Cmathbf%7Bv%7D%20%3D%200"></div>

<!-- $$
\frac{\rho}{\Gamma} \mathbf{v} \cdot \nabla \phi = \Delta \phi
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Crho%7D%7B%5CGamma%7D%20%5Cmathbf%7Bv%7D%20%5Ccdot%20%5Cnabla%20%5Cphi%20%3D%20%5CDelta%20%5Cphi"></div>

## Diagonal flow case

Let <!-- $\phi_\text{low} < \phi_\text{high}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cphi_%5Ctext%7Blow%7D%20%3C%20%5Cphi_%5Ctext%7Bhigh%7D"> be two given constants, and let <!-- $L > 0$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=L%20%3E%200"> be a fixed length giving the square domain <!-- $\Omega = (0,L) \times (0,L)$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5COmega%20%3D%20(0%2CL)%20%5Ctimes%20(0%2CL)">. The velocity field is given by <!-- $\mathbf{v} = \frac{v_0}{\sqrt{2}} \mathbf{i} + \frac{v_0}{\sqrt{2}} \mathbf{j}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bv%7D%20%3D%20%5Cfrac%7Bv_0%7D%7B%5Csqrt%7B2%7D%7D%20%5Cmathbf%7Bi%7D%20%2B%20%5Cfrac%7Bv_0%7D%7B%5Csqrt%7B2%7D%7D%20%5Cmathbf%7Bj%7D">, with <!-- $v_0 > 0$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=v_0%20%3E%200"> a known constant. The the CDE is

<!-- $$
\frac{\rho}{\Gamma} \mathbf{v} \cdot \nabla \phi = \frac{1}{\sqrt{2} L} \frac{\rho v_0 L}{\Gamma} \left( \frac{\partial \phi}{\partial x} + \frac{\partial \phi}{\partial y} \right) = \beta \, \mathrm{Pe} \left( \frac{\partial \phi}{\partial x} + \frac{\partial \phi}{\partial y} \right) = \Delta \phi
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Crho%7D%7B%5CGamma%7D%20%5Cmathbf%7Bv%7D%20%5Ccdot%20%5Cnabla%20%5Cphi%20%3D%20%5Cfrac%7B1%7D%7B%5Csqrt%7B2%7D%20L%7D%20%5Cfrac%7B%5Crho%20v_0%20L%7D%7B%5CGamma%7D%20%5Cleft(%20%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20x%7D%20%2B%20%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20y%7D%20%5Cright)%20%3D%20%5Cbeta%20%5C%2C%20%5Cmathrm%7BPe%7D%20%5Cleft(%20%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20x%7D%20%2B%20%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20y%7D%20%5Cright)%20%3D%20%5CDelta%20%5Cphi"></div>

where <!-- $\mathrm{Pe}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathrm%7BPe%7D"> is Péclet's number. Let <!-- $C_1 = [0,L) \times {0} \cup {L} \times [0,L)$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=C_1%20%3D%20%5B0%2CL)%20%5Ctimes%20%7B0%7D%20%5Ccup%20%7BL%7D%20%5Ctimes%20%5B0%2CL)"> and <!-- $C_2 = {0} \times (0,L] \cup (0,L] \times {L}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=C_2%20%3D%20%7B0%7D%20%5Ctimes%20(0%2CL%5D%20%5Ccup%20(0%2CL%5D%20%5Ctimes%20%7BL%7D"> be two curves in <!-- $\partial \Omega$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cpartial%20%5COmega"> and consider the function

<!-- $$
g(x,y) = 
\left\{
    \begin{aligned}
        &\phi_\text{low} & &\text{if } (x,y) \in C_1 \\
        &\phi_\text{high} & &\text{if } (x,y) \in C_2 \\
        &0 & &\text{otherwise}
    \end{aligned}
\right.
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=g(x%2Cy)%20%3D%20%0A%5Cleft%5C%7B%0A%20%20%20%20%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%26%5Cphi_%5Ctext%7Blow%7D%20%26%20%26%5Ctext%7Bif%20%7D%20(x%2Cy)%20%5Cin%20C_1%20%5C%5C%0A%20%20%20%20%20%20%20%20%26%5Cphi_%5Ctext%7Bhigh%7D%20%26%20%26%5Ctext%7Bif%20%7D%20(x%2Cy)%20%5Cin%20C_2%20%5C%5C%0A%20%20%20%20%20%20%20%20%260%20%26%20%26%5Ctext%7Botherwise%7D%0A%20%20%20%20%5Cend%7Baligned%7D%0A%5Cright."></div>

The diagonal flow case problem is the following boundary-value problem:

<!-- $$
\left\{
    \begin{aligned} 
        \Delta \phi - \left( \frac{\partial \phi}{\partial x} + \frac{\partial \phi}{\partial y} \right) \beta \, \mathrm{Pe} &= 0 & &\text{in } \Omega \\
        \phi &= g & &\text{on } \partial \Omega
    \end{aligned}
\right.
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cleft%5C%7B%0A%20%20%20%20%5Cbegin%7Baligned%7D%20%0A%20%20%20%20%20%20%20%20%5CDelta%20%5Cphi%20-%20%5Cleft(%20%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20x%7D%20%2B%20%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20y%7D%20%5Cright)%20%5Cbeta%20%5C%2C%20%5Cmathrm%7BPe%7D%20%26%3D%200%20%26%20%26%5Ctext%7Bin%20%7D%20%5COmega%20%5C%5C%0A%20%20%20%20%20%20%20%20%5Cphi%20%26%3D%20g%20%26%20%26%5Ctext%7Bon%20%7D%20%5Cpartial%20%5COmega%0A%20%20%20%20%5Cend%7Baligned%7D%0A%5Cright."></div>

### Analytical study

The boundary-value problem is studied for two values of <!-- $\mathrm{Pe}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathrm%7BPe%7D">:

- For <!-- $\mathrm{Pe} = \infty$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathrm%7BPe%7D%20%3D%20%5Cinfty">, the previous boundary-value problem becomes a transport problem. The analytical solution is found via the method of characteristics. However, it is not a solution in the classical sense. The weak solution is not studied.
- For <!-- $\mathrm{Pe} \in [0,\infty)$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathrm%7BPe%7D%20%5Cin%20%5B0%2C%5Cinfty)">, a second-order elliptic PDE is obtained. By means of energy methods in Sobolev spaces (explained in Chapter 6 of Lawrence C. Evan's excellent book ''Partial Differential Equations''), the existence of weak solution is proved, although it cannot be said whether it is unique or not. This solution turns out to be a <!-- $\mathcal{C}^\infty$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathcal%7BC%7D%5E%5Cinfty"> function in <!-- $\Omega$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5COmega">.

### Numerical study

For several values of Péclet's number in the range <!-- $[10^{-9}, 10^9]$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5B10%5E%7B-9%7D%2C%2010%5E9%5D">, the numerical solution is computed using the aforementioned C++ code. 

![Alt text](readme_images/diagonal_case_1.PNG?raw=true "Title")

![Alt text](readme_images/diagonal_case_2.PNG?raw=true "Title")

## Smith-Hutton case

Let <!-- $L > 0$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=L%20%3E%200"> be a fixed length and the domain <!-- $\Omega = (-L,L) \times (0,L)$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5COmega%20%3D%20(-L%2CL)%20%5Ctimes%20(0%2CL)">. Take the velocity field given by

<!-- $$
\mathbf{v} = 2 y (1 - x^2) \mathbf{i} - 2 x (1 - y^2) \mathbf{j}
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bv%7D%20%3D%202%20y%20(1%20-%20x%5E2)%20%5Cmathbf%7Bi%7D%20-%202%20x%20(1%20-%20y%5E2)%20%5Cmathbf%7Bj%7D"></div>

defined in <!-- $\overline{\Omega}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Coverline%7B%5COmega%7D">. Consider the curves $C_1 = [−L, 0] \times \{0\}$, $C_2 = (\{−L\} \times (0,L)) \cup ([−L,L] \times \{L\}) \cup (\{L\} \times [0,L))$ and $C_3 = (0,L) × \{0\}$, which give a partition of <!-- $\partial \Omega$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cpartial%20%5COmega">. Now let <!-- $g \colon C_1 \cup C_2 \to \mathbb{R}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=g%20%5Ccolon%20C_1%20%5Ccup%20C_2%20%5Cto%20%5Cmathbb%7BR%7D"> be the function given by

<!-- $$
g(x,y) = 
\left\{
    \begin{aligned}
        &1 + \tanh{(10(2x+1))} & &\text{if } (x,y) \in C_1 \\
        &1 - \tanh{(10)} & &\text{if } (x,y) \in C_2
    \end{aligned}
\right.
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=g(x%2Cy)%20%3D%20%0A%5Cleft%5C%7B%0A%20%20%20%20%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%261%20%2B%20%5Ctanh%7B(10(2x%2B1))%7D%20%26%20%26%5Ctext%7Bif%20%7D%20(x%2Cy)%20%5Cin%20C_1%20%5C%5C%0A%20%20%20%20%20%20%20%20%261%20-%20%5Ctanh%7B(10)%7D%20%26%20%26%5Ctext%7Bif%20%7D%20(x%2Cy)%20%5Cin%20C_2%0A%20%20%20%20%5Cend%7Baligned%7D%0A%5Cright."></div>

The boundary-value problem for the steady state Smith-Hutton problem is

<!-- $$
\left\{
    \begin{aligned}
        \Delta \phi - \frac{\rho}{\Gamma} \mathbf{v} \cdot \nabla \phi &= 0 & &\text{in } \Omega \\
        \phi &= g & &\text{on } C_1 \cup C_2 \\
        \frac{\partial \phi}{\partial y} &= 0 & &\text{on } C_3
    \end{aligned}
\right.
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cleft%5C%7B%0A%20%20%20%20%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%5CDelta%20%5Cphi%20-%20%5Cfrac%7B%5Crho%7D%7B%5CGamma%7D%20%5Cmathbf%7Bv%7D%20%5Ccdot%20%5Cnabla%20%5Cphi%20%26%3D%200%20%26%20%26%5Ctext%7Bin%20%7D%20%5COmega%20%5C%5C%0A%20%20%20%20%20%20%20%20%5Cphi%20%26%3D%20g%20%26%20%26%5Ctext%7Bon%20%7D%20C_1%20%5Ccup%20C_2%20%5C%5C%0A%20%20%20%20%20%20%20%20%5Cfrac%7B%5Cpartial%20%5Cphi%7D%7B%5Cpartial%20y%7D%20%26%3D%200%20%26%20%26%5Ctext%7Bon%20%7D%20C_3%0A%20%20%20%20%5Cend%7Baligned%7D%0A%5Cright."></div>

Hereinafter, we shall assume <!-- $L = 1 \ \mathrm{m}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=L%20%3D%201%20%5C%20%5Cmathrm%7Bm%7D">.

### Velocity field

Smith-Hutton's problem velocity field studied, in particular, the existence and uniqueness of streamlines is proved. In addition, these lines are computed numerically for some initial conditions, yielding the following plot:

![Alt text](readme_images/smith_hutton_velocity_field.PNG?raw=true "Title")

### Analytical solution

Since the value <!-- $\rho / \Gamma$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Crho%20%2F%20%5CGamma"> appears in the PDE, it seems reasonable to assume that the solution will depend on this constant:

- For <!-- $\rho / \Gamma = \infty$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Crho%20%2F%20%5CGamma%20%3D%20%5Cinfty">, the boundary-value problem turns into a transport-like problem. The characteristics of this problem are the streamlines found in the previous section. Due to the nature of these curves, the analytical solution is harder to find. The expression for this solution is given in terms of another function, although it is not proved whether it is unique. Nonetheless, it is demonstrated that if two classical solutions to the transport-like problem exist, then they must be equal. 
- For <!-- $\rho / \Gamma \in [0,\infty)$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Crho%20%2F%20%5CGamma%20%5Cin%20%5B0%2C%5Cinfty)">, a second-order elliptic PDE is obtained. The theory studied for the diagonal-flow case cannot be applied here (at least in a straightforward fashion), since on <!-- $C_3$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=C_3"> a Neumann boundary condition is prescribed (instead of a Dirichlet boundary condition).

### Numerical solution

For several values in <!-- $\rho / \Gamma$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Crho%20%2F%20%5CGamma"> in the range <!-- $[10^{-9}, 10^9]$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5B10%5E%7B-9%7D%2C%2010%5E9%5D">, the numerical solution to the Smith-Hutton case is computed using the C++ code developed. Below some examples are shown:

![Alt text](readme_images/smith_huttton_+1.PNG?raw=true "Title")

![Alt text](readme_images/smith_huttton_+9.PNG?raw=true "Title")

![Alt text](readme_images/smith_huttton_-9.PNG?raw=true "Title")







