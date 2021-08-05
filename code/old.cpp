// // Physical data
// Sizes
double x0 = 0;  // Lower left corner x coordinate for rectangular domain    [m]
double y0 = 0;  // Lower left corner y coordinate for rectangular domain    [m]
double lx = 1;  // Domain size in x axis                                    [m]
double ly = 1;  // Domain size in y axis                                    [m]
double lz = 1;  // Domain size in z axis                                    [m]
// Boundary conditions
const double phi_low = 10;      // Minimum value for phi
const double phi_high = 20;     // Maximum value for phi
// Flow
const double v0 = 1;                // Velocity modulus     [m/s]
const double alpha = 0.25 * M_PI;   // Velocity angle       [rad]
// Thermophysical properties for water at 20 ºC
const double lambda = 0.5861;       // Thermal conductivity                         [W/(k*m)]
const double cv = 4183;             // Specific heat at constant volume (pressure)  [J/(kg*K)]
const double rho = 998.2;           // Density                                      [kg/m^3]
const double gamma = lambda / cv;   // Diffusion coefficient

// // Numerical data
unsigned int nx = 5;      // Number of nodes in x axis
unsigned int ny = 5;      // Number of nodes in y axis
const double phi0 = 1;      // Initial value to fill phi vector for linear system resolution
const double tol = 1e-15;   // Tolerance to stop iteration

double* nodeX = (double*) malloc(nx * sizeof(double*));
double* nodeY = (double*) malloc(ny * sizeof(double*));

double* distX = (double*) malloc((nx - 1) * sizeof(double*));
double* distY = (double*) malloc((ny - 1) * sizeof(double*));

double* faceX = (double*) malloc((nx + 1) * sizeof(double*));
double* faceY = (double*) malloc((ny + 1) * sizeof(double*));

double* surfX = (double*) malloc(ny * sizeof(double*));
double* surfY = (double*) malloc(nx * sizeof(double*));

double* vol = (double*) malloc(nx * ny * sizeof(double*));

compute2DUniformRectangularMesh(x0, y0, nx, ny, lx, ly, lz, nodeX, nodeY, distX, distY, faceX, faceY, surfX, surfY, vol);
// printMeshInfo(x0, y0, nx, ny, lx, ly, lz, nodeX, nodeY, distX, distY, faceX, faceY, surfX, surfY, vol);





double* A = (double*) malloc(5 * nx * ny * sizeof(double*));
double* b = (double*) malloc(nx * ny * sizeof(double*));
double* phi_boundary = (double*) malloc(2 * sizeof(double*));
phi_boundary[0] = phi_low;
phi_boundary[1] = phi_high;

double* v = (double*) malloc(2 * sizeof(double*));
v[0] = v0 * cos(alpha);
v[1] = v0 * sin(alpha);

int scheme = 0;
computeDiscretizationCoefficientsDiagonalCase(nx, ny, nodeX, nodeY, distX, distY, faceX, faceY, surfX, surfY, vol, phi_boundary, v, rho, gamma, A, b, scheme);

checkSystemMatrix(nx, ny, tol, A);
double* phi = (double*) malloc(nx * ny * sizeof(double*));
std::fill_n(phi, nx*ny, phi0);

const int maxIt = 1e6;
solveSystem(nx, ny, tol, maxIt, A, b, phi);

// verification(lx, ly, lz, nx, ny, nodeX, nodeY, distX, distY, faceX, faceY, surfX, surfY, vol, phi_boundary, v, rho, gamma, A, b, phi, scheme);

// Free memory allocated
free(nodeX);
free(nodeY);
free(distX);
free(distY);
free(faceX);
free(faceY);
free(surfX);
free(surfY);
free(vol);

free(A);
free(b);
free(phi_boundary);
free(v);
free(phi);

return 0;
