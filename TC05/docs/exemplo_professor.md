# Modelagem e Simulação Multifísica - Aula 18

## Trabalho Computacional 05 - O Levitador Magnético Linear

No Trabalho Computacional 01, apresentamos as equações diferenciais ordinárias que descrevem a dinâmica de um sistema de levitação magnética:

\begin{eqnarray}
u(t) &=& RI(t) + L(x)\frac{dI(t)}{dt} + Bl\frac{dx(t)}{dt} \\
m \ddot x &=& mg - f_m
\end{eqnarray}

onde $x$ é a distância entre o objeto em suspensão e o eletroímã, $I$ é a corrente na bobina, $g$ é a aceleração da gravidade, $m$ é a massa do objeto e $f_m$ é a força magnética do eletroímã.

Neste trabalho, estaremos interessados em encontrar valores dos parâmetros usados no modelo dinâmico a partir da solução da equação diferencial parcial que governa o problema eletromagnético do levitador. 

## 1 - Conceitos básicos

Considere o eletroímã apresentado a seguir:  

Fig. 1 Eletroímã planar
<p>
    <img width="450" height="100" src="Figs/Aula16_eletroima.png">
</p>

Assumindo que a corrente no eletroímã é constante, o problema pode ser modelado como um problema magnetostático que obedece às seguintes equações diferenciais parciais em cada ponto do domínio:

\begin{eqnarray}
\nabla\cdot \textbf{B} &=& 0, \\%\label{eq_Gauss}\\
\nabla\times \textbf{H} &=& \textbf{J},%\label{eq_Ampere}
\end{eqnarray}

onde $\textbf{H}$ é Campo Magnético $(A/m)$, $\textbf{B}$ é a Densidade de Fluxo Magnético dada em Tesla $(T)$ e $\textbf{J}$ é a densidade de corrente elétrica $(A/m^2)$. Os campos $\textbf{B}$ e $\textbf{H}$ se relacionam com as características dos materiais através da equação constitutiva
 
\begin{equation}
\textbf{B} = \mu\textbf{H} %\label{eq_constitutiva}
\end{equation}

onde $\mu$ é a permeabilidade magnética do material.

O eletroímã, é composto por três materiais distintos sendo eles i) o cobre que compõe a bobina, ii) o material ferromagnético que compõe o eletroímã e a peça em levitação e iii) o ar que completa o domínio do problema.

Tanto o ar quanto o cobre não possuem características magnéticas relevantes e suas permeabilidades magnéticas podem ser consideradas iguais a do espaço livre ($\mu_{ar}=\mu_{Cu}=\mu_0=4\pi\cdot 10^{-7}H/m$). Já materiais ferromagnéticos, apresentam uma permeabilidade elevada quando comparada com a do espaço livre. Tipicamente, materiais ferromagnéticos tem características não lineares e a permeabilidade é geralmente obtida da curva de magnetização do material (Curva $\textbf{B}\times\textbf{H}$). Entretanto, nesse trabalho iremos considerar que o material ferromagnético é linear, podendo ser aproximado por $\mu_{fer} = \mu_r\mu_0$, sendo $\mu_r$ uma constante positiva. 

### 1.1 Aproximação analítica

<p>
    Uma <em> <a href="https://www.youtube.com/watch?v=A51jWH5Ar4s"> solução analítica </a> para esse problema pode ser obtida assumindo-se as seguintes considerações: </em>
</p>

i) O levitador irá operar fora da região de saturação da curva $\textbf{B}\times\textbf{H}$, logo o valor de $\mu$ no material ferromagnético será elevado e o campo magnético no interior desse material poderá ser desprezado.

ii) Como a distância entre a peça e o eletroímã é pequena, o espraiamento das linhas de fluxo magnético poderá ser desprezado e o fluxo magnético no entreferro (região do ar entre as pontas do eletroímã e a peça) será constante e normal à superfície $S$.

Nesse caso, aplicando a lei de Ampère em um caminho fechado que percorre os entreferros e os materiais ferromagnéticos temos:

\begin{eqnarray}
\oint_{L'(S')} \textbf{H}\cdot d\textbf{l}' &=& \int_{S'} \textbf{J}\cdot d\textbf{S}'.\nonumber
\end{eqnarray}

Aplicando as considerações anteriores, a equação acima se reduz a
\begin{eqnarray}
2Hx &=& NI,\nonumber
\end{eqnarray}

uma vez que o campo $\textbf{H}$ é diferente de zero apenas nos dois entreferros que possuem comprimento $x$.

Reescrevendo a equação, é possível obter o módulo do campo magnético 
\begin{eqnarray}
H &=& \frac{NI}{2x}\nonumber
\end{eqnarray}
onde $N$ é o número de espiras da bobina.


### 1.2 Cálculo da Indutância
<p>
<em>
Uma vez conhecido o campo magnético, podemos calcular a  <a href="https://www.youtube.com/watch?v=pq3QZRXfsjA"> indutância própria </a> do eletroímã. </em>
</p>

Se lembrarmos que a indutância própria é a capacidade de uma bobina de $N$ espiras em criar o fluxo magnético $\phi$ com uma determinada corrente $I$, podemos encontra a indutância da bobina por meio da seguinte expressão:
\begin{equation}
L = \frac{\lambda}{I} = \frac{N\phi}{I},
\end{equation}
onde o fluxo concatenado $\lambda = N\phi$ é o fluxo total enlaçado pelas $N$ espiras da bobina e o fluxo $\phi$ representa a integral da densidade de fluxo magnético que atravessa a seção transversal das bobinas. Matematicamente:
\begin{equation}
\phi = \int \mathbf{B}\cdot d\mathbf{s}.
\end{equation}

Desprezando-se o espraiamento no entreferro, o fluxo que atravessa as espiras será igual ao fluxo que que atravessa qualquer uma das pernas do entreferro. Como no entreferro o campo magnético é constante, a integral de superfície será igual ao produto da densidade de fluxo pela seção transversal $S$ do entreferro. 
\begin{equation}
\phi = BS = \mu_0 HS.
\end{equation}    
Dessa forma, a indutância se reduz a:
\begin{equation}
L = \frac{\mu_0 N^2S}{2x}
\end{equation}

### 1.3 Cálculo da força magnética

A força magnética pode ser obtida por meio do  <a href="https://www.youtube.com/watch?v=CpTi77y6yL8"> Trabalho Virtual </a> ou por meio do <a href="https://www.youtube.com/watch?v=dX4rY2P9_BI8"> Tensor de Maxwell </a>. Ambas as abordagens resultam em resultados similares, sendo apresentado aqui apenas os resultados para o Tensor de Maxwell.
    
Uma vez conhecido o campo magnético exterior à superfície da peça em levitação, a densidade de força por unidade de área exercida pelo campo pode ser obtida por meio do tensor de Maxwell:

\begin{eqnarray}
\frac{d\textbf{F}}{ds} &=& -\frac{\mu_0}{2}H^2\textbf{n}+\mu_0(\textbf{H}\cdot\textbf{n})\textbf{H} %\label{eq_tensorMaxwell}
\end{eqnarray}

Integrando a equação anterior sobre a superfície da peça, lembrando que o campo $\textbf{H}$ é normal à peça e diferente de zero apenas nos entreferros, obtemos a força exercida pelo eletroímã

\begin{eqnarray}
F &=& \frac{\mu_0N^2S}{4} \left(\frac{I}{x} \right)^2 \nonumber
\end{eqnarray}

### 1.4 Cálculo do fator de força ($Bl$) e $k$

O Fator de força é definido como a razão entre a força magnética e a corrente 
\begin{equation}
Bl = \frac{fm}{i} = \frac{\mu_0N^2S}{4} \left(\frac{I}{x^2} \right)
\end{equation}

Comparando as equações de indutância e força com as utilizadas na aula 03:

\begin{eqnarray}
L(x) &=& \frac{k}{x} \nonumber\\
f_m(i,x) &=& \frac{k}{2}\left(\frac{I}{x^2} \right) \nonumber
\end{eqnarray}

é possível encontrar uma forma fechada para $k$
\begin{equation}
k = \frac{\mu_0N^2S}{2}.
\end{equation}

### 1.5 Cálculo da resistência

Por fim, o último parâmetro a ser calculado é a resistência da bobina que pode ser obtida diretamente da condutividade do cobre $\sigma = 58.0 (MS/m)$, da seção transversal do fio $S_f$ e do comprimento total do fio $L_f$
\begin{equation}
R = \frac{L_f}{\sigma S_f}.
\end{equation}

Embora o modelo apresentado proporcione a obtenção dos parâmetros de forma analítica, a geometria utilizada difere consideravelmente da utilizada no TC05, sendo de pouca aplicação prática. Na seção seguinte, um modelo mais realista é apresentado. Entretanto, ele não possui solução analítica e será resolvido utilizando o método de elementos finitos.

## 2 - Um levitador axissimétrico

O exemplo apresenta uma bobina contendo N espiras e um núcleo ferromagnético linear de permeabilidade conhecida. Abaixo da bobina está localizada uma peça ferromagnética cilíndrica. 
Todo o problema é parametrizável de forma que diferentes configurações podem ser testadas. 
No caso do enrolamento da bobina, utiliza-se como referência a padronização <a href="https://en.wikipedia.org/wiki/American_wire_gauge"> AWG </a>.

O problema é resolvido usando simetria axissimétrica. Dessa forma, apenas a seção transversal do sólido de revolução precisa ser representada e o problema pode ser resolvido em um domínio bidimensional. 

Fig. 2 Eletroímã axissimétrico
<p>
    <img width="450" height="100" src="Figs/Aula16_levitadorGeoAxiss.png">
</p>

Esse exemplo pode ser validado como o software FEMM em: https://www.femm.info/wiki/MagneticsTutorial.

### 2.1 Plataforma FEniCS

O problema será resolvido com o auxílio da plataforma de computação científica <a href="https://fenicsproject.org/"> FEniCS </a>. As bibliotecas utilizadas são apresentadas a seguir


```python
import numpy as np # For vector and matrix manipulation
import gmsh # Interface with the software GMSH
from scipy import integrate # to post-process the results

# Visualization 
import pyvista
from mpi4py import MPI
import matplotlib.pyplot as plt

# FEniCS
from ufl import (SpatialCoordinate, TestFunction, TrialFunction, as_vector,
                 dx, ds, dS, dot, grad, inner)
from petsc4py import PETSc

from dolfinx import default_scalar_type
from dolfinx.mesh import create_mesh, compute_midpoints, locate_entities_boundary
from dolfinx.io import gmshio, XDMFFile, VTKFile
from dolfinx.io.gmshio import model_to_mesh
from dolfinx.fem import (Constant, dirichletbc, Function, Expression, functionspace, assemble_scalar,
                         form, locate_dofs_geometrical, locate_dofs_topological)
from dolfinx.fem.petsc import LinearProblem
from dolfinx.plot import vtk_mesh
from dolfinx import fem
```

### 2.2 Parâmetros de entrada

Valores comentados como default são os valores que foram usados para a comparação com o software FEMM.


```python
# --- mesh parameters ---
tol = 0 # mesh tolerance 

# --- Source parameter
I = 1.0 #(A) current (default 1.0)

# --- geometry parameters ---

# Core parameters
r_core = 1.0e-2 # (m) (default 1.0e-2)
z_core = 5.0e-2 # (m) (default 5.0e-2)

mu_core = 1200 # relative permeability (default 1200)

# Wire
class AWG:
    def __init__(self, type):
        self.AWG = type
        self.diameter = 0.127e-3*92**((36-type)/39) #(m)
        self.area = np.pi/4*self.diameter**2
        self.sigma = 58.0e6 # (S/m)
awg = AWG(18) # (default 18)

# Coil
N = 100 # number of turns (default 100)
Coil_thickness = N*awg.area/z_core # (m)
r_coil = r_core + Coil_thickness # (m)
print("Coil_thickness = ", Coil_thickness, "(m)")

# Cylinder
x_gap = .5e-2 # distance from the electromagnet (m) (default 0.5e-2)
r_cyl = r_core # (m) (default r_core)
z_cyl = r_core # (m) (default r_core)

mu_cyl = mu_core # relative permeability  (default mu_core)

# air box
r_air = 5*r_coil  # (default 5*r_coil)
z_air = 2.5*z_core # (default 2.5*z_core)
```

    Coil_thickness =  0.0016460936674626287 (m)


### 2.3 Construindo a malha para o Método de Elementos Finitos


```python
rank = MPI.COMM_WORLD.rank

gmsh.initialize()

gdim = 2  # Geometric dimension of the mesh
model_rank = 0
mesh_comm = MPI.COMM_WORLD

if mesh_comm.rank == model_rank:

    # Tags
    core_tag = 1 #inside core
    coil_tag = 2 #inside coil
    cyl_tag = 3  # Cylinder tag
    air_tag = 4
    
    # Creates the air rectangle
    air_rectangle = gmsh.model.occ.addRectangle(0, -z_air/2, 0, r_air, z_air, tag = air_tag)
    gmsh.model.occ.synchronize()

    holes = [] # List of all geometries inside the air_rectangle
    
    # Creates the magnetic core
    core_rectangle = gmsh.model.occ.addRectangle(0, 0, 0, r_core, z_core, tag = core_tag)
    gmsh.model.occ.synchronize()
    holes.append((gdim, core_tag))

    # Creates the coil using AWG 18 wire
    coil_rectangle = gmsh.model.occ.addRectangle(r_core, 0, 0, Coil_thickness, z_core, tag = coil_tag)
    gmsh.model.occ.synchronize()
    holes.append((gdim, coil_tag))
    
    # Magnetic cylinder in levitation
    magnetic_cylinder = gmsh.model.occ.addRectangle(0, -x_gap, 0, r_cyl, -z_cyl, tag = cyl_tag)
    gmsh.model.occ.synchronize()
    holes.append((gdim, cyl_tag))
    
    whole_domain = gmsh.model.occ.fragment([(2, air_rectangle)], holes)
    gmsh.model.occ.synchronize()

    # Creating physical groups
    gmsh.model.addPhysicalGroup(gdim, [core_tag], tag = 1, name = "core_tag") 
    gmsh.model.addPhysicalGroup(gdim, [coil_tag], tag = 2, name = "coil_tag") 
    gmsh.model.addPhysicalGroup(gdim, [cyl_tag], tag = 3, name = "cyl_tag") 
    gmsh.model.addPhysicalGroup(gdim, [air_tag], tag = 4, name = "air_tag") 
    gmsh.model.occ.synchronize()
    
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", r_core/8)
    gmsh.model.mesh.generate(gdim)  
    #gmsh.write("Aula18b.msh")

# Converts the mesh created in GMSH to the format used by Dolfinx.
mesh,cell_tags,facet_tags = gmshio.model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=gdim)
# Close GMSH
gmsh.finalize()

```

    Info    : Meshing 1D...                                                                                                  
    Info    : [  0%] Meshing curve 1 (Line)
    Info    : [ 10%] Meshing curve 2 (Line)
    Info    : [ 20%] Meshing curve 3 (Line)
    Info    : [ 20%] Meshing curve 4 (Line)
    Info    : [ 30%] Meshing curve 5 (Line)
    Info    : [ 30%] Meshing curve 6 (Line)
    Info    : [ 40%] Meshing curve 7 (Line)
    Info    : [ 50%] Meshing curve 8 (Line)
    Info    : [ 50%] Meshing curve 9 (Line)
    Info    : [ 60%] Meshing curve 10 (Line)
    Info    : [ 60%] Meshing curve 11 (Line)
    Info    : [ 70%] Meshing curve 12 (Line)
    Info    : [ 80%] Meshing curve 13 (Line)
    Info    : [ 80%] Meshing curve 14 (Line)
    Info    : [ 90%] Meshing curve 15 (Line)
    Info    : [ 90%] Meshing curve 16 (Line)
    Info    : [100%] Meshing curve 17 (Line)
    Info    : Done meshing 1D (Wall 0.00288746s, CPU 0.003786s)
    Info    : Meshing 2D...
    Info    : [  0%] Meshing surface 1 (Plane, Frontal-Delaunay)
    Info    : [ 30%] Meshing surface 2 (Plane, Frontal-Delaunay)
    Info    : [ 60%] Meshing surface 3 (Plane, Frontal-Delaunay)
    Info    : [ 80%] Meshing surface 4 (Plane, Frontal-Delaunay)
    Info    : Done meshing 2D (Wall 0.161981s, CPU 0.148265s)
    Info    : 5794 nodes 11724 elements


Mostrando a geometria e os PhysicalGroups


```python
# Sets up the visualization with PyVista
pyvista.start_xvfb()

plotter = pyvista.Plotter()
mesh.topology.create_connectivity(gdim, gdim)
grid = pyvista.UnstructuredGrid(*vtk_mesh(mesh, gdim))
num_local_cells = mesh.topology.index_map(gdim).size_local
grid.cell_data["Marker"] = cell_tags.values[cell_tags.indices < num_local_cells]
grid.set_active_scalars("Marker")

# Applies the "viridis" colormap to color the mesh based on "Cell Markers"
actor = plotter.add_mesh(grid, show_edges=True, cmap="viridis", edge_color="black")
actor = plotter.add_title(
    'Fig. 3 Geometria do problema axissimétrico', font='courier', color='k', font_size=7)
plotter.view_xy()
plotter.show();
print('Fig. 3 Geometria do problema axissimétrico')
```

    error: XDG_RUNTIME_DIR is invalid or not set in the environment.
    MESA: error: ZINK: failed to choose pdev
    glx: failed to create drisw screen



    
![png](exemplo_professor_files/exemplo_professor_11_1.png)
    


    Fig. 3 Geometria do problema axissimétrico


## 3 - Representação do problema como um problema de valor de contorno 

### 3.1 Forma forte do problema magnetostático

Em problemas que a densidade de corrente não é nula, definimos um potencial vetor $\textbf{A}$ tal que:

\begin{equation}
\textbf{B}= \nabla\times\textbf{A}\nonumber
\end{equation}

Essa relação é coerente com a equação de Gauss para o campo magnético $\nabla\cdot\textbf{B}=0$ pois o
$\nabla\cdot(\nabla\times\textbf{A})$ é nulo.

Podemos então aplicar essa relação na lei de Ampère da magnetostática:

\begin{eqnarray}
\nabla\times\textbf{H} &=& \textbf{J} \nonumber\\ 
\nabla\times(\frac{1}{\mu}\nabla\times\textbf{A}) &=& \textbf{J} \nonumber
\end{eqnarray}

Além das equações de Maxwell, é preciso garantir as condições de interface na passagem entre meios materiais

\begin{eqnarray}
\hat{n}\cdot(\textbf{B}_2-\textbf{B}_1) &=& 0;\; em\,\Gamma_{12} \nonumber\\
\hat{n}\times(\textbf{H}_2-\textbf{H}_1) &=& \textbf{k};\; em\,\Gamma_{12} \nonumber
\end{eqnarray}

onde $\Gamma_{12}$ representa a fronteira entre dois materiais, como por exemplo, a fronteira entre o núcleo e a região da bobina, ou a fronteira entre a bobina e o ar. Para o nosso exemplo, a densidade linear de corrente $\textbf{k}$ é nula, existindo apenas a densidade superficial de corrente $\textbf{J}$.

Reescrevendo as condições de interface em função do potencial vetor magnético e adicionando condições de contorno para garantir unicidade da solução, encontramos a forma forte do problema: 

Dado $\mu$ e $\textbf{J}$ encontre $A:\bar{\Omega}\rightarrow \mathbb{R}$ tal que:

\begin{eqnarray}
\nabla\times(\frac{1}{\mu}\nabla\times\textbf{A}) = \textbf{J}; \; \forall x \in \Omega \nonumber\\
\textbf{A} = \textbf{A}_q; \; \forall x \in \Gamma_q \nonumber\\
\hat{n}\times\frac{1}{\mu}\nabla\times\textbf{A} = \textbf{H}_t;\; \forall x \in \Gamma_h \nonumber\\
\hat{n}\times(\frac{1}{\mu_2}\nabla\times\textbf{A}_2 - \frac{1}{\mu_1}\nabla\times\textbf{A}_1) = \textbf{k};\; \forall x \in \Gamma_{12} \nonumber
\end{eqnarray}

onde $\Gamma_{q}$ é a região da fronteira onde é imposta uma condição de Dirichlet e $\Gamma_{h}$ representa uma região onde será imposta uma condição de Neumann. 

### 3.2 Condição de contorno

Para o levitador axissimétrico, apenas a condição de Dirichlet Homogênea ($A=0$) é necessária, sendo imposta nas quatro bordas do domínio. Essa condição impõe que a Densidade de Fluxo Magnético deva ser tangencial à fronteira. No FEniCS essa condição é imposta conforme mostrado a seguir:


```python
V = functionspace(mesh, ("Lagrange", 2))
tdim = mesh.topology.dim

# Define the radial coordinate r from x[0]
x = SpatialCoordinate(mesh)
r = x[0]

# Find the boundary nodes and define that all of them belong to the Dirichlet boundary (A=0).
facets = locate_entities_boundary(mesh, tdim - 1, lambda x: np.full(x.shape[1], True))
dofs = locate_dofs_topological(V, tdim - 1, facets)
bc = dirichletbc(default_scalar_type(0), dofs, V)
```

### 3.3 Permeabilidade Magnética ($\mu)$

No caso da permeabilidade dos materiais, tanto a permeabilidade magnética do cobre quanto a do ar são muito próximas da permeabilidade magnética do vácuo, $\mu_0 = 4\pi\cdot 10^{-7}(H/m)$. Para as partes magnéticas, consideraremos os materiais lineares e isotrópicos com $\mu_r = 1200$. Assim como no problema da bobina, usaremos os elementos descontínuos de ordem zero para definir o valor de $\mu$ em cada elemento do domínio.  


```python
Q = functionspace(mesh, ("DG", 0))
mu = Function(Q) #

# Permeability
mu0 = Constant(mesh, PETSc.ScalarType(4e-7*np.pi))
muf = Constant(mesh, PETSc.ScalarType(1200*4e-7*np.pi))

# List of elements in each domain
core_elements = cell_tags.find(1) # core_tag = 1 
coil_elements = cell_tags.find(2) # coil_tag = 2
cyl_elements = cell_tags.find(3) # cyl_tag = 3 
air_elements = cell_tags.find(4) # air_tag = 4

# Defining the permeability values
mu.x.array[core_elements] = np.full_like(core_elements, muf, dtype=default_scalar_type)
mu.x.array[coil_elements] = np.full_like(coil_elements, mu0, dtype=default_scalar_type)
mu.x.array[cyl_elements] = np.full_like(cyl_elements, muf, dtype=default_scalar_type)
mu.x.array[air_elements] = np.full_like(air_elements, mu0, dtype=default_scalar_type)
```

### 3.4 Encontrando a densidade de corrente (J)
Como a seção transversal do fio AWG é conhecida, a densidade de corrente pode ser calculada de duas formas distintas. 
i) corrente do fio sobre a área do fio ou 
ii) corrente total (NI) sobre a área da bobina. Os dois resultados devem ser equivalentes se a relação área total igual a N vezes a área do fio for mantida.



```python
# J obtained from the wire's area:
j = I/awg.area
print('current density, J =',j)

# J obtained from the total area:
j = N*I/(Coil_thickness * z_core)
print('current density, J =',j)

# Creates a function J in the space Q that is non-zero only inside the coil
J = Function(Q)
Js = Constant(mesh, PETSc.ScalarType(j))
J.x.array[coil_elements] = np.full_like(coil_elements, Js, dtype=default_scalar_type)
```

    current density, J = 1214997.6878793903
    current density, J = 1214997.6878793903


Repare que os resultados são iguais. No entanto, em problemas reais, o isolamento do fio e o fato de as bobinas não se encaixarem perfeitamente fazem com que a área da bobina seja um pouco maior do que $N$ vezes a área do fio. De modo geral, o fator de enchimento (ou taxa de ocupação), que representa a proporção do espaço útil efetivamente preenchido pelo fio de cobre nos enrolamentos de um motor ou transformador, varia entre 0,30 e 0,75, dependendo do método de enrolamento e do isolamento utilizados. Como exemplo, a Figura 4 ilustra a construção manual de uma bobina de 2000 espiras usada em um protótipo de levitador. O lado esquerdo representa o enrolamento das primeiras 404 espiras, enquanto o lado direito corresponde à bobina com 1000 espiras. É possível notar que, no processo manual, a taxa de ocupação piora à medida que mais espiras são adicionadas.

Fig. 4 Construção manual de um protótipo de Eletroímã de 2000 espiras
<p>
    <img width="600" height="100" src="Figs/Aula_16_Bobina_manual.png">
</p>

## 4 - Construção da forma variacional

A forma variacional para esse problema é exatamente a mesma do problema axissimétrico da bobina no espaço livre, a única diferença são os materiais envolvidos.


```python
import ufl
# Define Trial and Test functions
u = TrialFunction(V)
v = TestFunction(V)

dx = ufl.Measure("dx", domain=mesh, subdomain_data=cell_tags);

# Variational equation with cylindrical coordinates
a = (1 / mu) * (1 / r) * dot(grad(u), grad(v)) * dx

# Source term of the axially symmetric weak form
L = J * v * dx

# Solving the linear problem
A_ = Function(V) # A_ = r*A_alpha
problem = LinearProblem(a, L, u=A_, bcs=[bc])
problem.solve()
```




    Coefficient(FunctionSpace(Mesh(blocked element (Basix element (P, triangle, 1, equispaced, unset, False, float64, []), (2,)), 0), Basix element (P, triangle, 2, gll_warped, unset, False, float64, [])), 2)



## 5 - Pós-processamento dos resultados

### 5.1 Potencial Vetor Magnético


```python
# Creates the plotter
plotter = pyvista.Plotter()

# Converts the mesh to an UnstructuredGrid
A_grid = pyvista.UnstructuredGrid(*vtk_mesh(V))

# Adds the A_z field data to the grid
A_grid.point_data["A_"] = A_.x.array
A_grid.set_active_scalars("A_")

# Applies the deformation based on the A_z field
warp = A_grid.warp_by_scalar("A_")

# Adds the deformed mesh to the plotter
actor = plotter.add_mesh(warp, show_edges=False)
actor = plotter.add_title(
    'Fig. 5 Potencial vetor magnético', font='courier', color='k', font_size=7)

plotter.view_xy()
plotter.show()

```

    error: XDG_RUNTIME_DIR is invalid or not set in the environment.
    MESA: error: ZINK: failed to choose pdev
    glx: failed to create drisw screen



    
![png](exemplo_professor_files/exemplo_professor_21_1.png)
    


### 5.2 Densidade de Fluxo Magnético ($\mathbf{B}$) e Campo Magnético ($\mathbf{H}$)


```python
# Compute Magnetic flux density (B = curl A)
W = functionspace(mesh, ("DG", 0, (mesh.geometry.dim, )))
B = Function(W)
B_expr = Expression(as_vector((-(1/r)*A_.dx(1), (1/r)*A_.dx(0))), W.element.interpolation_points())
B.interpolate(B_expr)

# Compute Magnetic Field (H = B/mu)
H = Function(W)
H_expr = Expression(as_vector((-(1/r/mu)*A_.dx(1), (1/r/mu)*A_.dx(0))), W.element.interpolation_points())
H.interpolate(H_expr)
```

#### 5.2.1 Gráfico de  ($\mathbf{B}$)


```python
# Iterpolate B again to mach vtk_mesh DoF.
Wl = functionspace(mesh, ("Lagrange", 2, (mesh.geometry.dim, )))
Bl = Function(Wl)
Bl.interpolate(B)

pyvista.start_xvfb()
topology, cell_types, geo = vtk_mesh(V)
values = np.zeros((geo.shape[0], 3), dtype=np.float64)
values[:, :len(Bl)] = Bl.x.array.real.reshape((geo.shape[0], len(Bl)))

# Create a point cloud of glyphs
function_grid = pyvista.UnstructuredGrid(topology, cell_types, geo)
function_grid["Bl"] = values
glyphs = function_grid.glyph(orient="Bl", factor=.3)

# Create a pyvista-grid for the mesh
mesh.topology.create_connectivity(mesh.topology.dim, mesh.topology.dim)
grid = pyvista.UnstructuredGrid(*vtk_mesh(mesh, mesh.topology.dim))

# Create plotter
plotter = pyvista.Plotter()
#plotter.add_mesh(grid, style="wireframe", color="k")
plotter.add_mesh(glyphs)
plotter.view_xy()
actor = plotter.add_title(
    'Fig. 6 Densidade de fluxo magnético', font='courier', color='k', font_size=7)
#plotter.window_size = [1000, 250];
#plotter.camera.zoom(3)
plotter.show();
```

    error: XDG_RUNTIME_DIR is invalid or not set in the environment.
    MESA: error: ZINK: failed to choose pdev
    glx: failed to create drisw screen



    
![png](exemplo_professor_files/exemplo_professor_25_1.png)
    


#### 5.2.2 Gráfico de  ($\mathbf{H}$)

Como $\mathbf{H} = \mathbf{B}/\mu$, podemos observar que o campo magnético tende a zero nos materiais com $\mu$ elevado, validando a aproximação analítica apresentada anteriormente. 


```python
# Iterpolate H again to mach vtk_mesh DoF.
Wl = functionspace(mesh, ("Lagrange", 2, (mesh.geometry.dim, )))
Hl = Function(Wl)
Hl.interpolate(H)

pyvista.start_xvfb()
topology, cell_types, geo = vtk_mesh(V)
values = np.zeros((geo.shape[0], 3), dtype=np.float64)
values[:, :len(Hl)] = Hl.x.array.real.reshape((geo.shape[0], len(Hl)))

# Create a point cloud of glyphs
function_grid = pyvista.UnstructuredGrid(topology, cell_types, geo)
function_grid["Hl"] = values
glyphs = function_grid.glyph(orient="Hl", factor=1.5e-6)

# Create a pyvista-grid for the mesh
mesh.topology.create_connectivity(mesh.topology.dim, mesh.topology.dim)
grid = pyvista.UnstructuredGrid(*vtk_mesh(mesh, mesh.topology.dim))

# Create plotter
plotter = pyvista.Plotter()
#plotter.add_mesh(grid, style="wireframe", color="k")
plotter.add_mesh(glyphs)
plotter.view_xy()
actor = plotter.add_title(
    'Fig. 7 Campo magnético', font='courier', color='k', font_size=7)
#plotter.window_size = [1000, 250];
#plotter.camera.zoom(3)
plotter.show();
```

    error: XDG_RUNTIME_DIR is invalid or not set in the environment.
    MESA: error: ZINK: failed to choose pdev
    glx: failed to create drisw screen



    
![png](exemplo_professor_files/exemplo_professor_27_1.png)
    


### 5.3 Indutância

A indutância pode ser calculada de duas maneiras distintas: i) Por meio da relação entre o fluxo total e a corrente de entrada ou ii) pela relação entre a energia total do sistema e a corrente. 

#### 5.3.1 Utilizando o fluxo magnético

Conforme apresentado anteriormente, podemos encontrar a indutância da bobina por meio da seguinte expressão:
\begin{equation}
L = \frac{N\phi}{I},
\end{equation}

Entretanto, é possível observar a partir da fig. 5 que $\mathbf{B}$ não é constante ao longo das bobinas. Isso fica mais evidente se traçarmos o gráfico do $|B|$ sobre o eixo vertical de simetria do problema.


```python
import matplotlib.pyplot as plt 
from dolfinx import geometry

#Creating a balanced tree with the elements of the mesh
bb_tree = geometry.bb_tree(mesh, mesh.topology.dim) 

# sample points
tol = 1e-6
z_points = np.linspace(-z_air/2 + tol, z_air/2-tol, 101)
r_points = np.full_like(z_points, 0)
points = np.zeros((3, 101))
points[0] = r_points
points[1] = z_points

B_values = []
cells = []
points_on_proc = []

# Find cells whose bounding-box collide with the the points
cell_candidates = geometry.compute_collisions_points(bb_tree, points.T)
# Choose one of the cells that contains the point
colliding_cells = geometry.compute_colliding_cells(mesh, cell_candidates, points.T)
for i, point in enumerate(points.T):
    if len(colliding_cells.links(i)) > 0:
        points_on_proc.append(point)
        cells.append(colliding_cells.links(i)[0])

points_on_proc = np.array(points_on_proc, dtype=np.float64)
B_values = B.eval(points_on_proc, cells)

# plot horizontal line
fig = plt.figure(figsize=(8,5))
plt.plot(z_points, B_values[:,0], 'y', linewidth=2); # magnify w
plt.plot(z_points, B_values[:,1], 'g', linewidth=2); # magnify w
plt.grid(True);

# Core and Cylinder boundaries
plt.plot(np.array([0,0]),np.array([0,np.max(B_values[:,1])]),'-.b', linewidth=1.5);
plt.plot(np.array([-x_gap,-x_gap]),np.array([0,np.max(B_values[:,1])]),'-.g', linewidth=1.5);
plt.plot(np.array([z_core,z_core]),np.array([0,np.max(B_values[:,1])]),'-.b', linewidth=1.5);
plt.plot(np.array([-x_gap-r_cyl,-x_gap-r_cyl]),np.array([0,np.max(B_values[:,1])]),'-.g',linewidth=1.5);


plt.legend(['Br (T)','Bz (T)', 'Core boudaries', 'Cylinder boundaries'], loc='upper left');
fig.suptitle("Fig. 8 Indução magnética");
plt.grid(True)
plt.xlabel('$z$ (m), r = 0')
plt.ylabel('$|B|$')


```




    Text(0, 0.5, '$|B|$')




    
![png](exemplo_professor_files/exemplo_professor_29_1.png)
    


A figura 8 mostra que, à medida que andamos verticalmente sobre o eixo de simetria ($r=0$), a intensidade do fluxo magnético varia, sendo máxima no centro do núcleo ferromagnético da bobina. Isso tem impacto no cálculo da indutância uma vez que espiras diferentes são enlaçadas por quantidades de fluxo diferentes. 

A figura 9 apresenta o componente vertical da densidade de fluxo magnético $B_z$ em função do raio $r$. Para essa figura, a altura $z$ é mantida fixa no centro da bobina. 

A partir da figura 9, é possível observar que o campo é relativamente constante dentro do núcleo, havendo uma transição abrupta na interface com o enrolamento. Essa transição indica que as linhas de fluxo percorrem o núcleo no sentido positivo de $z$ e retornam pelo ar (no sentido negativo de $z$). 

Repare que foi definida apenas uma variável, *Mesh.CharacteristicLengthMax* (ver Célula [3]), que controla a qualidade global da malha. Dessa forma, poucos elementos são usados para representar o enrolamento, resultando na característica pontiaguda da curva. Recomenda-se que o estudante teste diferentes densidades de malha, bem como diferentes tamanhos de retângulos de ar para ver o impacto na qualidade e custo computacional do problema.



```python
# Plot B
r_line = np.linspace(tol , r_air, 101)
z_line = z_core/2*np.ones(np.shape(r_line))

points = np.zeros((3, 101))
points[0] = r_line
points[1] = z_line

B_values = []
cells = []
points_on_proc = []

# Find cells whose bounding-box collide with the the points
cell_candidates = geometry.compute_collisions_points(bb_tree, points.T)
# Choose one of the cells that contains the point
colliding_cells = geometry.compute_colliding_cells(mesh, cell_candidates, points.T)
for i, point in enumerate(points.T):
    if len(colliding_cells.links(i)) > 0:
        points_on_proc.append(point)
        cells.append(colliding_cells.links(i)[0])

points_on_proc = np.array(points_on_proc, dtype=np.float64)
B_values = B.eval(points_on_proc, cells)

# plot horizontal line
fig = plt.figure(figsize=(8,5))
plt.plot(r_line, B_values[:,0], 'y', linewidth=2); # magnify w
plt.plot(r_line, B_values[:,1], 'g', linewidth=2); # magnify w
plt.grid(True);


plt.legend(['Br (T)','Bz (T)'], loc='upper right');
fig.suptitle("Fig. 9 Indução magnética");
plt.grid(True)
plt.xlabel('$r$ (m), z = r_core/2')
plt.ylabel('$B$')
```




    Text(0, 0.5, '$B$')




    
![png](exemplo_professor_files/exemplo_professor_31_1.png)
    


Uma vez conhecida a densidade de fluxo no interior da bobina, podemos encontrar a indutância. Uma aproximação simples consiste em assumir que $\mathbf{B}$ é constante no interior da bobina e pode ser aproximado pelo valor no centro do núcleo. Nesse caso temos:
\begin{eqnarray}
\phi &\approx& \int_0^{2\pi}\int_0^{r_{core}} \mathbf{B}\cdot rd\alpha dr\hat{z} \nonumber \\ 
     &\approx& 2\pi\int_0^{r_{core}} B_z r dr \nonumber \\
     &\approx& \pi r_{core}^2 B_z, \\
\end{eqnarray}

e a indutância é dada por
\begin{equation}
L = \frac{N\pi r_{core}^2 B_z}{I}. 
\end{equation}


```python
# Calc L using phi
Bz = B_values[0,1] # np.max(norm_B_v_line)
print('Bz = ', Bz, '(T)')
print('phi = ', np.pi*r_core**2*Bz, '(Wb)')
L_phi = N*np.pi*r_core**2*Bz/I
print('L = ', L_phi, '(H)')
```

    Bz =  0.021284443822654903 (T)
    phi =  6.686705234899729e-06 (Wb)
    L =  0.0006686705234899731 (H)


#### Utilizando a energia armazenada no Campo Magnético

A energia armazenada em um indutor na forma de Campo Magnético é dada por:
\begin{equation}
W = \frac{1}{2}LI^2. 
\end{equation}

Alternativamente, essa energia pode ser obtida integrando-se a densidade de energia magnética distribuída no volume. Para meios lineares temos:

\begin{equation}
W = \int \frac{1}{2}\mathbf{B}\cdot\mathbf{H} dv = \int \frac{1}{2\mu}|B|^2dv. 
\end{equation}

Lembrando que para o problema axissimétrico:
\begin{equation}
W = 2\pi \int_{-z_{air}}^{z_{air}}\int_{0}^{r_{air}} \frac{1}{2\mu}|B|^2rdrdz,
\end{equation}
cuja integral pode ser facilmente computada pelo *FEniCS* por meio da função *assemble()*. Dessa forma, a indutância pode ser estimada como: 

\begin{equation}
L = \frac{2W}{I^2}. 
\end{equation}


```python
Energy = np.pi*fem.assemble_scalar(fem.form(r/mu*dot(B,B)*dx))
L_w = 2*Energy/I**2

print(" Energy = {0:e} (J)". format(Energy))
print(" L = {0:e} (H)". format(L_w))
```

     Energy = 2.722820e-04 (J)
     L = 5.445640e-04 (H)


#### Comparando as duas abordagens

Reparem que ambas as abordagens apresentam resultados semelhantes. Contudo, a abordagem por meio do cálculo da energia não apresenta a simplificação de campo constante ao longo do cilindro. Adicionalmente, quando o comprimento da bobina é muito maior que sua seção transversal, e na ausência da peça cilíndrica, ambos os resultados tendem à aproximação analítica. 
\begin{equation}
L = \frac{\mu\pi r^2N}{L}.
\end{equation}


```python
print('Results:')
print('L(phi) = ', L_phi, '(H)')
print('L(W) = ', L_w, '(H)\n')
print('L(phi)/L(W) = ', L_phi/L_w, '\n')
L_ana = 1200*(4e-7*np.pi)*np.pi*r_core**2*N/z_core
print('Analytical Solution L = ', L_ana, '(H)')
print('L(W)/L(ana) = ', L_w/L_ana)
```

    Results:
    L(phi) =  0.0006686705234899731 (H)
    L(W) =  0.0005445639779876547 (H)
    
    L(phi)/L(W) =  1.2279007619287148 
    
    Analytical Solution L =  0.0009474820225045784 (H)
    L(W)/L(ana) =  0.5747486126946786


### Força magnética

Uma vez conhecido o campo magnético sobre a superfície da peça em levitação, a densidade de força por unidade de área exercida pelo campo pode ser obtida por meio do Tensor de Maxwell.

Integrando a densidade de força $d\textbf{F}$ sobre a superfície da peça obtemos a força total exercida pelo eletroímã. No caso da peça cilíndrica representada com simetria axissimétrica, temos três superfícies de interesse: i) a superfície superior da peça, ii) a superfície inferior e iii) a superfície lateral. Para cada superfície, o vetor normal $\hat{n}$ é igual a $+\hat{z}$, $+\hat{r}$ e $-\hat{z}$ respectivamente. A figura 10 apresenta os valores das componentes $H_z$ e $H_r$ ao longo de cada superfície. Reparem que o campo $H_z$ é predominante na superfície superior da peça, sendo essa a principal componente responsável pela força.


```python
# Lines on surfaces
tol2 = 1e-8
r_line = np.linspace(tol , r_cyl, 101)

# Points on each surface
points_up = np.zeros((3, 101))
points_down = np.zeros((3, 101))
points_side = np.zeros((3, 101))

# suf_up
points_up[0] = r_line
points_up[1] = (-x_gap+tol2)*np.ones(101)
# suf_down
points_down[0] = r_line
points_down[1] = (-x_gap-r_cyl-tol2)*np.ones(101)
#suf_side
points_side[0] = (r_cyl+tol2)*np.ones(101)
points_side[1] = -r_line-x_gap

points = np.concatenate((points_up, points_down, points_side), axis=1)

# H on surfaces
H_values = []
H_up = []
H_down = []
H_side = []

cells = []
points_on_proc = []

# Find cells whose bounding-box collide with the the points
cell_candidates = geometry.compute_collisions_points(bb_tree, points.T)
# Choose one of the cells that contains the point
colliding_cells = geometry.compute_colliding_cells(mesh, cell_candidates, points.T)
for i, point in enumerate(points.T):
    if len(colliding_cells.links(i)) > 0:
        points_on_proc.append(point)
        cells.append(colliding_cells.links(i)[0])

points_on_proc = np.array(points_on_proc, dtype=np.float64)
H_values = H.eval(points_on_proc, cells)

print(np.shape(H_values))
[H_up, H_down, H_side] = np.split(H_values, 3, axis=0)

# plot
fig = plt.figure(figsize=(15,4))
#print('Fig. 6 - Campo H nas superfícies da peça em levitação \n')
# suf_up
plt.subplot(1, 3, 1);
plt.plot(r_line, H_up[:,0],'b')
plt.plot(r_line, H_up[:,1],'g')
plt.legend(['$H_r$ (A/m)','$H_z$ (A/m)'], loc='center left');
plt.grid(True)
plt.xlabel('$r$ (m)')
plt.ylabel('$H$')
plt.title('Fig. 10 - Campo H nas superfícies da peça em levitação \n a) superfície superior');

# suf_down
plt.subplot(1, 3, 2);
plt.plot(r_line, H_down[:,0],'b')
plt.plot(r_line, H_down[:,1],'g')
#plt.legend(['$H_r$ (A/m)','$H_z$ (A/m)'], loc='center left');
plt.grid(True)
plt.xlabel('$r$ (m)')
#plt.ylabel('$H$')
plt.title('b) superfície inferior');

# suf_side
plt.subplot(1, 3, 3);
plt.plot(-x_gap-r_line, H_side[:,0],'b')
plt.plot(-x_gap-r_line, H_side[:,1],'g')
#plt.legend(['$H_r$ (A/m)','$H_z$ (A/m)'], loc='center left');
plt.grid(True)
plt.xlabel('$z$ (m)')
#plt.ylabel('$H$')
plt.title('b) superfície lateral');
```

    (303, 2)



    
![png](exemplo_professor_files/exemplo_professor_39_1.png)
    


#### Integral na superfície superior

Na superfície superior $\hat{n} = \hat{z}$, resultando em: 
\begin{eqnarray}
\frac{d\textbf{F}}{ds} &=& -\frac{\mu_0}{2}\left(H_r^2+H_z^2\right)\hat{z}+\mu_0 H_z\left(H_r\hat{r}+H_z\hat{z}\right) \nonumber\\
&=& \mu_0 H_zH_r\hat{r} + \frac{\mu_0}{2}\left(H_z^2-H_r^2\right)\hat{z}.
\end{eqnarray}

Consequentemente, a parcela da força total referente à superfície superior será igual a:
\begin{eqnarray}
\mathbf{F_{sup}}&=& \int_0^{2\pi}\int_0^{r_{cyl}} \left[ \mu_0 H_z H_r\hat{r} + \frac{\mu_0}{2}\left(H_z^2-H_r^2\right) \hat{z}\right]  rd\alpha dr \nonumber \\ 
&=& 2\pi\mu_0 \int_0^{r_{cyl}} \frac{1}{2}(H_z^2-H_r^2) r dr \hat{z}
\end{eqnarray}
onde a componente $F_r$ é eliminada da integral devido à natureza axissimétrica do problema, ($H_r(\alpha)=-H_r(-\alpha))$. 

Avaliando a equação:


```python
# using scipy to integ
mi0 = 4e-7*np.pi
df_z = np.pi*mi0*(H_up[:,1]**2 -H_up[:,0]**2)*r_line
F_sup = integrate.simpson(df_z, x=r_line)
print('F_sup = ', F_sup ,'(N)')
```

    F_sup =  0.0027247634111716356 (N)


#### Integral na superfície inferior

Na superfície superior $\hat{n} = -\hat{z}$, logo:
\begin{eqnarray}
\frac{d\textbf{F}}{ds} &=& +\frac{\mu_0}{2}\left(H_r^2+H_z^2\right)\hat{z}-\mu_0 H_z\left(H_r\hat{r}+H_z\hat{z}\right) \nonumber\\
&=& -\mu_0 H_zH_r\hat{r} + \frac{\mu_0}{2}\left(H_r^2-H_z^2\right)\hat{z}.
\end{eqnarray}

Consequentemente, a parcela da força total referente à superfície inferior será igual a:
\begin{eqnarray}
\mathbf{F_{inf}} &=& -2\pi\mu_0 \int_0^{r_{cyl}} \frac{1}{2}(H_z^2-H_r^2) r dr\hat{z}
\end{eqnarray}

Avaliando a equação:


```python
df_z = -np.pi*mi0*(H_down[:,1]**2 -H_down[:,0]**2)*r_line
F_inf = integrate.simpson(df_z, x=r_line)
print('F_inf = ', F_inf ,'(N)')
```

    F_inf =  -0.0003794667817811139 (N)


#### Integral na superfície lateral

Na superfície superior $\hat{n} = \hat{r}$, resultando em: 
\begin{eqnarray}
\frac{d\textbf{F}}{ds} &=& -\frac{\mu_0}{2}\left(H_r^2+H_z^2\right)\hat{r}+\mu_0 H_r\left(H_r\hat{r}+H_z\hat{z}\right) \nonumber\\
&=& \frac{\mu_0}{2}\left(H_r^2-H_z^2\right)\hat{z} +\mu_0 H_rH_z\hat{z}.
\end{eqnarray}

Consequentemente, a parcela da força total referente à superfície lateral será igual a:
\begin{eqnarray}
\mathbf{F_{lat}}&=& \int_0^{2\pi}\int_{-x}^{-x-z_{cyl}} \left[\frac{\mu_0}{2}\left(H_r^2-H_z^2\right)\hat{r} + \mu_0 H_r H_z\hat{z} \right]r_{cyl}d\alpha dz  \nonumber \\ 
&=& 2\pi\mu_0 \int_{-x}^{-x-z_{cyl}} H_rH_z r_{cyl} dz \hat{z}
\end{eqnarray}
onde a componente $F_r$ é eliminada da integral devido à natureza axissimétrica do problema. 

Avaliando a equação:


```python
df_z = 2*np.pi*mi0*(H_side[:,0]*H_side[:,1])*r_cyl
F_lat = integrate.simpson(df_z, x=r_line)
print('F_lat = ', F_lat ,'(N)')
```

    F_lat =  -0.00010861654666246493 (N)


#### Força total

Como resultado, a força total é a soma algébrica das contribuições individuais das superfícies. Como esperado, o termo proveniente da superfície superior é predominante.


```python
F_tot = F_sup + F_inf  + F_lat 
print('F_total = ', F_tot ,'(N)')
F_femm = 0.00202562 # result using the FEMM software.
print('Erro (%) em relação ao software FEMM: ', 100*(F_tot-F_femm)/F_femm)
```

    F_total =  0.0022366800827280566 (N)
    Erro (%) em relação ao software FEMM:  10.419529957645386


#### Cálculo do tensor de Maxwell usando as funções *FacetNormal()* e *assemble()*

O cálculo da força total, conforme apresentado anteriormente, foi possível devido à geometria bem definida da peça em levitação. Entretanto, ao substituir a peça cilíndrica por uma esférica, seria necessário trocar as três integrais de superfície por uma integral ao longo do arco que representa a projeção bidimensional da esfera. Felizmente, existe um método automatizado para calcular integrais de componentes normais (ou tangenciais) de um campo vetorial nas interfaces do domínio. Esse método utiliza as funções *FacetNormal()* e *assemble()*.

A função *FacetNormal()* retorna o vetor normal às fronteiras de uma determinada malha, enquanto a função *assemble()* realiza integrações numéricas sobre um domínio ($dx$) ou sobre sua fronteira ($ds$).

Para ilustrar o uso da função *assemble()* e a definição da fronteira $ds$, utilizaremos essas ferramentas para calcular a soma dos perímetros das três superfícies de integração mencionadas na seção anterior: (i) a superfície superior da peça, (ii) a superfície inferior e (iii) a superfície lateral.

O primeiro passo para utilizar a função assemble() é construir o elemento de integração ($ds$), que define quais arestas da malha compõem a fronteira externa da peça cilíndrica. Existem diferentes abordagens para construir esse elemento. Neste exemplo, utilizaremos a função *locate_entities()* para identificar as arestas na superfície da peça. Essa informação é armazenada em um objeto do tipo *meshtag*, juntamente com uma tag de identificação. Alternativamente, o processo de identificação da região de interesse poderia ser realizado previamente no GMSH por meio dos *physical groups*.


```python
from dolfinx.mesh import locate_entities, meshtags

# Boolean functions to find the edges of the cylinder
def on_top(x): 
    return np.logical_and(x[0]<=r_cyl,np.isclose(x[1], -x_gap))
def on_bottom(x):
    return np.logical_and(x[0]<=r_cyl,np.isclose(x[1], -x_gap-z_cyl))
def on_side(x):
    return np.logical_and(np.isclose(x[0], r_cyl), np.logical_and(x[1]<=-x_gap,x[1]>=-x_gap-z_cyl))
    
# Finding the edges on the cylinder surfaces using the locate_entitiesy function
top_facets = locate_entities( # list of all edges  on the upper surface
    mesh,gdim -1, on_top) 
bottom_facets = locate_entities( # list of all edges on the bottom surface
    mesh,gdim -1, on_bottom)
side_facets = locate_entities( # list of all edges on the side surface
    mesh,gdim -1, on_side)

cyl_facets = top_facets
cyl_facets = np.append(cyl_facets, bottom_facets)
cyl_facets = np.append(cyl_facets, side_facets)

# Creating a label for each surface
# Defining a marker equal to "1" for the top edges, "2" for the bottom and "3" for the side ones. 
facet_markers = np.ones(np.shape(top_facets),dtype=np.int32)
facet_markers = np.append(facet_markers,2*np.ones(np.shape(bottom_facets),dtype=np.int32))
facet_markers = np.append(facet_markers,3*np.ones(np.shape(side_facets),dtype=np.int32))

# a dolfinx.meshMes class that associate data (markers) with a subset of mesh entities of a given dimension.
facet_tag = meshtags(mesh, gdim-1, cyl_facets, facet_markers)

# Defining a surface integration element that operates only on the defined edges
ds = ufl.Measure("dS", domain=mesh, subdomain_data=facet_tag)
```

O elemento $d_s$ possui a informação de todas as arestas selecionadas e suas tags. Como as arestas devidamente identificadas, podemos criar o elemento de integração e usar a função *assemble()* para calcular o perímetro. Reparem que usaremos apenas os elementos de $ds$ marcados com os índices 1-3.


```python
# Function assemble simple test
perimeter = fem.assemble_scalar(fem.form(1*ds(1))) 
perimeter2 = fem.assemble_scalar(fem.form(1*ds(2))) 
perimeter3 = fem.assemble_scalar(fem.form(1*ds(3))) 

print('perimeter = ', perimeter, '(m)')
print('perimeter2 = ', perimeter2, '(m)')
print('perimeter3 = ', perimeter3, '(m)')
```

    perimeter =  0.01 (m)
    perimeter2 =  0.01 (m)
    perimeter3 =  0.009999999999999998 (m)


#### Cálculo da Força

Como a função *assemble()* retorna um escalar, a densidade de força $d\textbf{F}$ deve ser decomposta em suas componentes $dF_r$ e $dF_z$. Por simetria, $dF_r=0$, sendo necessário apenas encontrar a componente $dF_z$. 

\begin{eqnarray}
%{dF_r} &=& 2\pi\left[ -\frac{\mu_0}{2}(H_r^2+H_z^2)n_r + \mu_0(H_r\cdot n_r+ H_z\cdot n_z)H_r\right]rds \\
{dF_z} &=& \left[ -\frac{\mu_0}{2}(H_r^2+H_z^2)n_z + \mu_0(\textbf{H}\cdot \hat{n})H_z\right]ds 
\end{eqnarray}
onde $n_r$ e $n_z$ são as projeções do vetor normal unitário, obtidas por meio da função *FacetNormal()*.


```python
n = ufl.FacetNormal(mesh) 

F_z = list()
for i in range(3):
    dF_z = fem.form(2*np.pi*mi0*r*( -0.5*( H[0]('-')**2+H[1]('-')**2 )*n[1]('+') + dot(H('-'),n('+'))*H[1]('-') )*ds(i+1) + Constant(mesh, PETSc.ScalarType(0.0))*dx(3));
    F_z.append(fem.assemble_scalar(dF_z));
    print('Fz',i+1,'= ', F_z[i] ,'(N)')

F_tot2 = F_z[0]+F_z[1]+F_z[2]
print('F_tot = ', F_tot2,'(N)\n')

print('Comparando com o resultado anterior:\n')
print('Erro (%) = ', 100*(F_tot2-F_tot)/F_tot)
```

    Fz 1 =  0.0027223210217848922 (N)
    Fz 2 =  -1.5622896071994275e-10 (N)
    Fz 3 =  2.143042949954794e-09 (N)
    F_tot =  0.0027223230085988816 (N)
    
    Comparando com o resultado anterior:
    
    Erro (%) =  21.71266823633048


Repare que a função de integração é a mesma para cada interface, o que muda é apenas o elemento de integração $ds$. Entretanto, o resultado está correto apenas para a primeira integral.

Alguns detalhes merecem ser observados: i) A função Measure() recebe o argumento 'dS' (S maiúsculo), que indica que as interfaces que estamos definindo estão dentro do domínio. Logo, a interface pode pertencer simultaneamente a dois elementos. Quando usamos o argumento 'ds' (s minúsculo), apenas a fronteira externa do domínio é selecionada, e o vetor normal é definido para fora da malha. ii) É necessário definir a orientação das grandezas vetoriais por meio dos argumentos ('+') e ('-'), uma vez que, dentro do domínio, não temos uma definição de dentro ou fora para o vetor normal. Os argumentos ('+') e ('-') também são usados no Campo Magnético para definir em qual elemento iremos avaliar o campo. Veja na Fig. 7 que o campo $H$ é descontínuo na interface entre a peça e o ar. No caso do tensor de Maxwell, estamos interessados no Campo Magnético exterior à peça em levitação. Dependendo da escolha dos argumentos de orientação, os resultados das integrais serão diferentes. Recomenda-se que o aluno avalie diferentes configurações.

Para resolver a ambiguidade do Campo Magnético e da orientação do vetor normal nas versões antigas do FEniCS, a função assemble() recebia um argumento adicional (*Constant(0)dx(3)) que indica que estamos integrando também uma constante nula no subdomínio 3 (peça cilíndrica). Embora o valor dessa integral seja nulo, ela era usada para padronizar a orientação do vetor normal. Dessa forma, a orientação ('+') é definida 'para fora' da peça cilíndrica e é usada juntamente com o vetor normal $n$. No caso do campo $\textbf{H}$, deveríamos usar a orientação ('-') para avaliar o campo externo à peça. Entretanto, esse tipo de solução não está disponível na versão 0.9.0.0 usada neste notebook.

Para essa versão, iremos construir um elemento de integração 'ds' (s minúsculo). Por ser um elemento de integração externo com cada aresta associada apenas a um único elemento. Nesse caso, teremos que criar essa associação manualmente, seguindo os passos descritos <em> <a href="https://fenicsproject.discourse.group/t/wrong-facetnormal-vector-on-internal-boundaries/12887/2"> aqui </a>.


```python
# Create connectivties required for defining integration entities
mesh.topology.create_entities(gdim-1)
mesh.topology.create_connectivity(gdim-1, gdim)
mesh.topology.create_connectivity(gdim, gdim-1)

# Give a set of facets, get a consistent orientation for an interior integral
facets_to_integrate = facet_tag.find(1)
facets_to_integrate = np.append(facets_to_integrate, facet_tag.find(2))
facets_to_integrate = np.append(facets_to_integrate, facet_tag.find(3))
#print('facets_to_integrate', facets_to_integrate)

# f_to_c generates, for each edge, the list of elements in which it appears.
f_to_c = mesh.topology.connectivity(gdim-1, gdim)
# c_to_f performs the opposite operation.
c_to_f = mesh.topology.connectivity(gdim, gdim-1)

#From f_to_c, we can see that each edge is shared by two elements, and we only want the ones belonging to the air region (air_tag = 4).

# Compute integration entities for a single facet of a cell.
# Each facet is represented as a tuple (cell_index, local_facet_index), where cell_index is local to process
# local_facet_index is the local indexing of a facet for a given cell
integration_entities = []
for i, facet in enumerate(facets_to_integrate):
    # Find cells connected to facet
    cells = f_to_c.links(facet)
    # Get value of cells
    marked_cells = cell_tags.values[cells]
    # Get the cell marked with 4
    correct_cell = np.flatnonzero(marked_cells == 4)
    
    assert len(correct_cell) == 1
    # Get local index of facet
    local_facets = c_to_f.links(cells[correct_cell[0]])
    local_index = np.flatnonzero(local_facets == facet)
    assert len(local_index) == 1

    # Append integration entities
    integration_entities.append(cells[correct_cell[0]])
    integration_entities.append(local_index[0])

#print('integration_entities', integration_entities)
# Create custom integration measure for one-sided integrals
ds2 = ufl.Measure("ds", domain=mesh, subdomain_data=[
                 (1, np.asarray(integration_entities, dtype=np.int32))])
```

Com o novo elemento de integração, a dualidade na definição dos termos é eliminada, e a integral se resume a:


```python
dF_z = fem.form(2*np.pi*mi0*r*( -0.5*( H[0]**2+H[1]**2 )*-n[1] + dot(H,-n)*H[1])*ds2(1))
F_tot3 = fem.assemble_scalar(dF_z)

print('F_tot = ', F_tot3,'(N)\n')
print('Comparando com o primeiro resultado:')
print('Erro (%) = ', 100*(F_tot3-F_tot)/F_tot)
```

    F_tot =  0.0022313614878044555 (N)
    
    Comparando com o primeiro resultado:
    Erro (%) =  -0.23778970290261833


## 6 - Comparação com o modelo circuital

### Indutância e força magnética

Assim como no exemplo analítico, a solução do problema possibilita o cálculo dos parâmetros usados no modelo ODE do problema. Uma vez parametrizado o problema, pode-se resolver o modelo para diferentes posições de $x$ a fim de se levantar as curvas $L(x)$ e $f_m(i,x)$. Reparem que como o modelo é linear, não é necessário varia a corrente durante as simulações. O mesmo não é verdade se considerarmos a curva $B\times H$ do material ferromagnético.  

### Constante $k$ 

Um ponto interessante a ser observado é que a relação entre os parâmetros circuitais e os parâmetros mecânicos foram definidos a parti de uma constante de proporcionalidade $k$ onde:
\begin{eqnarray}
L(x)&=& \frac{k}{x} \\ 
f_m(I,x) &=& \frac{k}{2}\left(\frac{I}{x}\right)^2
\end{eqnarray}

Reescrevendo as equações:
\begin{eqnarray}
k &=& xL(x) \\ 
k &=& 2f_m\left(\frac{x}{I}\right)^2
\end{eqnarray}

Todavia, usando os resultados para encontrar k obtemos:


```python
k1=L_w*x_gap
print('k(L) = ',k1)

k2 = 2*F_tot*x_gap**2/I**2
print('k(f_m) = ', k2)
print('k(L)/k(f_m) = ', k1/k2)
```

    k(L) =  2.7228198899382737e-06
    k(f_m) =  1.1183400413640283e-07
    k(L)/k(f_m) =  24.346976672830895


A justificativa para tal diferença encontra-se no fato que as relações acima foram obtidas considerando que todas as linhas de campo geradas pela bobina passarão pela peça em levitação. Isso foi verdade para o levitador analítico e tanto a indutância quando a força magnética puderam ser representadas por um $k$ constante. Contudo, para o levitador axissimétrico, uma parcela das linhas de campo circula o eletroímã sem entrar em contato com a peça em levitação. Essa parcela de fluxo contribui para a indutância, mas não impacta a força magnética. Isso justificaria as constantes diferentes.

### Resistência

Finalmente, a resitência do levitador axissimétrico pode ser obtida aproximando-se o comprimento médio das espiras da bobina. 


```python
# Average loop length
all = 2*np.pi*(r_core+Coil_thickness/2)
# Total length
L = N*all
R = L/awg.sigma/awg.area
print('R = ',R, '(ohms)')
```

    R =  0.1424547267492916 (ohms)


Alternativamente, a resitência pode ser calculada por meio da potência dissipada no enrolamento. Da teoria de circuitos, sabemos que a potência dissipada em um condutor é dada por $P = RI^2$. Onde a potência pode ser obtida por meio da integral da densidade de potência na região da bobina ($\Omega_c$):
\begin{equation}
P = \int_{\Omega_c} \frac{j^2}{\sigma} d\Omega_c = 2\pi z\int_{r_{int}}^{r_{ext}} \frac{j^2}{\sigma} rdr
\end{equation}
onde $z$ é a altura da bobina e $r_{int}$ e $r_{ext}$ são respectivamente os raios interno e externo da bobina.


```python
P = 2*np.pi*fem.assemble_scalar(fem.form(r*J**2/awg.sigma*dx(2)))
print('P = ',P, '(W)')
print('R = ',P/I**2, '(ohms)')
```

    P =  0.14245472674929152 (W)
    R =  0.14245472674929152 (ohms)


#### Comparação entre as aproximações para a resistência

Calculando o erro (%) entre a resistência obtida pela aproximação por um fio desenrolado $R_{fio}$ e a resistência obtida por meio da integral da densidade de potência $R_{int}$.
\begin{equation}
erro = \frac{R_{int} - R_{fio}}{R_{int}}\cdot100\%
\end{equation}




```python
print('erro = ',(P/I**2 - R)/(P/I**2)*100, '(%)')
```

    erro =  -5.845136117766682e-14 (%)


O erro obtido é da ordem da precisão da máquina, o que demonstra que os resultados são equivalentes. Esse comportamento já era esperado, uma vez que a geometria foi projetada para impor um fator de enchimento igual a um.

Uma vez obtidos a indutância, o fator de força e a resistência, é possível simular o modelo circuital do levitador e criar um gêmeo digital sem necessidade de construir o protótipo. Com o modelo, podem-se avaliar mecanismos de controle e prever o desempenho do dispositivo. As vantagens do gêmeo digital incluem redução de custos, agilidade no teste de controladores, possibilidade de otimização iterativa e eliminação de riscos físicos. Contudo, para que o gêmeo digital funcione corretamente, é essencial que as propriedades magnéticas dos materiais sejam conhecidas com precisão.
