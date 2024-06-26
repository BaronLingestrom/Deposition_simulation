# Deposition simulation - WIP
by Robin Kucsera and Zoltan Santha, 2024\
project wotk for BME Computer Simulation in Physics course (BMETE15MF74)

## Description

## Usage

## External libraries
- NumPy
- SciPy
- Matplotlib
- (copy module)

## Implementation

### Global variables

### Projection functions

### Classes

#### Particle
- Properties:
  - pos: &nbsp;
  postition in lattice matrix (3D)
  - E: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  kinetic energy
- Methods:
  ~~~python
  class Particle():
    def __init__(self, pos, E):
      ...
      return
    def __copy__(self):
      ...
      return
    def __repr__(self):
      ...
      return
  ~~~

#### Lattice
- Properties:
  - size: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  x,y dimensions of a layer (lattice size)
  - lattice_template: &nbsp;
  template array for adding new layer
  - neighbors: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  list of vectors pointing to (nearest) neighboring sites (matrix represention)
  - lattice:
  - heightmap:
  - E_bond:
  - E_bound:
  - E_substrate:
  - particles:
