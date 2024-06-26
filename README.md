# Deposition simulation
by Robin Kucsera and Zoltan Santha, 2024\
project wotk for BME Computer Simulation in Physics course (BMETE15MF74)

## Description

The project is about thin-film deposition using evaporation. We concider a small target and a vapor source located far away. This results in a uniform spatial distribution of atoms originating from the beam, simplifying the model. The evaporated atoms form a gas resulting in that the incoming particles have an energy distribution according to the kinetic theory of gases:

$f(E) = 2\sqrt{\frac{E}{\pi}}(\frac{1}{k_BT})^{3/2}\exp(-\frac{E}{k_BT})$,

which can be written equivalently as a gamma distribution, with $3/2$ shape factor and $k_BT$ scale parameter.

The implementation builds on the kinetic Monte-Carlo algorithm. The rate of the particle hopping in 3D is determined by

$w_{i\rightarrow j}=(\frac{E_{kin}}{T})^2\exp(\frac{E_{int,j}-E_{int,i}}{T})$,

where $E_{kin}$ is the kinetic energy of the particle, $E_{int,i}$ and $E_{int,j}$ is the interaction energy of the particle at position $i$ and $j$ respectively. The particle desorption is considered as a hopping to a site with no neighboring particles, so this effect is also taken into consideration. We also implemented a kinetic energy transfer for moving particles to model the thermal effects of a particle movement. Each particle movement results in the kinetic energy change of the nearby particles according to

$\Delta E_i=\frac{\bar{E}(\bar{E}-E_i)}{\sum_j|\bar{E}-E_j|}$,

where $\bar{E}$ is the average kinetic energy of the considered (nearby) particles.

The output of the simulation is GIF of four time evolving plots: the atomic deposition in 3D, the 2D heightmap of the depostied particles, the number of particles in each layer and the kinetic energy distribution of the deposited particles.

## Usage
The simulation is written in one single python script. After changing the global varibles and projection functions, the script can be simply run from a terminal using python. A GIF is saved after the simulation.

## External libraries
- NumPy
- SciPy
- Matplotlib
- (copy module)

## Implementation

### Global variables
  ~~~python
  SIZE                          # x,y matrix dimensions of a layer
  LATTICE                       # matrix for a layer
  NEIGHBORING_SITES             # list of vectors pointing to neighboring sites
  TEMPERATURES                  # list of temperatures
  HOPPING_PER_TEMPERATURE       # hoppings to perform at each temperature
  NEW_PARTICLE_PER_TEMPERATURE  # new particles to add at each temperature
  E_BOND                        # intralayer bonding energy (same layer)
  E_BOUND                       # interlayer bonding energy (diff. layers)
  E_SUBSTRATE                   # substrate bonding energy (lowest layer to substrate)
  FRAME_RATE                    # actions between each snapshot
  ~~~

### Projection functions
Transformation functions from matrix representation to real space representation.

Inputs: lists of coordinates for particles (particle i: x[i], y[i], z[i])

#### Example: (triangle lattice)
~~~python
def project_x(x,y,z):
    x_out = []
    for i in range(len(x)):
        x_out += [x[i]+y[i]/2]
    return x_out
def project_y(x,y,z):
    y_out = []
    for i in range(len(y)):
        y_out += [y[i]*np.sqrt(3)/2]
    return y_out
def project_z(x,y,z):
    z_out = []
    for i in range(len(z)):
        z_out += [z[i]]
    return z_out
~~~

### Main function
At each temperature value HOPPING_PER_TEMPERATURE number of hoppings and NEW_PARTICLE_PER_TEMPERATURE number of particle additions are performed in random sequence.

  ~~~python
  def main():
      # main function
      lattice = Lattice(SIZE,LATTICE,NEIGHBORING_SITES,E_BOND,E_BOUND,E_SUBSTRATE)

      particles = []; heightmap = [] # for snippets during simulation -> animation
      for i,T in enumerate(TEMPERATURES):
          print("Current temperature: T = {:.3f}\t ({}/{})".format(T,i+1,len(TEMPERATURES)))
          action_sequence = np.array([0]*HOPPING_PER_TEMPERATURE + [1]*NEW_PARTICLE_PER_TEMPERATURE)
          np.random.shuffle(action_sequence)
          # action: hopping(0) or incident particle(1) happens
          print("  Action sequence: {}".format(action_sequence))
          for j,action in enumerate(action_sequence):
              if action == 0:     # hopping action
                  lattice.hopping_action(T)
              else:               # new incident particle
                  lattice.new_particle_action()
              # animation sampling
              if ( i*(HOPPING_PER_TEMPERATURE+NEW_PARTICLE_PER_TEMPERATURE)+j )%FRAME_RATE == 0:
                  particles.append(copy.deepcopy(lattice.particles))
                  heightmap.append(copy.deepcopy(lattice.heightmap))
          print("  Current particle number: {}".format(len(lattice.particles)))

      animation(particles,heightmap)
      return
  ~~~

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
    def __copy__(self):
      return type(self)(self.pos,self.E)
    def __repr__(self):
      return str(self.pos)
  ~~~

#### Lattice
- Properties:
  - size: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  x,y dimensions of a layer (lattice size)
  - lattice_template: &nbsp;
  template array for adding new layer
  - neighbors: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  list of vectors pointing to (nearest) neighboring sites (matrix represention)
  - lattice: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  3D matrix of sites (occupation information)
  - heightmap: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  highest occuped layer number at each x,y lattice site
  - E_bond: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  bonding energy intralayer (within a layer)
  - E_bound: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  bonding energy interlayer (between layers)
  - E_substrate: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  bonding energy to the substrate
  - particles: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  list of particles (Particle() class objects)

- Methods:
  ~~~python
  class Lattice():
    def __init__(self,size,lattice,neighbors,E_bond,E_bound,E_substrate):
      ...
      return
    def extend_layers(self):
      # if the top layer has a particle, then adds a new empty layer on top
      ...
      return
    def get_max_height(self,pos_2d):
      # returns the highest occupied layer number at pos_2d
      ...
      return top_particle_layer_idx + 1
    def hopping_sites(self,pos):
      # returns the possible hopping destinations from pos
      ...
      return pos_to_hopp
    def intralayer_interaction_sites(self,pos):
      # returns the sites with intralayer(within layer) interaction to pos site
      ...
      return pos_of_sites
    def interlayer_interaction_sites(self,pos):
      # returns the sites with interlayer(between layers) interaction to pos site
      ...
      return pos_of_sites
    def all_interaction_sites(self,pos):
      # returns all sites with interaction to pos site
      return np.append(self.intralayer_interaction_sites(pos),self.interlayer_interaction_sites(pos),axis=0)
    def add_particle(self,particle):
      # adds particle to the lattice
      ...
      return particle.pos
    def sub_particle(self,index):
      # removes the particle with index (in particles) from the lattice
      ...
      return
    def get_particle_index(self,pos):
      # returns the index (in particles) of particle at pos
      ...
      return i
    def generate_new_particle(self):
      # generates a new particle with Maxwell-Boltzmann energy distribution and uniform spatial distribution
      ...
      return Particle(pos, E)     # position in 2D
    def calculate_interaction_energy(self,pos):
      # calculates the interaction energy of particle at pos
      ...
      return E
    def calculate_rate_factors(self, T):
      # calculates the hopping rates of the particles
      ...
      return origin_list,destination_list,rate_list
    def movement_energy_transfer(self,pos):
      # distributes energy around pos due to particle movement
      ...
      return
    def hopping_energy_transfer(self,origin,destination):
      # updates energies due hopping (changed bond structure)
      ...
      return
    def hopping_action(self,T):
      # performs a particle hopping at T temperature
      ...
      return
    def new_particle_action(self):
      # adds a new particle to the system
      ...
      return
  ~~~

### Animation functions
  ~~~python
  def animation(particles,heightmap):
    # creates animation from snippets made during simulation
    ...
    return
  def update(frame):
    # plot a single snippet
    ...
    return scat,hmap,ppl,en_dist,
  ~~~
