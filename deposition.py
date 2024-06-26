# Deposition simulation
# by Robin Kucsera and Zoltan Santha, 2024
# project work for BME Computer Simulation in Physics course (BMETE15MF74)

# EXTERNAL LIBRARIES
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import copy

# GLOBAL VARIABLES
SIZE = np.array([10,10])
LATTICE = np.zeros(SIZE)
#   vectors to neighboring sites
NEIGHBORING_SITES = [np.array([0,0,1]),np.array([0,0,-1]),
                     np.array([1,0,0]),np.array([0,1,0]),
                     np.array([-1,0,0]),np.array([0,-1,0]),
                     np.array([1,-1,0]),np.array([-1,1,0])] # triangle lattice
TEMPERATURES = 1/(np.linspace(0.1,10,10))
HOPPING_PER_TEMPERATURE = 3
NEW_PARTICLE_PER_TEMPERATURE = 10
E_BOND = 1                                  # intralayer bonding energy (same layer)
E_BOUND = 1                                 # interlayer bonding energy (diff. layers)
E_SUBSTRATE = 3                             # substrate bonding energy (lowest layer to substrate)

FRAME_RATE = 5  # actions between each snapshot
# plot projections
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

# SOURCE CODE
class Particle():
    def __init__(self, pos, E):
        self.pos = pos
        self.E = E
    def __copy__(self):
        return type(self)(self.pos,self.E)
    def __repr__(self):
        return str(self.pos)

class Lattice():
    def __init__(self,size,lattice,neighbors,E_bond,E_bound,E_substrate):
        self.size = size
        self.lattice_template = np.expand_dims(np.zeros_like(copy.deepcopy(lattice),dtype=bool),axis=2)
        self.neighbors = copy.deepcopy(neighbors)
        self.lattice = np.append(copy.deepcopy(self.lattice_template),copy.deepcopy(self.lattice_template),axis=2)
        self.heightmap = np.zeros(self.size,dtype=int)
        self.E_bond = E_bond
        self.E_bound = E_bound
        self.E_substrate = E_substrate
        self.particles = []
        return
    
    def extend_layers(self):
        if np.max(self.heightmap)==self.lattice.shape[2]-1:
            self.lattice = np.append(self.lattice,copy.deepcopy(self.lattice_template),axis=2)
        return
    
    def get_max_height(self,pos_2d):
        top_particle_layer_idx = max([i for i in range(self.lattice.shape[2]) if self.lattice[*pos_2d,i]],default=0) # 0 for 1st layer particle
        return top_particle_layer_idx + 1
    
    def hopping_sites(self,pos):
        pos_to_hopp = pos + np.array([v for v in self.neighbors if not np.array_equal(v,[0,0,0]) and pos[2]+v[2]>=0])
        pos_to_hopp[:,0] %= self.size[0]
        pos_to_hopp[:,1] %= self.size[1]
        return pos_to_hopp
    
    def intralayer_interaction_sites(self,pos):
        pos_of_sites = pos + np.array([v for v in self.neighbors if not np.array_equal(v,[0,0,0]) and v[2]==0])
        pos_of_sites[:,0] %= self.size[0]
        pos_of_sites[:,1] %= self.size[1]
        return pos_of_sites
    
    def interlayer_interaction_sites(self,pos):
        pos_of_sites = pos + np.array([v for v in self.neighbors if not np.array_equal(v,[0,0,0]) and v[2]!=0 and pos[2]+v[2]>=0])
        pos_of_sites[:,0] %= self.size[0]
        pos_of_sites[:,1] %= self.size[1]
        return pos_of_sites
    
    def all_interaction_sites(self,pos):
        return np.append(self.intralayer_interaction_sites(pos),self.interlayer_interaction_sites(pos),axis=0)
    
    def add_particle(self,particle):
        pos = particle.pos # 2D
        particle.pos = [*pos, self.heightmap[*pos]] # 3D
        self.heightmap[*pos]+=1
        self.lattice[*particle.pos] = True
        self.particles.append(particle)
        self.extend_layers()
        return particle.pos
    
    def sub_particle(self,index):
        self.lattice[*self.particles[index].pos] = False
        pos = self.particles[index].pos[:2] # 2D
        self.heightmap[*pos] = self.get_max_height(pos)
        self.particles.pop(index)
        return
    
    def get_particle_index(self,pos):
        for i,p in enumerate(self.particles):
            if np.array_equal(p.pos, list(pos)):
                return i
            
    def generate_new_particle(self):
        # generating a random incoming particle with energy distribution of 
        # Maxwell-Boltzmann and a random position with a uniform angular
        # distribution assuming the slit of the effusion cell is really far away
        E = sp.stats.gamma.rvs(1.5, loc=0, scale=10)  #  gamma.pdf(y, a) / scale    with    y = (x - loc) / scale
        pos = list(np.random.randint(0,self.size))
        return Particle(pos, E)     # position in 2D
    
    def calculate_interaction_energy(self,pos):
        E = sum([self.E_bond for site in self.intralayer_interaction_sites(pos) if self.lattice[*site]]\
                +[self.E_bound for site in self.interlayer_interaction_sites(pos) if self.lattice[*site]])
        if pos[2] == 0: # lowest layer -> correct for substrate contribution
            E += sum([self.E_substrate for site in self.interlayer_interaction_sites(np.array([0,0,1])) if site[2]==0])
        return E
    
    def calculate_rate_factors(self, T):
        occ_x,occ_y,occ_z = np.where(self.lattice==True)
        occ_sites = [[x,y,z] for x,y,z in zip(occ_x,occ_y,occ_z)]
        unocc_nextto_occ = [[site for site in self.all_interaction_sites(occ) if (not self.lattice[*site] and site[2]>=0)] for occ in occ_sites]
        origin_list = []
        destination_list = []
        rate_list = []
        for origin,destinations in zip(occ_sites,unocc_nextto_occ):
            E_old = self.calculate_interaction_energy(origin)
            for destination in destinations:
                origin_list += [origin]
                destination_list += [destination]
                E_hop = self.calculate_interaction_energy(destination)
                # correction of self-interaction
                delta = destination-origin
                if delta[2]==0: # intralayer hopping
                    E_hop -= self.E_bond
                else: # interlayer hopping
                    E_hop -= self.E_bound # (substrate double counting not possible)
                rate_list += [(self.particles[self.get_particle_index(origin)].E/T)**2 * np.exp((E_hop-E_old)/T)]
        return origin_list,destination_list,rate_list
    
    def movement_energy_transfer(self,pos):
        nearest_particle_pos_list = [pos] + [site for site in self.all_interaction_sites(pos) if self.lattice[*site]]
        E_mean = np.mean([self.particles[self.get_particle_index(site)].E for site in nearest_particle_pos_list])
        E_sum_abs_dev_from_mean = np.sum([np.abs(E_mean-self.particles[self.get_particle_index(site)].E) for site in nearest_particle_pos_list])
        if E_sum_abs_dev_from_mean==0:
            return # no transfer
        for site in nearest_particle_pos_list:
            self.particles[self.get_particle_index(site)].E += E_mean*(E_mean-self.particles[self.get_particle_index(site)].E)/E_sum_abs_dev_from_mean
        return
    
    def hopping_energy_transfer(self,origin,destination):
        # energy change due to bond structure change:
        #  around origin (exclude already moved hopping particle from update)
        for pos in [site for site in self.intralayer_interaction_sites(origin) if self.lattice[*site] and not np.array_equal(site,destination)]:
            self.particles[self.get_particle_index(pos)].E += self.E_bond
        for pos in [site for site in self.interlayer_interaction_sites(origin) if self.lattice[*site] and not np.array_equal(site,destination)]:
            self.particles[self.get_particle_index(pos)].E += self.E_bound
        #  around destination
        for pos in [site for site in self.intralayer_interaction_sites(destination) if self.lattice[*site]]:
            self.particles[self.get_particle_index(pos)].E -= self.E_bond
        for pos in [site for site in self.interlayer_interaction_sites(destination) if self.lattice[*site]]:
            self.particles[self.get_particle_index(pos)].E -= self.E_bound
        #  moving particle
        delta = destination-origin
        E_change = self.calculate_interaction_energy(origin)-self.calculate_interaction_energy(destination)
        # self-interaction correction
        if delta[2]==0: # intralayer hopping
            E_change -= self.E_bond
        else: # interlayer hopping
            E_change -= self.E_bound
        self.particles[self.get_particle_index(destination)].E += E_change
        # energy change due to hopping movement
        self.movement_energy_transfer(destination)
        return
    
    def hopping_action(self,T):
        origin_list,destination_list,rate_list = self.calculate_rate_factors(T)
        if not origin_list:
            return # no action to be done
        
        idx = np.random.choice(np.arange(len(origin_list)),p=rate_list/np.sum(rate_list))
        origin = origin_list[idx]
        destination = destination_list[idx]
        particle_idx = self.get_particle_index(origin)
        # hopping: update parameters
        self.lattice[*origin] = False
        self.lattice[*destination] = True
        self.particles[particle_idx].pos = destination
        self.heightmap[*origin[:2]] = self.get_max_height(origin[:2])
        self.heightmap[*destination[:2]] = self.get_max_height(destination[:2])
        neighbor_count = len([site for site in self.intralayer_interaction_sites(destination) if self.lattice[*site]]
                             +[site for site in self.interlayer_interaction_sites(destination) if self.lattice[*site]])
        if neighbor_count==0 and destination[2]!=0: # desorption
            self.sub_particle(particle_idx)
        else: # simple hopping
            self.hopping_energy_transfer(origin,destination)
        self.extend_layers() # add new layer if needed
        return
    
    def new_particle_action(self):
        pos = self.add_particle(self.generate_new_particle())
        self.movement_energy_transfer(pos)
        return
        
def animation(particles,heightmap):
    print("Starting animation...")
    fig = plt.figure(figsize=[24,20])
    ax_scatter = fig.add_subplot(2,2,1,projection='3d')
    ax_heightmap = fig.add_subplot(2,2,2)
    ax_particle_per_layer = fig.add_subplot(2,2,3)
    ax_energy_distribution = fig.add_subplot(2,2,4)
    
    scat = ax_scatter.scatter(0, 0, 0, s=90, c=0, cmap='plasma', vmin=0, vmax=10*(E_BOND+E_BOUND+E_SUBSTRATE))
    fig.colorbar(scat,label="Kinetic energy of particle")
    ax_scatter.view_init(elev=12, azim=-40)
    ax_scatter.set_title("3D plot")
    ax_scatter.set_xlabel("x")
    ax_scatter.set_ylabel("y")
    ax_scatter.set_zlabel("z")
    ax_scatter.clear()
    
    hmap = ax_heightmap.scatter(0,0,s=90,c=0,cmap="viridis",vmin=0,vmax=3)
    fig.colorbar(hmap,label="Max. height")
    ax_heightmap.set_title("Heightmap")
    ax_heightmap.set_xlabel("x")
    ax_heightmap.set_ylabel("y")
    ax_heightmap.clear()
    
    ppl = ax_particle_per_layer.hist([0],bins=10,range=(0.5,10.5))
    ax_particle_per_layer.set_title("Number of particles per layer")
    ax_particle_per_layer.set_xlabel("Layer")
    ax_particle_per_layer.set_ylabel("Number of particles")
    ax_particle_per_layer.clear()
    
    en_dist = ax_energy_distribution.hist([0],bins=20,density=True,range=(0,10*(E_BOND+E_BOUND+E_SUBSTRATE)))
    ax_energy_distribution.set_title("Distribution of particle kinetic energy")
    ax_energy_distribution.set_xlabel("Energy")
    ax_energy_distribution.set_ylabel("Relative number of particles")
    ax_energy_distribution.clear()

    def update(frame):
        print("Working on frame {}".format(frame))
        x = [];    y = [];    z = []
        E = []
        for par in particles[frame]:
            x += [par.pos[0]]; y += [par.pos[1]]; z += [par.pos[2]]
            E += [par.E]
        x_hmap = [];    y_hmap = [];    z_hmap = []
        for x_i in range(SIZE[0]):
            for y_i in range(SIZE[1]):
                if heightmap[frame][x_i,y_i] != 0:
                    x_hmap += [x_i];  y_hmap += [y_i]; z_hmap += [heightmap[frame][x_i,y_i]]
        
        ax_scatter.clear()
        ax_scatter.set_title("3D plot")
        ax_scatter.set_xlabel("x")
        ax_scatter.set_ylabel("y")
        ax_scatter.set_zlabel("z")
        scat = ax_scatter.scatter(project_x(x,y,z), project_y(x,y,z), project_z(x,y,z), s=90, c=E, cmap='plasma', vmin=0, vmax=10*(E_BOND+E_BOUND+E_SUBSTRATE))
        ax_scatter.set_xlim(project_x([0],[0],[0])[0], project_x([SIZE[0]],[SIZE[1]],[0])[0])
        ax_scatter.set_ylim(project_y([0],[0],[0])[0], project_y([SIZE[0]],[SIZE[1]],[0])[0])
        ax_scatter.set_zlim(0, 2)
        
        ax_heightmap.clear()
        ax_heightmap.set_title("Heightmap")
        ax_heightmap.set_xlabel("x")
        ax_heightmap.set_ylabel("y")
        hmap = ax_heightmap.scatter(project_x(x_hmap,y_hmap,z_hmap), project_y(x_hmap,y_hmap,z_hmap), s=200, c=z_hmap, cmap="viridis", vmin=0, vmax=3)
        ax_heightmap.set_xlim(project_x([0],[0],[0])[0]-0.5, project_x([SIZE[0]],[SIZE[1]],[0])[0]-0.5)
        ax_heightmap.set_ylim(project_y([0],[0],[0])[0]-0.5, project_y([SIZE[0]],[SIZE[1]],[0])[0]-0.5)
        
        ax_particle_per_layer.clear()
        ax_particle_per_layer.set_title("Number of particles per layer")
        ax_particle_per_layer.set_xlabel("Layer")
        ax_particle_per_layer.set_ylabel("Number of particles")
        ppl = ax_particle_per_layer.hist([z[i]+1 for i in range(len(z))],bins=10,range=(0.5,10.5))
        ax_particle_per_layer.set_ylim(0,SIZE[0]*SIZE[1])
        
        ax_energy_distribution.clear()
        ax_energy_distribution.set_title("Distribution of particle kinetic energy")
        ax_energy_distribution.set_xlabel("Energy")
        ax_energy_distribution.set_ylabel("Relative number of particles")
        en_dist = ax_energy_distribution.hist(E,bins=20,density=True,range=(0,10*(E_BOND+E_BOUND+E_SUBSTRATE)))
        return scat,hmap,ppl,en_dist,
    
    anim = ani.FuncAnimation(fig=fig, func=update, frames=len(particles), interval=300)
    print("Animation completed!")
    anim.save("deposition.gif", writer='pillow')
    print("Animation saved!")
    return

def main():
    lattice = Lattice(SIZE,LATTICE,NEIGHBORING_SITES,E_BOND,E_BOUND,E_SUBSTRATE)

    particles = []; heightmap = []
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

if __name__=="__main__":
    print("Running simulation...")
    main()
    print("Simulation finished!")