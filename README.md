# Deposition_simulation

Részecske hozzáadás:
  - külső függvény
  - energia generálás eloszlás (Maxwell)
  - random rácsbeli x-y hely
  - output részecske változóban

Classes:
  - Chamber: overall parameters
  - Lattice:
    - lattice bool(occupied, unoccupied)
    - particle list
    - adjecency matrix (getting the feeling that it does not give enough information to be worth having)
  - Particle:
    - position
    - energy
    - (maybe) type
    - incoming energy from maxwell-boltzmann distribution: Lambert W distr
    - layer-layer bonding energy different from in-layer bonding energy

Ráta faktorok:
  - szomszédok befolyása
  - globál T
  - kinetikus energia
  - w = E_kin/T * exp(beta*Delta_n), if goes vacuum then a layer bonding contribution gets added 

Termalizáció:
  - a szomszédok energiájával arányosan osztjuk szét a hopping utáni energia differneciát
  - termál fürdő: konstans kicsi kiegyenlítés T-vel

Rács 3D elrendezése - primitive cubic 





