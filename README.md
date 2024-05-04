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

Ráta faktorok:
  - szomszédok befolyása
  - globál T
  - kinetikus energia
  - w = E_kin * exp(beta*Delta_n)

Termalizáció:
  - a szomszédok energiájával arányosan osztjuk szét a hopping utáni energia differneciát
  - termál fürdő: konstans kicsi kiegyenlítés T-vel

Rács 3D elrendezése ???

hőmérséklet rendeződés hopping tagja ???

Ráta faktorra kinetikust jobban modellezni ???





