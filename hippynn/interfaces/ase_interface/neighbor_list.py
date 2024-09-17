from dataclasses import dataclass
import numpy as np
from typing import Union, Literal

class NeighborList:
    def __init__(self, coords: np.ndarray, cell: np.ndarray, pbc: np.ndarray, cutoff: float, 
                 buffer: float = 0, sorted=False, backend: Literal['ASE', 'matscipy'] = "ASE"):
        self.coords = coords.copy()
        self.pbc = pbc.copy()
        self.cell = cell.copy()
        self.cutoff = cutoff
        self.buffer = buffer
        self.sorted = sorted
        self.backend = backend
        self.sender = None
        self.receiver = None
        self.shift = None


    @classmethod
    def from_atoms(cls, atoms, cutoff, buffer=0,
                   sorted = False,
                   backend: Literal['ASE', 'matscipy'] = "ASE"):
        #cell = atoms.cell.array if any(atoms.pbc) else None
        return cls(atoms.get_positions().copy(),
                   atoms.cell.array.copy(),
                   atoms.pbc.copy(),
                   cutoff,
                   buffer,
                   sorted,
                   backend)

    def __post_init__(self):
       self.update(self.coords)

    def update(self, coords, cell, pbc):
        # if any atom moved more than "buffer", update the nbr list
        if (np.any(np.sum((self.coords - coords)**2, axis=0) > self.buffer**2)
            or (np.any(cell != self.cell))
            or (np.any(pbc != self.pbc))):
            if self.backend == "ASE":
                [self.sender, 
                self.receiver,
                self.shift] = generate_neighbor_list_ASE(self.coords, 
                                                        self.cell, 
                                                        self.cutoff + self.buffer)
            if self.backend == "matscipy":
                [self.sender, 
                self.receiver,
                self.shift] = generate_neighbor_list_matscipy(self.coords, 
                                                              self.cell, 
                                                              self.cutoff + self.buffer)
            
        if len(coords) > 0 and self.sorted:
            mask = np.argsort(self.sender * len(self.sender) +
                              self.receiver)
            self.sender = self.sender[mask]
            self.receiver = self.receiver[mask]
            self.shift = self.shift[mask]

        self.coords = coords.copy()
        self.cell = cell.copy()
        self.pbc = pbc.copy()
def generate_neighbor_list_matscipy(coords: np.ndarray, cell: Union[np.ndarray, None], 
                                    pbc: Union[np.ndarray, None],
                                    cutoff: float):
    from matscipy.neighbours import neighbour_list
    return neighbour_list(quantities="ijS",
                           cutoff=cutoff,
                           positions=coords,
                           cell=cell,
                           pbc=pbc)

def generate_neighbor_list_ASE(coords: np.ndarray, cell: Union[np.ndarray, None], 
                               pbc: Union[np.ndarray, None],
                               cutoff: float):
    from ase.neighborlist import primitive_neighbor_list
    
    return primitive_neighbor_list(
            'ijS', pbc, cell, coords, cutoff,
            self_interaction=False)