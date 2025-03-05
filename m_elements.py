from dolfin import *

class Elements():
    def __init__(self, MeshClass):
        self.mesh = MeshClass.mesh
        
        sDG0_elem = FiniteElement("DG", self.mesh.ufl_cell(), 0) 
        self.sDG0 = FunctionSpace(self.mesh, sDG0_elem)

        tDG0_elem = TensorElement("DG", self.mesh.ufl_cell(), 0) 
        self.tDG0 = FunctionSpace(self.mesh, tDG0_elem)

        vCG1_elem = VectorElement("Lagrange", self.mesh.ufl_cell(), 1)  
        self.vCG1 = FunctionSpace(self.mesh, vCG1_elem)

        vCG2_elem = VectorElement("Lagrange", self.mesh.ufl_cell(), 2)
        self.vCG2 = FunctionSpace(self.mesh, vCG2_elem)

        sCG1_elem = FiniteElement("Lagrange", self.mesh.ufl_cell(), 1)
        self.sCG1 = FunctionSpace(self.mesh, sCG1_elem)

        sCG2_elem = FiniteElement("Lagrange", self.mesh.ufl_cell(), 2)
        self.sCG2 = FunctionSpace(self.mesh, sCG2_elem)

        vB_elem  = VectorElement("B", self.mesh.ufl_cell(), self.mesh.topology().dim() + 1)
        self.V = FunctionSpace(self.mesh, MixedElement([sCG1_elem, vCG1_elem, vB_elem]))

        