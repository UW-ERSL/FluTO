import numpy as np
import numpy.matlib
import time
import torch
from torch_sparse_solve import solve
from FE_template import computeFluidMatrices
  #-----------------------#
def to_torch(x):
  return torch.tensor(x).float()
  #-----------------------#
class FE:
  #-----------------------#
  def __init__(self, mesh, matProp, bc):
    self.mesh = mesh
    self.bc=bc
    self.Amu, self.B, self.matArea, \
      self.AalphaX, self.AalphaY, \
      self.AalphaXY = computeFluidMatrices(\
                mesh.meshSpec['elemSize'][0],\
                mesh.meshSpec['elemSize'][1],\
                matProp['mu'])
    self.force=np.zeros(self.mesh.meshSpec['numDOFs']['total'])
    self.sol=np.zeros((self.mesh.meshSpec['numDOFs']['total']))
    
    V = np.zeros((self.mesh.meshSpec['numDOFs']['total'], self.mesh.meshSpec['numDOFs']['total']));
    V[self.mesh.meshSpec['bcDOFs']['velocity'],self.mesh.meshSpec['bcDOFs']['velocity']] = 1.
    V = torch.tensor(V[np.newaxis])
    indices = torch.nonzero(V).t()
    values = V[indices[0], indices[1], indices[2]] # modify this based on dimensionality    
    penal = 1e15;
    self.fixedBCPenaltyMatrix = \
        penal*torch.sparse_coo_tensor(indices, values, V.size())
    self.bc['force'][self.mesh.meshSpec['bcDOFs']['velocity']] = self.bc['force'][self.mesh.meshSpec['bcDOFs']['velocity']]*1e15
    self.f = torch.tensor(self.bc['force']).unsqueeze(0).unsqueeze(2)
    self.matProp = matProp
    self.objectiveHandle = (self.objective)
  def objective(self, C):
     
    def assembleK(C):
      Aelem = torch.einsum('e,jk->ejk', C['00'], self.AalphaX) +\
              torch.einsum('e,jk->ejk', C['11'], self.AalphaY) +\
              torch.einsum('e,jk->ejk', C['01'], self.AalphaXY) +\
              self.Amu[np.newaxis,:,:]
           
              
      Kelem = Aelem + self.B[np.newaxis,:,:] + self.matArea[np.newaxis,:,:]
  
      Aasm =  torch.sparse_coo_tensor(self.mesh.nodeIdx_A, Aelem[:,0:18,0:18].flatten(), \
                  (1, self.mesh.meshSpec['numDOFs']['velocity'], self.mesh.meshSpec['numDOFs']['velocity']))   
   
      Kasm = torch.sparse_coo_tensor(self.mesh.nodeIdx, Kelem.flatten(), \
                  (1, self.mesh.meshSpec['numDOFs']['total'], self.mesh.meshSpec['numDOFs']['total']))
      Kt = (Kasm + self.fixedBCPenaltyMatrix).coalesce()
  
      return Kt, Aasm
    #     #--------------------------#
  
    def solveFE( Kt):
      sol = solve(Kt, self.f).flatten()    
      U = sol[0:self.mesh.meshSpec['numDOFs']['velocity']];
      uVelocity = U[0::2];  vVelocity = U[1::2];
      Pressure = sol[self.mesh.meshSpec['numDOFs']['velocity']:\
                              self.mesh.meshSpec['numDOFs']['velocity']+self.mesh.meshSpec['numDOFs']['pressure']];
      return sol
    #     #--------------------------#
  
    def computeObjective(sol, Aasm):
      U = sol[0:self.mesh.meshSpec['numDOFs']['velocity']]
      F = torch.einsum('ij,j->i', Aasm.to_dense()[0,:,:], U.float())
      f = 0.5*torch.einsum('i,i->', F, U.float())
      return f
    #     #--------------------------#
  
    def computeElemVelPress( sol):
      UV = sol[0:self.mesh.meshSpec['numDOFs']['velocity']];
      P = sol[self.mesh.meshSpec['numDOFs']['velocity']:\
      self.mesh.meshSpec['numDOFs']['velocity']+ self.mesh.meshSpec['numDOFs']['pressure']];
      UVElem = UV[self.mesh.edofMat['velocity']].\
      reshape((self.mesh.meshSpec['numElems'],self.mesh.numVelocityDOFsPerElem))
      PatDofs = P[self.mesh.edofMat['pressure']].\
      reshape((self.mesh.meshSpec['numElems'],self.mesh.numPressureDOFsPerElem))
      U = torch.mean(UVElem[:,0::2], axis=1)
      V = torch.mean(UVElem[:,1::2], axis=1)
      Pressure = torch.mean(PatDofs, axis=1)
      return U, V, Pressure
    #     #--------------------------#
    Kt, Aasm = assembleK(C)
    sol=solveFE(Kt)
    f = computeObjective(sol, Aasm)
    Uvel, Vvel, Pressure = computeElemVelPress(sol)
    return f, Uvel, Vvel, Pressure
