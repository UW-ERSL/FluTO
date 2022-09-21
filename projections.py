import torch
import numpy as np

#-----------------------------#
def to_np(x):
  return x.detach().cpu().numpy()
#-----------------------------#

#-------FOURIER LENGTH SCALE-----------#
def computeFourierMap(mesh, fourierMap):
  # compute the map
  coordnMapSize = (mesh.meshSpec['ndim'], fourierMap['numTerms'])
  freqSign = np.random.choice([-1.,1.], coordnMapSize)
  stdUniform = np.random.uniform(0.,1., coordnMapSize) 
  wmin = 1./(2*fourierMap['maxRadius']*mesh.meshSpec['elemSize'][0])
  wmax = 1./(2*fourierMap['minRadius']*mesh.meshSpec['elemSize'][0]) # w~1/R
  wu = wmin +  (wmax - wmin)*stdUniform
  coordnMap = np.einsum('ij,ij->ij', freqSign, wu)
  coordnMap = torch.tensor(coordnMap).float()
  return coordnMap
#-----------------#
def applyFourierMap(xy, fourierMap):
  if(fourierMap['isOn']):
    c = torch.cos(2*np.pi*torch.einsum('ij,jk->ik', xy, fourierMap['map']))
    s = torch.sin(2*np.pi*torch.einsum('ij,jk->ik', xy, fourierMap['map']))
    xy = torch.cat((c,s), axis = 1);
  return xy

#-------DENSITY PROJECTION-----------#

def applyDensityProjection(x, densityProj):
  if(densityProj['isOn']):
    b = densityProj['sharpness']
    nmr = np.tanh(0.5*b) + torch.tanh(b*(x-0.5))
    x = 0.5*nmr/np.tanh(0.5*b)
  return x

#-------SYMMETRY-----------#
# def applySymmetry(x, symMap):
#   if(symMap['YAxis']['isOn']):
#     xv = index_update( x[:,0], index[:], symMap['YAxis']['midPt'] \
#                           + torch.abs(x[:,0] - symMap['YAxis']['midPt']) )
#   else:
#     xv = x[:, 0]
#   if(symMap['XAxis']['isOn']):
#     yv = index_update( x[:,1], index[:], symMap['XAxis']['midPt'] \
#                           + torch.abs(x[:,1] - symMap['XAxis']['midPt']) )
#   else:
#     yv = x[:, 1]
#   x = torch.transpose(torch.stack((xv,yv)),0,1);
#   return x
#-------REFLECTION-----------#
def applyReflection(x, symMap):
  signs = {}
  if(symMap['YAxis']['isOn']):
    signs['Y'] =  -torch.sign(x[:,0] - symMap['YAxis']['midPt'] + 1e-6)
    xv =( symMap['YAxis']['midPt'] + torch.abs( x[:,0] - symMap['YAxis']['midPt']));
  else:
    signs['Y'] = torch.ones((x.shape[0]))
    xv = x[:, 0]
  if(symMap['XAxis']['isOn']):
    signs['X'] = -torch.sign(x[:,1] - symMap['XAxis']['midPt'] + 1e-6)
    yv = (symMap['XAxis']['midPt'] + torch.abs( x[:,1] - symMap['XAxis']['midPt'])) ;
  else:
    signs['X'] = torch.ones((x.shape[0]))
    yv = x[:, 1]
  
  x = torch.transpose(torch.stack((xv,yv)), 0, 1);
  return x, signs
#--------------------------#

