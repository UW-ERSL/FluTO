import numpy as np

def computeBC(bb, profile, mesh, nodetype = 'velocity'):
  bcNodes = []
  bcValues = []
  for i in range(mesh.nodeXY[nodetype][:,0].shape[0]):
    if(mesh.nodeXY[nodetype][i,0] >= bb['xmin'] and \
       mesh.nodeXY[nodetype][i,0] < bb['xmax'] and \
       mesh.nodeXY[nodetype][i,1] >= bb['ymin'] and \
       mesh.nodeXY[nodetype][i,1] < bb['ymax']):
      
      bcNodes.append(i)
      bcValues.append(profile(mesh.nodeXY[nodetype][i,0], \
                              mesh.nodeXY[nodetype][i,1]))

  return bcNodes, np.array(bcValues)

def getBC(example, mesh):
  '''
  Parameters
  ----------
  example : Example no, choose from defined set of examples: integer
  meshSpec : Dictionary should contain nelx, nely

  Returns
  -------
  boundary condition dictionary. only defined for rect geom here

  '''
  if(example == 1):
    bc = diffuser(mesh)
  elif(example == 2):
    bc = bend(mesh)
  else:
    bc = double_pipe(mesh)
  return bc


def diffuser(mesh):
  bc = {'velocity':{'u':np.zeros((mesh.meshSpec['numNodes']['velocity'])),\
                    'v':np.zeros((mesh.meshSpec['numNodes']['velocity']))}, \
        'pressure':np.zeros((mesh.meshSpec['numNodes']['pressure'])),\
        'force':np.zeros((mesh.meshSpec['numDOFs']['total']))}
  tol=1e-8
  lx, ly = mesh.meshSpec['nelx']*mesh.meshSpec['elemSize'][0],\
            mesh.meshSpec['nely']*mesh.meshSpec['elemSize'][1]
  inletDomain = {'xmin':0., 'xmax':tol, 'ymin':0., 'ymax':ly+ tol}
  outletDomain = {'xmin':lx-tol,'xmax':lx+tol,'ymin':ly/3.-tol, 'ymax':2.*ly/3.+tol}
  topDomain = {'xmin':0.+0.5*mesh.meshSpec['elemSize'][0],\
                'xmax':lx-0.5*mesh.meshSpec['elemSize'][0]+tol, 'ymin':ly, 'ymax':ly+tol}
  bottomDomain = {'xmin':0.+0.5*mesh.meshSpec['elemSize'][0],\
                  'xmax':lx-0.5*mesh.meshSpec['elemSize'][0]+tol, 'ymin':0, 'ymax':tol}
  totaloutletDomain = {'xmin':lx-tol,'xmax':lx+tol,'ymin':0, 'ymax':ly+tol}
  # TODO: there could some tolerancing issue (minor) that needs fixing
  inletVelocityProfile = lambda x,y: 4*y*(ly-y)
  outletVelocityProfile = lambda x,y: (y-ly/3.)*(2.*ly/3.-y)
  topVelocityProfile = lambda x,y: 0
  bottomVelocityProfile = lambda x,y: 0
  totaloutletVelocityProfile = lambda x,y: 0
  
  inletNodes, inletBC = computeBC(inletDomain, inletVelocityProfile,\
                                                  mesh, 'velocity')
  
  inletBC = (inletBC-min(inletBC))/(max(inletBC)-min(inletBC)) # normalize inlet vel
  
  outletNodes, outletBC = computeBC(outletDomain, outletVelocityProfile, \
                                                  mesh, 'velocity')
  outletBC = np.sum(inletBC)*outletBC/np.sum(outletBC)
  
  topNodes, topBC = computeBC(topDomain, topVelocityProfile,\
                                                  mesh, 'velocity')
  bottomNodes, bottomBC = computeBC(bottomDomain, bottomVelocityProfile,\
                                                  mesh, 'velocity')
  totaloutletNodes, totaloutletBC = computeBC(totaloutletDomain, totaloutletVelocityProfile,\
                                                  mesh, 'velocity')
  
  bc['velocity']['u'][inletNodes] = inletBC
  bc['velocity']['u'][topNodes] = topBC
  bc['velocity']['u'][bottomNodes] = bottomBC
  bc['velocity']['u'][totaloutletNodes] = totaloutletBC
  bc['velocity']['u'][outletNodes] = outletBC
  force_nodes=np.array((inletNodes+outletNodes))
  force=(np.hstack((inletBC,outletBC)))
  
  mesh.meshSpec['bcNodes']['velocity'] = (np.array((inletNodes+bottomNodes+topNodes+totaloutletNodes)))
  mesh.meshSpec['bcDOFs']['velocity'] =np.zeros(2*mesh.meshSpec['bcNodes']['velocity'].shape,dtype=np.int)
  (mesh.meshSpec['bcDOFs']['velocity'][0: :2])=(2*(mesh.meshSpec['bcNodes']['velocity'][:]))
  (mesh.meshSpec['bcDOFs']['velocity'][1: :2])=(2*(mesh.meshSpec['bcNodes']['velocity'][:])+1)
  bc['force'][2*force_nodes[:]]=force[:]
  return bc

def bend(mesh):
  bc = {'velocity':{'u':np.zeros((mesh.meshSpec['numNodes']['velocity'])),\
                    'v':np.zeros((mesh.meshSpec['numNodes']['velocity']))}, \
        'pressure':np.zeros((mesh.meshSpec['numNodes']['pressure'])),\
        'force':np.zeros((mesh.meshSpec['numDOFs']['total']))}
  tol=1e-8
  lx, ly = mesh.meshSpec['nelx']*mesh.meshSpec['elemSize'][0],\
            mesh.meshSpec['nely']*mesh.meshSpec['elemSize'][1]       
  totalinletDomain = {'xmin':0-tol,'xmax':0+tol,'ymin':0-tol, 'ymax':ly+tol}
  inletDomain = {'xmin':0.-tol, 'xmax':tol, 'ymin':(ly*0.8)-tol, 'ymax':(ly*1)+tol}
  outletDomain = {'xmin':0-tol,'xmax':0+tol, 'ymin':(ly*0)-tol, 'ymax':(ly*0.2)+tol}
  bottomDomain = {'xmin':0.+0.5*mesh.meshSpec['elemSize'][0],\
                'xmax':lx-0.5*mesh.meshSpec['elemSize'][0]+tol, 'ymin':0, 'ymax':tol}
  rightDomain = {'xmin':lx-tol,'xmax':lx+tol,'ymin':0, 'ymax':ly+tol}
  # print('jjj',inletDomain2,outletDomain1,outletDomain2)
  topDomain = {'xmin':0.+0.5*mesh.meshSpec['elemSize'][0],\
                'xmax':lx-0.5*mesh.meshSpec['elemSize'][0]+tol, 'ymin':ly, 'ymax':ly+tol}
  # TODO: there could some tolerancing issue (minor) that needs fixing
  totalinletVelocityProfile = lambda x,y: 0
  inletVelocityProfile = lambda x,y: (y-ly*0.8)*(ly*1.-y)
  outletVelocityProfile = lambda x,y:(y-ly*0.)*(ly*0.2-y)
  topVelocityProfile = lambda x,y: 0
  rightVelocityProfile = lambda x,y: 0
  totaloutletVelocityProfile = lambda x,y: 0
  bottomVelocityProfile = lambda x,y: 0
  
  inletNodes, inletBC = computeBC(inletDomain, inletVelocityProfile,\
                                                  mesh, 'velocity')


  inletBC = (inletBC-min(inletBC))/(max(inletBC)-min(inletBC)) # normalize inlet vel  
  totalinletNodes, totalinletBC = computeBC(totalinletDomain, totalinletVelocityProfile,\
                                                  mesh, 'velocity')
  outletNodes, outletBC = computeBC(outletDomain, outletVelocityProfile, \
                                                  mesh, 'velocity')

  outletBC = np.sum(inletBC)*outletBC/np.sum(outletBC)
  
  topNodes, topBC = computeBC(topDomain, topVelocityProfile,\
                                                  mesh, 'velocity')
  rightNodes, rightBC = computeBC(rightDomain, rightVelocityProfile,\
                                                  mesh, 'velocity')
  bottomNodes, bottomBC = computeBC(bottomDomain, bottomVelocityProfile,\
                                                  mesh, 'velocity')
  
  bc['velocity']['u'][totalinletNodes] = totalinletBC
  bc['velocity']['u'][inletNodes] = inletBC
  bc['velocity']['u'][topNodes] = topBC
  bc['velocity']['u'][rightNodes] = rightBC

  bc['velocity']['v'][outletNodes] = outletBC
  bc['velocity']['u'][bottomNodes] = bottomBC
  force_inlet_nodes=np.array((inletNodes))
  force_outlet_nodes=np.array((outletNodes))
  force_inlet=(np.hstack((inletBC)))
  force_outlet=(np.hstack((outletBC)))
  
  mesh.meshSpec['bcNodes']['velocity'] = (np.array((totalinletNodes+rightNodes+topNodes+bottomNodes)))
  mesh.meshSpec['bcDOFs']['velocity'] =np.zeros(2*mesh.meshSpec['bcNodes']['velocity'].shape,dtype=np.int)
  (mesh.meshSpec['bcDOFs']['velocity'][0: :2])=(2*(mesh.meshSpec['bcNodes']['velocity'][:]))
  (mesh.meshSpec['bcDOFs']['velocity'][1: :2])=(2*(mesh.meshSpec['bcNodes']['velocity'][:])+1)
  bc['force'][2*force_inlet_nodes[:]]=force_inlet[:]
  bc['force'][2*force_outlet_nodes[:]]=-force_outlet[:]

  return bc



def double_pipe(mesh):
  bc = {'velocity':{'u':np.zeros((mesh.meshSpec['numNodes']['velocity'])),\
                    'v':np.zeros((mesh.meshSpec['numNodes']['velocity']))}, \
        'pressure':np.zeros((mesh.meshSpec['numNodes']['pressure'])),\
        'force':np.zeros((mesh.meshSpec['numDOFs']['total']))}
  tol=1e-8
  lx, ly = mesh.meshSpec['nelx']*mesh.meshSpec['elemSize'][0],\
               mesh.meshSpec['nely']*mesh.meshSpec['elemSize'][1]       
  totalinletDomain = {'xmin':0-tol,'xmax':0+tol,'ymin':0-tol, 'ymax':ly+tol}
  inletDomain1 = {'xmin':0., 'xmax':tol, 'ymin':(ly/6.)-tol, 'ymax':(ly/3)+tol}
  inletDomain2 = {'xmin':0., 'xmax':tol, 'ymin':(2*ly/3)-tol, 'ymax':(5*ly/6)+tol}
  outletDomain1 = {'xmin':lx-tol,'xmax':lx+tol,'ymin':(ly/6.)-tol, 'ymax':(ly/3)+tol}
  outletDomain2 = {'xmin':lx-tol,'xmax':lx+tol,'ymin':(2*ly/3.)-tol, 'ymax':(5*ly/6)+tol}
  topDomain = {'xmin':0.+0.5*mesh.meshSpec['elemSize'][0],\
               'xmax':lx-0.5*mesh.meshSpec['elemSize'][0]+tol, 'ymin':ly, 'ymax':ly+tol}
  bottomDomain = {'xmin':0.+0.5*mesh.meshSpec['elemSize'][0],\
                  'xmax':lx-0.5*mesh.meshSpec['elemSize'][0]+tol, 'ymin':0, 'ymax':tol}
  totaloutletDomain = {'xmin':lx-tol,'xmax':lx+tol,'ymin':0, 'ymax':ly+tol}
  # TODO: there could some tolerancing issue (minor) that needs fixing
  totalinletVelocityProfile = lambda x,y: 0
  inletVelocityProfile1 = lambda x,y: (y-ly/6.)*(ly/3.-y)
  inletVelocityProfile2 = lambda x,y: (y-2*ly/3.)*(5*ly/6.-y)
  outletVelocityProfile1 = lambda x,y: (y-ly/6.)*(ly/3.-y)
  outletVelocityProfile2 = lambda x,y: (y-2*ly/3.)*(5*ly/6.-y)
  topVelocityProfile = lambda x,y: 0
  bottomVelocityProfile = lambda x,y: 0
  totaloutletVelocityProfile = lambda x,y: 0
  
  inletNodes1, inletBC1 = computeBC(inletDomain1, inletVelocityProfile1,\
                                                  mesh, 'velocity')
  # print(inletNodes1, inletBC1)
  inletNodes2, inletBC2 = computeBC(inletDomain2, inletVelocityProfile2,\
                                                  mesh, 'velocity')
  inletBC1 = 1.3*(inletBC1-min(inletBC1))/(max(inletBC1)-min(inletBC1)) # normalize inlet vel
  inletBC2 = 1.3*(inletBC2-min(inletBC2))/(max(inletBC2)-min(inletBC2)) # normalize inlet vel
  totalinletNodes, totalinletBC = computeBC(totalinletDomain, totalinletVelocityProfile,\
                                                  mesh, 'velocity')
  outletNodes1, outletBC1 = computeBC(outletDomain1, outletVelocityProfile1, \
                                                  mesh, 'velocity')
  outletNodes2, outletBC2 = computeBC(outletDomain2, outletVelocityProfile2, \
                                                  mesh, 'velocity')
  outletBC1 = np.sum(inletBC1)*outletBC1/np.sum(outletBC1)
  outletBC2 = np.sum(inletBC2)*outletBC2/np.sum(outletBC2)
  
  topNodes, topBC = computeBC(topDomain, topVelocityProfile,\
                                                  mesh, 'velocity')
  bottomNodes, bottomBC = computeBC(bottomDomain, bottomVelocityProfile,\
                                                  mesh, 'velocity')
  totaloutletNodes, totaloutletBC = computeBC(totaloutletDomain, totaloutletVelocityProfile,\
                                                  mesh, 'velocity')
  
  bc['velocity']['u'][totalinletNodes] = totalinletBC
  bc['velocity']['u'][inletNodes1] = inletBC1
  bc['velocity']['u'][inletNodes2] = inletBC2
  bc['velocity']['u'][topNodes] = topBC
  bc['velocity']['u'][bottomNodes] = bottomBC
  bc['velocity']['u'][totaloutletNodes] = totaloutletBC
  bc['velocity']['u'][outletNodes1] = outletBC1
  bc['velocity']['u'][outletNodes2] = outletBC2
  force_nodes=np.array((inletNodes1+inletNodes2+outletNodes1+outletNodes2))
  force=(np.hstack((inletBC1,inletBC2,outletBC1,outletBC2)))
  
  mesh.meshSpec['bcNodes']['velocity'] = (np.array((totalinletNodes+bottomNodes+topNodes+totaloutletNodes)))
  mesh.meshSpec['bcDOFs']['velocity'] =np.zeros(2*mesh.meshSpec['bcNodes']['velocity'].shape,dtype=np.int)
  (mesh.meshSpec['bcDOFs']['velocity'][0: :2])=(2*(mesh.meshSpec['bcNodes']['velocity'][:]))
  (mesh.meshSpec['bcDOFs']['velocity'][1: :2])=(2*(mesh.meshSpec['bcNodes']['velocity'][:])+1)
  bc['force'][2*force_nodes[:]]=force[:]  
  return bc
    
    
    