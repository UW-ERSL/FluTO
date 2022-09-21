import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.collections

#--------------------------#
class Q2Q1Mesher: # only for grid mesh
    

###########################################################################################################
  def __init__(self, meshSpec):
    self.meshSpec = meshSpec
    
    
    self.meshSpec['numNodes'] = {}
    self.meshSpec['numNodes']['velocity'] = (2*meshSpec['nelx']+1)*(2*meshSpec['nely']+1)
    self.meshSpec['numNodes']['pressure'] = (meshSpec['nelx']+1)*(meshSpec['nely']+1)
            
    self.meshSpec['numDOFs'] = {}
    self.meshSpec['numDOFs']['velocity'] = 2*self.meshSpec['numNodes']['velocity'] 
    self.meshSpec['numDOFs']['pressure'] = self.meshSpec['numNodes']['pressure']
    self.meshSpec['numDOFs']['total'] = self.meshSpec['numDOFs']['velocity']+\
                                self.meshSpec['numNodes']['pressure']+1  
    bb = {}
    bb['xmin'], bb['xmax'], bb['ymin'], bb['ymax'] = 0.,\
      meshSpec['nelx']*meshSpec['elemSize'][0], 0.,\
      meshSpec['nely']*meshSpec['elemSize'][1]
    self.meshSpec['bb'] = bb
    self.meshSpec['bcNodes'] = {}
    self.meshSpec['bcDOFs'] = {}
    self.meshSpec['bcNodes']['velocity'] = []
    self.meshSpec['bcDOFs']['velocity'] = []
    
    self.nodes = {'velocity':np.arange(self.meshSpec['numNodes']['velocity']),\
                  'pressure':np.arange(self.meshSpec['numNodes']['pressure'])}
      
    self.meshSpec['dofs'] = {'velocity':np.arange(self.meshSpec['numDOFs']['velocity']),\
                              'pressure':np.arange(self.meshSpec['numDOFs']['pressure'])}
    ######################################################################
    self.meshSpec['globaldofs']=[]
    self.meshSpec['freedofs']=[]
    self.globaldofs=np.arange(self.meshSpec['numDOFs']['velocity'])
    arr2=(np.arange(self.meshSpec['numDOFs']['velocity'],\
                    self.meshSpec['numDOFs']['velocity']+1+ self.meshSpec['numDOFs']['pressure']))
    self.globaldofs=np.concatenate((self.globaldofs, arr2), axis=0) 
    self.meshSpec['globaldofs']=self.globaldofs.tolist()

    ##################################################################
    self.nodeXY = {}
    x = np.linspace(0,meshSpec['elemSize'][0]*meshSpec['nelx'], 2*meshSpec['nelx']+1)
    y = np.linspace(0,meshSpec['elemSize'][1]*meshSpec['nely'], 2*meshSpec['nely']+1)
    Y,X = np.meshgrid(y,x)
    self.nodeXY['velocity'] = np.vstack((X.reshape(-1), Y.reshape(-1))).T

    
    x = np.linspace(0,meshSpec['elemSize'][0]*meshSpec['nelx'], meshSpec['nelx']+1)
    y = np.linspace(0,meshSpec['elemSize'][1]*meshSpec['nely'], meshSpec['nely']+1)
    Y,X = np.meshgrid(y,x)
    self.nodeXY['pressure'] = np.vstack((X.reshape(-1), Y.reshape(-1))).T
    
    self.elemNodes = {}
    self.edofMat = {}
    self.numVelocityNodesPerElem = 9
    self.numVelocityDOFsPerElem = 18
    self.numPressureNodesPerElem = 4 # same as numNodesPerElem
    self.numPressureDOFsPerElem = 4 # same as numNodesPerElem
    self.totalDOFsPerElem = self.numVelocityDOFsPerElem + self.numPressureDOFsPerElem +1
    self.elemNodes['velocity'] = np.zeros((self.meshSpec['numElems'], \
                                           self.numVelocityNodesPerElem))
    self.elemNodes['pressure'] = np.zeros((self.meshSpec['numElems'], \
                                           self.numPressureNodesPerElem))
    self.edofMat['velocity'] = np.zeros((self.meshSpec['numElems'], \
                                            self.numVelocityDOFsPerElem))

    for elx in range(self.meshSpec['nelx']):
      for ely in range(self.meshSpec['nely']):
        el = ely+elx*self.meshSpec['nely']
        n1 = (2*self.meshSpec['nely']+1)*2*elx + 2*ely
        n2 = (2*self.meshSpec['nely']+1)*(2*elx+1) + 2*ely
        n3 = (2*self.meshSpec['nely']+1)*(2*elx+2) + 2*ely
        self.elemNodes['velocity'][el,:] = np.array([n1, n3, n3+2, n1+2, \
                                                   n2, n3+1, n2+2, n1+1, n2+1])

        self.edofMat['velocity'][el,:] = np.array([2*n1, 2*n1+1, 2*(n3), 2*(n3)+1,\
                                            2*(n3+2), 2*(n3+2)+1, 2*(n1+2), 2*(n1+2)+1,\
                                            2*(n2), 2*(n2)+1, 2*(n3+1), 2*(n3+1)+1,\
                                            2*(n2+2), 2*(n2+2)+1, 2*(n1+1), 2*(n1+1)+1,\
                                            2*(n2+1), 2*(n2+1)+1])

    self.elemNodes['velocity'] = self.elemNodes['velocity'].astype(int)
    self.edofMat['velocity'] = self.edofMat['velocity'].astype(int)

    self.edofMat['pressure'] = np.zeros((self.meshSpec['numElems'], \
                                            self.numPressureDOFsPerElem))
    ctr = 0
    for elx in range(self.meshSpec['nelx']):
      for ely in range(self.meshSpec['nely']):
        el = ely + elx*self.meshSpec['nely']
        n1 = (self.meshSpec['nely']+1)*elx+ely
        n2 = (self.meshSpec['nely']+1)*(elx+1)+ely
        self.elemNodes['pressure'][el,:] = np.array([n1, n2, n2+1, n1+1])
        self.edofMat['pressure'][el,:] = ctr + np.array([n1, n2, n2+1, n1+1])

    self.elemNodes['pressure'] = self.elemNodes['pressure'].astype(int)
    self.edofMat['pressure'] = self.edofMat['pressure'].astype(int)
    
    self.edofMat['global'] = np.zeros((self.meshSpec['numElems'], \
                            self.numVelocityDOFsPerElem+self.numPressureDOFsPerElem))
    onematrix = np.ones((self.meshSpec['numElems'],1))
    self.edofMat['global'] = np.hstack((self.edofMat['velocity'], \
                    self.meshSpec['numDOFs']['velocity']+self.edofMat['pressure'],\
                    self.meshSpec['numDOFs']['velocity']+\
                    self.meshSpec['numDOFs']['pressure']+onematrix-1))

    self.iK= np.kron(self.edofMat['global'],\
                        np.ones((self.totalDOFsPerElem ,1))).flatten().astype(int)
    self.jK = np.kron(self.edofMat['global'],\
                        np.ones((1,self.totalDOFsPerElem))).flatten().astype(int)
    bK = tuple(np.zeros((len(self.iK))).astype(int)) #batch values
    #self.nodeIdx = jax.ops.index[self.iK,self.jK]
    self.nodeIdx = [bK,self.iK,self.jK] # TODO: Getback batch dim!!!!
    
    self.iA= np.kron(self.edofMat['velocity'],\
                        np.ones((self.numVelocityDOFsPerElem ,1))).flatten().astype(int)
    self.jA = np.kron(self.edofMat['velocity'],\
                        np.ones((1,self.numVelocityDOFsPerElem))).flatten().astype(int)
    bA = tuple(np.zeros((len(self.iA))).astype(int)) #batch values
    self.nodeIdx_A = [bA,self.iA,self.jA] # TODO: Getback batch dim!!!!
    #self.nodeIdx_A = jax.ops.index[self.iA,self.jA]
    self.fig, self.ax = plt.subplots()
  #--------------------------#
  def plotMesh(self, field, bc, scatterVelocityNodes, scatterPressureNodes, \
               annotateVelocityNodes, annotatePressureNodes, plotBC):
    # Field is elemental field
    # return the figure axis
    
    fig, ax = plt.subplots()
    
    # scatter the pressure and velocity nodes
    if(scatterVelocityNodes):
      plt.scatter(self.nodeXY['velocity'][:,0], self.nodeXY['velocity'][:,1],\
                  marker = 'o', facecolors='none', edgecolors='b')
    if(scatterPressureNodes):
      plt.scatter(self.nodeXY['pressure'][:,0], self.nodeXY['pressure'][:,1],\
                  marker = 's', facecolors='none', edgecolors='r')
    
    # plot grid and field
    def quatplot(y,z, quatrangles, values, ax=None, **kwargs):
      if not ax: ax=plt.gca()
      yz = np.c_[y,z]
      verts= yz[quatrangles]
      pc = matplotlib.collections.PolyCollection(verts, **kwargs)
      pc.set_array(values)
      ax.add_collection(pc)
      ax.autoscale()
      return pc
    field = None # HARD CODED TEST
    pc = quatplot(self.nodeXY['pressure'][:,0], self.nodeXY['pressure'][:,1],\
                  np.asarray(self.elemNodes['pressure']), field, ax=ax, 
                  edgecolor="crimson", cmap="gray", alpha = 0.2)#, facecolor = 'none'
    plt.colorbar(pc)
    
    # annotate the pressure nodes
    if(annotatePressureNodes):
      for i in range(self.nodeXY['pressure'].shape[0]):
        ax.annotate(i, (self.nodeXY['pressure'][i,0], \
                        self.nodeXY['pressure'][i,1]), size = 13, c = 'r')
    
    # annotate the velocity nodes
    if(annotateVelocityNodes):
      for i in range(self.nodeXY['velocity'].shape[0]):
        ax.annotate(i, (self.nodeXY['velocity'][i,0], \
                        self.nodeXY['velocity'][i,1]), size = 13, c = 'b')
    
    if(plotBC):
      velMag = np.sqrt(bc['velocity']['u']**2 + bc['velocity']['v']**2)
      maxVelMag = np.max(velMag)
      uvel = bc['velocity']['u']/(1e-3+maxVelMag)
      vvel = bc['velocity']['v']/(1e-3+maxVelMag)
      # print(uvel)
      plt.quiver(self.nodeXY['velocity'][:,0],\
                 self.nodeXY['velocity'][:,1],\
          np.ma.array(uvel, mask = velMag <= 0.),\
          np.ma.array(vvel, mask = velMag <= 0.),\
          color = 'purple',scale = 3) # TODO: scaling?
      
    plt.axis('equal')
    plt.show()
    
    return fig, ax
  #--------------------------#
  def generatePoints(self,  resolution = 1): # generate points in elements
      ctr = 0;
      xy = np.zeros((resolution*self.meshSpec['nelx']*resolution*self.meshSpec['nely'],2));

      for i in range(resolution*self.meshSpec['nelx']):
          for j in range(resolution*self.meshSpec['nely']):
              xy[ctr,0] = self.meshSpec['elemSize'][0]*(i + 0.5)/resolution
              xy[ctr,1] = self.meshSpec['elemSize'][1]*(j + 0.5)/resolution
              ctr += 1;

      return xy
  #-----------------------#
  def plotField(self, field, titleStr, res = 1):
    plt.ion(); plt.clf()
    plt.imshow(field.reshape((res*self.meshSpec['nelx'], res*self.meshSpec['nely'])).T, cmap='gray',\
                interpolation='none', origin = 'lower')
    plt.axis('Equal')
    plt.colorbar()
    plt.title(titleStr)
    plt.grid(False)
    self.fig.canvas.draw()
    plt.pause(0.01)

