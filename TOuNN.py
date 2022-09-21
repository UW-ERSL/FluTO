import numpy as np
import torch
import torch.optim as optim
import matplotlib.pyplot as plt
from matplotlib import colors
from network import TopNet
from projections import applyFourierMap, \
    applyReflection, applyDensityProjection

from materialCoeffs import microStrs
import pickle
import scipy.io as sio
from scipy.ndimage import rotate
#-----------------------------#
def to_torch(x):
  return torch.tensor(x).float()
#-----------------------------#
def to_np(x):
  return x.detach().cpu().numpy()
#-----------------------------#

class TOuNN:
  def __init__(self, mesh, material, FE, nnSettings, fourierMap, symMap, bc, optimizationMethod, overrideGPU= True):
    
    self.device = self.setDevice(overrideGPU);
    self.FE = FE(mesh, material, bc)
    xy = self.FE.mesh.generatePoints()
    self.xy = torch.tensor(xy, requires_grad = True).\
                                    float().view(-1,2).to(self.device);
    self.fourierMap = fourierMap
    if(fourierMap['isOn']):
      nnSettings['inputDim'] = 2*fourierMap['numTerms']
    else:
      nnSettings['inputDim'] = self.FE.mesh.meshSpec['elemSize']
    self.topNet = TopNet(nnSettings).to(self.device)
    # # densProj to generate 0-1 designs
    self.densityProj = {'isOn':False, 'sharpness':5.}
    self.symMap = symMap
    self.mstrData = microStrs
    self.optimizationMethod = optimizationMethod
  #-----------------------------#
 
  def setDevice(self, overrideGPU):
    if(torch.cuda.is_available() and (overrideGPU == False) ):
      device = torch.device("cuda:0");
      print("GPU enabled")
    else:
      device = torch.device("cpu")
      print("Running on CPU")
    return device;
  #-----------------------------#
  def rotateCMatrix(self, C, theta):
    c = torch.cos(theta)
    s = torch.sin(theta)
    
    Cprime = {}
    Cprime['00'] = (1/C['00'])*c**2 + (1/C['11'])*s**2   #Taking the inverse first then doing rotation transformation
    Cprime['11'] = (1/C['11'])*c**2 + (1/C['00'])*s**2   #As inv(rot*C*rot.T) = rot*inv(C)*rot.T
    Cprime['01'] = (1/C['00'])*c*s - (1/C['11'])*c*s
    Cprime['area'] = C['area']
    Cprime['perim'] = C['perim']
    return Cprime
  # #-----------------------------#
  def optimizeDesign(self, optParams, savedNet):
    convgHistory = {'epoch':[], 'af':[], 'J':[]}
    xyR, signsReflection = applyReflection(self.xy, self.symMap)
    if(self.fourierMap['isOn']):
      self.xyF = applyFourierMap(xyR, self.fourierMap)
    else:
      self.xyF = xyR
    penal = 1
    # C Matrix
    def getCfromCPolynomial(sizePow, mstr):
      C = {}
      for c in ['00', '11', 'perim']: # init
        C[c] = torch.zeros((self.FE.mesh.meshSpec['numElems']))
      for c in ['00', '11', 'perim']: # c = Q0*a^0 + Q1*a^1 + ... + Qn*a^n
        for pw in range(mstr['order']+1):
            C[c] = C[c]+(mstr[c][str(pw)]*sizePow[str(pw)])
      return C
    

    def getCMatrix(mstrType, size, theta):
      sizePow = {} # compute the powers once to avoid repeated calc
      for pw in range(self.mstrData['squircle']['order']+1):# TODO: use the max order of all mstrs
        sizePow[str(pw)] = size**pw
      C = {}
      for c in ['00', '11', 'area', 'perim']:
        C[c] = torch.zeros((self.FE.mesh.meshSpec['numElems']))
      
      for mstrCtr, mstr in enumerate(self.mstrData):
        Cmstr = getCfromCPolynomial(sizePow, self.mstrData[mstr])
        Cmstr['area'] = 1. - self.mstrData[mstr]['sizeCoeff']*sizePow[str(2)]
        mstrPenal =  mstrType[:,mstrCtr]**penal # now select mstrs
        mstrWithoutPenal = mstrType[:,mstrCtr]
        for c in ['00', '11']:
          Cmstr[c] = (torch.maximum(to_torch(0), Cmstr[c])+1e-4)
          C[c] = C[c] + (torch.einsum('i,i->i', mstrPenal, Cmstr[c]))
          # C[c] = Cmstr[c]
        for c in ['perim', 'area']:
          Cmstr[c] = (torch.maximum(to_torch(0), Cmstr[c])+1e-4)
          C[c] = C[c] + (torch.einsum('i,i->i', mstrWithoutPenal, Cmstr[c]))
      Cprime = self.rotateCMatrix(C, theta)
      Cprime['perim'] *= self.FE.mesh.meshSpec['elemSize'][0] # perim generated assuming elemsize of 1... scale by element size
      return Cprime
    self.getCMatrix = getCMatrix
    #--------------------------#
    def saveTopologyNetwork(self, fileName):
      torch.save(self.topNet.state_dict(), fileName)
    #--------------------------#
    def loadTopologyNetwork(self, fileName):
      self.topNet.load_state_dict(torch.load(fileName))
    #--------------------------#
    # optimizer
  
    if(self.optimizationMethod == 'adam'):
        optParams['epochSnapShot'] = 20
        muVar = 0; # lagrange multiplier   
        self.optimizer = optim.Adam(self.topNet.parameters(), amsgrad=True,lr=optParams['learningRate']);  
    else:
        optParams['lossMethod']['alpha0'] = 10 # start
        optParams['lossMethod']['delAlpha'] = 0.15
        optParams['epochSnapShot'] = 1
        muVar = 0;
        self.optimizer = optim.LBFGS(self.topNet.parameters(),
                                     line_search_fn = 'strong_wolfe')
 
    
    #--------------------------#
    def computeLoss(objective, constraints):
      if(optParams['lossMethod']['type'] == 'penalty'):
        alpha = min(100.,optParams['lossMethod']['alpha0'] + \
                self.epoch*optParams['lossMethod']['delAlpha']) # penalty method
        loss = objective
        for c in constraints:
          loss = loss + alpha*c**2
      if(optParams['lossMethod']['type'] == 'logBarrier'):
        t = optParams['lossMethod']['t0']* \
                          optParams['lossMethod']['mu']**self.epoch
        loss = objective
        for c in constraints:
          if(c < (-1/t**2)):
            loss = loss - torch.log(-c)/t
          else:
            loss = loss + t*c - torch.log(1/t**2)/t + 1./t
      if(optParams['lossMethod']['type'] == 'augLag'):
        
        alpha = min(100,  optParams['lossMethod']['alpha0'] + \
                    optParams['lossMethod']['delAlpha']); # Augumented Lagrangian
        # for c in constraints:       
        loss = objective
        for c in constraints:
          loss = loss + alpha*c**2 + muVar*c
      return loss, alpha
    #--------------------------#
    def closure():
      convgHistory = {'epoch':[], 'af':[], 'J':[]}
      self.optimizer.zero_grad()
      self.mstrType, self.size, self.theta = self.topNet(self.xyF)
      # self.size = applyDensityProjection(self.size, self.densityProj)
      self.theta = torch.einsum('i,i->i', self.theta, signsReflection['X'])
      self.theta = torch.einsum('i,i->i', self.theta, signsReflection['Y']) #(180 - theta) and not -theta
      self.C = getCMatrix(self.mstrType, 0.01+self.size, self.theta)
      areaCons = (torch.mean(self.C['area'])/optParams['desiredAreaFraction'])- 1.
      perimCons = 1. - (torch.sum(self.C['perim'])/optParams['desiredPerimeter'])
      self.J, self.Uvel, self.Vvel,self.Pressure = self.FE.objectiveHandle(self.C)
      if(self.epoch == 0 or self.epoch == 10):
        self.J0 = self.J.item()
      self.constraints =  [perimCons] #areaCons , perimCons
      loss, self.alpha = computeLoss(self.J/self.J0,  self.constraints)  
      loss.backward(retain_graph = True);
      return loss# , areaCons , perimCons
    #--------------------------#
    loss_old = 100.
    loss = 0.
    
    for self.epoch in range(optParams['maxEpochs']):
      if abs(loss-loss_old)<0.00001:
        break
      loss_old = loss
      print(loss)
      penal = min(8.0, 1. + self.epoch*0.02)
      loss = closure()
      if(self.optimizationMethod == 'adam'):
          
          self.optimizer.step()
          
      else:

          self.optimizer.step(closure)
      for c in self.constraints:
        muVar = muVar +  self.alpha*2*c.item()

          
   

      self.mstrType, self.size, self.theta = self.topNet(self.xyF)
      # self.size = applyDensityProjection(self.size, self.densityProj)
      self.theta = torch.einsum('i,i->i', self.theta, signsReflection['X'])
      self.theta = torch.einsum('i,i->i', self.theta, signsReflection['Y']) #(180 - theta) and not -theta
      self.C = getCMatrix(self.mstrType, 0.01+self.size, self.theta)
      areaCons = (torch.mean(self.C['area'])/optParams['desiredAreaFraction'])- 1.
      perimCons = 1. - (torch.sum(self.C['perim'])/optParams['desiredPerimeter'])
      self.J, self.Uvel, self.Vvel,self.Pressure = self.FE.objectiveHandle(self.C)
      if(self.epoch == 0 or self.epoch == 10):
        self.J0 = self.J.item()
      self.constraints =  [ perimCons] 
      loss, self.alpha = computeLoss(self.J/self.J0,  self.constraints)  


      loss.backward(retain_graph = True);

      if(self.epoch% optParams['epochSnapShot'] == 0):
        convgHistory['epoch'].append(self.epoch)
        
        convgHistory['J'].append(self.J.item())
        areaf= torch.mean(self.C['area'])
        convgHistory['af'].append(areaf.item())
        exposed_perim = torch.sum(self.C['perim'])
        status = 'epoch {:d} J \t {:.2E} areaf {:.2F} perim {:.2F}'.\
                              format(self.epoch, self.J, areaf, exposed_perim)
        print(status)
        if(self.epoch%30 == 0):
          self.FE.mesh.plotField(to_torch(self.C['perim']), 'perimeter')
          # self.FE.mesh.plotField(to_torch(self.C['area']), 'area')
          # self.plotCompositeTopology(1) # comment this when not saving frames
      # if(savedNet['isDump']):
      #   # saveTopologyNetwork(fileName)
    self.plotVelocity(self.Uvel.detach().numpy(), self.Vvel.detach().numpy(), self.Pressure.detach().numpy())
    self.FE.mesh.plotField(to_torch(self.C['area']), 'area')
    NX = 1*int(np.ceil((self.FE.mesh.meshSpec['bb']['xmax']-\
        self.FE.mesh.meshSpec['bb']['xmin'])/self.FE.mesh.meshSpec['elemSize'][0]))
    NY = 1*int(np.ceil((self.FE.mesh.meshSpec['bb']['ymax']-\
        self.FE.mesh.meshSpec['bb']['ymin'])/self.FE.mesh.meshSpec['elemSize'][1]))
    plt.figure()
    plt.imshow(self.size.detach().numpy().reshape((NX, NY)))
    plt.title('size')
    return convgHistory
    
    
    

  #-----------------------#
  def plotVelocity(self, uVelocity, vVelocity, Pressure):
    plt.figure()
    P =plt.imshow((Pressure).\
    reshape((self.FE.mesh.meshSpec['nelx'], \
    self.FE.mesh.meshSpec['nely'])).T, origin = 'lower', cmap ='rainbow')
    v = np.linspace(np.amin(Pressure), np.amax(Pressure), 10, endpoint=True)
    plt.colorbar(ticks=v)
    plt.title('$Pressure$')
    plt.axis('equal')
  
    
    
    plt.figure()
    a =plt.imshow(uVelocity.\
    reshape((self.FE.mesh.meshSpec['nelx'], \
    self.FE.mesh.meshSpec['nely'])).T, origin = 'lower', cmap ='rainbow')
    plt.colorbar(a)
    plt.title('$u_X$')
    plt.axis('equal')
    
    plt.figure()
    b =plt.imshow(vVelocity.\
    reshape((self.FE.mesh.meshSpec['nelx'], \
    self.FE.mesh.meshSpec['nely'])).T, origin = 'lower', cmap ='rainbow')
    plt.colorbar(b)
    plt.title('$u_Y$')
    plt.axis('equal')

    plt.figure()
    netV = np.sqrt(uVelocity**2 + vVelocity**2)
    b =plt.imshow(netV.\
    reshape((self.FE.mesh.meshSpec['nelx'], \
    self.FE.mesh.meshSpec['nely'])).T, origin = 'lower', cmap ='rainbow')
    plt.colorbar(b)
    plt.title('||u||')
    plt.axis('equal')

   

    
    
    plt.show()
  #-----------------------#
  def plotCompositeTopology(self, res):
    xy = self.FE.mesh.generatePoints(res)
    xyR, signsReflection = applyReflection(to_torch(xy), self.symMap)
    if(self.fourierMap['isOn']):
      xyF = applyFourierMap(xyR, self.fourierMap)
    else:
      xyF = xyR
    mstrType, size, theta = self.topNet(xyF)
    C = self.getCMatrix(mstrType, 0.01+size, theta)
    C['area'] = 1.-C['area']
    theta = torch.einsum('i,i->i', theta, signsReflection['X'])
    theta = torch.einsum('i,i->i', theta, signsReflection['Y'])
    fillColors = ['white',  (0.5,0.5,0), (0,1,0), (0,0,1), (0,0,0),\
        (1,0,1), (0.5,0,0.5), (0,0.5,0.5), (0,0,0.5), (0,0.5,0),\
          (1,0,0), (0.5,0,0)]
    
    microstrImages = to_torch(sio.loadmat('./data_generation/fluid/microStrImages.mat')['microstructures'])
    
    NX = res*int(np.ceil((self.FE.mesh.meshSpec['bb']['xmax']-\
        self.FE.mesh.meshSpec['bb']['xmin'])/self.FE.mesh.meshSpec['elemSize'][0]))
    NY = res*int(np.ceil((self.FE.mesh.meshSpec['bb']['ymax']-\
        self.FE.mesh.meshSpec['bb']['ymin'])/self.FE.mesh.meshSpec['elemSize'][1]))
    nx, ny = microstrImages.shape[2], microstrImages.shape[3]
    compositeImg = torch.zeros((NX*nx, NY*ny))
    colorImg = torch.zeros((NX, NY))
    densityImg = torch.zeros((NX, NY))
    maxC = 0
    step = 0.01 # THIS IS THE STEP USED WHEN GEN THE MSTR IMAGES! HC FOR NOW
    cutOff = 10.98 # val above which its a dark square
    cutOffLower = -0.0725
    th = torch.zeros((theta.shape[0]))
    dens_np = torch.zeros((size.shape[0]))
    for elem in range(xy.shape[0]):
      cx = int((res*xy[elem,0])/self.FE.mesh.meshSpec['elemSize'][0])
      cy = int((res*xy[elem,1])/self.FE.mesh.meshSpec['elemSize'][1])
      dens_np[elem] = int(100.*(C['area'][elem]))
      densityImg[cx, cy] = dens_np[elem]
      if(size[elem] > cutOff):
        compositeImg[nx*cx:(cx+1)*nx, ny*cy:(cy+1)*ny] = torch.ones((nx, ny))
        colorImg[cx, cy] = 1
      elif(size[elem] < cutOffLower):
        compositeImg[nx*cx:(cx+1)*nx, ny*cy:(cy+1)*ny] = torch.zeros((nx, ny))
        colorImg[cx, cy] = 0
      else:
        mstrIdx = min(microstrImages.shape[1]-1, int((size[elem]//step)))
        densityImg[cx, cy] = mstrIdx
        mstrTypeIdx = torch.argmax(mstrType[elem,:])
        mstrimg = to_torch(microstrImages[mstrTypeIdx, mstrIdx,:,:])
        th[elem] = 90 + int(180*theta[elem]/torch.pi)
        mstrimg = to_torch(rotate(mstrimg,th[elem],reshape = False, order = 0, mode = 'nearest'))
        c = torch.argmax(mstrType[elem,:])+1
        if(c > maxC):
          maxC = c
        compositeImg[nx*cx:(cx+1)*nx, ny*cy:(cy+1)*ny] = (1.-mstrimg)*c
        colorImg[cx, cy] = c
    plt.figure()
    plt.imshow(compositeImg.T, cmap = colors.ListedColormap(fillColors[:maxC+1]),\
                interpolation='none',vmin=0, vmax=maxC, origin = 'lower') 
    
    ax = plt.gca()
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])
    # plt.legend(['squircle', 'fishscale', 'fishscale_half', 'square', 'traingle', 'circle', 'ellipse', 'plus', 'mucosa10', 'mucosa20'])
    plt.show()

    plt.savefig(f'./frames/top_{self.epoch:d}.pdf', dpi=300)

    plt.figure()
    plt.imshow(colorImg.T, cmap = colors.ListedColormap(fillColors[:maxC+1]),\
                interpolation='none',vmin=0, vmax=maxC, origin = 'lower') 
    ax = plt.gca()
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])
    plt.show()

    plt.figure()
    plt.imshow(100-densityImg.T, cmap = 'gray',\
                interpolation='none', origin = 'lower') 
    plt.colorbar()
    ax = plt.gca()
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])
    plt.show()

    plt.figure()
    xy = to_torch(xy)
    a = plt.quiver(xy[:,0].reshape((NX, NY)),\
          xy[:,1].reshape((NX, NY)), \
          np.sin(np.pi*th/180.),\
          np.cos(np.pi*th/180.),\
          dens_np,\
          cmap = 'rainbow', scale= 10, headwidth = 0.1)#plt.cm.gnuplot
    plt.show()

    plt.figure()
    a = plt.imshow((th).reshape((NX,NY)).T, cmap = 'jet', origin = 'lower')
    plt.title('orientation')
    plt.colorbar(a)
    plt.show()
    #-----------------------#
