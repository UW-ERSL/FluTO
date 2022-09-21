import numpy as np
from mesher import Q2Q1Mesher
import configparser
from bc import getBC
from projections import computeFourierMap
import matplotlib.pyplot as plt



from TOuNN import TOuNN
from FE_Fluid import FE
from plotUtil import plotConvergence

import time
start_time = time.time()



#%% read config file
configFile = './config.txt'
config = configparser.ConfigParser()
config.read(configFile)

#%% Mesh and BC
meshConfig = config['MESH']
nelx, nely = meshConfig.getint('nelx'), meshConfig.getint('nely')
elemSzX, elemSzY = meshConfig.getfloat('elemSizeX'), meshConfig.getfloat('elemSizeY')
meshSpec = {'ndim':2, 'nelx':nelx, 'nely':nely,\
            'elemSize':np.array([elemSzX,elemSzY]),\
            'numElems':nelx*nely, 'elemArea':elemSzX*elemSzY}

mesh = Q2Q1Mesher(meshSpec)


#%% BC
bc = getBC(example = 2, mesh = mesh) #1= Diffuser 2=bent pipe  3= Double Pipe

#%% Material
materialConfig = config['MATERIAL']
mu =  materialConfig.getfloat('mu')# mu is viscosity
matProp = {'physics':'fluid', 'mu':mu}


#%% NN
tounnConfig = config['TOUNN']
nnSettings = {'numLayers': tounnConfig.getint('numLayers'),\
              'numNeuronsPerLayer':tounnConfig.getint('hiddenDim'),\
              'outputDim':tounnConfig.getint('outputDim')}
  
fourierMap = {'isOn':tounnConfig.getboolean('fourier_isOn'),\
              'minRadius':tounnConfig.getfloat('fourier_minRadius'), \
              'maxRadius':tounnConfig.getfloat('fourier_maxRadius'),\
              'numTerms':tounnConfig.getint('fourier_numTerms')}

fourierMap['map'] = computeFourierMap(mesh, fourierMap)


#%% Optimization params
lossConfig = config['LOSS']
# lossMethod = {'type':'logBarrier', 't0':lossConfig.getfloat('t0'),\
#               'mu':lossConfig.getfloat('mu')}
lossMethod = {'type':'augLag', 'alpha0':lossConfig.getfloat('alpha0'), \
              'delAlpha':lossConfig.getfloat('delAlpha')}
          
optConfig = config['OPTIMIZATION']
optParams = {'maxEpochs':optConfig.getint('numEpochs'),\
              'lossMethod':lossMethod,\
              'learningRate':optConfig.getfloat('lr'),\
              'desiredAreaFraction':optConfig.getfloat('desiredAreaFraction'),\
              'desiredPerimeter':optConfig.getfloat('desiredPerimeter'),\
              'epochSnapShot':optConfig.getint('epochSnapShot'),\
              'gradclip':{'isOn':optConfig.getboolean('gradClip_isOn'),\
                          'thresh':optConfig.getfloat('gradClip_clipNorm')}}

#%% Run optimization
symMap = {'XAxis':{'isOn':True, \
          'midPt': 0.5*meshSpec['elemSize'][1]*meshSpec['nely']},\
          'YAxis':{'isOn':False, \
          'midPt': 0.5*meshSpec['elemSize'][0]*meshSpec['nelx']}}
plt.close('all')
savedNet = {'isAvailable':False, 'file':'./netWeights.pkl', 'isDump':False}
overrideGPU = True
tounn = TOuNN( mesh, matProp, FE, nnSettings, fourierMap, symMap, bc, 'lbfgs', overrideGPU)
convgHistory = tounn.optimizeDesign(optParams, savedNet)
tounn.plotCompositeTopology(1)
plotConvergence(convgHistory)
# tounn.plotCompositeTopology(1)
# plt.show(block=True)
print("--- %s seconds ---" % (time.time() - start_time))