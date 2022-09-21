import numpy as np
import random
import torch
import torch.nn as nn


def set_seed(manualSeed):
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.manual_seed(manualSeed)
    torch.cuda.manual_seed(manualSeed)
    torch.cuda.manual_seed_all(manualSeed)
    np.random.seed(manualSeed)
    random.seed(manualSeed)
#%% Neural network
class TopNet(nn.Module):
    def __init__(self, nnSettings):
        self.inputDim =  nnSettings['inputDim'] # x and y coordn of the point
        super().__init__()
        self.layers = nn.ModuleList()
        manualSeed = 77 # NN are seeded manually
        set_seed(manualSeed);
        current_dim = self.inputDim;
        for lyr in range(nnSettings['numLayers']): # define the layers
            l = nn.Linear(current_dim, nnSettings['numNeuronsPerLayer'])
            nn.init.xavier_normal_(l.weight)
            nn.init.zeros_(l.bias)
            self.layers.append(l)
            current_dim = nnSettings['numNeuronsPerLayer']
        self.layers.append(nn.Linear(current_dim, nnSettings['outputDim']))
        self.bnLayer = nn.ModuleList()
        for lyr in range(nnSettings['numLayers']): # batch norm
            self.bnLayer.append(nn.BatchNorm1d(nnSettings['numNeuronsPerLayer']))

    def forward(self, x):
        m = nn.LeakyReLU(); #nn.SiLU()#  SiLU is swish
        ctr = 0;
        
        for layer in self.layers[:-1]: # forward prop
          
            x =  m(self.bnLayer[ctr](layer(x)))
            ctr += 1;
        
        nnOut = (self.layers[-1](x))
        sizeLower, sizeUpper = 0.1, 1.# size = 0 means no solid  #for fix size fish_scale 0.83 to 0.87
        theta = np.pi*torch.sigmoid(nnOut[:,0]).reshape(-1)

        size = (sizeLower + (sizeUpper - sizeLower)*(torch.sigmoid(nnOut[:,1])).reshape(-1))
        mstrType = torch.softmax(nnOut[:,2:], dim=1, dtype=None)
        return  mstrType, size, theta