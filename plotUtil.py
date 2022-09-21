import matplotlib.pyplot as plt
import numpy as np

def plotConvergence( convg):
  x = np.array(convg['epoch'])
  for key in convg:    
    if(key == 'epoch'):
      continue # epoch is x axis for all plots
    plt.figure()
    y = np.array(convg[key])
    plt.semilogy(x, y, label = str(key))
    plt.xlabel('Iterations')
    plt.ylabel(str(key))
    plt.grid('True')


