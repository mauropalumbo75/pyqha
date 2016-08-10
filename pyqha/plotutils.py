# -*- coding: utf-8 -*-

################################################################################
#
# Some tests with matplotlib

import matplotlib.pyplot as plt
from read import read_alpha

def plot_x(T,x,fignum):
  plt.figure(fignum)
  plt.subplot(211)
  plt.xlabel('T (K)')
  plt.ylabel('a (a.u.) ')
  plt.plot(T, x[0], 'r')
  plt.subplot(212)
  plt.xlabel('T (K)')
  plt.ylabel('c/a ')
  plt.plot(T, x[1], 'r')

def plot_alpha(T,alpha,fignum):
  plt.figure(fignum)
  plt.subplot(211)
  plt.xlabel('T (K)')
  plt.ylabel(r'$\alpha_x$')
  plt.plot(T, alpha[0], 'r')
  plt.subplot(212)
  plt.xlabel('T (K)')
  plt.ylabel(r'$\alpha_z$')
  plt.plot(T, alpha[1], 'r')
  plt.show()

################################################################################
#   MAIN
################################################################################

if __name__ == "__main__":
  T, x, alpha = read_alpha("Os/output_anharm.dat.python.quad")
  
  plot_x(T,x,1)
  plot_alpha(T,alpha,2)
  



