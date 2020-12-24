import MDSplus as mds
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


def get_signal(connection, tree, shot, tag):
  # get 1D signal from MDSplus tree
  connection.openTree(tree, shot)
  data = connection.get(tag).value
  time = connection.get('dim_of(' + tag + ')').value
  connection.closeAllTrees()
  return time, data


shot = 204660
connection_path = 'skylark.pppl.gov'
connection = mds.Connection(connection_path)
tree = 'activespec'

# t,tstime = get_signal(connection, tree, shot, '.MPTS.OUTPUT_DATA.BEST.TS_TIMES')
t,fit_te = get_signal(connection, tree, shot, '.MPTS.OUTPUT_DATA.BEST:FIT_TE')


print(t)
print(fit_te.shape)
plt.plot(t,fit_te)
plt.show()