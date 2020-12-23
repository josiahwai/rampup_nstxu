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
tree = 'efit01'


signals = ['rcur', 'zcur']
tags = ['.RESULTS.AEQDSK:RCUR', '.RESULTS.AEQDSK:ZCUR']
save_fn = './rz_cur.mat'

signal_data = {}
for i, signal in enumerate(signals):
  t, signal_data[signal] = get_signal(connection, tree, shot, tags[i])
  plt.plot(t, signal_data[signal])

plt.legend(signals)
plt.show()

signal_data['time'] = t
sio.savemat(save_fn, signal_data)



'''
tree = 'activespec'
t,fit_te = get_signal(connection, tree, shot, '.MPTS.OUTPUT_DATA.BEST:FIT_TE')
plt.plot(t,fit_te.T)
plt.xlim([0,1])
plt.ylim([0,1.8])
plt.show()
'''


'''
tags = []
tags.append('.CHERS.ANALYS.CT1:TI' )  # ion temp
tags.append('.CHERS.ANALYS.CT1:ATI')  # apparent ion temp
tags.append('.CHERS.ANALYS.CT1:TIS')  # spline ion temp
tags.append('.CHERS.ANALYS.CT1:ZTI')  # zeeman corrected ion temp
tags.append('.CHERS.ANALYS.CT2:TI' )
tags.append('.CHERS.ANALYS.CT2:ATI')
tags.append('.CHERS.ANALYS.CT2:TIS')
tags.append('.CHERS.ANALYS.CT2:ZTI')

t,Ti = get_signal(connection, tree, shot, '.CHERS.ANALYSIS.CT1:TIS')
plt.plot(t,Ti.T)
plt.show()
'''