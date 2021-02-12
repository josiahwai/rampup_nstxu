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
  
  
def get_conductor_data(connection,tree,shot,ip_threshold=200000):
  t,ccefit = get_signal(connection, tree, shot, '.RESULTS.AEQDSK:CCBRSP')
  t,ecefit = get_signal(connection, tree, shot, '.RESULTS.AEQDSK:ECCURT')
  t, ip = get_signal(connection, tree, shot, '.RESULTS.AEQDSK:IPMEAS')
  mask = np.where(ip>ip_threshold)
  ip = ip[mask]
  t = t[mask]
  coil_currents = np.zeros((len(mask[0]),13))
  vessel_currents = np.zeros((len(mask[0]), 40))
  coil_currents[:, 0] = ecefit[mask]#OH
  coil_currents[:, 1] = ccefit[mask, 0]#PF1AU
  coil_currents[:, 2] = ccefit[mask, 1]#PF1BU
  coil_currents[:, 3] = ccefit[mask, 2]#PF1CU
  coil_currents[:, 4] = ccefit[mask, 3]#PF2U
  coil_currents[:, 5] = ccefit[mask, 4]#PF3U
  coil_currents[:, 6] = (ccefit[mask, 5]+ccefit[mask, 8])/2.0#PF4
  coil_currents[:, 7] = (ccefit[mask, 6]+ccefit[mask, 7])/2.0#PF5
  coil_currents[:, 8] = ccefit[mask, 9]#PF3L
  coil_currents[:, 9] = ccefit[mask, 10]#PF2L
  coil_currents[:, 10] = ccefit[mask, 11]#PF1CL
  coil_currents[:, 11] = ccefit[mask, 12]#PF1BL
  coil_currents[:, 12] = ccefit[mask, 13]#PF1AL
  vessel_currents[:,:] = ccefit[mask,14:]  
  
  return coil_currents, vessel_currents, t
  
shot = 204660

connection_path = 'skylark.pppl.gov'
connection = mds.Connection(connection_path)
tree = 'efit01'
ic, iv, times = get_conductor_data(connection,tree,shot)

save_times = np.arange(20,1001,1) # ms
isave = []
for i, t in enumerate(save_times):
  isave.append(np.abs(times - t/1000).argmin())

plt.plot(times, ic)
plt.show()

coils = {'ic':ic[isave,:], 'iv':iv[isave,:], 't':times[isave], 'shot':shot}
fn = './data/coils' + str(shot) + '.mat' 
sio.savemat(fn, {'coils':coils})








