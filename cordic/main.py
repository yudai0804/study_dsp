import os, subprocess, math
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def print_ndarray(nd):
    [print('{:.16f}'.format(_nd)) for _nd in nd]

mode = ['sin', 'cordic_sin']
cordic_n = 18
dx = 0.001
start = -math.pi / 2
end = math.pi / 2
# start = -0.9
# end = 0.9
# main.pyのあるディレクトリに移動
os.chdir(os.path.dirname(__file__))
# cmake&実行
subprocess.run('cmake -S . -B build', shell=True)
subprocess.run('cmake --build build', shell=True)
for i in range(len(mode)):
    subprocess.run(['./build/main', os.path.join('data', mode[i]+'.csv'), str(cordic_n), mode[i], str(start), str(end), str(dx)])

_y1 = np.loadtxt(os.path.join('data', mode[0]+'.csv'), delimiter=',', dtype='float64')
_y2 = np.loadtxt(os.path.join('data', mode[1]+'.csv'), delimiter=',', dtype='float64')

y1 = np.zeros(len(_y1), dtype='float64')
y2 = np.zeros(len(_y1), dtype='float64')
x = np.zeros(len(_y1), dtype='float64')

for i in range(len(x)):
    x[i] = _y1[i][0]
    y1[i] = _y1[i][1]
    y2[i] = _y2[i][1]

diff = np.zeros(len(y1), dtype='float64')
for i in range(len(y1)):
    diff[i] = (y1[i] - y2[i])
# print_ndarray(diff)

# y1 = y1 * 180 / np.pi
# y2 = y2 * 180 / np.pi

fig, ax = plt.subplots()
# ax.plot(x,y1)
# ax.plot(x,y2)
ax.plot(x, diff)
plt.show()