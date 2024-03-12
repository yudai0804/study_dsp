import os, subprocess, math
import matplotlib as plt
import numpy as np

def print_ndarray(nd):
    [print('{:.16f}'.format(_nd)) for _nd in nd]

mode = ['sin', 'cordic_sin']
cordic_n = 18
dx = 0.01
start = 0.0
end = math.pi / 2
# main.pyのあるディレクトリに移動
os.chdir(os.path.dirname(__file__))
# cmake&実行
subprocess.run('cmake -S . -B build', shell=True)
subprocess.run('cmake --build build', shell=True)
for i in range(len(mode)):
    subprocess.run(['./build/main', os.path.join('data', mode[i]+'.csv'), str(cordic_n), mode[i], str(start), str(end), str(dx)])

_sin_nd = np.loadtxt('data/sin.csv', delimiter=',', dtype='float64')
_cordic_sin_nd = np.loadtxt('data/cordic_sin.csv', delimiter=',', dtype='float64')

sin_nd = np.zeros(len(_sin_nd), dtype='float64')
cordic_sin_nd = np.zeros(len(_sin_nd), dtype='float64')
x = np.zeros(len(_sin_nd), dtype='float64')

for i in range(len(x)):
    x[i] = _sin_nd[i][0]
    sin_nd[i] = _sin_nd[i][1]
    cordic_sin_nd[i] = _cordic_sin_nd[i][1]

diff = np.zeros(len(sin_nd), dtype='float64')
for i in range(len(sin_nd)):
    diff[i] = (sin_nd[i] - cordic_sin_nd[i])
print_ndarray(diff)