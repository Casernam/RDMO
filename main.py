# Media Onda No Controlado con fuente sinusoidal
# vf (t) = p2  v  sen(wt)
# Carga tipo RL

# Entrada de datos
# V= input ('Tension efectiva de la fuente sinusoidal ');
# R= input ('Resistencia [Ohm] ');
# L= input ('Inductancia [H] ');
# f= input ('Frecuencia de la fuente [Hz] ');
from cProfile import label
import math
from turtle import color
import numpy
import matplotlib.pyplot as plt
from scipy import signal

import numpy as np

def cust_range(*args, rtol=1e-05, atol=1e-08, include=[True, False]):
    """
    Combines crange and numpy.isclose to mimic
    open, half-open and closed intervals.
    Avoids also floating point rounding errors as with
    >>> crange(1, 1.3, 0.1)
    array([1. , 1.1, 1.2, 1.3])

    args: [start, ]stop, [step, ]
        as in crange
    rtol, atol: floats
        floating point tolerance as in numpy.isclose
    include: boolean list-like, length 2
        if start and end point are included
    """
    # process arguments
    if len(args) == 1:
        start = 0
        stop = args[0]
        step = 1
    elif len(args) == 2:
        start, stop = args
        step = 1
    else:
        assert len(args) == 3
        start, stop, step = tuple(args)

    # determine number of segments
    n = (stop-start)/step + 1

    # do rounding for n
    if np.isclose(n, np.round(n), rtol=rtol, atol=atol):
        n = np.round(n)

    # correct for start/end is exluded
    if not include[0]:
        n -= 1
        start += step
    if not include[1]:
        n -= 1
        stop -= step

    return np.linspace(start, stop, int(n))

def crange(*args, **kwargs):
    return cust_range(*args, **kwargs, include=[True, True])

def orange(*args, **kwargs):
    return cust_range(*args, **kwargs, include=[True, False])


def applyFunction(array,function):
  toSin = numpy.vectorize(function)
  return toSin(array)


V= float(input ("Ingrese valor del voltaje (V): "))
R= float(input("Ingrese valor resistencia (Ohm): "))
L= float (input("Ingrese valor inductancia (H): "))
f= float(input("Ingrese valor frecuencia de la fuente (Hz):"))

while V <= 0:
  V= float(input ("Ingrese valor del voltaje (V): "))
while R <=0:
  R= float(input("Ingrese valor resistencia (Ohm): "))
while L <=0:
  L= float (input("Ingrese valor inductancia (H): "))
while f <=0:
  f= float(input("Ingrese valor frecuencia de la fuente (Hz):"))

 # Parametros
fi= math.atan (2* math.pi*f*L/R)         # ángulo de desfase
Z= math.sqrt ((2* math.pi*f*L)**2+R**2)    #calcular la impedancia 
# condición de regimen permanente 
I02= math.sqrt(2)*V*math.sin(fi)/Z*(1+math.exp(-math.pi/math.tan(fi)))/(math.exp(math.pi/math.tan(fi))-math.exp(-math.pi/math.tan(fi)))
I01= I02*math.exp(math.pi/math.tan(fi))
if I02 < 0:
  I02 =0
  I01= math.sqrt(2)*V/Z*(math.sin(math.pi -fi)+math.sin(fi)*math.exp(-( math.pi)/math.tan(fi)))

# Función en el tiempo
t1 = crange(0, math.pi, 0.001)

t2 = crange(math.pi, 2* math.pi,0.001)

def factorFi(fi, t1):
    return applyFunction(-(t1)/math.tan(fi),math.exp)

id1= math.sqrt (2)*V/Z*(applyFunction(t1-fi,math.sin) + 
  math.sin(fi)*factorFi(fi,t1)) + I02*factorFi(fi, t1)
id2=( I01 * applyFunction(-((t2 - math.pi) / math.tan(fi)),math.exp))
t=[t1,t2]
# Rizado
Rizado =( max(id1)-min(id2))/2
# Corriente media y efectiva de los diodos y la carga
Io_d1 =1/(2* math.pi)* numpy.trapz (id1, t1)
Io_d2 =1/(2* math.pi)* numpy.trapz (id2, t2)
Io= Io_d1 + Io_d2
Irms_d1 = math.sqrt (1/(2* math.pi)* numpy.trapz (id1**2,t1) )
Irms_d2 = math.sqrt (1/(2* math.pi)* numpy.trapz (id2**2, t2) )
Irms = math.sqrt ( Irms_d1 ** 2+ Irms_d2 ** 2)
# tension media y efectiva
Vo= math.sqrt (2)*V/math.pi
Vrms =V/math.sqrt(2)
# Factor de rizado
FR_i = math.sqrt ( Irms ** 2 - Io ** 2) / Io
FR_v = math.sqrt ( Vrms ** 2- Vo ** 2) / Vo
# Graficas
fig, ax = plt.subplots(1)
t2 = crange( 0, 2*math.pi,0.001)

vx= math.sqrt(2)*V*applyFunction(t2,math.sin)*(signal.square(t2))*0.5
vf= math.sqrt(2)*V*applyFunction(t2,math.sin)
ix = numpy.concatenate ((id1 ,id2))

plt.plot(t2 ,vf,linestyle="--",label="Fuente")
plt.xlabel('Tiempo')
plt.ylabel('Magnitud')
plt.plot(t2 ,vx,color="r",label="Carga")
ax.set_xticks(crange(0, 2*math.pi, math.pi/6))
ax.set_xticklabels(['0','T/12','T/6','T/4','T/3','5T/12','T/2','7T/12','2T/3','3T/4','5T/6','11T/12','T'])
ax.set_yticklabels([" "])



fig1, ax1 = plt.subplots(1)
plt.plot (t2 ,ix)
plt.xlabel('Tiempo')
plt.ylabel('Magnitud') 
ax1.set_xticks(crange(0, 2*math.pi, math.pi/6))
ax1.set_xticklabels(['0','T/12','T/6','T/4','T/3','5T/12','T/2','7T/12','2T/3','3T/4','5T/6','11T/12','T'])
ax1.set_yticklabels([" "])


fig2, ax2 = plt.subplots (1)
plt.plot (t2 ,vx , color='r')
plt.show()


#Desarrollado por César Serna