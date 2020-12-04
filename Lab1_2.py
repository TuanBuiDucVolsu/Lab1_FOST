import easygui
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.animation as animation


msg = "Enter velocity & angle"
title = "Input for projectile simulation"
fieldNames = ["Velocity x","Velocity y","Angle","k"]
fieldValues = []  
fieldValues = easygui.multenterbox(msg,title, fieldNames)
print(fieldValues)

k = float(fieldValues[3]) 
m = 0.5 
g = 9.81 
Vx0 = float(fieldValues[0]) 
Vy0 = float(fieldValues[1]) 
thetha = float(fieldValues[2])  * np.pi / 180.0  
V0 = np.sqrt(Vx0**2+Vy0**2) 

t0 = 0.0 
tmax = 2.0  
steps = 100  
tLF = np.linspace(t0, tmax, steps+1)

y0LF = [0.0, Vx0, 0.0, Vy0]  

def deriv(yLF,tF): #Creates the function which computes the derivatives 
    Vx = yLF[1]   # Identifies the velocity along x axis
    Vy = yLF[3]   #Identifies the velocity along y axis
    return [Vx, -k*Vx/m, Vy, -g-k*Vy/m]   # The function outputs [dx/dt, dVx/dt, dy/dt, dVy/dt]

yM = odeint(deriv, y0LF, tLF)  #  The 4-D array containing the solution for the differential equation

plt.plot(yM[:,0],yM[:,2],'.',label='Numerical solution') #Plots y over x numerically.

#Analytical Solution:
VT = m*g/k #Вычисляет конечную скорость
anal_x = ((V0*VT)/g)*np.cos(thetha)*(1-np.exp(-g*tLF/VT)) 
anal_y = (VT/g)*(V0*np.sin(thetha)+VT)*(1-np.exp(-g*tLF/VT))-VT*tLF 

plt.plot(anal_x,anal_y,label='Analytical solution') #Plots y over x analytically.
plt.grid()
plt.xlabel('Horizontal axis [meters]')
plt.ylabel('Vertical axis [meters]')
plt.legend()
plt.title('Projectile motion with air resistance.')
plt.savefig('ballistic.pdf')
#Computing answer of the questions:

i = np.abs(yM[1:steps,2]).argmin() #calculates the point from the y-array closest to 0 (after the initial point)and assigns it to i (=impact).

D = yM[i+1,0] #Вычисляет расстояние до точки удара.

H = np.amax(yM[:,2]) #Находит максимальное значение массива anal_y, представляющего наивысшую точку траектории, и присваивает его H.

TF = tLF[i+1]

Vxi = yM[i+1,1]
Vyi = yM[i+1,3]
Vi = np.sqrt(Vxi**2+Vyi**2)

#Creates the file .txt and saves the answers of the questions
f = open('Lab1_2.txt','w')
f.write(' The distance to the point of impact is {}.\n The highest point of the trajectory is at {} meters above ground.\n The time for flight is {} s. \n The impact velocity is {} m/s.'.format(str(D),str(H),str(TF),str(Vi)))
f.close()


#Extra: animation of the projectile motion.
fig, ax = plt.subplots() #sets the initial figure.
line, = ax.plot(anal_x,anal_y)   # Sets the values of 'line'

def animate(k): #Creates the animation function, which updates the data in each frame. 
    line.set_data(anal_x[:k],anal_y[:k])
    return line,

plt.axis([0.0, 20, 0.0, 10]) #Sets the range of the axis. 

ani = animation.FuncAnimation(fig, animate,100,interval=50,blit=True) #Creates the animation. 
plt.title('Animation of projectile motion with air resistance.')
#Adds an arrow pointing at the landing spot:
#plt.annotate('landing point', xy=(D, 0), xytext=(1.8, 0.8), arrowprops=dict(facecolor='blue', shrink=0.02),)
plt.show()
