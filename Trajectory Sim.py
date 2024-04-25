#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook


t, Fm=np.genfromtxt('RocketFileG1.csv', delimiter=',', skip_header=1, usecols=(0,1), unpack=True)
diameter, l, m, mp, Cd, Pdiameter, PCd,  h0,  z0,  ivw, wtht, tht=np.genfromtxt('RocketFileG1.csv', delimiter=',', skip_header=1,usecols=(2,3,4,5,6,7,8,9,10,11,12,13), unpack=True, max_rows=1)

p0=101325 #Air pressure
T0=288.15 #Temperature
g=9.81 #Gravity
L= 0.0065 #temperature lapse rate,
R=8.31 #Universal Gas Constant
M=0.0289652 #molar mass of dry air
rho=1.225 #Initial Air Density 

axlist=[0] #Lateral Acceleration
aylist=[0] #Vertical Acceleration
vxlist=[0] #Lateral Velocity
vylist=[0] #Vertical Velocity
hlist=[0] #Height (Vertical Distance)
dlist=[float(0)] #Lateral Distance
wtht=math.radians(wtht) #Wind Angle from lateral (East-South)
A=3.1415*(0.5*diameter)**2 #Nosecone Area
PA=3.1415*(0.5*Pdiameter)**2 #Parachute Area
LA=diameter*l #Cross Sectional Side-area of Rocket
t=t.tolist() #Converting from Array to List
tburn=t[-1] #Total Burn Time
t2=np.arange(t[-1],200,0.1) 
t.extend(t2) #Creating a list spanning the time of mission
Z=np.zeros(10000) 
Fm=np.concatenate((Fm,Z)) #Thrust is zero after burn time
Fmx=Fm*math.cos(math.radians(tht)) 
#Seperating into vertical and lateral components
Fmy=Fm*math.sin(math.radians(tht))

for n in range(len(t)-1):
    
    rho=((p0*M)/(R*T0))*((1-(L*(hlist[n]+h0))/T0)**((g*M)/(R*L)-1))
    #Air Density as a function of height
    vw=ivw*((np.log((hlist[n]+h0)/z0))/(np.log((h0)/z0)))
    #Wind speed as a function of height
    
    dt=t[n+1]-t[n] #Time Differential
    if t[n]==tburn: 
        m=m-mp #Mass of propellant is assumed to all be expelled
        
    if (hlist[n]-hlist[n-1])>=0: #Before Apogee
        Fdy=0.5*rho*vylist[n]**2*Cd*A
        #Drag Components in x and y direction
        Fdx=0.5*rho*vxlist[n]**2*Cd*A
    elif (hlist[n]-hlist[n-1])<0: #Detects Apogee for Parachute Deployment
        Fdy=0.5*rho*vylist[n]**2*PCd*PA
        Fdx=0.5*rho*vxlist[n]**2*Cd*LA
        
    if vxlist[n]>=0: #Finds appropriate direction for drag
        axlist.append((Fmx[n+1]-Fdx)/m)
        vxlist.append(vxlist[n]+axlist[n]*dt)
    else:
        axlist.append((Fmx[n+1]+Fdx)/m)
        vxlist.append(vxlist[n]+axlist[n]*dt)
    if vylist[n]>=0:
        aylist.append((Fmy[n+1]-Fdy)/m-g)
        vylist.append(vylist[n]+aylist[n]*dt)
    else:
        aylist.append((Fmy[n+1]+Fdy)/m-g)
        vylist.append(vylist[n]+aylist[n]*dt)
        
    hlist.append(hlist[n]+vylist[n]*dt+0.5*aylist[n]*dt**2) 
    #Finite Difference Equations
    dlist.append(dlist[n]+vxlist[n]*dt-vw*dt+0.5*axlist[n]*dt**2)
    
    if hlist[n]<0:
        break

print("Apogee reached:", np.max(hlist), "m" )
print("Landing Site Distance:", dlist[-1], "m")
print(np.max(dlist))

plt.plot(t[:len(hlist)],hlist, 'o', label="Vertical Distance") 
plt.plot(t[:len(dlist)],dlist, 'o', label="Lateral Distance")
plt.xlabel("Time (s)")
plt.ylabel("Distance (m)")

plt.legend()
plt.show()

    
    
    
    
    
    
    
    
    
    
    


# In[5]:


fig = plt.figure()
ax = plt.axes(projection='3d')
dmod=np.array(dlist)
dmod=dmod.tolist()
N=math.cos(wtht)
E=math.sin(wtht)

z=hlist
y=[x*N for x in dmod]
x=[x*E for x in dmod]

fig=ax.plot3D(x, y, z, 'r')
ax.set_xlabel('Lateral Distance (m)', fontsize=10)
ax.set_ylabel('Longitudinal Distance (m)',  fontsize=10)
ax.set_zlabel('Height (m)')
ax.set_box_aspect(aspect=None, zoom=0.8)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
ax.zaxis.set_tick_params(labelsize=8)



# In[6]:


with cbook.get_sample_data(r"C:\Users\Knifi\Jupyter\Coding\ExampleMap2 (2) (1).png") as image_data:
    image = plt.imread(image_data)


fig, ax = plt.subplot_mosaic([
    ['Map']
], figsize=(7, 3.5))

ax['Map'].imshow(image)




markers = [(x[0]+500, y[0]+500), (x[-1]+500,y[-1]+500)]
s0,sf = zip(*markers)
ax['Map'].plot(s0, sf, 'rx')
ax['Map'].set_xlabel('Latitudinal Distance (m)')
ax['Map'].set_ylabel('Longitudinal Distance (m)')
ax['Map'].set_xlim(0,1001)
ax['Map'].set_ylim(0,1001)

plt.show()


# In[4]:


thtlist=[]
dmaxlist=[0]
dlist=[0]
hlist=[0]
vxlist=[0]
vylist=[0]
axlist=[0]
aylist=[0]
m=0.824
for tht in np.arange(91,70,-0.1):
    Fmx=Fm*math.cos(math.radians(tht))
    Fmy=Fm*math.sin(math.radians(tht))
    for n in range(len(t)-1):
        rho=((p0*M)/(R*T0))*((1-(L*(hlist[n]+h0))/T0)**((g*M)/(R*L)-1))
        vw=ivw*((np.log((hlist[n]+h0)/z0))/(np.log((h0)/z0)))
        dt=t[n+1]-t[n]
        if t[n]==tburn:
            m=m-mp
        if (hlist[n]-hlist[n-1])>=0:
            Fdy=0.5*rho*vylist[n]**2*Cd*A
            Fdx=0.5*rho*vxlist[n]**2*Cd*A
        elif (hlist[n]-hlist[n-1])<0:
            Fdy=0.5*rho*vylist[n]**2*PCd*PA
            Fdx=0.5*rho*vxlist[n]**2*Cd*LA
        if vxlist[n]>=0:
            axlist.append((Fmx[n+1]-Fdx)/m)
            vxlist.append(vxlist[n]+axlist[n]*dt/2)
        else:
            axlist.append((Fmx[n+1]+Fdx)/m)
            vxlist.append(vxlist[n]+axlist[n]*dt/2)
        if vylist[n]>=0:
            aylist.append((Fmy[n+1]-Fdy)/m-g)
            vylist.append(vylist[n]+aylist[n]*dt/2)
        else:
            aylist.append((Fmy[n+1]+Fdy)/m-g)
            vylist.append(vylist[n]+aylist[n]*dt/2)
        hlist.append(hlist[n]+vylist[n]*dt+0.5*aylist[n]*dt**2)
        dlist.append(dlist[n]+vxlist[n]*dt-vw*dt+0.5*axlist[n]*dt**2)
        if hlist[n]<0:
            thtlist.append(tht)
            dmaxlist.append(dlist[-1])
            break
    if dmaxlist[-1]>0:
        break
    else:
        dlist=[0]
        hlist=[0]
        vxlist=[0]
        vylist=[0]
        axlist=[0]
        aylist=[0]
        m=0.824
print("Angle to Minimise Lateral Distance:", round(thtlist[-1],2), "degrees")


# In[ ]:





# In[ ]:




