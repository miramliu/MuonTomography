# functions used in Geometric Trajectory and Detector Placement
#generally in order of use in notebooks
import numpy as np
import math

#global coordinates
#pyramid coordinates
p1 = np.linspace(-115.15, 115.15,100)
p2 = np.empty(100)
p2.fill(115.15)
p3 = np.empty(100)
p3.fill(-115.15)


#global functions
#  plotting average angles as unit vectors along the z direction (rather than along the y direction)
def unit(x,y,z):
    M = np.sqrt(x**2 + y**2 + z**2)
    #M=1
    return x/M, z/M,y/M #flipping y and z to project vertically

# SinogramSpace: takes the location of the far detector pixel's i and j (in ij space i.e. coordinate axis of the detector), 
#takes the dx, dy, and dz between the two pixels in question, and returns the location (xi_0 and phi) in sinogram space.
def SinogramSpace(i,j,dx,dy,dz):
    if dx == 0:
        psi = np.pi/2
    else:
        psi = math.atan(dz/dx)
        if psi<=0:
            psi = math.atan(dz/dx)+np.pi
    phi = psi-(np.pi/2)
    x = (i-(960/2))/100 # to be in meters not cm??
    y = -((230.33/2)+25) # in meters
    xi0 = x*math.cos(phi)+y*math.sin(phi)
    return phi, xi0, psi, x, y

# ok plotting this given a point and angle as a vector of length 300
def VectorPlot(x1,y1,psi):
    # given a point and an angle in radians, get the equation. 
    r = 300 #length
    if psi == np.pi/2:
        X = np.empty(50)
        X.fill(x1)
        Y = np.linspace(y1,r+y1,50)
    else:
        Y = []
        X = []
        m = math.tan(psi)
        xlen = r*math.cos(psi)
        for x in np.linspace(x1,xlen+x1,50):
            y = m*(x-x1)+y1
            Y.append(y)
            X.append(x)
    return (X,Y)

#rotate the radon coordinate system by phi (with respect to the coordinate system of the pyramid) into x-y space!
def Rotatexy(x,y,phi):
    xi = x*math.cos(phi)- y*math.sin(phi)
    eta = x*math.sin(phi) + y*math.cos(phi)
    return xi,eta
def RadonCoordinates(phi): #create the radon coordinates rotated
    x = np.linspace(-250,250,100)
    p = np.zeros(100)
    y = np.linspace(-250,250,100)
    xi_x, xi_y = Rotatexy(x,p,phi)
    eta_x, eta_y = Rotatexy(p,y,phi)
    return xi_x, xi_y, eta_x, eta_y #these are the radon transfer coordinates IN x-y spaace 

# 3d plot of trajectories. converts from detector coordinate to pyramid coordinate
def VectorPlot3D(i,j,dx,dy,dz,L):
    #L = 300 #length
    x1 = (i-(960/2))/100 #x in pyramid coordinate (meters)
    y1 = -((230.33/2)+25)# y in pyramid coordinate (meters)
    z1 = j # z in pyramid coordinate
    if dy==0:
        M = L/r #but when would it be 0??
    else: 
        M = L/dy # scale factor
    if dx ==0: #straight forward
        x = np.empty(50)
        x.fill(x1)
        y = np.linspace(y1,M*dy+y1,50)
        z = np.linspace(z1,M*dz,50)
    if dz ==0: #straight forward
        z = np.empty(50)
        z.fill(j)
        x = np.linspace(x1,M*dx,50)
        y = np.linspace(y1,M*dy+y1,50)
    else:  
        x = np.linspace(x1,M*dx,50)
        y = np.linspace(y1,M*dy+y1,50)
        z = np.linspace(z1,M*dz,50)
    return (x,y,z)

#find index of closest value in an array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#distance travelled by a (3D) vector to reach a plane. converts from detector coordinate to pyramid coordinate
#discarded for now. using similar triangles instead.
def rDistance(i,j,dx,dy,dz,plane):
    x1 = (i-(960/2))/100 #x in pyramid coordinate (meters)
    y1 = -((230.33/2)+25)# y in pyramid coordinate (meters)
    L =plane-y1 #distance between y and vertical plane 
    z1 = j # z in pyramid coordinate
    if dy==0:
        M = L/r #but when would it be 0?? gotta check this. 
    else: 
        M = L/dy # scale factor
    if dx ==0: #straight forward
        x = np.empty(50)
        x.fill(x1)
        y = np.linspace(y1,M*dy+y1,50)
        z = np.linspace(z1,M*dz,50)
        DY = (y1-M*dy+y1) #should always be 300
        DZ = (z1-M*dz)
        DX = 0
    if dz ==0: #straight forward
        z = np.empty(50)
        z.fill(j)
        x = np.linspace(x1,M*dx,50)
        y = np.linspace(y1,M*dy+y1,50)
        DX = (x1-M*dx)
        DY = y1-M*dy+y1 #should still be 300 still
        DZ = 0
    else:  
        DX = (x1-M*dx)
        DY = y1-M*dy+y1
        DZ = z1=M*dz
    r = np.sqrt(DX**2 + DY**2 + DZ**2)
    return r

#given location and solid angle, make the lines of the pixel.
def SAPixels(xC,zC,SR):
    L = .5*np.sqrt(SR)
    z1 = zC - L
    z2 = zC + L
    x1 = xC - L
    x2 = xC + L
    zL1 = np.empty(10)
    zL1.fill(z1)
    zL2 = np.empty(10)
    zL2.fill(z2)
    xL1 = np.empty(10)
    xL1.fill(x1)
    xL2 = np.empty(10)
    xL2.fill(x2)
    x = np.linspace(min(x1,x2),max(x1,x2),10)
    z = np.linspace(min(z1,z2),max(z1,z2),10)
    Lines = [xL1,xL2,zL1,zL2]
    return x,z,Lines


# "ShiftedRotated" =SR, and SRSinogramSpace: takes the location of the FAR detector pixel's i and j (in ij space i.e. 2d coordinate axis of the detector with bottom left as origin and cm), 
#takes rotation of detector with respect to the Y axis in DEGREES. (i.e. facing straight ahead is 0, facing along negative x-axis is -90 degrees)
#takes the dx, dy, and dz between the two pixels in question, and the location of the CENTER of the detector's FURTHEST panel in PYRAMID space (i.e. 3d coordinate axis with center of pyramid as origin and in m)
#and returns the location (xi_0 and phi) in sinogram space.

#this only does BOTTOM HORIZONTAL SLICE OF PYRAMID. (meaning it is one row of pixels: no tilt or j>0)
def SRSinogramSpace(i,j,X,Y,theta,dx,dy,dz):
    theta = theta*np.pi/180 #translation to theta in radians
    if dx == 0:
        psi = np.pi/2 + theta #directly forward in detector system
    else:
        psi = math.atan(dz/dx)
        if psi<=0:
            psi = math.atan(dz/dx)+np.pi
        psi = psi- theta
    phi = psi-(np.pi/2)
    x = (i-(960/2))/100 + X # in pyramid space
    y = Y # in pyramid space
    z = 0 #in pyramid space, with j = 0
    xi0 = x*math.cos(phi)+y*math.sin(phi)
    return phi, xi0, psi, x, y

#calculates sinogram space for one location XY of the detector
def OneSinogramSpace(i,j,X,Y,theta):
    phis = []
    xis = []
    psis = []
    x1s = []
    y1s = []
    for Hpix in i: #for all 480 pixels along the x axis
        for H in i: #for every pixel on the second detector x axis
            # for every combination of pixels
            di = Hpix-H #horizontal relative displacement (in cm)
            dj = 0 #zero vertical relative displacement
            Vpix = 0#by default
            phi_n, xi0_n,psi,x1,y1 = SRSinogramSpace(Hpix, Vpix,X,Y, theta,di,dj,200)
            phis.append(phi_n)
            xis.append(xi0_n)
            psis.append(psi)
            x1s.append(x1)
            y1s.append(y1)
    return phis,xis,psis,x1s,y1s
