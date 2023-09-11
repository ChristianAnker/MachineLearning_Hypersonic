import yt
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import csv

class dataset():
    def __init__(self, shape, case, alpha, dmax):
        self.shape = shape
        self.case  = case
        self.alpha = alpha
        self.dmax = dmax

# Show data
def plot(ds, type:str):
    slc = yt.SlicePlot(ds,'z',type)
    slc.show()


def load_data(filename):
    return yt.load(filename)

# Load shape
def load_shape(filename, position= np.array([0,0]), num = 200, flip = False):
    x_list = list()
    y_list = list()
    # Open the file in write mode
    with open(filename, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')

        # Write the data rows
        for row in reader:
            x_list.append( row[0] )
            y_list.append( row[1] )

    x_list.append( x_list[0] )
    y_list.append( y_list[0] )

    coordinates = np.array([x_list,y_list],dtype=float).T 
    coordinates -= coordinates.mean(axis=0)
    coordinates += np.array(position) #np.array([1.2, 1.5])

    # compute s = distance
    s = [0]
    for i in range(1,len(coordinates)):
        dis = np.linalg.norm(coordinates[i,:]-coordinates[i-1,:])
        s.append( s[i-1] + dis )
    coordinates = np.append(coordinates,np.array([s]).T,axis=1)
    num = len(coordinates)

    co = coordinates[:(num+1)//2,:]

    if flip:
        co[:,0] = np.flip(co[:,0])

    return co




#Shape = pd.DataFrame(coordinates+np.array([1.482,1.5]),columns=['x','y'])

def unfold(ds, coordinates,distance = 0.01, n=100, m = 100, l = 10, mirror = False):
    # Distance is: points used at a distance normal to the surface [m] 
    # n is number of zones in x, m is number of zones in y, l is number of zones away from the coordinate, points will be looked for 

    # Transform to dataset
    ad = ds.all_data()
    AllData  = pd.DataFrame({'x':ad['x'], 'y':ad['y'], 'T':ad['nvtp'], 'I':ad['nvki']})

    # Sort by a rectangle around the shape
    dis = 0.15
    maxx = max(coordinates[:,0])+dis
    minx = min(coordinates[:,0])-dis
    disx = maxx - minx
    maxy = max(coordinates[:,1])+dis
    miny = min(coordinates[:,1])-dis
    disy = maxy - miny
    Data = AllData[AllData['x']<maxx]
    Data = Data[Data['x']>minx]
    Data = Data[Data['y']<maxy]
    Data = Data[Data['y']>miny]
    Data = Data[Data['I']==0]
    points = np.array([Data['x'],Data['y'],Data['T']]).T
    
    # Map the remainer points into a grid
    map = np.empty((n,m),dtype = list)
    for i in range(n):
        for j in range(m):
            map[i,j] = []

    # id in x is: (x-xmin) / disX * n
    normpoints = (points - np.array([minx, miny, 0]) ) / np.array([disx/n, disy/m, 1])
    for i, normpoint in enumerate(normpoints):
        ix = int(normpoint[0])
        iy = int(normpoint[1])
        map[ix,iy] += [i]

    s = 0
    p_list = list()
    for k in range(coordinates.shape[0]-1):
        
        co = coordinates[k,:2]
        ny = coordinates[k+1,:2] - co
        dis = np.linalg.norm(ny)
        ny = ny/dis
        nx = np.array([ny[1],-ny[0]]) if mirror else np.array([-ny[1],+ny[0]])

        # Find mapping close to
        normco = (co - np.array([minx, miny]) ) / np.array([disx/n, disy/m])
        ix = int(normco[0])
        iy = int(normco[1])

        # List all points within the grid close to the coordinate
        ip_list = list()
        for i in range(ix-l, ix+l+1):
            for j in range(iy-l, iy+l+1):
                ip_list += map[i,j]
        
        # Compute the normal distance and tangent distance for each point found
        for ip in ip_list:
            point = points[ip,]
            vec = point[:2] - co
            vecx = np.dot(vec, nx)
            vecy = np.dot(vec, ny)
            if vecx < distance and vecx > 0 and vecy < dis and vecy > 0:
                p_list.append(np.array([vecx, vecy+s, point[0], point[1], point[2]]))
        s+= dis

    return np.array(p_list)

def gradient_line(boundary, dis_tangent, dis_normal, polynomial, num=300, remove = 0.0):
    y_list = np.linspace(0,max(boundary[:,1]), num)
    g_list = list()

    # Add reversed points
    revpoints = boundary[boundary[:,1] < dis_normal/2]
    revpoints[:,1] = -revpoints[:,1]
    boundary = np.concatenate([boundary, revpoints])

    revpoints = boundary[boundary[:,0] > max(boundary[:,1])-dis_normal/2]
    revpoints[:,1] = 2*max(boundary[:,1]) - revpoints[:,1]
    boundary = np.concatenate([boundary, revpoints])

    for y in y_list:
        points = boundary[boundary[:,1] > y-dis_normal/2, :]
        points = points[points[:,1] < y+dis_normal/2, :]
        points = points[points[:,0] < dis_tangent, :]
        points = points[points[:,0] > remove, :]
        g_list.append(compute_grad(points,polynomial))

    return np.array([y_list, g_list]).T

# Compute the gradient
def compute_grad(points:np.ndarray, polynomial:int, plot = False, id=0)->float:

    coeffs = np.polyfit(points[:,0],  points[:,4], polynomial)
    poly = np.poly1d(coeffs)
    # Generate x values for plotting the regression line
    x_fit = np.linspace(0, max(points[:,0]), 100)
    T_fit = poly(x_fit)
    
    if plot:
        plot.plot(points[:,0],points[:,4],'.',c=f'C{id}')

        plot.plot(x_fit,T_fit,c=f'C{id}')

    # Compute the gradient at the given x-value
    derivative_coeffs = np.polyder(coeffs)
    gradient = np.polyval(derivative_coeffs, min(points[:,0]))


    return gradient


# Make the shapes
def make_shape(filename, X, Y):
    num = len(X)
    # Mirror
    x_list = list()
    y_list = list()
    for i in range(num):
        x_list.append( X[i] )
        y_list.append( Y[i] )
    for i in range(num-2,0,-1):
        x_list.append( X[i] )
        y_list.append( -Y[i])

    # Open the file in write mode
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')

        # Write the data rows
        for x, y in zip(x_list, y_list):
            writer.writerow([x,y])

    print(f"REMEMBER! points must be anticlockwise and describe a symmetric geometry\nData has been written to {filename}")


def curvature_line(coordinate, dis_tangent, degree, num=100):
    y_list = np.linspace(0,coordinate[-1,2], num)
    c_list = list()
    for y in y_list:
        c = curvature(coordinate, y, dis_tangent, degree)
        c_list.append( np.concatenate(( [y], np.ravel(c)) ))
    
    c_list[0][2] = 0
    return np.array(c_list)

def curvature(coordinate, position, length, degree):

    # Last point before hitting s
    ids = len(coordinate[coordinate[:,2]-position <= 0])-1
    # Last point before hitting s+l
    idmax = len(coordinate[coordinate[:,2]-(position+length/2) < 0])
    # First point after hitting s-l
    idmin = len(coordinate[coordinate[:,2]-(position-length/2) < 0])

    if ids == len(coordinate)-1: # In the case the last point is the ids
        vec = coordinate[ids] - coordinate[ids-1]
    else:
        vec = coordinate[ids+1] - coordinate[ids]

    # vec[2] is the distance between points
    n_t = vec[:2] / vec[2]
    angle = (np.arctan2(n_t[0], n_t[1])) #% (2*np.pi) 
    n_p = np.array([-n_t[1],n_t[0]])
    rot = np.array([n_t,n_p]).T # Rotate

    mpoint = coordinate[ids,:2] + (position-coordinate[ids,2]) * n_t
    
    # Move and rotate all points
    points = coordinate[idmin:idmax,:2] - mpoint
    points = np.matmul(points, rot)

    coeffs = np.polyfit(points[:,0],  points[:,1], degree+1)

    # Compute the second gradient (ie curvature) at x=0 
    curvature = np.polyder(coeffs,2)
    
    return np.concatenate(([vec[2]],[angle], np.ravel(np.flip(curvature))))


def distance_angle(curve, position, angles):

    ids = len(curve[curve[:,0]-position <= 0])-1
    if ids == len(curve)-1: # In the case the last point is the ids
        return 0

    vec = curve[ids+1] - curve[ids]
    p = curve[ids] + vec/vec[0] * (position-curve[ids,0])
    curvature = list()
    curv = np.insert(curve[ids+1:], 0, p, axis=0)
    curv[:,0] -= position
    for angle in angles:
        a = angle + p[2]
        # find angle:
        i = 1
        while a > max(curv[i,2], curv[i-1,2]) or a < min(curv[i,2], curv[i-1,2]):
            i += 1
            if i == len(curv):
                curvature.append(1e6)
                break
        if i == len(curv): continue
        # interpolate between points
        
        dis0 = curv[i-1,0]
        a0 = curv[i-1,2]
        dis1 = curv[i, 0]
        a1 = curv[i,2]
        
        dis = (a - a0)/(a1 - a0) * (dis1 - dis0) + dis0
        curvature.append( dis )

    return np.array(curvature)

def extract_gradient(filename, etamax = 0.005, etamin = 0.0015, Dzeta = 0.01, p=1, Dcurve = 0.1, degree = 1, num=100, a = 1):
    coor = load_shape(filename+'/shape', [2,0])
    ds = load_data(filename+'/plt')

    l = 0.12
    boundary = unfold(ds, coor,l)
    gradient = gradient_line(boundary, etamax, Dzeta, p, remove=etamin, num=num)
    curv = curvature_line(coor, Dcurve, degree, num=num)

    if a == 1: # First version
        input = np.zeros((num, 3))
        input[0,0] = gradient[0,1] # set the initial gradient  as the first input
        input[1:,:] = curv[1:,1:]
        output = gradient[1:,1]
    elif a == 2: # Second version
        input = curv[:,1:]
        output = gradient[:,1]

    return input, output


import os
def make_visitheader(directory_path):

    # Get all folder names in the directory
    folder_names = [folder for folder in os.listdir(directory_path) if os.path.isdir(os.path.join(directory_path, folder))]

    # Create a file to store the folder names
    output_file_path = directory_path+'/movie.visit'

    with open(output_file_path, 'w') as file:
        # Write each folder name to the file
        for folder_name in folder_names:
            file.write(folder_name + '\Header\n')

    print(f"Folder names have been saved to '{output_file_path}'.")


def boundary_scatter(boundary, coordinates):
    fig,ax = plt.subplots()
    fig.set_figwidth(8)
    fig.set_figheight(6)
    ax.scatter(boundary[:,2], boundary[:,3], s=50, c=boundary[:,4], label='Data',marker = '.')
    ax.plot(coordinates[:,0], coordinates[:,1],'r')
    # Add a zo

    # Customize the plot
    ax.set_title('Scatter Plot')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    #ax.axis('equal')
    ax.legend()