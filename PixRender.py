from PIL import Image
import numpy as np
import math




from ctypes import cdll, c_float, c_double, c_int, byref, Structure, POINTER
accel = cdll.LoadLibrary('/home/peter/Desktop/Pixie/build/lib.linux-x86_64-3.5/RayCastAccel.cpython-35m-x86_64-linux-gnu.so')

accel.setCamera.argtypes = [c_int, c_int, c_float, c_float, c_float, c_float,
                            c_float, c_float, c_float, c_float]
accel.setCamera.restype = c_int

accel.setVolumeWithPos.argtypes = [c_float, c_float, c_float, c_float, c_float, c_float,
                                   c_int, c_int, c_int, POINTER(c_float), c_float]
accel.setVolumeWithPos.restype = c_int

accel.setVolume.argtypes = [c_float, c_float, c_float, c_float, c_float, c_float,
                            c_int, c_int, c_int, POINTER(c_float), c_float]
accel.setVolume.restype = c_int

accel.setLight.argtypes = [c_float, c_float, c_float, c_float]
accel.setLight.restype = c_int

accel.renderVolume.argtypes = [POINTER(c_float)]
accel.renderVolume.restype = c_int



class Ray:
    def __init__(self, orig, dirc):
        self.orig = orig
        self.dirc = dirc
        self.invdirc = [1/x if x != 0 else float("inf") for x in dirc]
        self.sign = [1 if x < 0 else 0 for x in self.invdirc]

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.matrix([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def rayBoxIntersection(ray, box):
    bnd = box.bounds
    
    # c implementation is about 2.5 times slower
    """rox = c_float(ray.orig[0])
    roy = c_float(ray.orig[1])
    roz = c_float(ray.orig[2])
    rdx = c_float(ray.dirc[0])
    rdy = c_float(ray.dirc[1])
    rdz = c_float(ray.dirc[2])
    blx = c_float(bnd[0][0])
    bly = c_float(bnd[0][1])
    blz = c_float(bnd[0][2])
    bux = c_float(bnd[1][0])
    buy = c_float(bnd[1][1])
    buz = c_float(bnd[1][2])
    tmin = c_float(0)
    tmax = c_float(0)

    result = accel.rayBoxIntersection(rox, roy, roz, rdx, rdy, rdz,
                                      blx, bly, blz, bux, buy, buz,
                                      byref(tmin), byref(tmax))
    
    tmin = max(tmin.value, 0)
    tmax = max(tmax.value, 0)
    
    return tmin, tmax"""
    
    tmin = (bnd[ray.sign[0]][0] - ray.orig[0]) * ray.invdirc[0];
    tmax = (bnd[1-ray.sign[0]][0] - ray.orig[0]) * ray.invdirc[0];
    tymin = (bnd[ray.sign[1]][1] - ray.orig[1]) * ray.invdirc[1];
    tymax = (bnd[1-ray.sign[1]][1] - ray.orig[1]) * ray.invdirc[1];

    if ((tmin > tymax) or (tymin > tmax)):
        return False
    if (tymin > tmin):
        tmin = tymin
    if (tymax < tmax):
        tmax = tymax

    tzmin = (bnd[ray.sign[2]][2] - ray.orig[2]) * ray.invdirc[2];
    tzmax = (bnd[1-ray.sign[2]][2] - ray.orig[2]) * ray.invdirc[2];

    if ((tmin > tzmax) or (tzmin > tmax)):
        return False
    if (tzmin > tmin):
        tmin = tzmin
    if (tzmax < tmax):
        tmax = tzmax

    return tmin, tmax


def castRay(ray, scene):
    intersection = rayBoxIntersection(ray, scene.volume)
    if intersection and intersection[1] > 0:
        tmin, tmax = intersection
        tmin = max(tmin, 0)
        return (tmin, tmax, 0.5, .5)
    else:
        return (0.0, 0.0, 0.0, 0.0)



def renderScene(scene):
    camera = scene.camera
    pos = camera.pos
    rot = camera.rot
    res = camera.resolution
    ang = np.deg2rad(camera.angle)
    foc = camera.focalLength
    swd = camera.sensorWidth
    lgt = scene.light

    im = Image.new("RGBA", (res[0], res[1]))
    pix = im.load()
    
    width = c_int(res[0])
    height = c_int(res[1])
    
    #print(pos)    
    #print([round(tmp, 1) for tmp in pos])
    #camPosx = c_float(round(pos[0], 1))
    #camPosy = c_float(round(pos[1], 1))
    #camPosz = c_float(round(pos[2], 1))
    #print(camPosx, camPosy, camPosz)
    camPosx = c_float(pos[0])
    camPosy = c_float(pos[1])
    camPosz = c_float(pos[2])
    
    camRotx = c_float(rot[0])
    camRoty = c_float(rot[1])
    camRotz = c_float(rot[2])
    
    camFocal = c_float(foc)
    camSensorWidth = c_float(swd)
    
    accel.setCamera(width, height, camPosx, camPosy, camPosz, camRotx, camRoty,
                    camRotz, camFocal, camSensorWidth)
    
    vol = scene.volume
    
    volPosx = c_float(vol.pos[0])
    volPosy = c_float(vol.pos[1])
    volPosz = c_float(vol.pos[2])
    
    volScalex = c_float(vol.scale[0])
    volScaley = c_float(vol.scale[1])
    volScalez = c_float(vol.scale[2])
    
    volDimx = c_int(vol.dim[0])
    volDimy = c_int(vol.dim[1])
    volDimz = c_int(vol.dim[2])
    
    c_float_p = POINTER(c_float)
    grid = scene.volume.grid.ctypes.data_as(c_float_p)
    
    densityMul = c_float(vol.densityMul)
    
    accel.setVolumeWithPos(volPosx, volPosy, volPosz, volScalex, volScaley, 
                           volScalez, volDimx, volDimy, volDimz, grid, densityMul)
    
    lightx = c_float(lgt.pos[0])
    lighty = c_float(lgt.pos[1])
    lightz = c_float(lgt.pos[2])
    
    lightIntensity = c_float(lgt.intensity)
    
    accel.setLight(lightx, lighty, lightz, lightIntensity)
    
    results = (c_float * (res[0] * res[1] * 4))()
    
    succ = accel.renderVolume(results)
    
    for w in range(res[0]):
        for h in range(res[1]):
            p = (int(255*results[w*res[1]*4+h*4+0]),
                 int(255*results[w*res[1]*4+h*4+1]),
                 int(255*results[w*res[1]*4+h*4+2]),
                 int(255*results[w*res[1]*4+h*4+3]))
            pix[w, h] = p
    
    return im
    
    """unit = np.array([0, 1, 0])
    
    Mx = rotation_matrix((1,0,0), np.deg2rad(rot[0]))
    My = rotation_matrix((0,1,0), np.deg2rad(rot[1]))
    Mz = rotation_matrix((0,0,1), np.deg2rad(rot[2]))
    Mrot = Mz*My*Mx
    camDir = unit*Mrot
    
    size = 2/res[0]
    xOffset = -1+size/2
    yOffset = -(res[1]-1)*size/2
    for w in range(res[0]):
        for h in range(res[1]):
            x = xOffset+w*size
            y = pln
            z = yOffset+h*size
            mod = math.sqrt(x**2 + y**2 + z**2)
            x /= mod
            y /= mod
            z /= mod
            #print((x, y, z))
            direction = np.array((x, y, z)) * Mrot
            
            ray = Ray(pos, (direction.item((0,0)), direction.item((0,1)),
                            direction.item((0,2))))
            result = castRay(ray, scene)
            p = tuple([int(128*result[a]) for a in range(3)])
            pix[w, h] = p"""
    
    
    """ays = [-ang*res[1]/res[0]/2 + y*ang/(res[0]-1) for y in range(res[1])]
    camMyCache = [camDir*rotation_matrix((1,0,0), ay) for ay in ays]
    
    for x in range(res[0]):
        ax = -ang/2 + x*(ang/(res[0]-1))
        Mx = rotation_matrix((0,0,1), ax)
        for y in range(res[1]):
            camMy = camMyCache[y]
            
            direction = camMy*Mx
            
            ray = Ray(pos, (direction.item((0,0)), direction.item((0,1)),
                            direction.item((0,2))))
            #result = castRay(ray, scene)
            
            print
            
            #p = tuple([int(128*result[a]) for a in range(3)])
            p = (int(255/2*(direction.item((0, 0))+1)),
                 int(255/2*(direction.item((0, 1))+1)),
                 int(255/2*(direction.item((0, 2))+1)))
            pix[x, y] = p"""
    
    """dim = scene.volume.dim

    for x in range(dim[2]):
        for y in range(dim[1]):
            total = 0
            for z in range(dim[0]):
                #total += nb[z*dim[0]*dim[1] + y*dim[0] + x]
                total += dat[z][y][x]
            v = int(255*total/dim[2])
            pix[x, y] = (v, v, v)"""
    
    return im
