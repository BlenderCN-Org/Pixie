#define false 0
#define true 1

#define LOOK(grid, xx, yy, zz) grid[zz*gridDim.y*gridDim.x+yy*gridDim.x+xx]

#include <stdio.h>
#include <math.h>

#include "PixVector.c"

typedef struct Pixels {
    float r;
    float g;
    float b;
    float a;
} Pixel;


struct Camera {
    VecF pos;
    VecF rot;
    int dimx;
    int dimy;
    float focalLength;
    float sensorWidth;
};

struct Volume {
    VecF boundLower;
    VecF boundUpper;
    VecI gridDim;
    float* grid;
    float densityMul;
};

struct Light {
    VecF pos;
    float intensity;
};


struct Camera pCam;
struct Volume pVol;
struct Light pLgt;


int setCamera(int resx, int resy, float posx, float posy, float posz,
              float rotx, float roty, float rotz, float focalLength,
              float sensorWidth) {
    pCam.dimx = resx;
    pCam.dimy = resy;
    
    pCam.pos.x = posx;
    pCam.pos.y = posy;
    pCam.pos.z = posz;
    
    pCam.rot.x = rotx;
    pCam.rot.y = roty;
    pCam.rot.z = rotz;
    
    pCam.focalLength = focalLength;
    pCam.sensorWidth = sensorWidth;
    
    return true;
}

int setVolumeWithPos(float posx, float posy, float posz,
                     float scalex, float scaley, float scalez,
                     int dimx, int dimy, int dimz, float grid[], float densityMul) {
    pVol.boundLower.x = posx - scalex/2;
    pVol.boundLower.y = posy - scaley/2;
    pVol.boundLower.z = posz - scalez/2;
    
    pVol.boundUpper.x = posx + scalex/2;
    pVol.boundUpper.y = posy + scaley/2;
    pVol.boundUpper.z = posz + scalez/2;
    
    pVol.gridDim.x = dimx;
    pVol.gridDim.y = dimy;
    pVol.gridDim.z = dimz;
    
    pVol.grid = grid;
    pVol.densityMul = densityMul;
    
    return true;
}

int setVolume(float blx, float bly, float blz, float bux, float buy, float buz,
              int dimx, int dimy, int dimz, float grid[], float densityMul) {
    pVol.boundLower.x = blx;
    pVol.boundLower.y = bly;
    pVol.boundLower.z = blz;
    
    pVol.boundUpper.x = bux;
    pVol.boundUpper.y = buy;
    pVol.boundUpper.z = buz;
    
    pVol.gridDim.x = dimx;
    pVol.gridDim.y = dimy;
    pVol.gridDim.z = dimz;
    
    pVol.grid = grid;
    pVol.densityMul = densityMul;
    
    return true;
}

int setLight(float posx, float posy, float posz, float intensity) {
    pLgt.pos.x = posx;
    pLgt.pos.y = posy;
    pLgt.pos.z = posz;
    
    pLgt.intensity = intensity;
    
    return true;
}


int rayBoxIntersection(VecF rayOrig, VecF rayDirc, VecF bndLower, VecF bndUpper,
                       float* tminR, float* tmaxR) {
    // rox/y/z - ray origin x/y/z
    // rdx/y/z - ray direction x/y/z
    // blx/y/z - box lower x/y/z
    // bux/y/z - box upper x/y/z
    // tminR, tmaxR - The distance to entering the box and leaving the box
    //                   in terms of the VecF rdx/y/z.
    
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    
    VecF bounds[2];
    bounds[0] = bndLower;
    bounds[1] = bndUpper;
    
    VecF invdir;
    invdir.x = 1/rayDirc.x;
    invdir.y = 1/rayDirc.y;
    invdir.z = 1/rayDirc.z;
    
    int sign[3];
    sign[0] = (invdir.x < 0);
    sign[1] = (invdir.y < 0);
    sign[2] = (invdir.z < 0); 
    
    tmin = (bounds[sign[0]].x - rayOrig.x) * invdir.x;
    tmax = (bounds[1-sign[0]].x - rayOrig.x) * invdir.x;
    tymin = (bounds[sign[1]].y - rayOrig.y) * invdir.y;
    tymax = (bounds[1-sign[1]].y - rayOrig.y) * invdir.y;

    if ((tmin > tymax) || (tymin > tmax))
    return false;
    if (tymin > tmin)
    tmin = tymin;
    if (tymax < tmax)
    tmax = tymax;

    tzmin = (bounds[sign[2]].z - rayOrig.z) * invdir.z;
    tzmax = (bounds[1-sign[2]].z - rayOrig.z) * invdir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
    return false;
    if (tzmin > tmin)
    tmin = tzmin;
    if (tzmax < tmax)
    tmax = tzmax;

    *tminR = tmin;
    *tmaxR = tmax;

    return true;
}

int rotationMatrix(float x, float y, float z, float result[3][3]) {
    //printf("%f, %f %f\n", x, y, z);
    float cx = cos(x);
    float cy = cos(y);
    float cz = cos(z);
    float sx = sin(x);
    float sy = sin(y);
    float sz = sin(z);
    result[0][0] = cy*cz;
    result[0][1] = cz*sx*sy-cx*sz;
    result[0][2] = cx*cz*sy+sx*sz;
    result[1][0] = cy*sz;
    result[1][1] = cx*cz+sx*sy*sz;
    result[1][2] = -cz*sx+cx*sy*sz;
    result[2][0] = -sy;
    result[2][1] = cy*sx;
    result[2][2] = cx*cy;
    //printf("R %f,%f,%f,%f,%f,%f,%f,%f,%f", result[0][0], result[0][1], result[0][2], result[1][0], result[1][1], result[1][2], result[2][0], result[2][1], result[2][2]); 
    return true;
}

Pixel castRayInVolumeOld(VecF start, VecF direc, VecF boundLower, VecF boundUpper,
                      VecI gridDim, float grid[], float stepSize, float tDiff) {
    VecF direction = VmulF(direc, stepSize);
    Pixel result = {1.0, 1.0, 1.0, 0.0};
    
    int samples = (int) (tDiff/stepSize) + 1;
    for (int i=0; i<samples; i++) {
        VecF p = VdivV(VminV(VaddV(VmulF(direction, i), start), boundLower), VminV(boundUpper, boundLower));
        
        VecI co = VmulVi(p, gridDim);
        
        float density = grid[co.z*pVol.gridDim.y*pVol.gridDim.x+co.y*pVol.gridDim.x+co.x];
        
        result.r += (1-result.a) * stepSize * density;
        result.g += (1-result.a) * stepSize * density;
        result.b += (1-result.a) * stepSize * density;
        
        result.a += (1-result.a) * stepSize * density;
    }
    return result;
}


Pixel castRayInVolume(VecF start, float startDist, VecF direc, VecF boundLower, VecF boundUpper,
                      VecI gridDim, float grid[], float stepSize,
                      int lightRay) {
    float dMul = pVol.densityMul;
    if (lightRay) {
        stepSize *= 50;
        dMul *= 1;
    }
    VecF direction = VmulF(direc, stepSize);
    //printf("direction: %f, %f, %f\n", direction.x, direction.y, direction.z);
    Pixel result = {0.0, 0.0, 0.0, 1.0};
    VecF current; current.x = start.x; current.y = start.y; current.z = start.z;
    
    // start will be on the very edge of the volume and there can be numerical errors
    if (current.x > boundUpper.x) current.x = boundUpper.x;
    if (current.y > boundUpper.y) current.y = boundUpper.y;
    if (current.z > boundUpper.z) current.z = boundUpper.z;
    
    if (current.x < boundLower.x) current.x = boundLower.x;
    if (current.y < boundLower.y) current.y = boundLower.y;
    if (current.z < boundLower.z) current.z = boundLower.z;
    
    while (current.x <= boundUpper.x && current.x >= boundLower.x &&
           current.y <= boundUpper.y && current.y >= boundLower.y &&
           current.z <= boundUpper.z && current.z >= boundLower.z &&
           result.a > 0.0) {
        //if (firstSample) printf("c: %f, %f, %f; l: %f, %f, %f; u: %f, %f, %f;\n", current.x, current.y, current.z, boundLower.x, boundLower.y, boundLower.z, boundUpper.x, boundUpper.y, boundUpper.z);
        VecF p = VdivV(VminV(current, boundLower), VminV(boundUpper, boundLower));
        VecF fco = VimulV(gridDim, p);
        VecF rco = VaddF(fco, .5);
        
        float density;
        
        if (rco.x <= 1.0 || rco.x >= (gridDim.x-1) || rco.y <= 1.0 || rco.y >= (gridDim.y-1) ||
            rco.z <= 1.0 || rco.z >= (gridDim.z-1)) {
            // For sample that are on the very edge and so can't use
            //   trilinear interpolation.
            density = 0;
        } else {
            VecI co = VfToVi(rco);
            // Trilinear interpolation
            VecF cellCo = VminVi(rco, co);
            float c000 = LOOK(grid, co.x, co.y, co.z);
            float c001 = LOOK(grid, co.x, co.y, (co.z+1));
            float c010 = LOOK(grid, co.x, (co.y+1), co.z);
            float c011 = LOOK(grid, co.x, (co.y+1), (co.z+1));
            float c100 = LOOK(grid, (co.x+1), co.y, co.z);
            float c101 = LOOK(grid, (co.x+1), co.y, (co.z+1));
            float c110 = LOOK(grid, (co.x+1), (co.y+1), co.z);
            float c111 = LOOK(grid, (co.x+1), (co.y+1), (co.z+1));
            
            float c00 = c000*(1-cellCo.x) + c100*cellCo.x;
            float c01 = c001*(1-cellCo.x) + c101*cellCo.x;
            float c10 = c010*(1-cellCo.x) + c110*cellCo.x;
            float c11 = c011*(1-cellCo.x) + c111*cellCo.x;
            
            float c0 = c00*(1-cellCo.y) + c10*cellCo.y;
            float c1 = c01*(1-cellCo.y) + c11*cellCo.y;
            
            float c = c0*(1-cellCo.z) + c1*cellCo.z;
            
            density = c * dMul;
            if (density > 1.0) density = 1.0;
            
            //Nearest (ie. no interpolation)
            //density = grid[co.z*pVol.gridDim.y*pVol.gridDim.x+co.y*pVol.gridDim.x+co.x];
        }
        
        //printf("%f, %f, %f; \n", p.x, p.y, p.z);
        //if (p.x == 1.0) co.x--;
        //if (p.y == 1.0) co.y--;
        //if (p.z == 1.0) co.z--;
        //printf("%i %i %i\n", co.x, co.y, co.z);
        //float density = grid[co.z*pVol.gridDim.y*pVol.gridDim.x+co.y*pVol.gridDim.x+co.x];
        
        //printf("one\n");
        if (!lightRay && density > 0.0) {
            //printf("not lightray\n");
            float distToLight;
            VecF lRay = VnormWithLen(VminV(pLgt.pos, current), &distToLight);
            Pixel secondary = castRayInVolume(current, 0, lRay, boundLower, boundUpper,
                                              gridDim, grid, stepSize, true);
            float intensity = pLgt.intensity * secondary.a / (distToLight*distToLight);
            
            result.r += result.a * density * intensity;
            result.g += result.a * density * intensity;
            result.b += result.a * density * intensity;
        }
        
        result.a *= (1 - density);

        current = VaddV(current, direction);
    }
    //printf("r %f; g %f; b %f; a %f\n", result.r, result.g, result.b, result.a);
    //result.a = pow(1-result.a, stepSize);
    if (!lightRay) {
        result.a = 1 - result.a;
    }
    return result;
}


/*int renderVolume(int width, int height, double d_camPosx, double d_camPosy, double d_camPosz,
                 double d_camRotx, double d_camRoty, double d_camRotz, float camPlane,
                 float blx, float bly, float blz, float bux, float buy, float buz,
                 int dimx, int dimy, int dimz, float lightx, float lighty,
                 float lightz, float lightIntens, float grid[], float results[]) {*/
int renderVolume(float results[]) {
    float MCamRot[3][3];
    //printf("rot: %f, %f, %f\n", pCam.rot.x, pCam.rot.y, pCam.rot.z);
    rotationMatrix(pCam.rot.x, pCam.rot.y, pCam.rot.z, MCamRot);
    
    float pixelSize = pCam.sensorWidth/((float) pCam.dimx);
    float xOffset = -(pCam.sensorWidth/2) + pixelSize/2;
    float yOffset = -(pCam.dimy-1)*pixelSize/2;
    for (int w=0; w<pCam.dimx; w++) {
        float x = xOffset + w*pixelSize;
        for (int h=0; h<pCam.dimy; h++) {
            float z = yOffset + h*pixelSize;
            float mod = sqrt(x*x + pCam.focalLength*pCam.focalLength + z*z);
            float normX = x/mod;
            float normY = pCam.focalLength/mod;
            float normZ = z/mod;
            
            //printf("M %f,%f,%f,%f,%f,%f,%f,%f,%f", MCamRot[0][0], MCamRot[0][1], MCamRot[0][2], MCamRot[1][0], MCamRot[1][1], MCamRot[1][2], MCamRot[2][0], MCamRot[2][1], MCamRot[2][2]); 
            
            float tranX = MCamRot[0][0]*normX + MCamRot[0][1]*normY + MCamRot[0][2]*normZ;
            float tranY = MCamRot[1][0]*normX + MCamRot[1][1]*normY + MCamRot[1][2]*normZ;
            float tranZ = MCamRot[2][0]*normX + MCamRot[2][1]*normY + MCamRot[2][2]*normZ;
            
            float tminR, tmaxR;
            VecF direction = {tranX, tranY, tranZ};
            //printf("s: %f, %f, %f, t: %f, %f, %f\n", normX, normY, normZ, tranX, tranY, tranZ);
            int intersection = rayBoxIntersection(pCam.pos, direction,
                                                  pVol.boundLower, pVol.boundUpper,
                                                  &tminR, &tmaxR);

            if (intersection) {
                VecF firstIntersect = VaddV(pCam.pos, VmulF(direction, tminR));
                Pixel result = castRayInVolume(firstIntersect, tminR, direction,
                                               pVol.boundLower, pVol.boundUpper,
                                               pVol.gridDim, pVol.grid, 0.008,
                                               false);
                /*Pixel result = castRayInVolumeOld(firstIntersect, direction,
                                               pVol.boundLower, pVol.boundUpper,
                                               pVol.gridDim, pVol.grid, 0.1, tmaxR-tminR);*/
                results[w*pCam.dimy*4+h*4+0] = result.r;
                results[w*pCam.dimy*4+h*4+1] = result.g;
                results[w*pCam.dimy*4+h*4+2] = result.b;
                results[w*pCam.dimy*4+h*4+3] = result.a;
                //results[w*pCam.dimy*3+h*3+0] = result.r * result.a;
                //results[w*pCam.dimy*3+h*3+1] = result.g * result.a;
                //results[w*pCam.dimy*3+h*3+2] = result.b * result.a;
            } else {
                results[w*pCam.dimy*4+h*4+0] = 0.0;
                results[w*pCam.dimy*4+h*4+1] = 0.0;
                results[w*pCam.dimy*4+h*4+2] = 0.0; 
                results[w*pCam.dimy*4+h*4+3] = 0.0;                
            }
        }
    }
    return true;
}
                 
