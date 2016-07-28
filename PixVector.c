#include <math.h>

typedef struct VectorFloat {
    float x, y, z;
} VecF;

typedef struct VectorInteger {
    int x, y, z;
} VecI;

VecF ViToVf(VecI i) {
    VecF result;
    result.x = i.x;
    result.y = i.y;
    result.z = i.z;
    return result;
}

VecI VfToVi(VecF v) {
    VecI result;
    result.x = v.x;
    result.y = v.y;
    result.z = v.z;
    return result;
}

VecF VmulF(VecF v, float f) {
    VecF result;
    result.x = v.x * f;
    result.y = v.y * f;
    result.z = v.z * f;
    return result;
}

VecF VmulV(VecF a, VecF b) {
    VecF result;
    result.x = a.x * b.x;
    result.y = a.y * b.y;
    result.z = a.z * b.z;
    return result;
}

VecI VmulVi(VecF f, VecI i) {
    VecI result;
    result.x = f.x * i.x;
    result.y = f.y * i.y;
    result.z = f.z * i.z;
    return result;
}

VecF VimulV(VecI i, VecF f) {
    VecF result;
    result.x = f.x * i.x;
    result.y = f.y * i.y;
    result.z = f.z * i.z;
    return result;
}

VecF VdivV(VecF a, VecF b) {
    VecF result;
    result.x = a.x / b.x;
    result.y = a.y / b.y;
    result.z = a.z / b.z;
    return result;
}

VecF VdivVi(VecF a, VecI b) {
    VecF result;
    result.x = a.x / b.x;
    result.y = a.y / b.y;
    result.z = a.z / b.z;
    return result;
}

VecF VaddF(VecF v, float f) {
    VecF result;
    result.x = v.x + f;
    result.y = v.y + f;
    result.z = v.z + f;
    return result;
}

VecF VaddV(VecF a, VecF b) {
    VecF result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

VecF VminV(VecF a, VecF b) {
    VecF result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}

VecF VminVi(VecF v, VecI i) {
    VecF result;
    result.x = v.x - i.x;
    result.y = v.y - i.y;
    result.z = v.z - i.z;
    return result;
}

VecF Vnorm(VecF v) {
    VecF result;
    float mod = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    result.x = v.x/mod;
    result.y = v.y/mod;
    result.z = v.y/mod;
    return result;
}

VecF VnormWithLen(VecF v, float* l) {
    VecF result;
    float mod = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    result.x = v.x/mod;
    result.y = v.y/mod;
    result.z = v.y/mod;
    *l = mod;
    return result;
}

float Vlen(VecF v) {
    return sqrt(v.x*v.x + v.y*v.y + v.z+v.z);
}
