import ctypes
import numpy
import pathlib
import math

libfile = pathlib.Path(__file__).parent / "p3plib.so"
p3pclib = ctypes.CDLL(str(libfile))

type_pvec2 = ctypes.POINTER(ctypes.c_double * 2)
type_vec3 = ctypes.c_double * 3
type_vec2 = ctypes.c_double * 2
type_pvec3 = ctypes.POINTER(type_vec3)
type_mat33 = ((ctypes.c_double * 3) * 3)

p3pclib.p3p_lambdatwist.restype = ctypes.c_int
p3pclib.p3p_lambdatwist.argtypes = [type_vec2, type_vec2, type_vec2, type_vec3, type_vec3, type_vec3, type_mat33*4, type_vec3*4, ctypes.c_int]

p3pclib.toHomography.argtypes = [ctypes.POINTER(ctypes.c_double * 8),
    ctypes.POINTER(type_mat33), type_pvec3]

p3pclib.applyHomography.argtypes = [type_vec2, ctypes.c_double * 8, type_vec2]

p3pclib.applyRT.argtypes = [type_vec2, type_mat33, type_vec3, type_vec3]

def lambdatwist(y1, y2, y3, x1, x2, x3, iterations=5):
    Rs = (type_mat33 * 4)()
    Ts = (type_vec3 * 4)()
    _x1 = (type_vec3)(*x1) if len(x1)>=3 else (type_vec3)(x1[0],x1[1],0)
    _x2 = (type_vec3)(*x2) if len(x1)>=3 else (type_vec3)(x2[0],x2[1],0)
    _x3 = (type_vec3)(*x3) if len(x1)>=3 else (type_vec3)(x3[0],x3[1],0)
    valid = p3pclib.p3p_lambdatwist(type_vec2(*y1),type_vec2(*y2),type_vec2(*y3),
        _x1,_x2,_x3,Rs,Ts,iterations)
    return Rs[:valid],Ts[:valid]
    
def toHomography(R,T):
    h = (ctypes.c_double * 8)()
    p3pclib.toHomography(h, R, T)
    return h

def homographies(y1,y2,y3,x1,x2,x3,iterations=5):
    Rs,Ts = lambdatwist(y1,y2,y3,x1,x2,x3,iterations=iterations)
    out = []
    for i in range(len(Rs)):
        out.append(toHomography(Rs[i], Ts[i]))
    return out

def applyRT(R,T,v):
    _v = (type_vec3)(*v) if len(v)>=3 else (type_vec3)(v[0],v[1],0)
    out = (type_vec2)()
    p3pclib.applyRT(out, R, T, _v)
    return out

def p4p(dest,source,iterations=5):
    Rs,Ts = lambdatwist(*dest[:3],*source[:3],iterations=iterations)
    if len(Rs) == 0:
        return None,None
    bestI = -1
    bestD = float("inf")
    for i in range(len(Rs)):
        out = applyRT(Rs[i],Ts[i],source[3])
        d = math.hypot(out[0]-dest[3][0],out[1]-dest[3][1])
        if d < bestD:
            bestI = i
            bestD = d
    if bestI is None:
        return None,None
    return Rs[bestI],Ts[bestI]
    
def applyHomography(h,x):
    _x = (type_vec2)(*x[:2])
    out = (type_vec2)()
    p3pclib.applyHomography(out,h,_x)
    return out
        
def p3p_test():
    h = homographies([0,0],[0,1],[1,0],[0,0,0],[0,2,0],[1.5,0,0], 5)
    print(len(h),list(h[0]))
    print(list(applyHomography(h[0],[0,2,0])))
