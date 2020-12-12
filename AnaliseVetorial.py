import numpy as np
import sympy as sym

#retangulares
x, y, z = sym.symbols('x y z')

#cilindricas
rho, phi = sym.symbols('rho, phi')

#esféricas
r, theta = sym.symbols('r, theta')


def DivCar(V):
    """V = [ax, ay, az]"""
    V = np.array(V)
    ans = np.array([sym.diff(V[0], x),
        sym.diff(V[1], y),
        sym.diff(V[2], z)])
    return np.sum(ans)

def DivCil(V):
    """V = [a_rho, a_phi, a_z]
    Ângulos em redianos!!!"""
    V = np.array(V) 
    ans = np.array([(1/rho)*sym.diff(rho*V[0], rho),
        (1/rho)*sym.diff(V[1], phi),
        sym.diff(V[2], z)])
    return np.sum(ans)

def DivEsf(V):
    """V = [a_r, a_theta, a_phi]
    Ângulos em redianos!!!"""
    V = np.array(V) 
    ans = np.array([(1/r**2)*sym.diff(r**2*V[0], r),
                    (1/(r*sym.sin(theta)))*sym.diff(sym.sin(theta)*V[1], theta),
                    (1/(r*sym.sin(theta)))*sym.diff(V[2], phi)])
    return np.sum(ans)

    
def GradCar(f):
    """função"""
    ans = np.array([sym.diff(f, x),
                    sym.diff(f, y),
                    sym.diff(f, z)])
    return ans

def GradCil(f):
    """função (rad)"""
    ans = np.array([sym.diff(f, rho),
                    (1/rho)*sym.diff(f, phi),
                    sym.diff(f, z)])
    return ans

def GradEsf(f):
    """função (rad)"""
    ans = np.array([sym.diff(f, r),
                    (1/r)*sym.diff(f, theta),
                    (1/(r*sym.sin(theta))*sym.diff(f, phi))])
    return ans


def RotCar(V):
    """V = [ax, ay, ay]"""
    V = np.array(V)
    rot = np.array([sym.diff(V[2], y) - sym.diff(V[1], z) ,
         sym.diff(V[0], z) - sym.diff(V[2], x) ,
         sym.diff(V[1], x) - sym.diff(V[0], y)])
    return rot

def RotCil(V):
    """V = [a_rho, a_phi, a_z]
    Ângulos em redianos!!!"""
    V = np.array(V)
    rot = np.array([sym.diff(V[2], phi)/rho - sym.diff(V[1], z) ,
         sym.diff(V[0], z) - sym.diff(V[2], rho) ,
         (sym.diff(rho*V[1], rho) - sym.diff(V[0], phi))/rho])
    return rot

def RotEsf(V):
    """V = [a_r, a_theta, a_phi]
    Ângulos em redianos!!!"""
    V = np.array(V)
    rot = np.array([(sym.diff(V[2]*sym.sin(theta), theta) - sym.diff(V[1], phi))/(r*sym.sin(theta)) ,
        (sym.diff(V[0], phi)/sym.sin(theta) - sym.diff(r*V[2], r))/r ,
        (sym.diff(r*V[1], r) - sym.diff(V[0], theta))/r])
    return rot



def ModSym(V):
    mod = V/sym.sqrt(np.sum(V**2))
    return mod

def ProdVetCar(U, V):
    ans = np.array([U[1]*V[2] - U[2]*V[1], U[2]*V[0] - U[0]*V[2], U[0]*V[1] - U[1]*V[0]])
    return ans

    
