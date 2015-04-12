import subprocess
import struct
import math

class Particle:
    fmtstr = struct.Struct("dddddddddddd")

    def __init__(self, string="", x=0.0, y=0.0, z=0.0,
                 vx=0.0, vy=0.0, vz=0.0,
                 ax=0.0, ay=0.0, az=0.0,
                 m=0.0, r=0.0, lastcollision=0.0):
        if len(string) > 0:
            self.string = string
            data = self.fmtstr.unpack(string)
            (self.x, self.y, self.z,
                self.vx, self.vy, self.vz,
                self.ax, self.ay, self.az,
                self.m, self.r, self.lastcollision) = data
        else:
            self.string = self.fmtstr.pack(
                x, y, z, vx, vy, vz, ax, ay, az, m, r, lastcollision)

            (self.x, self.y, self.z,
                self.vx, self.vy, self.vz,
                self.ax, self.ay, self.az,
                self.m, self.r, self.lastcollision) = (x, y, z,
                                                       vx, vy, vz,
                                                       ax, ay, az,
                                                       m, r, lastcollision)
    def __str__(self):
        return "<rebound.Particle object, x=%f y=%f z=%f vx=%f vy=%f vz=%f m=%f r=%f>"%(self.x,self.y,self.z,self.vx,self.vy,self.vz,self.m,self.r)

class Orbit:
    fmtstr = struct.Struct("dddddddddd")

    def __init__(self, string="",
                 a=0.0, r=0.0, h=0.0, P=0.0, l=0.0, e=0.0, inc=0.0,
                 Omega=0.0, omega=0.0, f=0.0):
        if len(string) > 0:
            self.string = string
            try:
                data = self.fmtstr.unpack(string)
            except Exception, e:
                print string
                print e
            (self.a, self.r, self.h, self.P, self.l, self.e, self.inc,
             self.Omega, self.omega, self.f) = data
        else:
            self.string = self.fmtstr.pack(
                a, r, h, P, l, e, inc,
                Omega, omega, f)

            (self.a, self.r, self.h, self.P, self.l, self.e, self.inc,
             self.Omega, self.omega, self.f) = (a, r, h, P, l, e, inc,
                                                                Omega, omega, f)
    def __str__(self):
        return "<rebound.Orbit object, a=%f r=%f e=%f Omega=%f omega=%f>"%(self.a, self.r, self.e, self.Omega, self.omega)
    def toJson(self):
        return json.dumps({
            'ma': this.f,
            'a': this.a,
            'e': this.e,
            'i': this.inc,
            'w': this.omega,
            'om': this.Omega,
            'P': this.P
        })


TWOPI = 2.*math.pi
def mod2pi(f):
    """Returns the angle f modulo 2 pi."""
    while f < 0.:
        f += TWOPI
    while f > TWOPI:
        f -= TWOPI
    return f

def notNone(a):
    """Returns True if array a contains at least one element that is not None. Returns False otherwise."""
    return a.count(None) != len(a)

def eccentricAnomaly(e, M):
    """Returns the eccentric anomaly given the eccentricity and mean anomaly of a Keplerian orbit.

    Keyword arguments:
    e -- the eccentricity
    M -- the mean anomaly
    """
    E = M if e < 0.8 else math.pi
    
    F = E - e*math.sin(M) - M
    for i in range(100):
        E = E - F/(1.0-e*math.cos(E))
        F = E - e*math.sin(E) - M
        if math.fabs(F)<1e-16:
            break
    E = mod2pi(E)
    return E

k = 0.01720209895       # Gaussian constant 
G = k**2

def kepler_particle(m,
                    primary,    # central body (rebound.Particle object)
                    a,          # semimajor axis
                    anom=0.,  # anomaly
                    e=0.,     # eccentricity
                    omega=0., # argument of pericenter
                    inc=0.,   # inclination
                    Omega=0., # longitude of ascending node
                    MEAN=False):    # mean anomaly
    """Returns a particle structure initialized with the passed set of
        orbital elements. Mass (m), primary and 'a' are required (see Parameters
        below, and any orbital mechanics text, e.g., Murray & Dermott
        Solar System Dynamics for definitions). Other values default to zero.
        All angles should be passed in radians. Units are set by the
        gravitational constant G (default = 1.). If MEAN is set to True, anom is
        taken as the mean anomaly, rather than the true anomaly.
        
        Usage
        _____
        primary = rebound.Particle(m=1.) # particle with unit mass at origin & v=0
        
        # test particle (m=0) with specified elements using mean anomaly
        p = kepler_particle(m=0.,primary=primary,a=2.5, anom=math.pi/2,e=0.3,
        omega=math.pi/6,inc=math.pi/3,Omega=0.,MEAN=True)
        
        # m=0.1 particle on circular orbit math.pi/4 from x axis in xy plane
        p = kepler_particle(0.1,primary,2.5,math.pi/4)
        
        Parameters
        __________
        m       : (float)            Mass of the particle
        primary : (rebound.Particle) Particle structure for the central body
        a       : (float)            Semimajor axis
        anom    : (float)            True anomaly (default).
        Mean anomaly if MEAN is set to True
        e       : (float)            Eccentricity
        omega   : (float)            Argument of pericenter
        inc     : (float)            Inclination (to xy plane)
        Omega   : (float)            Longitude of the ascending node
        MEAN    : (boolean)          If False (default), anom = true anomaly
        If True, anom = mean anomaly
        
        Returns
        _______
        A rebound.Particle structure initialized with the given orbital parameters
        """
    
    if not(0.<=e<1.): raise ValueError('e must be in range [0,1)')
    # not sure if these equations work for parabolic/hyperbolic obits
    if not(0.<=inc<=math.pi): raise ValueError('inc must be in range [0,pi]')
    
    if MEAN is True: # need to calculate f
        E = eccentricAnomaly(e,anom)
        f = mod2pi(2.*math.atan(math.sqrt((1.+ e)/(1. - e))*math.tan(0.5*E)))
    else:
        f = anom
    
    cO = math.cos(Omega)
    sO = math.sin(Omega)
    co = math.cos(omega)
    so = math.sin(omega)
    cf = math.cos(f)
    sf = math.sin(f)
    ci = math.cos(inc)
    si = math.sin(inc)
    
    r = a*(1.-e**2)/(1.+e*cf)
    
    # Murray & Dermott Eq. 2.122
    x  = primary.x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
    y  = primary.y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
    z  = primary.z + r*(so*cf+co*sf)*si
    
    n = math.sqrt(G*(primary.m+m)/(a**3))
    v0 = n*a/math.sqrt(1.-e**2)
    
    # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
    vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
    vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
    vz = primary.vz + v0*((e+cf)*co*si - sf*si*so)
    
    return Particle(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)

sun = Particle(m=1.00000597682, x=-4.06428567034226e-3, y=-6.08813756435987e-3, z=-1.66162304225834e-6, vx=+6.69048890636161e-6, vy=-6.33922479583593e-6, vz=-3.13202145590767e-9)

def convert_to_coord(dic):
    dic['GM'] = dic.get('GM', 0.000001);
    if dic['GM'] == '':
        dic['GM'] = 0
    return kepler_particle(float(dic['GM']), sun, dic['a'],
                           anom=math.radians(dic['ma']), MEAN=True, e=dic['e'],
                           omega=math.radians(dic['w']), inc=math.radians(dic['i']), Omega=math.radians(dic['om'])
                           )


def simulate(rocket, objects):
    particles = []
    particles.append(rocket)

    for i in objects:
        i['r'] = 0.0000001
        i['GM'] = 0.000001
        particles.append(convert_to_coord(i))

    p = subprocess.Popen(["./nbody"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    for i in particles:
        p.stdin.write(i.string)
    p.stdin.close()
    out = p.stdout.read()
    for i in range(len(out)):
        print "thing found!"
        piece = out[i:i+10]
        orbit = Orbit(string=piece)
        print orbit.toJson()
