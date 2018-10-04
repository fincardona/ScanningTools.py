"""
This module implements the quaternions as a Python's class.
"""
import numpy as np

class Quaternion(object):
    """
    This class implements the quaternions.
    """
    def __init__(self, w=None, x=0, y=0, z=0):
        """
        It initializes a quaternion object: q = (w, x, y, z). The quaternion should be normalized. 
        If not, it will be done automatically.

        Parameters
        ----------
        w : float, default = 0
            If None the quaternion will bw initialized with w = 0.
        x : float, default = 0
            The x coordinate.
        y : float, default = 0
            The y coordinate.
        z : float, default = 0
            The z coordinate.
        """
        if x==0 and y==0 and z==0:
            self.x = 0
            self.y = 0
            self.z = 0
            self.w = 1
        elif w==0 or w==None:
            module = np.sqrt(x**2 + y**2 + z**2)
            self.w = 0
            self.x = x / module
            self.y = y / module
            self.z = z / module
        else:
            module = np.sqrt(w**2 + x**2 + y**2 + z**2)
            self.w = w / module
            self.x = x / module
            self.y = y / module
            self.z = z / module
        
    def __str__(self):
        return np.str([self.w, self.x, self.y, self.z])

    def __call__(self):
        return np.array([self.w, self.x, self.y, self.z])

    def __add__(self, q):
        return Quaternion(self.w + q.w, self.x + q.x, self.y + q.y, self.z + q.z)

    def __sub__(self, q):
        return Quaternion(self.w - q.w, self.x - q.x, self.y - q.y, self.z - q.z)

    def __mul__(self, q):
        return Quaternion(
            self.w*q.w - self.x*q.x - self.y*q.y - self.z*q.z,
            self.w*q.x + self.x*q.w + self.y*q.z - self.z*q.y,
            self.w*q.y - self.x*q.z + self.y*q.w + self.z*q.x,
            self.w*q.z + self.x*q.y - self.y*q.x + self.z*q.w)

    def __abs__(self):
        return np.sqrt(self.w*self.w + self.x*self.x + self.y*self.y + self.z*self.z)

    def __neg__(self):
        return Quaternion(-self.w, -self.x, -self.y, -self.z)

    def __invert__(self):
        return self.conjugate()

    def conjugate(self):
        """
        It returns the conjugate of a quaternion.

        Example
        -------
        >>> q = Quaternion(2, 2, 2, 2)
        >>> print q 
        '[2.0, 2.0, 2.0, 2.0]'
        >>> print q.conjugate()
        '[2.0, -2.0, -2.0, -2.0]'
        """
        return Quaternion(self.w, -self.x, -self.y, -self.z)
   
    def get_array(self):
        """
        It returns a numpy array of shape (4,) whose components are the components of the 
        quaternion.
        """
        return np.array([self.w, self.x, self.y, self.z])

    def get_matrix(self):
        """
        It returns a numpy array of shape (4, 4) whose components are related to the components of 
        the quaternion, according to:
        [[w, -x, -y, -z],
         [x,  w, -z,  y],
         [y,  z,  w, -x], 
         [z, -y,  x,  w]]
        """
        return np.array([[self.w, - self.x, - self.y, - self.z],
                         [self.x,   self.w, - self.z,   self.y],
                         [self.y,   self.z,   self.w, - self.x],
                         [self.z, - self.y,   self.x,   self.w]])
    
    def get_versor(self):
        """
        It returns a numpy array of shape (3,) whose components are the spatial coordinates of the 
        quaternion: [x, y, z].
        """
        return self.get_array()[1:]
        
    def rotate_vector_by_quaternion(self, v, direct=True):
        """
        It rotates the vector v according to rotation described by the quaternion. It returns a 
        quaternion whose (x, y, z) components are the spatial coordinates of the rotated vector and
        w = 0.

        Parameters
        ----------
        v      : numpy array of shape (3,)
                 It is the vector that should be rotated.
        direct : bool, default = True
                 If True, the counterclockwise rotation is performed: q * v * q.conjugate. 
                 Otherwise, the clockwise rotation is performed: q.conjugate * v * q.

        Returns
        -------
        out    : quaternion object.
        """
        if direct:
            return (self * Quaternion.get_quaternion_from_vector(v) * self.conjugate())
        else:
            return (self.conjugate() * Quaternion.get_quaternion_from_vector(v) * self)
    
    @staticmethod
    def get_interpolation(v0, v1, t, t0, t1):
        """
        Given two vectors, v0 and v1, and some instants of time t, included between t0 and t1, it 
        returns the quaternion included between the quaternions corrisponding to v0 and v1. The 
        Gram-Schmidt method for the linear interpolation (LERP) is used.
        
        Parameters
        ----------
        v0, v1 : numpy arrays of shape (3,)
                 Are the Cartesian coordinates (x, y, z) of a given pointing.
        t0, t1 : integers, seconds
                 Are two integer instants of time.
        t      : numpy array, seconds
                 They are not integer instants of time included between t0 and t1.

        Returns
        -------
        out    : numpy array of shape (len(t), 3).

        """
        time = (t - t0) / (t1 - t0)
        mul = Quaternion.get_quaternion_from_vector(v0) * Quaternion.get_quaternion_from_vector(v1)
        phi = np.arccos(mul.w)
        return (v0 * (np.sin((1 - time) * phi) / np.sin(phi))[..., None] +
                v1 * (np.sin(time * phi) / np.sin(phi))[..., None])
               
    @staticmethod
    def get_quaternion_from_vector(v):
        """
        It returns the quaternion corresponding to the vector v with w = 0. 

        Parameters
        ----------
        v : numpy array of shape (3,)

        Example
        -------
        >>> v = np.array([0, 0, 1])
        >>> V = Quaternion.get_quaternion_from_vector(v)
        >>> print V
        '[0, 0, 0, 1]'
        """
        return Quaternion(x=v[0], y=v[1], z=v[2])
        
    @staticmethod
    def get_quaternion_from_rotation(theta, versor=np.array([0, 0, 1])):
        """
        It returns the quaternion correspondig to a rotation about a given direction (versor).

        Parameters
        ----------
        theta  : float, radians
                 The rotation angle.
        versor : numpy array of shape (3,)
                 They are the versor's components about which the rotation is performed.
        """
        cos = np.cos(theta / 2)
        sin = np.sin(theta / 2)
        return Quaternion(cos, versor[0] * sin, versor[1] * sin, versor[2] * sin)
    
    @staticmethod
    def get_quaternion_from_euler(alpha=0, beta=0, gamma=0):
        """
        Given two reference systems(x, y, z) and (X, Y, Z) with the following Euler angles: 
        
        alpha, the angle between the x axis and the line of nodes; 
        beta , the angle between the z and Z axes;
        gamma, the angle between the X axis and the line of nodes;
        
        it returns the quaternion correspondig to the rotation of alpha, beta, gamma about 
        respectively the axes (z, x, z). So that, the following rotation is performed:

        Rz(alpha) * Rx(beta) * Rz(gamma)

        Parameters
        ----------
        alpha : float, radians, default = 0
                Angle between the x axis and the line of nodes.
        beta  : float, radians, default = 0
                Angle between the z and Z axes.
        gamma : float, radians, default = 0
                Angle between the X axis and the line of nodes.
        """
        return(
            Quaternion.get_quaternion_from_rotation(alpha, versor=np.array([0, 0, 1])) *
            Quaternion.get_quaternion_from_rotation(beta, versor=np.array([1, 0, 0])) * 
            Quaternion.get_quaternion_from_rotation(gamma, versor=np.array([0, 0, 1])))
