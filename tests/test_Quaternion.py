"""QUATERNION TEST MODULE"""

import unittest
import numpy as np
import matplotlib.pyplot as mp

from ScanningTools import ScanningTools, Quaternion

def rot_matrix(theta, u=np.array([0, 0, 1])):
    cos = np.cos(theta)
    sin = np.sin(theta)
    cos_1 = 1 - cos
    return np.array([
        [cos + u[0]**2 * cos_1, u[0] * u[1] * cos_1 - u[2] * sin, u[0] * u[2] * cos_1 + u[1] * sin],
        [u[1] * u[0] * cos_1 + u[2] * sin, cos + u[1]**2 * cos_1, u[1] * u[2] * cos_1 - u[0] * sin],
        [u[2] * u[0] * cos_1 - u[1] * sin, u[2] * u[1] * cos_1 + u[0] * sin, cos + u[2]**2 * cos_1]
        ])

def euler_rotation_matrix(alpha, beta, gamma):
    c1, c2, c3 = (np.cos(alpha), np.cos(beta), np.cos(gamma))
    s1, s2, s3 = (np.sin(alpha), np.sin(beta), np.sin(gamma))
    return np.array([
        [c1*c3 - c2*s1*s3, -c1*s3 - c2*c3*s1, s1*s2],
        [c3*s1 + c1*c2*s3, c1*c2*c3 - s1*s3, -c1*s2],
        [s2*s3, c3*s2, c2]])

class TestQuaternions(unittest.TestCase):
    def test_multiplication(self):
        q = Quaternion(1, 2, 3, 4)
        p = Quaternion(5, 6, 7, 8)
        qp = np.array([-0.83045479853739979, 0.16609095970747995, 0.41522739926869989,
                       0.33218191941495989])
        pq = np.array([-0.83045479853739979, 0.27681826617913324, 0.19377278632539327,
                       0.44290922588661319])
        self.assertTrue(np.allclose((q*p)(), qp))
        self.assertTrue(np.allclose((p*q)(), pq))

    def test_rotation(self):
        x_versor = np.array([1, 0, 0])
        rot = Quaternion.get_quaternion_from_rotation(np.pi/2, versor=(0, 0, 1))
        y_versor = rot.rotate_vector_by_quaternion(x_versor)
        minus_y_versor = rot.rotate_vector_by_quaternion(x_versor, direct=False)
        self.assertTrue(np.allclose(y_versor(), np.array([0, 0, 1, 0])))
        self.assertTrue(np.allclose(minus_y_versor(), np.array([0, 0, -1, 0])))
        u = np.array([np.sqrt(1/2), 0, np.sqrt(1/2)])
        matrix = rot_matrix(np.pi, u)
        z = np.dot(matrix, x_versor)
        qrot = Quaternion.get_quaternion_from_rotation(np.pi, versor=u)
        z_versor = qrot.rotate_vector_by_quaternion(x_versor)
        self.assertTrue(np.allclose(z_versor(), Quaternion.get_quaternion_from_vector(z)()))
        
    def test_interpolation(self):
        x_versor = np.array([1, 0, 0])
        y_versor = np.array([0, 1, 0])
        t0 = 0
        t1 = 1
        t = np.linspace(t0, t1, 5)
        interpolated_vectors = Quaternion.get_interpolation(x_versor, y_versor, t, t0, t1)
        # mp.figure()
        # mp.plot(interpolated_vectors[..., 0], interpolated_vectors[..., 1])
        # mp.title('SLERP interpolation')
        # mp.show()   
        self.assertEqual(np.sum(np.sqrt(np.sum(interpolated_vectors**2, axis=1))) /
                         len(interpolated_vectors), 1)

    def test_euler_angles_rotation(self):
        x_versor = np.array([1, 0, 0])
        y_versor = np.array([0, 1, 0])
        alpha, beta, gamma = np.radians((10, 10, 50))
        test_vector = np.array([3, 2, 1])
        euler_matrix = euler_rotation_matrix(alpha, beta, gamma)
        t_vector = np.dot(euler_matrix, test_vector) 
        t_vector /= np.sqrt(np.sum(t_vector**2))
        qeulerrot = Quaternion.get_quaternion_from_euler(alpha, beta, gamma)
        t_quat_vector = qeulerrot.rotate_vector_by_quaternion(test_vector)
        self.assertTrue(np.allclose(t_vector, t_quat_vector.get_versor()))
        qrot2 = Quaternion.get_quaternion_from_euler(np.pi/2, 0, 0)
        y_versor2 = qrot2.rotate_vector_by_quaternion(x_versor)
        self.assertTrue(np.allclose(y_versor, y_versor2.get_versor()))
        
if __name__ == '__main__':
    unittest.main()

