# basic quaternion algebra
# expects normalized quaternions
# rule-of-thumb: Don't give any rotations more than pi (there is actually no need)
# q=[1,0,0,0] do not rotate the vector

import numpy as np
import scipy.linalg

def qq_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return np.array([w, x, y, z])


def q_conjugate(q):
    w, x, y, z = q
    return np.array([w, -x, -y, -z])


def qv_mult(q1, v):
    q2 = [0.0] + list(v)
    return qq_mult(qq_mult(q1, q2), q_conjugate(q1))[1:]


def axisangle_to_q(v, theta):
    if not np.allclose(v, [0,0,0]):
        v = np.array(v)/np.linalg.norm(v)
    x, y, z = v
    theta /= 2
    w = np.cos(theta)
    x = x * np.sin(theta)
    y = y * np.sin(theta)
    z = z * np.sin(theta)
    return np.array([w, x, y, z])


def q_to_axisangle(q):
    q = np.array(q)
    w, v = q[0], q[1:]
    theta = np.arccos(w) * 2.0
    if not np.allclose(v, [0,0,0]):
        v = np.array(v)/np.linalg.norm(v)
    return v, theta


def q_between_vectors(v1, v2, direction_if_antiparallel=np.array([0,0,0])):
    # q that will rotate v1 to v2
    # if v1 is antiparallel to v2, then choose a simple v to rotate v1
    v1=np.array(v1)/np.linalg.norm(v1)
    v2=np.array(v2)/np.linalg.norm(v2)
    dot = np.dot(v1, v2)
    if np.isclose(dot,-1):
        if np.isclose(np.array(direction_if_antiparallel),np.array([0,0,0])).all():
            v = scipy.linalg.null_space(np.array([v1]))[:,0]
        else: #check if direction_if_antiparallel is actually in the null_space
            null_space = scipy.linalg.null_space(np.array([v1]))
            if np.isclose(np.linalg.norm(direction_if_antiparallel),np.linalg.norm([np.dot(null_space[:,0],direction_if_antiparallel),np.dot(null_space[:,1],direction_if_antiparallel)])):
                v = direction_if_antiparallel
            else:
                raise ValueError('direction_if_antiparallel is not perpendicular to v1')
    else:
        v = np.cross(v1, v2)
    costheta = dot/(np.linalg.norm(v1)*np.linalg.norm(v2))
    if np.isclose(costheta,1):
        costheta= 1
    elif np.isclose(costheta,-1):
        costheta= -1
    theta = np.arccos(costheta)
    q = axisangle_to_q(v, theta)
    
    return q


def angle_v1tov2(v1, v2):
    v1=np.array(v1)/np.linalg.norm(v1)
    v2=np.array(v2)/np.linalg.norm(v2)

    costheta = v1.dot(v2)
    sintheta = np.cross(v1,v2)[2]
    if sintheta >= 0:
        return np.arccos(costheta)
    else:
        return -np.arccos(costheta)
