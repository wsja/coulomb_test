# Try rotations

import numpy as num
import math

pi = num.pi
r2d = 180./pi
d2r = pi/180.

def ned_to_normal_shear_tangential(tensor, strike, dip, rake):
    """
    Rotate symmetric (stress or strain) 9-element tensor to coordinate system
    of the earthquake defined by strike, dip and rake
    tensor: (float) nn, ne, nd, en, ee, ed, nd, ed, dd
    strike: (float) stike of the fault [degree]
    dip: (float) dip of the fault [degree]
    rake: (float) sense of shearing on the fault plane [degree]
    """

    def euler_to_matrix(alpha, beta, gamma):
        '''Given euler angle triplet, create rotation matrix

        Given coordinate system `(x,y,z)` and rotated system `(xs,ys,zs)`
        the line of nodes is the intersection between the `x,y` and the `xs,ys`
        planes.

        :param alpha: is the angle between the `z`-axis and the `zs`-axis [rad]
        :param beta:  is the angle between the `x`-axis and the line of nodes [rad]
        :param gamma: is the angle between the line of nodes and the `xs`-axis
            [rad]

        Usage for moment tensors::

            m_unrot = numpy.matrix([[0,0,-1],[0,0,0],[-1,0,0]])
            euler_to_matrix(dip,strike,-rake, rotmat)
            m = rotmat.T * m_unrot * rotmat

        '''

        ca = math.cos(alpha)
        cb = math.cos(beta)
        cg = math.cos(gamma)
        sa = math.sin(alpha)
        sb = math.sin(beta)
        sg = math.sin(gamma)

        mat = num.matrix([[cb*cg-ca*sb*sg,  sb*cg+ca*cb*sg,  sa*sg],
                          [-cb*sg-ca*sb*cg, -sb*sg+ca*cb*cg, sa*cg],
                          [sa*sb,           -sa*cb,          ca]],
                          dtype=num.float)
        return mat

    tensor = num.matrix(num.reshape(tensor, (3, 3))).A
    print(tensor)

    # north-east-down to normal strike rotation matrix
    # rotmat = euler_to_matrix((90+dip)*d2r, -strike*d2r, 90*d2r).A
    #rotmat = euler_to_matrix(dip*d2r, strike*d2r, -rake*d2r).A
    #rotmat = euler_to_matrix((dip+180)*d2r, strike*d2r, 0).A

    # tensor_nsd = rotmat @ tensor @ rotmat.T  # normal, strike, dip
    # print(tensor_nsd)

    # normal = -tensor_nsd[0,0]
    #slip = -tensor_nsd[0,1]
    #tang = -tensor_nsd[0,2]

    #  Shear stress on the fault plane, in strike and dip (positive up)
    #  coordinates
    # plane_sd = num.matrix([[tensor_nsd[0, 1], 0],
                           # [0, tensor_nsd[0, 2]]]).A
    # print(plane_sd)

    # cr = math.cos(rake*d2r)
    # sr = math.sin(rake*d2r)
    # rotmat_plane = num.matrix([[cr, -sr],
                               # [sr, cr]]).A

    # plane_st = rotmat_plane.T @ plane_sd @ rotmat_plane  # slip, tangential
    # print(plane_st)
    # slip = plane_st[0,0]
    # tang = plane_st[1,1]
    # print(normal, slip, tang)

    rotmat = euler_to_matrix((dip + 180.) * d2r, strike * d2r, rake * d2r).A
    tensor_stn = num.linalg.multi_dot(
        [rotmat, tensor, rotmat.T])  # slip, tangential, normal
    print(tensor_stn)

    slip, tang, normal = tensor_stn[0, 0], tensor_stn[1, 1], -tensor_stn[2, 2]

    return normal, slip, tang
