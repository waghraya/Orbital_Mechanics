import numpy as np

def PQW2ECI(vec_PQW, Ohm, omega, i):
    rotO = np.array([
        [np.cos(-Ohm), np.sin(-Ohm), 0],
        [-np.sin(-Ohm), np.cos(-Ohm), 0],
        [0,0,1]
    ])
    roto = np.array([
        [np.cos(-omega), np.sin(-omega), 0],
        [-np.sin(-omega),np.cos(-omega),0],
        [0,0,1]
    ])
    roti = np.array([
        [1,0,0],
        [0,np.cos(-i),np.sin(-i)],
        [0,-np.sin(-i),np.cos(-i)]
    ])
    vec_ECI = rotO@roti@roto@vec_PQW
    return vec_ECI


def ECI2PQW(vec_ECI,Ohm,omega,i):
    rotOT = np.array([
        [np.cos(-Ohm), np.sin(-Ohm), 0],
        [-np.sin(-Ohm), np.cos(-Ohm), 0],
        [0,0,1]
    ]).T
    print(rotOT.shape)
    rotoT = np.array([
        [np.cos(-omega), np.sin(-omega), 0],
        [-np.sin(-omega),np.cos(-omega),0],
        [0,0,1]
    ]).T
    rotiT = np.array([
        [1,0,0],
        [0,np.cos(-i),np.sin(-i)],
        [0,-np.sin(-i),np.cos(-i)]
    ]).T
    vec_SEZ = rotoT@rotiT@rotOT@vec_ECI
    print('vec_SEZ:', vec_SEZ)
    return vec_SEZ