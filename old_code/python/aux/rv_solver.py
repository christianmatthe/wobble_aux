import numpy as np
from scipy.optimize import curve_fit,fsolve

def mean_anomaly(t,T0,P, M0 = 0):
    #M0 in deg, t and T0 in JD, P in JD
    if M0 == 0:
        M = 2*np.pi*((t-T0)%P)/P
    else:
        M = (M0*(np.pi/180) + 2*np.pi*((t-T0))/P)%(2*np.pi) #accounts for startig M0 in fitted model
    return M

def solve_kep_eqn(M,e):
    """ Solve Keplers equation E - e*sin(E) = M for x"""
    # Calculates eccentric anonaly E from mean anomaly M
    #try:
        #M[0]
        #res = np.zeros(M.shape)
        #for i,Mi in enumerate(l):
            #tmp,= fsolve(lambda E: E-e*np.sin(E) - Mi,Mi)
            #res[i] = tmp
    #except IndexError:
        #res, = fsolve(lambda E: E - e*np.sin(E)-M,M)
    res = fsolve(lambda E: E - e*np.sin(E)-M,M)
    return res

def true_anomaly(E,e):
    theta = 2* np.arctan(np.sqrt((1+e)/(1-e)) * np.tan(E/2)) #this form can do a full 2 pi rotation
    return theta

def radial_velocity(t, K, P, e, omega, T0):
    # K in ms^-1
    #P ,t , T0 in days
    #omega in deg
    omega_rad = np.pi * omega/180
    M = mean_anomaly(t, T0, P)
    E = solve_kep_eqn(M,e)
    theta = true_anomaly(E,e)
    v = K * (np.cos(theta + omega_rad) + e * np.cos(omega_rad))
    return v

def radial_velocity_M0(t, K, P, e, omega, T0, M0):
    # K in ms^-1
    #P ,t , T0 in days
    #omega in deg
    #M0 in  deg
    omega_rad = np.pi * omega/180
    M = mean_anomaly(t, T0, P, M0)
    E = solve_kep_eqn(M,e)
    theta = true_anomaly(E,e)
    v = K * (np.cos(theta + omega_rad) + e * np.cos(omega_rad))
    return v
    
if __name__=="__main__":
    test = [radial_velocity(t/10, 17.38, 2.644, 0.152, 325.8, 1552.077) for t in range(15520,15560)] 
    print(test)
