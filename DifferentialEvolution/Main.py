import numpy as np
from scipy.integrate import odeint
from scipy.optimize import differential_evolution
from math import sin,cos, isnan, isfinite
from re import findall
import time

'''
here we use parameters of a male human with 1.8m height and 80kg weight

loc = 0.43m
lcb = 0.4m
lba = 0.08m

from the book we obtain the following data:

moc_max = 12.13kg
moc_min = 10.53kg
mcb_max = 3.7kg
mcb_min = 3.22kg
mba_max = 1.22
mba_min = 0.97

Joc_max = 0.17kg*m^2
Joc_min = 0.12kg*m^2
Jcb_max = 0.49kg*m^2
Jcb_min = 0.38kg*m^2
Jba_max = 0.0006kg*m^2
Jba_min = 0.0003kg*m^2

koc_max = 0.2m
koc_min = 0.15m
kcb_max = 0.18m
kcb_min = 0.14m
kba_max = 0.04m
kba_min = 0.032m
'''
#system parameters

#exoskeleton part length
loc = 0.43
lcb = 0.4
lba = 0.08

#gravity
g = 9.81


#read the file and init a list of 3 (6) np.array-s
# 3 positional coordinates and 3 velocities (if needed)
#just like the result of scipy.integrate.odeint
class Reader:
    def __init__(self, hip_pos_filename, shin_pos_filename, foot_pos_filename, hip_vel_filename, shin_vel_filename, foot_vel_filename):
        # self.data = [self.read_hip_pos(hip_pos_filename), self.read_shin_pos(shin_pos_filename), self.read_foot_pos(foot_pos_filename),
        #             self.read_hip_vel(hip_vel_filename), self.read_shin_vel(shin_vel_filename), self.read_foot_vel(foot_vel_filename)]
        self.data = solve_equations(moc=22,mcb=7,mba=2,Joc=0.8,Jcb=0.3,Jba=0.0008,koc=0.18,kcb=0.15,kba=0.035)
        pass

    def read_hip_pos(self, hip_pos_filename):
        hip_stream = open(hip_pos_filename, 'r')
        values = []
        for line in hip_stream:
            data = line
            result = findall(r"[-+]?\d*\.\d+|\d+", data)
            values.append(float(result[-1]))
        return np.asarray(values)

    def read_shin_pos(self, shin_pos_filename):
        shin_stream = open(shin_pos_filename, 'r')
        values = []
        for line in shin_stream:
            data = line
            result = findall(r"[-+]?\d*\.\d+|\d+", data)
            values.append(float(result[-1]))
        return np.asarray(values)

    def read_foot_pos(self, foot_pos_filename):
        foot_stream = open(foot_pos_filename, 'r')
        values = []
        for line in foot_stream:
            data = line
            result = findall(r"[-+]?\d*\.\d+|\d+", data)
            values.append(float(result[-1]))
        return np.asarray(values)

    def read_hip_vel(self, hip_vel_filename):
        hip_stream = open(hip_vel_filename, 'r')
        values = []
        for line in hip_stream:
            data = line
            result = findall(r"[-+]?\d*\.\d+|\d+", data)
            values.append(float(result[-1]))
        return np.asarray(values)

    def read_shin_vel(self, shin_vel_filename):
        shin_stream = open(shin_vel_filename, 'r')
        values = []
        for line in shin_stream:
            data = line
            result = findall(r"[-+]?\d*\.\d+|\d+", data)
            values.append(float(result[-1]))
        return np.asarray(values)

    def read_foot_vel(self, foot_vel_filename):
        foot_stream = open(foot_vel_filename, 'r')
        values = []
        for line in foot_stream:
            data = line
            result = findall(r"[-+]?\d*\.\d+|\d+", data)
            values.append(float(result[-1]))
        return np.asarray(values)

# making motion eqns with numeriacl coefficients
# and solving them with Adams and BDF method
def solve_equations(moc, mcb, mba, Joc, Jcb, Jba, koc, kcb, kba):

    def equations(y, t):
        alpha, beta, gamma, alpha1, beta1, gamma1 = y

        global loc
        global lcb
        global lba

        # friction torques in joints
        q = 100 - 50 * alpha1
        u = 15 - 50 * (beta1 - alpha1)
        s = 0.5 - 50 * (gamma1 - beta1)

        #gravity
        global g

        eqn = [alpha1, beta1, gamma1,
               (Jba * (2 * sin(alpha - beta) * loc * mba * beta1 ** 2 *
               lcb ** 3 + (sin(2 * (alpha - beta)) * loc ** 2 * alpha1 ** 2 +
               2 * kba * loc *
               (sin(alpha - gamma) * mba -
               cos(alpha - beta) * sin(beta - gamma)) * gamma1 ** 2 +
               2 * loc * mba * ((-g) * cos(alpha - beta) * sin(beta) +
               g * sin(alpha) * mba +
               (sin(alpha - beta) * kcb * beta1 ** 2 + g * sin(alpha)) *
               mcb) + 2 * mba * (-q + u + g * sin(alpha) * koc * moc)) * lcb ** 2 -
               2 * cos(alpha - beta) *
               loc * (s - u +
               kcb * (-2 * sin(alpha - beta) * loc * alpha1 ** 2 +
               sin(beta - gamma) * kba * gamma1 ** 2 +
               g * sin(beta) * (mba + 1)) * mcb) * lcb +
               2 * cos(alpha - beta) * kcb * loc *
               mcb * (-s + u +
               kcb * (alpha1 ** 2 * sin(alpha - beta) * loc - g * sin(beta)) * mcb) +
               2 * Jcb * (-q + u +
               loc * (sin(alpha - beta) * (lcb + kcb * mcb) * beta1 ** 2 +
               sin(alpha - gamma) * kba * gamma1 ** 2 +
               g * sin(alpha) * (mba + mcb)) +
               g * sin(alpha) * koc * moc)) -
               kba * (2 * cos(beta - gamma) * sin(alpha - beta) * kba ** 2 * lcb ** 2 * loc *
               gamma1 ** 2 +
               s * loc * (lcb * (lcb * (-2 * mba * cos(alpha - gamma) +
               cos(alpha - gamma) + cos(alpha - 2 * beta + gamma)) +
               2 * cos(alpha - beta) * cos(beta - gamma) * kcb * mcb) -
               2 * cos(alpha - gamma) * Jcb) +
               kba * (2 * loc * (cos(beta - gamma) * sin(alpha - gamma) -
               cos(alpha - gamma) * sin(beta - gamma) * mba) * beta1 ** 2 *
               lcb ** 3 +
               (2 * (-q + u + g * sin(alpha) * koc * moc) * cos(beta - gamma) ** 2 +
               loc **
               2 * (sin(2 * (alpha - beta)) + sin(2 * (alpha - gamma)) -
               sin(2 * (alpha - gamma)) * mba) * alpha1 ** 2 +
               loc * (g * mba * (sin(alpha - 2 * beta) +
               sin(alpha - 2 * gamma) +
               2 * cos(alpha - gamma) * sin(gamma) * mba) +
               2 * cos(beta -
               gamma) * (sin(alpha - gamma) * kcb * beta1 ** 2 +
               g * cos(beta - gamma) * sin(alpha)) *
               mcb)) * lcb ** 2 +
               loc * (kcb * (alpha1 **
               2 * (sin(2 * (alpha - beta)) + sin(2 * (alpha - gamma))) *
               loc - 2 * g *
               cos(beta - gamma) * (cos(alpha - gamma) * sin(beta) +
               cos(alpha - beta) * sin(gamma) * mba)) * mcb -
               2 * cos(alpha -
               gamma) * (sin(beta - gamma) * Jcb * beta1 ** 2 + (s - u) *
               cos(beta - gamma))) * lcb +
               2 * cos(alpha - gamma) * Jcb *
               loc * (g * sin(gamma) * mba -
               alpha1 ** 2 * sin(alpha - gamma) * loc)))) /
               ((2 * cos(beta - gamma) ** 2 * Joc * lcb ** 2 +
               loc ** 2 * (2 * Jcb * cos(alpha - gamma) ** 2 +
               lcb * (lcb * ((cos(2 * (alpha - gamma)) +
               cos(2 * (beta - gamma)) + 2) * mba -
               4 * cos(alpha - beta) * cos(alpha - gamma) *
               cos(beta - gamma)) +
               2 * cos(beta - gamma) * (cos(beta - gamma) * lcb -
               2 * cos(alpha - beta) * cos(
               alpha - gamma) * kcb) * mcb))) *
               kba ** 2 +
               2 * Jba * (cos(alpha - beta) ** 2 * lcb ** 2 * loc ** 2 +
               cos(alpha - beta) ** 2 * kcb ** 2 * mcb ** 2 * loc ** 2 -
               (Jcb + lcb * (lcb * mba - 2 * cos(alpha - beta) ** 2 * kcb)) * mcb *
               loc ** 2 - (mba * lcb ** 2 + Jcb) *
               (mba * loc ** 2 + Joc))),





               (Jba * (2 * Joc * (s - u +
               lcb * ((-sin(alpha - beta)) * loc * alpha1 ** 2 +
               sin(beta - gamma) * kba * gamma1 ** 2 + g * sin(beta) * mba) +
               kcb * (g * sin(beta) - alpha1 ** 2 * sin(alpha - beta) * loc) *
               mcb) - loc * (sin(2 * (alpha - beta)) * lcb ** 2 * loc * beta1 ** 2 +
               2 *
               lcb * (sin(alpha - beta) * loc ** 2 * (mba + mcb) * alpha1 ** 2 +
               kba * loc * (cos(alpha - beta) * sin(alpha - gamma) -
               sin(beta - gamma) * (mba + mcb)) *
               gamma1 ** 2 +
               loc * (g * mba * (cos(alpha - beta) * sin(alpha) -
               sin(beta) * mba) + (
               sin(2 * (alpha - beta)) * kcb *
               beta1 ** 2 + g * cos(alpha - beta) * sin(
               alpha) -
               g * sin(beta) * mba) * mcb) +
               cos(alpha - beta) * (-q + u + g * sin(alpha) * koc * moc)) +
               2 * (sin(alpha - beta) * kcb * loc ** 2 * mcb * (mba + mcb) *
               alpha1 ** 2 +
               loc * (mba * (-s + u +
               g * cos(alpha) * sin(alpha - beta) * kcb * mcb) +
               mcb * (-s + u +
               kcb * (cos(alpha - beta) * sin(alpha - gamma) * kba *
               gamma1 ** 2 +
               sin(alpha -
               beta) * (cos(
               alpha - beta) * kcb * beta1 ** 2 +
               g * cos(alpha)) * mcb))) +
               cos(alpha - beta) * kcb *
               mcb * (-q + u + g * sin(alpha) * koc * moc)))) +
               kba * (2 * cos(beta - gamma) * Joc *
               lcb * (s +
               kba * (sin(alpha - gamma) * loc * alpha1 ** 2 +
               sin(beta - gamma) * lcb * beta1 ** 2 - g * sin(gamma) * mba)) +
               loc * (2 * cos(alpha - gamma) * sin(alpha - beta) * kba ** 2 * lcb * loc *
               gamma1 ** 2 + 2 * s * loc *
               (lcb * (cos(beta - gamma) * (mba + mcb) -
               cos(alpha - beta) * cos(alpha - gamma)) -
               cos(alpha - beta) * cos(alpha - gamma) * kcb * mcb) +
               kba * (lcb ** 2 *
               loc * (2 * cos(alpha - gamma) * sin(alpha - 2 * beta + gamma) +
               sin(2 * (beta - gamma)) * (mba + mcb)) * beta1 ** 2 +
               2 * cos(alpha - gamma) *
               loc * ((u - s) * cos(alpha - gamma) -
               kcb * (sin(beta - gamma) * loc * alpha1 ** 2 +
               g * cos(alpha - gamma) * sin(beta) -
               g * cos(alpha - beta) * sin(gamma) * mba) * mcb) +
               2 * lcb * (loc **
               2 * (cos(beta - gamma) *
               sin(alpha - gamma) * (mba + mcb) -
               cos(alpha - gamma) * sin(beta - gamma)) * alpha1 ** 2 +
               loc * (g * mba * (cos(alpha - gamma) *
               sin(alpha - beta + gamma) -
               cos(beta - gamma) * sin(gamma) *
               mba) + (cos(alpha - gamma) *
               sin(
               alpha - 2 * beta + gamma) * kcb * beta1 ** 2 +
               g * cos(beta -
               gamma) * (
               cos(alpha - gamma) * sin(alpha) -
               sin(gamma) * mba)) * mcb) +
               cos(alpha - gamma) *
               cos(beta - gamma) * (-q + u +
               g * sin(alpha) * koc * moc)))))) /
               ((2 * cos(beta - gamma) ** 2 * Joc * lcb ** 2 +
               loc ** 2 * (2 * Jcb * cos(alpha - gamma) ** 2 +
               lcb * (lcb * ((cos(2 * (alpha - gamma)) +
               cos(2 * (beta - gamma)) + 2) * mba -
               4 * cos(alpha - beta) * cos(alpha - gamma) *
               cos(beta - gamma)) +
               2 * cos(beta - gamma) * (cos(beta - gamma) * lcb -
               2 * cos(alpha - beta) * cos(
               alpha - gamma) * kcb) * mcb))) *
               kba ** 2 +
               2 * Jba * (cos(alpha - beta) ** 2 * lcb ** 2 * loc ** 2 +
               cos(alpha - beta) ** 2 * kcb ** 2 * mcb ** 2 * loc ** 2 -
               (Jcb + lcb * (lcb * mba - 2 * cos(alpha - beta) ** 2 * kcb)) * mcb *
               loc ** 2 - (mba * lcb ** 2 + Jcb) *
               (mba * loc ** 2 + Joc))),





               (Joc * lcb * (-2 * s * lcb * mba +
               2 * kba * (lcb * (loc * (cos(beta - gamma) * sin(alpha - beta) -
               sin(alpha - gamma) * mba) * alpha1 ** 2 +
               g * mba * (sin(gamma) * mba - cos(beta - gamma) * sin(beta))) +
               cos(beta - gamma) * (-s + u +
               kcb * (alpha1 ** 2 * sin(alpha - beta) * loc - g * sin(
               beta)) *
               mcb) +
               sin(beta - gamma) * lcb ** 2 * mba * (-beta1 ** 2)) +
               sin(2 * (beta - gamma)) * kba ** 2 * lcb * (-gamma1 ** 2)) -
               Jcb * (2 *
               Joc * (s +
               kba * (sin(alpha - gamma) * loc * alpha1 ** 2 +
               sin(beta - gamma) * lcb * beta1 ** 2 - g * sin(gamma) * mba)) +
               loc * (sin(2 * (alpha - gamma)) * kba ** 2 * loc * gamma1 ** 2 +
               2 * s * loc * (mba + mcb) +
               2 *
               kba * (sin(alpha - gamma) * loc ** 2 * (mba + mcb) * alpha1 ** 2 +
               lcb * loc * (cos(alpha - gamma) * sin(alpha - beta) +
               sin(beta - gamma) * (mba + mcb)) *
               beta1 ** 2 +
               loc * (g * mba * (cos(alpha - gamma) * sin(alpha) -
               sin(gamma) *
               mba) + (cos(alpha -
               gamma) * (sin(alpha - beta) * kcb * beta1 ** 2 +
               g * sin(alpha)) - g * sin(
               gamma) * mba) * mcb) +
               cos(alpha - gamma) * (-q + u + g * sin(alpha) * koc * moc)))) +
               loc * (kba ** 2 * lcb *
               loc * (2 * cos(alpha - beta) * sin(alpha + beta - 2 * gamma) * kcb * mcb +
               lcb * (-2 * cos(alpha - beta) *
               sin(alpha + beta - 2 * gamma) * (mba - 1) -
               sin(2 * (beta - gamma)) * mcb)) * gamma1 ** 2 +
               2 * s *
               loc * (kcb ** 2 * mcb ** 2 * cos(alpha - beta) ** 2 +
               2 * kcb * lcb * mcb * cos(alpha - beta) ** 2 +
               lcb ** 2 * (cos(alpha - beta) ** 2 - mba * (mba + mcb))) +
               kba * (2 *
               loc * (cos(alpha - beta) * sin(alpha - gamma) -
               mba * (cos(alpha - gamma) * sin(alpha - beta) +
               sin(beta - gamma) * (mba + mcb))) * beta1 ** 2 * lcb ** 3 +
               (loc **
               2 * (2 * cos(alpha - beta) *
               sin(beta - gamma) + (-2 * mba * sin(alpha - gamma) +
               sin(alpha - gamma) +
               sin(alpha - 2 * beta + gamma)) * (mba + mcb)) *
               alpha1 ** 2 +
               loc * (2 * g * (mba - 1) *
               mba * (sin(gamma) * mba -
               cos(alpha - beta) * sin(alpha + beta - gamma)) + (2 *
               g * sin(
               gamma) * mba ** 2 -
               (2 * cos(
               alpha - gamma) * sin(
               alpha - beta) * kcb *
               beta1 ** 2 +
               g * (sin(
               2 * alpha - gamma) + sin(
               2 * beta - gamma) +
               2 * sin(
               gamma))) * mba +
               2 * cos(alpha -
               beta) * (2 * sin(
               alpha - gamma) * kcb * beta1 ** 2 +
               g * cos(
               beta - gamma) * sin(alpha))) * mcb) + (-2 * mba *
               cos(
               alpha - gamma) + cos(
               alpha - gamma) +
               cos(
               alpha - 2 * beta + gamma)) * (
               -q + u +
               g * sin(alpha) * koc * moc)) * lcb ** 2 +
               (2 * kcb * loc ** 2 *
                mcb * (2 * cos(alpha - beta) * sin(beta - gamma) +
                cos(beta - gamma) * sin(alpha - beta) * (mba + mcb)) *
                alpha1 ** 2 +
                loc * (2 * kcb * (cos(alpha - beta) * sin(alpha - gamma) * kcb *
                beta1 ** 2 +
                g * cos(alpha) * cos(beta - gamma) *
                sin(alpha - beta)) *
                mcb ** 2 + (2 * (u - s) * cos(beta - gamma) +
                g * kcb * (2 * cos(alpha - beta) * cos(alpha - gamma) *
                sin(beta) + (sin(2 * alpha -
                gamma) - (cos(2 * (alpha - beta)) + 2) *
                sin(gamma)) * mba)) * mcb +
                2 * (s - u) * (cos(alpha - beta) * cos(alpha - gamma) -
                cos(beta - gamma) * mba)) +
                2 * cos(alpha - beta) * cos(beta - gamma) * kcb *
                mcb * (-q + u + g * sin(alpha) * koc * moc)) * lcb +
                2 * cos(alpha - beta) * kcb * loc *
                mcb * ((s - u) * cos(alpha - gamma) +
                kcb * (sin(beta - gamma) * loc * alpha1 ** 2 +
                g * cos(alpha - gamma) * sin(beta) -
                g * cos(alpha - beta) * sin(gamma) * mba) * mcb)))) /
                ((2 * cos(beta - gamma) ** 2 * Joc * lcb ** 2 +
                loc ** 2 * (2 * Jcb * cos(alpha - gamma) ** 2 +
                lcb * (lcb * ((cos(2 * (alpha - gamma)) +
                cos(2 * (beta - gamma)) + 2) * mba -
                4 * cos(alpha - beta) * cos(alpha - gamma) *
                cos(beta - gamma)) +
                2 * cos(beta - gamma) * (cos(beta - gamma) * lcb -
                2 * cos(alpha - beta) * cos(
                alpha - gamma) * kcb) * mcb))) *
                kba ** 2 +
                2 * Jba * (cos(alpha - beta) ** 2 * lcb ** 2 * loc ** 2 +
                cos(alpha - beta) ** 2 * kcb ** 2 * mcb ** 2 * loc ** 2 -
                (Jcb + lcb * (lcb * mba - 2 * cos(alpha - beta) ** 2 * kcb)) * mcb *
                loc ** 2 - (mba * lcb ** 2 + Jcb) *
                (mba * loc ** 2 + Joc)))

               ]

        return eqn

    #time linspace
    t = np.linspace(0, 10, 501)

    #initial conditions
    y0 = [0, 0, 0, 0, 0, 0]

    w = odeint(equations, y0, t)

    #positions
    #alpha
    ans0 = w[:, 0]
    #beta
    ans1 = w[:, 1]
    #gamma
    ans2 = w[:, 2]
    #velocities
    #alpha
    ans3 = w[:, 3]
    #beta
    ans4 = w[:, 4]
    #gamma
    ans5 = w[:, 5]

    return [ans0, ans1, ans2, ans3, ans4, ans5]

#quality function (have we achieved good enough answer?)
def qual_func(x): #x[0] = moc, x[1] = mcb, x[2] = mba, x[3] = Joc, x[4] = Jcb, x[5] = Jba

    global reader

    ans = solve_equations(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8])

    diff = [(reader.data[0]-ans[0])**2, (reader.data[1]-ans[1])**2, (reader.data[2]-ans[2])**2,
            (reader.data[3]-ans[3])**2, (reader.data[4]-ans[4])**2, (reader.data[5]-ans[5])**2]

    d0max=diff[0].max()
    d1max=diff[1].max()
    d2max=diff[2].max()
    d3max=diff[3].max()
    d4max=diff[4].max()
    d5max=diff[5].max()

    result = d0max + d1max + d2max

    #scipy de is dumb and does not understand that NaN is bad.
    #so we give a little help...
    if(isnan(result) or not isfinite(result)):
        result = 1.e+10
    #print(result)
    return result

#differential evolution function
def diff_evolution(bounds, tol, atol):

    #atol gives unpredictable results, so it remains default
    res = differential_evolution(qual_func, bounds, workers=-1, strategy='currenttobest1exp', tol = tol, atol = atol, popsize = 2, mutation = 1.0, recombination= 1)
    return res

#inside name==main it is a LOCAL variable, but i need global
reader = Reader(hip_pos_filename="hip_pos.txt", shin_pos_filename="shin_pos.txt", foot_pos_filename="foot_pos.txt",
                hip_vel_filename="hip_vel.txt", shin_vel_filename="shin_vel.txt", foot_vel_filename="foot_vel.txt")


#enter point
if __name__ == '__main__':

    #search area bounds
    #obtained from the book
    bounds = [(21.5, 23), (6.5, 7.5), (1.5, 2.5),
              (0.2, 0.4), (0.07, 0.09), (0.0007, 0.0009),
              (0.17,0.19), (0.14,0.16), (0.03, 0.04)]

    t = time.time()
    t2 = time.clock()
    #tolerance is for precision tolerance. the lower value - the greater precision. it stands somewhere in [0,1]
    res = diff_evolution(bounds, tol=0.95, atol = 0.6) #t=0.9
    print(time.time() - t)
    print(time.clock()-t2)


    print("result")
    print(res)
    #print(qual_func([21.9,7.2, 2.3, 0.6, 0.35, 0.001,  0.16, 0.14, 0.031]))
    #print(qual_func([22,  7,   2,   0.8, 0.3,  0.0008, 0.18, 0.15, 0.035]))