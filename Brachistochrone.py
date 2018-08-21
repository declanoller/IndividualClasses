from random import shuffle,randint,sample,random
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt,sin,cos
from scipy.optimize import fsolve

class Brachistochrone:

    def __init__(self,N=20,height=1.0):
        #height is as a ratio of the width.
        #N is the number of *segments*, that is, "blocks", not points (which will be N+1).
        #The state is the height of each point, so will be length N+1.
        #To make things easy, I'll also have an xpos list.
        self.N_segments = N
        self.N_pts = N+1
        self.width = 1.0
        self.height = height
        self.state = [self.height]
        self.xpos = [0]
        self.delta_x = self.width/self.N_segments

        for i in range(self.N_segments-1):
            self.xpos.append((i+1)*self.delta_x)
            #self.state.append(random()*self.height)
            frac = 1
            self.state.append(-self.height/frac + random()*(self.height - (-self.height/frac)))

        self.xpos.append(1.0)
        self.state.append(0)

        #Maybe at some point, try something where it starts with very few points, solves as good
        #as it can with them, and then doubles the number of points, in between

        self.sol = None

        #This is just assuming it dropping vertically and then going horizonally at that speed; it doesn't have to be perfect.
        g = 9.8
        t1 = sqrt(2*self.height/g)
        v1 = g*t1
        t2 = self.width/v1
        self.max_FF = 40*(t1 + t2)

        self.getBrachistochroneSol()

    def copyState(self,other_state):
        self.state = other_state.state

    def getBrachistochroneSol(self):
        w = self.width
        h = self.height

        #This is all solved with the assumption that the starting point,
        #where the bead is dropped, is (0,0), meaning that the ending point is
        #(w,-h). See https://math.stackexchange.com/questions/889187/finding-the-equation-for-a-inverted-cycloid-given-two-points
        #Importantly, this means that what you'll inevitably get will be in that coord. system.
        #So, to match it up with the coords we've been using (dropped at (0,h), ending at (w,0)), simply add h to y in the end.

        f_t = lambda t: np.cos(t)-1+ (-h/w)*(np.sin(t)-t)
        t = fsolve(f_t,3.14)[0]

        a = w/(t-sin(t))

        '''print('a:',a)
        print('t:',t)'''

        self.t_range = np.linspace(0,t,self.N_pts)

        self.x = lambda t: a*(t-np.sin(t))
        self.y = lambda t: h + a*(np.cos(t)-1)

        self.sol = (self.t_range,self.x,self.y)

        sol_numeric_y = []

        for x_pt in self.xpos:
            f = lambda t: self.x(t)-x_pt
            tval = fsolve(f,3.14)[0]
            sol_numeric_y.append(self.y(tval))

        temp_state = self.state
        self.state = sol_numeric_y
        self.sol_numeric_y = sol_numeric_y

        self.t_ideal = self.fitnessFunction()
        #print('theoretical best time:',self.t_ideal)

        self.state = temp_state

    def isSameState(self,other_state):

        #return(False)
        max_diff = ((self.N_pts-2)*self.height)**2
        diff = np.array(self.state) - np.array(other_state.state)
        abs_diff = sum(diff**2)
        thresh = 1*10.0**(-5)
        if abs_diff/max_diff < thresh:
            return(True)
        else:
            return(False)

    def solFound(self):

        max_diff = ((self.N_pts-2)*self.height)**2
        diff = np.array(self.state) - np.array(self.sol_numeric_y)
        abs_diff = sum(diff**2)
        thresh = 1*10.0**(-5)
        if abs_diff/max_diff < thresh:
            print('abs_diff/max_diff:',abs_diff/max_diff)
            return(True)
        else:
            return(False)

        '''if abs(self.fitnessFunction() - self.t_ideal)/self.t_ideal < 0.005:
            return(True)
        else:
            return(False)'''


    def plotState(self,plot_axis=None,show=False):

        if plot_axis is None:
            ax = plt.gca()
        else:
            ax = plot_axis

        #ax.clear()

        if self.sol is not None:
            t = self.sol[0]
            x = self.sol[1]
            y = self.sol[2]
            ax.plot(x(t),y(t),'-',color='cornflowerblue')
            ax.text(.8*self.width,.9*self.height,'ideal: {:.3f}'.format(self.t_ideal))


        ax.text(.8*self.width,.8*self.height,'actual: {:.3f}'.format(self.fitnessFunction()))
        #ax.plot(self.xpos,self.state,'o-',color='darkred')
        ax.plot(self.xpos,self.state,'o-')

        if show:
            plt.show()
        #return(plt)




    def mutate(self):
        #inclusive,inclusive



        index = randint(1,self.N_segments-1)
        sway = self.height/20.0
        new_height = self.state[index] + ((-sway) + random()*2*sway)
        #new_height = ((1.0-sway) + random()*2*sway)*self.state[index]
        if new_height < self.height:
            self.state[index] = new_height


        '''index = randint(1,self.N_segments-1)
        sway = .2
        new_height = ((1.0-sway) + random()*2*sway)*self.state[index]
        if new_height < self.height:
            self.state[index] = new_height'''
        '''index = randint(1,self.N_segments-1)
        self.state[index] = random()*self.height'''
        '''M = 3
        index = randint(1,self.N_segments-1-M)
        rand = random()*self.height
        self.state[index:(index+M)] = [rand+i*.0001 for i in range(M)]'''
        '''N_switch = randint(1,self.N_pts-3)
        switch_indices = sample(list(range(1,self.N_pts-2)),N_switch)
        for index in switch_indices:
            self.state[index] = random()*self.height'''


    def fitnessFunction(self):

        g = 9.8

        #So if the next point is lower than the previous one, d will be *positive* (i.e., the y axis is down, opposite with the plot axis.)
        d = -np.array([self.state[i+1] - self.state[i] for i in range(self.N_segments)])

        #Be careful with signs and indices!
        v = sqrt(2*g)*np.sqrt([0] + [sum(d[:(i+1)]) for i in range(len(d))])

        if np.isnan(v).any():
            print('\n\nbad v:',v)
            print('\nbad d sum:',[sum(d[:(i+1)]) for i in range(len(d))])
            print('\nstate',self.state)
            plt.savefig('test_bad_np.png')
            exit(0)


        #v = np.sqrt([0] + [sum(d[:(i+1)]) for i in range(len(d))])
        v = v[:-1]
        t = (np.sqrt(v**2 + 2*g*d) - v)/(g*d/np.sqrt(d**2 + self.delta_x**2))
        '''print('\n\n')
        print('state',self.state)
        print('d',d)
        print([sum(d[:(i+1)]) for i in range(len(d))])
        print('v',v)
        print('t',t)'''


        return(sum(t))





    #*************************************** GA stuff

    def mate(self,other_individ):


        #return(self.mateRandomIndices(other_individ))
        return(self.mateCrossover(other_individ))


    def mateRandomIndices(self,other_individ):

        newindivid_1 = deepcopy(self)
        newindivid_2 = deepcopy(other_individ)
        #inclusive, inclusive
        N_switch = randint(1,self.N_segments-1)
        switch_indices = sample(list(range(1,self.N_segments)),N_switch)

        for index in switch_indices:
            temp = newindivid_1.state[index]
            newindivid_1.state[index] = newindivid_2.state[index]
            newindivid_2.state[index] = temp

        return(newindivid_1,newindivid_2)



    def mateCrossover(self,other_individ):

        newindivid_1 = deepcopy(self)
        newindivid_2 = deepcopy(other_individ)
        #inclusive, inclusive
        #index = randint(2,self.N_pts-2)

        r1 = randint(1,self.N_pts-3)
        r2 = randint(r1+1,self.N_pts-2)

        temp = newindivid_1.state[r1:r2]
        newindivid_1.state[r1:r2] = newindivid_2.state[r1:r2]
        newindivid_2.state[r1:r2] = temp

        '''temp = newindivid_1.state[:index]
        newindivid_1.state[:index] = newindivid_2.state[:index]
        newindivid_2.state[:index] = temp'''

        return(newindivid_1,newindivid_2)


    def mateAvg(self,other_individ):

        newindivid_1 = deepcopy(self)
        newindivid_2 = deepcopy(other_individ)

        newindivid_1.state = ((np.array(newindivid_1.state) + np.array(newindivid_2.state))/2).tolist()


        return(newindivid_1,newindivid_1)

#
