
"""

@author: Laura Malo

Cours: Quantum simulators

PROBLEM SET 3


"""

from __future__ import division
import numpy as np
from itertools import permutations
from math import pi
from numpy import linalg as LA
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------                      
# ----------------------------------------------------------------------------------------      
#                                        FUNCTIONS DEFINITION                          
# ---------------------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------------------  

def calcul_L(state):
    """Computes the total angular momentum of a state"""
    L = np.sum([i*N for i,N in enumerate(state)])
    return L

# -----------------------------------------------------------------------------------
#                                  Creates the base
# -----------------------------------------------------------------------------------
""" The following functions are used to create the Fock basis
"""
def state_generator(state):
  tuple_possible_states = {s for s in permutations(state)}
  #Change tuple states into lists
  possible_states = [list(s) for s in tuple_possible_states]
  return possible_states

# BASIC STATES CREATION

"""Funtion to create all possible states
"""
def basic_states_generator(num_particles, L):
  basic_state = [0]*num_particles
  basic_state[0] = num_particles

  basic_states_list = [basic_state.copy()]
  num_moviments = int(num_particles*(num_particles-1) / 2)
  for n in range(num_moviments):
    basic_state = next_state( basic_state.copy())
   
    basic_states_list.append( basic_state)
    
  states_list = delete_duplicated( basic_states_list)
  
  #Give length associated to L
  #L = num_particles*( num_particles - 1)
  if L > num_particles:
    zeros_extres = [0]*(L-num_particles + 1)
  
    for state in states_list:
      state.extend(zeros_extres)
  
  return states_list
  
  
def next_state(state):
  index = 0
  STATUS = True
  
  while STATUS == True and index < len(state):
    
    if state[index] > 1:
      state[index] -= 1
      state[index+1] += 1
      STATUS = False
      
    else:
      index += 1

  return state
  
def delete_duplicated(states):
  for i, state in enumerate(states):
    states[i] = sorted(state, reverse = True)
    
  states_set = set(map(tuple,states))  #need to convert the inner lists to tuples so they are hashable
  states =[list(s) for s in states_set]   
  return  states

#BASIS GENERATION

def basis_generator(num_particles, L):
  """Function to generate the basis asociated to the number of particles we have"""
  #L = num_particles*(num_particles-1)
  
  basic_states = basic_states_generator(num_particles, L)
  
  #Generation of all possible states
  possible_states = [state_generator(s) for s in basic_states]
  
  #Validation of the states
  states = []
  for p_states in possible_states:
    for s in p_states:
      if calcul_L(s) <= L:
        states.append(s)
      
  #Sort the states depending on L
  sorted_states= sorted(states, key = calcul_L)
  
  return sorted_states

# -----------------------------------------------------------------------------------
#                          creates the blocs of momentum L  
# -----------------------------------------------------------------------------------

def bloc_L(basis,L):
    """ Given the full Fock basis and an angular momentum L, it returns the Fock basis corresponding to
    the subspace of the Fock space with L.
    """
    bloc = []
    for b in basis:
        if (calcul_L(b) == L):
            bloc.append(b)
    return bloc

# -----------------------------------------------------------------------------------
#                          single particle hamiltonian
# -----------------------------------------------------------------------------------
    
def Hsp_element(state,omega):
    """ Computes the single particle energy of a state for a given value of the ratio Omega/omega"""
    res= 0.0
    for i in range(len(state)):
        e = 1+i*(1-omega)
        res=res+state[i]*e
    return res

def H1(bloc,omega):
    """Constructs the matrix of the single particle term of the Hamiltonian"""
    mat=np.zeros([len(bloc),len(bloc)])
    for i in range(len(bloc)):
        mat[i,i]=Hsp_element(bloc[i],omega)
    return mat

# -----------------------------------------------------------------------------------
#                          interaction hamiltonian
# -----------------------------------------------------------------------------------
        
def anihilation(i,state_in):
    """Anihilation operator over a state_in for the position i"""
    if not (state_in[i] == 0):
        coef = np.sqrt(state_in[i])
        state_out=state_in.copy()
        state_out[i]=state_out[i]-1
        stop = False
        return state_out,coef,stop
    else:
        #print('This state cant be lowered at', i,'!', )
        stop = True 
        state_out= []
        coef=0
        return state_out,coef,stop
    
def creation(i,state_in):
    """Creation operator over a state state_in for the position i"""
    coef = np.sqrt(state_in[i]+1)
    state_out=state_in.copy()
    state_out[i] = state_out[i]+1
    return state_out,coef  

def factorial(k):
    """Computes the factorial of a number: k!"""
    i=k
    res=1
    while i>0:
        res = res*i
        i=i-1
    return res
        

def I(k,l,p,q):
    """Computes the integral I for four indexes"""
    if (k+l-p-q==0):
        return (1/(2*pi))*(1/(2**(k+l)))*(1/np.sqrt(factorial(k)*factorial(l)*
            factorial(p)*factorial(q))) *factorial(k+l)        
    else:
        return 0

def Hint(state_ini):
    """Computes how the interaction hamiltonian acts over a given state"""
    states = []
    coefs = []
    for k in range(len(state_ini)):
        for l in range(len(state_ini)):
            for p in range(len(state_ini)):
                for q in range(len(state_ini)):
                    if not (I(k,l,p,q)==0):
                        state1,coef1,stop1 = anihilation(p,state_ini)
                        if not stop1: 
                            state2,coef2,stop2 = anihilation(q,state1)
                            if not stop2: 
                                state3,coef3 = creation(l,state2)
                                state4,coef4 = creation(k,state3)
                                states.append(state4)
                                coefs.append(I(k,l,p,q)*coef1*coef2*coef3*coef4)
                    
    return states,coefs

def construct_H(state_bra,state_ket):
    """Constructs the matrix elements of the interaction term by computing the braket
    over all of the states of the Fock basis"""
    states,coefs=Hint(state_ket)
    res = 0
    for i in range(len(states)):
        if (state_bra==states[i]):
            res = res+ coefs[i]
    return 0.5*res 

def H2(bloc):
    """Returns the matrix corresponding to the interaction term"""
    mat=np.zeros([len(bloc),len(bloc)])
    for i in range(len(bloc)):
        for j in range(len(bloc)):
            mat[i,j]=construct_H(bloc[i],bloc[j])
    return mat
            

# -----------------------------------------------------------------------------------
#                          full hamiltonian 
# -----------------------------------------------------------------------------------  
    
def hamiltonian_matrix(bloc,omega):
    """Constructs the full Hamiltonina by adding the single particle and the interaction
    term"""
    H_1 = H1(bloc, omega)
    H_2 = H2(bloc)
    return H_1+H_2

# -----------------------------------------------------------------------------------
#                                Ground state 
# -----------------------------------------------------------------------------------  

def gs_ene(bloc, omega):
    """Computes the energy of ground state for a subspace of angular momentun given by 
    the bloc we are considering"""
    hamiltonian = hamiltonian_matrix(bloc,omega)
    w,v = LA.eig(hamiltonian)
    return np.min(w)

def gs_omega(blocs,omega):
    """Computes the ground state for all the angular momentums for a given omega. This is
    used to compute the Yrast line."""
    gs_results = []
    for i in range(len(blocs)):
        gs_results.append(gs_ene(bloc = blocs[i],omega = omega))
    return gs_results

def gs(bloc,omega):
    """Computes the ground state. Retuns the energy of the ground state and the
    Fock coefficients of the gs."""
    hamiltonian = hamiltonian_matrix(bloc,omega)
    w, v = LA.eig(hamiltonian)
    gs_ene = np.min(w)
    for i in range(len(w)):
        if w[i]==gs_ene:
            gs_state = v[i]
    return gs_ene, gs_state   
    
# -----------------------------------------------------------------------------------
#                                     Density
# -----------------------------------------------------------------------------------  
    
def psi(z,k):
    """Returns the wavefunction for a given position and angular momentum"""
    return (1/np.sqrt(pi*factorial(k)))*z**k * np.exp(-(np.abs(z)**2)/2)

def braket(i,j,bra,ket):
    """Computes the bracket that appears in the expression of the density"""
    res = 0
    state1,coef1,stop1 = anihilation(j,ket)
    if not stop1: 
        state2,coef2 = creation(i,state1)
        if np.all(state2 == bra):
            res = coef1*coef2
    return res
        
def density(bloc,omega,z):
    """Computes the density for a given position z"""
    _,coefs_fock = gs(bloc, omega)
    den = 0
    states_fock = bloc.copy()
    for x in range(len(coefs_fock)):
        for y in range(len(coefs_fock)):
                        for i in range(len(states_fock[x])):
                            for j in range(len(states_fock[x])):
                                add = coefs_fock[x]*coefs_fock[y]*psi(z,i)*psi(z,j)*braket(i,j,states_fock[x],states_fock[y])
                                den = den + add
    return den

                      
def den_z(bloc, omega, z_values):
    """Density for a range of positions of z"""
    den_values = []
    for z in z_values:
        den_values.append(density(bloc,omega,z))
    return den_values

# -----------------------------------------------------------------------------------
#                                       Pair correlation
# -----------------------------------------------------------------------------------  
    
def braket_pair(i,j,k,l,bra,ket):
    """Computes the braket that appears in the expression of the pair correlation"""
    res = 0
    state1,coef1,stop1 = anihilation(l,ket)
    if not stop1:
        state2,coef2,stop2 = anihilation(k,state1)
        if not stop2:
            state3,coef3 = creation(j,state2)
            state4,coef4 = creation(i, state3)
            if np.all(state4 == bra):
                res = coef1*coef2*coef3*coef4
    return res 


def pair_correlation(bloc,omega,z1,z2):
    """Computes the pair correlation of a block of angular momentum for the postion
    z1 and z2"""
    _,coefs_fock = gs(bloc, omega)
    pair = 0
    states_fock = bloc.copy()
    for x in range(len(coefs_fock)):
        for y in range(len(coefs_fock)):
            for i in range(len(states_fock[0])):
                for j in range(len(states_fock[0])):
                    for k in range(len(states_fock[0])):
                        for l in range(len(states_fock[0])):
                            add = coefs_fock[x]*coefs_fock[y]*psi(z1,i)*psi(z2,j)*psi(z1,k)*psi(z2,l)*braket_pair(i,j,k,l,states_fock[x],states_fock[y])
                            pair = pair + add
    return pair

def paircorrelations_z(bloc, omega, z1_values, z2_values):
    """Pair correlation for a range of values of z1 and z2"""
    mat = np.zeros([len(z1_values),len(z2_values)])
    for i in range(len(z1_values)):
        for j in range(len(z2_values)):
            mat[i,j] =  pair_correlation(bloc = bloc ,omega = omega,z1 = z1_values[i],z2 = z2_values[j]) 
    return mat
                           
# ----------------------------------------------------------------------------------------                      
# ----------------------------------------------------------------------------------------      
#                                        MAIN CODE                          
# ---------------------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------------------  

# Fock basis             
B =  basis_generator(num_particles = 3,L = 6) #all states

# boxes of fixed angular momentum       
B0 = bloc_L(basis = B, L= 0)
B1 = bloc_L(basis = B, L= 1)
B2 = bloc_L(basis = B, L= 2)
B3 = bloc_L(basis = B, L= 3)
B4 = bloc_L(basis = B, L= 4) 
B5 = bloc_L(basis = B, L= 5)
B6 = bloc_L(basis = B, L= 6) 

#crate lists with the basis of each bloc and its corresponding angular momentum
blocs = np.array([B0,B1,B2,B3,B4,B5,B6])
L = np.array([0,1,2,3,4,5,6])

# ----------------------------------------------------------------------------------------                      
# ----------------------------------------------------------------------------------------      
#                                        PLOTS                        
# ---------------------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------------------  
"""
This plot computes the Yrast line for the different values of Omega/omega
"""

plt.figure()
plt.plot(L,gs_omega(blocs = blocs, omega = 0.70),linewidth=1,linestyle='dashed',label='$\Omega$/$\omega$=0.70')
plt.plot(L,gs_omega(blocs = blocs, omega = 0.75),linewidth=1,linestyle='dashed',label='$\Omega$/$\omega$=0.75')
plt.plot(L,gs_omega(blocs = blocs, omega = 0.80),linewidth=1,linestyle='dashed',label='$\Omega$/$\omega$=0.80')
plt.plot(L,gs_omega(blocs = blocs, omega = 0.85),linewidth=1,linestyle='dashed',label='$\Omega$/$\omega$=0.85')
plt.plot(L,gs_omega(blocs = blocs, omega = 0.90),linewidth=1,linestyle='dashed',label='$\Omega$/$\omega$=0.90')
plt.plot(L,gs_omega(blocs = blocs, omega = 0.95),linewidth=1,linestyle='dashed',label='$\Omega$/$\omega$=0.95')
plt.plot(L,gs_omega(blocs = blocs, omega = 1.0),linewidth=1,linestyle='dashed',label='$\Omega$/$\omega$=1.0')
plt.xlabel('Angular momentum, L')
plt.legend(loc="best")
plt.ylabel('E = <GS|H|GS>')
plt.xlim(0,6)
plt.ylim(3,5) 

# ----------------------------------------------------------------------------------------  

"""
This plot shows the density profile for a particular value of the ratio Omega/omega (the value is 
specified for each case) for a range of position going from 0 up to 4. 
"""
z_values = np.linspace(0.0,4.0, num = 50)
plt.figure()
plt.plot(z_values,den_z(bloc = B0, omega =0.80, z_values = z_values), linewidth = 0.8, label = 'L=0')
plt.plot(z_values,den_z(bloc = B1, omega =0.80, z_values = z_values), linewidth = 0.8, label = 'L=1')
plt.plot(z_values,den_z(bloc = B2, omega =0.90, z_values = z_values), linewidth = 0.8, label = 'L=2')
plt.plot(z_values,den_z(bloc = B3, omega =0.90, z_values = z_values), linewidth = 0.8, label = 'L=3')
plt.plot(z_values,den_z(bloc = B4, omega =0.90, z_values = z_values), linewidth = 0.8, label = 'L=4')
plt.plot(z_values,den_z(bloc = B5, omega =0.98, z_values = z_values), linewidth = 0.8, label = 'L=5')
plt.plot(z_values,den_z(bloc = B6, omega =0.98, z_values = z_values), linewidth = 0.8, label = 'L=6')
plt.xlim(0,4)
plt.legend(loc = 'best')
plt.xlabel('z')
plt.ylabel('Density')


# ---------------------------------------------------------------------------------------- 

"""
This plot shows the pair correlation for different angular momentum. The value of Omgega/omega is the same 
in each case =0.8. As the computation time is very high for this plot, it is commented. To obtain the plot, 
remove the semmicolons before and after the lines of the code.
"""

"""
z_values = np.linspace(-3.0,3.0, num = 40)
plt.figure()
plt.subplot(2,3,1)
plt.title('L=0')
plt.imshow(paircorrelations_z(bloc = B0, omega = 0.8, z1_values = z_values, z2_values = z_values))
plt.axis('off')
plt.colorbar()
plt.subplot(2,3,2)
plt.title('L=1')
plt.imshow(paircorrelations_z(bloc = B1, omega = 0.8, z1_values = z_values, z2_values = z_values))
plt.axis('off')
plt.colorbar()
plt.subplot(2,3,3)
plt.title('L=2')
plt.imshow(paircorrelations_z(bloc = B2, omega = 0.8, z1_values = z_values, z2_values = z_values))
plt.axis('off')
plt.colorbar()
plt.subplot(2,3,4)
plt.title('L=3')
plt.imshow(paircorrelations_z(bloc = B3, omega = 0.8, z1_values = z_values, z2_values = z_values))
plt.axis('off')
plt.colorbar()
plt.subplot(2,3,5)
plt.title('L=4')
plt.imshow(paircorrelations_z(bloc = B4, omega = 0.8, z1_values = z_values, z2_values = z_values))
plt.axis('off')
plt.colorbar()
plt.subplot(2,3,6)
plt.title('L=6')
plt.imshow(paircorrelations_z(bloc = B6, omega = 0.8, z1_values = z_values, z2_values = z_values))
plt.axis('off')
plt.colorbar()"""