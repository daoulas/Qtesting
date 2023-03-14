
# This program tests how D-Wave handles the minimization of a square-gradient term defined via a simple forward difference scheme:
# F = (Kappa*D/2)*\Sum_{i=0}^{nodes-1}[(p(i+1)-p(i))/D]^2
# Periodic boundary conditions apply: p(nodes) = p(0) and p(-1) = p(nodes-1)

from dwave.system import DWaveSampler,EmbeddingComposite
  
import dimod, math

import dwave.inspector

Q = {}

V = {}

density_result = {}

# number of variables. We use 3-bit a approximation, so each lattice node has 3 binary variables for our functional. 
# q_k[i] k = 1 .. 3 and i = 0 .. n indicates the node.  In total we have n + 1 nodes. We map q_k[i] on a 1-D matrix
# q[k+i*3]

nodes = 3 

runs = 500

#largest value of p

h = 1

#step of lattice in space 
D = 1

#

u = 5.

Kappa = 5.*0

Coeff = (Kappa/2.)*(1/D)  

po = 1

po2 = po*po

chempot = 0.

# units of p[i]: ~1/Length
# units of Kappa: energy*Length^3
# units of Coeff: energy*Length^2 

print("Coeff = ", Coeff)

b1 = h/7
b2 = 2*h/7
b3 = 4*h/7

#define powers 

b1p2 = math.pow(b1,2)
b1p3 = math.pow(b1,3)
b1p4 = math.pow(b1,4)

b2p2 = math.pow(b2,2)
b2p3 = math.pow(b2,3)
b2p4 = math.pow(b2,4)

b3p2 = math.pow(b3,2)
b3p3 = math.pow(b3,3)
b3p4 = math.pow(b3,4)

#Weight coeffs for constraints

Lan2 = 1.5

for i in range(0,nodes+1):
      Q[1+i*5,2+i*5] = (2*b1p3*b2*u + 3*b1p2*b2p2*u - 3*b1p2*b2*u*po + 2*b1*b2p3*u - 3*b1*b2p2*u*po + b1*b2*po2*u)*D + Lan2*2 + 4*b1*b2*Coeff
      Q[1+i*5,3+i*5] = (2*b1p3*b3*u + 3*b1p2*b3p2*u - 3*b1p2*b3*u*po + 2*b1*b3p3*u - 3*b1*b3p2*u*po + b1*b3*po2*u)*D + Lan2*2 + 4*b1*b3*Coeff
      Q[2+i*5,3+i*5] = (2*b2p3*b3*u + 3*b2p2*b3p2*u - 3*b2p2*b3*u*po + 2*b2*b3p3*u - 3*b2*b3p2*u*po + b2*b3*po2*u)*D + Lan2*2 + 4*b2*b3*Coeff
      Q[1+i*5,4+i*5] = -Lan2*2
      Q[1+i*5,5+i*5] = -Lan2*4
      Q[2+i*5,4+i*5] = -Lan2*2
      Q[2+i*5,5+i*5] = -Lan2*4
      Q[3+i*5,4+i*5] = -Lan2*2
      Q[3+i*5,5+i*5] = -Lan2*4
      Q[4+i*5,5+i*5] = (6*b1p2*b2*b3*u + 6*b1*b2p2*b3*u + 6*b1*b2*b3p2*u - 6*b1*b2*b3*po*u)*D + Lan2*4


for i in range(0,nodes+1):
        V[1+i*5] = (-b1p3*u*po + 0.5*b1p2*u*po2 + 0.5*b1p4*u) + Lan2 - chempot*D*b1 + Coeff*2*b1p2
        V[2+i*5] = (-b2p3*u*po + 0.5*b2p2*u*po2 + 0.5*b2p4*u) + Lan2 - chempot*D*b2 + Coeff*2*b2p2
        V[3+i*5] = (-b3p3*u*po + 0.5*b3p2*u*po2 + 0.5*b3p4*u) + Lan2 - chempot*D*b3 + Coeff*2*b3p2
        V[4+i*5] = Lan2
        V[5+i*5] = 4*Lan2


for i in range(0,nodes): 
     Q[1+i*5,1+(i+1)*5] = -2*b1p2*Coeff
     Q[1+i*5,2+(i+1)*5] = -2*b1*b2*Coeff
     Q[2+i*5,1+(i+1)*5] = -2*b1*b2*Coeff
     Q[1+i*5,3+(i+1)*5] = -2*b1*b3*Coeff
     Q[3+i*5,1+(i+1)*5] = -2*b1*b3*Coeff
     Q[2+i*5,2+(i+1)*5] = -2*b2p2*Coeff
     Q[2+i*5,3+(i+1)*5] = -2*b2*b3*Coeff
     Q[3+i*5,2+(i+1)*5] = -2*b2*b3*Coeff 
     Q[3+i*5,3+(i+1)*5] = -2*b3p2*Coeff

Q[1+0*5,1+nodes*5] = -2*b1p2*Coeff
Q[1+0*5,2+nodes*5] = -2*b1*b2*Coeff
Q[2+0*5,1+nodes*5] = -2*b1*b2*Coeff
Q[1+0*5,3+nodes*5] = -2*b1*b3*Coeff
Q[3+0*5,1+nodes*5] = -2*b1*b3*Coeff
Q[2+0*5,2+nodes*5] = -2*b2p2*Coeff
Q[2+0*5,3+nodes*5] = -2*b2*b3*Coeff
Q[3+0*5,2+nodes*5] = -2*b2*b3*Coeff 
Q[3+0*5,3+nodes*5] = -2*b3p2*Coeff


offset = 0

vartype = dimod.BINARY

model = dimod.BinaryQuadraticModel(V, Q, offset, vartype)

#This is the old way I used:

#sampler = EmbeddingComposite(DWaveSampler())

#New way for calling the sampler

qpu_advantage = DWaveSampler()

sampler = EmbeddingComposite(qpu_advantage)

print("Maximum available annealing time:")
print(qpu_advantage.properties["annealing_time_range"][1]) 

long_time = qpu_advantage.properties["annealing_time_range"][1]*0.4

long_time = 40

sampleset = sampler.sample(model, num_reads = runs, chain_strength = 10, annealing_time=long_time)

#sampleset = sampler.sample(model, num_reads = runs, chain_strength = Kappa)

#ANALYSIS OF SOLUTIONS STARTS HERE 

#This prints out the energies of all solutions 

all_solutions = len(sampleset.record)

print("The total number of solutions found is:", all_solutions) 

#for i in range(0,all_solutions):
#    print(sampleset.record[i][1])

#Establish the smallest energy in the sample of the solutions. 
#This smallest energy can be achieved by different solutions, i.e. several equal minima.

# Set some very high upper boundary for the energy 

Energymax = 100000

numoptsolutions = 0

for i in range(0,all_solutions):
    if(sampleset.record[i][1]<=Energymax):
         Energymax = sampleset.record[i][1]

print("The minimum energy found in these runs is")

print(Energymax)
    
#Output all solutions that correspond to the minimum energy

# Tol is the tolerance for differnces between minimum solutions due to numerics

Tol = 1e-08
for i in range(0,all_solutions):
    if(abs(sampleset.record[i][1] - Energymax)<Tol):
        print("Solution", numoptsolutions)
        numoptsolutions = numoptsolutions + 1

        print("The energy of this solution is", sampleset.record[i][1]) 
        print("Density distribution")
        for j in range(0,nodes+1):    
          density_result[j] = b1*sampleset.record[i][0][j*3+0] + b2*sampleset.record[i][0][j*3+1] + b3*sampleset.record[i][0][j*3+2]
          print(0.5*D + j*D, density_result[j])


print("=======================")

#print(sampleset.record)

#print("=======================")

# Cross-check the previous output. 
# The output below gives one minimum solution that must be part of the previous output

# for i in range(0,nodes+1):
# density_result[i] = b1*sampleset.first[0][i*3+1] + b2*sampleset.first[0][i*3+2] + b3*sampleset.first[0][i*3+3]
# print(0.5*D + i*D, density_result[i])

print("Calculate total mass in box")

Mass = 0.
#Fenergy = 0

for i in range(0,nodes+1):
    Mass = Mass + density_result[i]*D
#    Fenergy = Fenergy +0.5*u*density_result[i]*(density_result[i]-po)*density_result[i]*(density_result[i]-po)*D

print("Mass = ", Mass)

#print("Free energy = ", Fenergy - chempot*Mass)

#  print(i,sampleset.first[0][i*5+4]+sampleset.first[0][i*5+5])
 
print("======") 

print(sampleset.info["timing"])

dwave.inspector.show(sampleset)

#for i in range(1,(nodes+1)*5+1):
#    print(sampleset.first[0][i])
    

#result = sampleset.first.sample

#print(result)



