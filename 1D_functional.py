from dwave.system import DWaveSampler,EmbeddingComposite
  
import dimod, math

# import greedy

import dwave.inspector

Q = {}

V = {}

density_result = {}

# number of variables. We use 3-bit a approximation, so each lattice node has 5 binary variables for our functional. 
# q_k[i] k = 1 .. 5 and i = 0 .. n indicates the node.  In total we have n + 1 nodes. We map q_k[i] on a 1-D matrix
# q[k+i*5]

nodes = 0

#largest value of p

h = 1

#reference value of p

po = 1

#chemical potential

chempot = 0 

#step of lattice in space 
D = 1

u = 5.

D2 = D*D 
#
po2 = po*po

b1 = h/7
b2 = 2*h/7
b3 = 4*h/7

#Weight coeffs for constraints
Lan2 = 3.5

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

for i in range(0,nodes+1):
      Q[1+i*5,2+i*5] = (2*b1p3*b2*u + 3*b1p2*b2p2*u - 3*b1p2*b2*u*po + 2*b1*b2p3*u - 3*b1*b2p2*u*po + b1*b2*po2*u)*D + Lan2*2
      Q[1+i*5,3+i*5] = (2*b1p3*b3*u + 3*b1p2*b3p2*u - 3*b1p2*b3*u*po + 2*b1*b3p3*u - 3*b1*b3p2*u*po + b1*b3*po2*u)*D + Lan2*2
      Q[2+i*5,3+i*5] = (2*b2p3*b3*u + 3*b2p2*b3p2*u - 3*b2p2*b3*u*po + 2*b2*b3p3*u - 3*b2*b3p2*u*po + b2*b3*po2*u)*D + Lan2*2
      Q[1+i*5,4+i*5] = -Lan2*2
      Q[1+i*5,5+i*5] = -Lan2*4
      Q[2+i*5,4+i*5] = -Lan2*2
      Q[2+i*5,5+i*5] = -Lan2*4
      Q[3+i*5,4+i*5] = -Lan2*2
      Q[3+i*5,5+i*5] = -Lan2*4
      Q[4+i*5,5+i*5] = (6*b1p2*b2*b3*u + 6*b1*b2p2*b3*u + 6*b1*b2*b3p2*u - 6*b1*b2*b3*po*u)*D + Lan2*4


for i in range(0,nodes+1):
        V[1+i*5] = (-b1p3*u*po + 0.5*b1p2*u*po2 + 0.5*b1p4*u) + Lan2 - chempot*D*b1
        V[2+i*5] = (-b2p3*u*po + 0.5*b2p2*u*po2 + 0.5*b2p4*u) + Lan2 - chempot*D*b2
        V[3+i*5] = (-b3p3*u*po + 0.5*b3p2*u*po2 + 0.5*b3p4*u) + Lan2 - chempot*D*b3
        V[4+i*5] = Lan2
        V[5+i*5] = 4*Lan2

print("Coupling coefficients")

for i in range(1,6):
    for j in range(i+1,6):
        print(i, j, Q[i,j])
 
for i in range(1,6):
    print(i,V[i])

offset = 0

vartype = dimod.BINARY

model = dimod.BinaryQuadraticModel(V, Q, offset, vartype)

sampler = EmbeddingComposite(DWaveSampler())

#sampler = greedy.SteepestDescentSolver()

sampleset = sampler.sample(model, num_reads = 3000, chain_strength = 15)

print(sampleset.variables)

print("=======================")

#total number of variables (nodes+1)*5

for i in range(0,nodes+1):
    density_result[i] = b1*sampleset.first[0][i*5+1] + b2*sampleset.first[0][i*5+2] + b3*sampleset.first[0][i*5+3]
    print(0.5*D + i*D, density_result[i])

print("Calculate total mass in box")

Mass = 0.
Fenergy = 0

for i in range(0,nodes+1):
    Mass = Mass + density_result[i]*D
    Fenergy = Fenergy +0.5*u*density_result[i]*(density_result[i]-po)*density_result[i]*(density_result[i]-po)*D

print("Mass = ", Mass)

print("Free energy = ", Fenergy - chempot*Mass)

print("Testing constraint between q1,q2,q3 and q4,q5")

for i in range(0,nodes+1):
    print("node:", i, "constraint value:", sampleset.first[0][i*5+1]+sampleset.first[0][i*5+2]+sampleset.first[0][i*5+3]-
          sampleset.first[0][i*5+4]-2*sampleset.first[0][i*5+5])
        
#  print(i,sampleset.first[0][i*5+4]+sampleset.first[0][i*5+5])
 
print("======") 

dwave.inspector.show(sampleset)

#for i in range(1,(nodes+1)*5+1):
#    print(sampleset.first[0][i])
    

#result = sampleset.first.sample

#print(result)



