from dwave.system import DWaveSampler,EmbeddingComposite
  
import dimod, math

import dwave.inspector

Q = {}

V = {}

density_result = {}

# Minimize the quadratic free energy F = u*p*(po-p)/2

# number of variables. We use 3-bit a approximation, so we have 3 binary variables for our functional. 
# q_k k = 1,2,3. 

#largest value of p

h = 1

#reference value of p

po = 1

u = 5.

b1 = h/7
b2 = 2*h/7
b3 = 4*h/7

#define powers 

b1p2 = math.pow(b1,2)

b2p2 = math.pow(b2,2)

b3p2 = math.pow(b3,2)


Q[1,2] = -b1*b2*u
Q[1,3] = -b1*b3*u
Q[2,3] = -b2*b3*u

V[1] = -0.5*b1p2*u + 0.5*b1*po*u
V[2] = -0.5*b2p2*u + 0.5*b2*po*u
V[3] = -0.5*b3p2*u + 0.5*b3*po*u

print("Coupling coefficients")

for i in range(1,4):
    for j in range(i+1,4):
        print(i, j, Q[i,j])
 
for i in range(1,4):
    print(i,V[i])

offset = 0

vartype = dimod.BINARY

model = dimod.BinaryQuadraticModel(V, Q, offset, vartype)

sampler = EmbeddingComposite(DWaveSampler())

sampleset = sampler.sample(model, num_reads = 1000)

print(sampleset.variables)

print("=======================")

density_result = b1*sampleset.first[0][1] + b2*sampleset.first[0][2] + b3*sampleset.first[0][3]
print(density_result)

print("Free energy = ", 0.5*u*density_result*(po-density_result))

print("======") 

dwave.inspector.show(sampleset)

#for i in range(1,(nodes+1)*5+1):
#    print(sampleset.first[0][i])
    

#result = sampleset.first.sample

#print(result)



