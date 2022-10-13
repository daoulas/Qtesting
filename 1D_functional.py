from dwave.system import DWaveSampler,EmbeddingComposite
  
import dimod, math

Q = {}

V = {}

# number of variables. We use 3-bit a approximation, so each lattice node has 5 binary variables for our functional. 
# q_k[i] k = 1 .. 5 and i = 0 .. n indicates the node.  In total we have n + 1 nodes. We map q_k[i] on a 1-D matrix
# q[k+i*5]

nodes = 5 

h = 1
D = 1

u = 5
M = 3

D2 = D*D 

b1 = h/7
b2 = 2*h/7
b3 = 4*h/7

Lan1 = 100

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
      Q[1+i*5,2+i*5] = (2*b1p3*b2*u + 3*b1p2*b2p2*u - 3*b1p2*b2*u + 2*b1*b2p3*u - 3*b1*b2p2*u + b1*b2*u)*D + Lan1*2*D2*b1*b2
      Q[1+i*5,3+i*5] = (2*b1p3*b3*u + 3*b1p2*b3p2*u - 3*b1p2*b3*u + 2*b1*b3p3*u - 3*b1*b3p2*u + b1*b3*u)*D + Lan1*2*D2*b1*b3
      Q[2+i*5,3+i*5] = (2*b2p3*b3*u + 3*b2p2*b3p2*u - 3*b2p2*b3*u + 2*b2*b3p3*u - 3*b2*b3p2*u + b2*b3*u)*D + Lan1*2*D2*b2*b3
      Q[4+i*5,5+i*5] = (6*b1p2*b2*b3*u + 6*b1*b2p2*b3*u + 6*b1*b2*b3p2*u - 6*b1*b2*b3*u)*D


for i in range(0,nodes+1):
        V[1+i*5] = (-b1p3*u + 0.5*b1p2*u + 0.5*b1p4*u) + Lan1*(D2*b1p2 - 2*D*M*b1)
        V[2+i*5] = (-b2p3*u + 0.5*b2p2*u + 0.5*b2p4*u) + Lan1*(D2*b2p2 - 2*D*M*b2)
        V[3+i*5] = (-b3p3*u + 0.5*b3p2*u + 0.5*b3p4*u) + Lan1*(D2*b3p2 - 2*D*M*b3)

offset = Lan1*M*M

vartype = dimod.BINARY

model = dimod.BinaryQuadraticModel(V, Q, offset, vartype)

sampler = EmbeddingComposite(DWaveSampler())

sampleset = sampler.sample(model, num_reads = 1)

result = sampleset.first.sample

print(result)

#first use the exact solver which enumerates all possible states

# exactsolver = dimod.ExactSolver()

#results = exactsolver.sample(model)

#length = len(results)

#l = 0

#unsorted output 

#while l<length:
#   print(results.record[l][0], results.record[l][1], results.record[l][2])
#   l = l+1

#print(results)




