

from dwave.system import DWaveSampler,EmbeddingComposite

import dimod

import dwave.inspector

Q = {(1,2): 52, (1,3): 16, (1,4): 32, (2,3): 32, (2,4): 64, 
     (3,4): 20, (1,1):-45, (2,2):-64, (3,3): -31, (4,4): -52}

model = dimod.BinaryQuadraticModel.from_qubo(Q, offset = 0.0)

#exactsolver = dimod.ExactSolver()

#results = exactsolver.sample(model)

#print(results)

#Checking various posibilities for out put

sampler = dimod.ExactSolver()

sampleset = sampler.sample(model)

print(sampleset)

print("=======")

print(sampleset.variables)

print("=======")

print(sampleset.record)

print("=======")

print(sampleset.first)

print("=======")

#This writes out the values of the four variables of the minimum solution
for i in range(1,5):
  print(sampleset.first[0][i])

#Now use the quantum sampler

sampler = EmbeddingComposite(DWaveSampler())

sampleset = sampler.sample(model, num_reads = 15)

print("Now the results from quantum sampler")

print(sampleset)

print("========")

print(sampleset.first)

for i in range (1,5):
 print(sampleset.first[0][i])

dwave.inspector.show(sampleset)

