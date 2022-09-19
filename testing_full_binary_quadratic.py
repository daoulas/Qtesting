
from dwave.system import DWaveSampler,EmbeddingComposite
  
import dimod

Q = {}

V = {}

N = 5 

J = 10

for i in range(1,N):
        Q[i,i+1] = J

V[1] = -J
V[N] = -J

for i in range(2,N):
        V[i] = -2*J

offset = J*(N-1)

print(Q)

print(V)

print(offset)

vartype = dimod.BINARY

model = dimod.BinaryQuadraticModel(V, Q, offset, vartype)

exactsolver = dimod.ExactSolver()

results = exactsolver.sample(model)

print(results)
