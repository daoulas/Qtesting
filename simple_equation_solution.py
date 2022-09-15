from dwave.system import DWaveSampler,EmbeddingComposite

import dimod

Q = {(1,2): 52, (1,3): 16, (1,4): 32, (2,3): 32, (2,4): 32, 
     (3,4): 20, (1,1):-45, (2,2):-64, (3,3): -31, (4,4): -52}

model = dimod.BinaryQuadraticModel.from_qubo(Q, offset = 0.0)

exactsolver = dimod.ExactSolver()

results = exactsolver.sample(model)

print(results)
