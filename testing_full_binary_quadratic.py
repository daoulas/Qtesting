# This exercise searches for a minimum of a binary quadratic model obtained by "expanding" the compact form 
# J*sum_i={1,N-1} (1-s(i))*(1-s(i+1)), where s(i), i=1...N are binary variables. The goal of the exercise 
# is to understand how to use  the definition via dimod.BinaryQuadraticModel. 


# Uses exact solver

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

#print(Q)

#print(V)

#print(offset)

vartype = dimod.BINARY

model = dimod.BinaryQuadraticModel(V, Q, offset, vartype)

#first use the exact solver which enumerates all possible states

exactsolver = dimod.ExactSolver()

results = exactsolver.sample(model)

length = len(results)

l = 0

#unsorted output 

while l<length:
   print(results.record[l][0], results.record[l][1], results.record[l][2])
   l = l+1

print(results)

# We can directly consider the "compact" form J*sum_i={1,N-1} t(i)*t(i+1), where we substituted t(i) = 1 - s(i).
# The results should be the same to those we obtained from the expanded form (keeping  in mind that s(i) = 1 -t(i)

for i in range(1,N):
       Q[i,i+1] = J


for i in range(1,N+1):
       V[i] = 0

offset = 0

model = dimod.BinaryQuadraticModel(V,Q,offset,vartype)

exactsolver = dimod.ExactSolver()

results = exactsolver.sample(model)

print ('=========================')

length = len(results)

l = 0

#unsorted output

while l<length:
    print(1 - results.record[l][0], results.record[l][1], results.record[l][2])
    l = l+1

#print(results)



