import random
import sys

## Generates a random dense matrix
## Use: python randomMatrix.py M N

if (len(sys.argv) > 2):

  M = int(sys.argv[1])
  N = int(sys.argv[2])

  m = sys.argv[1]
  n = sys.argv[2]

  f = open('matrix','w');


  table= [ [ 0 for i in range(M) ] for j in range(N) ]
  table2= [ [ 0 for i in range(N) ] for j in range(M) ]


  for i in range(N):
    for j in range(M):

      table[i][j] = random.random()
      table2[j][i] = table[i][j]
      f.write(str(table[i][j]))
      f.write('\t')

    f.write('\n')

  f.close()


  if((len(sys.argv) > 3) and (sys.argv[3] == "-py")):
  
    for i in range(M):
     print table2[i]



else:

  print "please input two numbers"

