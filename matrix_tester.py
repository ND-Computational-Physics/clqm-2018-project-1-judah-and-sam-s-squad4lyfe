import numpy as np
import schrodinger_matrix_solver as sms
import schrodinger_to_matrix as stm

my_matrix = np.array(([1,0,3],[0,3,6],[0,2,2]))
print(sms.eigensolver(my_matrix)[0])
print('')
print(sms.eigensolver(my_matrix)[1])
