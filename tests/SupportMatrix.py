#%%

import sympy as sp
import numpy as np

from IPython.display import display


#creates the d variable
d = 2

M = 5

doesNotExist = 'DNEx'


#creates the list of values
supportArrayInitial = []
#creates the list of values
for i in range(d+M):
    tempArray = []
    for j in range(d+M):

        #obtains the k and l variables
        k = d - i
        l = d - j
        #gets the k,l array
        k_l_array = np.array([k,l])

        #gets the current valid range
        currentRange = [-np.min(k_l_array), d+1-np.max(k_l_array)]
        
        
        if currentRange[0] >= currentRange[1]:
            currentRange = doesNotExist
        
        tempArray.append(currentRange)

    
    supportArrayInitial.append(tempArray)




display(supportArrayInitial)

supportArray = []
#goes through and obtains the actual support array
for i in range(d+M):
    tempArray = []
    for j in range(d+M):

        #gets the temp
        tempCondition = supportArrayInitial[i][j]

        if isinstance(tempCondition, str):
            tempArray.append(doesNotExist)
        else:
            #gets the lower condition
            lowerCondition = tempCondition[0]
            #gets the upper condition
            upperCondition =  tempCondition[1]
            if lowerCondition < 0:
                newLowerCondition = 0
            else:
                newLowerCondition = lowerCondition

            if upperCondition > M:
                newUpperCondition = M
            else:
                newUpperCondition = upperCondition

            tempArray.append([newLowerCondition, newUpperCondition])
    supportArray.append(tempArray)

display(supportArray)

# %%
