#creates the scripts to display the Sections and integrations for an example basis function

import numpy as np
import matplotlib.pyplot as plt



#creates the function definition 
def degree2Basis(t: np.ndarray): #time variable
    

    outputList = []

    for i in range(len(t)):

        currentTime = t.item(i)
        evaluation = 0.0

        if 0.0 <= currentTime and currentTime < 1.0:
            evaluation = (1.0/2.0)*(currentTime**2)
        elif 1.0 <= currentTime and currentTime < 2.0:
            evaluation = -(currentTime**2) + 3*currentTime - (3.0/2.0)
        elif 2.0 <= currentTime and currentTime < 3.0:
            evaluation = ((3.0-currentTime)**2)/(2.0)
        else:
            evaluation = 0.0

        #appends the evaluation to the list
        outputList.append(evaluation)


    outputList = np.array(outputList)


    return outputList

        


startTime = 0.0
endTime = 1.0
numSamples = 1000

#creates a linspace 
t = np.linspace(startTime, endTime, numSamples)


#gets the evaluation for first, second, and third sections

firstSection = degree2Basis(t=t)

secondSection = degree2Basis(t=(t+1.0))

thirdSection = degree2Basis(t=(t+2.0))


plt.figure(figsize=(10,3))
plt.plot(t, firstSection, color='blue')
plt.plot(t, secondSection, color='blue')
plt.plot(t, thirdSection, color='blue')
plt.xticks(np.arange(0.0, 1.1, 0.5))
plt.yticks(np.arange(0.0, 1.1, 0.5))
plt.title("Degree 2 Polynomial Basis sections")
plt.xlabel('time')
plt.ylabel('f(t)')
plt.show()





potatoSalad = 0
