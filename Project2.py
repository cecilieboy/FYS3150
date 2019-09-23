#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
from numpy import linalg as LA


# In[3]:


n = 2
A = 2*np.eye(n,n) + -1 * np.eye(n,n,k=-1) + -1 * np.eye(n,n,k=1)

print(A)

tau = (A[0,0] - A[1,1]) / (2 * A[1,0])
tplus = - tau + np.sqrt(1 + tau**2)
tminus = -tau - np.sqrt(1 + tau**2)

c = 1 / (np.sqrt(1 + tau**2))
s = c * tminus

theta = np.arcsin(s)

print(theta)

B = np.zeros((2,2))

B[0,0] = A[0,0] * c**2 - 2 * A[0,1] * c * s + A[1,1] * s**2

B[1,1] = A[0,0] * c**2  + 2 * A[0,1] * c * s + A[0,0] * s**2

B[0,1] = B[1,0] = (A[0,0] - A[1,1]) * c * s + A[0,1] * (c**2 - s**2)

print(B)


# In[32]:



def code(n):  
   
   
###########################   funktion for rotationen ###########################
   
   
   def rotation(A,l,k,v): 
# A er matricen der skal diagonaliseres
# l, k er positionerne for rotationsledende og positionerne i A der bliver minimeret
# v er en mastrix bestående af basisvektorene
       

       B = np.copy(A)
       
       w = np.copy(v)

       tau = (A[l,l] - A[k,k]) / (2 * A[k,l])
#tau er defineret sådan så vi får en andengradsligning

       tplus  = -tau + np.sqrt(1 + tau**2)
       tminus = -tau - np.sqrt(1 + tau**2)
#Udfra betingelsen at b_kl = 0, får vi en andengradsligning med disse løsninger for t        

       if abs(tplus) < abs(tminus):
           t = tplus
       else:
           t = tminus
#Vi skal åbenbart altid vælge den mindre t?
           
       c = 1 / (np.sqrt(1 + t**2))

       s = t * c

       for j in range(len(v)):
           print(k,l)
           print(j)
           print()
           w[j,k] =   c * v[j,k] - s * v[j,l]
           w[j,l] =   s * v[j,k] + c * v[j,l]
           #print(w[k,j])
           #print(w[l,j])
           
       print(w)
#Vi fandt ud af ved testing at, hvis vi bruger rotatione matricen, sker der kun noget med k'te og l'te række i vær emkel søjle
           

       B[k,k] = A[k,k] * c**2 - 2 * A[k,l] * c * s + A[l,l] * s**2
       B[l,l] = A[l,l] * c**2 + 2 * A[k,l] * c * s + A[k,k] * s**2
       B[k,l] = B[l,k] = 0
# Her sker selve rotationen

       for i in range(len(A)):


           if i != k and i != l:
               #B[i,i] = A[i,i]
               B[i,k] = B[k,i] = A[i,k] * c - A[i,l] * s
               B[i,l] = B[l,i] = A[i,l] * c + A[i,k] * s
#Alle andre led i k,l -række/søjle ændrer også værdi


       return B,w











   def condition(A,max):
#condition som stopper alogrithmen (returner False) ligeså snart den ikke finder nogen værdi som ligger over max


       for i in range(len(A)):
           for j in range(len(A)):


               if i !=j:
                   if np.abs(A[i,j]) > max:
                       return True

       return False






#Definere A udfra vores andengradslignings problem, definere v som matrix med enhedsvektorer, definere max:    
   A =          2*np.eye(n,n) + -1 * np.eye(n,n,k=-1) + -1 * np.eye(n,n,k=1)
   Aunchanged = 2*np.eye(n,n) + -1 * np.eye(n,n,k=-1) + -1 * np.eye(n,n,k=1)
   v = np.eye(n,n)
   print(len(v))
   max = 10**(-10)   
   count = 0
   max_count = 100
   
   
   
   
   
   while condition(A,max) and count < max_count:    


       highestmatrixentry = 0  


       for i in range(len(A)):
           for j in range(len(A)):

               if i != j: 

                   if np.abs(A[i,j]) > highestmatrixentry:

                       highestmatrixentry = np.abs(A[i,j])
                       l, k = i, j
#Finder positionen af den højeste matrice indgang, med et loop som går igennem alle indgangene og opdaterer den højeste indgangs vœrdi/position

       
       A,v = rotation(A,l,k,v)
       #print("===============")
       #print(A)
       print(k,l)


       count += 1


       
       
       
       
       
       
       
       

   print('diagonalized matrix:')
   print(A)
   #print('number of rotations =', count)
   print('Matrix of eigvectors:')
   print(v)
   
   print()
   
   #print(np.dot(np.transpose(v),v))
   
   eigvalues = np.diag(A)
   #print(eigvalues)
   
   a1 = np.dot(Aunchanged,v[:,0]) 
   a2 = eigvalues[0] * v[:,0]
   
   print('own test')
   print(a1)
   print(a2)
   
   #print(a1,a2)
   
   print('pythons solutions:')
   pythoneigvalues, pythoneigvectors = LA.eig(Aunchanged)
   print(pythoneigvalues)
   print(pythoneigvectors)
   print(np.dot(Aunchanged,pythoneigvectors[:,0]))
   print(pythoneigvalues[0] * pythoneigvectors[:,0])
   
   
   
code(3)  

           
           
   
   



           
               

               
                   
                   
                   
                   
        

   


# In[ ]:


A = [[1,4],[1,3]]
s = np.max (A)
print(s)
p = len(A)
print(p)


# In[ ]:





# In[ ]:




