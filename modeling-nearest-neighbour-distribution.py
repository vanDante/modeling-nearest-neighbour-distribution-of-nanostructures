import numpy as np
import seaborn as sns
from scipy.special import k0

#--------------------------------------------------------------------------
# INTRODUCTION
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
# INPUT PARAMETERS
#--------------------------------------------------------------------------

# set trageted NW density per cm^2
density=1.1*10**9      

# set substrate spacings in nm
length=10000   
width=10000   

# set list diffusion length L
L_list=[120]

# set list for m (m is related to the number of atoms in the critical nucleus) 
m_list=[4]

# set experimentally observed diameter of NWs
rNW=20  



#--------------------------------------------------------------------------
# PREPROCESSING AND INITIALIZATTION OF TARGET PARAMETERS
#--------------------------------------------------------------------------
 
# number of NWs coressponding to density and substrate spcing 
N=round(density*length*width*10**(-14))

# number of sets resulting from L_list and m_list
sets=len(L_list)*len(m_list)


# create list to store 
NWlist=np.zeros((N,4)) #x value | y value | nearest neighbor distance | next nearest neighbor distance

allData=np.zeros((N,sets))
allParameters=np.zeros((2,sets))
        
maxi=int(np.sqrt(length**2+width**2))
plist=np.zeros((maxi,2))

NNdist=10**5*np.ones((N))



#--------------------------------------------------------------------------
# MODELING
#--------------------------------------------------------------------------

#CREATE NWs
counter=0
for L in L_list:
    for m in m_list:

        #solve besselfunction
        for k in range(maxi):
          plist[k,0]=k
          #prob function according to settings
          c=(1-k0(k/L)/k0(rNW/L))
          plist[k,1]=c**m
                
        #create first NW
        NWlist[0,0]=int(np.random.uniform()*length)
        NWlist[0,1]=int(np.random.uniform()*width)

        #create further NWs
        i=1
        while i<N:
          x=int(np.random.uniform()*length)
          y=int(np.random.uniform()*width)
          p=1
          for j in range(i):
            d=int(np.sqrt((x-NWlist[j,0])**2 + (y-NWlist[j,1])**2))
            if d==0:
              d=2*rNW
            p=p*plist[d,1]
          c=np.random.uniform()
          if p>=c: 
            NWlist[i,0]=x
            NWlist[i,1]=y
            i=i+1

        # NN DISTRIBUTION
        for i in range(N):
            for j in range(N):
                NNdist[j]=np.sqrt((NWlist[i,0]-NWlist[j,0])**2 + (NWlist[i,1]-NWlist[j,1])**2)
          
            NNdist=np.sort(NNdist)
            NWlist[i,2]=int(NNdist[1])
            NWlist[i,3]=int(NNdist[2])

        allData[:,counter]=NWlist[:,2]
        allParameters[0,counter]=L
        allParameters[1,counter]=m

        counter=counter+1
        # show progress in console window
        print('set '+str(counter)+'/'+str(sets))
        


#--------------------------------------------------------------------------
# PLOTTING
#--------------------------------------------------------------------------        
for i in range(counter): 
    plot = sns.kdeplot(allData[:,i], gridsize=50)
plot.set_xlabel("Nearest neighbor distance")

        

#--------------------------------------------------------------------------
# SAVE DATA
#--------------------------------------------------------------------------
