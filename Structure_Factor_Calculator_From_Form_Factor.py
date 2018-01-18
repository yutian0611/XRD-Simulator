#!/usr/bin/env python

import math
import numpy
import cmath

#=========================================== Simulation Constants Set Ups ============================================
Lattice_Constant=4.19
# Imaginary_Number is a complex number that only has the imaginary part, 1i, can be used as i in calculation.
Imaginary_Number=complex(0,1)
# A unit vector in the direction of the normal vector of the 001 planes
Normal_Unit=numpy.array([0,0,1])
# Wavelength of X-ray
Lambda_XR=1.540598
# Starting point of theta
Theta_In_Radians=0
#================== A class that contains element's all the information used in the calculation of Form Factor ==================
class Element_Class:
    def __init__(self,E_Name,a_1,a_2,a_3,a_4,b_1,b_2,b_3,b_4,E_c):
        self.Element_Name=E_Name      
        self.a=[a_1,a_2,a_3,a_4]
        self.b=[b_1,b_2,b_3,b_4]
        self.c=E_c
      
        '''
    def Form_Factor_Calculator(self):
        temp_F=self.c
        for i in range(4):
            temp_F=temp_F+self.a[i]*math.exp(-self.b[i]*(math.sin(Theta_In_Radians)/Lambda_XR)*(math.sin(Theta_In_Radians)/Lambda_XR))
        return temp_F
        '''
# Here is 4 element objects used in our experiment

Ba_Element_2_Positive=Element_Class("Ba",20.1807,19.1136,10.9054,0.77634,3.21367,0.28331,20.0558,51.746,3.02902)
Zr_Element_4_Positive=Element_Class("Zr",18.1668,10.0562,1.01118,-2.6479,1.2148,10.1483,21.6054,-0.10276,9.41454)
O_Element_2_Negative=Element_Class("O",3.0485,2.2868,1.5463,0.867,13.2771,5.7011,0.3239,32.9089,0.2508)
Mg_Element_2_Positive=Element_Class("Mg",3.4988,3.8378,1.3284,0.8497,2.1676,4.7542,0.185,10.1411,0.4853)


#======================= A class that contains all the useful information of an atom in XRD =======================



class Atom_Class:
    def __init__(self, Atom_F,Atom_x,Atom_y,Atom_z):
        self.F=Atom_F
        self.x=Atom_x
        self.y=Atom_y
        self.z=Atom_z
        self.E=self.Element_Selection()
        self.R=numpy.array([self.x,self.y,self.z])
        print (self.F)
    # A function to calculate the form factor
    def Form_Factor_Calculator(self):
        # self.E is the element of THIS atom
        temp_E=self.E
        temp_F=temp_E.c
        for i in range(4):
            temp_F=temp_F+temp_E.a[i]*math.exp(-temp_E.b[i]*(math.sin(Theta_In_Radians)/Lambda_XR)*(math.sin(Theta_In_Radians)/Lambda_XR))
        if(int(self.F)==8):
            print (temp_F)
        return temp_F   
    
    def Structure_Factor_Calculator(self):
        #Calculate the dot product by the normal unit vector and the location vector of this atom (n*r)
        dot_product=Normal_Unit.dot(self.R)
        #Calculate the exponent used in the structure calculator with the dot product
        #Here the exponent of e is an imaginary number, i.e. an imaginary part of a complex number.
        e_power=Imaginary_Number*(dot_product*4*math.pi*math.sin(Theta_In_Radians)/Lambda_XR)
        #Calculate the structure factor of this atom with the exponent and the form factor calculated with function above.
        return (self.Form_Factor_Calculator())*(cmath.exp(e_power))
       # return (self.Atom_Element.Form_Factor_Calculator())*(cmath.exp(e_power))
       
    #==============================================================================================================================       
    # This is a function to let the program know which element does THIS atom object belong to, with the Z number in the data.
    def Element_Selection(self):
        if int(self.F)==56:
            return Ba_Element_2_Positive
        elif int(self.F)==40:
            return Zr_Element_4_Positive
        elif int(self.F)==8:
            return O_Element_2_Negative
        elif int(self.F)==12:
            return Mg_Element_2_Positive       
    
#========================= An function calculates relative intensity from structure factor =========================
def Intensity_Calculator(Comp_Temp):
    Intensity_Temp=Comp_Temp*(numpy.conjugate(Comp_Temp))
    return Intensity_Temp


#======================================================= Main Program ========================================================
#========================================== Create atoms with data from config file ==========================================
Atoms_Config=open('30_Lattices_MgO_BaZrO3.txt','r')
Lattice_Cell=[]
Number_Of_Atoms=int(Atoms_Config.readline())

#================ A loop to create an array, to store all the atoms objects with the parameters from config file ================
for j in range(0,Number_Of_Atoms):
    Temp_Atom_F=float(Atoms_Config.readline())
    Temp_Atom_x=float(Atoms_Config.readline())
    Temp_Atom_y=float(Atoms_Config.readline())
    Temp_Atom_z=float(Atoms_Config.readline())
    
    Temp_Atom=Atom_Class(Temp_Atom_F,Temp_Atom_x,Temp_Atom_y,Temp_Atom_z)
    
    Lattice_Cell.append(Temp_Atom)

# A file manipulater to write intensities in a file
Intensity_With_Theta=open("30_Lattices_MgO_BaZrO3_Result_With_F.txt","w")
#===================== A loop to go through from theta=0 to theta=180, in degrees, with steplength=0.1 degree =====================
for Theta in numpy.arange(0,180,0.1):
    Theta_In_Radians=math.radians(Theta)
    Total_Structure_Factor_At_This_Theta=0
    #===================== A loop to calculate the total structure factor of the lattice cell =====================
    for n in range(0,Number_Of_Atoms):
        Total_Structure_Factor_At_This_Theta=Total_Structure_Factor_At_This_Theta+Lattice_Cell[n].Structure_Factor_Calculator()
    print ("The Total Structure Factor At ",Theta,"Degree is ",Total_Structure_Factor_At_This_Theta)
    Intensity_At_This_Theta=Intensity_Calculator(Total_Structure_Factor_At_This_Theta)
    print ("The Relative Intensity At ",Theta,"Degree is ",Intensity_At_This_Theta)
    print (Theta_In_Radians,Intensity_At_This_Theta.real,"\n")
        
    Intensity_With_Theta.write(str(Theta))
    Intensity_With_Theta.write("  ")
    Intensity_With_Theta.write(str(Intensity_At_This_Theta.real))
    Intensity_With_Theta.write("\n")

Intensity_With_Theta.close()
print ("=======================================================================")

input()
