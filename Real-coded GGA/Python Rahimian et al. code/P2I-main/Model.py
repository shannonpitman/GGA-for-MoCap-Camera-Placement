from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import math
from numpy import linalg as LA
from Qualitymetric import define_rotation_translation, matrice_translation_3D,matrice_rotation_3D_X,matrice_rotation_3D_Y,matrice_rotation_3D_Z,Camera,angle_between_vectors,intersection_between_vectors

print("Debut du programme")
# canonical frame of reference
vecteurmondex = [1, 0, 0,0]
vecteurmondey = [0, 1, 0,0]
vecteurmondez = [0, 0, 1,0]
#####
# Create pairs of cameras wich have a rotation along all the different axis
#####

##
# Along Z axis
##
# pairs of cameras have only different x positions
translation_one = [5,5,5]
translation_two = [8,5,5]
# Camera 1
referential_one = define_rotation_translation(0,translation_one,"Z")
camera_one = Camera(referential_one[0],referential_one[1],referential_one[2],referential_one[3],"camera_one")
# Camera 2
referential_two = define_rotation_translation(90,translation_two,"Z")
camera_two = Camera(referential_two[0],referential_two[1],referential_two[2],referential_two[3],"camera_two")
##
# Along Y axis
##
translation_three = [2,2,5]
translation_four = [2,2,8]
# Camera 3
referential_three = define_rotation_translation(0,translation_three,"Y")
camera_three = Camera(referential_three[0],referential_three[1],referential_three[2],referential_three[3],"camera_three")
# Camera 4
referential_four = define_rotation_translation(90,translation_four,"Y")
camera_four = Camera(referential_four[0],referential_four[1],referential_four[2],referential_four[3],"camera_four")
##
# Along X axis
##
translation_five = [4,2,4]
translation_six = [4,5,4]
# Camera 5
referential_five = define_rotation_translation(0,translation_five,"X")
camera_five = Camera(referential_five[0],referential_five[1],referential_five[2],referential_five[3],"camera_five")
# Camera 6
referential_six= define_rotation_translation(90,translation_six,"X")
camera_six = Camera(referential_six[0],referential_six[1],referential_six[2],referential_six[3],"camera_six")

# Test of position between pairs of cameras
# C1 and C2
angle_between_vectors(camera_one.vector,camera_two.vector)
intersection_between_vectors(camera_one.vector,camera_two.vector,camera_one.origine,camera_two.origine)
# C3 and C4
angle_between_vectors(camera_three.vector,camera_four.vector)
intersection_between_vectors(camera_three.vector,camera_four.vector,camera_three.origine,camera_four.origine)
# C5 and C6
angle_between_vectors(camera_five.vector,camera_six.vector)
intersection_between_vectors(camera_five.vector,camera_six.vector,camera_five.origine,camera_six.origine)

##########
# Cost function
##########
# Point
P=[5.5,3.5,5.5]
def cost_function():

   
    Q=0
    cameras = [camera_one,camera_two,camera_three,camera_four,camera_five,camera_six]
    
    for cam in cameras:
        print("cam",cam)
        new_list = cameras[:] #Slicing pour créer une copie de la liste cameras
        new_list.remove(cam)
        print("new_list",new_list)
        for cam2 in new_list:
            print("cam2", cam2)
            if angle_between_vectors(cam.vector,cam2.vector)!=False and intersection_between_vectors(cam.vector,cam2.vector,cam.origine,cam2.origine)!=False:
                if P == intersection_between_vectors(cam.vector,cam2.vector,cam.origine,cam2.origine):
                    print(f"Les deux caméras  {cam.name} et {cam2.name} respectent la condition de triangulabilité")
                else:

                    Q+=angle_between_vectors(cam.vector,cam2.vector)
                    print("Q",Q)
            else:
                print("Les deux caméras ne respectent pas les conditions de quality metric")
          
        
        
    # Avoid sum twice for one section 
    print("Q",Q/2)
    return Q/2
cost_function()
##########
# Display
##########
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111,projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
# set limit
ax.set(xlim=(-10, 10), ylim=(-10, 10), zlim=(0, 20))
# For camera 1 
X1,Y1,Z1 = [camera_one.origine, camera_one.origine, camera_one.origine]
U1,V1,W1 = [referential_one[0][0:3],referential_one[1][0:3],referential_one[2][0:3]]
# plot the new 3 axis with ax.quiver
ax.quiver(X1[0], Y1[1], Z1[2], U1[0], U1[1], U1[2], color="red",
          normalize = True, length = 1, label = 'x')
ax.quiver(X1[0], Y1[1], Z1[2], V1[0], V1[1], V1[2], color="red",
          normalize = True, length = 1, label = 'y')
ax.quiver(X1[0], Y1[1], Z1[2], W1[0], W1[1], W1[2], color="red",
          normalize = True, length = 1, label = 'z')
ax.plot([X1[0],camera_one.vector[0]],[Y1[1],camera_one.vector[1]],[Z1[2],camera_one.vector[2]],"r-")
# For camera 2 
X2, Y2, Z2 = [camera_two.origine,camera_two.origine,camera_two.origine]

U2, V2, W2 = [referential_two[0][0:3],referential_two[1][0:3], referential_two[2][0:3]]

ax.quiver(X2[0], Y2[1], Z2[2], U2[0], U2[1], U2[2], color="red",
          normalize = True, length = 1, label = 'x')
ax.quiver(X2[0], Y2[1], Z2[2], V2[0], V2[1], V2[2], color="red",
          normalize = True, length = 1, label = 'y')
ax.quiver(X2[0], Y2[1], Z2[2], W2[0], W2[1], W2[2], color="red",
          normalize = True, length = 1, label = 'z')

ax.plot([X2[0],camera_two.vector[0]],[Y2[1],camera_two.vector[1]],[Z2[2],camera_two.vector[2]],"r-")
# For camera 3
X3, Y3, Z3 = [camera_three.origine,camera_three.origine,camera_three.origine]

U3, V3, W3 = [referential_three[0][0:3],referential_three[1][0:3], referential_three[2][0:3]]

ax.quiver(X3[0], Y3[1], Z3[2], U3[0], U3[1], U3[2], color="blue",
          normalize = True, length = 1, label = 'x')
ax.quiver(X3[0], Y3[1], Z3[2], V3[0], V3[1], V3[2], color="blue",
          normalize = True, length = 1, label = 'y')
ax.quiver(X3[0], Y3[1], Z3[2], W3[0], W3[1], W3[2], color="blue",
          normalize = True, length = 1, label = 'z')

ax.plot([X3[0],camera_three.vector[0]],[Y3[1],camera_three.vector[1]],[Z3[2],camera_three.vector[2]],"b-")

# For camera 4
X4, Y4, Z4 = [camera_four.origine,camera_four.origine,camera_four.origine]

U4, V4, W4 = [referential_four[0][0:3],referential_four[1][0:3], referential_four[2][0:3]]

ax.quiver(X4[0], Y4[1], Z4[2], U4[0], U4[1], U4[2], color="blue",
          normalize = True, length = 1, label = 'x')
ax.quiver(X4[0], Y4[1], Z4[2], V4[0], V4[1], V4[2], color="blue",
          normalize = True, length = 1, label = 'y')
ax.quiver(X4[0], Y4[1], Z4[2], W4[0], W4[1], W4[2], color="blue",
          normalize = True, length = 1, label = 'z')

ax.plot([X4[0],camera_four.vector[0]],[Y4[1],camera_four.vector[1]],[Z4[2],camera_four.vector[2]],"b-")
# For camera 5
X5, Y5, Z5 = [camera_five.origine,camera_five.origine,camera_five.origine]

U5, V5, W5 = [referential_five[0][0:3],referential_five[1][0:3], referential_five[2][0:3]]

ax.quiver(X5[0], Y5[1], Z5[2], U5[0], U5[1], U5[2], color="green",
          normalize = True, length = 1, label = 'x')
ax.quiver(X5[0], Y5[1], Z5[2], V5[0], V5[1], V5[2], color="green",
          normalize = True, length = 1, label = 'y')
ax.quiver(X5[0], Y5[1], Z5[2], W5[0], W5[1], W5[2], color="green",
          normalize = True, length = 1, label = 'z')

ax.plot([X5[0],camera_five.vector[0]],[Y5[1],camera_five.vector[1]],[Z5[2],camera_five.vector[2]],"g-")
# For camera 6
X6, Y6, Z6 = [camera_six.origine,camera_six.origine,camera_six.origine]

U6, V6, W6 = [referential_six[0][0:3],referential_six[1][0:3], referential_six[2][0:3]]

ax.quiver(X6[0], Y6[1], Z6[2], U6[0], U6[1], U6[2], color="green",
          normalize = True, length = 1, label = 'x')
ax.quiver(X6[0], Y6[1], Z6[2], V6[0], V6[1], V6[2], color="green",
          normalize = True, length = 1, label = 'y')
ax.quiver(X6[0], Y6[1], Z6[2], W6[0], W6[1], W6[2], color="green",
          normalize = True, length = 1, label = 'z')

ax.plot([X6[0],camera_six.vector[0]],[Y6[1],camera_six.vector[1]],[Z6[2],camera_six.vector[2]],"g-")
ax.scatter(P[0],P[1],P[2],'o')
plt.show()
