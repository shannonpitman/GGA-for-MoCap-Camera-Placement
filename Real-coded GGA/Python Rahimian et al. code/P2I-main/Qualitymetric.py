from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import math
from numpy import linalg as LA
origine = [0,0,0,1]
vecteurmondex = [1, 0, 0,0]
vecteurmondey = [0, 1, 0,0]
vecteurmondez = [0, 0, 1,0]
p0 = [1, 0, 0]
p1 = [0, 1, 0]
p2 = [0, 0, 1]

# Matrice de rotation 3D autour de l'axe x
def matrice_rotation_3D_X(teta):
     teta = math.radians(teta)
     return np.array([[1, 0, 0, 0], [0, math.cos(teta), -math.sin(teta), 0], [0, math.sin(teta), math.cos(teta), 0], [0, 0, 0, 1]])
# Matrice de rotation autour axe y
def matrice_rotation_3D_Y(teta):
    teta = math.radians(teta)
    return np.array([[math.cos(teta), 0, math.sin(teta), 0], [0, 1, 0, 0], [-math.sin(teta), 0, math.cos(teta), 0], [0, 0, 0, 1]])
# Matrice de rotation autour axe z
def matrice_rotation_3D_Z(teta):
     teta = math.radians(teta)
     return np.array([[math.cos(teta), -math.sin(teta), 0, 0], [math.sin(teta), math.cos(teta), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
# Matrice de Translation
def matrice_translation_3D(Tx, Ty, Tz):
    return np.array([[1, 0, 0, Tx], [0, 1, 0, Ty], [0, 0, 1, Tz], [0, 0, 0, 1]])
def point_on_line(p, v,origine):
    """
    The function takes in three parameters, a point, a vector, and an origin, and returns whether the
    point lies on the line defined by the vector and origin.
    
    :param p: A point on the line (as a tuple or list of coordinates)
    :param v: A list or tuple representing the vector in 2D or 3D space
    :param origine: A list or tuple representing the origin of the vector
    """
    # On cherche à résoudre
    # p = q + t*v
    # Coordonnées du point
    x, y, z = p
    # Coordonnées du vecteur
    u, v, w = v
    # Coordonnées de l'origine du vecteur par rapport à la base canonique
    uO,vO,wO = origine
    # Vérification des valeurs nulles
    if u != 0:
        t = (x - uO) / u
    elif v != 0:
        t = (y - vO) / v
    elif w != 0:
        t = (z - wO) / w
    else:
        # Vecteur nul
        return False

    # Si le point p appartient à la droite engendrée par le vecteur
    # Si la valeur de t vérifie les 3 équations 
    print("coefficient V*t",t)
    if t <= 1:
        res = [t*u+uO, t*v+vO, t*w+wO]
        return res == p
    else:
        return False
def angle_between_vectors(v1,v2):
    """
    The function takes in two vectors and returns the angle between them.
    
    :param v1: A list or tuple representing the first vector in 2D or 3D space
    :param v2: A list or tuple representing the first vector in 2D or 3D space
    """
# V1 et V2 coordoonées du vecteur dans leur bases respectives
# Calcul de l'angle entre ces deux vecteurs
    produitscalaire = np.dot(v1,v2)
    print(produitscalaire)
    normev1 = LA.norm(v1,2)
    normev2 = LA.norm(v2,2)
    cos = produitscalaire/(normev1*normev2)
    angle = math.degrees(math.acos(cos))
    if angle>=40 and angle <= 140:
        print("angle between vectors in degree",angle)
        return angle
    else:
        print("L'angle ne respecte pas les conditions sur l'angle de convergence")
        return False
    ##
def intersection_between_vectors(v1,v2,origine1,origine2):
    # Coordoonées du vecteur V1
    u1,v1,w1 = v1
    # Coordonnées de l'origine du vecteur 1
    u1_O,v1_O,w1_O = origine1
    # Coordonnées du vecteur V2
    u2,v2,w2 = v2
   # Coordonnées de l'origine du vecteur 2 
    u2_O,v2_O,w2_O = origine2
    # Chaque vecteur peut engendrer une droite représenté sous la forme d'un système paramétrique
    # Ex x = t*u1 + u1_O
    # Avec t un coeff entre 0 et 1
    # Ainsi on détermine si un point est commun avec les deux droites 
    if u1!=u2:
        t = abs((u1_O-u2_O)/(u1-u2))
    elif v1!=v2:
        t = abs((v1_O-v2_O)/(v1-v2))
    elif w1!=w2:
        t = abs((w1_O-w2_O)/(w1-w2))
    else:
        return False
    if t<=1:
        res1 = [t*u1+u1_O,t*v1+v1_O,t*w1+w1_O]
        res2 = [t*u2+u2_O,t*v2+v2_O,t*w2+w2_O]
        if res1 == res2:
            print("Intersection entre les vecteurs",res1)
            return res1
        else:
            print("Pas d'intersection")
            return False
# Fonction permettant de réaliser une transformation géométrique par rapport au référentiel canonique
# Retourne les coordonnées de l'origine de la nouvelle base ainsi que les coordonnées des vecteurs unitaires de la base
# rotation en degrés et translation en array 3D
# axis en str pour déterminer selon quels axes
def define_rotation_translation(angle,translation,axis):
    """
    This function is intended to define a rotation and translation operation in three-dimensional space
    given an angle, translation vector, and axis of rotation.
    
    :param angle: The angle of rotation in degrees
    :param translation: A list or tuple of three values representing the translation vector in x, y, and
    z directions respectively
    :param axis: The axis parameter specifies the axis of rotation or translation. It can be either "x",
    "y", or "z" for the x-axis, y-axis, or z-axis respectively
    Returns the coordinates of the origin of the frame of reference and the coordinates of the unit vectors
    """
    if axis == "X":
            monderotationx = np.dot(matrice_rotation_3D_X(angle),vecteurmondex)
            monderotationy = np.dot(matrice_rotation_3D_X(angle),vecteurmondey)
            monderotationz = np.dot(matrice_rotation_3D_X(angle),vecteurmondez)
    elif axis == "Y":
            monderotationx = np.dot(matrice_rotation_3D_Y(angle),vecteurmondex)
            monderotationy = np.dot(matrice_rotation_3D_Y(angle),vecteurmondey)
            monderotationz = np.dot(matrice_rotation_3D_Y(angle),vecteurmondez)
    elif axis == "Z":
            monderotationx = np.dot(matrice_rotation_3D_Z(angle),vecteurmondex)
            monderotationy = np.dot(matrice_rotation_3D_Z(angle),vecteurmondey)
            monderotationz = np.dot(matrice_rotation_3D_Z(angle),vecteurmondez)
    else:
        print("Veuillez saisir X, Y ou Z dans axis.")
        # Normalisation des vecteurs
    monderotationx = monderotationx/np.linalg.norm(monderotationx)
    monderotationy =monderotationy/np.linalg.norm(monderotationy)
    monderotationz = monderotationz/np.linalg.norm(monderotationz)
    origine = np.dot(matrice_translation_3D(translation[0],translation[1],translation[2]),[0,0,0,1])
    return np.array([monderotationx,monderotationy,monderotationz,origine])

# Création d'une classe représentant les caméras
class Camera: 
    def __init__(self,monderotationx,monderotationy,monderotationz,origine,name):
        self.origine = origine[0:3]
        self.vector = np.dot(365,monderotationx[0:3])+np.dot(365,monderotationy[0:3])+np.dot(365,monderotationz[0:3])
        self.name = name

