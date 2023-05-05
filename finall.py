# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import numpy as np
import sys

sys.setrecursionlimit(10**6)

#codes by hossein armanpoor

#############*****************************************************************************************#############
#*************The following codes are the algoritm for generating non_overlap spheres in a finite volume**********#
#############*****************************************************************************************#############

#np.random.seed(1)
#user should change these values in order to achive desire foam
alpha, beta, gama, zeta, resize_factor, gam_zaviye = 1, 2, 1, 3, 1, np.pi/12
length_x_box, length_y_box, length_z_box = 100, 100, 100
tedad_halaty_ke_baraye_har_koreMarja_test_shavad = 144
x_center_spheres, y_center_spheres, z_center_spheres = [], [], []
r_spheres, all_coordinates, list_ref_spheres, refrences= [], [], [],[]
volume_sphere = []
r_min, r_max, r_ref = 3, 7, 0

#this function has to benefists for us: 1) can investigate in the L list in order to check for overlap 2)check for if somthing exists in a list
def investigate_all_elements(list, exist=False):
    All_is_possitive = True
    Most_Negative_L = None
    Already_exist = False
    for item in list:
        try:
            if item == exist:
                Already_exist = True
            if item < 0:
                All_is_possitive = False
                if not Most_Negative_L:
                    Most_Negative_L = item
                else:
                    if item < Most_Negative_L:
                        Most_Negative_L = item
        except Exception:
            if item == exist:
                Already_exist = True
    return All_is_possitive, Most_Negative_L , Already_exist

#append each refrence sphere cordinate in a list
def choose_RefSphere_randomly(all_coordinates,list_ref_spheres):
    iterates = 1
    index_tasadofy = np.random.randint(0,len(all_coordinates))
    exists = investigate_all_elements(list_ref_spheres, all_coordinates[index_tasadofy])[2]
    if exists == False:
        Current_RefSphere = all_coordinates[index_tasadofy]
        list_ref_spheres.append(Current_RefSphere)
        return Current_RefSphere
    else:
        if iterates == 15:
            Current_RefSphere = all_coordinates[index_tasadofy]
            list_ref_spheres.append(Current_RefSphere)
            return Current_RefSphere
        else:
            return (choose_RefSphere_randomly(all_coordinates,list_ref_spheres))

#this function generate random numbers for x,y,z and radius for new spheres
def find_coordinate(length_x_box, length_y_box, length_z_box, r_max):
    new_x = np.random.uniform(0, length_x_box)
    new_y = np.random.uniform(0, length_y_box)
    new_z = np.random.uniform(0, length_z_box)
    new_r = np.random.uniform(0.5*r_max, r_max)
    return new_x, new_y, new_z, new_r

#this function generate random numbers for x,y,z and radius for refrence spheres
def find_coordinate2(length_x_box, length_y_box, length_z_box, x_ref, y_ref, r_min, r_max, alpha, beta, zeta):
    new_r = np.random.uniform(r_min, 0.5*r_max)
    dr_min = (new_r + r_ref) * alpha
    dr_max = (new_r + r_ref) * beta
    dr = (dr_max - dr_min) * np.random.uniform(1, zeta) + dr_min
    global my_bool
    my_bool = True
    while my_bool:
        teta = np.random.uniform(0, np.pi)
        phi = np.random.uniform(0, np.pi)
        new_x = abs(dr*np.sin(phi)*np.cos(teta)) + x_ref
        new_y = abs(dr * np.sin(phi) * np.sin(teta)) + y_ref
        new_z = abs(dr * np.cos(phi) + z_ref)
        if new_x<length_x_box and new_y < length_y_box and new_z < length_z_box:
            my_bool = False
            break
        else:
            continue
    return new_x, new_y, new_z, new_r

#gerenerate L as the formule mentioned in the doc
def Produce_list_L(new_x, new_y, new_z, new_r, all_coordinates, gama):
    for k in range(len(all_coordinates)):
        L = np.sqrt((new_x - all_coordinates[k][0]) ** 2 + (new_y - all_coordinates[k][1]) ** 2
                    + (new_z - all_coordinates[k][2]) ** 2) - (new_r + all_coordinates[k][3]) * gama
        list_L.append(L)
    return list_L

#checks overlap depends on the method mentioned in document
def overlap_cheker(new_x, new_y, new_z, new_r, all_coordinates, gama, r_min):
    global list_L
    list_L = []
    list_L = Produce_list_L(new_x, new_y, new_z, new_r, all_coordinates, gama)
    All_positive, Most_Negative_L, v = investigate_all_elements(list_L)
    if All_positive:
        x_center_spheres.append(new_x)
        y_center_spheres.append(new_y)
        z_center_spheres.append(new_z)
        r_spheres.append(new_r)
        volume_sphere.append(4 / 3 * np.pi * new_r ** 3)
        all_coordinates.append([new_x, new_y, new_z, new_r])

        list_L = []
    else:
        rabete10_gostare_mojaz = new_r - abs(Most_Negative_L) > r_min
        if rabete10_gostare_mojaz:
            new_r = new_r - abs(Most_Negative_L)
            list_L = []
            overlap_cheker(new_x, new_y, new_z, new_r, all_coordinates, gama, r_min)

    return  x_center_spheres, y_center_spheres, z_center_spheres, r_spheres, volume_sphere, all_coordinates

#here the major algoritm starts
tedad_kore_marja = 30
print("The times that the Algoritm will choose refrence spheres =", tedad_kore_marja)
for i in range(tedad_kore_marja):
    #generate first sphere coordinate
    if not all_coordinates:
        r = np.random.uniform(r_min, r_max)
        x,y,z = 0,0,0
        x_center_spheres.append(x)
        y_center_spheres.append(y)
        z_center_spheres.append(z)
        r_spheres.append(r)
        volume_sphere.append(4/3*np.pi*r**3)
        x_ref,y_ref,z_ref,r_ref = x,y,z,r
        Current_RefSphere = [x_ref, y_ref, z_ref, r_ref]
        list_ref_spheres.append(Current_RefSphere)
        all_coordinates.append([x_center_spheres[0] , y_center_spheres[0] , z_center_spheres[0] , r])

    else:
        #generate possible spheres around the refrence sphere
        for j in range (tedad_halaty_ke_baraye_har_koreMarja_test_shavad):
            if j == tedad_halaty_ke_baraye_har_koreMarja_test_shavad - 1:
                print i*100/tedad_kore_marja, "% of spheres generated"
                Current_RefSphere = choose_RefSphere_randomly(all_coordinates, list_ref_spheres)
                x_ref ,y_ref , z_ref , r_ref= Current_RefSphere[0], Current_RefSphere[1], Current_RefSphere[2], Current_RefSphere[3]
                break

            if i<tedad_kore_marja/2:
                new_x, new_y, new_z, new_r = find_coordinate(length_x_box, length_y_box, length_z_box, r_max)

            else:
                new_x, new_y, new_z, new_r = find_coordinate2(length_x_box, length_y_box, length_z_box, x_ref, y_ref,
                                                      r_min, r_max, alpha, beta, zeta)

            x_center_spheres, y_center_spheres, z_center_spheres, r_spheres, volume_sphere, all_coordinates = \
                overlap_cheker(new_x, new_y, new_z, new_r, all_coordinates, gama, r_min)
            list_L = []

del all_coordinates
print("Algoritm for generating spheres has finished seccussfully")

#############*******************************************************#############
#*************The following codes are for link the algoritm and abaqus**********#
#############*******************************************************#############

print("Lets go for visualize them in abaqus")

# this code blongs to assembely module
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)

# create cube
def create_cube (length_x_box, length_y_box, length_z_box):
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0),
        point2=(length_x_box, length_y_box))
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='cube', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['cube'].BaseSolidExtrude(depth=length_z_box, sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']
create_cube(length_x_box, length_y_box, length_z_box)

def create_sphere_1 (radius, i):
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0,
        -100.0), point2=(0.0, 100.0))
    mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0)
        , direction=CLOCKWISE, point1=(0.0, radius), point2=(0.0, -radius))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, radius), point2=(
        0.0, -radius))
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='sphere%d'%i, type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['sphere%d'%i].BaseSolidRevolve(angle=360.0,
        flipRevolveDirection=OFF, sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']

if len(r_spheres)%2 != 0:
    tedad_koreha = len(r_spheres) + 1
else:
    tedad_koreha = len(r_spheres)

b = 5
# create sphere and translate them to there coordinate
for i in range (len(r_spheres)):
    #create
    a = i*100/tedad_koreha
    if int(a) %b == 0:
        print(b, "% of spheres visualized in abaqus")
        b += 5
    create_sphere_1(r_spheres[i], i)
    #translate
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='sphere%d-1' % i,
                                                part=mdb.models['Model-1'].parts['sphere%d' % i])
    mdb.models['Model-1'].rootAssembly.translate(instanceList=('sphere%d-1' % i,),
                                                 vector=(x_center_spheres[i], y_center_spheres[i], z_center_spheres[i]))


def function_density(volume_sphere, length_x_box , length_y_box, length_z_box):
    volume_box = length_x_box*length_y_box*length_z_box
    total_volume_sphere = round(np.sum(volume_sphere))
    density = round((((volume_box-total_volume_sphere) / volume_box) *100),2)
    return volume_box, total_volume_sphere, density

volume_box, total_volume_sphere, density = function_density(volume_sphere, length_x_box , length_y_box, length_z_box)
print "volume_box", volume_box
print "total_volume_sphere", total_volume_sphere
print "density", density, "%"
