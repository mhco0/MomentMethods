import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import sys

CONST_VACUUM_PERMITTIBITY = 8.854187812813e-12 
#F ⋅ m−1

def estimated_full_charge(dist_consts, delta):
    qpa = sum(dist_consts)

    area = delta ** 2

    return qpa * area

def plot_surface(central_points, dist_consts):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    xs = []
    ys = []
    zs = []

    for i in range(len(central_points)):
        xs.append(central_points[i][0])
        ys.append(central_points[i][1])
        zs.append(dist_consts[i])

    ax.scatter(xs, ys, zs)
    plt.show()

def find_dist_consts(impedance_matrix, tension_arr):
    # Z * a = V

    dist_consts = np.linalg.solve(impedance_matrix, tension_arr)

    return dist_consts

def build_impedance_matrix(central_points, delta):
    impedance_matrix = np.zeros((len(central_points), len(central_points)))
    
    for i in range(len(central_points)):
        for j in range(len(central_points)):
            if i != j:
                dx2 = (central_points[i][0] - central_points[j][0]) ** 2
                dy2 = (central_points[i][1] - central_points[j][1]) ** 2

                impedance_matrix[i][j] =  (delta * delta) / (4.0 * np.pi * CONST_VACUUM_PERMITTIBITY * np.sqrt(dx2 + dy2))
            else:
                impedance_matrix[i][j] = (delta * np.log(1.0 + np.sqrt(2))) / (np.pi * CONST_VACUUM_PERMITTIBITY)

    return impedance_matrix

def discretize(plate_lenght, delta):
    central_points = []

    for j in np.arange(0, plate_lenght, delta):
        for i in np.arange(0, plate_lenght, delta):

            x = i + (delta / 2.0)
            y = j + (delta / 2.0)
            
            if (x > plate_lenght / 2.0 and y > plate_lenght / 2.0):
                continue
            central_points.append([x, y])


    return central_points   

def mm_help():
    print("HELP GUIDE : ")
    print("")
    print(".Usage :")
    print("") 
    print("\tpython3 MomentMethods.py L V0 N [-q [-nplot]]")
    print("")
    print(".Details :")
    print("")
    print("L := The L length used for the metal plate (in m)")
    print("V0 := The constant tension used in the metal plate (in Volts)")
    print("N := A even number, represents the division in the metal plate")
    print("-q := Optional argument, with this set the code will print the code will print the estimated full charge.")
    print("-nplot := Optional argument, if this is set with the '-q' arg, the program won't plot the graph.")
    print("")
    print("--------------METAL PLATE--------------")
    print("")
    print("L ^ --------        ")
    print("  | --------        ")
    print("  | --------        =====> V0")
    print("  | --------'-------")
    print("  | ----------------")
    print("  | ----------------")
    print("0 | ________________>")
    print("  0        L/2      L")
    print(". L = N * delta")
    print(". ' = (L/2, L/2)")
    print("---------------------------------------")


def main():
    args_len = len(sys.argv)

    if args_len < 4:
        mm_help()
    else:
        print_plot = True
        print_q_total = False

        if args_len > 5:
            if(sys.argv[4] == "-q"):
                print_q_total = True
            else:
                print("Option " + sys.argv[4] + " not found.")

            if args_len == 6:
                if(sys.argv[5] == "-nplot"):
                    print_plot = False
                else:
                    print("Option " + sys.argv[5] + " not found.")

        plate_lenght = float(sys.argv[1])
        tension = float(sys.argv[2])
        sub_divisions = int(sys.argv[3])

        delta = plate_lenght / float(sub_divisions)

        central_points = discretize(plate_lenght, delta)

        impedance_matrix = build_impedance_matrix(central_points, delta)

        tension_arr = np.full(len(central_points), tension)

        dist_consts = find_dist_consts(impedance_matrix, tension_arr)

        if print_plot:
            plot_surface(central_points, dist_consts)

        if print_q_total:
            print("Qtotal : ", estimated_full_charge(dist_consts, delta))

if __name__ == "__main__":
    main()