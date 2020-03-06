import numpy as np
import matplotlib.pyplot as plt

BOUND_ANGLE = 250

def convex_angle(point_1, point_2, point_3):
    crossp = ((point_2[0] - point_1[0]) * (point_3[1] - point_1[1]) -
              (point_2[1] - point_1[1]) * (point_3[0] - point_1[0]))
    return crossp > BOUND_ANGLE

def inside_triangle(point_1, point_2, point_3, point_to_analise):
    column_with_ones = np.array([1, 1, 1])

    combinations = [[point_1, point_2], [point_2, point_3], [point_3, point_1]]
    determinants = []

    for combination in combinations:
        point_a, point_b = combination

        matrix = np.matrix([column_with_ones,
                               np.array([point_a[0], point_b[0], point_to_analise[0]]),
                               np.array([point_a[1], point_b[1], point_to_analise[1]])])

        determinant = np.linalg.det(matrix.T)
        determinants.append(determinant)
    
    min_determinant = min(determinants)
    
    is_inside_triangle = False
    if min_determinant >= 0:
        is_inside_triangle = True
    
    return is_inside_triangle

def is_ear(point_1, point_2, point_3, index_to_disconsider, polygon):
    points = [point_1, point_2, point_3]

    is_ear = True

    if convex_angle(point_1, point_2, point_3):
        for i, point_to_analise in enumerate(polygon):
            # To desconsider the point that is the ear
            if i == index_to_disconsider:
                continue
            repeated_point = False
            # To don't try repeated points of my triangle
            for point in points:
                if np.array_equal(point, point_to_analise):
                    repeated_point = True
                    break
            if repeated_point:
                continue
            
            if inside_triangle(point_1, point_2, point_3, point_to_analise):
                is_ear = False
                break

        return is_ear

def get_ears(polygon):
    index_ears_list = []
    
    for i in range(len(polygon)):
        point_1 = polygon[i - 1]
        point_2 = polygon[i]
        point_3 = polygon[(i + 1) % len(polygon)]

        points = [point_1, point_2, point_3]
        
        is_ear = True

        if convex_angle(point_1, point_2, point_3):
            for point_to_analise in polygon:
                repeated_point = False
                # To don't try repeated points of my triangle
                for point in points:
                    if np.array_equal(point, point_to_analise):
                        repeated_point = True
                        break
                if repeated_point:
                    continue
                
                if inside_triangle(point_1, point_2, point_3, point_to_analise):
                    is_ear = False
                    break

            if is_ear:        
                index_ears_list.append(i)

    return index_ears_list

def ears_clipping(polygon):
    copy_polygon = np.copy(polygon)
    
    # Get all the first ears of solution
    index_ears_list = get_ears(polygon)
    
    edges_solution = []
    size_polygon = len(polygon)
    
    while index_ears_list and size_polygon > 3:
        # Ear that will be used
        index_ear = index_ears_list.pop(0)

        # Neibourhood of ear to test if the vertexes used to construct the diagonal become ear 
        index_ear_plus_one = (index_ear + 1) % len(polygon)
        index_ear_plus_two = (index_ear + 2) % len(polygon)
        index_ear_minus_one = (index_ear - 1) % len(polygon)
        index_ear_minus_two = (index_ear - 2) % len(polygon)
        
        # Add vertexes that construct a diagonal in the solution
        edges_solution.append((polygon[index_ear_minus_one], polygon[index_ear_plus_one]))

        # Delete used vertexes to construct diagonal from the list of ears if they are there
        if index_ear_minus_one in index_ears_list:
            index_ears_list.remove(index_ear_minus_one)
        if index_ear_plus_one in index_ears_list:
            index_ears_list.remove(index_ear_plus_one)

        # Test for the vertexes used to costruct diagonal if they are ears, if yes add only one, because we don't have to ears consecutive
        last_point_is_ear = False
        if is_ear(polygon[index_ear_minus_one], polygon[index_ear_plus_one], polygon[index_ear_plus_two], index_ear, polygon):
            last_point_is_ear = True
            index_ears_list.append(index_ear_plus_one)

        if not last_point_is_ear and is_ear(polygon[index_ear_minus_two], polygon[index_ear_minus_one], polygon[index_ear_plus_one], index_ear, polygon):
            index_ears_list.append(index_ear_minus_one)

        # Delete used ear
        polygon = np.delete(polygon, index_ear, 0)
        size_polygon -= 1

        # Decrement indexes bigger than of the element used as ear
        i = 0
        while i < len(index_ears_list):
            if index_ears_list[i] > index_ear:
                index_ears_list[i] -= 1
            
            i += 1

    # Construct a list of the indexes of points used in diagonals
    index_edges = []
    edges_solution = np.array(edges_solution)
    for e in edges_solution:
        u = e[0, :]
        v = e[1, :]
        index_edges.append((np.where(np.all(u == copy_polygon, axis=1))[0][0],
                            np.where(np.all(v == copy_polygon, axis=1))[0][0]))

    return index_edges

def print_solution(file_name, diagonals):       
    with open(file_name, mode='w+') as fp:
        for diagonal in diagonals:
            print(str(diagonal[0]), str(diagonal[1]), file=fp)
        fp.close()

def draw_polygon(polygon, new_edges, datasetname):
    polygon = np.append(polygon, [polygon[0]], axis=0)
    x_polygon = polygon[:,0]
    y_polygon = polygon[:,1]

    plt.plot(x_polygon, y_polygon, linewidth=1, color="black")

    for e in new_edges:
        u = polygon[e[0]]
        v = polygon[e[1]]
        uv = np.array([u, v])
        plt.plot(uv[:, 0], uv[:, 1], linewidth=1, color="green")
    
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.suptitle(str(datasetname))
    plt.savefig("fig/"+ str(datasetname) +".png")
    plt.close()

def run():
    datasets = ["dataset/polygon1.txt", "dataset/polygon2.txt"]

    for dataset in datasets:
        # Dataset path
        INPUT_PATH = dataset
        # Load dataset
        polygon = np.loadtxt(INPUT_PATH).astype(np.float)
        # Call function to construct solution
        index_edges = ears_clipping(polygon)
        # Name of solution file
        output_file_name = 'diagonals{}.txt'.format(INPUT_PATH[-5])
        # Function to save solution file
        print_solution("solution/{}".format(output_file_name), index_edges)
        # Dataset name
        datasetname = output_file_name.split(".")[0]
        # Draw solution
        draw_polygon(polygon, index_edges, datasetname)

run()