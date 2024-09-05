import numpy as np, trimesh
from scipy.spatial.transform import Rotation as R
from shapely.geometry import Polygon
from shapely.affinity import rotate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Hull(trimesh.Trimesh):
    def __init__(self, *args, **kwargs):
        '''
        :param args:
        :param kwargs:
        Class Renaming: The class CustomTrimesh is renamed to Hull.
        Union (+):
            Hull1 + Hull2: Combines two Hull objects (Union).
            Hull1 + [x, y, z]: translate Hull1 the displacement of [x, y, z]
        Difference (-):
            Hull1 - Hull2: Computes the difference between two Hull objects.
        Intersection (*):
            Hull1 * Hull2: Computes the intersection of two Hull objects.
        Rotation (@):
            Hull1 @ [angle_x, angle_y, angle_z]: Rotates a Hull object by specified angles.
                angle_x, angle_y, angle_z are rotation angles around x, y, z in radians
        Twist (^):
            Hull1 ^ [angle_x, angle_y, angle_z]: Twists a Hull object by specified angles per unit length along the x, y, and z axes.
            angle_x, angle_y, angle_z are rotation angles per unit length around x, y, z in radians
        Shear (%):
            Hull1 % [shear_zx, shear_yz, shear_zy, shear_xy, shear_yz, shear_xz]: Shears a Hull object by specified shear factors.
            shear_zx, shear_yz, shear_zy, shear_xy, shear_yz, shear_xz are the shear factors (The distance a point moves due to shear divided by the perpendicular distance of a point from the invariant line.) on the specific planes respectively.
                e.g. shear_xy: shear around x-axis to the direction of y
        '''


        if len(args) == 1 and isinstance(args[0], trimesh.Trimesh):
            original_mesh = args[0]
            super().__init__(
                vertices=original_mesh.vertices,
                faces=original_mesh.faces,
                vertex_normals=original_mesh.vertex_normals if original_mesh.vertex_normals is not None else None,
                face_normals=original_mesh.face_normals if original_mesh.face_normals is not None else None,
                metadata=original_mesh.metadata if original_mesh.metadata is not None else None
            )
        else:
            vertices = kwargs.get('vertices', None)
            faces = kwargs.get('faces', None)
            if vertices is not None and faces is None:
                convex_hull = trimesh.convex.convex_hull(vertices)
                kwargs['vertices'] = convex_hull.vertices
                kwargs['faces'] = convex_hull.faces
            super().__init__(*args, **kwargs)

    def __add__(self, other):
        if isinstance(other, trimesh.Trimesh):
            if not self.is_watertight or not other.is_watertight:
                raise ValueError("Both meshes must be watertight (closed volumes) to perform union.")
            return Hull(**trimesh.boolean.union([self, other]).to_dict())
        elif isinstance(other, (np.ndarray, list, tuple)) and len(other) == 3:
            translation_matrix = trimesh.transformations.translation_matrix(other)
            translated_mesh = self.copy()
            translated_mesh.apply_transform(translation_matrix)
            return Hull(**translated_mesh.to_dict())
        else:
            raise ValueError("Unsupported operand type for +: 'Hull' and '{}'".format(type(other)))

    def __sub__(self, other):
        if isinstance(other, trimesh.Trimesh):
            if not self.is_watertight or not other.is_watertight:
                raise ValueError("Both meshes must be watertight (closed volumes) to perform difference.")
            return Hull(**trimesh.boolean.difference([self, other]).to_dict())
        elif isinstance(other, (np.ndarray, list, tuple)) and len(other) == 3:
            translation_matrix = trimesh.transformations.translation_matrix(-np.array(other))
            translated_mesh = self.copy()
            translated_mesh.apply_transform(translation_matrix)
            return Hull(**translated_mesh.to_dict())
        else:
            raise ValueError("Unsupported operand type for -: 'Hull' and '{}'".format(type(other)))

    def __mul__(self, other):
        if isinstance(other, trimesh.Trimesh):
            if not self.is_watertight or not other.is_watertight:
                raise ValueError("Both meshes must be watertight (closed volumes) to perform intersection.")
            return Hull(**trimesh.boolean.intersection([self, other]).to_dict())
        else:
            raise ValueError("Unsupported operand type for *: 'Hull' and '{}'".format(type(other)))

    def __matmul__(self, other):
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 3:
            angle_x, angle_y, angle_z = other
            rotation_matrix = trimesh.transformations.euler_matrix(angle_x, angle_y, angle_z)
            rotated_mesh = self.copy()
            rotated_mesh.apply_transform(rotation_matrix)
            return Hull(rotated_mesh)
        else:
            raise ValueError("Unsupported operand type for @: 'Hull' and '{}'".format(type(other)))

    def __xor__(self, other):

        '''
        Twist an object
        :param angle_x: Twist angle per unit length along the x-axis
        :param angle_y: Twist angle per unit length along the y-axis
        :param angle_z: Twist angle per unit length along the z-axis
        :return: twisted object in self.hull
        '''
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 3:
            angle_x, angle_y, angle_z = other
        else: angle_x, angle_y, angle_z = (0, 0, 0)

        # Get the vertices of the mesh
        vertices = self.vertices.copy()

        # Apply the twist transformation
        for i, vertex in enumerate(vertices):
            x, y, z = vertex

            # Calculate the twist angles based on the position along the axes
            twist_angle_x = angle_x * x
            twist_angle_y = angle_y * y
            twist_angle_z = angle_z * z

            # Rotation matrices for twisting around the x, y, and z axes
            rotation_matrix_x = np.array([
                [1, 0, 0],
                [0, np.cos(twist_angle_x), -np.sin(twist_angle_x)],
                [0, np.sin(twist_angle_x), np.cos(twist_angle_x)]
            ])

            rotation_matrix_y = np.array([
                [np.cos(twist_angle_y), 0, np.sin(twist_angle_y)],
                [0, 1, 0],
                [-np.sin(twist_angle_y), 0, np.cos(twist_angle_y)]
            ])

            rotation_matrix_z = np.array([
                [np.cos(twist_angle_z), -np.sin(twist_angle_z), 0],
                [np.sin(twist_angle_z), np.cos(twist_angle_z), 0],
                [0, 0, 1]
            ])

            # Combined rotation matrix
            combined_rotation_matrix = rotation_matrix_x @ rotation_matrix_y @ rotation_matrix_z

            # Apply the combined rotation to the vertex
            new_vertex = combined_rotation_matrix @ np.array([x, y, z])
            vertices[i] = new_vertex

        return Hull(vertices = vertices)

    def __mod__(self, other):
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 6:
            shear_matrix = np.eye(4)
            shear_matrix[0, 1] = other[0]
            shear_matrix[0, 2] = other[1]
            shear_matrix[1, 0] = other[2]
            shear_matrix[1, 2] = other[3]
            shear_matrix[2, 0] = other[4]
            shear_matrix[2, 1] = other[5]
            sheared_mesh = self.copy()
            sheared_mesh.apply_transform(shear_matrix)
            return Hull(sheared_mesh)
        else:
            raise ValueError("Unsupported operand type for %: 'Hull' and '{}'".format(type(other)))

    def __str__(self):
        l = self.vertices.tolist()
        s = 'Hull(vertices = [' + str(l[0])
        for row in l[1:]:
            s += ',\n' + str(row)

        s += '],'
        l = self.faces.tolist()
        s += '\n   faces = [' + str(l[0])
        for row in l[1:]:
            s += ',\n' + str(row)

        return s + '])'

    def __repr__(self):
        return self.__str__()


def showx(Hull):
    '''
    Display a Hull object.
    :param Hull: Hull object to display
    :return: a window to display the hull
    '''
    scene = trimesh.Scene(Hull)
    scene.add_geometry(__create_custom_axis(length=1.0, divisions=10, radius=0.01))
    scene.show()

def __create_custom_axis(length=1.0, divisions=10, radius=0.01):
    axis = trimesh.creation.axis(origin_size=radius, origin_color=[1, 1, 1, 1], axis_length=length, axis_radius=radius)
    for i in range(divisions + 1):
        tick_pos = length * i / divisions
        x_tick = trimesh.creation.cylinder(radius=radius / 2, height=radius * 2, sections=8)
        x_tick.apply_translation([tick_pos, 0, 0])
        y_tick = trimesh.creation.cylinder(radius=radius / 2, height=radius * 2, sections=8)
        y_tick.apply_translation([0, tick_pos, 0])
        z_tick = trimesh.creation.cylinder(radius=radius / 2, height=radius * 2, sections=8)
        z_tick.apply_translation([0, 0, tick_pos])
        axis = axis + x_tick + y_tick + z_tick
    return axis

def cuboid(a=1, b=1, c=1):
    '''
    Create a cuboid of sides a, b, c
    :param a: side a
    :param b: side b
    :param c: side c
    :return: a cuboid in self.hull
    '''
    obj = trimesh.creation.box(extents=[a, b, c])
    return Hull(obj)


def cube(a=1):
    '''
    Create a cube of side a
    :param a: side length of the cube
    :return: a cube in self.hull
    '''
    return cuboid(a=a, b=a, c=a)

def __get_num_edges(radius=1, division = 0.1):
    '''
    Get the number of edges of a circle from self.division.
    The minimum number of sides is 6.
    :param radius: radius of the circle
    :return: the number of sides
    '''
    return max(int(2 * np.pi * radius / division), 6)

def cylinder(radius=1, height=1):
    '''
    Create a cylinder of radius and height along a
    :param radius: radius of the bottom circle
    :param height: height along z
    :return: a cylinder along z
    '''
    obj = trimesh.creation.cylinder(radius=radius, height=height, sections=__get_num_edges(radius=radius))
    return Hull(vertices=obj.vertices, faces=obj.faces)

def cone(radius=1, height=1):
    '''
    Create a cone in radius a and height along z
    :param radius: radius of the bottom of cone
    :param height: height of the cone
    :return: a cone in self.hull
    '''
    obj = trimesh.creation.cone(radius=radius, height=height, sections=__get_num_edges(radius=radius))
    return Hull(vertices=obj.vertices, faces=obj.faces)

def ellipsoid(a=1, b=1, c=1, division = 0.1):
    '''
    Create an ellipsoid (an oval object)
        (x/a)**2 + (y/b)**2 + (z/c)**2 = 1
    :param a: radius along x
    :param b: radius along y
    :param c: radius along z
    :param division: the distance between neighboring x, y points
    :return: an ellipsoid in self.hull
    '''
    vertices = np.empty((0, 3))
    for x in np.arange(-a, a + division, division):
        y_limit = b * (1 - (x / a) ** 2) ** 0.5
        if not np.isnan(y_limit) and not np.isinf(y_limit):
            for y in np.arange(-y_limit, y_limit + division, division):
                z = c * (1 - (x / a) ** 2 - (y / b) ** 2) ** 0.5
                if not np.isnan(z):
                    p = np.array([[x, y, z], [x, y, -z]])
                    vertices = np.append(vertices, p, axis=0)

    return Hull(vertices = vertices)

def sphere(radius=1, face_index = 3):
    '''
    Create a sphere of radius
    :param radius: spherical radius
           face_index: 3 (4**3 = 64 faces) or 4 (4**4 = 256 faces)
    :return: a sphere in self.hull
    '''
    #num = __get_num_edges(radius=radius, division = division))
    obj = trimesh.creation.icosphere(subdivisions= face_index, radius = radius)
    return Hull(vertices=obj.vertices, faces=obj.faces)

def __polygon_vertices(radius, num_sides):
    '''
    Internal function to generate the vertices of a regular polygon from raidus and num_sides
    :param radius: radius of the regular polygon
    :param num_sides: number of sides of the regular polygon
    :return: ndarray of vertices of the regular polygon
    '''
    vertices = []
    angle_increment = 2 * np.pi / num_sides

    for i in np.arange(num_sides):
        x = radius * np.cos(i * angle_increment)
        y = radius * np.sin(i * angle_increment)
        vertices.append([x, y])

    return np.array(vertices)

def pyramid(radius=1, n=3, height=1, division = 0.1):
    '''
    A pryramid of the bottom of n sides of regular polygon in radius and height along z.
    :param radius: radius of the regular polygon
    :param n: number of sides of polygon
    :param height: height along z
    :return: a regular pyramid in self.hull
    '''

    vertices = np.empty((0, 3))
    dr = radius * division / height  # reduction of radius per division height
    for m, z in enumerate(np.arange(0, height + division, division)):
        r = radius - m * dr
        if r >= 0:
            polygon_vertices = __polygon_vertices(radius=r, num_sides=n)

            # Create a column vector with the number to append
            column_to_append = np.full((polygon_vertices.shape[0], 1), r)

            # Append the column to the original array
            vertices_3d = np.hstack((polygon_vertices, column_to_append))

            vertices = np.append(vertices, vertices_3d, axis=0)
        else:
            break

    return Hull(vertices = vertices)

def prism(radius=1, n=3, height=1):  # polygon prism
    '''
    Create a prism with the base of a regular polygon of number of sides n or a polygon defined by an array of vertices.
    A prism is an object extended from the regular polygon (on xy-plane) of sides n along z
    :param radius: radius of the polygon
    :param n: number of sides of the polygon
    :param height: extension height along z
    :return: prism in self.hull
    '''

    vertices_2d = __polygon_vertices(radius, n)
    bottom_1d = np.zeros(vertices_2d.shape[0]).reshape(-1, 1)
    top_1d = np.full(vertices_2d.shape[0], height).reshape(-1, 1)

    # Concatenate the 2D array and the reshaped 1D array
    bottom = np.hstack((vertices_2d, bottom_1d))
    top = np.hstack((vertices_2d, top_1d))

    return Hull(vertices=np.concatenate([bottom, top], axis = 0))

def torus(major_radius=2, minor_radius=1, division = 0.1):
    '''
    Create a torus (tire shape)
    :param major_radius: bigger radius from center to the ring of the tire
    :param minor_radius: ring radius of the tire
    :return: torus in self.hull
    '''
    major_sections = __get_num_edges(radius=major_radius, division = division)
    minor_sections = __get_num_edges(radius=minor_radius, division = division)
    obj = trimesh.creation.torus(major_radius=major_radius, minor_radius=minor_radius,
                                       major_sections=major_sections, minor_sections=minor_sections)
    return Hull(vertices=obj.vertices, faces=obj.faces)

def tetrahedron(a=1):
    '''
    Create a tetrahedron (4 triangular faces)
    :param a: radius
    :return: tetrahedron in self.hull
    '''
    # a tetrahedron of vertices (0, 0, 0), (a, 0, 0), (a/2, 3**0.5/2*a, 0) and (a/2, a/2/3**0.5, 2**0.5*a/3**0.5)
    vertices = np.array([
        [0, 0, 0],
        [a, 0, 0],
        [0.5 * a, 0.5 * 3 ** 0.5 * a, 0],
        [0.5 * a, 0.5 * a / 3 ** 0.5, 2 ** 0.5 / 3 ** 0.5 * a]
    ])

    return Hull(vertices = vertices)

def octahedron(a=1):
    '''
    Create an octahedron (8 faces)
    :param a: radius
    :return: octahedron in self.hull
    '''
    x = a / 2 ** 0.5
    vertices = np.array([
        [x, 0, 0],
        [-x, 0, 0],
        [0, x, 0],
        [0, -x, 0],
        [0, 0, x],
        [0, 0, -x]
    ])

    return Hull(vertices = vertices)

def dodecahedron(a=1):
    '''
    Create a dodecahedron (12 faces)
    :param a: radius
    :return: dodecahedron in self.hull
    '''
    phi = 0.5 * (1 + 5 ** 0.5)

    c1 = 1 / (3 - phi) ** 0.5
    x1 = a * c1

    c2 = 1 / phi * c1
    c3 = phi * c1
    x2 = a * c2
    x3 = a * c3

    vertices = np.array([
        [x1, x1, x1],
        [-x1, x1, x1],
        [x1, -x1, x1],
        [x1, x1, -x1],
        [-x1, -x1, x1],
        [-x1, x1, -x1],
        [x1, -x1, -x1],
        [-x1, -x1, -x1],
        [0, x2, x3],
        [0, -x2, x3],
        [0, x2, -x3],
        [0, -x2, -x3],
        [x2, x3, 0],
        [-x2, x3, 0],
        [x2, -x3, 0],
        [-x2, -x3, 0],
        [x3, 0, x2],
        [-x3, 0, x2],
        [x3, 0, -x2],
        [-x3, 0, -x2]
    ])

    return Hull(vertices = vertices)

def icosahedron(a=1):
    '''
    Create an icosahedron (20 faces)
    :param a: radius of a
    :return: icosahedron in self.hull
    '''
    phi = 0.5 * (1 + 5 ** 0.5)

    x1 = a * 0.5
    x2 = x1 * phi

    vertices = np.array([
        [x1, x2, 0],
        [-x1, x2, 0],
        [x1, -x2, 0],
        [-x1, -x2, 0],
        [x2, 0, x1],
        [-x2, 0, x1],
        [x2, 0, -x1],
        [-x2, 0, -x1],
        [0, x1, x2],
        [0, -x1, x2],
        [0, x1, -x2],
        [0, -x1, -x2]
    ])

    return Hull(vertices = vertices)

def extremes(hull):
    # Combine all vertices from all meshes
    vertices = hull.vertices

    # Compute extremes from the combined set of vertices
    min_x = np.min(vertices[:, 0])
    max_x = np.max(vertices[:, 0])
    min_y = np.min(vertices[:, 1])
    max_y = np.max(vertices[:, 1])
    min_z = np.min(vertices[:, 2])
    max_z = np.max(vertices[:, 2])

    return min_x, max_x, min_y, max_y, min_z, max_z

def plothull(hull, axis = True, origin_size = 0.05, axis_length = None):
    '''
    Simple hull display
    :param hull: trimesh object to display: Hull.hull
    :param axis: red: x; green: y; blue: z
            True: show with axis
            False: show without axis
    :param origin_size: size of the origin dot
    :param axis_length: length of the axis
            if None, the length is 1.5 times of the extreme vertices.
    :return: display the trimesh object
    '''
    obj = hull
    if axis:
        if axis_length is None:
           axis_length = max([abs(extreme) for extreme in extremes(hull)])*1.5

        axes = trimesh.creation.axis(origin_size=origin_size, axis_length=axis_length)

        obj = [hull, axes]

    # Combine the mesh and the coordinate axes
    scene = trimesh.Scene(obj)

    # Show the scene
    scene.show()

############### 2d shapes #################
class Shape2d:
    '''
    A 2d shape object
    '''
    def __init__(self, vertices):
        self.vertices = np.array(vertices)

    def extremes(self):
        '''
        show extreme points
        :return: leftmost, rightmost, top, bottom
        '''
        xs = self.vertices[:, 0]
        ys = self.vertices[:, 1]

        leftmost = xs.min()
        rightmost = xs.max()
        top = ys.max()
        bottom = ys.min()

        return leftmost, rightmost, top, bottom

    # rotate
    def __rotate(self, alpha):
        """
        Rotate an array of points in 3D space around the x, y, and z axes by the given angles.
        Points should be a numpy array of shape (n, 3).
        Angles should be in radians.
        """
        shape = Polygon(self.vertices)
        result = rotate(shape, angle=alpha, origin='centroid', use_radians=True)
        vertices = np.array(result.exterior.coords)
        return vertices

    # shift
    def __shift(self, displacement):
        return self.vertices + np.array([displacement])

    # make_closed/open
    def __toggle(self):
        if np.equal(self.vertices[0], self.vertices[-1]).all(): return self.vertices[:-1, :] # open the curve
        else: return np.concatenate((self.vertices, np.array([self.vertices[0, :]]))) # close the curve

    # stretch
    def __stretch(self, mag = np.array([1, 1])):
        return self.vertices*np.array(mag).reshape(1, 2)

    # shear
    def __shear(self, alpha, beta):
        # Create shear matrices
        Sx = np.array([[1, np.tan(alpha)], [0, 1]])
        Sy = np.array([[1, 0], [np.tan(beta), 1]])

        # Combine the shear matrices
        S = np.dot(Sx, Sy)

        # Apply the shear transformation to the vertices
        sheared_vertices = np.dot(self.vertices, S.T)

        return sheared_vertices

    # union 2d
    def __union(self, other):
        polygon1 = Polygon(self.vertices)
        polygon2 = Polygon(other.vertices)
        result = polygon1.union(polygon2)
        vertices = np.array(result.exterior.coords)
        return vertices

    # difference
    def __difference(self, other):
        polygon1 = Polygon(self.vertices)
        polygon2 = Polygon(other.vertices)
        result = polygon1.difference(polygon2)
        vertices = np.array(result.exterior.coords)
        return vertices

    # intersect
    def __intersection(self, other):
        polygon1 = Polygon(self.vertices)
        polygon2 = Polygon(other.vertices)
        result = polygon1.intersection(polygon2)
        vertices = np.array(result.exterior.coords)
        return vertices

    def __add__(self, other):
        if isinstance(other, self.__class__):
            vertices = self.__union(other)
            return Shape2d(vertices)
        elif isinstance(other, (np.ndarray, list, tuple)) and len(other) == 2:
            vertices = self.__shift(other)
            return Shape2d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape2d' and '{}'".format(type(other)))

    def __sub__(self, other):
        if isinstance(other, self.__class__):
            vertices = self.__difference(other)
            return Shape2d(vertices)
        elif isinstance(other, (np.ndarray, list, tuple)) and len(other) == 2:
            vertices = self.__shift(-np.array(other))
            return Shape2d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape2d' and '{}'".format(type(other)))

    def __mul__(self, other):
        if isinstance(other, self.__class__):
            vertices = self.__intersection(other)
            return Shape2d(vertices)
        elif isinstance(other, (np.ndarray, list, tuple)) and len(other) == 2:
            vertices = self.__stretch(other)
            return Shape2d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape2d' and '{}'".format(type(other)))

    def __matmul__(self, other):
        if isinstance(other, (float, int)):
            vertices = self.__rotate(other)
            return Shape2d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape2d' and '{}'".format(type(other)))

    def __mod__(self, other):
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 2:
            vertices = self.__shear(other[0], other[1])
            return Shape2d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape2d' and '{}'".format(type(other)))

    def __invert__(self):
        vertices = self.__toggle()
        return Shape2d(vertices)

    def __str__(self):
        l = self.vertices.tolist()
        s = 'Shape2d([' + str(l[0])
        for row in l[1:]:
            s += ',\n' + str(row)

        return s + '])'

    def __repr__(self):
        return self.__str__()




def plot2d(shape2d, marker_size = 5, grid = False):
    if isinstance(shape2d, Shape2d):
        shape2d = [shape2d]

    for shape in shape2d:
        # Extract x and y coordinates from the points array
        x_coords = shape.vertices[:, 0]
        y_coords = shape.vertices[:, 1]

        # Plotting
        plt.figure(figsize=(8, 8))  # Optional: Adjust figure size

        # Plot lines connecting the points
        plt.plot(x_coords, y_coords, color='red', linestyle='-', linewidth=1, marker='o', markersize=marker_size)

    # Adding labels and title
    plt.xlabel('X')
    plt.ylabel('Y')

    # Display the plot
    plt.grid(grid)  # Optional: Add grid
    plt.tight_layout()  # Optional: Improve spacing
    plt.show()

# regular polygon
def polygon(R = 1, n = 6):
    '''
    Create an n-side polygon with the radius of R on the xy plane. The center of the polygon is (0, 0, 0)
    :param R: radius of the polygon (default 1)
    :param n: number of sides (default 6)
    :return: (n, 2) 2darray of vertices
    '''
    # Generate the angles for the vertices
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    # Calculate the x and y coordinates
    x = R * np.cos(angles)
    y = R * np.sin(angles)
    vertices = np.array([x, y]).T
    return Shape2d(vertices)

# circle
def circle(R = 1, step = 0.1):
    '''
    Create a circle made up of n sides centered at (0, 0, 0) and the radius is R on the xy plane.
    :param R: radius of circle
    :param step: determines the number of sides:  n = max(circumference/step, 20)
    :return: (n, 2) 2darray of vertices
    '''
    n = int(max(np.ceil(2*np.pi*R/step), 20)) # number of sides of the circle
    return polygon(R, n)

# rectangle
def rectangle(a = 2, b = 1):
    '''
    Create a counterclockwise rectangle whose corners are:
    [
            [0, 0],
            [a, 0],
            [a, b],
            [0, b]
        ]
    :param a: horizontal side length
    :param b: vertical side length
    :return: (n, 2) 2darray of vertices
    '''

    vertices = np.array(
        [
            [0, 0],
            [a, 0],
            [a, b],
            [0, b]
        ]
    )
    return Shape2d(vertices)

# square
def square(a=1):
    '''
    Create a counterclockwise square whose corners are:
    [
            [0, 0],
            [a, 0],
            [a, a],
            [0, a]
        ]
    :param a: side length
    :return: (n, 2) 2darray of vertices
    '''

    return rectangle(a, a)

# parallelogram
def parallelogram(a=1, b=1, angle = np.pi/4):
    '''
    Create a counterclockwise parallelogram whose left bottom angle is angle in radians
    :param a: horizontal side length
    :param b: vertical side length
    :param angle: the left bottom corner
    :return: (n, 2) 2darray of vertices
    '''
    x_shift = b*np.cos(angle)
    y = b*np.sin(angle)
    vertices = np.array(
        [
            [0, 0],
            [a, 0],
            [a + x_shift, y],
            [x_shift, y]
        ]
    )

    return Shape2d(vertices)

# rhombus
def rhombus(a=1, angle = np.pi/4):
    '''
    Create a counterclockwise rhombus whose left bottom angle is angle in radians
    :param a: side length
    :param angle: the left bottom corner
    :return: (n, 3) 2darray of vertices
    '''
    x_shift = a*np.cos(angle)
    y = a*np.sin(angle)
    vertices = np.array(
        [
            [0, 0],
            [a, 0],
            [a + x_shift, y],
            [x_shift, y]
        ]
    )

    return Shape2d(vertices)

# isosceles trapezoid
def trapezoid(upper_side = 1, lower_side = 2, height = 1):
    '''
    Create an isosceles trapezoid whose upper_side < lower_side and the height is height. The left bottom corner is at the origin. The parallel sides are parallel to the x-axis on the xy plane.
    :param upper_side: upper side (shorter side)
    :param lower_side: lower side (longer side)
    :param height: height of the trapezoid
    :return: (n, 3) 2darray of vertices
    '''

    x_shift = (lower_side - upper_side)*0.5
    vertices = np.array(
        [
            [0, 0],
            [lower_side, 0],
            [lower_side - x_shift, height],
            [x_shift, height]
        ]
    )
    return Shape2d(vertices)

# pentagon
def pentagon(R = 1):
    '''
    Create a pentagon
    :param R: Radius of the pentagon
    :return: (n, 2) 2darray of vertices
    '''
    return polygon(R, 5)

# hexagon
def hexagon(R = 1):
    '''
    Create a hexagon
    :param R: Radius of the hexagon
    :return: (n, 2) 2darray of vertices
    '''
    return polygon(R, 6)

# heptagon
def heptagon(R = 1):
    '''
    Create a heptagon
    :param R: Radius of the heptagon
    :return: (n, 2) 2darray of vertices
    '''
    return polygon(R, 7)

# octagon
def octagon(R = 1):
    '''
    Create a octagon
    :param R: Radius of the octagon
    :return: (n, 2) 2darray of vertices
    '''
    return polygon(R, 8)

# nonagon
def nonagon(R = 1):
    '''
    Create a nonagon
    :param R: Radius of the nonagon
    :return: (n, 2) 2darray of vertices
    '''
    return polygon(R, 9)

# decagon
def decagon(R = 1):
    '''
    Create a decagon
    :param R: Radius of the decagon
    :return: (n, 2) 2darray of vertices
    '''
    return polygon(R, 10)

# ellipse
def ellipse(a = 2, b = 1, step = 0.1):
    '''
    Create an ellipse (oval) and center is at the origin
    :param a: semi-major axis
    :param b: semi-minor axis
    :param step: spacing between consecutive points.
    :return: (n, 2) 2darray of points
    '''
    # Calculate the circumference of the ellipse using an approximation formula
    h = ((a - b) ** 2) / ((a + b) ** 2)
    circumference = np.pi * (a + b) * (1 + 3 * h / (10 + np.sqrt(4 - 3 * h)))

    # Calculate the number of points needed
    num_points = int(np.ceil(circumference / step))

    # Generate the angles for the vertices
    angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)

    # Calculate the x and y coordinates of the ellipse
    x = a * np.cos(angles)
    y = b * np.sin(angles)

    vertices = np.array([x, y]).T
    return Shape2d(vertices)

# parabola
def parabola(a = 1, b = 2, w = 1, step = 0.1):
    '''
    Generate
    :param a: Parameter of the parabola
    :param b: x range of the parabola
    :param w: Width of the parabola
    :param step: step length of consecutive points along x
    :return: (n, 2) 2darray of points
    '''
    # Calculate the length of the parabola using an approximation formula
    # length = 8 * a / 3  # Approximate length of the parabola

    # Calculate the number of points needed
    num_points = int(np.ceil(w / step))

    # Generate the x coordinates for the vertices
    x = np.linspace(-b, b, num_points)

    # Calculate the corresponding y coordinates of the parabola
    y = x ** 2 / (4 * a)

    vertices = np.array([x, y]).T
    return Shape2d(vertices)

# hyperbola
def hyperbola(a = 1, b = 1, w = 1, step = 0.1):
    '''
    Create a hyperbola.
    :param a: Parameter along x
    :param b: Parameter along y
    :param w: width factor
    :param step: point step
    :return: (n, 2) 2darray of points
    '''

    # Calculate the approximate length of the hyperbola
    length = 4 * a * np.arcsinh(b / a)

    # Calculate the number of points needed
    num_points = int(np.ceil(length / step))

    final_angle = w*np.pi

    angles = np.linspace(-final_angle, final_angle, num_points)

    # Calculate the x and y coordinates of the hyperbola
    x = a * np.cosh(angles)
    y = b * np.sinh(angles)

    vertices = np.array([x, y]).T
    return Shape2d(vertices)

# spiral
def spiral(a = 1, b = 1, num_turns=2, step = 0.1):
    '''
    Create a spiral
    :param a: Parameter controlling growth in x-direction
    :param b: Parameter controlling growth in y-direction
    :param num_turns: number of turns of spiral
    :param step: point step
    :return: (n, 2) 2darray of points
    '''
    # Calculate the approximate length of the spiral
    length = num_turns * (np.sqrt(a ** 2 + b ** 2) + a * np.arcsinh(b / a))

    # Calculate the number of points needed
    num_points = int(np.ceil(length / step))

    # Generate the angles for the vertices
    angles = np.linspace(0, 2 * np.pi * num_turns, num_points)

    # Calculate the x and y coordinates of the spiral
    x = a * angles * np.cos(angles)
    y = b * angles * np.sin(angles)

    vertices = np.array([x, y])
    return Shape2d(vertices.T)

# cardioid
def cardioid(a = 1, step = 0.1):
    '''
    Create a cardioid shape
    :param a: size parameter
    :param step: point step
    :return: (n, 2) 2darray of points
    '''
    # Calculate the length of the cardioid using an approximation formula
    length = 8 * a

    # Calculate the number of points needed
    num_points = int(np.ceil(length / step))

    # Generate the angles for the vertices
    angles = np.linspace(0, 2 * np.pi, num_points)

    # Calculate the x and y coordinates of the cardioid
    x = a * (2 * np.cos(angles) - np.cos(2 * angles))
    y = a * (2 * np.sin(angles) - np.sin(2 * angles))

    vertices = np.array([x, y]).T
    return Shape2d(vertices)

# cycloid
def cycloid(a = 1, num_cycles=2, step = 0.1):
    '''
    Create a cycloid
    :param a: Radius of the generating circle
    :param num_cycles: Number of cycles of the cycloid
    :param step: point step

    :return: (n, 2) 2darray of points
    '''
    # Calculate the length of the cycloid using an approximation formula
    length = 8 * a * num_cycles

    # Calculate the number of points needed
    num_points = int(np.ceil(length / step))

    # Generate the parameter t for the vertices
    t = np.linspace(0, 2 * np.pi * num_cycles, num_points)

    # Calculate x and y coordinates of the cycloid
    x = a * (t - np.sin(t))
    y = a * (1 - np.cos(t))

    vertices = np.array([x, y])
    return Shape2d(vertices.T)

# astroid
def astroid(a = 1, step = 0.1):
    '''
    Create an astroid
    :param a: Radius of the astroid
    :param step: Point step
    :return: (n, 2) 2darray of points
    '''
    # Calculate the length of the astroid using an approximation formula
    length = 6 * a

    # Calculate the number of points needed
    num_points = int(np.ceil(length / step))

    # Generate the parameter t for the vertices
    t = np.linspace(0, 2 * np.pi, num_points)

    # Calculate x and y coordinates of the astroid
    x = a * np.cos(t) ** 3
    y = a * np.sin(t) ** 3

    vertices = np.array([x, y]).T
    return Shape2d(vertices)

# catenary
def catenary(a = 1, step = 0.1):
    '''
    Create a catenary
    :param a: size of the catenary
    :param step: point step
    :return: (n, 2) 2darray of points
    '''

    num_points = int(a/step)
    # Generate parameter t for the vertices
    t = np.linspace(-1, 1, num_points)

    # Calculate x and y coordinates of the catenary
    x = a * (np.cosh(t) - 1)
    y = a * np.sinh(t)

    vertices = np.array([x, y])
    return Shape2d(vertices.T)

# epicycloid
def epicycloid(R = 10, r = 1, step = 0.1):
    '''
    Create an epicycloid
    :param R: Radius of the fixed circle
    :param r: Radius of the rolling circle
    :param step: point step
    :return: (n, 2) 2darray of points
    '''
    # Calculate the length of one loop of the epicycloid using an approximation formula
    length = 4 * (R + r)

    # Calculate the number of points needed
    num_points = int(np.ceil(length / step))

    # Generate the parameter theta for the vertices
    theta = np.linspace(0, 2 * np.pi, num_points)

    # Calculate x and y coordinates of the epicycloid
    x = (R + r) * np.cos(theta) - r * np.cos((R + r) / r * theta)
    y = (R + r) * np.sin(theta) - r * np.sin((R + r) / r * theta)

    vertices = np.array([x, y]).T
    return Shape2d(vertices)

# hypocycloid
def hypocycloid(R = 6, r = 1, step = 0.1):
    '''
    Create a hypocycloid
    :param R: Radius of the fixed circle
    :param r: Radius of the rolling circle
    :param step: point step
    :return: (n, 2) 2darray of points
    '''
    # Calculate the length of one loop of the hypocycloid using an approximation formula
    length = 4 * (R - r)

    # Calculate the number of points needed
    num_points = int(np.ceil(length / step))

    # Generate the parameter theta for the vertices
    theta = np.linspace(0, 2 * np.pi, num_points)

    # Calculate x and y coordinates of the hypocycloid
    x = (R - r) * np.cos(theta) + r * np.cos((R - r) / r * theta)
    y = (R - r) * np.sin(theta) - r * np.sin((R - r) / r * theta)

    vertices = np.array([x, y])
    return Shape2d(vertices.T)

# sine curve
def sine(amplitude = 1, frequency = 1, x_start = 0, x_end = 2*np.pi,  step = 0.1):
    '''
    Create a sine curve
    :param amplitude: Amplitude of the curve
    :param frequency: Frequency of the curve
    :param x_start: beginning x
    :param x_end: ending x
    :param step: point step
    :return: (n, 2) 2darray of points
    '''
    # Calculate the total length of the sine curve
    x_dense = np.linspace(x_start, x_end, int((x_end - x_start)/step))
    y_dense = amplitude * np.sin(frequency * x_dense)

    # Compute the arc length using dense sampling points
    dx = np.diff(x_dense)
    dy = np.diff(y_dense)
    arc_length = np.sum(np.sqrt(dx ** 2 + dy ** 2))

    # Calculate the number of points needed based on the step size
    num_points = int(arc_length // step) + 1

    # Generate the vertices of the sine curve
    x_vertices = np.linspace(x_start, x_end, num_points)
    y_vertices = amplitude * np.sin(frequency * x_vertices)

    vertices = np.column_stack((x_vertices, y_vertices))
    return Shape2d(vertices)

# tangent curve
def tangent(x_start = -0.5*np.pi, x_end = 0.5*np.pi, step = 0.1,  height = 100):
    '''
    Create a tangent curve
    :param x_start: beginning x
    :param x_end: ending x
    :param step: point step
    :param height: height of the curve
    :return: (n, 2) 2darray of points
    '''
    x_dense = np.linspace(x_start, x_end, int((x_end - x_start)/step))
    y_dense = np.tan(x_dense)

    # Avoiding infinite values near asymptotes by masking
    mask = np.abs(y_dense) < height
    x_dense = x_dense[mask]
    y_dense = y_dense[mask]

    # Compute the arc length using dense sampling points
    dx = np.diff(x_dense)
    dy = np.diff(y_dense)
    arc_length = np.sum(np.sqrt(dx ** 2 + dy ** 2))

    # Calculate the number of points needed based on the step size
    num_points = int(arc_length // step) + 1

    # Generate the vertices of the tan curve
    x_vertices = np.linspace(x_dense[0], x_dense[-1], num_points)
    y_vertices = np.tan(x_vertices)

    # Avoiding infinite values near asymptotes again
    mask_vertices = np.abs(y_vertices) < height
    x_vertices = x_vertices[mask_vertices]
    y_vertices = y_vertices[mask_vertices]

    vertices = np.column_stack((x_vertices, y_vertices))
    return Shape2d(vertices)


class Shape3d:
    '''
    Convert a Shape2d (2d) object to a 3d object (Shape3d) on the xy-plane at z
    '''
    def __init__(self, shape2d, z = 0):
        if isinstance(shape2d, Shape2d):
            vertices = np.array(shape2d.vertices)
            self.vertices = np.hstack((vertices, np.full((vertices.shape[0], 1), z)))

        elif isinstance(shape2d, np.ndarray):
            if shape2d.shape[1] == 2:
                vertices = shape2d
                self.vertices = np.hstack((vertices, np.full((vertices.shape[0], 1), z)))
            elif shape2d.shape[1] == 3:
                self.vertices = shape2d

        elif isinstance(shape2d, (tuple, list)):
            if len(shape2d[0]) == 2:
                vertices = np.array(shape2d)
                self.vertices = np.hstack((vertices, np.full((vertices.shape[0], 1), z)))
            elif len(shape2d[0]) == 3:
                self.vertices = np.array(shape2d)

    # rotate
    def __rotation_matrix(self, axis, angle):
        """
        Return the rotation matrix for a rotation around the given axis by the given angle.
        Axis should be a string: 'x', 'y', or 'z'.
        Angle should be in radians.
        """
        if axis == 'x':
            return np.array([
                [1, 0, 0],
                [0, np.cos(angle), -np.sin(angle)],
                [0, np.sin(angle), np.cos(angle)]
            ])
        elif axis == 'y':
            return np.array([
                [np.cos(angle), 0, np.sin(angle)],
                [0, 1, 0],
                [-np.sin(angle), 0, np.cos(angle)]
            ])
        elif axis == 'z':
            return np.array([
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle), 0],
                [0, 0, 1]
            ])
        else:
            raise ValueError("The axis must be 'x', 'y' or 'z'.")

    def __rotate(self, alpha, beta, gamma):
        """
        Rotate an array of points in 3D space around the x, y, and z axes by the given angles.
        Points should be a numpy array of shape (n, 3).
        Angles should be in radians.
        """
        Rx = self.__rotation_matrix('x', alpha)
        Ry = self.__rotation_matrix('y', beta)
        Rz = self.__rotation_matrix('z', gamma)

        # Combine the rotations: first z, then y, then x
        R = np.dot(Rz, np.dot(Ry, Rx))

        return np.dot(self.vertices, R.T)

    def __shear(self, s_xy, s_xz, s_yz, s_yx, s_zx, s_zy):
        shear_matrix = np.array([
            [1, s_xy, s_xz],
            [s_yx, 1, s_yz],
            [s_zx, s_zy, 1]
        ])
        # Apply the shear matrix to each point
        return self.vertices.dot(shear_matrix.T)

    def __twist(self, other):
        '''
                Twist an object
                :param angle_x: Twist angle per unit length along the x-axis
                :param angle_y: Twist angle per unit length along the y-axis
                :param angle_z: Twist angle per unit length along the z-axis
                :return: twisted object in Shpe3d
                '''
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 3:
            angle_x, angle_y, angle_z = other
        else:
            angle_x, angle_y, angle_z = (0, 0, 0)

        # Get the vertices of the mesh
        vertices = self.vertices.copy()

        # Apply the twist transformation
        for i, vertex in enumerate(vertices):
            x, y, z = vertex

            # Calculate the twist angles based on the position along the axes
            twist_angle_x = angle_x * x
            twist_angle_y = angle_y * y
            twist_angle_z = angle_z * z

            # Rotation matrices for twisting around the x, y, and z axes
            rotation_matrix_x = np.array([
                [1, 0, 0],
                [0, np.cos(twist_angle_x), -np.sin(twist_angle_x)],
                [0, np.sin(twist_angle_x), np.cos(twist_angle_x)]
            ])

            rotation_matrix_y = np.array([
                [np.cos(twist_angle_y), 0, np.sin(twist_angle_y)],
                [0, 1, 0],
                [-np.sin(twist_angle_y), 0, np.cos(twist_angle_y)]
            ])

            rotation_matrix_z = np.array([
                [np.cos(twist_angle_z), -np.sin(twist_angle_z), 0],
                [np.sin(twist_angle_z), np.cos(twist_angle_z), 0],
                [0, 0, 1]
            ])

            # Combined rotation matrix
            combined_rotation_matrix = rotation_matrix_x @ rotation_matrix_y @ rotation_matrix_z

            # Apply the combined rotation to the vertex
            new_vertex = combined_rotation_matrix @ np.array([x, y, z])
            vertices[i] = new_vertex

        return vertices

    def __shift(self, displacement):
        return self.vertices + np.array([displacement])

    def __stretch(self, mag=np.array([1, 1, 1])):
        return self.vertices * np.array(mag).reshape(1, 3)

    # make_closed/open
    def __toggle(self):
        if np.equal(self.vertices[0, :], self.vertices[-1, :]).all():
            return self.vertices[:-1, :]  # open the curve
        else:
            return np.concatenate((self.vertices, np.array([self.vertices[0, :]])))  # close the curve

    def __add__(self, other):
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 3:
            vertices = self.__shift(other)
            return Shape3d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape3d' and '{}'".format(type(other)))

    def __sub__(self, other):
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 3:
            vertices = self.__shift(-np.array(other))
            return Shape3d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape3d' and '{}'".format(type(other)))

    def __mul__(self, other):
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 3:
            vertices = self.__stretch(other)
            return Shape3d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape3d' and '{}'".format(type(other)))

    def __matmul__(self, other):
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 3:
            vertices = self.__rotate(other[0], other[1], other[2])
            return Shape3d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape3d' and '{}'".format(type(other)))

    def __mod__(self, other):
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 6:
            vertices = self.__shear(other[0], other[1], other[2], other[3], other[4], other[5])
            return Shape3d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape3d' and '{}'".format(type(other)))

    def __xor__(self, other):
        if isinstance(other, (np.ndarray, list, tuple)) and len(other) == 3:
            vertices = self.__twist(other)
            return Shape3d(vertices)
        else:
            raise ValueError("Unsupported operand type for +: 'Shape3d' and '{}'".format(type(other)))

    def __invert__(self):
        vertices = self.__toggle()
        return Shape3d(vertices)

    def __str__(self):
        l = self.vertices.tolist()
        s = 'Shape3d([' + str(l[0])
        for row in l[1:]:
            s += ',\n' + str(row)

        return s + '])'

    def __repr__(self):
        return self.__str__()


def solenoid(cross_section, dz = 0.1, n_coils = 10):
    coils = []
    for i in range(n_coils):
        coil = Shape3d(cross_section, i*dz)
        coils.append(coil.vertices)

    return Shape3d(np.vstack(coils))

def plot3d(shape3d):
    if isinstance(shape3d, Shape3d): shape3d = [shape3d]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    for shape in shape3d:
        # Extract x, y, z coordinates
        x = shape.vertices[:, 0]
        y = shape.vertices[:, 1]
        z = shape.vertices[:, 2]

        # Plot the first 3D curve
        ax.plot(x, y, z, color='blue', marker='o')

    # Set labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Show the plot
    plt.show()


