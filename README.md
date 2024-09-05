# 3D and 2D Model Creation Python Module

## Overview

This Python module provides a set of tools for creating, manipulating, and visualizing 3D and 2D models. It leverages the power of libraries like `trimesh`, `numpy`, `matplotlib`, and `shapely` to handle various geometric operations on 3D and 2D shapes.

## Installation

To use this module, ensure you have the following dependencies installed:

```bash
pip install numpy trimesh shapely matplotlib scipy
```

## Usage

### Classes and Functions


#### `Hull`
A class for 3D mesh representation and manipulation, inheriting from `trimesh.Trimesh`. Supports operations like union, difference, intersection, rotation, twisting, and shearing.

- **Constructor**:
  - `Hull(*args, **kwargs)`

- **Operators**:
  - `+` : Union of two Hull objects or translation by a vector.
            Hull1 + Hull2: Combines two Hull objects (Union).
            hull1 + (x, y, z) or [x, y, z] or 1darray([x, y, z]) translates the hull1 by the displacement of (x, y, z).

  - `-` : Difference between two Hull objects or translation by a negative vector.
            Hull1 - Hull2: Computes the difference between two Hull objects.
            hull1 - (x, y, z) or [x, y, z] or 1darray([x, y, z]) translates the hull1 by the displacement of -(x, y, z).

  - `*` : Intersection of two Hull objects.
            Hull1 * Hull2: Computes the intersection of two Hull objects.

  - `@` : Rotation by angles `[angle_x, angle_y, angle_z]`.
            Hull1 @ [angle_x, angle_y, angle_z]: Rotates a Hull object by specified angles.
            angle_x, angle_y, angle_z are rotation angles around x, y, z in radians

  - `^` : Twisting by angles `[angle_x, angle_y, angle_z]` per unit length.
            Hull1 ^ [angle_x, angle_y, angle_z]: Twists a Hull object by specified angles per unit length along the x, y, and z axes.
            angle_x, angle_y, angle_z are rotation angles per unit length around x, y, z in radians

  - `%` : Shearing by factors `[shear_zx, shear_yz, shear_zy, shear_xy, shear_yz, shear_xz]`.
            Hull1 % [shear_zx, shear_yz, shear_zy, shear_xy, shear_yz, shear_xz]: Shears a Hull object by specified shear factors.
            shear_zx, shear_yz, shear_zy, shear_xy, shear_yz, shear_xz are the shear factors (The distance a point moves due to shear divided by the perpendicular distance of a point from the invariant line.) on the specific planes respectively.
                e.g. shear_xy: shear around x-axis to the direction of y

- **Methods**:
  - `extremes(hull)`: Find out the extreme points (vertices) of a hull object (hull = Hull object). The output is a tuple of 6 elements to provide minimum extent of x, maximum extent of x, minimum extent of y, maximum extent of y, minimum extent of z and maximum extent of z (min_x, max_x, min_y, max_y, min_z, max_z) respectively.
  - `showx(Hull)`: Plot a hull object (Hull) with coordinate system. X-axis: red, Y-axis: green and Z-axis: blue.
  - `plothull(hull, axis = True, origin_size = 0.05, axis_length = None)`: Plot a hull with options. hull: the Hull object; axis: True: show axes. False: do not show axes; Axis colors: red: x; green: y; blue: z; origin_size: size of the origin dot; axis_length: length of the axes; if None, the length is 1.5 times of the extreme vertices.




#### `Shape2d`
Represents and handles 2D shapes and provides geometric transformations such as rotation, translation, union, difference, intersection, and shearing.

- **Constructor**:
  - `Shape2d(vertices)`

- **Operators**:
  - `+` : Union of two 2D shapes or translation by a vector.
  - `-` : Difference between two 2D shapes or translation by a negative vector.
  - `*` : Intersection of two 2D shapes or stretching.
  - `@` : Rotation by an angle.
  - `%` : Shearing by factors.
  - `~` : Toggle between open and closed shapes.

- **Methods**:
  - `plot2d(shape2d, marker_size=5, grid=False)`: Plot the 2D shape using `matplotlib`.

#### `Shape3d`
Transforms a 2D shape into a 3D object on the xy-plane at a given z coordinate. It forms a 2d curve/shape on a plane in the 3d space. Supports geometric transformations like rotation, translation, stretching, twisting, and shearing.

- **Constructor**:
  - `Shape3d(shape2d, z=0)`

- **Operators**:
  - `+` : Translation by a vector `[x, y, z]`.
  - `-` : Translation by a negative vector `[x, y, z]`.
  - `*` : Stretching by factors `[a, b, c]`. It stretches the shape a times along x, b times along y and z times along z.
  - `@` : Rotation by angles `[alpha, beta, gamma]`.
  - `%` : Shearing by factors of `[shear_zx, shear_yz, shear_zy, shear_xy, shear_yz, shear_xz]`.
  - `^` : Twisting by angles `[angle_x, angle_y, angle_z]`.
  - `~` : Toggle between open and closed shapes.

- **Methods**:
  - `plot3d(shape3d)`: Plot the 3D shape using `matplotlib`.

#### `solenoid`
The `solenoid` function creates a 3D solenoid shape by extruding a 2D cross-section along the z-axis in a helical manner.

- **Function Signature**:
  - `solenoid(cross_section, dz=0.1, n_coils=10)`

  **Parameters**:
  - `cross_section`: A 2D shape (e.g., `Shape2d`) that will be extruded to form the solenoid.
  - `dz`: The vertical distance between consecutive coils. Default is `0.1`.
  - `n_coils`: The number of coils in the solenoid. Default is `10`.

  **Returns**:
  - A `Shape3d` object representing the solenoid.

  **Example**:
  ```python
  from module_name import circle, solenoid, plot3d

  # Create a circular cross-section
  cross_section = circle(R=1)

  # Create a solenoid with the circular cross-section
  my_solenoid = solenoid(cross_section, dz=0.1, n_coils=20)

  # Plot the solenoid
  plot3d(my_solenoid)
  ```

### Creating Shapes

#### 3D Shapes
- **Cuboid**: `cuboid(a=1, b=1, c=1)`: Creates a cuboid of sides a, b and c, centered at the origin.
- **Cube**: `cube(a=1)`: Creates a cube of side a, centered at the origin.
- **Cylinder**: `cylinder(radius=1, height=1)`: Creates a cylinder of the radius of the cross sectional circle and the height centered at the origin. The cylinder extends along z.
- **Cone**: `cone(radius=1, height=1)`: Creates a cone of the radius of the cross sectional circle and the height. The center of the base circle lies at the origin. The cone extends along +z.
- **Ellipsoid**: `ellipsoid(a=1, b=1, c=1, division=0.1)`: Creates an ellipsoid as the equation of (x/a)**2 + (y/b)**2 + (z/c)**2 = 1. The center of the ellipsoid lies at the origin. The division is the spacing between adjacent griding points on x-y plane.
- **Sphere**: `sphere(radius=1, face_index=3)`: Creates a sphere of radius centered at the origin. The recommended face_index can be 3 (4**3 = 64 faces) or 4 (4**4 = 256 faces).
- **Pyramid**: `pyramid(radius=1, n=3, height=1, division=0.1)`: Creates a pyramid of the bottom of n sides of a regular polygon in radius and height along z. The top of the pyramid lies at the origin. The polygon base of the pyramid lies at (0, 0, height). The division is the spacing of z points.
- **Prism**: `prism(radius=1, n=3, height=1)`: Creates a prism with the base of a regular polygon of number of sides n. The center of the polygon lies at the origin and the prism extends along +z. A prism is an object extended from the regular polygon (on xy-plane) of sides n or a polygon defined by an array of vertices along z.
- **Torus**: `torus(major_radius=2, minor_radius=1, division=0.1)`: Creates a torus that is a tire-like shape with the major radius (major_raduis, the surrounding circle) and the minor_radius (minor_radius, the cross sectional circle). The division is the spacing of points along the circumference of the minor circle and the spacing between adjacent cross sectional circles. The center of the major circle lies at the origin.
- **Tetrahedron**: `tetrahedron(a=1)`: Creates a tetrahedron (4 faces) of radius of a, one vertex of which lies at the origin, one side lies along x, and one face lies on x-y plane, extending along +z.
- **Octahedron**: `octahedron(a=1)`: Creates an octahedron (8 faces) of radius of a whose center is at the origin and the vertices are at the x, y and z axis symmetric to the origin.
- **Dodecahedron**: `dodecahedron(a=1)`: Creates a dodecahedron (12 faces) of radius of a whose center is at the origin.
- **icosahedron**: `icosahedron(a=1)`: Creates an icosahedron (20 faces) of radius of a whose center is at the origin.

#### 2D Shapes
- **Polygon**: `polygon(R=1, n=6)`: Creates a polygon with the number of sides of n in the radius of R on the xy-plane. The center of the polygon is (0, 0, 0).
- **Circle**: `circle(R=1, step=0.1)`: Creates a circle made up of n sides centered at (0, 0, 0) and the radius is R on the xy-plane. The number of sides is determined by n = circumference/step but the minimum number of sides is 20.
- **Rectangle**: `rectangle(a=2, b=1)`: Creates a rectangle with the horizontal length of a and vertical length of b, whose corners are:
                                            [
                                                    [0, 0],
                                                    [a, 0],
                                                    [a, b],
                                                    [0, b]
                                                ]
- **Square**: `square(a=1)`: Creates a square of side length of a, whose corners are:
                                            [
                                                    [0, 0],
                                                    [a, 0],
                                                    [a, a],
                                                    [0, a]
                                                ]
- **Parallelogram**: `parallelogram(a=1, b=1, angle=np.pi/4)`: Creates a counterclockwise parallelogram whose left bottom angle is angle in radians with a being the horizontal length and b being the vertical length.
- **Rhombus**: `rhombus(a=1, angle=np.pi/4)`: Creates a counterclockwise rhombus whose left bottom angle is angle in radians and the side length is a.
- **Trapezoid**: `trapezoid(upper_side=1, lower_side=2, height=1)`: Creates an isosceles trapezoid whose upper_side is shorter than the lower_side and the height is height. The left bottom corner is at the origin. The parallel sides are parallel to the x-axis on the xy-plane.
- **Pentagon**: `pentagon(R=1)`: 5-side polygon centered at the origin in radius of R.
- **Hexagon**: `hexagon(R=1)`: 6-side polygon centered at the origin in radius of R.
- **Heptagon**: `heptagon(R=1)`: 7-side polygon centered at the origin in radius of R.
- **Octagon**: `octagon(R=1)`: 8-side polygon centered at the origin in radius of R.
- **Nonagon**: `nonagon(R=1)`: 9-side polygon centered at the origin in radius of R.
- **Decagon**: `decagon(R=1)`: 10-side polygon centered at the origin in radius of R.
- **Ellipse**: `ellipse(a=2, b=1, step=0.1)`: Creates an ellipse (oval) centered at the origin for the semi-major axis as a and the semi-minor axis as b, where step is the spacing between adjacent points.
- **Parabola**: `parabola(a=1, b=2, w=1, step=0.1)`: Creates a convex parabola whose vertex is at the origin for the equation is y = x**2/(4*a) and 2*b is the width over x-axis. step is the spacing between adjacent points.
- **Hyperbola**: `hyperbola(a=1, b=1, w=1, step=0.1)`: Creates a hyperbola open to the right, whose vertex is at (1, 0). The step is the spacing between adjacent points.
- **Spiral**: `spiral(a=1, b=1, num_turns=2, step=0.1)`: Creates a spiral beginning at the origin. The a is the growth along x, b is the growth along y and num_turns is the number of turns of the spiral. The step is the spacing between adjacent points.
- **Cardioid**: `cardioid(a=1, step=0.1)`: Creates a cardioid (heart shape) symmetric to the x-axis. The a is the size of the shape. The step is the spacing between adjacent points.
- **Cycloid**: `cycloid(a=1, num_cycles=2, step=0.1)`: Creates a cycloid along x from the origin. The a is the radius of the generating circle and num_cycles is the number of rounds of the circle to rotate. The step is the spacing between adjacent points.
- **Astroid**: `astroid(a=1, step=0.1)`: Creates an asteroid centered at the origin. It is a hypocycloid with four cusps. Specifically, it is the locus of a point on a circle as it rolls inside a fixed circle with four times the radius. The a is radius of the asteroid. The step is the spacing between adjacent points.
- **Catenary**: `catenary(a=1, step=0.1)`: Creates a catenary for the vertex at the origin to open to the right of the size of a. The step is the spacing between adjacent points.
- **Epicycloid**: `epicycloid(R=10, r=1, step=0.1)`: Creates an epicycloid centered at the origin for the radius of the fixed circle as a and that of the outside rolling circle as b. The outside rolling circle rolls on the fixed circle.The step is the spacing between adjacent points.
- **Hypocycloid**: `hypocycloid(R=6, r=1, step=0.1)`: Creates a hypocycloid centered at the origin for the radius of the fixed circle as a and that of the inside rolling circle as b. The inside rolling circle rolls on the fixed circle.The step is the spacing between adjacent points.
- **Sine Curve**: `sine(amplitude=1, frequency=1, x_start=0, x_end=2*np.pi, step=0.1)`: Creates a sine curve. The step is the spacing between adjacent points.
- **Tangent Curve**: `tangent(x_start=-0.5*np.pi, x_end=0.5*np.pi, step=0.1, height=100)`: Creates a tangent curve. The height is the height of the curve (the maximum of the curve). The step is the spacing between adjacent points.

### Example

#### Create and Display a 3D Sphere

```python
from module_name import sphere, plot3d

# Create a sphere with radius 2
my_sphere = sphere(radius=2)

# Plot the 3D sphere
plot3d(my_sphere)
```

#### Create and Display a 2D Circle

```python
from module_name import circle, plot2d

# Create a circle with radius 3
my_circle = circle(R=3)

# Plot the 2D circle
plot2d(my_circle)
```

---

This manual provides an overview of the module's capabilities. For more advanced operations, consult the docstrings provided in the module's code.
