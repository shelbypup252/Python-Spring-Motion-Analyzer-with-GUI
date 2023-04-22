import math

def vector_components(magnitude, angle_degrees):
    """
    Calculate the magnitude and direction of a 2D vector.

    Args:
        x_component (float): The x-component of the vector.
        y_component (float): The y-component of the vector.

    Returns:
        tuple: A tuple containing the magnitude and direction of the vector.
               The magnitude is a float and the direction is in degrees.
    """
    angle_radians = math.radians(angle_degrees)
    x_component = magnitude * math.cos(angle_radians)
    y_component = magnitude * math.sin(angle_radians)
    return x_component, y_component

def vector_magnitude_and_direction(x_component, y_component):
    """
        Convert a 2D vector's x and y components into a vector.

        Args:
            x (float): The x-component of the vector.
            y (float): The y-component of the vector.

        Returns:
            tuple: A tuple containing the x and y components of the vector.
    """
    magnitude = math.sqrt(x_component**2 + y_component**2)
    direction_radians = math.atan2(y_component, x_component)
    direction_degrees = math.degrees(direction_radians)
    return magnitude, direction_degrees

def spherical_to_cartesian(radius, inclination, azimuth):
    """
    Convert spherical coordinates to Cartesian coordinates.

    Args:
        radius (float): The radial distance from the origin in meters.
        inclination (float): The angle between the positive z-axis and the point in radians.
        azimuth (float): The angle between the positive x-axis and the projection of the
                         point onto the xy-plane in radians.

    Returns:
        tuple: A tuple containing the x, y, and z components of the point in meters.
    """
    # Convert spherical coordinates to Cartesian coordinates
    x = radius * math.sin(inclination) * math.cos(azimuth)
    y = radius * math.sin(inclination) * math.sin(azimuth)
    z = radius * math.cos(inclination)

    return (x, y, z)

def cartesian_to_spherical(x, y, z):
    """
    Convert a 3D vector's cartesian coordinates to spherical coordinates.

    Args:
        x (float): The x-coordinate of the vector.
        y (float): The y-coordinate of the vector.
        z (float): The z-coordinate of the vector.

    Returns:
        tuple: A tuple containing the radius, inclination, and azimuth of the vector
               in spherical coordinates. The radius is a float, and the inclination
               and azimuth are in radians.
    """
    radius = math.sqrt(x**2 + y**2 + z**2)
    inclination = math.acos(z / radius)
    azimuth = math.atan2(y, x)
    return radius, inclination, azimuth

def momentum(mass, velocity):
    """
    Calculate the momentum of an object.

    Args:
        mass (float): The mass of the object in kilograms.
        velocity (float): The velocity of the object in meters per second.

    Returns:
        float: The momentum of the object in kilogram-meters per second.
    """
    momentum = mass * velocity
    return momentum

def moment_of_inertia(shape, axis, mass, radius=None, length=None):
    """
    Calculate the moment of inertia of a shape around a specified axis.

    Args:
        shape (str): The shape of the object. Possible values are 'solid_disc',
                     'hollow_disc', 'solid_sphere', 'hollow_sphere', and 'rod'.
        axis (str): The axis of rotation. Possible values are 'center', 'end',
                    and 'edge'. Only applicable to the 'rod' shape.
        mass (float): The mass of the object in kilograms.
        radius (float): The radius of the object in meters. Only applicable to
                        'solid_disc', 'hollow_disc', 'solid_sphere', and 'hollow_sphere'.
        length (float): The length of the object in meters. Only applicable to 'rod'.

    Returns:
        float: The moment of inertia of the object around the specified axis in
               kilogram-meters squared.
    """
    if shape == 'solid_disc':
        moment = (1 / 2) * mass * radius ** 2
    elif shape == 'hollow_disc':
        moment = (1 / 2) * mass * radius ** 2
    elif shape == 'solid_sphere':
        moment = (2 / 5) * mass * radius ** 2
    elif shape == 'hollow_sphere':
        moment = (2 / 3) * mass * radius ** 2
    elif shape == 'rod':
        if axis == 'center':
            moment = (1 / 12) * mass * length ** 2
        elif axis == 'end':
            moment = (1 / 3) * mass * length ** 2
        elif axis == 'edge':
            moment = (1 / 4) * mass * length ** 2
        else:
            raise ValueError("Invalid axis: '{}'".format(axis))
    else:
        raise ValueError("Invalid shape: '{}'".format(shape))

    return moment

def torque_from_moment(moment_of_inertia, angular_acceleration):
    """
    Calculate the torque on an object given its moment of inertia and angular acceleration.

    Args:
        moment_of_inertia (float): The moment of inertia of the object in
                                   kilograms * meter^2.
        angular_acceleration (float): The angular acceleration of the object in
                                      radians per second squared.

    Returns:
        float: The torque on the object in newton-meters.
    """
    torque = moment_of_inertia * angular_acceleration
    return torque

def torque(shape, axis, mass, radius=None, length=None, angular_acceleration=None):
    """
    Calculate the torque on an object given its shape, axis of rotation, mass,
    radius or length, and angular acceleration.

    Args:
        shape (str): The shape of the object. Possible values are 'solid_disc',
                     'hollow_disc', 'solid_sphere', 'hollow_sphere', and 'rod'.
        axis (str): The axis of rotation. Possible values are 'center', 'end',
                    and 'edge'. Only applicable to the 'rod' shape.
        mass (float): The mass of the object in kilograms.
        radius (float): The radius of the object in meters. Only applicable to
                        'solid_disc', 'hollow_disc', 'solid_sphere', and 'hollow_sphere'.
        length (float): The length of the object in meters. Only applicable to 'rod'.
        angular_acceleration (float): The angular acceleration of the object in
                                      radians per second squared.

    Returns:
        float: The torque on the object in newton-meters.
    """
    moment = moment_of_inertia(shape, axis, mass, radius, length)
    torque = moment * angular_acceleration
    return torque


def torque_from_vectors_cross(displacement, force):
    """
    Calculate the torque on an object given its displacement vector and force vector.

    Args:
        displacement (tuple): A tuple containing the x, y, and z components of the
                              displacement vector in meters.
        force (tuple): A tuple containing the x, y, and z components of the
                       force vector in newtons.

    Returns:
        tuple: A tuple containing the x, y, and z components of the torque vector
               in newton-meters.
    """
    # Calculate the torque using the cross product of the displacement and force vectors
    torque_x = displacement[1] * force[2] - displacement[2] * force[1]
    torque_y = displacement[2] * force[0] - displacement[0] * force[2]
    torque_z = displacement[0] * force[1] - displacement[1] * force[0]

    return (torque_x, torque_y, torque_z)



def spring_motion(mass, k, x0):
    """
    Calculates the amplitude, period, and equation of motion for an object attached to a spring.

    Args:
        mass (float): The mass of the object in kilograms.
        k (float): The spring constant in newtons per meter.
        x0 (float): The initial displacement of the object from the equilibrium position in meters.

    Returns:
        tuple: A tuple containing the amplitude, period, and equation of motion for the object.
    """
    omega = math.sqrt(k/mass)  # Calculate the angular frequency
    T = 2*math.pi/omega        # Calculate the period
    A = abs(x0)                # Calculate the amplitude

    # Construct the equation of motion string
    if x0 > 0:
        equation = f"f(x) = {A:.2f} * sin({omega:.2f}t)"
    elif x0 < 0:
        equation = f"f(x) = -{A:.2f} * sin({omega:.2f}t)"
    else:
        equation = f"f(x) = 0"

    return (A, T, equation)

