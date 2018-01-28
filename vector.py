from math import sqrt, acos, pi
from decimal import Decimal, getcontext
from collections import Iterator

getcontext().prec = 30
class Vector(object):

    CANNOT_NORMALIZED_ZERO_VECTOR_MSG = "Cannot normalize the zero vector"
    ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG = "Only defined in two/three dims"

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
#            self.coordinates = tuple(coordinates)
            self.dimension = len(coordinates)
            self.start = 0

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def __iter__(self):
        self.start = 0
        return iter(self)


    def next(self):
        if self.start >= len(self.coordinates):
            raise StopIteration
        next_ve = self.coordinates[self.start]
        self.start += 1
        return next_ve



    def plus(self, v):  # Plus Function
        new_coordinates = [x+y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def minus(self, v):  # Minus Function
        new_coordinates = [x-y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def times_scalar(self, c):  # Scalar Function
        new_coordinates = [Decimal(c) * x for x in self.coordinates]
#        new_coordinates = [c*x for x in self.coordinates]
        return Vector(new_coordinates)


    def magnitude(self):  # Calculate the Magnitude
        coordinates_squared = [x**2 for x in self.coordinates]
        return (sum(coordinates_squared)).sqrt()

    def normalized(self):  # Normalize the vector
        try:
            magnitude = self.magnitude()
            return self.times_scalar(Decimal('1.0')/magnitude)
#            return self.times_scalar(1.0/magnitude)

        except ZeroDivisionError:  # Error when zero vector
            raise Exception(self.CANNOT_NORMALIZED_ZERO_VECTOR_MSG)


    def dot(self, v):  # Dot Product
        return sum([x*y for x,y in zip(self.coordinates, v.coordinates)])


    def angle_with(self, v, in_degrees=False):  # Calculate angle
        try:
            u1 = self.normalized()
            u2 = v.normalized()
            angle_in_radians = acos(u1.dot(u2))

            if in_degrees:  # Different Return Result Format
                degrees_per_radian = 180./pi
                return angle_in_radians * degrees_per_radian
            else:
                return angle_in_radians

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZED_ZERO_VECTOR_MSG:
                raise Exception("Cannot compute an angle with the zero vector")
            else:
                raise e

    def project_parr(self, v):
        u1 = v.normalized()
        t = self.dot(u1)
        return u1.times_scalar(t)

    def project_cross(self, v):
        u1 = self.project_parr(v)
        u2 = self.minus(u1)
        return u2

    def is_parallel_to(self, v):
        return(self.is_zero() or v.is_zero() or
               self.angle_with(v) == 0 or
               self.angle_with(v) == pi)
    def is_orthogonal_to(self, v, tolerance=1e-10):
        return abs(self.dot(v)) < tolerance

    def is_zero(self, tolerance = 1e-10):
        return self.magnitude() < tolerance

    def component_paraller_to(self,basis):
        u = basis.normalized()
        weight = self.dot(u)
        return u.times_scalar(weight)

    def component_orthogonal_to(self,basis):
        try:
            projection = self.component_paraller_to(basis)
            return self.minus(projection)

        except Exception as e:
            raise e


    def cross(self, v):
        try:
            x_1, y_1, z_1 = self.coordinates
            x_2, y_2, z_2 = v.coordinates
            new_coordinates = [  y_1*z_2 - y_2*z_1,
                               -(x_1*z_2 - x_2*z_1),
                                 x_1*y_2 - x_2*y_1  ]
            return new_coordinates

        except ValueError as e:
            msg = str(e)
            if msg == 'need more than 2 values to unpadck':
                self_embedded_in_R3 = Vector(self.coordinates + ('0',))
                v_embedded_in_R3 = Vector(v.coordinates + ('0',))
                return self_embedded_in_R3.cross(v_embedded_in_R3)
            elif (msg == 'too many values to unpack' or
                  msg == 'need more than 1 value to unpack'):
                raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)
            else:
                raise e
