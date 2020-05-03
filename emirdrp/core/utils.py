
import numpy
import math


def create_rot2d(angle):
    """Create 2D rotation matrix"""
    ca = math.cos(angle)
    sa = math.sin(angle)
    return numpy.array([[ca, -sa], [sa, ca]])