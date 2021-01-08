# -*- coding: utf-8 -*-
from math import *
import sys

# Ellipsoids
ELLIPSOID_GRS80 = 1 # GRS80
ELLIPSOID_WGS84 = 2 # WGS84

# Major axis radius and oblateness by ellipsoid
GEODETIC_DATUM = {
    ELLIPSOID_GRS80: [
        6378137.0,         # [GRS80]長軸半径
        1 / 298.257222101, # [GRS80]扁平率
    ],
    ELLIPSOID_WGS84: [
        6378137.0,         # [WGS84]長軸半径
        1 / 298.257223563, # [WGS84]扁平率
    ],
}

# Upper limit for iterative calculations
ITERATION_LIMIT = 1000

'''
Vincenty method(forward method)
Find the coordinates and azimuth angle of the end point from the coordinates (latitude and longitude), azimuth angle, and distance of the start point
:param lat: latitude
:param lon: longitude
:param azimuth: azimuth
:param distance: distance
:param ellipsoid: ellipsoid
:return: Coordinates of the endpoind and azimuth
'''
def vincenty_direct(lat, lon, azimuth, distance, ellipsoid=None):

    # Obtain the major axis radius (a) and flatness (ƒ) required for the calculation from the constants, and calculate the minor axis radius (b)
    # If ellipsoid is not specified, the value of GRS80 is used.
    a, ƒ = GEODETIC_DATUM.get(ellipsoid, GEODETIC_DATUM.get(ELLIPSOID_GRS80))
    b = (1 - ƒ) * a

    # Convert to radians (except distance)
    φ1 = radians(lat)
    λ1 = radians(lon)
    α1 = radians(azimuth)
    s = distance

    sinα1 = sin(α1)
    cosα1 = cos(α1)

    # latitude on the auxiliary sphere
    U1 = atan((1 - ƒ) * tan(φ1))

    sinU1 = sin(U1)
    cosU1 = cos(U1)
    tanU1 = tan(U1)

    σ1 = atan2(tanU1, cosα1)
    sinα = cosU1 * sinα1
    cos2α = 1 - sinα ** 2
    u2 = cos2α * (a ** 2 - b ** 2) / (b ** 2)
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))

    # Initialize σ with s/(b*A)
    σ = s / (b * A)

    # Iterate the following calculations until σ converges
    # Set an upper limit on the number of iterations, since convergence may not be possible depending on the location.
    for i in range(ITERATION_LIMIT):
        cos2σm = cos(2 * σ1 + σ)
        sinσ = sin(σ)
        cosσ = cos(σ)
        Δσ = B * sinσ * (cos2σm + B / 4 * (cosσ * (-1 + 2 * cos2σm ** 2) - B / 6 * cos2σm * (-3 + 4 * sinσ ** 2) * (-3 + 4 * cos2σm ** 2)))
        σʹ = σ
        σ = s / (b * A) + Δσ

        # If the deviation is less than .0000000001, break
        if abs(σ - σʹ) <= 1e-12:
            break
        else:
            # If the calculation does not converge, return None.
            return None

    # Once σ has converged to the desired accuracy, perform the following calculations
    x = sinU1 * sinσ - cosU1 * cosσ * cosα1
    φ2 = atan2(sinU1 * cosσ + cosU1 * sinσ * cosα1, (1 - ƒ) * sqrt(sinα ** 2 + x ** 2))
    λ = atan2(sinσ * sinα1, cosU1 * cosσ - sinU1 * sinσ * cosα1)
    C = ƒ / 16 * cos2α * (4 + ƒ * (4 - 3 * cos2α))
    L = λ - (1 - C) * ƒ * sinα * (σ + C * sinσ * (cos2σm + C * cosσ * (-1 + 2 * cos2σm ** 2)))
    λ2 = L + λ1

    α2 = atan2(sinα, -x) + pi

    return {
        'lat': degrees(φ2),     # latitude of the endpoint
        'lon': degrees(λ2),     # longitude  of the endpoint
        'azimuth': degrees(α2), # azimuth of the endpoint
    }


if __name__ == '__main__':
    args = sys.argv
    print(args)
    print("1st arg - latitude of the start point: " + args[1])
    print("2nd arg - longitude of the start point: " + args[2])
    print("3rd arg - azimuth: " + args[3])
    print("4th arg - distance: " + args[4])
    print("5th arg - ellipsoid: " + args[5])
    vincenty_direct(args[1], args[2], args[3], args[4], args[5])
