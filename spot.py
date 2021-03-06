
import numpy as np



##################################################
# Auxiliary functions

#Great circle distance, i.e. angle between values on top of a sphere
def great_circle_distance(lon1, lon2, col1, col2):

    lat1 = np.pi/2 - col1
    lat2 = np.pi/2 - col2

    dlon = np.abs(lon2 - lon1)
    dlat = np.abs(lat2 - lat1)


    xx = (np.cos(lat2)*np.sin(dlon))**2 + (np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(dlon))**2
    yy = np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(dlon)

    return np.arctan2(np.sqrt(xx), yy)


#Spot class that defines limits/shape for the surface emission
class Spot:

    star_time = 0.0 #internal time in the surface

    def __init__(self, colat, rho, angvel):

        self.colat  = np.deg2rad(colat)
        self.rho    = np.deg2rad(rho)
        #self.angvel = freq*2.0*np.pi
        self.angvel = angvel


    def circular_spot(self, rot_phi, phi, theta):
        #ang_distance = great_circle_distance(rot_phi, phi, self.colat, theta)
        ang_distance = great_circle_distance(phi, rot_phi, self.colat, theta)

        if ang_distance <= self.rho:
            return True
        else:
            return False

    def hit(self, coords):

        time  = coords[0] #note that this is negative
        theta = coords[1]
        phi   = coords[2]

        phi_rot = (self.star_time + time)*self.angvel

        inside = self.circular_spot(phi_rot, phi, theta)
        #inside = self.circular_spot(0.0, phi_rot, theta)

        #If we are inside, lets construct an emission class
        if inside:
            return True

        else: #not inside, hence return empty emission
            return False






