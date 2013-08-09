#!/usr/bin/env python

import sys
import math
from PIL import Image, ImageDraw
import numpy as np
import tsase
from tsase.data import *
import time


class TEM_image_comparison:
	def __init__(self, *args):
		np.seterr(all='ignore')
		#POSCAR file
		data = None
		try:
			data = tsase.io.read(args[0])
			if type(data) is not list:
				data = data
		except:
			print "failed to load", args[0]
			return

		#TEM image
		tem_image = None
		try:
			tem_image = Image.open(args[1])
		except:
			print "failed to load", args[1]
			return

		tem = 1-(np.array(tem_image.getdata())/255.0)
		self.tem_data = tem[:,0]

		#variable setup
		self.atoms     = data.get_positions()
		self.w, self.h = tem_image.size

		if len(args) == 2:
			self.bg        = 0
			self.intensity = 0
			self.sigma     = 0
			self.scale     = 0
			self.x_trans   = 0
			self.y_trans   = 0
			self.x_rot     = 0
			self.y_rot     = 0
			self.z_rot     = 0
		else: #only used in GUI, needed to follow gradient for a single variable. 
			self.bg        = args[2]
			self.intensity = args[3]
			self.sigma     = args[4]
			self.scale     = args[5]
			self.x_trans   = args[6]
			self.y_trans   = args[7]
			self.x_rot     = args[8]
			self.y_rot     = args[9]
			self.z_rot     = args[10]

		self.bgi_bool       = False
		self.bg_bool        = False
		self.intensity_bool = False
		self.sigma_bool     = False
		self.scale_bool     = False
		self.x_trans_bool   = False
		self.y_trans_bool   = False
		self.x_rot_bool     = False
		self.y_rot_bool     = False
		self.z_rot_bool     = False

		self.bg_slope        = None
		self.intensity_slope = None
		self.sigma_slope     = None
		self.scale_slope     = None
		self.x_trans_slope   = None
		self.y_trans_slope   = None
		self.x_rot_slope     = None
		self.y_rot_slope     = None
		self.z_rot_slope     = None

		self.bgi_o = None
		#self.score(self.bg, self.intensity, self.sigma, self.scale, self.x_trans, self.y_trans, self.x_rot, self.y_rot, self.z_rot)



###################################################################################################
# Score
###################################################################################################

	def score(self, bg, intensity, sigma, scale, x_trans, y_trans, x_rot, y_rot, z_rot):
        #setup variables
		w = self.w
		h = self.h
		atoms = np.array(self.atoms)

		#constants that do not need to be caclulated for every pixel (used in matrix creation functions)
		c1 = 2*sigma**2
		c2 = sigma**2
		c3 = sigma**3
		c4 = intensity/c2
		c6 = intensity/c1

		minx = min(atoms[:, 0])   
		miny = min(atoms[:, 1])
		minz = min(atoms[:, 2])
		maxx = max(atoms[:, 0])
		maxy = max(atoms[:, 1])
		maxz = max(atoms[:, 2])
		midx = (minx + maxx)/2
		midy = (miny + maxy)/2
		midz = (minz + maxz)/2

		#matrix creation function defenitions, all but fake_tem are used for finding slopes
		fake_tem = lambda index_array, x, y: (intensity*(np.e**(-(((index_array[0]- x)**2 + (index_array[1]-y)**2)/c1))))
		scale_helper = lambda index_array, x, y: c4 * (np.e**(-(((index_array[0]- x)**2 + (index_array[1]-y)**2)/c1))) * (((x/scale)*(index_array[0]-x)) + ((y/scale)*(index_array[1]-y)))
		sigma_helper = lambda index_array, x, y: intensity*(np.e**(-(((index_array[0]- x)**2 + (index_array[1]-y)**2)/c1)))*(((index_array[0]- x)**2 + (index_array[1]-y)**2)/c3)
		x_trans_helper = lambda index_array, x, y: (c4*((index_array[0] - x))*(np.e**(-(((index_array[0]- x)**2 + (index_array[1]-y)**2)/c1))))
		y_trans_helper = lambda index_array, x, y: (c4*((index_array[1] - y))*(np.e**(-(((index_array[0]- x)**2 + (index_array[1]-y)**2)/c1))))
		x_rot_helper = lambda index_array, x, y, x_o, z_o: c4 * (np.e**(-(((index_array[0]- x)**2 + (index_array[1]-y)**2)/c1))) * ((index_array[1] - y)*(scale*(-(y_o-midy)*np.sin(x_rot) + (z_o-midz)*np.cos(x_rot))))
		y_rot_helper = lambda index_array, x, y, y_o, z_o: -c4 * (np.e**(-(((index_array[0]- x)**2 + (index_array[1]-y)**2)/c1))) * ((index_array[0] - x)*(scale*( (x_o-midx)*np.sin(y_rot) + (z_o-midz)*np.cos(y_rot))))
		z_rot_helper = lambda index_array, x, y, x_o, y_o: -c6 * (np.e**(-(((index_array[0]- x)**2 + (index_array[1]-y)**2)/c1))) * ((2*(index_array[0] - x) * -scale*(-np.sin(z_rot)*(x_o-midx) + np.cos(z_rot)*(y_o-midy))) + (2*(index_array[1] - y) * -scale*(-np.cos(z_rot)*(x_o-midx) - np.sin(z_rot)*(y_o-midy))))
		
		#apply rotations
		atoms_original = np.array(atoms)
		atoms = self.rotate(atoms, x_rot, "x")
		atoms = self.rotate(atoms, y_rot, "y")
		atoms = self.rotate(atoms, z_rot, "z")

		#apply matrix creation functions to all atoms
		indices = np.indices((w,h))
		xyz = None

		for i in range(len(atoms)):
			x = ((atoms[i][0] + x_trans) * scale)
			y = ((atoms[i][1] + y_trans) * scale)
			x_o = atoms_original[i][0]
			y_o = atoms_original[i][1]
			z_o = atoms_original[i][2]

			if xyz == None:
				xyz = fake_tem(indices, x, y) + bg*255
				if self.scale_bool:
					scale_matrix = scale_helper(indices, x, y)
				if self.sigma_bool:
					sigma_matrix = sigma_helper(indices, x, y)
				if self.x_trans_bool: 
					x_trans_matrix = x_trans_helper(indices, x, y)
				if self.y_trans_bool:
					y_trans_matrix = y_trans_helper(indices, x, y)
				if self.x_rot_bool:
					x_rot_matrix = x_rot_helper(indices, x, y, x_o, z_o)
				if self.y_rot_bool:
					y_rot_matrix = y_rot_helper(indices, x, y, y_o, z_o)
				if self.z_rot_bool:
					z_rot_matrix = z_rot_helper(indices, x, y, x_o, y_o)
			else:
				xyz = xyz + fake_tem(indices, x, y)
				if self.scale_bool:
					scale_matrix = scale_matrix + scale_helper(indices, x, y)
				if self.sigma_bool:
					sigma_matrix = sigma_matrix + sigma_helper(indices, x, y)
				if self.x_trans_bool:
					x_trans_matrix = x_trans_matrix + x_trans_helper(indices, x, y)
				if self.y_trans_bool:
					y_trans_matrix = y_trans_matrix + y_trans_helper(indices, x, y)
				if self.x_rot_bool:
					x_rot_matrix = x_rot_matrix + x_rot_helper(indices, x, y, x_o, z_o)
				if self.y_rot_bool:
					y_rot_matrix = y_rot_matrix + y_rot_helper(indices, x, y, y_o, z_o)
				if self.z_rot_bool:
					z_rot_matrix = z_rot_matrix + z_rot_helper(indices, x, y, x_o, y_o)

		#score 
		tem = self.tem_data
		self.xyz = xyz #self.xyz is used to create the PIL image
		xyz = xyz.reshape(w*h)
		N = float(len(tem))
		newarray = (tem - (xyz/255.0))**2
		score = 1 - (np.sum(newarray)/N)
		print "score:", score
		self.xyz_data = xyz

		#gradients
		x = xyz/255.0
		t = tem

		if self.bgi_bool: #background and intensity (analytically solved for hence no slope)
			g = (N*np.sum(x*t) - np.sum(x)*np.sum(t))/(N*np.sum(x*x) - np.sum(x)*np.sum(x))
			b = np.sum(t-(g*x))/N
			self.bgi_o = [g, b]

		if self.bg_bool:
			self.bg_slope = (2/N)*np.sum((t-x))

		if self.intensity_bool:
			self.intensity_slope = (2/N)*np.sum((t-x)*(-x/intensity))

		if self.sigma_bool: #width of all atoms (sigma)
			sm = (sigma_matrix/255.0).reshape(w*h)
			self.sigma_slope = (2/N)*np.sum((t-x)*sm)

		if self.scale_bool: #scale of all atoms
			s = (scale_matrix/255.0).reshape(w*h)
			self.scale_slope = (2/N)*np.sum((t-x)*s)

		if self.x_trans_bool: #x translation of all atoms
			xt = (x_trans_matrix/255.0).reshape(w*h)
			self.x_trans_slope = (2/N)*np.sum((t-x)*xt)

		if self.y_trans_bool: #y translation of all atoms
			yt = (y_trans_matrix/255.0).reshape(w*h)
			self.y_trans_slope = (2/N)*np.sum((t-x)*yt)

		if self.x_rot_bool: #x rotation of all atoms
			xr = (x_rot_matrix/255.0).reshape(w*h)
			self.x_rot_slope = (2/N)*np.sum((t-x)*xr)

		if self.y_rot_bool: #y rotation of all atoms
			yr = (y_rot_matrix/255.0).reshape(w*h)
			self.y_rot_slope = (2/N)*np.sum((t-x)*yr)

		if self.z_rot_bool: #z rotation of all atoms
			zr = (z_rot_matrix/255.0).reshape(w*h)
			self.z_rot_slope = (2/N)*np.sum((t-x)*zr)

		return score

		

###################################################################################################
# slopes
###################################################################################################
	
	def slope_sigma(self, sigma): #working
		self.sigma_bool = True
		self.score(self.bg, self.intensity, sigma, self.scale, self.x_trans, self.y_trans, self.x_rot, self.y_rot, self.z_rot)
		self.sigma_bool = False
		return self.sigma_slope

	def slope_scale(self, scale): #working
		self.scale_bool = True
		self.score(self.bg, self.intensity, self.sigma, scale, self.x_trans, self.y_trans, self.x_rot, self.y_rot, self.z_rot)
		self.scale_bool = False
		return self.scale_slope

	def slope_xt(self, x_trans): #working
		self.x_trans_bool = True
		self.score(self.bg, self.intensity, self.sigma, self.scale, x_trans, self.y_trans, self.x_rot, self.y_rot, self.z_rot)
		self.x_trans_bool = False
		return self.x_trans_slope

	def slope_yt(self, y_trans): #working
		self.y_trans_bool = True
		self.score(self.bg, self.intensity, self.sigma, self.scale, self.x_trans, y_trans, self.x_rot, self.y_rot, self.z_rot)
		self.y_trans_bool = False
		return self.y_trans_slope

	def slope_xr(self, x_rot): #working
		self.x_rot_bool = True
		self.score(self.bg, self.intensity, self.sigma, self.scale, self.x_trans, self.y_trans, x_rot, self.y_rot, self.z_rot)
		self.x_rot_bool = False
		return self.x_rot_slope

	def slope_yr(self, y_rot): #working
		self.y_rot_bool = True
		self.score(self.bg, self.intensity, self.sigma, self.scale, self.x_trans, self.y_trans, self.x_rot, y_rot, self.z_rot)
		self.y_rot_bool = False
		return self.y_rot_slope

	def slope_zr(self, z_rot): #working
		self.z_rot_bool = True
		self.score(self.bg, self.intensity, self.sigma, self.scale, self.x_trans, self.y_trans, self.x_rot, self.y_rot, z_rot)
		self.z_rot_bool = False
		return self.z_rot_slope


		
###################################################################################################
# Follow Gradient Functions
###################################################################################################
	
	#working as intended, note background is between 0-1 not 0-255.
	def get_bgi(self, *args):
		self.bgi_bool = True
		self.score(0, 1, self.sigma, self.scale, self.x_trans, self.y_trans, self.x_rot, self.y_rot, self.z_rot)
		self.bgi_bool = False
		g, b = self.bgi_o
		return self.bgi_o



	def get_scale_o(self, *args):#working
		s = self.scale
		a = 5
		slope = 1
		while(math.fabs(slope)> .00001):
			slope = self.slope_scale(s)
			print "slope:", slope
			print ""
			s = s + (slope*a)
		return s



	def get_sigma_o(self, *args): #working
		w = self.sigma
		a = 50
		slope = 1
		while(math.fabs(slope) > .00001):
			slope = self.slope_sigma(w)
			w = w + (slope*a)
		return w



	def get_xt_o(self, *args): #working
		x = self.x_trans
		a = 250
		slope = 1
		while (math.fabs(slope) > .00001):
			slope = self.slope_xt(x)
			x = x + (slope*a)
		return x



	def get_yt_o(self, *args): #working, stop if translation is less than 1 pixel. 
		y = self.y_trans
		a = 250
		prev = -100
		slope = 1
		while (math.fabs(slope) > .00001):
			slope = self.slope_yt(y)
			prev = y
			y = y + (slope*a)
			if math.fabs(y-prev) < 1/self.scale:
				break
		return y



	def get_xr_o(self, *args): #working, might want to play with 'a' however
		xr = self.x_rot
		a = .1
		slope = 1
		while(math.fabs(slope) > .00001):
			slope = self.slope_xr(xr)
			xr = xr + (slope*a)
		return xr



	def get_yr_o(self, *args): #working, might want to play with 'a' aswell.
		yr = self.y_rot
		prev = -100
		a = 1
		slope = 1
		while(math.fabs(slope) > .00001):
			slope = self.slope_yr(yr)
			prev = yr
			print "slope:", slope
			print "y:", yr
			print ""
			yr = yr + (slope*a)
		return yr



	def get_zr_o(self, *args): #working, again 'a' needs to be tinkered with
		zr = self.z_rot
		prev = -100
		a = 1
		slope = 1
		while(math.fabs(slope) > .00001):
			slope = self.slope_zr(zr)
			prev = zr
			zr = zr + (slope*a)
		return zr



	def get_optimal_score(self, *args): #sigma always wants to 0 if initial guess is far from optimal position, NOTE: takes a very long time. 
		while(True):
			nums = self.get_bgi()
			x_trans = self.x_trans + 250*self.slope_xt(self.x_trans)
			y_trans = self.y_trans + 250*self.slope_yt(self.y_trans)
			sigma = self.sigma + 50*self.slope_sigma(self.sigma)
			scale = self.scale + 5*self.slope_scale(self.scale)
			x_rot = self.x_rot + 1*self.slope_xr(self.x_rot)
			y_rot = self.y_rot + 1*self.slope_yr(self.y_rot)
			z_rot = self.z_rot + 1*self.slope_zr(self.z_rot)
			
			self.x_trans = x_trans
			self.y_trans = y_trans
			self.sigma = sigma
			if self.sigma < 6:
				self.sigma = 6
			self.intesnity = nums[0]
			self.bg = nums[1]
			self.scale = scale
			self.x_rot = x_rot
			self.y_rot = y_rot
			self.z_rot = z_rot

			if math.fabs(self.sigma_slope) < .001 and math.fabs(self.scale_slope) < .001 and math.fabs(self.x_trans_slope) < .001 and math.fabs(self.y_trans_slope) < .001 and math.fabs(self.x_rot_slope) < .001 and math.fabs(self.y_rot_slope) < .001 and math.fabs(self.z_rot_slope) < .001:
				break

		return [self.bg, self.intensity, self.sigma, self.scale, self.x_trans, self.y_trans, self.x_rot, self.y_rot, self.z_rot]





###################################################################################################
# extras
###################################################################################################		
	def draw_image(self, *args): #creates image of the current state of the fake TEM image. 
		im = Image.new("RGBA", (self.w, self.h), (0, 0, 0, 0))
		draw = ImageDraw.Draw(im)
		self.score(self.bg, self.intensity, self.sigma, self.scale, self.x_trans, self.y_trans, self.x_rot, self.y_rot, self.z_rot)
		matrix = self.xyz
		for i in range(self.w):
			for j in range(self.h):
				value = int(matrix[i,j])
				draw.point((i,j), fill=(0,0,0, value))
		return im



	def rotate(self, atoms, theta, string):
		rot = None
		ct = math.cos(theta)
		st = math.sin(theta)
		if string == "x":
			rot = np.array([[1, 0, 0], [0, ct, -st], [0, st, ct]])
		if string == "y":
			rot = np.array([[ct, 0, st], [0, 1, 0], [-st, 0, ct]])
		if string == "z":
			rot = np.array([[ct, -st, 0], [st, ct, 0], [0, 0, 1]])
		
		minx = min(atoms[:, 0])   
		miny = min(atoms[:, 1])
		minz = min(atoms[:, 2])
		maxx = max(atoms[:, 0])
		maxy = max(atoms[:, 1])
		maxz = max(atoms[:, 2])
		midx = (minx + maxx)/2
		midy = (miny + maxy)/2
		midz = (minz + maxz)/2

		for i in atoms:
			i[0] = i[0] - midx
			i[1] = i[1] - midy
			i[2] = i[2] - midz

		atoms = atoms.dot(rot)
		
		for i in atoms:
			i[0] = i[0] + midx
			i[1] = i[1] + midy
			i[2] = i[2] + midz

		return atoms



	def get_score(self, vector):
		b  = vector[0]
		i  = vector[1]
		w  = vector[2]
		s  = vector[3]
		xt = vector[4]
		yt = vector[5]
		xr = vector[6]
		yr = vector[7]
		zr = vector[8]

		return self.score(b, i, w, s, xt, yt, xr, yr, zr)



	def get_gradient(self, vector):
		b  = vector[0]
		i  = vector[1]
		w  = vector[2]
		s  = vector[3]
		xt = vector[4]
		yt = vector[5]
		xr = vector[6]
		yr = vector[7]
		zr = vector[8]

		self.bgi_bool       = True
		self.bg_bool        = True
		self.intensity_bool = True
		self.sigma_bool     = True
		self.scale_bool     = True
		self.x_trans_bool   = True
		self.y_trans_bool   = True
		self.x_rot_bool     = True
		self.y_rot_bool     = True
		self.z_rot_bool     = True

		self.score(b, i, w, s, xt, yt, xr, yr, zr)

		self.bgi_bool       = False
		self.bg_bool        = False
		self.intensity_bool = False
		self.sigma_bool     = False
		self.scale_bool     = False
		self.x_trans_bool   = False
		self.y_trans_bool   = False
		self.x_rot_bool     = False
		self.y_rot_bool     = False
		self.z_rot_bool     = False

		return [self.bg_slope, self.intensity_slope, self.sigma_slope, self.scale_slope, self.x_trans_slope, self.y_trans_slope, self.x_rot_slope, self.y_rot_slope, self.z_rot_slope]



###################################################################################################
# Main
###################################################################################################
if __name__ == "__main__":
	print len(sys.argv)
	if len(sys.argv) == 3:
		TEM_image_comparison(sys.argv[1], sys.argv[2])
	else:
		print "imageEval.py needs a POSCAR file and a TEM image."
