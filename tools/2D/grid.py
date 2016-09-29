#!/usr/bin/env python
from __future__ import print_function
import sys
import csv
import itertools
import math

def dist(p1, p2):
  return math.sqrt( (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 )

class circle(object):
  def __init__(self, center, radius):
    self.center = center
    self.__act_center = [center[0] + 0.5, center[1] + 0.5]
    self.radius = radius
    self.indices = []
    
    for y in range(self.center[1] - self.radius, self.center[1] + self.radius + 2):
      for x in range(self.center[0] - self.radius, self.center[0] + self.radius + 2):
	if dist([self.center[0] + 0.5, self.center[1] + 0.5], [x, y]) <= radius:
	  self.indices.append([x, y])
	      
  def __repr__(self):
    return "circle"
  
  def __eq__(self, other):
    return self.center == other.center and self.radius == other.radius
    
  def distance_to(self, other):
    if type(other) is circle:
      distance = dist(self.center, other.center) 
      return distance - self.radius - other.radius if distance > self.radius + other.radius else 0
    elif type(other) is rectangle:
      #left of rectangle
      distance = 0
      if self.__act_center[0] < other.bottomleft[0] and other.bottomleft[1] <= self.__act_center[1] <= other.topright[1]:
	distance = other.bottomleft[0] - self.__act_center[0]
      #right
      elif self.__act_center[0] > other.topright[0] and other.bottomleft[1] <= self.__act_center[1] <= other.topright[1]:
	distance = self.__act_center[0] -  other.topright[0]
      #bottom
      elif self.__act_center[1] < other.bottomleft[1] and other.bottomleft[0] <= self.__act_center[0] <= other.topright[0]:
	distance = other.bottomleft[1] - self.__act_center[1]
      #top
      elif self.__act_center[1] > other.topright[1] and other.bottomleft[0] <= self.__act_center[0] <= other.topright[0]:
	distance = self.__act_center[1] - other.topright[1]
      #bottomleft
      elif self.__act_center[0] < other.bottomleft[0] and self.__act_center[1] < other.bottomleft[1]:
	distance = dist(self.__act_center, other.bottomleft)
      #bottomright
      elif self.__act_center[0] > other.topright[0] and self.__act_center[1] < other.bottomleft[1]:
	distance = dist(self.__act_center, [other.topright[0], other.bottomleft[1]])
      #topleft
      elif self.__act_center[0] < other.bottomleft[0] and self.__act_center[1] > other.topright[1]:
	distance = dist(self.__act_center, [other.bottomleft[0], other.topright[1]])
      #topright
      elif self.__act_center[0] > other.topright[0] and self.__act_center[1] > other.topright[1]:
	distance = dist(self.__act_center, other.topright)
      else:
	distance = 1
	
      return distance - self.radius if distance > self.radius else 0
	  
      
      
class rectangle(object):
  def __init__(self, bottomleft, topright):
    self.bottomleft = bottomleft
    self.topright = topright
    self.indices = itertools.product(range(bottomleft[0], topright[0] + 1), range(bottomleft[1], topright[1] + 1))
    
  def __repr__(self):
    return "rectangle"
  
  def __eq__(self, other):
    return self.bottomleft == other.bottomleft and self.topright == other.topright
    
  def distance_to(self, other):
    if type(other) is circle:
      return other.distance_to(self)
    elif type(other) is rectangle:
      distance = 0
            
      #left
      if self.topright[0] < other.bottomleft[0] and self.topright[1] >= other.bottomleft[1] and self.bottomleft[1] <= other.topright[1]:
	distance = other.bottomleft[0] - self.topright[0]
      #right
      elif self.bottomleft[0] > other.topright[0] and self.topright[1] >= other.bottomleft[1] and self.bottomleft[1] <= other.topright[1]:
	distance = self.bottomleft[0] - other.topright[0]
      #bottom
      elif self.topright[1] < other.bottomleft[1] and self.bottomleft[0] <= other.topright[0] and self.topright[0] >= other.bottomleft[0]:
	distance = other.bottomleft[1] - self.topright[1]
      #top
      elif self.bottomleft[1] > other.topright[1] and self.bottomleft[0] <= other.topright[0] and self.topright[0] >= other.bottomleft[0]:
	distance = self.bottomleft[1] - other.topright[1]
      #bottomleft
      elif self.topright[0] < other.bottomleft[0] and self.topright[1] < other.bottomleft[1]:
	distance = dist(self.topright, other.bottomleft)
      #bottomright
      elif self.bottomleft[0] > other.topright[0] and self.topright[1] < other.bottomleft[1]:
	distance = dist([self.bottomleft[0], self.topright[1]], [other.topright[0], other.bottomleft[1]])
      #topleft
      elif self.topright[0] < other.bottomleft[0] and self.bottomleft[1] > other.topright[1]:
	distance = dist([self.topright[0], self.bottomleft[1]], [other.bottomleft[0], other.topright[1]])
      #topright
      elif self.bottomleft[0] > other.topright[0] and self.bottomleft[1] > other.topright[1]:
	distance = dist(self.bottomleft, other.topright)
      else:
	distance = 0

      return distance - 1 if distance != 0 else 0

class zigzag(object):
  def __init__(self, bottomleft, width, distance):
    self.bottomleft = bottomleft
    self.width = width
    self.distance = distance
    self.indices = []

    start_i = bottomleft[0] if width > 0 else bottomleft[0] + abs(width) - 1

    if distance > 0:
      i_order = 1 if width > 0 else -1    

      width = abs(width)
    
      j = bottomleft[1]

      count = 0
      while(count < width + distance):
        for i in range(start_i, start_i + i_order * min(bottomleft[1] + distance - j + width, width), i_order):
          self.indices.append([i, j])
        
        count += 1
        start_i += i_order
        j += 1

    if distance < 0:
      i_order = 1 if width > 0 else -1    

      width = abs(width)
      abs_distance = abs(distance)
    
      j = bottomleft[1]

      count = 0
      while(count < width + abs_distance):
        for i in range(start_i, start_i + i_order * min(j - (bottomleft[1] - abs_distance - width), width), i_order):
          self.indices.append([i, j])
        
        count += 1
        start_i += i_order
        j -= 1		
    
class grid:
  def __init__(self, cols, rows, inverted = False, eps = 1e-14):
    self.__setval = 1 if not inverted else 0

    self.__data = [[1] * cols]    
    for j in range(1, rows - 1):
      self.__data.append([1] + [1 - self.__setval] * (cols - 2) + [1])
    self.__data.append([1] * cols)
    
    self.cols = cols
    self.rows = rows
    self.eps = 1e-14
    self.circles = []
    self.rectangles = []
    self.zigzags = []
    self.boundaries = []
    #left
    self.boundaries.append(rectangle([0, 0], [0, rows - 1]))
    #right
    self.boundaries.append(rectangle([cols - 1, 0], [cols - 1, rows - 1]))
    #bottom
    self.boundaries.append(rectangle([0, 0], [cols - 1, 0]))
    #top
    self.boundaries.append(rectangle([0, rows - 1], [cols - 1, rows - 1]))
  
  def __getitem__(self, tup):
    x, y = tup
    return self.__data[self.rows - 1 - y][x]
  
  def __setitem__(self, tup, val):
    x, y = tup
    self.__data[self.rows - 1 - y][x] = val
    
  def __repr__(self):
    return "grid()"
    
  def __str__(self):
    output = ""
    for y, row in enumerate(self.__data):
      output += str(row)
      if y != self.rows - 1:
	output += "\n"
	
    return output
  
  def __iter__(self):
    return iter(self.__data)
  
  def reset(self):
    self.circles = []
    self.rectangles = []
	
  def apply_shapes(self):
    for shape in itertools.chain(self.circles, itertools.chain(self.rectangles, self.zigzags)):
      self.set_indices(shape.indices)
    
  def write_to(self, file_path):
    with open(file_path, 'wb') as csvfile:
      gridwriter = csv.writer(csvfile, delimiter=',', quotechar='"')
	
      for row in self.__data:
	gridwriter.writerow(row)
	
  def is_viable(self, shape, distance, include_boundary):
    if type(shape) is zigzag and distance != 0:
      for other in self.zigzags:
        for idx in shape.indices:
          for idy in other.indices:
            if idx[0] == idy[0] and idx[1] == idy[1]:
              return False      

    else:
      for other in itertools.chain(self.circles, self.rectangles):
        if distance != 0 and shape.distance_to(other) < distance:
          return False
      
      if include_boundary:
        for bnd in self.boundaries:
          if distance != 0 and shape.distance_to(bnd) < distance:
            return False
    
    return True
  
  def set_indices(self, indices):
    for index in indices:
      x, y = index
      if 0 <= x < self.cols and 0 <= y < self.rows:
	      self[x, y] = self.__setval
	
  def add_shape(self, shape, distance = 0, include_boundary = False):
    if self.is_viable(shape, distance, include_boundary):
      if type(shape) is circle:
	      self.circles.append(shape)
      elif type(shape) is rectangle:
	      self.rectangles.append(shape)
      elif type(shape) is zigzag:
        self.zigzags.append(shape)      
      
      return True
   
    return False
  
if __name__ == "__main__":

  b = set([1])

  if b:
    print("blabla")
  else:
    print("bubub")

   
  cols = 30
  rows = 30
  g = grid(cols, rows, True)
  #print(g.add_shape(rectangle([2, 2], [5, 5])))
  #print(g.add_shape(rectangle([8, 2], [10, 5]), 1))
  #print(g.add_shape(circle([10, 10], 2)))
  #print(g.add_shape(rectangle([8, 16], [10, 23]), 1))
  #print(g.add_shape(rectangle([3, 8], [6, 23]), 1))
  #z = zigzag([10, 10], [20, 20], 4)
  #print(z.indices)
 # g.add_shape(z)
#
 # g.apply_shapes()
  #print(g)
 
