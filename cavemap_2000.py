#!/usr/bin/env python
# coding: utf-8

import csv
import math

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np


#%%
# helper parsing functions
def parse_float(string):
    if string == '':
        return None
    else:
        return float(string)

def parse_multi_float(string):
    """ If string format is float/float, returns a (float,float) tuple.
        Otherwise, returns a float. Throws an exception if neither works."""
    try:
        return parse_float(string)
    except ValueError:
        return tuple([float(x) for x in string.split('/')[:2]])

def csv_to_shots(filepath):
    """ Converts a CSV file to shot_info dicts.
    Doesn't do the work of verifying that values are sensible"""
    
    # read in data
    with open(filepath, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        labels = reader.__next__()
        rows = []
        for row in reader:
            rows.append(row)
            
    # use labels to make row to dict parser
    columns = {}
    for i, label in enumerate(labels):
        columns[label] = i
     
    string_fields = ['from', 'name', 'note']
    float_fields = ['distance']
    union_fields = ['azimuth', 'inclination', 'left', 'right', 'up', 'down']

    def parse_row(row):
        info = {}
        try:  # different treatment for different fields
            for field in string_fields:
                info[field] = row[columns[field]]
            for field in float_fields:
                info[field] = parse_float(row[columns[field]])
            for field in union_fields:
                info[field] = parse_multi_float(row[columns[field]])
        except ValueError:
            print(row, field)
            raise
        return info
        
    return list(map(parse_row, rows))
        


#%%
class LinePlot():
    def __init__(self, title='Cave', dist_units='feet', angle_tol=2.0):
        self.shots = []
        self.shot_names = [] # extra list for faster lookup
        self.title = title
        self.dist_units = dist_units
        self.angle_tol = angle_tol
        self.origin_name = None
        
        
    def add_shot(self, shot_info):
        """Adds a shot to the line plot. 
        Args:
            shot_info (dict):
                name (str): name of target station
                from (str): name of reference station
                distance (float): shot length
                azimuth (tuple(float, float) / float): foresight or fore and backsight
                inclination (tuple(float, float) / float): foresight or fore and backsight
                left (tuple(float, float) / float): distance or immediate and farthest distance
                right (tuple(float, float) / float): distance or immediate and farthest distance
                up (tuple(float, float) / float): distance or immediate and farthest distance
                down (tuple(float, float) / float): distance or immediate and farthest distance
                note (str): any comments
"""
        if len(self.shots) > 0:
            assert shot_info['from'] in self.shot_names
        for field in ['azimuth', 'inclination']:
            value = shot_info[field]
            if isinstance(value, tuple): 
                if field=='inclination':
                    val_1 = math.copysign(value[1], value[0])
                else:
                    val_1 = (value[1]+180)%360
                if abs(value[0]-val_1) > self.angle_tol:
                    name = shot_info['from'] + '-->' + shot_info['name']
                    print(f"WARNING: {value} in shot {name}  has a value of {abs(value[0]-val_1)} and is not within bounds of {self.angle_tol}")
        assert shot_info['distance'] > 0
        self.shots.append(shot_info)
        self.shot_names.append(shot_info['name'])
        
        
    def process(self):
        """ Modifies shots to prepare for plotting.
        Takes average of angle coords
        Adds new global position field of type tuple(x,y,z)"""

        def calc_pos(prev, shot): # TODO - put this elsewhere?
            """Calculates the absolute position of a new shot"""
            x0, y0, z0 = prev['position']
            dist, azi, incl = shot['distance'], shot['azimuth'], shot['inclination']
            x = x0 + dist * math.cos(incl*2*np.pi/360) * math.sin(azi*2*np.pi/360)
            y = y0 + dist * math.cos(incl*2*np.pi/360) * math.cos(azi*2*np.pi/360)
            z = z0 + dist * math.sin(incl*2*np.pi/360)
            return (x,y,z)


        # take angle averages where applicable
        for shot in self.shots:
            azi = shot['azimuth']
            if isinstance(azi, tuple):
                shot['azimuth'] = (azi[0] + (azi[1]+180)%360)/2.0
            incl = shot['inclination']
            if isinstance(incl, tuple):
                shot['inclination'] = (incl[0] + incl[1])/2.0
        
        unfinished = self.shots[:]
        done = []
        done_names = lambda: [shot['name'] for shot in done]
        
        # we arbitratily call the first point the origin
        first = unfinished.pop(0)
        first['position'] = calc_pos({'position':(0,0,0)}, first)
        self.origin_name = first['from']
        done.append(first)
        

        
        # look for shots we can now link
        # TODO - make this a pretty algorithm
        looped_twice = False
        while len(unfinished) > 0:
            for shot in unfinished:
                if shot['from'] in done_names():
                    prev = next(s for s in done if s['name']==shot['from'])
                    # calculate new global position
                    shot['position'] = calc_pos(prev, shot)
                    # finished, move between lists
                    done.append(shot)
                    unfinished.remove(shot)
                    looped_twice = False
                    continue
            if looped_twice:
                raise ValueError('Provided shots do not connect from origin')
            else:
                looped_twice = True
        self.shots = done
                    
        #TODO - check for unfilled values earlier on
        

    def plot(self, view='3d', angle=0):
        """ Plots a profile view at angle relative to North."""
        
        # prepare segments and labels to plot
        points = [shot['position'] for shot in self.shots]
        labels = [shot['name'] for shot in self.shots]
        # add origin
        points.append((0,0,0))
        labels.append(self.origin_name)
        segments = []
        for shot in self.shots:
            try:
                prev = next(s for s in self.shots if s['name']==shot['from'])
                p0 = prev['position']
            except StopIteration:
                if shot['from'] == self.origin_name:
                    p0 = 0, 0, 0
                else:
                    print(self.origin_name)
                    raise
            p1 = shot['position']
            segments.append(tuple(zip(p0,p1)))

        # prep            
        fig = plt.figure()
                
        if view == '3d':
            ax = fig.gca(projection='3d')
        elif view == 'plan':
            ax = plt.gca()
            ax.set_aspect('equal', adjustable='box')
            # reduce points and lines to 2d
            points = [p[:2] for p in points]
            segments = [s[:2] for s in segments]
        else:
            raise KeyError(f'Invalid view: {view}')
        
        ax.scatter(*zip(*points))
        for seg in segments:
            ax.plot(*seg)
            if view == '3d':
                ax.plot(*seg, lw=10, alpha=0.2)
        annotations = []
        for i, txt in enumerate(labels):
            annotations.append(ax.text(*points[i], txt))
        
        plt.show()


#%%
lineplot = LinePlot()
infos = csv_to_shots('cave.csv')
for info in infos:
    lineplot.add_shot(info)
lineplot.process()
lineplot.plot(view='3d')



