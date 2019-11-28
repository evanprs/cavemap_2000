#!/usr/bin/env python
# coding: utf-8

import argparse
import csv
import math

import matplotlib.pyplot as plt
import numpy as np

import helpers

parser = argparse.ArgumentParser(description='cavemap 2000')
parser.add_argument('files', nargs='*', type=str)
parser.add_argument('--plan', action='store_true', default=False,
                    dest='plan',
                    help='Plot a plan view')
parser.add_argument('--profile', action='store_true', default=False,
                    dest='profile',
                    help='Plot a profile view')
parser.add_argument('--flat', action='store_true', default=False,
                    dest='flat',
                    help='Plot a flat profile view')
parser.add_argument('--3d', action='store_true', default=False,
                    dest='three_d',
                    help='Plot a 3d view')

args = parser.parse_args()
paths = args.files
path = paths[0]
plan = args.plan
profile = args.profile
flat = args.flat
three_d = args.three_d

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
    """ Object for storing and plotting survey data.
    Args:
        title (str): Name of the cave
        dist_units (str): Distance unit label used in survey
        angle_tol (float): Allowable difference between fore and backsight
    """

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
                    val_1 = -value[1]
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

        def calc_pos_flat(prev, shot): # TODO - put this elsewhere?
            """Calculates the position of a new shot in a flattened profile view"""
            w0, z0 = prev['flat_position']
            dist, incl = shot['distance'], shot['inclination']
            w = w0 + dist * math.cos(incl*2*np.pi/360)
            z = z0 + dist * math.sin(incl*2*np.pi/360)
            return (w,z)



        # take angle averages where applicable
        for shot in self.shots:
            azi = shot['azimuth']
            if isinstance(azi, tuple):
                shot['azimuth'] = (azi[0] + (azi[1]+180)%360)/2.0
            incl = shot['inclination']
            if isinstance(incl, tuple):
                shot['inclination'] = (incl[0] - incl[1])/2.0

        unfinished = self.shots[:]
        done = []
        done_names = lambda: [shot['name'] for shot in done]

        # we arbitratily call the first point the origin
        first = unfinished.pop(0)
        first['position'] = calc_pos({'position':(0,0,0)}, first)
        first['flat_position'] = calc_pos_flat({'flat_position':(0,0)}, first)
        self.origin_name = first['from']
        done.append(first)


        # look for shots we can now link
        # TODO - check for unfilled values earlier on
        # TODO - make this a pretty algorithm
        looped_twice = False
        while len(unfinished) > 0:
            for shot in unfinished:
                if shot['from'] in done_names():
                    prev = next(s for s in done if s['name']==shot['from'])
                    # calculate new global position
                    shot['position'] = calc_pos(prev, shot)
                    shot['flat_position'] = calc_pos_flat(prev, shot)
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


    def plot(self, view='3d', angle=0):
        """ Plots a cave.
        Args:
            view (str): 3d, plan, profile, or flat_profile
            angle (float): if profile is selected, plot at angle relative to N
            """
        # TODO - use angle

        # prepare segments and labels to plot
        pos_type = 'flat_position' if view == 'flat_profile' else 'position'
        points = [shot[pos_type] for shot in self.shots]
        labels = [shot['name'] for shot in self.shots]
        # add origin
        points.append((0, 0, 0))
        labels.append(self.origin_name)


        segments = []
        for shot in self.shots:
            try:
                prev = next(s for s in self.shots if s['name']==shot['from'])
                p0 = prev[pos_type]
            except StopIteration:
                if shot['from'] == self.origin_name:
                    p0 = (0, 0, 0)
                else:
                    print(self.origin_name)
                    raise
            p1 = shot[pos_type]
            segments.append(tuple(zip(p0,p1)))

        if view in ['plan', 'flat_profile']:
            # reduce points and lines to 2d. for flat_prof this just trims the origin point
            points = [p[:2] for p in points]
            segments = [s[:2] for s in segments]
        elif view == 'profile':
            # reduce points and lines to 2d
            points = [p[1:] for p in points]
            segments = [s[1:] for s in segments]
        elif view == '3d':
            pass
        else:
            raise KeyError(f'Invalid view: {view}')

        # prepare axes
        fig = plt.figure()
        if view == '3d':
            ax = fig.gca(projection='3d')
        else:
            ax = fig.gca()
            ax.set_aspect('equal', adjustable='box')
            plt.axis('off')

        # create plot
        ax.scatter(*zip(*points), marker='^')
        for seg in segments:
            ax.plot(*seg, c='black')
            if view == '3d':
                ax.plot(*seg, lw=10, alpha=0.2)
                helpers.set_axes_equal(ax)
        annotations = []
        for i, txt in enumerate(labels):
            annotations.append(ax.text(*points[i], txt))
        ax.set_title(self.title)
        plt.show()


#%%
lineplot = LinePlot()
infos = csv_to_shots(path)
for info in infos:
    lineplot.add_shot(info)
lineplot.process()
if profile:
    lineplot.plot(view='profile')
if plan:
    lineplot.plot(view='plan')
if flat:
    lineplot.plot(view='flat_profile')
if three_d:
    lineplot.plot(view='3d')




