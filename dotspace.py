#!/usr/bin/env python

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed # This is just for de-buggin'
import time
import os

def timeit(method):

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print ('%r (%r, %r) %2.2f sec' % \
                (method.__name__, args, kw, te-ts))
        return result
    return timed

def defect_layer(args,defects):
    """ Creates defect layer"""
    # Place random defects on positions which have not been specified in 
    # "defects".
    for position in args["random_numbers"]:
        defects[position] = 1 # This may change later; 
        # currently marking defect position
    defects.shape = (args["height"], args["width"])
    return defects
 
def get_randoms(args):
    defects = np.zeros(args["width"]*args["height"], dtype=int)
    for position in args["defects"]:
        defects[position] = 1 # This may change later; 
        # currently marking defect position
    
    zeros = np.where(defects == 0)[0]
    if bool(args["random"]) is True:
        random_numbers = np.random.choice(zeros, args["random"], replace=False) 
    else:
        random_numbers = []
    return random_numbers,defects

def create(args, defects, orientation):
    """ Creates pattern with basis broken at defect points"""
    
    pattern = np.zeros(args["width"]*args["height"], dtype=int)
    
    if orientation == 90:
        steps = 1
        basis = args["basis90_digital"]
    elif orientation == 180:
        steps = 2
        basis = args["basis180_digital"]
    elif orientation == 270:
        steps = 3
        basis = args["basis270_digital"]
    elif orientation == 0:
        steps = 0
        basis = args["basis_digital"]

    basis_length = len(basis)
    
    defects = np.rot90(defects,steps)
    
    defects = np.ravel(defects)

    # Working in 1D
    counter = 0 
    for position in np.arange(len(pattern)):
        if defects[position] == 1:
            pattern[position] = basis[0] # re-start basis at defect
            counter = 1
        else:
            if basis[counter] == 1:
                pattern[position] = 1
            counter = ((counter + 1) % basis_length) # clock maths!
    
    # 2D to the eye
    pattern.shape = (args["height"], args["width"])
    pattern = np.rot90(pattern,4-steps)    
    return pattern

def overlay(pattern):
    """generator expression used"""
    
    length = pattern[0].shape[0] # We know they are square
    pattern_overlaid = np.zeros((length, length))
    
    for i in np.arange(length):
        for j in np.arange(length):
            if any(x[i,j] == 1 for x in pattern):
               pattern_overlaid[i,j] = 1
               
    return pattern_overlaid

def visualise(args, defects, patterns):
    """Matplotlib magic: Generates PNG files for defects, bottom layer, 
    top layer and overlaid. 1's (dots) are black, 0's (spaces) are white."""
    image = plt.imshow(defects, interpolation='nearest')
    image.set_cmap('Greys')
    plt.axis('off')
    plt.savefig(args["identity"]+"-"+"defects"+".png", 
                    bbox_inches='tight',dpi=500)
    plt.close()
    
    if np.all(patterns[-1]):
        print ("WARNING: YOU HAVE A FULLY DOTTED DOTSPACE") 
    image = plt.imshow(patterns[-1], interpolation='nearest')
    image.set_cmap('Greys')
    plt.axis('off')
    plt.savefig(args["identity"]+"-"+"overlaid"+".png", 
                    bbox_inches='tight',dpi=500)
    plt.close()


    for index,x in enumerate(args["orientations"]):
        image = plt.imshow(patterns[index], interpolation='nearest')
        image.set_cmap('Greys')
        plt.axis('off')
        plt.savefig(args["identity"]+"-layer_"+str(x)+".png", 
                    bbox_inches='tight',dpi=500)
        plt.close()

def jitterbug(pattern,X,Y):
    
    # Histogram: index is distance dx^2+dy^2 
    total = np.zeros((X-1)*(X-1)+(Y-1)*(Y-1)+1)
    count = np.zeros((X-1)*(X-1)+(Y-1)*(Y-1)+1)
    
    dX = np.arange(-X+1, X)
    dY = np.arange(-Y+1, Y)
    
    pattern[pattern == 0] = -1
    
    store=np.array([np.empty(4)])
    for x in np.arange(X):
        for y in np.arange(Y):
            for dx in dX:
                # Pesky negative indices could be calculated 
                # as [-1]==[last element in array] etc!
                if x+dx < 0 or x+dx >= pattern.shape[0]:
                    pass
                else:
                    for dy in dY:
                        if y+dy < 0 or y+dy >= pattern.shape[1]:
                            pass
                        else:    
                            # Check if it's been done t'other way around
                            if [x+dx,y+dy,x,y] in store:
                                pass
                            else:
                                # You may now proceed 
                                np.vstack((store,[x,y,dx,dy]))
                                product = pattern[x,y]*pattern[x+dx,y+dy]
                                # Update histograms
                                total[dx*dx+dy*dy]=total[dx*dx+dy*dy]+product
                                count[dx*dx+dy*dy]=count[dx*dx+dy*dy]+1
    return total, count
   
@timeit
def correlate(args, pattern):
    """Calculates radial correlation function of pattern."""
    
    # See Ziman, Models of disorder p.? 
    # Thanks to Andrew Goodwin and Jarvist Frost for discussion

    # Can't jit dictionaries
    if "cutoff" not in args:
        X = pattern.shape[0]
        Y = pattern.shape[1]
    else:
        X = args["cutoff"]+1
        Y = args["cutoff"]+1

    # Ugly, slow  - but working! -  loops
    # (It really is slow...TODO: Parallelise?)
    # Store points so that there's no double counting

    # Jit-terbug magic
    total, count = jitterbug(pattern,X,Y)
    
    # Fudge to stop nan error. Everything's so easy with Numpy.
    count = np.where(count==0,1,count)
    
    # Weight total by number of partners
    rho = np.divide(total,count)
    
    distance_squared = np.arange((X-1)*(X-1) + (Y-1)*(Y-1) + 1)
    distance = [np.sqrt(r2) for r2 in distance_squared]
    
    # Plot 
    fig, ax = plt.subplots()
    ax.stem(distance, rho, markerfmt=' ')
    plt.savefig(args["identity"]+"-corr.png", 
                bbox_inches='tight')
    plt.close()
    
    with open(args["identity"]+"-corr.txt","w") as f:
        f.write(str(np.vstack([np.multiply(distance,distance), 
                               distance, rho, count]).T))
        
    # TODO: calculate correlation length
    # dummy number, too ridiculous to be mistaken as truth...
    correlation = 301249135 
    data = { 'correlation' : correlation }
    print ("Correlation length is: "+ str(correlation))
    
    return data

def save(args, data):
    """ Saves inputs to file so that pattern can be reproduced """
    
    with open(args["identity"]+".txt", "w") as f:
        for keys,values in args.items():
            f.write(str(keys)+" : " + str(values))
            f.write('\n') 
        if data:
            for keys, values in data.items():
                f.write(str(keys)+" : " + str(values))
                f.write('\n')

def main(args):
    os.mkdir(args["identity"])
    os.chdir(args["identity"])
    
    print ("Dot-Space ID: " + args['identity'])
    
    print ("Creating defect layer")
    args["random_numbers"], defect_initial = get_randoms(args)
    defects = defect_layer(args,defect_initial)
    
    patterns = []
    for x in args["orientations"]:
        print ("Creating layer %i" % (x))
        patterns.append(create(args, defects, x))
    
    if args["width"] == args["height"]: 
        print ("Overlaying patterns...")
        patterns.append(overlay(patterns))
    
    print ("Plotting visual representation of pattern...")
    visualise(args, defects, patterns)
   
    if args["correlation"] == True:
        print ("Calculating correlation lengths...")
        data = correlate(args,patterns[-1])
    else: 
        data = False

    print ("Saving data...")
    save(args, data)
    
    print ("All done.")

def parse_basis_input(line,name):
        args[name+"_digital"] = []
        digits_list = list(map(int,line.split("=")[1].strip()))
        for x in digits_list:
            args[name+"_digital"] = args[name+"_digital"] + [1] + [0]*(x-1)
        args[name] = line.split("=")[1].strip()

def import_input(input_filename):

    with open(input_filename) as f:
        for line in f:
            if line.startswith("#"):
                pass
            else:
                if line.split("=")[0].strip() == "IDENTITY":
                    args["identity"] = line.split("=")[1].strip()
                if line.split("=")[0].strip() == "WIDTH":    
                    args["width"] = int(line.split("=")[1].strip())
                if line.split("=")[0].strip() == "HEIGHT":           
                    args["height"] = int(line.split("=")[1].strip())
                if line.split("=")[0].strip() == "BASIS":
                    parse_basis_input(line,"basis")
                if line.split("=")[0].strip() == "BASIS90":
                    parse_basis_input(line,"basis90")
                if line.split("=")[0].strip() == "BASIS180":
                    parse_basis_input(line, "basis180")
                if line.split("=")[0].strip() == "BASIS270":
                    parse_basis_input(line, "basis270")
                if line.split("=")[0].strip() == "DEFECTS":
                    if "x" in line.split("=")[1]:
            # WIP: not sure how to specify input
                        args[defects] = []
                        entries = line.split("=")[1].split()
                        
                        for e in entries:
                            for i in range(int(e.split("x")[0])):
                                args[name].append(int(e.split("x")[1]))
                    else:
                        args["defects"] = [int(x) for x in line.split("=")[1].split()]
                if line.split("=")[0].strip() == "CORRELATION":
                    args["correlation"] = line.split("=")[1].strip().lower() in ("yes","true")
                if line.split("=")[0].strip() == "RANDOM_DEFECTS":
                    args["random"] = int(line.split("=")[1].strip())
                if line.split("=")[0].strip() == "CUTOFF":
                    args["cutoff"] = int(line.split("=")[1].strip())
                if line.split("=")[0].strip() == "ULAM":
                    args["ulam"] = line.split("=")[1].strip().lower() in ("yes","true")
                if line.split("=")[0].strip() == "ORIENTATIONS":
                    args["orientations"] = [int(x) for x in line.split("=")[1].split()]
    return args

if __name__=='__main__':
    philosophy = """
   -------------------------------------------------------------------
   |  Dotty for space code.                                          |
   |  Written by Lucy Whalley (github: lucydot) to test-run and      |
   |   analyse the Dot-Space drawings create by Richard Scott         | 
   |  (richard-md-scott.tumblr.com).                                 |
   |                                                                 |
   |  Remember! Indexing starts from zero....                        |
   |                                                                 |
   |  "Whether obvious or subtle, visible or elusive, these          |
   |  symmetries are inherent and real."                             |
   |      Agnes Denes, Notes on a visual philosophy                  |
   -------------------------------------------------------------------
     """ 
    if len(sys.argv) > 1:
    	input_file = sys.argv[1]
    else:
    	input_file = "INPUT"

    print (philosophy)
    args = dict()        
    args = import_input(input_file)
    if 'random' not in args:
        args["random"] = 0
    if 'correlation' not in args:
        args["correlation"] = False
    if 'cutoff' not in args:
        args["cutoff"] = None
    if 'defects' not in args:
        args["defects"] = []
    if 'identity' not in args:
        raise ValueError("You must supply identity")
    if 'basis' not in args:
        raise ValueError("You must supply a basis")
    if 'width' not in args:
        raise ValueError("You must supply a width")
    if 'height' not in args:
        raise ValueError("You must supply a height")
    for item in args["defects"]: 
        if item > (args["width"]*args["height"])-1:
            raise ValueError("Defects specified beyond pattern size")
    if len(args["basis"]) > args["width"]*args["height"]:
        raise ValueError("Basis specified beyond pattern size")
    if args["random"]  > (args["width"]*args["height"]):
        raise ValueError("Number of defects cannot exceed pattern size")
    if ("basis90" not in args) and (90 in args["orientations"]):
        args["basis90"] = args["basis"]
    if ("basis180" not in args) and (180 in args["orientations"]):
        args["basis180"] = args["basis"]
    if ("basis270" not in args) and (270 in args["orientations"]):
        args["basis270"] = args["basis"]
       
    main(args) 


