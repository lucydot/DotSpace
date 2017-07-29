#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed # This is just for de-buggin'
import time

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
        basis = args["basis90"]
    elif orientation == 180:
        steps = 2
        basis = args["basis180"]
    elif orientation == 270:
        steps = 3
        basis = args["basis270"]
    elif orientation == 0:
        steps = 0
        basis = args["basis"]

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
    if args["cutoff"] is None:
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
   
    if args["no_correlation"] is not True:
        print ("Calculating correlation lengths...")
        data = correlate(args,patterns[-1])
    else: 
        data = False

    print ("Saving data...")
    save(args, data)
    
    print ("All done.")


if __name__=='__main__':
    philosophy = """
   -------------------------------------------------------------------
   |  Dotty for space code.                                          |
   |  Written by Lucy Whalley (github: lucydot) to analyse Dot-Space | 
   |  drawings by Richard Scott (richard-md-scott.tumblr.com).       |
   |                                                                 |
   |  Remember! Indexing starts from zero....                        |
   |                                                                 |
   |  "Whether obvious or subtle, visible or elusive, these          |
   |  symmetries are inherent and real."                             |
   |      Agnes Denes, Notes on a visual philosophy                  |
   -------------------------------------------------------------------
     """ 
    print (philosophy)

    parser=argparse.ArgumentParser(description="Dotty for space code",
                                   add_help=False)
    parser.add_argument("-i", "--identity", required=True, type=str,
                        help="identity for filenames etc")
    parser.add_argument("-w", "--width", required=True, type=int, 
                        help="width of pattern")
    parser.add_argument("-h", "--height", required=True, type=int, 
                        help="height of pattern")
    parser.add_argument("-b", "--basis", required=True, type=int, nargs='+',
                        help="repeating basis: 1=dot, 0=empty")
    parser.add_argument("-b90", "--basis90", type=int, nargs='+',
                        help="repeating basis: 1=dot, 0=empty")
    parser.add_argument("-b180", "--basis180", type=int, nargs='+',
                        help="repeating basis: 1=dot, 0=empty")
    parser.add_argument("-b270", "--basis270", type=int, nargs='+',
                        help="repeating basis: 1=dot, 0=empty")
    parser.add_argument("-d", "--defects", type=int, nargs='+', default=[],
                        help="position of defects "
                        "(in 1D, remember count from 0)")
    parser.add_argument("-a", "--args", action='help', help="print args")
    parser.add_argument("-nc", "--no_correlation", action="store_true",
                        help="if selected then correlation length will not be"
                             " calculated")
    parser.add_argument("-r", "--random", type=int, default=0,
                        help="number of defects to place randomly throughout"
                             "the pattern. Default is 0.")
    parser.add_argument("-x", "--cutoff", type=int, default=None,
                        help="largest horizontal and vertical distance for "
                             "the correlation calculation")
    parser.add_argument("-nn","--nonumba", action='store_true', default=False,
                        help="if selected, then numba will not be used")
    parser.add_argument("-o","--orientations", type=int, nargs='+', default=[0,90],
                        help="orientation of layers")
    args = vars(parser.parse_args())
    
    for item in args["defects"]: 
        if item > (args["width"]*args["height"])-1:
            raise ValueError("Defects specified beyond pattern size")
    if len(args["basis"]) > args["width"]*args["height"]:
        raise ValueError("Basis specified beyond pattern size")
    if args["random"]  > (args["width"]*args["height"]):
        raise ValueError("Number of defects cannot exceed pattern size")
    if (args["basis90"] == None) and (90 in args["orientations"]):
        args["basis90"] = args["basis"]
    if (args["basis180"] == None) and (180 in args["orientations"]):
        args["basis180"] = args["basis"]
    if (args["basis270"] == None) and (270 in args["orientations"]):
        args["basis270"] = args["basis"]
    
    main(args) 


