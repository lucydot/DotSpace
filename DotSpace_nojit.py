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

def create(args):
    """ Creates pattern with basis broken at defect points"""
    
    basis_length = len(args["basis"])
    pattern = np.zeros(args["width"]*args["height"], dtype=int)
    defects = np.zeros(args["width"]*args["height"], dtype=int)
    
    for position in args["defects"]:
        defects[position] = 1 # This may change later; 
        # currently marking defect position
    
    # Place random defects on positions which have not been specified in 
    # "defects".
    zeros = np.where(defects == 0)[0]
    random_numbers = np.random.choice(zeros, args["random"], replace=False) 
    args["random_numbers"] = random_numbers

    for position in random_numbers:
        defects[position] = 1 # This may change later; 
        # currently marking defect position
        
    # Working in 1D
    counter = 0 
    for position in np.arange(len(pattern)):
        if defects[position] == 1:
            pattern[position] = args["basis"][0] # re-start basis at defect
            counter = 1
        else:
            if args["basis"][counter] == 1:
                pattern[position] = 1
            counter = ((counter + 1) % basis_length) # clock maths!
    
    # 2D to the eye
    pattern.shape = (args["height"], args["width"])
    defects.shape = (args["height"], args["width"])
    
    return args, defects, [pattern]

def overlay(pattern_bottom, pattern_top):
    """OR-logic applied to two arrays element wise"""
    
    length = pattern_bottom.shape[0] # We know they are square
    pattern_overlaid = np.zeros((length, length))
    for i in np.arange(length):
        for j in np.arange(length):
            if pattern_bottom[i,j] == 1 or pattern_top[i,j] == 1:
               pattern_overlaid[i,j] = 1
               
    return pattern_overlaid

def visualise(args, defects, patterns):
    """Matplotlib magic: Generates PNG files for defects, bottom layer, 
    top layer and overlaid. 1's (dots) are black, 0's (spaces) are white."""
    
    patterns.insert(0,defects)
    names = ["defects", "bottom", "top", "overlaid"]
    for index in np.arange(len(patterns)):
        image = plt.imshow(patterns[index], interpolation='nearest')
        image.set_cmap('Greys')
        plt.axis('off')
        plt.savefig(args["identity"]+"-"+names[index]+".png", 
                    bbox_inches='tight')
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
    
    distance=[]
    for distance_squared in np.arange(len(rho)):
        distance.append(np.sqrt(distance_squared)) 
    
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
    
    print ("Creating pattern...")
    args, defects, patterns = create(args)
    
    if args["width"] == args["height"] and args["no_overlay"] is False: 
        patterns.append(np.rot90(patterns[0]))
        print ("Overlaying patterns...")
        patterns.append(overlay(patterns[0], patterns[1]))
    
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
    parser.add_argument("-d", "--defects", type=int, nargs='+', default=[],
                        help="position of defects "
                        "(in 1D, remember count from 0)")
    parser.add_argument("-a", "--args", action='help', help="print args")
    parser.add_argument("-no", "--no_overlay", action="store_true", 
                        help="if selected then only a single pattern is "
                             "created")
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
    args = vars(parser.parse_args())
    
    for item in args["defects"]: 
        if item > (args["width"]*args["height"])-1:
            raise ValueError("Defects specified beyond pattern size")
    if len(args["basis"]) > args["width"]*args["height"]:
            raise ValueError("Basis specified beyond pattern size")
    if args["random"]  > (args["width"]*args["height"]):
            raise ValueError("Number of defects cannot exceed pattern size")
    
    main(args) 


