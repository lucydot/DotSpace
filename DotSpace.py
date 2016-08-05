#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed

def create(args):
    basis_length = len(args["basis"])
    array = np.zeros(args["width"]*args["height"], dtype=int)
    defects = np.zeros(args["width"]*args["height"], dtype=int)
    
    for position in args["defects"]:
        defects[position] = 1
    
    zeros = np.where(defects == 0)[0]
    random_numbers = np.random.choice(zeros, args["random"], replace=False)
    args["random_numbers"] = random_numbers

    for position in random_numbers:
        defects[position] = 1
    
    counter = 0 
    for position in np.arange(len(array)):
        if defects[position] == 1:
            array[position] = 1
            counter = 1
        else:
            if args["basis"][counter] == 1:
                array[position] = 1
            counter = ((counter + 1) % basis_length)
    
    array.shape = (args["height"], args["width"])
    defects.shape = (args["height"], args["width"])
    return args, defects, [array]

def overlay(array_one, array_two):
    dim = array_one.shape[0]
    array = np.zeros((dim, dim))
    for i in np.arange(dim):
        for j in np.arange(dim):
            if array_one[i,j] == 1 or array_two[i,j] == 1:
               array[i,j] = 1
    return array

def visualise(args, defects, arrays):
    arrays.insert(0,defects)
    array_names = ["defects", "bottom", "top", "both"]
    for index in np.arange(len(arrays)):
        image = plt.imshow(arrays[index], interpolation='nearest')
        image.set_cmap('Greys')
        plt.axis('off')
        plt.savefig(args["identity"]+"-"+array_names[index]+".png", bbox_inches='tight')
        plt.close()

def correlate(args, array):
    # See Ziman, models of disorder. 
    # Thanks to Andrew Goodwin for discussion
    array[array == 0] = -1

    if args["cutoff"] is None:
        X = array.shape[0]
        Y = array.shape[1]
    else:
        X = args["cutoff"]+1
        Y = args["cutoff"]+1

    dX = np.arange(-X+1, X)
    dY = np.arange(-Y+1, Y)

    total = np.zeros((X-1)*(X-1)+(Y-1)*(Y-1)+1)
    count = np.zeros((X-1)*(X-1)+(Y-1)*(Y-1)+1)
    store=[]
    for x in np.arange(X):
        for y in np.arange(Y):
            for dx in dX:
                if x+dx < 0:
                    pass
                else:
                    for dy in dY:
                        if y+dy < 0:
                            pass
                        else:    
                            if [x+dx,y+dy,x,y] in store:
                                pass
                            else:
                                store.append([x,y,dx,dy])
                                try:
                                    product = array[x,y]*array[x+dx,y+dy]
                                    total[dx*dx+dy*dy]=total[dx*dx+dy*dy]+product
                                    count[dx*dx+dy*dy]=count[dx*dx+dy*dy]+1
                                except IndexError:
                                    pass
    rho = np.divide(total,count)
    r=[]
    for index in np.arange(len(rho)):
        r.append(np.sqrt(index)) 
    fig, ax = plt.subplots()
    ax.stem(r, rho, markerfmt=' ')
    plt.savefig(args["identity"]+"-corr.png", 
                bbox_inches='tight')
    plt.close()
    with open(args["identity"]+"-corr.txt","w") as f:
        f.write(str(np.vstack([np.multiply(r,r),r,rho, count]).T))
        
    
    correlation = 3 
    data = { 'correlation' : correlation }
    print ("Correlation length is: "+ str(correlation))
    return data

def save(args, data):
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
    
    print ("Creating arrays...")
    args, defects, arrays = create(args)
    
    if args["width"] == args["height"] and args["no_overlay"] is False: 
        arrays.append(np.rot90(arrays[0]))
        print ("Overlaying arrays...")
        arrays.append(overlay(arrays[0], arrays[1]))
    
    print ("Plotting visual representation of array...")
    visualise(args, defects, arrays)
   
    if args["no_correlation"] is not True:
        print ("Calculating correlation lengths...")
        data = correlate(args,arrays[-1])
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
                        help="width of array")
    parser.add_argument("-h", "--height", required=True, type=int, 
                        help="height of arary")
    parser.add_argument("-b", "--basis", required=True, type=int, nargs='+',
                        help="repeating basis pattern: 1=dot, 0=empty")
    parser.add_argument("-d", "--defects", type=int, nargs='+', default=[],
                        help="position of defects (in 1D, count from 0)")
    parser.add_argument("-a", "--args", action='help', help="print args")
    parser.add_argument("-no", "--no_overlay", action="store_true", 
                        help="if selected then only a single array is created")
    parser.add_argument("-nc", "--no_correlation", action="store_true",
                        help="if selected then correlation length will not be"
                        " calculated")
    parser.add_argument("-r", "--random", type=int, default=0,
                        help="number of defects to place randomly throughout"
                        "the array. Default is 0.")
    parser.add_argument("-x", "--cutoff", type=int, default=None,
                        help="largest horizontal/vertical distance for "
                        "correlation calculation")
    args = vars(parser.parse_args())
    for item in args["defects"]: 
        if item > (args["width"]*args["height"])-1:
            raise ValueError("Defects specified beyond array size")
    if len(args["basis"]) > args["width"]*args["height"]:
            raise ValueError("Basis beyond array size")
    if args["random"]  > (args["width"]*args["height"]):
            raise ValueError("Number of defects cannot exceed array size")
    main(args) 


