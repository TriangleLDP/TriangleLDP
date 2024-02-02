#!/usr/bin/env python3
import csv
import sys
import numpy as np
from scipy.sparse import lil_matrix

################################# Parameters ##################################
if len(sys.argv) < 4:
    print("Usage:",sys.argv[0],"[GPlusFile (in)] [EdgeFile (out)] [DegFile (out)]")
    sys.exit(0)

# GPlus File (input)
GPlusFile = sys.argv[1]
# Or Edge File (output)
OrEdgeFile = sys.argv[2]
# Degree File (output)
DegFile = sys.argv[3]

UserNum = 107614

#################################### Main #####################################

edges_lil = lil_matrix((UserNum, UserNum))
deg = np.zeros(UserNum)

# Read edges from the GPlus file --> edges_lil
f = open(GPlusFile, "r")
nodeid2index = {}
node_num = 0
print("Reading the GPlus file")
for i, line in enumerate(f):
    if i % 1000000 == 0:
        print(i)
    lst = line.rstrip("\n").split(" ")
    nodeid1 = int(lst[0])
    nodeid2 = int(lst[1])
    if nodeid1 not in nodeid2index:
        nodeid2index[nodeid1] = node_num
        node_num += 1
    if nodeid2 not in nodeid2index:
        nodeid2index[nodeid2] = node_num
        node_num += 1
    index1 = nodeid2index[nodeid1]
    index2 = nodeid2index[nodeid2]
    edges_lil[index1, index2] = 1
f.close()

a1, a2 = edges_lil.nonzero()

# Output edge information
print("Outputting edge information.")
f = open(OrEdgeFile, "w")
print("#nodes", file=f)
print(UserNum, file=f)
print("node,node", file=f)
writer = csv.writer(f, lineterminator="\n")
for i in range(len(a1)):
    # user_ids --> user_id1, user_id2
    user_id1 = a1[i]
    user_id2 = a2[i]
    if edges_lil[user_id2, user_id1] == 1:
        if user_id1 < user_id2:
            lst = [user_id1, user_id2]
            writer.writerow(lst)
            deg[user_id1] += 1
            deg[user_id2] += 1
    else:
        if user_id1 < user_id2:
            lst = [user_id1, user_id2]
        else:
            lst = [user_id2, user_id1]
        writer.writerow(lst)
        deg[user_id1] += 1
        deg[user_id2] += 1
f.close()

# Output degree information
print("Outputting degree information.")
f = open(DegFile, "w")
print("node,deg", file=f)
writer = csv.writer(f, lineterminator="\n")
for user1 in range(UserNum):
    # user index and her degree --> lst
    lst = [user1, int(deg[user1])]
    writer.writerow(lst)
f.close()
