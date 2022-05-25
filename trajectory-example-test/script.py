import mdtraj as md
t = md.load('frame0.h5')
#print(t)
t.save_xyz('frame0.xyz')