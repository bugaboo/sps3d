from sys import argv

fname1=argv[1]
file1 = open(fname1, 'r')
list1 = []
#for a in map(lambda x: x[0]+1j*x[1], list(map(lambda x: list(map(float, x[:-1].split()))[:2], open(fname1, 'r').readlines()))):
for line in file1:
    tmp = line[:-1].split()[:2]
    list1.append(complex(float(tmp[0]), float(tmp[1])))
file1.close()
fname2 = argv[2]
file2 = open(fname2, 'r')
found = False
for x in list1:
    for line in file2:
        tmp = line[:-1].split()[:2]
        z = complex(float(tmp[0]), float(tmp[1]))
        if abs(x-z) < 1e-10:
            print (abs(x-z))
            found = True
            break
    if not found:
        print ("Some number is not found")
        break
file2.close()
