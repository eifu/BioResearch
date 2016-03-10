import histone
from model4 import nextGen, trackingHist, bitvec

histoneList1 = histone.createRandomHistoneList(50, 1, 81, 40)


histoneList2 = nextGen(histoneList1,0,1)

for i in range(81):
    histoneList2[i].display()

print(bitvec(histoneList2))

# tr = trackingHist(histoneList1,10,0,1)
# 
# for i in range(10):
#     print(tr[i][0])

#     histoneList2 = histone.createRandomHistoneList(50, 0, NUM_OF_HISTONE, BEFORE_PROMOTER)
#     x = np.array([histoneList1])
#     print(len(x))
#     print(x.shape)
#     print(x.dtype)
#     print(type(x))

#     print(histoneList)
#     bitvec(histoneList)
#     plt.bar(histoneList,[1 for _ in range(NUM_OF_HISTONE)])
#     plt.show()
#     v = bitvec(histoneList1)
#     print(v)
#     print(v.shape)
#     print(v.dtype)
#     print(type(v))