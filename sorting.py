import random
array = [1,93,42,33,39,25252,50,99999,7,5,4,3,1,2,444,555,666,5,4,3,2,1]

number = 0
IsSorted = False
length = len(array)
while(not IsSorted):
    IsSorted = True
    i = random.randint(0,length-1)
    j = random.randint(0,length-1)
    if(i > j):
        temp = i
        i = j
        j = temp
    if(array[i] > array[j]):
        IsSorted = False
        temp = array[i]
        array[i] = array[j]
        array[j] = temp
        number += 1
    for k in range(0,length-1):
        if(array[k] > array[k+1]):
            IsSorted = False
        
    number += 1 # counter

print("Sorted number list") 
print(array)
print("Total sorting number: %d" % number)
