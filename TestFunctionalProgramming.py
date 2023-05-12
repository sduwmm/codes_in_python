from functools import reduce

items = [1, 2, 3, 4, 5]
for item in items:
    item + 1


def inc(x):
    return x+1
print(list(map(inc,items)))
