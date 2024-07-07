def qselect(k, arr):
    return qselect_indexed(k-1, arr)

def qselect_indexed(k, arr):
    left = [x for x in arr if x < arr[0]]

    if len(left) > k:
      return qselect_indexed(k, left)
    elif len(left) == k:
      return arr[0]

    right = [x for x in arr[1:] if x >= arr[0]]
    return qselect_indexed(k-len(left)-1, right)