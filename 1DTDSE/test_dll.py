from ctypes import *


##########################
### TEST DLL
##########################

DLL = "/Users/tadeasnemec/Programming/Git/CUPRAD_TDSE_Hankel/1DTDSE/test_dll.so"
DLL_Func = CDLL(DLL)

class nested(Structure):
    _fields_ = [("test", c_int)]

class str(Structure):
    _fields_ = [("nested", nested),
                ("arr", POINTER(c_double)),
                ("arr2", POINTER(c_double)),
                ("arr_size", c_int),
                ("number", c_int),
                ("float_ptr", POINTER(c_double))]
    
### Array size
arr_size = 10
### Init double array of size 10
test_array = c_double * arr_size
### Declare array + unwrap list with the '*' operator
test_array = test_array(*[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
### Equivalent way: test_array = test_array(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
### Declare pointer to double
float_ptr = pointer(c_double(3.14))
### Number
number = c_int(1)
### nested
n = nested()
n.test = c_int(160000)

### Init structure
s = str(arr = test_array, arr_size = arr_size, number = number, float_ptr = float_ptr)
s.nested.test = 160000
result = DLL_Func.test_struct(byref(s))