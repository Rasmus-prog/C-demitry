# Compare the number of integrand evaluations with the python/numpy's integration routines. 
import numpy as np
from scipy import integrate


# Counter for function evaluations
class FunctionCounter:
    def __init__(self, func):
        self.func = func
        self.calls = 0
    
    def __call__(self, x):
        self.calls += 1
        return self.func(x)
    
    def reset(self):
        self.calls = 0


#  ∫01 dx 1/√(x) = 2 ,
def f2(x):
    return 1.0/np.sqrt(x)

#  ∫01 dx ln(x)/√(x) = -4 .
def f4(x):
    return np.log(x)/np.sqrt(x)

# Test your implementation on some (converging)
# ∫0∞ dx e^(-x) = 1
def f5(x):
    return np.exp(-x)


# Wrap functions with counters
f2 = FunctionCounter(f2)
f4 = FunctionCounter(f4)
f5 = FunctionCounter(f5)

result2, err2 = integrate.quad(f2, 0.0, 1.0)
result4, err4 = integrate.quad(f4, 0.0, 1.0)
result5, err5 = integrate.quad(f5, 0.0, np.inf)

print("f2 ∫01 dx 1/√(x) = 2")
print("result = ", result2)
print("result is within 1e-3 of 2: ", np.abs(result2 - 2.0) < 1e-3)
print("number of iterations (function evaluations): ", f2.calls)

print("\nf4 ∫01 dx ln(x)/√(x) = -4")
print("result = ", result4)
print("result is within 1e-3 of -4: ", np.abs(result4 + 4.0) < 1e-3)
print("number of iterations (function evaluations): ", f4.calls)

print("\nf5 ∫0∞ dx e^(-x) = 1")
print("result = ", result5)
print("result is within 1e-3 of 1: ", np.abs(result5 - 1.0) < 1e-3)
print("number of iterations (function evaluations): ", f5.calls)