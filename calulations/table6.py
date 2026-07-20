import math
import pandas as pd


# Minimum polynomial degrees (6 decimal places)

degrees = {
    "sin": {
        0.2: 5,
        0.4: 6,
        0.6: 7,
        0.8: 8
    },
    "exp": {
        0.2: 5,
        0.4: 6,
        0.6: 8,
        0.8: 9
    },
    "ln": {
        0.2: 7,
        0.4: 13,
        0.6: 22,
        0.8: 47
    }
}


# Taylor polynomial functions

def taylor_sin(x, n):
    total = 0
    k = 0
    while 2*k + 1 <= n:
        total += ((-1)**k) * x**(2*k+1) / math.factorial(2*k+1)
        k += 1
    return total


def taylor_exp(x, n):
    total = 0
    for k in range(n+1):
        total += x**k / math.factorial(k)
    return total


def taylor_ln(x, n):
    total = 0
    for k in range(1, n+1):
        total += ((-1)**(k+1)) * x**k / k
    return total


# Exact functions

functions = {
    "sin": math.sin,
    "exp": math.exp,
    "ln": lambda x: math.log(1+x)
}

taylor = {
    "sin": taylor_sin,
    "exp": taylor_exp,
    "ln": taylor_ln
}


# Lagrange remainder bounds

def lagrange_bound(func, x, n):

    if func == "sin":
        M = 1
        return M * x**(n+1) / math.factorial(n+1)

    elif func == "exp":
        M = math.exp(x)
        return M * x**(n+1) / math.factorial(n+1)

    elif func == "ln":
        # Maximum derivative occurs at x=0
        return x**(n+1)/(n+1)



# Generating Table 3.5

rows = []

for func in ["sin","exp","ln"]:

    for x,n in degrees[func].items():

        true_value = functions[func](x)

        approx = taylor[func](x,n)

        actual_error = abs(true_value-approx)

        bound = lagrange_bound(func,x,n)

        ratio = bound/actual_error

        rows.append({
            "Function":func,
            "x":x,
            "Degree n":n,
            "True Value":true_value,
            "Taylor Approximation":approx,
            "Actual Error":actual_error,
            "Lagrange Bound":bound,
            "Conservativeness Ratio":ratio
        })

df = pd.DataFrame(rows)

pd.set_option("display.precision",12)

print(df)


# Save to CSV
df.to_csv("Table3_5.csv",index=False)