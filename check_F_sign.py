import math
expr_path = r"C:\Users\98kim\Desktop\Acta-Astronautica\Funcs_ISS_expr\Funcs_ISS10_expr.txt"
with open(expr_path, 'r', encoding='utf-8') as f:
    expr = f.read().strip()

def F(x,y,z):
    return eval(expr, {"__builtins__":None, "sqrt": math.sqrt}, {"x": x, "y": y, "z": z})

def classify(x,y,z):
    val = F(x,y,z)
    h = val - 1.0
    if abs(h) < 1e-6:
        tag = 'ON (F?1)'
    elif h > 0:
        tag = 'OUTSIDE (F>1)'
    else:
        tag = 'INSIDE (F<1)'
    return val, h, tag

samples = [
    (0.0, 0.0, 0.0),
    (0.0, 50.0, 0.0),
    (0.0, -50.0, 0.0),
    (100.0, 0.0, 0.0),
    (0.0, 0.0, 100.0),
    (0.0, 100.0, 0.0),
    (0.0, -100.0, 0.0),
]
print(f"Testing F from: {expr_path}")
for p in samples:
    val,h,tag = classify(*p)
    print(f"p={p} -> F={val:.6e}, F-1={h:.6e}, {tag}")
