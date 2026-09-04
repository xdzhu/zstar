"""Exact symbolic checks of shared-response definitions (research dependency)."""

import sympy as s

u = s.Matrix(s.symbols('u0:6'))
field = s.Matrix(s.symbols('E0:3'))
charge = s.symbols('e', positive=True)
z = s.Matrix(3, 6, s.symbols('z0:18'))
h = s.Matrix(6, 6, s.symbols('h0:36'))
phi = (h + h.T) / 2
f0 = s.Matrix(s.symbols('f0:6'))
energy = (u.T * phi * u)[0] / 2 - (f0.T * u)[0] - charge * (field.T * z * u)[0]
force = -s.Matrix([s.diff(energy, q) for q in u])
dipole = -s.Matrix([s.diff(energy, e) for e in field])
assert s.simplify(force.jacobian(u) + phi) == s.zeros(6)
assert s.simplify(dipole.jacobian(u) / charge - z) == s.zeros(3, 6)
assert s.simplify(force.jacobian(field) / charge - z.T) == s.zeros(6, 3)

# Neither Cartesian seeds nor a symmetric Born tensor are required.
x = s.Matrix([[1, 2, 3], [2, -1, 1], [3, 1, -1]])
b = s.Matrix(9, 3, s.symbols('b0:27'))
assert x.det() != 0
assert s.simplify((b * x) * x.inv() - b) == s.zeros(9, 3)

t, linear, quadratic, cubic = s.symbols('delta L Q C')
response = linear*t + quadratic*t**2 + cubic*t**3
central = s.expand((response - response.subs(t, -t)) / (2*t))
assert central == linear + cubic*t**2
print('PASS: force Hessian, Born/Maxwell derivatives, mixed-direction recovery, central truncation order')
