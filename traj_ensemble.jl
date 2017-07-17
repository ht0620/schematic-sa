module Matrix
	include("matrix.jl")
end

# constants
B = 1 / float(ARGS[1])
tobs = float(ARGS[2])
kobs = float(ARGS[3])

# matrix
H = Matrix.Hamiltonian()
w = Matrix.Rate(H, B)
r = Matrix.Escape(w)

p0 = zeros(Matrix.N)
pj = zeros(Matrix.N)

# vector
p0[1:Matrix.N] = 1 / Matrix.N
pj[1:Matrix.N] = 1

#iteration
@printf("# 1 bias, 2 Zs (inf), 3 k (inf), 4 Zx (inf), 5 t (inf), 6 k (finite), 7 Es (finite), 8 Cs (finite), 9 t (finite), 10 Ex (finite), 11 Cx (finite)\n")

for ib = -500:500
	s = ib / 5000
	x = ib / 5000

	ds = 1 / 1000000
	dx = 1 / 1000000

	# continous time
	W = Matrix.ContinousTime(w, r, s)
	dW= Matrix.ContinousTime(w, r, s + ds)

	# discrete time
	T = Matrix.DiscreteTime(w, r, x)
	dT= Matrix.DiscreteTime(w, r, x + dx)

	# generating function (finite)
	zs = dot(pj, expm(tobs * W) * p0)
	zd = dot(pj, expm(tobs *dW) * p0)
	qx = dot(pj, T ^ kobs * p0)
	qd = dot(pj,dT ^ kobs * p0)

	# energy, susceptibility (finite)
	Et = dot(pj, H * expm(tobs * W) * p0) / zs
 	Ct = dot(pj, H * H * expm(tobs * W) * p0) / zs - Et * Et
	Ek = dot(pj, H * T ^ kobs * p0) / qx
	Ck = dot(pj, H * H * T ^ kobs * p0) / qx - Ek * Ek

	# large deviation function (finite)
	zs = log(zs) / tobs
	zd = log(zd) / tobs
	qx = log(qx) / kobs
	qd = log(qd) / kobs

	# generating function (infinite)
	U = eigmax(W)
	dU= eigmax(dW)
	V = eigmax(T)
	dV= eigmax(dT)

	# activity, time (infinite)
	K = (U - dU) / ds
	t = (V - dV) / dx

	# output (infinite)
	@printf("%.8f, %.8f, %.8f, %.8f, %.8f, ", s, U, K, V, t)

	#activity, time (finite)
	K = -(zd - zs) / ds
	t = -(qd - qx) / dx

	# output (finite)
	@printf("%.8f, %.8f, %.8f, %.8f, %.8f, %.8f\n", K, Et, Ct, t, Ek, Ck)
end
