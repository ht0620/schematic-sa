module Matrix
	include("matrix.jl")
end

# constants
B = 1 / float(ARGS[1])
tobs = float(ARGS[2])

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
@printf("0 bias,1 K(inf),2 t(inf),3 K(finite),4 t(fin)\n")

for ib = -500:500
	s = ib / 5000
	ds = 1 / 1000000

	# continous time
	W = Matrix.ContinousTime(w, r, s)
	dW= Matrix.ContinousTime(w, r, s + ds)

	# generating function (finite)
	zs = dot(pj, expm(tobs * W) * p0)
	zd = dot(pj, expm(tobs *dW) * p0)

	# large deviation function (finite)
	zs = log(zs) / tobs
	zd = log(zd) / tobs

	# generating function (infinite)
	U = eigmax(W)
	dU= eigmax(dW)

	# activity, time (infinite)
	K = (U - dU) / ds
	t = 1 / K

	# output (infinite)
	@printf("%.16f,", s)
	@printf("%.16f,%.16f,", K, t)

	#activity, time (finite)
	K = -(zd - zs) / ds
	t = 1 / K

	# output (finite)
	@printf("%.16f,%.16f\n", K, t)
end
