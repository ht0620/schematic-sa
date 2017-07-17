module Matrix
	include("matrix.jl")
end

# constant variables
tobs = float(ARGS[1])
ds = 1 / 1000000

# matrix
H = Matrix.Hamiltonian()

# vector
p0 = zeros(Matrix.N)
pj = zeros(Matrix.N)
p0[1:Matrix.N] = 1 / Matrix.N
pj[1:Matrix.N] = 1

# iteration
tobs = 1.0

@printf("# 1 temp, 2 t_obs / <K>, 3 P_bound, 4 E(t), 5 C(t)\n")

temp = 0.1

while temp < 1
	B = 1 / temp

	w = Matrix.Rate(H, B)
	r = Matrix.Escape(w)
	
	W = Matrix.ContinousTime(w, r, 0)
	Wp= Matrix.ContinousTime(w, r, ds)
	Wm= Matrix.ContinousTime(w, r,-ds)

	pt = zeros(Matrix.N)
	pt = expm(W * tobs) * p0

	Et = dot(pj, H * pt)
	Ct = dot(pj, H * H * pt) - Et * Et

	zp = log(dot(pj, expm(Wp * tobs) * p0)) / tobs
	zm = log(dot(pj, expm(Wm * tobs) * p0)) / tobs
	k = - (zp - zm) / (2 * ds)

	println(p0)
	println(pt)

#	@printf("%.16f %.16f ", temp, 1 / k)
#	@printf("%.16f ", pt[Matrix.N - 1] + pt[Matrix.N])
#	@printf("%.16f %.16f\n", Et, Ct)

	temp += 0.1
end
