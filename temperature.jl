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
p0[4:6] = 1 / 3
pj[1:Matrix.N] = 1

# iteration
@printf("# 1 temp, 2 t_obs / <K>, 3 P_bound, 4 E(t), 5 C(t)\n")

temp = 0.01

while temp < 2
	B = 1 / temp

	w = Matrix.Rate(H, B)
	r = Matrix.Escape(w)

	W = Matrix.ContinousTime(w, r, 0)
	Wp= Matrix.ContinousTime(w, r, ds)
	Wm= Matrix.ContinousTime(w, r,-ds)

	pt = expm(tobs * W) * p0
	
	Et = dot(pj, H * pt)
	Ct = dot(pj, H * H * pt) - Et * Et

	zp = log(dot(pj, expm(Wp * tobs) * p0)) / tobs
	zm = log(dot(pj, expm(Wm * tobs) * p0)) / tobs
	k = - (zp - zm) / (2 * ds)

	@printf("%.16f %.16f ", temp, 1 / k)
	@printf("%.16f ", pt[Matrix.N - 1] + pt[Matrix.N])
	@printf("%.16f %.16f\n", Et, Ct)

	temp += 0.01
end
