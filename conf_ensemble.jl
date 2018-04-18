module Matrix
	include("matrix.jl")
end

# constant variables
B = 1 / float(ARGS[1])
Tmax = float(ARGS[2])

dt = 1 / 1000
ds = 1 / 1000000
dx = 1 / 1000000

# matrix
H = Matrix.Hamiltonian()
w = Matrix.Rate(H, B)
r = Matrix.Escape(w)

W = Matrix.ContinousTime(w, r, 0)
Wp= Matrix.ContinousTime(w, r, ds)
Wm= Matrix.ContinousTime(w, r,-ds)
T = Matrix.DiscreteTime(w, r, 0)

Tp= Matrix.DiscreteTime(w, r, dx)
Tm= Matrix.DiscreteTime(w, r,-dx)

# vector
p0 = zeros(Matrix.N)
pj = zeros(Matrix.N)
p0[3] = 1
pj[1:Matrix.N] = 1

# iteration
tobs = 0.01

@printf("0 t_obs,1 t_obs/<K>,2 P_bound,3 E(t),4 C(t)\n")

while tobs < Tmax
	pt = expm(W * tobs) * p0

	Et = dot(pj, H * expm(Wp * tobs) * p0)
	Ct = dot(pj, H * H * expm(Wp * tobs) * p0) - Et * Et

	zp = log(dot(pj, expm(Wp * tobs) * p0)) / tobs
	zm = log(dot(pj, expm(Wm * tobs) * p0)) / tobs
	k = - (zp - zm) / (2 * ds)

	@printf("%.16f,%.16f,", tobs, 1 / k)
	@printf("%.16f,", pt[Matrix.N])
	@printf("%.16f,%.16f\n", Et, Ct)

	tobs *= 1.1
end
