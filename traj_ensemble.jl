N = 8
B = 1 / float(ARGS[1])
tobs = float(ARGS[2])
kobs = float(ARGS[3])

I = eye(N)

function ContinousMatrix(w, r, s)
	exp(-s) * w - r
end

function DiscreteMatrix(w, r, x)
	q = zeros(N, N)

	for i = 1:N
		q[i, i] = 1 / (x + r[i, i])
	end

	return w * q
end

H = zeros(N, N)
w = zeros(N, N)
r = zeros(N, N)

H[1, 1] = 1.0
H[2, 2] = 1.0
H[3, 3] = 1.0
H[4, 4] = 2.0
H[5, 5] = 2.0
H[6, 6] = 2.0
H[7, 7] = 0.0
H[8, 8] = 0.0

# Hamiltonian
for i = 1:N, j = 1:N
	w[i, j] = min(1, exp(-(H[i, i] - H[j, j]) * B))
end

# Random walk
for i = 1:N, j = 1:N
	if(1 != abs(i - j))
		w[i, j] = 0
	end
end

# Escaping ratio
for i = 1:N, j = 1:N
	r[j, j] += w[i, j]
end

W = zeros(N, N)
dW= zeros(N, N)

T = zeros(N, N)
dT= zeros(N, N)

p0 = zeros(N)
pj = zeros(N)

p0[1:N] = 1 / N
pj[1:N] = 1

@printf("# 1 bias, 2 Zs (inf), 3 k (inf), 4 Zx (inf), 5 t (inf), 6 k (finite), 7 Es (finite), 8 Cs (finite), 9 t (finite), 10 Ex (finite), 11 Cx (finite)\n")

for ib = -500:500
	s = ib / 5000
	x = ib / 5000

	ds = 1 / 1000000
	dx = 1 / 1000000

	W = ContinousMatrix(w, r, s)
	dW= ContinousMatrix(w, r, s + ds)

	T = DiscreteMatrix(w, r, x)
	dT= DiscreteMatrix(w, r, x + dx)

	###

	zs = dot(pj, expm(tobs * W) * p0)
	zd = dot(pj, expm(tobs *dW) * p0)

	qx = dot(pj, T ^ kobs * p0)
	qd = dot(pj,dT ^ kobs * p0)

	Et = dot(pj, H * expm(tobs * W) * p0)
	Et/= zs

 	Ct = dot(pj, H * H * expm(tobs * W) * p0)
	Ct/= zs
	Ct-= Et * Et

	Ek = dot(pj, H * T ^ kobs * p0)
	Ek/= qx

	Ck = dot(pj, H * H * T ^ kobs * p0)
	Ck/= qx
	Ck-= Ek * Ek

	zs = log(zs) / tobs
	zd = log(zd) / tobs
	qx = log(qx) / kobs
	qd = log(qd) / kobs

	###

	U = eigmax(W)
	dU= eigmax(dW)

	V = eigmax(T)
	dV= eigmax(dT)

	K = (U - dU) / ds
	t = (V - dV) / dx

	@printf("%.8f, %.8f, %.8f, %.8f, %.8f, ", s, U, K, V, t)

	K = -(zd - zs) / ds
	t = -(qd - qx) / dx

	@printf("%.8f, %.8f, %.8f, %.8f, %.8f, %.8f\n", K, Et, Ct, t, Ek, Ck)
end

