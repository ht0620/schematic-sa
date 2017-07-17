N = 8
T = float(ARGS[1])

dt = 1 / 1000
ds = 1 / 1000000

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

H[1, 1]	= 1.0
H[2, 2]	= 1.0
H[3, 3]	= 1.0
H[4, 4] = 2.0
H[5, 5] = 2.0
H[6, 6] = 2.0
H[7, 7] = 0.0
H[8, 8] = 0.0

for i = 1:N, j = 1:N
	w[i, j] = min(1, exp(-(H[i, i] - H[j, j]) / T))
end

for i = 1:N, j = 1:N
	if(1 != abs(i - j))
		w[i, j] = 0
	end
end

for i = 1:N, j = 1:N
	r[j, j] += w[i, j]
end

W = zeros(N, N)
W = ContinousMatrix(w, r, 0)

Wp = zeros(N, N)
Wm = zeros(N, N)

Wp = ContinousMatrix(w, r, ds)
Wm = ContinousMatrix(w, r,-ds)

T = zeros(N, N)
T = DiscreteMatrix(w, r, 0)

p0 = zeros(N)
pj = zeros(N)

pt = zeros(N)
pk = zeros(N)

p0[4:6] = 1 / 3
pj[1:N] = 1

t = 1.0

@printf("# 1 t_obs, 2 t_obs / <K>, 3 P_bound, 4 E(t), 5 C(t)\n")

while t < 100000000
	pt = expm(W * t) * p0

	Et = dot(pj, H * expm(Wp * t) * p0)
	Ct = dot(pj, H * H * expm(Wp * t) * p0)
	Ct-= Et * Et

	zp = dot(pj, expm(Wp * t) * p0)
	zm = dot(pj, expm(Wm * t) * p0)
	zp = log(zp) / t
	zm = log(zm) / t
	k = - (zp - zm) / (2 * ds)

	@printf("%.16f %.16f %.16f %.16f %.16f\n", t, 1 / k, pt[N - 1] + pt[N], Et, Ct)

	t *= 1.1
end
