N = 8

function Hamiltonian()
	H = zeros(N, N)
	H[1, 1] = 1.0
	H[2, 2] = 1.0
	H[3, 3] = 1.0
	H[4, 4] = 2.0
	H[5, 5] = 2.0
	H[6, 6] = 2.0
	H[7, 7] = 0.0
	H[8, 8] = 0.0
	return H
end

function Rate(H, B)
	w = zeros(N, N)
	for i = 1:N, j = 1:N
		if(1 == abs(i - j))
			w[i, j] = min(1, exp(-(H[i, i] - H[j, j]) * B))
		end
	end
	return w
end

function Escape(w)
	r = zeros(N, N)
	for i = 1:N, j = 1:N
		r[j, j] += w[i, j]
	end
	return r
end

function ContinousTime(w, r, s)
	exp(-s) * w - r
end

function DiscreteTime(w, r, x)
	w * inv(x * eye(N) + r)
end
