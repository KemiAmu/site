def convergent_fraction(x: float, y: float, step: int = 10) -> list:
	"""
	Computes the continued fraction convergents of x/y (or y/x if x>y) up to 'step' terms.
	Uses the recurrence relation p_n = a_n*p_{n-1} + p_{n-2} (and similarly for q_n),
	where a_n are the continued fraction coefficients, to generate successive approximations.
	"""
	
	r = 1 if x == y else x/y if x < y else y/x

	p0, q0 = int(r), 1
	p1, q1 =      1, 0

	result = []
	while r != int(r) and len(result) < step:
		r = 1 / (r - int(r))
		a = int(r)
		
		p0, p1 = a * p0 + p1, p0
		q0, q1 = a * q0 + q1, q0

		result.append((p0, q0))

	return result

if __name__ == "__main__":
    print(convergent_fraction(1.414213562373095, 1)) # sqrt(2)
    print(convergent_fraction(3.141592653589793, 1)) # pi
    print(convergent_fraction(1.618033988749895, 1)) # golden ratio
