def plot_line(
	p0: tuple[float, float],
	p1: tuple[float, float]
) -> list[tuple[int, int]]:
	"""
	Implements a modified Bresenham's algorithm to generate integer pixel coordinates
	approximating a line between two floating-point endpoints. The algorithm operates
	by iteratively selecting adjacent pixels that minimize the perpendicular distance
	to the ideal mathematical line. It handles fractional starting positions through
	error term initialization that accounts for subpixel offsets.

	The core mathematical mechanism tracks a floating-point error term representing
	the signed vertical (or horizontal) distance from the true line. At each iteration,
	the algorithm moves either horizontally or vertically based on which direction
	minimizes accumulated error. The error term is updated using the line's slope
	properties (dx and dy) to maintain accuracy without requiring floating-point
	division operations.
	"""

	x = round(p0[0])
	y = round(p0[1])
	dx = abs(p0[0] - p1[0])
	dy = abs(p0[1] - p1[1])
	sx = 1 if p0[0] < p1[0] else -1
	sy = 1 if p0[1] < p1[1] else -1
	
	err = .5 * (dx-dy) + dy * sx * (p0[0]-x) - dx * sy * (p0[1]-y)

	points = []
	for _ in range(abs(x-round(p1[0])) + abs(y-round(p1[1]))):
		if err < 0: y += sy; err += dx
		else:       x += sx; err -= dy
		points.append((x, y))

	return points

def adaptive_bezier_path(
	p0: tuple[float, float],
	p1: tuple[float, float],
	p2: tuple[float, float],
	p3: tuple[float, float]
) -> list[tuple[float, float]]:
	"""
	Generates a smooth cubic Bezier curve through adaptive subdivision. The algorithm
	recursively splits the curve segments until the midpoint deviation from the
	chord approximation falls below a threshold (0.001), ensuring accurate rendering
	with minimal points. Returns a list of floating-point coordinates.
	"""

	def _bezier_point(t):
		u = 1.0 - t
		return (
			u**3*p0[0] + 3*u**2*t*p1[0] + 3*u*t**2*p2[0] + t**3*p3[0],
			u**3*p0[1] + 3*u**2*t*p1[1] + 3*u*t**2*p2[1] + t**3*p3[1]
		)

	curve = [p0]
	stack = [(0.5, 1.0), (0.0, 0.5)]

	while stack:
		t0, t1 = stack.pop()
		tm = (t0 + t1) / 2
		p_start = curve[-1]
		p_end = _bezier_point(t1)
		p_mid = _bezier_point(tm)

		mid = ((p_start[0] + p_end[0]) / 2, (p_start[1] + p_end[1]) / 2)
		if (p_mid[0] - mid[0])**2 + (p_mid[1] - mid[1])**2 < .001:
			curve.append(p_end)
		else:
			stack.append((tm, t1))
			stack.append((t0, tm))
	
	return curve

def cumulative_bezier(
	p0: tuple[float, float],
	p1: tuple[float, float],
	p2: tuple[float, float],
	p3: tuple[float, float]
) -> list[int]:
	"""
	Computes cumulative movement vectors along a Bezier curve path. The function
	first generates the curve's pixel coordinates, then calculates directional
	movement vectors by aggregating consecutive pixel steps with matching slopes.
	Returns a list of integers representing alternating x and y components of
	the cumulative movement vectors.
	"""
	
	bezier = adaptive_bezier_path(p0, p1, p2, p3)
	
	points = [(round(p0[0]), round(p0[1]))] + [
		point
		for start, end in zip(bezier[:-1], bezier[1:])
		for point in plot_line((start[0], start[1]), (end[0], end[1]))
	]

	delta = [[0, 0]]
	for (x0, y0), (x1, y1) in zip(points, points[1:]):
		dx, dy = x1 - x0, y1 - y0
		
		if dx * delta[-1][1] != dy * delta[-1][0]:
			delta.append([0, 0])
			
		delta[-1][0] += dx
		delta[-1][1] += dy

	return [d[0] or d[1] for d in delta]



if __name__ == "__main__":

	import matplotlib.pyplot as plt

	val = [
		(0, 100),
		(199, 0),
		(-99, 0),
		(100, 100)
	]
	
	# val = [
	# 	(0, 0),
	# 	(200, 50),
	# 	(-100, 50),
	# 	(100, 100)
	# ]

	curve = adaptive_bezier_path(*val)
	cumulative = cumulative_bezier(*val)

	# Plot the adaptive bezier curve in blue
	plt.plot([p[0] for p in curve], [p[1] for p in curve], 'b-')

	# Plot the cumulative vectors in red
	x, y = curve[0][0], curve[0][1]
	for dx, dy in zip(cumulative[::2], cumulative[1::2]):
		plt.plot([x, x + dx], [y, y + dy], 'r-')
		x += dx
		y += dy

	plt.gca().set_aspect('equal', adjustable='box')
	plt.grid(True)
	plt.show()
