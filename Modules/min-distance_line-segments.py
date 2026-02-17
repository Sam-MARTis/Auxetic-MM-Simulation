import matplotlib.pyplot as plt

class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def dot(self, other: 'Vector') -> float:
        return self.x * other.x + self.y * other.y
    def magnitude(self) -> float:
        return (self.x ** 2 + self.y ** 2) ** 0.5
    def __sub__(self, other: 'Vector') -> 'Vector':
        return Vector(self.x - other.x, self.y - other.y)
    def __add__(self, other: 'Vector') -> 'Vector':
        return Vector(self.x + other.x, self.y + other.y)
    def __truediv__(self, scalar: float) -> 'Vector':
        a = 1/scalar
        return Vector(self.x * a, self.y * a)
    def __mul__(self, scalar: float) -> 'Vector':
        return Vector(self.x * scalar, self.y * scalar)
    def __repr__(self):
        return f"v̄({self.x:.5g}, {self.y:.5g})"


"""
-- Source - https://stackoverflow.com/a/67102941
-- Posted by Trey Reynolds, modified by community. See post 'Timeline' for change history
-- Retrieved 2026-02-17, License - CC BY-SA 4.0

local eta = 1e-6
local function nearestPointsOnLineSegments(a0, a1, b0, b1)
    local r = b0 - a0
    local u = a1 - a0
    local v = b1 - b0

    local ru = r:Dot(u)
    local rv = r:Dot(v)
    local uu = u:Dot(u)
    local uv = u:Dot(v)
    local vv = v:Dot(v)

    local det = uu*vv - uv*uv
    local s, t

    if det < eta*uu*vv then
        s = math.clamp(ru/uu, 0, 1)
        t = 0
    else
        s = math.clamp((ru*vv - rv*uv)/det, 0, 1)
        t = math.clamp((ru*uv - rv*uu)/det, 0, 1)
    end

    local S = math.clamp((t*uv + ru)/uu, 0, 1)
    local T = math.clamp((s*uv - rv)/vv, 0, 1)

    local A = a0 + S*u
    local B = b0 + T*v
    return A, B, (B - A):Length()
end


"""

def nearest_points_on_line_segments(a0: Vector, a1: Vector, b0: Vector, b1: Vector):
    r = b0 - a0
    u = a1 - a0
    v = b1 - b0

    ru = r.dot(u)
    rv = r.dot(v)
    uu = u.dot(u)
    uv = u.dot(v)
    vv = v.dot(v)

    det = uu * vv - uv * uv
    s, t = 0, 0

    if det < 1e-6 * uu * vv:
        s = max(0, min(1, ru / uu))
        t = 0
    else:
        s = max(0, min(1, (ru * vv - rv * uv) / det))
        t = max(0, min(1, (ru * uv - rv * uu) / det))

    S = max(0, min(1, (t * uv + ru) / uu))
    T = max(0, min(1, (s * uv - rv) / vv))

    A:Vector = a0 + u * S
    B:Vector = b0 + v * T
    return A, B, (B - A).magnitude()


if __name__ == "__main__":
    a0 = Vector(1, 2)
    a1 = Vector(3, 4)
    b0 = Vector(1, 3)
    b1 = Vector(3, 4.5)
    A, B, dist = nearest_points_on_line_segments(a0, a1, b0, b1)
    print(f"Nearest points are {A} and {B} with distance {dist}")
    # Plot the line and point


    plt.plot([a0.x, a1.x], [a0.y, a1.y], 'b-', label='Line Segment A')
    plt.plot([b0.x, b1.x], [b0.y, b1.y], 'r-', label='Line Segment B')
    plt.plot([A.x, B.x], [A.y, B.y], 'g--', label='Nearest Points')
    plt.scatter([A.x, B.x], [A.y, B.y], color='black')
    plt.legend()
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Nearest Points on Line Segments')
    plt.grid()
    plt.axis('equal')
    plt.show()