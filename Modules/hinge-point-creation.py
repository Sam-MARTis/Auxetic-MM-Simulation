from minDistanceLineSegments import Vector, nearest_points_on_line_segments
import ezdxf as dxf
import matplotlib.pyplot as plt
def load_dxf_file(file_path: str):
    doc = dxf.readfile(file_path)
    msp = doc.modelspace()
    return msp, doc

# hinge_points = []
def get_bounds_and_hinge_points(filename:str, min_distance:float, max_distance:float, center:bool, size:float):
    msp, doc = load_dxf_file(filename)
    lines_vec:list[Vector] = []
    lines = []
    hinge_points:list[Vector] = []

    for e in msp.query("INSERT"):
        block = e.block()
        m = e.matrix44()
        for line in block.query("LINE"):
            start = line.dxf.start
            end = line.dxf.end
            start = m.transform(start)
            end = m.transform(end)
            lines.append((start, end))
            lines_vec.append((Vector(start.x, start.y), Vector(end.x, end.y)))
    for e in msp.query("LINE"):
        start = e.dxf.start
        end = e.dxf.end
        lines.append((start, end))
        lines_vec.append((Vector(start.x, start.y), Vector(end.x, end.y)))
    
    for i in range(len(lines_vec)):
        for j in range(i, len(lines_vec)):
            A, B, dist = nearest_points_on_line_segments(lines_vec[i][0], lines_vec[i][1], lines_vec[j][0], lines_vec[j][1])
            if min_distance < dist < max_distance:
                hinge_point = (A + B) / 2
                hinge_points.append(hinge_point)
    minx = float("inf")
    maxx = -float("inf")
    miny = float("inf")
    maxy = -float("inf")

    for i in lines_vec:
        if i[0].x < minx: minx = i[0].x
        if i[0].x > maxx: maxx = i[0].x
        if i[0].y < miny: miny = i[0].y
        if i[0].y > maxy: maxy = i[0].y

    bounds = (minx, maxx, miny, maxy)
    if center:
        cx = (minx + maxx) / 2
        cy = (miny + maxy) / 2
        bounds = (-cx, cx, -cy, cy)
        for hp in hinge_points:
            hp.x -= cx
            hp.y -= cy
        for i in lines_vec:
            i[0].x -= cx
            i[0].y -= cy
            i[1].x -= cx
            i[1].y -= cy
    if size:
        scaling = size / (maxx - minx)
        for hp in hinge_points:
            hp.x *= scaling
            hp.y *= scaling
        for i in lines_vec:
            i[0].x *= scaling
            i[0].y *= scaling
            i[1].x *= scaling
            i[1].y *= scaling
        bounds = (bounds[0] * scaling, bounds[1] * scaling, bounds[2] * scaling, bounds[3] * scaling)
    return lines_vec, hinge_points, bounds

if __name__ == "__main__":
    lines_vec, hinge_points, bounds = get_bounds_and_hinge_points("mm_pattern_01.dxf", 1e-6, 5e-2, True, None)
    for (start, end) in lines_vec:
        plt.plot([start.x, end.x], [start.y, end.y], 'b-')
    print(f"Bounds: {bounds}")
    # for i in range(len(lines_vec)):
    #     for j in range(i, len(lines_vec)):
    #         A, B, dist = nearest_points_on_line_segments(lines_vec[i][0], lines_vec[i][1], lines_vec[j][0], lines_vec[j][1])
    #         if 1e-6<dist < 5e-2:
    #             hinge_point = (A + B) / 2
    #             hinge_points.append(hinge_point)
    #             plt.plot([hinge_point.x], [hinge_point.y], 'ro')
    plt.plot([hp.x for hp in hinge_points], [hp.y for hp in hinge_points], 'ro')
    plt.show()

