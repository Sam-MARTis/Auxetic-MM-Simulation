from minDistanceLineSegments import Vector, nearest_points_on_line_segments
import ezdxf as dxf
import matplotlib.pyplot as plt
def load_dxf_file(file_path: str):
    doc = dxf.readfile(file_path)
    msp = doc.modelspace()
    return msp, doc


lines = []
lines_vec:list[Vector] = []

hinge_points = []
if __name__ == "__main__":
    msp, doc = load_dxf_file("mm_pattern_01.dxf")
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

        # print(start, end)
    # Plot all lines
    for (start, end) in lines_vec:
        plt.plot([start.x, end.x], [start.y, end.y], 'b-')
        
    for i in range(len(lines_vec)):
        for j in range(i, len(lines_vec)):
            A, B, dist = nearest_points_on_line_segments(lines_vec[i][0], lines_vec[i][1], lines_vec[j][0], lines_vec[j][1])
            if 1e-6<dist < 5e-2:
                hinge_point = (A + B) / 2
                hinge_points.append(hinge_point)
                plt.plot([hinge_point.x], [hinge_point.y], 'ro')
    plt.show()

