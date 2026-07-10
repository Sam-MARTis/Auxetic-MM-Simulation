import ezdxf
import numpy as np
from numpy import cos, sin

class Vec2D:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def normalize(self):
        mag = (self.x*self.x + self.y*self.y)**0.5
        return Vec2D(self.x / mag, self.y / mag)
    def __add__(self, other):
        return Vec2D(self.x + other.x, self.y + other.y)
    def __mul__(self, other):
        return Vec2D(self.x*other, self.y *other)
    def __sub__(self, other):
        return Vec2D(self.x - other.x, self.y - other.y)
    
    

class Unit_Cell: 
    def __init__(self, s, t, θ, a):
        self.internalCuts = []
        self.boundaryCuts = []
        
        deg_to_rad = 3.14159265/180
        
        cornerPoint1 = Vec2D(0, 0)
        cornerPoint2 = Vec2D(s*cos(60*deg_to_rad), s*sin(60*deg_to_rad))
        cornerPoint3 = Vec2D(s, 0)
        c12 = (cornerPoint2 - cornerPoint1).normalize()
        c23 = (cornerPoint3 - cornerPoint2).normalize()
        c31 = (cornerPoint1 - cornerPoint3).normalize()
        self.boundaryCuts.append([cornerPoint1, cornerPoint2])
        self.boundaryCuts.append([cornerPoint2, cornerPoint3])
        self.boundaryCuts.append([cornerPoint3, cornerPoint1])

        hp1 = cornerPoint1 + c12*t
        hp2 = cornerPoint2 + c23*t
        hp3 = cornerPoint3 + c31*t 
        
        d1 = Vec2D(1, np.tan(θ*deg_to_rad)).normalize()
        d2 = Vec2D(-1, -np.tan((60 + θ)*deg_to_rad)).normalize()
        d3 = Vec2D(-1, np.tan((60-θ)*deg_to_rad)).normalize()
        ep1 = hp1 + d1*a
        ep2 = hp2 + d2*a
        ep3 = hp3 + d3*a 

        self.internalCuts.append([hp1, ep1])
        self.internalCuts.append([hp2, ep2]) 
        self.internalCuts.append([hp3, ep3])


        
        # self.cut_size = cut_size
        # self.θ_in_radians = θ_in_radians
        # self.edge_gap = edge_gap
        # self.side_gap = side_gap
        # self.a = cut_size + edge_gap
        # self.l = (self.side_gap + self.a * cos(θ_in_radians)) / sin(θ_in_radians)

        # self.dx = self.side_gap + self.l * cos(θ_in_radians) + self.a * sin(θ_in_radians)
        # self.dy = self.dx 
        
        # self.xlb = 0
        # self.ylb = 0
        # self.xrb = self.dx
        # self.yrb = 0
        # self.xlt = 0
        # self.xrt = self.xrb
        # self.yrt = self.dy
        # self.ylt = self.yrt

        # self.lb_cut_x0 = self.xlb + self.side_gap
        # self.lb_cut_y0 = self.ylb
        # self.lb_cut_x1 = self.lb_cut_x0 + self.cut_size * cos(θ_in_radians)
        # self.lb_cut_y1 = self.lb_cut_y0 + self.cut_size * sin(θ_in_radians)

        # # self.lb_cut = self.cell.add_line((self.lb_cut_x0, self.lb_cut_y0), (self.lb_cut_x1, self.lb_cut_y1))
        
        
        # self.rb_cut_x0 = self.xrb
        # self.rb_cut_y0 = self.yrb + self.side_gap
        # self.rb_cut_x1 = self.rb_cut_x0 - self.cut_size * sin(θ_in_radians)
        # self.rb_cut_y1 = self.rb_cut_y0 + self.cut_size * cos(θ_in_radians)

        # # self.rb_cut = self.cell.add_line((self.rb_cut_x0, self.rb_cut_y0), (self.rb_cut_x1, self.rb_cut_y1))

        
        # self.rt_cut_x0 = self.xrt - self.side_gap
        # self.rt_cut_y0 = self.yrt
        # self.rt_cut_x1 = self.rt_cut_x0 - self.cut_size * cos(θ_in_radians)
        # self.rt_cut_y1 = self.rt_cut_y0 - self.cut_size * sin(θ_in_radians)

        # # self.rt_cut = self.cell.add_line((self.rt_cut_x0, self.rt_cut_y0), (self.rt_cut_x1, self.rt_cut_y1))


        # self.lt_cut_x0 = self.xlt
        # self.lt_cut_y0 = self.ylt - self.side_gap
        # self.lt_cut_x1 = self.lt_cut_x0 + self.cut_size * sin(θ_in_radians)
        # self.lt_cut_y1 = self.lt_cut_y0 - self.cut_size * cos(θ_in_radians)

        # # self.lt_cut = self.cell.add_line((self.lt_cut_x0, self.lt_cut_y0), (self.lt_cut_x1, self.lt_cut_y1))
        # self.lines = [
        #     ((self.lb_cut_x0, self.lb_cut_y0), (self.lb_cut_x1, self.lb_cut_y1)),
        #     ((self.rb_cut_x0, self.rb_cut_y0), (self.rb_cut_x1, self.rb_cut_y1)),
        #     ((self.rt_cut_x0, self.rt_cut_y0), (self.rt_cut_x1, self.rt_cut_y1)),
        #     ((self.lt_cut_x0, self.lt_cut_y0), (self.lt_cut_x1, self.lt_cut_y1)),
        # ]

    def draw(self, msp, mirror=False):
        """
        # Draw the unit cell into modelspace.

        # mirror=True mirrors the cell about its local vertical axis.
        # """
        for i in range(len(self.boundaryCuts)):
            line = self.boundaryCuts[i]
            sp = line[0]
            ep = line[1]
            msp.add_line((sp.x, sp.y), (ep.x, ep.y), dxfattribs={"layer": "BOUNDARY"})
        
        for i in range(len(self.internalCuts)):
            line = self.internalCuts[i]
            sp = line[0]
            ep = line[1]
            msp.add_line((sp.x, sp.y), (ep.x, ep.y), dxfattribs={"layer": "CUTS"})
        

        # for (p0, p1) in self.lines:
        #     x0, y0 = p0
        #     x1, y1 = p1

        #     if mirror:
        #         x0 = self.dx - x0
        #         x1 = self.dx - x1

        #     x0 += x_offset
        #     y0 += y_offset
        #     x1 += x_offset
        #     y1 += y_offset

        #     msp.add_line((x0, y0), (x1, y1), dxfattribs={"layer": "CUTS"})

        
def is_odd(n):
    return n % 2 == 1

if __name__ == "__main__":
    
    doc = ezdxf.new(dxfversion="R2010")
    doc.layers.new(name="CUTS", dxfattribs={"color": 1})
    doc.layers.new(name="BOUNDARY", dxfattribs={"color": 2})
    msp = doc.modelspace()

    s = 1
    edge_gap = 0.04
    t = 0.1
    θ = 10
    a = 0.60
    # θ_in_radians = np.radians(θ_in_degrees)

    # def __init__(self, s, t, θ, a):
    unit_cell = Unit_Cell(s, t, θ, a)
    
    # Generate pattern array
    n = 1  # number of repetitions in each direction
    for i in range(n):
        for j in range(n):
            # x = (i + is_odd(i+j)) * unit_cell.dx
            # y = j * unit_cell.dy
            unit_cell.draw(
                msp,
   
                mirror=is_odd(i + j),
            )
    # Draw lines around
    # msp.add_line((0, 0), (unit_cell.dx*n, 0), dxfattribs={"layer": "BOUNDARY"})
    # msp.add_line((unit_cell.dx*n, 0), (unit_cell.dx*n,unit_cell.dy*n), dxfattribs={"layer": "BOUNDARY"})
    # msp.add_line((unit_cell.dx*n, unit_cell.dy*n), (0, unit_cell.dy*n), dxfattribs={"layer": "BOUNDARY"})
    # msp.add_line((0, unit_cell.dy*n), (0, 0), dxfattribs={"layer": "BOUNDARY"})

    name = f"cell2.dxf"
    dir1 = "/home/odin/Academic/Research/Metamaterials/Auxetic_Mechanical_MM/Bistable/Computational_Exploration/Modules/DXF_Making/02_triangular"
    dir2 = "/home/odin/Academic/Research/Metamaterials/Auxetic_Mechanical_MM/Bistable/Computational_Exploration/CPP_Version/resources/DXFs"
    doc.saveas(f"{dir2}/{name}")
    doc.saveas(f"{dir1}/{name}")
    print(f"Pattern saved as {name}")
