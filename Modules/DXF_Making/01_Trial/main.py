import ezdxf
import numpy as np
from numpy import cos, sin


class Unit_Cell: 
    def __init__(self, cut_size, θ_in_radians, edge_gap, side_gap):
        self.cut_size = cut_size
        self.θ_in_radians = θ_in_radians
        self.edge_gap = edge_gap
        self.side_gap = side_gap
        self.a = cut_size + edge_gap
        self.l = (self.side_gap + self.a * cos(θ_in_radians)) / sin(θ_in_radians)

        self.dx = self.side_gap + self.l * cos(θ_in_radians) + self.a * sin(θ_in_radians)
        self.dy = self.dx 
        
        self.xlb = 0
        self.ylb = 0
        self.xrb = self.dx
        self.yrb = 0
        self.xlt = 0
        self.xrt = self.xrb
        self.yrt = self.dy
        self.ylt = self.yrt

        self.lb_cut_x0 = self.xlb + self.side_gap
        self.lb_cut_y0 = self.ylb
        self.lb_cut_x1 = self.lb_cut_x0 + self.cut_size * cos(θ_in_radians)
        self.lb_cut_y1 = self.lb_cut_y0 + self.cut_size * sin(θ_in_radians)

        # self.lb_cut = self.cell.add_line((self.lb_cut_x0, self.lb_cut_y0), (self.lb_cut_x1, self.lb_cut_y1))
        
        
        self.rb_cut_x0 = self.xrb
        self.rb_cut_y0 = self.yrb + self.side_gap
        self.rb_cut_x1 = self.rb_cut_x0 - self.cut_size * sin(θ_in_radians)
        self.rb_cut_y1 = self.rb_cut_y0 + self.cut_size * cos(θ_in_radians)

        # self.rb_cut = self.cell.add_line((self.rb_cut_x0, self.rb_cut_y0), (self.rb_cut_x1, self.rb_cut_y1))

        
        self.rt_cut_x0 = self.xrt - self.side_gap
        self.rt_cut_y0 = self.yrt
        self.rt_cut_x1 = self.rt_cut_x0 - self.cut_size * cos(θ_in_radians)
        self.rt_cut_y1 = self.rt_cut_y0 - self.cut_size * sin(θ_in_radians)

        # self.rt_cut = self.cell.add_line((self.rt_cut_x0, self.rt_cut_y0), (self.rt_cut_x1, self.rt_cut_y1))


        self.lt_cut_x0 = self.xlt
        self.lt_cut_y0 = self.ylt - self.side_gap
        self.lt_cut_x1 = self.lt_cut_x0 + self.cut_size * sin(θ_in_radians)
        self.lt_cut_y1 = self.lt_cut_y0 - self.cut_size * cos(θ_in_radians)

        # self.lt_cut = self.cell.add_line((self.lt_cut_x0, self.lt_cut_y0), (self.lt_cut_x1, self.lt_cut_y1))
        self.lines = [
            ((self.lb_cut_x0, self.lb_cut_y0), (self.lb_cut_x1, self.lb_cut_y1)),
            ((self.rb_cut_x0, self.rb_cut_y0), (self.rb_cut_x1, self.rb_cut_y1)),
            ((self.rt_cut_x0, self.rt_cut_y0), (self.rt_cut_x1, self.rt_cut_y1)),
            ((self.lt_cut_x0, self.lt_cut_y0), (self.lt_cut_x1, self.lt_cut_y1)),
        ]

    def draw(self, msp, x_offset=0.0, y_offset=0.0, mirror=False):
        """
        Draw the unit cell into modelspace.

        mirror=True mirrors the cell about its local vertical axis.
        """

        for (p0, p1) in self.lines:
            x0, y0 = p0
            x1, y1 = p1

            if mirror:
                x0 = self.dx - x0
                x1 = self.dx - x1

            x0 += x_offset
            y0 += y_offset
            x1 += x_offset
            y1 += y_offset

            msp.add_line((x0, y0), (x1, y1), dxfattribs={"layer": "CUTS"})
        
def is_odd(n):
    return n % 2 == 1

if __name__ == "__main__":
    
    doc = ezdxf.new(dxfversion="R2010")
    doc.layers.new(name="CUTS", dxfattribs={"color": 1})
    doc.layers.new(name="BOUNDARY", dxfattribs={"color": 2})
    msp = doc.modelspace()

    cut_size = 1
    edge_gap = 0.02
    side_gap = 0.1
    θ_in_degrees = 65
    θ_in_radians = np.radians(θ_in_degrees)

    
    unit_cell = Unit_Cell(cut_size, θ_in_radians, edge_gap, side_gap)
    
    # Generate pattern array
    n = 1  # number of repetitions in each direction
    for i in range(n):
        for j in range(n):
            x = (i + is_odd(i+j)) * unit_cell.dx
            y = j * unit_cell.dy
            unit_cell.draw(
                msp,
                x_offset=x,
                y_offset=y,
                mirror=is_odd(i + j),
            )
    # Draw lines around
    msp.add_line((0, 0),(0, unit_cell.dy*n),  dxfattribs={"layer": "BOUNDARY"})
    msp.add_line((0, unit_cell.dy*n), (unit_cell.dx*n, unit_cell.dy*n), dxfattribs={"layer": "BOUNDARY"})
    msp.add_line((unit_cell.dx*n,unit_cell.dy*n),(unit_cell.dx*n, 0),  dxfattribs={"layer": "BOUNDARY"})
    msp.add_line((unit_cell.dx*n, 0), (0, 0), dxfattribs={"layer": "BOUNDARY"})

    name = f"cell1.dxf"
    dir1 = "/home/odin/Academic/Research/Metamaterials/Auxetic_Mechanical_MM/Bistable/Computational_Exploration/Modules/DXF_Making/01_Trial"
    dir2 = "/home/odin/Academic/Research/Metamaterials/Auxetic_Mechanical_MM/Bistable/Computational_Exploration/CPP_Version/resources/DXFs"
    doc.saveas(f"{dir2}/{name}")
    doc.saveas(f"{dir1}/{name}")
    print(f"Pattern saved as {name}")
