import ezdxf
import numpy as np
from numpy import cos, sin


class Unit_Cell: 
    def __init__(self, doc, msp, name, cut_size, θ_in_radians, edge_gap, side_gap):
        self.doc = doc
        self.msp = msp
        # self.cell = self.doc.blocks.new(name=name)
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

        msp.add_line((x0, y0), (x1, y1))
        
def is_odd(n):
    return n % 2 == 1

if __name__ == "__main__":
    
    doc_unit = ezdxf.new(dxfversion="R2010")
    msp_unit = doc_unit.modelspace()
    
    doc_pattern = ezdxf.new(dxfversion="R2010")
    msp_pattern = doc_pattern.modelspace()

    cut_size = 1
    edge_gap = 0.04
    side_gap = 0.1
    θ_in_degrees = 65
    θ_in_radians = np.radians(θ_in_degrees)

    
    unit_cell = Unit_Cell(doc_unit, msp_unit, "UnitCell", cut_size, θ_in_radians, edge_gap, side_gap)
    
    msp_unit.add_blockref("UnitCell", insert=(0, 0))
    
    # doc_unit.saveas("unit_cell.dxf")
    
    # Create pattern - create the block in pattern document as well
    pattern_unit_cell = Unit_Cell(doc_pattern, msp_pattern, "UnitCell", cut_size, θ_in_radians, edge_gap, side_gap)
    
    # Generate pattern array
    n = 1  # number of repetitions in each direction
    for i in range(n):
        for j in range(n):
            x = (i + is_odd(i+j)) * unit_cell.dx
            y = j * unit_cell.dy
            msp_pattern.add_blockref("UnitCell", insert=(x, y),    dxfattribs={
        "xscale": -1 if is_odd(i+j) else 1,   # mirror in Y-axis
        "yscale": 1,
    })
    # Draw lines around
    msp_pattern.add_line((0, 0), (unit_cell.dx*n, 0))
    msp_pattern.add_line((unit_cell.dx*n, 0), (unit_cell.dx*n,unit_cell.dy*n))
    msp_pattern.add_line((unit_cell.dx*n, unit_cell.dy*n), (0, unit_cell.dy*n))
    msp_pattern.add_line((0, unit_cell.dy*n), (0, 0))

    doc_pattern.saveas(f"/home/odin/Academic/Research/Metamaterials/Auxetic_Mechanical_MM/Bistable/Computational_Exploration/Physics-Engine_Python_V3_02_Varying-Stiffness/dxf_files/pattern_{n}x{n}_edge-gap_{edge_gap}.dxf")
    print(f"Pattern saved as pattern_{n}x{n}_edge-gap_{edge_gap}.dxf")
