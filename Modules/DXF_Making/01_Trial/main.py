import ezdxf
import numpy as np
from numpy import cos, sin


# for i in range(n):
#     for j in range(n):
#         x, y = i * a, j * a
#         msp.add_circle(center=(x, y), radius=2.0)

# doc.saveas("pattern.dxf")

class Unit_Cell: 
    def __init__(self, doc, msp, name, cut_size, θ_in_radians, edge_gap, side_gap):
        self.doc = doc
        self.msp = msp
        self.cell = self.doc.blocks.new(name=name)
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
        
        # self.right_edge = self.cell.add_line((self.xrb, self.yrb), (self.xrt, self.yrt))
        # self.top_edge = self.cell.add_line((self.xlt, self.ylt), (self.xrt, self.yrt))
        # self.left_edge = self.cell.add_line((self.xlb, self.ylb), (self.xlt, self.ylt))
        # self.bottom_edge = self.cell.add_line((self.xlb, self.ylb), (self.xrb, self.yrb))
        
        self.lb_cut_x0 = self.xlb + self.side_gap
        self.lb_cut_y0 = self.ylb
        self.lb_cut_x1 = self.lb_cut_x0 + self.cut_size * cos(θ_in_radians)
        self.lb_cut_y1 = self.lb_cut_y0 + self.cut_size * sin(θ_in_radians)

        self.lb_cut = self.cell.add_line((self.lb_cut_x0, self.lb_cut_y0), (self.lb_cut_x1, self.lb_cut_y1))
        
        
        self.rb_cut_x0 = self.xrb
        self.rb_cut_y0 = self.yrb + self.side_gap
        self.rb_cut_x1 = self.rb_cut_x0 - self.cut_size * sin(θ_in_radians)
        self.rb_cut_y1 = self.rb_cut_y0 + self.cut_size * cos(θ_in_radians)

        self.rb_cut = self.cell.add_line((self.rb_cut_x0, self.rb_cut_y0), (self.rb_cut_x1, self.rb_cut_y1))

        
        self.rt_cut_x0 = self.xrt - self.side_gap
        self.rt_cut_y0 = self.yrt
        self.rt_cut_x1 = self.rt_cut_x0 - self.cut_size * cos(θ_in_radians)
        self.rt_cut_y1 = self.rt_cut_y0 - self.cut_size * sin(θ_in_radians)

        self.rt_cut = self.cell.add_line((self.rt_cut_x0, self.rt_cut_y0), (self.rt_cut_x1, self.rt_cut_y1))


        self.lt_cut_x0 = self.xlt
        self.lt_cut_y0 = self.ylt - self.side_gap
        self.lt_cut_x1 = self.lt_cut_x0 + self.cut_size * sin(θ_in_radians)
        self.lt_cut_y1 = self.lt_cut_y0 - self.cut_size * cos(θ_in_radians)

        self.lt_cut = self.cell.add_line((self.lt_cut_x0, self.lt_cut_y0), (self.lt_cut_x1, self.lt_cut_y1))

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
    θ_in_degrees = 60
    θ_in_radians = np.radians(θ_in_degrees)

    
    unit_cell = Unit_Cell(doc_unit, msp_unit, "UnitCell", cut_size, θ_in_radians, edge_gap, side_gap)
    
    msp_unit.add_blockref("UnitCell", insert=(0, 0))
    
    doc_unit.saveas("unit_cell.dxf")
    
    # Create pattern - create the block in pattern document as well
    pattern_unit_cell = Unit_Cell(doc_pattern, msp_pattern, "UnitCell", cut_size, θ_in_radians, edge_gap, side_gap)
    
    # Generate pattern array
    n = 5  # number of repetitions in each direction
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

    doc_pattern.saveas("pattern2.dxf")
