import ezdxf

def verify_dxf_file(filename):
    """Verify that a DXF file contains the expected geometry"""
    try:
        doc = ezdxf.readfile(filename)
        msp = doc.modelspace()
        
        print(f"\n=== Verifying {filename} ===")
        print(f"Entities in modelspace: {len(msp)}")
        
        # Check for block references
        block_refs = [e for e in msp if e.dxftype() == "INSERT"]
        print(f"Block references: {len(block_refs)}")
        
        # Check blocks
        blocks = list(doc.blocks)
        print(f"Total blocks: {len(blocks)}")
        
        for block in blocks:
            if not block.name.startswith('*'):  # Skip anonymous blocks
                entities = list(block)
                print(f"  Block '{block.name}': {len(entities)} entities")
                for entity in entities:
                    print(f"    - {entity.dxftype()}")
        
        # Check all entities in modelspace
        print("Entities in modelspace:")
        for entity in msp:
            print(f"  - {entity.dxftype()}")
            
        return True
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return False

if __name__ == "__main__":
    verify_dxf_file("unit_cell.dxf")
    verify_dxf_file("pattern.dxf")
