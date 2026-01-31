from Objects import Node, Vector


class Helper:
    @staticmethod
    def find_lo_between_nodes(node_a: Node, node_b: Node) -> float:
        dx = node_a.position.x - node_b.position.x
        dy = node_a.position.y - node_b.position.y
        return (dx ** 2 + dy ** 2) ** 0.5
