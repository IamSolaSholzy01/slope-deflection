import numpy as np

class Node:
    def __init__(self, x_coordinates, type_of_support, rotation_coefficient: int, deflection: float = 0.0):
        self.x_coordinates = x_coordinates
        self.type_of_support = type_of_support
        self.rotation = None
        self.rotation_coefficient = rotation_coefficient
        self.deflection = deflection

    def set_rotation(self, rotation: float):
        self.rotation = rotation

class Member:
    def __init__(self, node_a: Node, node_b: Node, length: float, E: float = None, I: float = None, EI = None):
        self.node_a = node_a
        self.node_b = node_b
        self.fEMa = 0.0
        self.fEMb = 0.0
        self.length = length
        self.EI = EI if EI else E * I
    
    def set_fems(self, fEMa: float, fEMb: float):
        self.fEMa = fEMa
        self.fEMb = fEMb

    def get_moments(self):
        m1 = self.fEMa + (self.EI / self.length) * (2 * self.node_a.rotation + self.node_b.rotation) # + (3 * self.node_a.deflection / self.length)
        m2 = self.fEMb + (self.EI / self.length) * (self.node_a.rotation + 2 * self.node_b.rotation) # - (3 * self.node_b.deflection / self.length)
        return m1, m2
    
    def __str__(self) -> str:
        return f'Member: {self.node_a.x_coordinates} -> {self.node_b.x_coordinates}, \
            Length: {self.length}, EI: {self.EI}, FEMa: {self.fEMa}, FEMb: {self.fEMb}, \
                Node A Rotation: {self.node_a.rotation}, Node B Rotation: {self.node_b.rotation}\
                    1st Moment: {self.get_moments()[0]}, 2nd Moment: {self.get_moments()[1]}'
    
class UniformDistributedLoad:
    def __init__(self, magnitude: float, start_end_x_coordinates: list[float]):
        self.magnitude = magnitude
        self.start_end_x_coordinates = start_end_x_coordinates

class PointLoad:
    def __init__(self, coordinates: float, magnitude: float):
        self.coordinates = coordinates
        self.magnitude = magnitude

support_rotations = {
    'fixed': 0,
    'roller': 1
}

def get_rotations(members: list[Member]): 
    if (len(members) == 1):
        c = members[0].EI / members[0].length
        print(c)
        coeff_matrix = np.array([[2 * c, c], [c, 2 * c]])
        res_matrix = np.array([[-members[0].fEMa], [-members[0].fEMb]])
        print(coeff_matrix)
        print(res_matrix)
        results = np.dot(np.linalg.inv(coeff_matrix), res_matrix)
        print(results)
        members[0].node_a.set_rotation(results[0][0])
        members[0].node_b.set_rotation(results[1][0])
        print(members[0])
    
    else:
        matrix = np.zeros((len(members) + 1, len(members) + 1))
        res_matrix = np.zeros((len(members) + 1, 1))
        for m in range(len(members)):
            member = members[m]
            c = member.EI / member.length
            matrix[m][m] += 2 * c
            matrix[m][m + 1] += c
            matrix[m + 1][m] += c
            matrix[m + 1][m + 1] += 2 * c
            res_matrix[m][0] += -member.fEMa
            res_matrix[m + 1][0] += -member.fEMb

        results = np.dot(np.linalg.inv(matrix), res_matrix)
        for i in range(len(members)):
            members[i].node_a.set_rotation(results[i][0])
            members[i].node_b.set_rotation(results[i + 1][0])
            print(members[i])


def calculate_fem_for_member():
    nodes = [
        Node(x_coordinates=0, type_of_support='roller', rotation_coefficient=support_rotations['roller']),
        # Node(x_coordinates=10, type_of_support='roller', rotation_coefficient=support_rotations['roller']),
        Node(x_coordinates=20, type_of_support='roller', rotation_coefficient=support_rotations['roller']),
        # Node(x_coordinates=30, type_of_support='roller', rotation_coefficient=support_rotations['roller'])
    ]

    members = []
    loads = [
        UniformDistributedLoad(magnitude = 20, start_end_x_coordinates = [2, 8]),
        PointLoad(coordinates=10, magnitude=10),
        UniformDistributedLoad(magnitude=15, start_end_x_coordinates=[10, 14]),
        PointLoad(coordinates=7, magnitude=80)
    ]

    number_of_members: int = len(nodes) - 1

    for i in range(number_of_members):
        members.append(Member(nodes[i], nodes[i + 1], nodes[i + 1].x_coordinates - nodes[i].x_coordinates, E=1, I=1))

    for i in range(len(members)):
        a = members[i].node_a.x_coordinates
        b = members[i].node_b.x_coordinates
        l = b - a

        for j in range(len(loads)):
            print(j)
            if isinstance(loads[j], PointLoad) and a <= loads[j].coordinates < b:
                point_load = loads[j]
                x1_distance = point_load.coordinates - a
                x2_distance = b - (point_load.coordinates - a)
                members[i].fEMa -= point_load.magnitude * x2_distance ** 2 * x1_distance / l ** 2
                members[i].fEMb += point_load.magnitude * x1_distance ** 2 * x2_distance / l ** 2

            elif isinstance(loads[j], UniformDistributedLoad):
                if loads[j].start_end_x_coordinates[0] <= a and \
                    loads[j].start_end_x_coordinates[1] >= b:
                    udl = loads[j]
                    members[i].fEMa -= udl.magnitude * l ** 2 / 12
                    members[i].fEMb += udl.magnitude * l ** 2 / 12

                elif loads[j].start_end_x_coordinates[0] <= a and loads[j].start_end_x_coordinates[1] < b:
                    udl = loads[j]
                    length = loads[j].start_end_x_coordinates[1] - a
                    members[i].fEMa -= (udl.magnitude * length ** 2) * ((3 * (length ** 2)) - (8 * length * l) + (6 * (l ** 2))) / (12 * l ** 2)
                    members[i].fEMb += (udl.magnitude * length ** 3) * ((4 * l) - (3 * length)) / (12 * l ** 2)

                elif loads[j].start_end_x_coordinates[0] > a and loads[j].start_end_x_coordinates[1] >= b:
                    udl = loads[j]
                    length = b - loads[j].start_end_x_coordinates[0]
                    members[i].fEMb -= (udl.magnitude * length ** 2) * ((3 * (length ** 2)) - (8 * length * l) + (6 * (l ** 2))) / (12 * l ** 2)
                    members[i].fEMa += (udl.magnitude * length ** 3) * ((4 * l) - (3 * length)) / (12 * l ** 2)

                elif loads[j].start_end_x_coordinates[0] > a and loads[j].start_end_x_coordinates[1] < b:
                    udl = loads[j]
                    length = loads[j].start_end_x_coordinates[1] - loads[j].start_end_x_coordinates[0]
                    l1 = loads[j].start_end_x_coordinates[1] - a
                    l2 = b - loads[j].start_end_x_coordinates[0]
                    members[i].fEMa -= (udl.magnitude * l1 ** 2) * ((3 * (l1 ** 2)) - (8 * l1 * l) + (6 * (l ** 2))) / (12 * l ** 2) + (udl.magnitude * l2 ** 3) * ((4 * l) - (3 * l2)) / (12 * l ** 2) - (udl.magnitude * l ** 2 / 12)
                    members[i].fEMb += (udl.magnitude * l1 ** 3) * ((4 * l) - (3 * l1)) / (12 * l ** 2) + (udl.magnitude * l2 ** 2) * ((3 * (l2 ** 2)) - (8 * l2 * l) + (6 * (l ** 2))) / (12 * l ** 2) - (udl.magnitude * l ** 2 / 12)

        print(f'FEM 1st: {members[i].fEMa}')
        print(f'FEM 2nd: {members[i].fEMb}')
        members[i].set_fems(members[i].fEMa, members[i].fEMb)

    get_rotations(members)

calculate_fem_for_member()

