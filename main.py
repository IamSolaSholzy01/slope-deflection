import numpy as np

### Current Limitations
# 1. Supports fixed, hinge or roller supports at ends
# 2. Only supports point loads and uniform distributed loads
# 3. Only supports 1D members
# 4. Only supports 1D members with 1D loads
# 5. Only supports 1D members with 1D loads and 1D supports
# 6. Does not cater for free ended members


support_rotations = {
    'fixed': 0,
    'roller': 1,
    'hinge': 1,
    # 'free': 1
}

class Node:
    def __init__(self, x_coordinates, type_of_support, deflection: float = 0.0):
        self.x_coordinates = x_coordinates
        self.type_of_support = type_of_support
        self.rotation = None
        self.rotation_coefficient = support_rotations[type_of_support]
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

def get_rotations(members: list[Member]): 

    if (len(members) == 1):
        c = members[0].EI / members[0].length
        if (members[0].node_a.type_of_support == 'fixed' and members[0].node_b.type_of_support == 'fixed'):
            results = np.array([0, 0])
            members[0].node_a.set_rotation(results[0])
            members[0].node_b.set_rotation(results[1])
            return
        elif (members[0].node_a.type_of_support == 'fixed'):
            results = np.array([[0], [-members[0].fEMb / (2 * c)]])
            members[0].node_a.set_rotation(0)
            members[0].node_b.set_rotation(-members[0].fEMb / (2 * c))
            return
        elif (members[0].node_b.type_of_support == 'fixed'):
            results = np.array([[-members[0].fEMa / (2 * c)], [0]])
            members[0].node_a.set_rotation(-members[0].fEMa / (2 * c))
            members[0].node_b.set_rotation(0)
            return
        
        coeff_matrix = np.array([
            [members[0].node_a.rotation_coefficient * 2 * c, members[0].node_b.rotation_coefficient * c], 
            [members[0].node_a.rotation_coefficient * c, members[0].node_b.rotation_coefficient * 2 * c]
            ])
        res_matrix = np.array([[-members[0].fEMa], [-members[0].fEMb]])
        results = np.dot(np.linalg.inv(coeff_matrix), res_matrix)
        members[0].node_a.set_rotation(results[0][0])
        members[0].node_b.set_rotation(results[1][0])
    
    else:
        if (members[0].node_a.type_of_support == 'fixed' or members[-1].node_b.type_of_support == 'fixed'):
            if (members[0].node_a.type_of_support == 'fixed'):
                members[0].node_a.set_rotation(0)
                if (members[-1].node_b.type_of_support == 'fixed'):
                    members[-1].node_b.set_rotation(0)
                    matrix = np.zeros((len(members) - 1, len(members) - 1))
                    res_matrix = np.zeros((len(members) -  1, 1))

                    for i in range(len(members)):
                        member = members[i]
                        c = member.EI / member.length

                        if (i == 0):
                            matrix[i][i] += member.node_b.rotation_coefficient * 2 * c
                            res_matrix[i][0] += -member.fEMb

                        if (i == len(members) - 1):
                            matrix[i-1][i-1] += member.node_a.rotation_coefficient * 2 * c
                            res_matrix[i - 1][0] += -member.fEMa

                        if i != 0 and i != len(members) - 1:
                            print(i)
                            matrix[i - 1][i - 1] += member.node_a.rotation_coefficient * 2 * c
                            matrix[i - 1][i] += member.node_b.rotation_coefficient * c
                            matrix[i][i - 1] += member.node_a.rotation_coefficient * c 
                            matrix[i][i] = member.node_b.rotation_coefficient * 2 * c
                            res_matrix[i - 1][0] += -member.fEMa
                            res_matrix[i][0] += -member.fEMb

                    results = np.dot(np.linalg.inv(matrix), res_matrix)

                    for i in range(len(members)):
                        if (i == 0 or i == len(members) - 1):
                            continue
                        
                        members[i].node_a.set_rotation(results[i - 1][0])
                        members[i].node_b.set_rotation(results[i][0])

                else:
                    matrix = np.zeros((len(members), len(members)))
                    res_matrix = np.zeros((len(members), 1))

                    for i in range(len(members)):
                        member = members[i]
                        c = member.EI / member.length

                        if (i == 0):
                            matrix[i][i] += member.node_b.rotation_coefficient * 2 * c
                            res_matrix[i][0] += -member.fEMb

                        if (i != 0):
                            matrix[i][i] += member.node_a.rotation_coefficient * 2 * c
                            matrix[i][i - 1] += member.node_b.rotation_coefficient * c
                            matrix[i - 1][i] += member.node_a.rotation_coefficient * c
                            matrix[i - 1][i - 1] += member.node_b.rotation_coefficient * 2 * c
                            res_matrix[i][0] += -member.fEMa

                    print(matrix)
                    print(res_matrix)

                    results = np.dot(np.linalg.inv(matrix), res_matrix)

                    for i in range(len(members)):
                        if (i == 0):
                            continue
                        members[i].node_a.set_rotation(results[i][0])
                        members[i].node_b.set_rotation(results[i - 1][0])
            
            else:
                members[-1].node_b.set_rotation(0)
                matrix = np.zeros((len(members), len(members)))
                res_matrix = np.zeros((len(members), 1))

                for i in range(len(members)):
                    member = members[i]
                    c = member.EI / member.length

                    if (i == len(members) - 1):
                            matrix[i-1][i-1] += member.node_a.rotation_coefficient * 2 * c
                            res_matrix[i - 1][0] += -member.fEMa

                    if (i != 0):
                        matrix[i][i] += member.node_a.rotation_coefficient * 2 * c
                        matrix[i][i - 1] += member.node_b.rotation_coefficient * c
                        matrix[i - 1][i] += member.node_a.rotation_coefficient * c
                        matrix[i - 1][i - 1] += member.node_b.rotation_coefficient * 2 * c
                        res_matrix[i][0] += -member.fEMa

                print(matrix)
                print(res_matrix)

                results = np.dot(np.linalg.inv(matrix), res_matrix)

                for i in range(len(members)):
                    if (i == len(members) - 1):
                        continue
                    members[i].node_a.set_rotation(results[i][0])
                    members[i].node_b.set_rotation(results[i - 1][0])
            return

        matrix = np.zeros((len(members) + 1, len(members) + 1))
        res_matrix = np.zeros((len(members) + 1, 1))
        for m in range(len(members)):
            member = members[m]
            c = member.EI / member.length
            matrix[m][m] += member.node_a.rotation_coefficient * 2 * c
            matrix[m][m + 1] += member.node_b.rotation_coefficient * c
            matrix[m + 1][m] += member.node_a.rotation_coefficient * c
            matrix[m + 1][m + 1] += member.node_b.rotation_coefficient * 2 * c
            res_matrix[m][0] += -member.fEMa
            res_matrix[m + 1][0] += -member.fEMb

        results = np.dot(np.linalg.inv(matrix), res_matrix)
        for i in range(len(members)):
            members[i].node_a.set_rotation(results[i][0])
            members[i].node_b.set_rotation(results[i + 1][0])
            print(members[i])


def calculate_fem_for_member():
    nodes = [
        Node(x_coordinates=0, type_of_support='roller'),
        Node(x_coordinates=8, type_of_support='roller'),
        Node(x_coordinates=16, type_of_support='roller'),
        Node(x_coordinates=24, type_of_support='fixed')
    ]

    members = []
    loads = [
        UniformDistributedLoad(magnitude = 20, start_end_x_coordinates = [0, 16]),
        PointLoad(coordinates=20, magnitude = 60),
        # UniformDistributedLoad(magnitude=15, start_end_x_coordinates=[10, 14]),
        # PointLoad(coordinates=7, magnitude=80)
    ]

    number_of_members: int = len(nodes) - 1

    for i in range(number_of_members):
        members.append(Member(nodes[i], nodes[i + 1], nodes[i + 1].x_coordinates - nodes[i].x_coordinates, E=1, I=1))

    for i in range(len(members)):
        a = members[i].node_a.x_coordinates
        b = members[i].node_b.x_coordinates
        l = members[i].length

        for j in range(len(loads)):
            if isinstance(loads[j], PointLoad) and a <= loads[j].coordinates < b:
                point_load = loads[j]
                x1_distance = point_load.coordinates - a
                x2_distance = b - point_load.coordinates
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
    for member in members:
        print(member)

calculate_fem_for_member()

