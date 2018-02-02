from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):
    """Class describe the main info of movie
    The __init__ method may be documented in either the class level
    docstring, or as a docstring on the __init__ method itself.
    Either form is acceptable, but the two should not be mixed. Choose one
    convention to document the __init__ method and be consistent with it.
    Note:
        Do not include the `self` parameter in the ``Args`` section.
    Args:
        planes(plane): Planes
    Attributes:
        dimension: Linsys's dimension
    Methods:
        swap_rows(row1, row2): swap row1 and row2
        multiply_coefficient_and_row(coefficient, row): Multiply row with coefficient
        add_multiple_times_row_to_row(coefficient, row_to_add, row_to_be_added_to):
                         Add coefficient times of row_to_add to row_to_be_added_to
        compute_triangular_form(): change self to trianular form
        compute_rref(): Reduced row echelon form
        compute_solution(): compute solution
        compute_solution_parametrize(): Compute parametrize solution
    """    
    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):  # swap row1 and row2
        self[row1], self[row2] = self[row2], self[row1]


    def multiply_coefficient_and_row(self, coefficient, row):
        # Multiply row with coefficient
        n = self[row].normal_vector
        k = self[row].constant_term
        new_normal_vector = n.times_scalar(coefficient)
        new_constant_term = k * coefficient
        self[row] = Plane(normal_vector = new_normal_vector,
                          constant_term = new_constant_term)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        # Add coefficient times of row_to_add to row_to_be_added_to
        n1 = self[row_to_add].normal_vector
        n2 = self[row_to_be_added_to].normal_vector
        k1 = self[row_to_add].constant_term
        k2 = self[row_to_be_added_to].constant_term
        new_normal_vector = n1.times_scalar(coefficient).plus(n2)
        new_constant_term = k1 * coefficient + k2
        self[row_to_be_added_to] = Plane(normal_vector = new_normal_vector,
                                         constant_term = new_constant_term)

    def compute_triangular_form(self):  # change self to trianular form
        system = deepcopy(self)
        num_equations = len(system)
        num_variables = system.dimension
        j = 0
        for i in range(num_equations):
            while j < num_variables:
                c = MyDecimal(system[i].normal_vector.coordinates[j])
                if c.is_near_zero():
                    swap_succeeded = system.swap_with_row_below_for_nonzero_coefficient_if_able(i,j)
                    if not swap_succeeded:
                        j += 1
                        continue
                system.clear_coefficients_below(i,j)
                j += 1
                break
        ''' MY WAY TO ACHIEVE
        for i in range(num_equations):
            tmp1 = system.indices_of_first_nonzero_terms_in_each_row()
            if i != tmp1[i]:
                for j in range(i+1, len(system.planes)):
                    if tmp1[j] == i :
                        system.swap_rows(i,j)
                        break
            for j in range(i+1, len(system.planes)):
                xj = system[j].normal_vector.coordinates[i]
                xi = system[i].normal_vector.coordinates[i]
                system.add_multiple_times_row_to_row(-1 * xj / xi, i, j)
        '''
        return system

    def swap_with_row_below_for_nonzero_coefficient_if_able(self,row,col):
        num_equations = len(self)
        for k in range(row+1, num_equations):
            coefficient = MyDecimal(self[k].normal_vector.coordinates[col])
            if not coefficient.is_near_zero():
                self.swap_rows(row, k)
                return True
        return False

    def clear_coefficients_below(self, row, col):
        num_equations = len(self)
        beta = MyDecimal(self[row].normal_vector.coordinates[col])
        for k in range(row+1, num_equations):
            n = self[k].normal_vector
            gamma = n.coordinates[col]
            alpha = -gamma/beta
            self.add_multiple_times_row_to_row(alpha,row,k)

    def clear_coefficients_above(self, row, col):
        for k in range(row)[::-1]:
            n = self[k].normal_vector
            alpha = -1 * (n.coordinates[col])
            self.add_multiple_times_row_to_row(alpha,row,k)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension
        indices = [-1] * num_equations
        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e
        return indices

    def compute_rref(self):  # Reduced row echelon form
        tf = self.compute_triangular_form()
        num_equations = len(tf)
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()
        for i in range(num_equations)[::-1]:
            j = pivot_indices[i]
            if j < 0:
                continue
            tf.scale_row_to_make_coefficient_equal_one(i, j)
            tf.clear_coefficients_above(i, j)
        return tf

    def scale_row_to_make_coefficient_equal_one(self, row, col):
        n = self[row].normal_vector
        beta = Decimal('1.0') / n.coordinates[col]
        self.multiply_coefficient_and_row(beta, row)

    def compute_solution(self):  # compute solution
        try:
            return self.do_gaussian_elimination_and_extract_solution()
        except Exception as e:
            if(str(e) == self.NO_SOLUTIONS_MSG or
               str(e) == self.INF_SOLUTIONS_MSG):
               return str(e)
            else:
               raise e

    def do_gaussian_elimination_and_extract_solution(self):
        rref = self.compute_rref()
        rref.raise_exception_if_cotradictory_equation()
        rref.raise_exception_if_too_few_pivots()
        num_variables = rref.dimension
        solution_coodinates = [rref.planes[i].constant_term for i in
                               range(num_variables)]
        return Vector(solution_coodinates)

    def raise_exception_if_cotradictory_equation(self):
        for p in self.planes:
            try:
                p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == 'No nonzero elements found':
                    constant_term = MyDecimal(p.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)
                else:
                    raise e

    def raise_exception_if_too_few_pivots(self):
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0  else 0 for index in pivot_indices])
        num_variables = self.dimension
        if num_pivots < num_variables:
            raise Exception(self.INF_SOLUTIONS_MSG)

    def compute_solution_parametrize(self):  # Compute parametrize solution
        try:
            return self.do_gaussian_elimination_and_parametrize_solution()
        except Exception as e:
            if(str(e) == self.NO_SOLUTIONS_MSG):
               return str(e)
            else:
               raise e

    def do_gaussian_elimination_and_parametrize_solution(self):
        rref = self.compute_rref()
        rref.raise_exception_if_cotradictory_equation()
        direction_vectors = rref.extract_direction_vectores_for_parametrization()
        basepoint = rref.extract_basepoint_for_parametrization()
        return Parametrization(basepoint, direction_vectors)

    def extract_direction_vectores_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        free_variable_indices = set(range(num_variables)) - set(pivot_indices)
        direction_vectors = []
        for free_var in free_variable_indices:
            vector_coords = [0] * num_variables
            vector_coords[free_var] = 1
            for i,p in enumerate(self.planes):
                pivot_var = pivot_indices[i]
                if pivot_var < 0:
                    break
                vector_coords[pivot_var] = ip.normal_vector[free_var]
            direction_vectors.append(Vector(vector_coords))

        return direction_vectors

    def extract_basepoint_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        basepoint_coords = [0] * num_variables
        for i,p in enumerate(self.planes):
            pivot_var = pivot_indices[i]
            if pivot_var < 0:
                break
            basepoint_coords[pivot_var] = p.constant_term
        return Vector(basepoint_coords)

    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps
