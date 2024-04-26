import sympy
from tqdm import tqdm

class RiemannGeometry:
    """
    A class for calculating Riemannian geometry tensors including the Christoffel symbols,
    Riemann curvature tensor, Ricci tensor, and Ricci scalar from a given metric tensor.

    Attributes:
        g / metric          (sympy.Matrix): The metric tensor as a sympy Matrix. Represents the metric tensor g_{ab}.
        g_inv / metric_inv  (sympy.Matrix): The inverse of the metric tensor as a sympy Matrix. Represents the inverse metric tensor g^{ab}.
        x / coordinates     (sympy.Matrix): A sympy Matrix of coordinates used in the metric tensor.
        Christoffel_symbols (list of sympy.Matrix): A list of matrices where each matrix represents the Christoffel symbols Gamma^a_{bc} for a fixed upper index a.
        Riemann_tensor      (list of lists of sympy.Matrix): A nested list structure where each inner list contains matrices representing the Riemann curvature tensor components R^a_{bcd}.
        Ricci_tensor        (sympy.Matrix): A matrix representing the Ricci tensor R_{ab}, derived from contracting the Riemann tensor.
        Ricci_scalar        (sympy expression): A scalar representing the Ricci scalar R, derived from further contracting the Ricci tensor with the metric tensor.
    
    Convention:
        g / metric:          Denoted as g_{ab}, this is a symmetric tensor that describes the geometry of the manifold.
        g_inv / metric_inv:  Denoted as g^{ab}, this is the matrix inverse of the metric tensor.
        Christoffel_symbols: Denoted as Gamma^a_{bc}, these symbols represent the connection coefficients in the Levi-Civita connection derived from the metric.
        Riemann_tensor:      Denoted as R^a_{bcd}, it measures the non-commutativity of the covariant derivative and is computed based on the Christoffel symbols.
        Ricci_tensor:        Denoted as R_{ab}, this tensor is obtained by contracting the first and the third indices of the Riemann tensor.
        Ricci_scalar:        Denoted as R, this is obtained by contracting an index of the Ricci tensor with the inverse of the metric tensor.

    Methods:
    calculate(progress=True):
        Computes the Christoffel symbols, Riemann tensor, Ricci tensor, and Ricci scalar for the metric tensor provided during initialization.
        Optionally displays a progress bar if `progress` is set to True. If `progress` is False, the calculation proceeds without visual progress indication.
        This method efficiently organizes calculations to reuse results and minimize computational overhead.

    latex(markdown=False):
        Returns a string in LaTeX format representing all the calculated tensors. If `markdown` is True, the string is formatted to be used directly in markdown environments like Jupyter notebooks.
        This method utilizes the calculated tensor data to generate a formatted LaTeX representation, suitable for academic papers, reports, and presentations.
    
    Example (3D-sphere):
        >>> r, theta, phi = sympy.symbols('r theta phi')
        >>> R = sympy.symbols('R')
        >>> x = sympy.Matrix([r, theta, phi])
        >>> g = sympy.Matrix([
        >>>     [1, 0, 0],
        >>>     [0, R**2, 0],
        >>>     [0, 0, R**2 * (sympy.sin(theta))**2]
        >>> ])
        >>> geo = RiemannGeometry(g, x)
        >>> geo.calculate()
        >>> print(geo)
    """
    def __init__(self, g, x):
        """
        Initializes the RiemannGeometry with a metric tensor and coordinates.
        
        Parameters:
            g (sympy.Matrix): The metric tensor.
            x (sympy.Matrix): The coordinates used in the metric tensor.
        """
        self.n = len(x)
        self.g = g
        try:
            self.g_inv = g.inv()  # Attempt to calculate the inverse of the metric tensor
        except sympy.MatrixError as e:
            raise ValueError("The provided metric tensor g is not invertible.") from e
        self.x = x
        self.calculated = False
        self.Christoffel_symbols = None
        self.Riemann_tensor = None
        self.Ricci_tensor = None
        self.Ricci_scalar = None

        self.coordinate = x
        self.metric = self.g
        self.metric_inv = self.g_inv
        
    def calculate(self, progress=True):
        """
        Computes the Christoffel symbols, Riemann tensor, Ricci tensor, and Ricci scalar
        based on the metric tensor provided during initialization. This method
        facilitates the optional display of a progress bar to monitor the calculation progress.
    
        Parameters:
            progress (bool): If True, displays a progress bar during the calculation process.
                             If False, calculations proceed without progress indication.
    
        Returns:
            None. Updates class attributes with computed values.
    
        Raises:
            ValueError: If the tensors are already calculated and `self.calculated` is not reset.
    
        Example:
            >>> geo = RiemannGeometry(g, x)
            >>> geo.calculate(progress=True)  # Calculation with a progress bar.
            >>> geo.calculate(progress=False) # Calculation without a progress bar.
            >>> geo.calculate()               # Defaults to calculation with a progress bar.
            >>> geo.calculated = False        # Resetting to recalculate.
            >>> geo.calculate()               # Recalculates the tensors.
        """
        if self.calculated and self.Christoffel_symbols and self.Riemann_tensor: 
            raise ValueError("Tensors already calculated. To recalculate, set 'calculated' to False (e.g., object_name.calculated = False).")
    
        if progress:
            self.calculate_with_progress()
        else: 
            self.calculate_without_progress()

    
    def calculate_with_progress(self):
    
        self.Christoffel_symbols = []
        for a in tqdm(range(self.n), desc="Calculating Christoffel Symbols"):
            _Christoffel = self.function_to_sympyMatrix(lambda b,c: self.calculate_Christoffel(a, b, c))
            self.Christoffel_symbols.append(_Christoffel)
    
        self.Riemann_tensor = []
        for a in tqdm(range(self.n), desc="Calculating Riemann Tensor"):
            _Riemann_row = []
            for b in range(self.n):
                _Riemann = self.function_to_sympyMatrix(lambda c, d: self.calculate_Reimann(a, b, c, d))
                _Riemann_row.append(_Riemann)
            self.Riemann_tensor.append(_Riemann_row)
    
        self.Ricci_tensor = []
        for a in tqdm(range(self.n), desc="Calculating Ricci Tensor"):
            _Ricci_row = []
            for b in range(self.n):
                _Ricci = self.calculate_Ricci_tensor(a, b)
                _Ricci_row.append(_Ricci)
            self.Ricci_tensor.append(_Ricci_row)
        self.Ricci_tensor = sympy.Matrix(self.Ricci_tensor)
    
        self.Ricci_scalar = 0
        for a in tqdm(range(self.n), desc="Calculating Ricci Scalar"):
            for b in range(self.n):
                self.Ricci_scalar += self.g_inv[a, b] * self.get_Ricci_tensor(a, b)
    
        self.calculated = True
        
    def calculate_without_progress(self):

        self.Christoffel_symbols = []
        for a in range(self.n):
            _Christoffel = self.function_to_sympyMatrix(lambda b,c: self.calculate_Christoffel(a, b, c))
            self.Christoffel_symbols.append(_Christoffel)

        self.Riemann_tensor = []
        for a in range(self.n):
            _Riemann_row = []
            for b in range(self.n):
                _Riemann = self.function_to_sympyMatrix(lambda c, d: self.calculate_Reimann(a, b, c, d))
                _Riemann_row.append(_Riemann)
            self.Riemann_tensor.append(_Riemann_row)

        self.Ricci_tensor = []
        for a in range(self.n):
            _Ricci_row = []
            for b in range(self.n):
                _Ricci = self.calculate_Ricci_tensor(a,b)
                _Ricci_row.append(_Ricci)
            self.Ricci_tensor.append(_Ricci_row)
        self.Ricci_tensor = sympy.Matrix(self.Ricci_tensor)

        self.Ricci_scalar = 0
        for a in range(self.n):
            for b in range(self.n):
                self.Ricci_scalar += self.g_inv[a,b] * self.get_Ricci_tensor(a,b)
        
        self.calculated = True
    
    def function_to_sympyMatrix(self, m):
        if not hasattr(m, '__call__'):
            raise TypeError("Input must be a function.")

        matrixList = []
        for a in range(self.n):
            _row = []
            for b in range(self.n):
                _row.append(m(a,b))
            matrixList.append(_row)
        return sympy.Matrix(matrixList)

    def sympyMatrix_to_list(self, m):
        matrixList = []
        for a in range(self.n):
            _row = []
            for b in range(self.n):
                _row.append(m[a,b])
            matrixList.append(_row)
        return matrixList
        
    def calculate_Christoffel(self, a, b, c):
        Gamma = 0
        for d in range(self.n):
            Gamma += self.g_inv[a,d] * (g[d,c].diff(x[b]) + g[d,b].diff(x[c]) - g[b,c].diff(x[d]))
        return Gamma/2
    
    def calculate_Reimann(self, a, b, c, d):
        R = self.get_Christoffel(a,b,d).diff(x[c]) - self.get_Christoffel(a,b,c).diff(x[d])
        for e in range(self.n):
            R += self.get_Christoffel(a,e,c) * self.get_Christoffel(e,b,d) 
            R -= self.get_Christoffel(a,e,d) * self.get_Christoffel(e,b,c)
        return R

    def calculate_Ricci_tensor(self, a, b):
        R = 0
        for c in range(self.n):
            R += self.get_Riemann(c, a, c, b)
        return R

    def calculate_Ricci(self, a, b):
        R = 0
        for c in range(self.n):
            R += self.get_Riemann(c, a, c, b)
        return R
    
    # ------------------------------
    def get_Christoffel(self, a, b, c):
        return self.Christoffel_symbols[a][b,c]

    def get_Riemann(self, a, b, c, d):
        return self.Riemann_tensor[a][b][c,d]

    def get_Ricci_tensor(self, a, b):
        return self.Ricci_tensor[a,b]

    # ------------------------------
    def x_latex(self):
        x_list = [sympy.latex(self.x[i]) for i in range(self.n)]
        return r"\begin{pmatrix}"+r"\\".join(x_list)+"\\end{pmatrix}"
    
    def coordinate_latex(self):
        return self.x_latex()

    def g_latex(self):
        return self.matrix_to_latex(self.g)

    def metric_latex(self):
        return self.g_latex()

    def g_inv_latex(self):
        return self.matrix_to_latex(self.g_inv)

    def metric_inv_latex(self):
        return self.g_inv_latex()
        
    def Christoffel_symbols_latex(self):
        Gamma_latex = ' & '.join(map(self.matrix_to_latex, self.Christoffel_symbols))
        return f"\\begin{{pmatrix}} {Gamma_latex} \\end{{pmatrix}}"
        
    def Riemann_tensor_latex(self):
        matrixList = []
        for a in range(self.n):
            _row = []
            for b in range(self.n):
                _row.append(self.matrix_to_latex(self.Riemann_tensor[a][b]))
            matrixList.append(_row)
        return self.matrix_to_latex(matrixList)

    def Ricci_tensor_latex(self):
        return self.matrix_to_latex(self.Ricci_tensor)
    
    def Ricci_scalar_latex(self):
        return sympy.latex(self.Ricci_scalar)

    def matrix_to_latex(self, m):
        if isinstance(m, list):
            matrix_rows = [" & ".join(str(m[i][j]) for j in range(self.n)) for i in range(self.n)]
            matrix_latex = "\\\\ ".join(matrix_rows)
            return f"\\begin{{pmatrix}} {matrix_latex} \\end{{pmatrix}}"
        
        elif isinstance(m, sympy.Matrix):
            matrix_rows = [" & ".join(sympy.latex(m[i, j]) for j in range(self.n)) for i in range(self.n)]
            matrix_latex = "\\\\ ".join(matrix_rows)
            return f"\\begin{{pmatrix}} {matrix_latex} \\end{{pmatrix}}"
        
        else:
            raise TypeError("Unsupported type for LaTeX matrix generation.")

    def string(self):
        if not self.calculated:
            return r"Riemann geometry not calculated. Please call the calculate() method."
        coordinate_list = [self.x[i] for i in range(self.n)]
        metric_list = self.sympyMatrix_to_list(self.g)
        metric_inv_list = self.sympyMatrix_to_list(self.g_inv)
        christoffel_list  = [self.sympyMatrix_to_list(self.Christoffel_symbols[a]) for a in range(self.n)]
        riemann_list = [[self.sympyMatrix_to_list(self.Riemann_tensor[a][b]) for b in range(self.n)] for a in range(self.n)]
        ricci_list = self.sympyMatrix_to_list(self.Ricci_tensor)
        ricci = self.Ricci_scalar
        output = "# Riemann geometry\n"
        output += f"## x: {coordinate_list}\n"
        output += f"## metric: {metric_list}\n"
        output += f"## inverse metric: {metric_inv_list}\n"
        output += f"## Christoffel symbols: {christoffel_list}\n"
        output += f"## Riemann tensor: {riemann_list}\n"
        output += f"## Ricci tensor: {ricci_list}\n"
        output += f"## Ricci scalar: {ricci}"
        return output
    
    def latex(self, markdown=False):
        if not self.calculated:
            return r"Riemann geometry not calculated. Please call the calculate() method."
        coordinate_latex = self.x_latex()
        metic_latex = self.metric_latex()
        metic_inv_latex = self.metric_inv_latex()
        christoffel_latex = self.Christoffel_symbols_latex()
        riemann_latex = self.Riemann_tensor_latex()
        ricci_tensor_latex = self.Ricci_tensor_latex()
        ricci_latex = self.Ricci_scalar_latex()
        output = ""
        if markdown:
            output += f"$$x^{{a}}={coordinate_latex},\quad g_{{ab}}={metic_latex},\quad g^{{ab}}={metic_inv_latex},$$\n"
            output += f"$${{\\Gamma^{{a}}}}_{{bc}}={christoffel_latex},$$\n"
            output += f"$${{R^{{a}}}}_{{bcd}}={riemann_latex},$$\n"
            output += f"$$R_{{ab}}={ricci_tensor_latex},\quad R={ricci_latex}.$$"
            return output
        output += "\\begin{aligned}"
        output += f"x^{{a}}&={coordinate_latex},\quad g_{{ab}}={metic_latex},\quad g^{{ab}}={metic_inv_latex},\\\\ "
        output += f"{{\\Gamma^{{a}}}}_{{bc}}&={christoffel_latex},\\\\ "
        output += f"{{R^{{a}}}}_{{bcd}}&={riemann_latex},\\\\ "
        output += f"R_{{ab}}&={ricci_tensor_latex},\quad R={ricci_latex}."
        output += "\\end{aligned}"
        return output
    
    def __repr__(self):
        return self.string()

    def _repr_latex_(self):
        return self.latex(markdown=True)