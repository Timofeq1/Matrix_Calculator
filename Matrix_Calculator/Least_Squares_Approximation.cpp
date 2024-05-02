// Timofey Ivlev
#include<iostream>
#include<vector>
#include<iomanip>
#include<cmath>
using namespace std;

class Matrix {
public:
    int row = 0;
    int column = 0;
    vector<vector<double>> table;
    Matrix(int rows, int columns) : row(rows), column(columns), table(rows, vector<double>(columns)) {
    };

    Matrix() = default;

    virtual ~Matrix() {}

    void swapRows(int a, int b) {
        swap(table[a], table[b]);
    }

    friend istream& operator>>(istream& in, Matrix& input) {
        in >> input.row;
        in >> input.column;
        input.table.resize(input.row, vector<double>(input.column));
        for (int i = 0; i < input.row; i++) {
            for (int j = 0; j < input.column; j++) {
                in >> input.table[i][j];
            }
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, Matrix& mat) {
        if (mat.row > 0 && mat.column > 0) {
            for (int i = 0; i < mat.row; i++) {
                for (int j = 0; j < mat.column; j++) {
                    out << fixed << setprecision(4) << mat.table[i][j] << " ";

                }
                out << endl;
            }
        }
        return out;
    }
    Matrix transpose() const {
        Matrix result(column, row);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                result.table[j][i] = table[i][j];
            }
        }
        return result;
    }
    virtual Matrix& operator=(Matrix& right) {
        row = right.row;
        column = right.column;
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                this->table[i][j] = right.table[i][j];
            }
        }
        return *this;
    }

    virtual Matrix operator+(Matrix& right) const {
        if (row != right.row || column != right.column) {
            cout << "Error: the dimensional problem occurred" << endl;
            return Matrix();
        }
        Matrix result(row, right.column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < right.column; j++) {
                result.table[i][j] = table[i][j] + right.table[i][j];
            }
        }
        return result;
    }

    virtual Matrix operator-(Matrix& right) const {
        if (row != right.row || column != right.column) {
            cout << "Error: the dimensional problem occurred" << endl;
            return Matrix();
        }
        Matrix result(row, right.column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < right.column; j++) {
                result.table[i][j] = table[i][j] - right.table[i][j];
            }
        }
        return result;
    }
    virtual Matrix operator*(Matrix& right) const {
        if (column != right.row) {
            cout << "Error: the dimensional problem occurred" << endl;
            return Matrix();
        }
        Matrix result(row, right.column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < right.column; j++) {
                double sum = 0;
                for (int r = 0; r < column; r++) {
                    double leftOp = table[i][r];
                    double rightOp = right.table[r][j];
                    sum += table[i][r] * right.table[r][j];
                }
                result.table[i][j] = sum;
            }
        }
        return result;
    }
};


class SquareMatrix : public Matrix {
public:

    int _step = 1;
    double _detSign = 1.0;
    SquareMatrix(int size) : Matrix(size, size) {}

    SquareMatrix() : Matrix() {}

    SquareMatrix(const Matrix& mat) : Matrix(mat) {}

    friend istream& operator>>(istream& in, SquareMatrix& input) {
        int size;
        in >> size;
        input.column = size;
        input.row = size;
        input.table.resize(size, vector<double>(size));
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                in >> input.table[i][j];
            }
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, const SquareMatrix& mat) {
        if (mat.row > 0 && mat.row > 0) {
            for (int i = 0; i < mat.row; i++) {
                for (int j = 0; j < mat.row; j++) {
                    out << fixed << setprecision(4) << mat.table[i][j] << " ";
                }
                out << endl;
            }
        }
        return out;
    }

    SquareMatrix operator+(SquareMatrix& right) {
        Matrix sum = Matrix::operator+(right);
        return static_cast<SquareMatrix>(sum);
    }
    SquareMatrix operator=(SquareMatrix& right) {
        Matrix newMat = Matrix::operator=(right);
        return static_cast<SquareMatrix>(newMat);
    }
    SquareMatrix operator*(SquareMatrix& right) {
        Matrix newMat = Matrix::operator*(right);
        return static_cast<SquareMatrix>(newMat);
    }
    SquareMatrix operator-(SquareMatrix& right) {
        Matrix newMat = Matrix::operator-(right);
        return static_cast<SquareMatrix>(newMat);
    }
    SquareMatrix transpose() {
        Matrix newMat = Matrix::transpose();
        return static_cast<SquareMatrix>(newMat);
    }

    int findPivot(int _row, int _size, bool print) {
        int pivotRow = _row;

        for (int j = _row + 1; j < _size; j++) {
            if (abs(this->table[j][_row]) > abs(this->table[pivotRow][_row])) {
                pivotRow = j;
            }
        }

        if (pivotRow != _row) {
            this->swapRows(_row, pivotRow);
            if (print) {
                cout << "step #" << _step << ": permutation" << endl;
                cout << *this;
            }
            _step++;
            _detSign *= -1.0;
        }
        return pivotRow;

    }

    double findSum() {
        double sum = 0.0;
        int size = this->column;
        for (int i = 1; i < size; i++) {
            for (int j = 0; j < i; j++) {
                sum += this->table[i][j];
            }
            if (sum != 0.0) { break; }
        }
        return sum;
    }

    double findDeterminant(bool print) {
        int size = this->column;
        int i = 0;
        int pivotRow = 0;

        double determinant = 1.0;

        double sum = findSum();

        while (i < column && sum != 0.0) {
            if (i + 1 >= size) { break; }

            pivotRow = findPivot(i, size, print);

            sum = findSum();
            if (sum == 0.0) { break; }

            for (int j = i + 1; j < size; j++) {
                double fac = this->table[j][i] / this->table[i][i];
                if (fac == 0) { continue; }
                for (int k = i; k < size; k++) {
                    if (fac > 0) { this->table[j][k] -= fac * this->table[i][k]; }
                    else {
                        this->table[j][k] += (fac * -1.0) * this->table[i][k];
                    }
                }
                if (print) {
                    cout << "step #" << _step << ": elimination" << endl;
                    cout << *this;
                }
                _step++;
                pivotRow = findPivot(i, size, print);
                sum = findSum();
                if (sum == 0.0) { break; }
            }
            i++;
        }

        for (int i = 0; i < size; i++) {
            determinant *= this->table[i][i];
        }
        determinant *= _detSign;
        if (print) {
            cout << "result:" << endl;
            cout << determinant << endl;
        }
        return determinant;
    }
};

class AugmentedMatrix : public Matrix {
public:
    int step = 1;

    friend ostream& operator<<(ostream& out, const AugmentedMatrix& mat) {
        if (mat.column == (2 * mat.row)) {
            for (int i = 0; i < mat.row; i++) {
                for (int j = 0; j < mat.column; j++) {
                    out << fixed << setprecision(4) << mat.table[i][j] << " ";
                }
                out << endl;
            }
            return out;
        }
        else {
            for (int i = 0; i < mat.row; i++) {
                for (int j = 0; j < mat.row; j++) {
                    out << fixed << setprecision(4) << mat.table[i][j] << " ";
                }
                out << endl;
            }
            for (int i = 0; i < mat.row; i++) {
                out << fixed << setprecision(4) << mat.table[i][mat.column - 1] << endl;
            }
            return out;
        }

    }

    AugmentedMatrix() : Matrix() {}
    AugmentedMatrix(SquareMatrix& mat) : Matrix(mat.row, mat.row * 2) {
        int size = mat.row;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                this->table[i][j] = mat.table[i][j];
            }

            for (int j = size - 1; j < size * 2; j++) {
                if (j - size == i) {
                    this->table[i][j] = 1.0;
                }
            }
        }
    }
    int findLowerPivot(int _row, int _size, bool print) {
        int pivotRow = _row;

        for (int j = _row + 1; j < _size; j++) {
            if (abs(this->table[j][_row]) > abs(this->table[pivotRow][_row])) {
                pivotRow = j;
            }
        }

        if (pivotRow != _row) {
            this->swapRows(_row, pivotRow);
            if (print) {
                cout << "step #" << step << ": permutation" << endl;
                cout << *this;
            }
            step++;
        }
        return pivotRow;
    }

    double findLowerSum() {
        double sum = 0.0;
        int size = this->row;
        for (int i = 1; i < size; i++) {
            for (int j = 0; j < i; j++) {
                sum += this->table[i][j];
            }
            if (sum != 0.0) { break; }
        }
        return sum;
    }

    void lowerElimination(bool print) {
        int size = this->row;
        int i = 0;
        int pivotRow = 0;

        double sum = findLowerSum();

        while (i < row && sum != 0.0) {
            if (i + 1 >= size) { break; }

            pivotRow = findLowerPivot(i, size, print);

            sum = findLowerSum();
            if (sum == 0.0) { break; }

            for (int j = i + 1; j < size; j++) {
                double fac = this->table[j][i] / this->table[i][i];
                if (fac == 0) { continue; }
                for (int k = i; k < column; k++) {
                    if (fac > 0) { this->table[j][k] -= fac * this->table[i][k]; }
                    else {
                        this->table[j][k] += (fac * -1.0) * this->table[i][k];
                    }
                }
                if (print) {
                    cout << "step #" << step << ": elimination" << endl;
                    cout << *this;
                }
                step++;
                pivotRow = findLowerPivot(i, size, print);
                sum = findLowerSum();
                if (sum == 0.0) { break; }
            }
            i++;
        }
    }

    double findUpperSum() {
        double sum = 0.0;
        int size = this->row; // size by row (matrix size)
        for (int i = 0; i < size - 1; i++) {
            for (int j = i + 1; j < size; j++) {
                sum += this->table[i][j];
                if (sum != 0.0) { return sum; }
            }

        }
        return sum;
    }
    void upperElimination(bool print) {

        int i = row - 1;

        double sum = findUpperSum();

        while (i >= 0 && sum != 0.0) {

            sum = findUpperSum();
            if (sum == 0.0) { break; }

            for (int j = i - 1; j >= 0; j--) {
                double fac = this->table[j][i] / this->table[i][i];
                for (int k = j + 1; k < column; k++) {
                    if (fac > 0) { this->table[j][k] -= fac * this->table[i][k]; }
                    else {
                        this->table[j][k] += (fac * -1.0) * this->table[i][k];
                    }
                }
                if (print) {
                    cout << "step #" << step << ": elimination" << endl;
                    cout << *this;
                }
                step++;

                sum = findUpperSum();
                if (sum == 0.0) { break; }
            }
            i--;
        }
    }

    void diagonalNormalization(bool print) {
        for (int i = 0; i < row; i++) {
            double fac = this->table[i][i];
            for (int j = i; j < column; j++) {
                this->table[i][j] /= fac;
            }
        }
        if (print) {
            cout << "Diagonal normalization:" << endl;
            cout << *this;
        }
    }

    void gaussianElimination(bool print) {
        lowerElimination(print);
        upperElimination(print);
        diagonalNormalization(print);
    }

    SquareMatrix returnInverse() {
        SquareMatrix out(row);

        for (int i = 0; i < row; i++) {
            for (int j = row; j < column; j++) {
                out.table[i][j - row] = this->table[i][j];
            }
        }
        return out;
    }
};


int main() {
    // reading length of dataset
    int m;
    cin >> m;
    vector<int> b_i(m, 0);
    vector<int> t_i(m, 0);
    // reading data to arrays
    for (int i = 0; i < m; i++) {
        int a;
        int b;
        cin >> a >> b;
        t_i[i] = a;
        b_i[i] = b;
    }
    // approximational polynomial
    int n;
    cin >> n;

    // A matrix
    Matrix A(m, 1 + n);
    // b column vector
    Matrix b(m, 1);

    // writing data to matrix A where each new column has degree of n
    for (int i = 0; i < m; i++) {
        A.table[i][0] = 1;
        for (int j = 1; j <= n; j++) {
            A.table[i][j] = pow(t_i[i], j);
        }
        b.table[i][0] = b_i[i];
    }

    cout << "A:" << endl;
    cout << A;

    Matrix A_transposed(A.transpose());
    Matrix A_t_A = (A_transposed * A);
    cout << "A_T*A:" << endl;
    cout << A_t_A;

    SquareMatrix A_t_A_inv(A_t_A);
    AugmentedMatrix Augmented(A_t_A_inv);
    Augmented.gaussianElimination(false);
    SquareMatrix A_t_A_inv2 = Augmented.returnInverse();
    cout << "(A_T*A)^-1:" << endl;
    cout << A_t_A_inv2;

    cout << "A_T*b:" << endl;
    Matrix A_transposed2(A.transpose());
    Matrix A_t_b = (A_transposed2 * b);
    cout << A_t_b;

    cout << "x~:" << endl;
    Matrix solution = ((Matrix)A_t_A_inv2 * A_t_b);
    cout << solution;

    return 0;
}