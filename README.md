#include <iostream>

#include <vector>

#define DEBUG 0

class Vector {

public:

    unsigned int size = 0;


    // Constructors
    explicit Vector(int size){
        if (size > 0) this->size = size;
        else throw std::runtime_error("Vector can't have size lower then 1");

        elements.reserve(size);
    }

    explicit Vector(const std::vector<float>& vector){
        if (vector.empty()) throw std::runtime_error("Vector can't have size lower then 1");

        size = vector.size();
        elements = vector;
    }

    // Subscription operators and methods
    float &operator [] (const std::size_t index) {
        return elements[index];
    }

    const float &operator [] (const std::size_t index) const {
        return elements[index];
    }

    float &at(const std::size_t index) {
        return elements.at(index);
    }

    const float &at(const std::size_t index) const {
        return elements.at(index);
    }

    // +, -, *
    Vector operator + (const Vector &vector) const {
        if (vector.size != size) throw std::runtime_error("You can't add vectors with different sizes");

        std::vector<float> copy = elements;
        for (int i = 0; i < size; ++i)
            copy[i] += vector.elements[i];
        return Vector(copy);
    }

    Vector operator - (const Vector &vector) const {
        if (vector.size != size) throw std::runtime_error("You can't subtract vectors with different sizes");

        std::vector<float> copy = elements;
        for (int i = 0; i < size; ++i)
            copy[i] -= vector.elements[i];
        return Vector(copy);
    }

    float operator * (const Vector &vector) const {
        if (vector.size != size) throw std::runtime_error("You can't multiply vectors with different sizes");

        float sum = 0;
        for (int i = 0; i < size; ++i)
            sum += elements[i] * vector.elements[i];
        return sum;
    }

    // +=, -=, *=
    void operator += (const Vector &vector) {
        if (vector.size != size) throw std::runtime_error("You can't add vectors with different sizes");

        for (int i = 0; i < size; ++i)
            elements[i] += vector.elements[i];
    }

    void operator -= (const Vector &vector) {
        if (vector.size != size) throw std::runtime_error("You can't subtract vectors with different sizes");

        for (int i = 0; i < size; ++i)
            elements[i] -= vector.elements[i];
    }

    void operator *= (const float multiplier) {
        for (int i = 0; i < size; ++i)
            elements[i] *= multiplier;
    }

private:

    std::vector<float> elements;

};





class Matrix {

public:

    unsigned int rows = 0, columns = 0;


    // Constructors
    Matrix(int rows, int columns){
        if (rows > 0 and columns > 0) {
            this->rows = rows;
            this->columns = columns;
        }
        else throw std::runtime_error("Matrix must have size greater than 0");

        for (int row = 0; row < rows; ++row)
            matrix.emplace_back(columns);
    }

    explicit Matrix(std::vector<std::vector<float>> matrix){
        if (matrix.empty() or matrix[0].empty()) throw std::runtime_error("Matrix must have size greater than 0");

        rows = matrix.size();
        columns = matrix[0].size();

        for (const std::vector<float>& vector : matrix){
            if (vector.size() != columns) throw std::runtime_error("Incorrect matrix");
            this->matrix.emplace_back(vector);
        }
    }

    // Subscription operators and methods
    Vector &operator [] (const std::size_t index) {
        return matrix[index];
    }

    const Vector &operator [] (const std::size_t index) const {
        return matrix[index];
    }

    Vector &at(const std::size_t index) {
        return matrix.at(index);
    }

    const Vector &at(const std::size_t index) const {
        return matrix.at(index);
    }

    // +, -, *
    Matrix operator + (const Matrix &othermatrix) const {
        if (rows != othermatrix.rows or columns != othermatrix.columns) throw std::runtime_error("You can't add matrix with different sizes");

        Matrix result((int)rows, (int)columns);
        for (int row = 0; row < rows; ++row)
            for (int col = 0; col < columns; ++col)
                result[row][col] = matrix[row][col] + othermatrix[row][col];
        return result;
    }

    Matrix operator - (const Matrix &othermatrix) const {
        if (rows != othermatrix.rows or columns != othermatrix.columns) throw std::runtime_error("You can't add matrix with different sizes");

        Matrix result((int)rows, (int)columns);
        for (int row = 0; row < rows; ++row)
            for (int col = 0; col < columns; ++col)
                result[row][col] = matrix[row][col] - othermatrix[row][col];
        return result;
    }

    Matrix operator * (const float multiplier) const {
        Matrix result((int)rows, (int)columns);
        for (int row = 0; row < rows; ++row)
            for (int col = 0; col < columns; ++col)
                result[row][col] = matrix[row][col] * multiplier;
        return result;
    }

    // +=, -=, *=
    void operator += (const Matrix &othermatrix) {
        if (rows != othermatrix.rows or columns != othermatrix.columns) throw std::runtime_error("You can't add matrix with different sizes");

        for (int row = 0; row < rows; ++row)
            for (int col = 0; col < columns; ++col)
                matrix[row][col] += othermatrix[row][col];
    }

    void operator -= (const Matrix &othermatrix) {
        if (rows != othermatrix.rows or columns != othermatrix.columns) throw std::runtime_error("You can't add matrix with different sizes");

        for (int row = 0; row < rows; ++row)
            for (int col = 0; col < columns; ++col)
                matrix[row][col] -= othermatrix[row][col];
    }

    void operator *= (const float multiplier) {
        for (int row = 0; row < rows; ++row)
            for (int col = 0; col < columns; ++col)
                matrix[row][col] *= multiplier;
    }

    // Matrix methods
    Matrix transpose() const {
        Matrix transposed((int)columns, (int)rows);
        for (int row = 0; row < rows; ++row)
            for (int col = 0; col < columns; ++col)
                transposed[col][row] = matrix[row][col];
        return transposed;
    }

    Matrix dot(const Matrix &othermatrix) {
        if (rows != othermatrix.columns)
            throw std::runtime_error("You cannot calculate the product of matrices if the number of rows of the first matrix does not match the number of columns of the second matrix");

        Matrix result((int)rows, (int)othermatrix.columns);

        Matrix transposed = othermatrix.transpose();

        for (int row = 0; row < rows; ++row)
            for (int column = 0; column < othermatrix.columns; ++column)
                result.matrix[row][column] = matrix[row] * transposed.matrix[column];
        return result;
    }

    float determinant(int n) const {
        if (rows != columns) throw std::runtime_error("You can't find determinant of not square matrix");

        float determinant = 0;

        if (n == 1)
            return matrix[0][0];

        if (n == 2)
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

        Matrix temp(n, n);
        int sign = 1;
        for (int i = 0; i < n; i++) {
            temp = createSubMatrix(0, i, n);
            determinant += sign * matrix[0][i] * temp.determinant(n - 1);
            sign = -sign;
        }
        return determinant;
    }

    Matrix inverseMatrix() const {
        if (rows != columns) throw std::runtime_error("You can't find inverse of not square matrix");
        if (determinant((int)rows) == 0) throw std::runtime_error("You can't find inverse of matrix with determinant that equals 0");

        return *this * (1 / determinant((int)rows));
    }

    // Static methods
    static Matrix identityMatrix(int size) {
        Matrix matrix(size, size);
        for (int row = 0; row < size; ++row){
            for (int column = 0; column < size; ++column){
                matrix[row][column] = (row == column) ? 1 : 0;
            }
        }
        return matrix;
    }

    static void print(const Matrix& matrix){
        std::cout << std::endl;
        for (int row = 0; row < matrix.rows; ++row){
            for (int column = 0; column < matrix.columns; ++column){
                std::cout << matrix[row][column] << " ";
            }
            std::cout << std::endl;
        }
    }

private:

    std::vector<Vector> matrix;


    Matrix createSubMatrix(int p, int q, int n) const {
        if (rows != columns) throw std::runtime_error("You can't create submatrix for not square matrix");

        Matrix temp(n, n);

        int i = 0, j = 0;
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                if (row != p and col != q) {
                    temp[i][j++] = matrix[row][col];
                    if (j == n - 1) {
                        j = 0;
                        i++;
                    }
                }
            }
        }
        return temp;
    }

};


void tests() {
    Matrix matrix1 ({{0, 0, 1},
                     {0, 1, 0},
                     {1, 0, 0}});
    Matrix matrix2 ({{3, 2, 1},
                     {2, 3, 1},
                     {1, 2, 3}});

    try {
        Matrix::print(matrix1 + matrix2);
        Matrix::print(matrix1 - matrix2);
        Matrix::print(matrix1.dot(matrix2));

        matrix1 = matrix1.transpose();
        Matrix::print(matrix1);
        Matrix::print(matrix1 * 5);

        Matrix::print(matrix1.inverseMatrix());
        Matrix::print(Matrix::identityMatrix(3));
        Matrix::print(matrix1.dot(matrix1));
    }
    catch (std::runtime_error &e) {
        std::cout << e.what();
    }
}




int main(){
    if (DEBUG){
        tests();
        return 0;
    }

    Matrix *matrix1 = nullptr, *matrix2 = nullptr;
    int rows = 0, columns = 0;
    float element = 0;
    std::string answer;
    std::cout << "calculator matrix" << std::endl << std::endl<< std::endl;
    while (true){
        std::cout << "1. Enter first matrix" << std::endl;
        std::cout << "2. Enter second matrix" << std::endl;
        if (matrix1){
            std::cout << "3. Find determinant of first matrix" << std::endl;
            std::cout << "4. Transpose first matrix" << std::endl;
            std::cout << "5. Find inverse matrix of first matrix" << std::endl;
            std::cout << "6. Multiply first matrix by a number" << std::endl;
        }
        if (matrix2){
            std::cout << "7. Find determinant of second matrix" << std::endl;
            std::cout << "8. Transpose second matrix" << std::endl;
            std::cout << "9. Find inverse matrix of second matrix" << std::endl;
            std::cout << "10. Multiply second matrix by a number" << std::endl;
        }
        if (matrix1 and matrix2){
            std::cout << "11. Add matrices" << std::endl;
            std::cout << "12. Subtract matrices" << std::endl;
            std::cout << "13. Dot production of matrices" << std::endl;
        }
        std::cout << "14. Exit" << std::endl;
        std::cin >> answer;
        try {
            int ans = std::stoi(answer);
            if (ans == 1){
                if (matrix1)
                    delete matrix1;
                std::cout << "Enter matrix sizes" << std::endl;
                std::cin >> rows >> columns;
                std::cout << "Enter matrix elements" << std::endl;
                std::vector<std::vector<float>> elements;
                for (int r = 0; r < rows; r++){
                    elements.emplace_back();
                    for (int c = 0; c < columns; c++){
                        std::cin >> element;
                        elements[r].push_back(element);
                    }
                }
                matrix1 = new Matrix(elements);
            }
            else if (ans == 2){
                if (matrix2)
                    delete matrix2;
                std::cout << "Enter matrix sizes" << std::endl;
                std::cin >> rows >> columns;
                std::cout << "Enter matrix elements" << std::endl;
                std::vector<std::vector<float>> elements;
                elements.reserve(rows);
                for (int r = 0; r < rows; r++){
                    elements.emplace_back();
                    for (int c = 0; c < columns; c++){
                        std::cin >> element;
                        elements[r].push_back(element);
                    }
                }
                matrix2 = new Matrix(elements);
            }
            else if (ans == 3){
                std::cout << matrix1->determinant(rows) << std::endl;
            }
            else if (ans == 4){
                Matrix::print(matrix1->transpose());
            }
            else if (ans == 5){
                Matrix::print(matrix1->inverseMatrix());
            }
            else if (ans == 6){
                std::cout << "Enter multiplier" << std::endl;
                std::cin >> ans;
                Matrix::print(*matrix1 * ans);
            }
            else if (ans == 7){
                std::cout << matrix2->determinant(rows) << std::endl;
            }
            else if (ans == 8){
                Matrix::print(matrix2->transpose());
            }
            else if (ans == 9){
                Matrix::print(matrix2->inverseMatrix());
            }
            else if (ans == 10){
                std::cout << "Enter multiplier" << std::endl;
                std::cin >> ans;
                Matrix::print(*matrix2 * ans);
            }
            else if (ans == 11){
                Matrix::print(*matrix1 + *matrix2);
            }
            else if (ans == 12){
                Matrix::print(*matrix1 - *matrix2);
            }
            else if (ans == 13){
                Matrix::print(matrix1->dot(*matrix2));
            }
            else if (ans == 14){
                break;
            }
        }
        catch (const std::runtime_error &err) {
            std::cout << err.what() << std::endl;
            break;
        }
    }
    return 0;
}
